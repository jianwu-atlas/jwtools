#' 保存当前环境所有变量到一个 .qs2 文件
#'
#' 对应你原来用的 \code{save(list = ls(all = TRUE), file = "xxx.RData")}，
#' 但用 \pkg{qs2} 压缩/读写，速度更快、体积更小，尤其适合装着 Seurat 对象
#' 的大型工作环境。
#'
#' @param file 输出文件路径，建议以 \code{.qs2} 结尾。
#' @param envir 要保存的环境，默认 \code{.GlobalEnv}（你当前的工作环境）。
#' @param nthreads 压缩/写入线程数，默认 6。
#' @param exclude 字符向量，按变量名精确排除（例如已经单独存过的大对象，
#'   避免重复占用磁盘）。默认不排除。
#' @param exclude_pattern 正则表达式字符串，按变量名模式排除（例如
#'   \code{"^p_|^tmp_"} 排除所有以 \code{p_} 或 \code{tmp_} 开头的变量）。
#'   默认不排除。
#' @param verbose 是否打印排除信息、变量大小排行和保存结果。默认 \code{TRUE}。
#'
#' @return 不可见地返回被保存的变量名字符向量（\code{invisible()}）。
#' @export
#'
#' @examples
#' \dontrun{
#' # 保存整个工作环境
#' qs_save_workspace("core_WT_YMO_workspace.qs2", nthreads = 14)
#'
#' # 排除已经单独存过的大对象，避免重复保存
#' qs_save_workspace(
#'   "core_WT_YMO_workspace.qs2",
#'   nthreads = 14,
#'   exclude = c("immune.combined")
#' )
#'
#' # 按正则排除所有 ggplot / 临时对象
#' qs_save_workspace(
#'   "core_WT_YMO_workspace.qs2",
#'   nthreads = 14,
#'   exclude_pattern = "^p_|^tmp_|^fig[0-9]"
#' )
#' }
qs_save_workspace <- function(
    file,
    envir = .GlobalEnv,
    nthreads = 6,
    exclude = NULL,
    exclude_pattern = NULL,
    verbose = TRUE
) {

  if (!requireNamespace("qs2", quietly = TRUE)) {
    stop("需要先安装 qs2 包: install.packages(\"qs2\")", call. = FALSE)
  }

  # 1. 获取当前环境所有变量名（包括隐藏变量，对应 ls(all=TRUE)）
  all_names <- ls(envir = envir, all.names = TRUE)

  # 2. 按变量名排除
  if (!is.null(exclude)) {
    removed_by_name <- intersect(all_names, exclude)
    all_names <- setdiff(all_names, exclude)
    if (verbose && length(removed_by_name) > 0) {
      message("排除变量(按名称): ", paste(removed_by_name, collapse = ", "))
    }
  }

  # 3. 按正则表达式排除
  if (!is.null(exclude_pattern)) {
    match_idx <- grepl(exclude_pattern, all_names)
    if (verbose && any(match_idx)) {
      message("排除变量(按正则'", exclude_pattern, "'): ",
              paste(all_names[match_idx], collapse = ", "))
    }
    all_names <- all_names[!match_idx]
  }

  if (length(all_names) == 0) {
    stop("没有可保存的变量", call. = FALSE)
  }

  # 4. 打包成list
  objs <- mget(all_names, envir = envir)

  # 5. 保存前展示大小信息（可选，方便确认）
  if (verbose) {
    obj_sizes_mb <- sapply(objs, function(x) as.numeric(utils::object.size(x)) / 1e6)
    size_df <- data.frame(
      name = names(objs),
      size_MB = round(obj_sizes_mb, 2)
    )
    size_df <- size_df[order(-size_df$size_MB), ]
    message("▶ 即将保存 ", length(objs), " 个变量，总大小约 ",
            round(sum(obj_sizes_mb), 1), " MB")
    print(utils::head(size_df, 10))  # 只显示最大的10个
  }

  # 6. 用qs2保存
  qs2::qs_save(objs, file, nthreads = nthreads)

  file_size_mb <- file.size(file) / 1e6
  message("✅ 已保存到: ", file)
  message("   文件大小: ", round(file_size_mb, 1), " MB (",
          round(file_size_mb / 1024, 2), " GB)")

  invisible(names(objs))
}


#' 从 .qs2 文件恢复所有变量到当前环境
#'
#' 对应你原来用的 \code{load("xxx.RData")}，用来读回
#' \code{\link{qs_save_workspace}} 保存的多变量文件。
#'
#' @param file 要读取的 \code{.qs2} 文件路径。
#' @param nthreads 读取线程数，默认 6。
#' @param envir 变量要恢复到哪个环境，默认 \code{.GlobalEnv}。
#' @param verbose 是否打印恢复了哪些变量。默认 \code{TRUE}。
#'
#' @return 不可见地返回被恢复的变量名字符向量（\code{invisible()}）。
#' @export
#'
#' @examples
#' \dontrun{
#' qs_load_workspace("core_WT_YMO_workspace.qs2", nthreads = 14)
#' }
qs_load_workspace <- function(
    file,
    nthreads = 6,
    envir = .GlobalEnv,
    verbose = TRUE
) {

  if (!requireNamespace("qs2", quietly = TRUE)) {
    stop("需要先安装 qs2 包: install.packages(\"qs2\")", call. = FALSE)
  }

  if (!file.exists(file)) {
    stop("文件不存在: ", file, call. = FALSE)
  }

  objs <- qs2::qs_read(file, nthreads = nthreads)

  if (!is.list(objs)) {
    stop("该文件不是通过 qs_save_workspace() 保存的多变量格式，",
         "请用 qs2::qs_read() 直接读取单个对象", call. = FALSE)
  }

  list2env(objs, envir = envir)

  if (verbose) {
    message("✅ 已恢复 ", length(objs), " 个变量到环境: ",
            paste(names(objs), collapse = ", "))
  }

  invisible(names(objs))
}
