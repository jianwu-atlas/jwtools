#' 加载并校验 TE Family/Class 注释表
#'
#' 读取转座元件(TE)家族与类别的映射注释表(scName -> TE_family -> TE_class),
#' 并对数据质量进行全面诊断(重复检测、冲突检测、NA值检查、分布概览)。
#' 该函数适用于 所有 Atlas 等项目中
#' TE/ERV Family富集分析(compareCluster/超几何检验)前的注释表准备步骤。
#'
#' @param path 注释表文件路径(CSV或TSV均可, 自动根据后缀识别分隔符)。
#'   默认为 NULL, 此时使用项目统一维护的清洗后版本(见函数内部代码)。
#' @param required_cols 必须存在的列名, 默认为 c("scName", "TE_family", "TE_class")。
#' @param verbose 逻辑值, 是否打印诊断信息。默认为 TRUE。
#' @param check_conflict 逻辑值, 是否检查scName对应多个不同family/class的数据冲突情况。
#'   默认为 TRUE。
#'
#' @return 一个 tibble, 包含三列: scName, TE_family, TE_class。
#'   附加属性 attr(x, "class_distribution") 和 attr(x, "family_distribution")
#'   保存了分布统计, 便于设定富集分析阈值参考。
#'
#' @examples
#' \dontrun{
#' te_family_map <- load_te_family_annotation()
#' head(te_family_map)
#' attr(te_family_map, "class_distribution")
#'
#' deg_df_annotated <- all_TE_degs |>
#'   dplyr::left_join(te_family_map, by = c("gene" = "scName"))
#' }
#'
#' @importFrom readr read_csv read_tsv
#' @importFrom dplyr group_by summarise filter n_distinct count arrange desc select slice ungroup n all_of
#' @importFrom stringr str_detect
#' @importFrom rlang .data
#' @export
load_te_family_annotation <- function(
    path = NULL,
    required_cols = c("scName", "TE_family", "TE_class"),
    verbose = TRUE,
    check_conflict = TRUE
) {

  # ---- 0. 默认路径 (放在函数体内, 避免 \usage 行过长的 Rd NOTE) ----
  if (is.null(path)) {
    path <- paste0(
      "/Users/Jack/jian_wu/ann/1_human/01_TEs/08_scTEname/",
      "z6_36_TEfamily_CleanedAnnotation_1019.csv"
    )
  }

  # ---- 1. 文件存在性检查 ----
  if (!file.exists(path)) {
    stop("TE family/class 注释表不存在: ", path,
         "\n请检查路径是否正确, 或该文件是否已被移动/重命名。")
  }

  # ---- 2. 根据后缀自动选择读取方式 ----
  ext <- tolower(tools::file_ext(path))
  raw_tab <- switch(
    ext,
    "csv" = readr::read_csv(path, show_col_types = FALSE),
    "tsv" = readr::read_tsv(path, show_col_types = FALSE),
    stop("不支持的文件格式: .", ext, " (仅支持 .csv 或 .tsv)")
  )

  if (verbose) {
    message("原始文件读取成功: ", nrow(raw_tab), " 行 x ", ncol(raw_tab), " 列")
  }

  # ---- 3. 必要列存在性检查 ----
  missing_cols <- setdiff(required_cols, colnames(raw_tab))
  if (length(missing_cols) > 0) {
    stop("注释表缺少必要列: ", paste(missing_cols, collapse = ", "),
         "\n实际列名为: ", paste(colnames(raw_tab), collapse = ", "))
  }

  # ---- 4. 重复/冲突检测 (对已清洗表通常返回0冲突, 作为二次质检保险) ----
  if (check_conflict) {
    dup_check <- raw_tab |>
      dplyr::group_by(.data$scName) |>
      dplyr::summarise(
        n_rows = dplyr::n(),
        n_unique_family = dplyr::n_distinct(.data$TE_family),
        n_unique_class  = dplyr::n_distinct(.data$TE_class),
        .groups = "drop"
      ) |>
      dplyr::filter(.data$n_rows > 1)

    n_conflict <- sum(dup_check$n_unique_family > 1 | dup_check$n_unique_class > 1)

    if (verbose) {
      message("重复/冲突诊断: 重复scName数 = ", nrow(dup_check),
              " | 真实冲突数 = ", n_conflict)
    }

    if (n_conflict > 0) {
      warning("发现 ", n_conflict, " 个 scName 存在 family/class 冲突, ",
              "建议检查原始注释源文件。当前优先保留含'-'的细粒度subfamily名称。")
    }
  }

  # ---- 5. 去重: 优先保留细粒度(含"-")的family名称, 确保结果可复现 ----
  te_family_map <- raw_tab |>
    dplyr::group_by(.data$scName) |>
    dplyr::arrange(.data$scName, dplyr::desc(stringr::str_detect(.data$TE_family, "-")),
                   .by_group = TRUE) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::select(dplyr::all_of(required_cols))

  # ---- 6. NA值检查 ----
  n_na_family <- sum(is.na(te_family_map$TE_family) | te_family_map$TE_family == "")
  n_na_class  <- sum(is.na(te_family_map$TE_class)  | te_family_map$TE_class  == "")

  if (verbose && (n_na_family > 0 | n_na_class > 0)) {
    message("注释表存在空值: TE_family为空 = ", n_na_family,
            " 条 | TE_class为空 = ", n_na_class, " 条")
  }

  # ---- 7. 分布概览 (附加为attribute, 不污染console, 需要时用attr()调取) ----
  class_dist <- te_family_map |>
    dplyr::count(.data$TE_class, name = "n_loci") |>
    dplyr::arrange(dplyr::desc(.data$n_loci))

  family_dist <- te_family_map |>
    dplyr::count(.data$TE_family, name = "n_loci") |>
    dplyr::arrange(dplyr::desc(.data$n_loci))

  attr(te_family_map, "class_distribution")  <- class_dist
  attr(te_family_map, "family_distribution") <- family_dist

  if (verbose) {
    message("最终映射表: ", nrow(te_family_map), " 个唯一 scName, ",
            dplyr::n_distinct(te_family_map$TE_family), " 个 family, ",
            dplyr::n_distinct(te_family_map$TE_class), " 个 class")
    message("(使用 attr(x, 'class_distribution') / attr(x, 'family_distribution') 查看分布概览)")
  }

  te_family_map
}
