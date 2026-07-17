#' Rename dataset_origin factor labels to standardized "FirstAuthor_Year" format
#'
#' @description
#' 将 Seurat 对象 metadata 中的 `dataset_origin` 列标准化为学界通用的
#' "第一作者_年份" 引用格式（如 "Habermann et al. 2020"）。
#' 由于这是永久性修改 metadata 的关键步骤，会影响所有下游分析和图表 legend，
#' 函数默认在修改前自动备份原始对象，并对未被映射覆盖的原始值给出警告，
#' 避免 `factor()` 重新赋值时静默产生 NA 导致数据丢失。
#'
#' @param seurat_obj A Seurat object，其 `@meta.data` 中需包含 `dataset_origin` 列。
#' @param old_levels Character vector. 需要被重命名的原始 factor levels。
#'   默认对应 LungAgingERV 项目的三个数据来源。
#' @param new_labels Character vector. 新的标准化标签，按位置与 `old_levels` 一一对应。
#'   长度必须与 `old_levels` 一致。
#' @param keep_unmapped Logical. 若为 `TRUE`，`dataset_origin` 中未出现在
#'   `old_levels` 里的原始值会被保留（原样保留原值），而不是被强制转为 `NA`。
#'   默认为 `FALSE`（保持与旧版本行为一致，但会先给出明确警告）。
#' @param backup Logical. 若为 `TRUE`（默认），在修改前保存一份 `.rds` 备份。
#' @param backup_path Character. 备份文件路径，仅在 `backup = TRUE` 时生效。
#'   默认: `"backup_before_rename_dataset_origin.rds"`。
#' @param verbose Logical. 若为 `TRUE`（默认），打印重命名后的摘要表。
#'
#' @return 修改后的 Seurat 对象，`dataset_origin` 已标准化。
#'
#' @details
#' 默认重命名方案（LungAgingERV 项目）：
#' \itemize{
#'   \item `Kaminski_2020_Adams` -> `Adams_2020`（Sci Adv, Kaminski lab, Yale）
#'   \item `Banovich_2019_Habermann` -> `Habermann_2020`（Sci Adv, Banovich/Kropski lab；
#'         注：年份已订正，实际发表于2020年7月，非2019）
#'   \item `Banovich_2023_Natri` -> `Natri_2024`（Nat Genet, Banovich lab, TGen；
#'         注：年份已订正，正式发表于2024年3月，非2023）
#' }
#'
#' 安全机制：
#' \itemize{
#'   \item 若 `old_levels` 中某个值在数据里不存在 -> 给出 warning（原有行为）。
#'   \item 若数据中存在未被 `old_levels` 覆盖的值 -> 给出 warning，
#'         并根据 `keep_unmapped` 决定是保留原值还是转为 NA。
#' }
#'
#' @examples
#' \dontrun{
#' # 使用默认映射
#' immune.combined <- rename_dataset_origin(immune.combined)
#'
#' # 自定义映射（用于其他项目，如 HeartERVmap / PH-ERVAtlas）
#' immune.combined <- rename_dataset_origin(
#'   immune.combined,
#'   old_levels = c("A", "B"),
#'   new_labels = c("X_2020", "Y_2021"),
#'   keep_unmapped = TRUE
#' )
#' }
#'
#' @export
rename_dataset_origin <- function(seurat_obj,
                                  old_levels = c("Kaminski_2020_Adams",
                                                 "Banovich_2019_Habermann",
                                                 "Banovich_2023_Natri"),
                                  new_labels = c("Adams_2020",
                                                 "Habermann_2020",
                                                 "Natri_2024"),
                                  keep_unmapped = FALSE,
                                  backup = TRUE,
                                  backup_path = "backup_before_rename_dataset_origin.rds",
                                  verbose = TRUE) {

  # --- 输入检查 ---
  if (!inherits(seurat_obj, "Seurat")) {
    stop("`seurat_obj` 必须是一个 Seurat 对象。")
  }
  if (!"dataset_origin" %in% colnames(seurat_obj@meta.data)) {
    stop("在 seurat_obj@meta.data 中未找到 `dataset_origin` 列。")
  }
  if (length(old_levels) != length(new_labels)) {
    stop("`old_levels` 和 `new_labels` 长度必须一致。")
  }

  current_levels <- unique(as.character(seurat_obj@meta.data$dataset_origin))

  # 正向检查：old_levels 里有没有在数据中缺失的
  missing_levels <- setdiff(old_levels, current_levels)
  if (length(missing_levels) > 0) {
    warning("以下 `old_levels` 未在 dataset_origin 中找到: ",
            paste(missing_levels, collapse = ", "))
  }

  # 反向检查：数据中是否存在未被 old_levels 覆盖的值
  unmapped <- setdiff(current_levels, old_levels)
  if (length(unmapped) > 0) {
    if (keep_unmapped) {
      warning("以下 dataset_origin 值未在 old_levels 中定义，将保留原值: ",
              paste(unmapped, collapse = ", "))
    } else {
      warning("以下 dataset_origin 值未在 old_levels 中定义，将变为 NA: ",
              paste(unmapped, collapse = ", "),
              "\n如需保留原值，请设置 keep_unmapped = TRUE")
    }
  }

  # --- 自动备份 ---
  if (backup) {
    if (verbose) message("正在保存备份至: ", backup_path)
    saveRDS(seurat_obj, backup_path)
  }

  # --- 重命名 ---
  original_char <- as.character(seurat_obj@meta.data$dataset_origin)

  if (keep_unmapped) {
    # 只替换 old_levels 中定义的值，其余保留原值
    idx <- match(original_char, old_levels)
    new_char <- ifelse(!is.na(idx), new_labels[idx], original_char)
    seurat_obj@meta.data$dataset_origin <- factor(new_char)
  } else {
    # 原有行为：未匹配的值转为 NA
    seurat_obj@meta.data$dataset_origin <- factor(
      original_char,
      levels = old_levels,
      labels = new_labels
    )
  }

  # --- 输出摘要 ---
  if (verbose) {
    message("dataset_origin 重命名完成：")
    print(table(seurat_obj@meta.data$dataset_origin, useNA = "ifany"))
  }

  return(seurat_obj)
}
