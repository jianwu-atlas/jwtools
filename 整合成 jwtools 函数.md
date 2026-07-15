## 整合成 jwtools 函数

以下是可以直接放入您 `jwtools` 包的完整函数，遵循 roxygen2 文档规范，方便后续 `devtools::document()` 自动生成帮助文档：

```r
#' Rename dataset_origin factor labels to standardized "FirstAuthor_Year" format
#'
#' @description
#' 将 Seurat 对象 metadata 中的 `dataset_origin` 列标准化为学界通用的
#' "第一作者_年份" 引用格式（如 "Habermann et al. 2020"）。
#' 由于这是永久性修改 metadata 的关键步骤，会影响所有下游分析和图表 legend，
#' 函数默认在修改前自动备份原始对象。
#'
#' @param seurat_obj A Seurat object，其 `@meta.data` 中需包含 `dataset_origin` 列。
#' @param old_levels Character vector. 需要被重命名的原始 factor levels。
#'   默认对应 LungAgingERV 项目的三个数据来源。
#' @param new_labels Character vector. 新的标准化标签，按位置与 `old_levels` 一一对应。
#'   长度必须与 `old_levels` 一致。
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
#' @examples
#' \dontrun{
#' # 使用默认映射
#' immune.combined <- rename_dataset_origin(immune.combined)
#'
#' # 自定义映射（用于其他项目）
#' immune.combined <- rename_dataset_origin(
#'   immune.combined,
#'   old_levels = c("A", "B"),
#'   new_labels = c("X_2020", "Y_2021")
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
  missing_levels <- setdiff(old_levels, current_levels)
  if (length(missing_levels) > 0) {
    warning("以下 `old_levels` 未在 dataset_origin 中找到: ",
            paste(missing_levels, collapse = ", "))
  }

  # --- 自动备份 ---
  if (backup) {
    if (verbose) message("正在保存备份至: ", backup_path)
    saveRDS(seurat_obj, backup_path)
  }

  # --- 重命名 ---
  seurat_obj@meta.data$dataset_origin <- factor(
    seurat_obj@meta.data$dataset_origin,
    levels = old_levels,
    labels = new_labels
  )

  # --- 输出摘要 ---
  if (verbose) {
    message("dataset_origin 重命名完成：")
    print(table(seurat_obj@meta.data$dataset_origin))
  }

  return(seurat_obj)
}
```

---

## 📦 整合进 jwtools 包的具体步骤

### 1. 创建脚本文件（放入 `R/` 目录）

```r
usethis::use_r("rename_dataset_origin")
```
这会自动在 `jwtools/R/rename_dataset_origin.R` 创建文件，把上面代码粘贴进去即可。

### 2. 声明依赖（Seurat）

由于函数内部用了 `inherits(seurat_obj, "Seurat")`，需要在 `DESCRIPTION` 中声明 Seurat 为依赖（Imports 而非 Depends，避免污染用户命名空间）：

```r
usethis::use_package("Seurat", type = "Imports")
```

### 3. 生成文档 + NAMESPACE

```r
devtools::document()
```
这一步会根据 roxygen2 注释自动生成 `man/rename_dataset_origin.Rd` 帮助文档，并在 `NAMESPACE` 中添加 `export(rename_dataset_origin)`。

### 4. 本地测试

```r
devtools::load_all()
?rename_dataset_origin   # 查看帮助文档
```

### 5. 检查包完整性（可选但推荐）

```r
devtools::check()
```

---

## 💡 使用效果预览

集成完成后，您在任何项目中只需：

```r
library(jwtools)

immune.combined <- rename_dataset_origin(immune.combined)
#> 正在保存备份至: backup_before_rename_dataset_origin.rds
#> dataset_origin 重命名完成：
#> 
#>     Adams_2020 Habermann_2020     Natri_2024 
#>          XXXXX          XXXXX          XXXXX
```

一行代码即可完成"备份 + 重命名 + 验证"三步操作，且这个函数具有通用性——未来若您的其他项目（如 HeartERVmap、PH-ERVAtlas）也需要类似的数据集来源重命名逻辑，只需修改 `old_levels`/`new_labels` 参数即可复用，无需重写代码。

