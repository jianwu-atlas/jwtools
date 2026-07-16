# 整合为 jwtools 函数：`ct_proportion_analysis()`

理解你的需求。我把整个 Extended Data Fig1F 脚本重构成一个**通用化、参数化的函数**，可以直接放进你的 `jwtools` 包里。设计原则：

1. **泛化所有硬编码**：celltype顺序、颜色、分组变量名、比较列表都改为参数，日后可用于任何项目（LungERVmap / HeartERVmap / PH-ERVAtlas 等）而不需要复制粘贴改代码
2. **dataset overlay 变为可选**（`dataset_col = NULL` 时自动跳过，退化为普通 jitter）
3. **保留完整的统计输出**（Kruskal-Wallis + pairwise Wilcoxon + outlier检测）
4. **`save_dual()` 作为包内辅助函数一并提供**，避免每个脚本都要 `if(!exists(...))` 判断

---

## 第一步：新建 `R/ct_proportion_analysis.R`

```r
# ============================================================
# File: R/ct_proportion_analysis.R
# Package: jwtools
# Purpose: 通用化的"细胞类型比例 × 分组变量"统计可视化函数
#          可复用于 LungAgingERV / LungERVIPF / PH-ERVAtlas 等
#          任何需要"donor-level proportion + pairwise stats +
#          batch/dataset overlay"的项目
# Author: Jian Wu
# ============================================================

#' 双格式保存 ggplot 对象（PDF + PNG）
#'
#' @param filename_prefix 不含扩展名的文件名前缀
#' @param plot ggplot / patchwork 对象
#' @param width,height 图形尺寸（inch）
#' @param dpi PNG分辨率，默认300
#' @param bg 背景色，默认white
#' @export
save_dual <- function(filename_prefix, plot, width, height,
                      dpi = 300, bg = "white") {
  ggplot2::ggsave(paste0(filename_prefix, ".pdf"), plot, width = width, height = height)
  ggplot2::ggsave(paste0(filename_prefix, ".png"), plot, width = width, height = height,
                  dpi = dpi, bg = bg)
  message("  \u2705 \u5df2\u4fdd\u5b58: ", filename_prefix, ".pdf / .png")
}


#' 细胞类型比例的分组统计分析与可视化（支持dataset/batch overlay）
#'
#' 计算 donor-level 的细胞类型比例，在指定分组变量（如 AgeGroup_short）
#' 间进行 Kruskal-Wallis 总检验 + 两两 Wilcoxon 检验（BH校正），
#' 并生成 boxplot + jitter 图；可选地将样本点按 dataset/batch 来源
#' 标注 shape+color，用于核查"组间差异是否被单一队列驱动"。
#'
#' @param seurat_obj Seurat对象，或直接传入包含相应列的 data.frame
#' @param celltype_col 细胞类型列名，默认 "Manuscript_Identity"
#' @param sample_col   样本(donor)ID列名，默认 "sample"
#' @param group_col    分组变量列名（如年龄组），默认 "AgeGroup_short"
#' @param group_levels 分组变量的顺序，NULL则自动从factor levels/排序推断
#' @param group_labels 分组变量的图例显示文字（named vector），NULL则用原始值
#' @param group_colors 分组变量的boxplot填充色（named vector），NULL则自动配色
#' @param dataset_col  数据集/批次来源列名，NULL则不绘制dataset overlay
#' @param dataset_colors dataset的散点颜色（named vector），NULL则自动生成
#' @param dataset_shapes dataset的散点形状（named vector），NULL则自动生成
#' @param celltype_order celltype在facet中的排列顺序，NULL则自动推断
#' @param celltype_colors celltype对应的strip文字颜色（named vector），NULL则灰色
#' @param comparisons 两两比较列表，如 list(c("Y","M"), c("M","O"), c("Y","O"))；
#'                    NULL则自动生成 group_levels 的全部两两组合
#' @param p_adjust_method p值校正方法，默认"BH"
#' @param facet_nrow  facet的行数，默认2
#' @param output_prefix 输出文件名前缀，默认"CTprop"
#' @param save_plots  是否自动保存PDF+PNG，默认TRUE
#' @param fig_width,fig_height 保存图形的尺寸
#' @param label_outliers 是否额外生成"标注离群样本ID"的诊断版本图，默认FALSE
#' @param outlier_coef IQR倍数阈值，默认1.5
#' @param verbose 是否打印中间信息，默认TRUE
#'
#' @return list，包含：
#'   \item{prop_df}{donor-level细胞比例明细表}
#'   \item{kw_test}{每个celltype的Kruskal-Wallis总检验结果}
#'   \item{pairwise_test}{每个celltype的两两Wilcoxon检验结果（BH校正）}
#'   \item{stat_test}{用于画bracket的完整统计表（含坐标）}
#'   \item{plot}{主图 ggplot对象}
#'   \item{plot_labeled}{（可选）含离群样本标注的诊断图}
#'   \item{outlier_df}{离群样本明细表}
#'
#' @examples
#' \dontrun{
#' res <- ct_proportion_analysis(
#'   seurat_obj    = immune.combined,
#'   celltype_col  = "Manuscript_Identity",
#'   group_col     = "AgeGroup_short",
#'   group_levels  = c("Y", "M", "O"),
#'   dataset_col   = "dataset_origin",
#'   celltype_order = c("EC","Epi","B","DC","Mac","Mono","NK","T","Fib","SMC"),
#'   celltype_colors = cols,   # 复用你主UMAP的配色向量
#'   output_prefix = "28_ExtFig1F"
#' )
#' res$plot
#' }
#' @export
ct_proportion_analysis <- function(
    seurat_obj,
    celltype_col     = "Manuscript_Identity",
    sample_col       = "sample",
    group_col        = "AgeGroup_short",
    group_levels     = NULL,
    group_labels     = NULL,
    group_colors     = NULL,
    dataset_col      = NULL,
    dataset_colors   = NULL,
    dataset_shapes   = NULL,
    celltype_order   = NULL,
    celltype_colors  = NULL,
    comparisons      = NULL,
    p_adjust_method  = "BH",
    facet_nrow       = 2,
    output_prefix    = "CTprop",
    save_plots       = TRUE,
    fig_width        = 15,
    fig_height       = 7.5,
    label_outliers   = FALSE,
    outlier_coef     = 1.5,
    verbose          = TRUE
) {

  # ---- 依赖包检查 ----
  req_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "ggpubr",
               "rstatix", "ggh4x", "ggrepel", "scales")
  missing_pkgs <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("\u7f3a\u5c11\u4f9d\u8d56\u5305\uff0c\u8bf7\u5148\u5b89\u88c5: ",
         paste(missing_pkgs, collapse = ", "))
  }

  `%>%` <- dplyr::`%>%`

  # ============================================================
  # Step 0. 提取 meta.data（兼容 Seurat 对象 或 已有 data.frame）
  # ============================================================
  if (inherits(seurat_obj, "Seurat")) {
    meta_full <- seurat_obj@meta.data
  } else if (is.data.frame(seurat_obj)) {
    meta_full <- seurat_obj
  } else {
    stop("seurat_obj \u5fc5\u987b\u662f Seurat \u5bf9\u8c61\u6216 data.frame")
  }

  needed_cols <- c(sample_col, group_col, celltype_col)
  if (!is.null(dataset_col)) needed_cols <- c(needed_cols, dataset_col)
  missing_cols <- setdiff(needed_cols, colnames(meta_full))
  if (length(missing_cols) > 0) {
    stop("meta.data \u4e2d\u7f3a\u5c11\u4ee5\u4e0b\u5217: ", paste(missing_cols, collapse = ", "))
  }

  meta <- meta_full %>%
    tibble::as_tibble(rownames = "cellid") %>%
    dplyr::select(dplyr::all_of(c(sample_col, group_col, celltype_col,
                                  if (!is.null(dataset_col)) dataset_col)))

  # 统一重命名为内部标准列名，方便后续处理
  rename_map <- c(sample = sample_col, group = group_col, celltype = celltype_col)
  if (!is.null(dataset_col)) rename_map <- c(rename_map, dataset = dataset_col)
  meta <- meta %>% dplyr::rename(!!!rename_map)

  # ============================================================
  # Step 1. 自动推断 group_levels / celltype_order / 默认配色
  # ============================================================
  if (is.null(group_levels)) {
    group_levels <- if (is.factor(meta$group)) levels(meta$group) else sort(unique(meta$group))
  }
  if (is.null(group_labels)) {
    group_labels <- stats::setNames(group_levels, group_levels)
  }
  if (is.null(group_colors)) {
    default_pal  <- scales::hue_pal()(length(group_levels))
    group_colors <- stats::setNames(default_pal, group_levels)
  }

  if (is.null(celltype_order)) {
    celltype_order <- if (is.factor(meta$celltype)) levels(meta$celltype) else sort(unique(meta$celltype))
  }
  stopifnot(all(celltype_order %in% unique(meta$celltype)))

  if (!is.null(dataset_col)) {
    ds_levels <- if (is.factor(meta$dataset)) levels(meta$dataset) else sort(unique(meta$dataset))
    if (is.null(dataset_colors)) {
      ds_pal <- RColorBrewer::brewer.pal(max(3, length(ds_levels)), "Dark2")[seq_along(ds_levels)]
      dataset_colors <- stats::setNames(ds_pal, ds_levels)
    }
    if (is.null(dataset_shapes)) {
      shape_pool <- c(16, 17, 15, 18, 3, 4, 8)
      dataset_shapes <- stats::setNames(shape_pool[seq_along(ds_levels)], ds_levels)
    }
  }

  if (is.null(celltype_colors)) {
    celltype_colors <- stats::setNames(rep("black", length(celltype_order)), celltype_order)
  }

  if (verbose) {
    message("group_levels: ", paste(group_levels, collapse = ", "))
    message("celltype_order: ", paste(celltype_order, collapse = ", "))
    if (!is.null(dataset_col)) message("dataset levels: ", paste(ds_levels, collapse = ", "))
  }

  # ============================================================
  # Step 2. donor-level 细胞比例计算（补零，构建完整网格）
  # ============================================================
  group_by_cols_total  <- c("sample", "group", if (!is.null(dataset_col)) "dataset")
  group_by_cols_counts <- c(group_by_cols_total, "celltype")

  sample_total <- meta %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_by_cols_total))) %>%
    dplyr::summarise(total_cells = dplyr::n(), .groups = "drop")

  sample_counts <- meta %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_by_cols_counts))) %>%
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop")

  full_grid <- sample_total %>%
    tidyr::crossing(celltype = celltype_order)

  join_cols <- group_by_cols_counts

  prop_df <- full_grid %>%
    dplyr::left_join(sample_counts, by = join_cols) %>%
    dplyr::mutate(
      n_cells = tidyr::replace_na(n_cells, 0),
      percent = n_cells / total_cells * 100
    )

  prop_df$celltype <- factor(prop_df$celltype, levels = celltype_order)
  prop_df$group    <- factor(prop_df$group,    levels = group_levels)
  if (!is.null(dataset_col)) prop_df$dataset <- factor(prop_df$dataset, levels = ds_levels)

  if (verbose) {
    message("\n\u6bcf\u4e2acelltype\u7684sample\u6570 (\u6309group):")
    print(prop_df %>% dplyr::count(celltype, group) %>%
           tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0))
  }

  if (save_plots) {
    utils::write.csv(prop_df, paste0(output_prefix, "_proportion_perSample.csv"), row.names = FALSE)
  }

  # ============================================================
  # Step 3. 统计检验：Kruskal-Wallis(整体) + 两两Wilcoxon(BH校正)
  # ============================================================
  if (is.null(comparisons)) {
    combo_mat   <- utils::combn(group_levels, 2)
    comparisons <- lapply(seq_len(ncol(combo_mat)), function(i) combo_mat[, i])
  }

  kw_test <- prop_df %>%
    dplyr::group_by(celltype) %>%
    rstatix::kruskal_test(percent ~ group) %>%
    dplyr::ungroup()

  pairwise_test <- prop_df %>%
    dplyr::group_by(celltype) %>%
    rstatix::wilcox_test(
      percent ~ group,
      comparisons     = comparisons,
      p.adjust.method = p_adjust_method
    ) %>%
    dplyr::ungroup() %>%
    rstatix::add_significance("p.adj")

  if (save_plots) {
    utils::write.csv(kw_test, paste0(output_prefix, "_KruskalWallis_overall.csv"), row.names = FALSE)
    utils::write.csv(pairwise_test, paste0(output_prefix, "_PairwiseWilcoxon_stats.csv"), row.names = FALSE)
  }

  # ---- 计算bracket坐标：按跨度(span)升序叠加y.position，避免重叠 ----
  y_max_df   <- prop_df %>%
    dplyr::group_by(celltype) %>%
    dplyr::summarise(y_max = max(percent, na.rm = TRUE), .groups = "drop")

  x_pos_map <- stats::setNames(seq_along(group_levels), group_levels)

  stat_test <- pairwise_test %>%
    dplyr::left_join(y_max_df, by = "celltype") %>%
    dplyr::mutate(
      xmin = x_pos_map[group1],
      xmax = x_pos_map[group2],
      span = abs(xmax - xmin)
    ) %>%
    dplyr::group_by(celltype) %>%
    dplyr::arrange(span, .by_group = TRUE) %>%
    dplyr::mutate(
      bracket_rank = dplyr::row_number(),
      y.position   = y_max * (1 + 0.14 * bracket_rank)
    ) %>%
    dplyr::ungroup()

  if (verbose) {
    message("\n=== Pairwise Wilcoxon (BH-adjusted) ===")
    print(stat_test %>% dplyr::select(celltype, group1, group2, p, p.adj, p.adj.signif))
  }

  # ============================================================
  # Step 4. 绘图
  # ============================================================
  strip_text_colors <- celltype_colors[celltype_order]

  p_main <- ggplot2::ggplot(prop_df, ggplot2::aes(x = group, y = percent)) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = group),
      outlier.shape = NA, width = 0.6, linewidth = 0.4,
      color = "grey20", alpha = 0.85
    )

  if (!is.null(dataset_col)) {
    p_main <- p_main +
      ggplot2::geom_jitter(
        ggplot2::aes(shape = dataset, color = dataset),
        width = 0.15, size = 1.6, alpha = 0.8, stroke = 0.6
      ) +
      ggplot2::scale_color_manual(values = dataset_colors, name = "Dataset of origin") +
      ggplot2::scale_shape_manual(values = dataset_shapes, name = "Dataset of origin")
  } else {
    p_main <- p_main +
      ggplot2::geom_jitter(width = 0.15, size = 1.4, alpha = 0.6, color = "black")
  }

  p_main <- p_main +
    ggpubr::stat_pvalue_manual(
      stat_test, label = "p.adj.signif",
      xmin = "xmin", xmax = "xmax", y.position = "y.position",
      tip.length = 0.01, size = 3, bracket.size = 0.35,
      hide.ns = FALSE
    ) +
    ggplot2::scale_fill_manual(values = group_colors, labels = group_labels, name = "Group") +
    ggplot2::scale_x_discrete(labels = group_labels) +
    ggh4x::facet_wrap2(
      ~ celltype, nrow = facet_nrow, scales = "free_y",
      strip = ggh4x::strip_themed(
        background_x = lapply(celltype_order, function(ct) ggplot2::element_rect(fill = "grey85", color = NA)),
        text_x = lapply(celltype_order, function(ct)
          ggplot2::element_text(color = strip_text_colors[ct], size = 10, face = "bold"))
      )
    ) +
    ggplot2::labs(
      x = NULL, y = "Percentage from total cells",
      caption = if (!is.null(dataset_col))
        "Each point represents one donor sample; point shape/color denotes the source cohort."
      else "Each point represents one donor sample."
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      legend.position    = "bottom",
      legend.box         = "vertical",
      legend.title       = ggplot2::element_text(size = 9, face = "bold"),
      legend.text        = ggplot2::element_text(size = 8),
      legend.key.size    = grid::unit(0.4, "cm"),
      axis.text.x        = ggplot2::element_text(color = "black", size = 9),
      axis.text.y        = ggplot2::element_text(color = "black", size = 9),
      axis.title.y       = ggplot2::element_text(size = 10),
      panel.grid.major.y = ggplot2::element_line(color = "grey92", linewidth = 0.3),
      strip.clip         = "off",
      panel.spacing      = grid::unit(0.5, "cm"),
      plot.caption       = ggplot2::element_text(size = 8, color = "grey30", hjust = 0),
      plot.margin        = ggplot2::margin(8, 10, 8, 8)
    ) +
    ggplot2::guides(
      fill  = ggplot2::guide_legend(order = 1, override.aes = list(shape = NA))
    )

  if (save_plots) {
    save_dual(paste0(output_prefix, "_main"), p_main, width = fig_width, height = fig_height)
  }

  # ============================================================
  # Step 5.（可选）离群样本标注诊断图
  # ============================================================
  outlier_df   <- NULL
  p_labeled    <- NULL

  if (label_outliers) {
    outlier_df <- prop_df %>%
      dplyr::group_by(celltype, group) %>%
      dplyr::mutate(
        q1 = stats::quantile(percent, 0.25, na.rm = TRUE),
        q3 = stats::quantile(percent, 0.75, na.rm = TRUE),
        iqr = q3 - q1,
        is_outlier = percent > q3 + outlier_coef * iqr | percent < q1 - outlier_coef * iqr
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(is_outlier)

    if (verbose) {
      message("\u68c0\u6d4b\u5230 ", nrow(outlier_df), " \u4e2a\u79bb\u7fa4\u6837\u672c\u70b9")
    }

    if (save_plots) {
      utils::write.csv(outlier_df, paste0(output_prefix, "_outlier_samples.csv"), row.names = FALSE)
    }

    p_labeled <- p_main +
      ggrepel::geom_text_repel(
        data = outlier_df,
        ggplot2::aes(x = group, y = percent, label = sample),
        size = 2.2, color = "black", max.overlaps = 15,
        segment.size = 0.2, box.padding = 0.3
      )

    if (save_plots) {
      save_dual(paste0(output_prefix, "_labeled"), p_labeled, width = fig_width, height = fig_height)
    }
  }

  # ============================================================
  # 返回结果
  # ============================================================
  invisible(list(
    prop_df       = prop_df,
    kw_test       = kw_test,
    pairwise_test = pairwise_test,
    stat_test     = stat_test,
    plot          = p_main,
    plot_labeled  = p_labeled,
    outlier_df    = outlier_df
  ))
}
```

---

## 第二步：注册到包并安装

```r
# 在 jwtools 项目根目录下执行：
usethis::use_package("dplyr")
usethis::use_package("tidyr")
usethis::use_package("tibble")
usethis::use_package("ggplot2")
usethis::use_package("ggpubr")
usethis::use_package("rstatix")
usethis::use_package("ggh4x")
usethis::use_package("ggrepel")
usethis::use_package("scales")
usethis::use_package("RColorBrewer")

devtools::document()   # 生成 man/*.Rd 和 NAMESPACE
devtools::install()    # 重新安装到本地，覆盖旧版本
```

---

## 第三步：调用方式（原本300行 → 现在10行）

```r
library(jwtools)

res <- ct_proportion_analysis(
  seurat_obj      = immune.combined,
  celltype_col    = "Manuscript_Identity",
  sample_col      = "sample",
  group_col       = "AgeGroup_short",
  group_levels    = c("Y", "M", "O"),
  group_labels    = c("Y" = "Y (\u226435 y.o)", "M" = "M (35-50 y.o)", "O" = "O (>50 y.o)"),
  group_colors    = c("Y" = "#F4A9A8", "M" = "#8FCB8F", "O" = "#7FB3D5"),
  dataset_col     = "dataset_origin",          # NULL的话则不画dataset overlay
  celltype_order  = c("EC","Epi","B","DC","Mac","Mono","NK","T","Fib","SMC"),
  celltype_colors = cols,                       # 复用你主UMAP的10色配色向量
  output_prefix   = "28_ExtFig1F",
  label_outliers  = TRUE                        # 需要溯源诊断图时打开
)

res$plot            # 主图
res$plot_labeled     # 含离群样本ID的诊断图
res$pairwise_test    # 两两Wilcoxon检验结果表
```

---

## 复用场景举例（体现"泛化"的价值）

同一个函数，改几个参数即可直接用于其他项目：

```r
# 例：PH-ERVAtlas 项目，分组变量换成 "Disease" (Control vs PAH)
ct_proportion_analysis(
  seurat_obj     = ph_combined,
  group_col      = "Disease",
  group_levels   = c("Control", "PAH"),
  dataset_col    = NULL,            # 单一数据集，不需要overlay
  celltype_order = ph_celltype_order,
  output_prefix  = "PH_CTprop"
)

# 例：LungERVIPF 项目，比较 Control vs IPF，同时保留数据集溯源
ct_proportion_analysis(
  seurat_obj     = ipf_combined,
  group_col      = "Diagnosis",
  group_levels   = c("Control", "IPF"),
  dataset_col    = "cohort",
  output_prefix  = "IPF_CTprop"
)
```

---

## ⚠️ 需要你确认的两点

1. **`sample_col` 的处理**：如果同一个 sample 在你的 meta.data 里确实只对应唯一一个 `dataset`（即不存在一个sample跨两个dataset的情况），当前的 `group_by(sample, group, dataset)` 逻辑是安全的。但如果存在数据不一致（比如同一sample ID在不同队列间重复使用导致冲突），会导致 `crossing()` 后网格重复行。建议先跑一次这行做健全性检查：

```r
immune.combined@meta.data %>%
  dplyr::distinct(sample, dataset_origin) %>%
  dplyr::count(sample) %>%
  dplyr::filter(n > 1)   # 若非空，说明有sample对应多个dataset，需要处理
```

2. **`ggh4x::strip_themed()` 中 `celltype_colors` 未提供时全部显示黑色**——如果你希望默认也能有区分度（比如复用 `scales::hue_pal()`自动生成），我可以再加一个 fallback 逻辑，告诉我是否需要。
