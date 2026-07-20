# R/globals.R
# ==============================================================================
# 目的: 消除 R CMD check 中因 dplyr/ggplot2 tidy-eval 语法产生的
#       "no visible binding for global variable" NOTE。
#       这些变量名均为数据框列名(在 mutate()/aes() 等函数体内部使用),
#       并非真正的全局变量, 静态代码分析无法识别, 此处显式声明以消除误报。
# ==============================================================================
utils::globalVariables(c(
  "bracket_rank", "celltype", "dataset", "group", "group1", "group2",
  "iqr", "is_outlier", "n_cells", "p", "p.adj", "p.adj.signif",
  "percent", "q1", "q3", "span", "total_cells", "xmax", "xmin", "y_max"
))
