# ============================================================
# File: R/save_dual.R
# Package: jwtools
# Purpose: 双格式(PDF+PNG)保存 ggplot 对象的通用工具函数
# Author: Jian Wu
# ============================================================

#' 双格式保存 ggplot 对象（PDF + PNG）
#'
#' 同时保存 ggplot / patchwork 对象为矢量图(PDF)和位图(PNG)两种格式,
#' 适用于分析流程中需要"出版级矢量图 + 快速预览位图"并存的场景。
#'
#' @param filename_prefix 不含扩展名的文件名前缀
#' @param plot ggplot / patchwork 对象
#' @param width,height 图形尺寸（inch）
#' @param dpi PNG分辨率，默认300
#' @param bg 背景色，默认white
#'
#' @return 无返回值, 调用后在磁盘上生成对应的 .pdf 和 .png 文件(invisible NULL)
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(mpg, wt)) + geom_point()
#' save_dual("results/mtcars_scatter", p, width = 6, height = 5)
#' }
#'
#' @importFrom ggplot2 ggsave
#' @export
save_dual <- function(filename_prefix, plot, width, height,
                      dpi = 300, bg = "white") {
  ggplot2::ggsave(paste0(filename_prefix, ".pdf"), plot, width = width, height = height)
  ggplot2::ggsave(paste0(filename_prefix, ".png"), plot, width = width, height = height,
                  dpi = dpi, bg = bg)
  message("  \u2705 \u5df2\u4fdd\u5b58: ", filename_prefix, ".pdf / .png")
}
