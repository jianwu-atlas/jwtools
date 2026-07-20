#' 双格式保存 base graphics 绘图 (PDF + PNG)
#'
#' 用于保存 pheatmap / ComplexHeatmap::draw() / plot() / scatterplot3d() 等
#' "命令式"绘图函数产生的图形。与 \code{\link{save_dual}} (专用于ggplot/patchwork对象)
#' 互补, 二者配合可覆盖分析流程中几乎所有绘图场景。
#'
#' @param draw_fun 一个无参函数, 调用时会把绘图命令执行到当前打开的图形设备上,
#'   例如 \code{function() pheatmap::pheatmap(mat)} 或
#'   \code{function() ComplexHeatmap::draw(ht_obj)}。
#' @param filename_prefix 不含扩展名的文件名前缀
#' @param width,height 图形尺寸(inch)
#' @param res PNG分辨率(dpi), 默认300
#'
#' @return 逻辑值, 表示PDF是否绘制成功(隐式返回, invisible)
#'
#' @examples
#' \dontrun{
#' save_dual_base(
#'   draw_fun = function() pheatmap::pheatmap(mat, cluster_rows = TRUE),
#'   filename_prefix = "results/heatmap_TE",
#'   width = 8, height = 9
#' )
#' }
#' @importFrom grDevices pdf png dev.off
#' @importFrom graphics plot.new text
#' @export
save_dual_base <- function(draw_fun, filename_prefix, width = 8, height = 6, res = 300) {
  dir.create(dirname(filename_prefix), showWarnings = FALSE, recursive = TRUE)

  pdf(paste0(filename_prefix, ".pdf"), width = width, height = height)
  ok_pdf <- tryCatch({
    draw_fun()
    TRUE
  }, error = function(e) {
    message("  \u274c \u7ed8\u56fe\u5931\u8d25(PDF): ", conditionMessage(e))
    FALSE
  }, finally = { dev.off() })

  if (ok_pdf) {
    png(paste0(filename_prefix, ".png"), width = width, height = height,
        units = "in", res = res, bg = "white")
    tryCatch({
      draw_fun()
    }, error = function(e) {
      message("  \u274c \u7ed8\u56fe\u5931\u8d25(PNG): ", conditionMessage(e))
    }, finally = { dev.off() })
    message("  \u2705 \u5df2\u4fdd\u5b58: ", filename_prefix, ".pdf / .png")
  } else {
    png(paste0(filename_prefix, ".png"), width = width, height = height,
        units = "in", res = res, bg = "white")
    plot.new()
    text(0.5, 0.5, "No plottable data", cex = 1.2)
    dev.off()
    message("  \u26a0\ufe0f  \u4ec5\u751f\u6210\u5360\u4f4d\u56fePNG: ", filename_prefix, ".png")
  }

  invisible(ok_pdf)
}
