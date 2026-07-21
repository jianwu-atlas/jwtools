#' 保存 ggplot 对象为 PDF + PNG 双格式
#'
#' @param plot_obj ggplot 对象
#' @param filename_base 不含后缀的输出路径前缀
#' @param width 图宽 (inch)
#' @param height 图高 (inch)
#' @param dpi PNG分辨率
#'
#' @return invisible(TRUE/FALSE)，标记保存是否成功
#' @export
save_dual_gg <- function(plot_obj, filename_base, width = 8, height = 6, dpi = 300) {
  dir.create(dirname(filename_base), showWarnings = FALSE, recursive = TRUE)
  ok <- tryCatch({
    ggsave(paste0(filename_base, ".pdf"), plot_obj, width = width, height = height,
           device = "pdf", useDingbats = FALSE)
    ggsave(paste0(filename_base, ".png"), plot_obj, width = width, height = height,
           dpi = dpi, device = "png", bg = "white")
    TRUE
  }, error = function(e) {
    message("  ❗ 保存失败: ", filename_base, " | ", conditionMessage(e))
    FALSE
  })
  if (ok) message("  ✅ 已保存: ", filename_base, ".pdf / .png")
  invisible(ok)
}
