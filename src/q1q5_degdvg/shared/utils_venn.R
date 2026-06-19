# utils_venn.R
# Generalised wrapper around VennDiagram::venn.diagram.
#
# Replaces the ad-hoc call pattern used in:
#   src/scripts/15analysis/figures/plot_venn_mad_shift_normal_vs_permute.R
#
# Usage:
#   source("shared/utils_venn.R")
#   draw_named_venn(
#     sets_list   = list(limma = c("g1","g2"), DESeq2 = c("g2","g3")),
#     output_path = "venn_updeg.png",
#     title       = "Up-DEG overlap"
#   )

draw_named_venn <- function(sets_list,
                            output_path,
                            title  = "",
                            colors = NULL,
                            ...) {
  n <- length(sets_list)
  if (n < 2L || n > 4L)
    stop("draw_named_venn: sets_list must have 2-4 named elements (got ", n, ")")
  if (is.null(names(sets_list)) || any(nchar(names(sets_list)) == 0L))
    stop("draw_named_venn: all elements of sets_list must be named")

  if (is.null(colors))
    colors <- RColorBrewer::brewer.pal(max(3L, n), "Pastel2")[seq_len(n)]

  log_path <- VennDiagram::venn.diagram(
    x        = sets_list,
    filename = output_path,
    fill     = colors,
    main     = title,
    output   = TRUE,
    ...
  )
  invisible(log_path)
}
