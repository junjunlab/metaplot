#' motif metaplot
#'
#' @param relpos a list of relative position from CalculatePeaksRelativeDist function.
#' @param exp sample name.
#' @param linewidth line width.
#'
#' @return ggplot
#' @export
plotRelMotif <- function(relpos = NULL,
                         exp = NULL,
                         linewidth = 0.5){

  # x = 1
  lapply(seq_along(relpos), function(x){
    data.frame(rel_pos = relpos[[x]],exp = exp[x])
  }) %>% do.call("rbind",.) -> motif_df

  # ============================================================================
  # plot

  ggplot(motif_df) +
    geom_density(aes(x = rel_pos,color = exp),linewidth = 0.5) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_color_brewer(palette = "Set1") +
    xlab("m6A(RRACH) motif relative to peak center") +
    ylab("peak density")
}
