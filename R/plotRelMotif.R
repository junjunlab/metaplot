#' motif metaplot
#'
#' @param relpos a list of relative position from CalculatePeaksRelativeDist function.
#' A single numeric vector is also accepted and treated as one sample.
#' @param exp sample name. Default is `NULL`, which uses the names of `relpos`
#' or `sample1`, `sample2`, ... when the list is unnamed.
#' @param linewidth line width.
#'
#' @return ggplot
#' @export
plotRelMotif <- function(relpos = NULL,
                         exp = NULL,
                         linewidth = 0.5){

  # ==============================================================================
  # validate inputs
  # ==============================================================================
  if(is.null(relpos)){
    stop("`relpos` must be supplied, e.g. the output of ",
         "`CalculatePeaksRelativeDist()`.",call. = FALSE)
  }

  if(is.numeric(relpos)){
    relpos <- list(relpos)
  }

  if(!is.list(relpos) || length(relpos) == 0){
    stop("`relpos` must be a non-empty list of numeric vectors.",call. = FALSE)
  }

  is_numeric_element <- vapply(relpos,is.numeric,logical(1))
  if(!all(is_numeric_element)){
    stop("Every element of `relpos` must be numeric; element(s) ",
         paste(which(!is_numeric_element),collapse = ", "),"are not.",
         call. = FALSE)
  }

  if(is.null(exp)){
    exp <- names(relpos)
  }

  if(is.null(exp)){
    exp <- paste0("sample",seq_along(relpos))
  }

  if(length(exp) != length(relpos)){
    stop("`exp` has length ",length(exp)," but `relpos` has ",
         length(relpos)," element(s); they must match.",call. = FALSE)
  }

  lapply(seq_along(relpos), function(x){
    data.frame(rel_pos = relpos[[x]],exp = exp[x])
  }) %>% do.call("rbind",.) -> motif_df

  if(nrow(motif_df) == 0){
    stop("No relative position was supplied, nothing to plot.",call. = FALSE)
  }

  # ============================================================================
  # plot

  p <- ggplot(motif_df) +
    geom_density(aes(x = rel_pos,color = exp),linewidth = linewidth) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    xlab("m6A(RRACH) motif relative to peak center") +
    ylab("peak density")

  # Set1 only provides 9 colours, fall back to the default hue palette beyond that
  if(length(unique(motif_df$exp)) <= 9){
    p <- p + scale_color_brewer(palette = "Set1")
  }

  p
}
