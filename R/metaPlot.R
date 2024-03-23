globalVariables(c("end", "end.1", "f_len", "gene_id", "gene_name",
                  "p_mid", "rel_pos", "seqnames", "start", "start.1", ".",
                  "strand", "transcript_id", "type", "utr3.sf", "utr5.sf", "width"))

#' metaplot
#'
#' @param bed_file bed format of peaks files.
#' @param group_names vector of group names.
#' @param features_anno features_anno object from prepare_features functions.
#' @param scale_region whether do scale for each feature length, default FALSE.
#' @param cut_ratio whether use custom ratio to plot feature regions,
#' for eaxample: c(0.1,0.8,0.1).
#' @param facet whether dispaly plot for each group, default FALSE.
#'
#' @return A list with data and plot
#'
#' @import ggplot2
#'
#' @export
metaPlot <- function(bed_file = NULL,
                     group_names = NULL,
                     features_anno = NULL,
                     scale_region = FALSE,
                     cut_ratio = NULL,
                     facet = FALSE){
  # ==============================================================================
  # overlap peaks and anno
  # ==============================================================================
  if(is.null(group_names)){
    group_names <- bed_file
  }

  lapply(seq_along(bed_file),function(x){
    df <- calculate_relative_position(bed_file = bed_file[[x]],
                                      features_anno = features_anno,
                                      scale_region = scale_region,
                                      cut_ratio = cut_ratio)

    # add group name
    df$group_names <- group_names[x]

    return(df)
  }) %>% do.call("rbind",.) -> plot_df

  # ==============================================================================
  # plot
  # ==============================================================================
  utr5.sf <- attr(plot_df,"utr5.sf")
  utr3.sf <- attr(plot_df,"utr3.sf")

  if(scale_region == TRUE){
    label.x <- c((1 - utr5.sf + 1)/2,1.5,(2 + utr3.sf + 2)/2)
  }else{
    label.x <- c(0.5,1.5,2.5)
  }

  # check samples number
  if(length(bed_file) > 1){
    density_layer <- geom_density(aes(x = rel_pos,fill = group_names),color = NA,alpha = 0.5)
  }else{
    density_layer <- geom_density(aes(x = rel_pos),fill = "#CC0033",color = NA)
  }

  pdensity <-
    ggplot(plot_df) +
    density_layer +
    geom_vline(xintercept = c(1,2),lty = "dashed") +
    theme_classic(base_size = 12) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    scale_x_continuous(breaks = label.x,
                       labels = c("5'UTR","CDS","3'UTR")) +
    theme(axis.line = element_line(arrow = arrow(length = unit(0.25,"cm"),
                                                 type = "closed")),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(face = "bold")) +
    ylab("Peaks density") + xlab("")

  if(facet == TRUE){
    plot <- pdensity +
      facet_wrap(~group_names,scales = "free")
  }else{
    plot <- pdensity
  }

  # output
  res <- list(data = plot_df,
              plot = plot)

  return(res)
}
