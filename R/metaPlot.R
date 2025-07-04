globalVariables(c("end", "end.1", "f_len", "gene_id", "gene_name",
                  "p_mid", "rel_pos", "seqnames", "start", "start.1", ".",
                  "strand", "transcript_id", "type", "utr3.sf", "utr5.sf", "width"))

#' MetaPlot Function
#'
#' Generate a meta plot of peak densities across genomic features such as UTRs and CDS using BED files.
#'
#' @param bed_file A list of BED file paths or data frames containing regions of interest. Default is `NULL`.
#' @param group_names A vector of names corresponding to each BED file for grouping in the plot. Default is `NULL`, automatically set to `bed_file`.
#' @param features_anno A features_anno object or file path used to map features. Default is `NULL`.
#' @param scale_region Logical, whether to scale regions relative to feature lengths. Default is `FALSE`.
#' @param cut_ratio A numeric value indicating the ratio to cut the regions for scaling. Default is `NULL`.
#' @param facet Logical, whether to use faceting in the plot for separate groups. Default is `FALSE`.
#' @param geom_density_bw Numeric, specifying the bandwidth for `geom_density`. Default is `NULL`, uses automatic bandwidth selection.
#' @param geom_density_n Integer, number of equally spaced points at which the density is estimated. Default is `512`.
#' @param col_alpha Numeric, alpha value for color transparency in the density plot. Default is `0.5`.
#' @param line_color A color code for lines in the density plot if there is only one bed file. Default is `"#CC0033"`.
#'
#' @return A list with two components:
#' \describe{
#'   \item{data}{A combined data frame with relative positions of peaks against annotated features.}
#'   \item{plot}{A ggplot2 object representing the density plot of peaks across the specified genomic features.}
#' }
#'
#'
#' @import ggplot2
#'
#' @export
metaPlot <- function(bed_file = NULL,
                     group_names = NULL,
                     features_anno = NULL,
                     scale_region = FALSE,
                     cut_ratio = NULL,
                     facet = FALSE,
                     geom_density_bw = NULL,
                     geom_density_n = 512,
                     col_alpha = 0.5,
                     line_color = "#CC0033"){
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
    if(is.null(geom_density_bw)){
      density_layer <- geom_density(aes(x = rel_pos,fill = group_names,color = group_names),alpha = col_alpha,
                                    linewidth = line_width,n = geom_density_n)
    }else{
      density_layer <- geom_density(aes(x = rel_pos,fill = group_names,color = group_names),alpha = col_alpha,
                                    linewidth = line_width,
                                    bw = geom_density_bw, n = geom_density_n)
    }
  }else{
    if(is.null(geom_density_bw)){
      density_layer <- geom_density(aes(x = rel_pos),fill = line_color,color = line_color,alpha = col_alpha,
                                    linewidth = line_width,n = geom_density_n)
    }else{
      density_layer <- geom_density(aes(x = rel_pos),fill = line_color,color = line_color,alpha = col_alpha,
                                    linewidth = line_width,
                                    bw = geom_density_bw, n = geom_density_n)
    }

  }



  pdensity <-
    ggplot(plot_df) +
    density_layer +
    geom_vline(xintercept = c(1,2),lty = "dashed") +
    theme_classic(base_size = 12) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    scale_x_continuous(breaks = label.x,
                       labels = c("5'UTR","CDS","3'UTR")) +
    theme_bw() +
    theme(axis.line = element_line(arrow = arrow(length = unit(0.25,"cm"),
                                                 type = "closed")),
          axis.text = element_text(colour = "black"),
          strip.text = element_text(face = "bold",size = rel(1)),
          panel.grid = element_blank(),
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
