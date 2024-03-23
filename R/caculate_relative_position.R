#' calculate relative position for each features
#'
#' @param bed_file bed format of peaks files.
#' @param features_anno features_anno object from prepare_features functions.
#' @param scale_region whether do scale for each feature length, default FALSE.
#' @param cut_ratio whether use custom ratio to plot feature regions,
#' for eaxample: c(0.1,0.8,0.1), default NULL.
#'
#' @return data frame
#' @export
calculate_relative_position <- function(bed_file = NULL,
                                        features_anno = NULL,
                                        scale_region = FALSE,
                                        cut_ratio = NULL){
  # load peaks data
  result <- try(rtracklayer::import.bed(bed_file),silent = T)

  # whether is bed file
  if(inherits(result, "try-error")){
    peaks <- utils::read.delim(bed_file,header = F)
    colnames(peaks)[1:6] <- c("seqnames","start","end","name","score","strand")
    peaks <- GenomicRanges::GRanges(peaks)
  }else{
    peaks <- rtracklayer::import.bed(bed_file)
  }

  # overlap
  ov <- IRanges::findOverlaps(query = peaks,subject = features_anno)

  lo <- cbind(as.data.frame(peaks[S4Vectors::queryHits(ov)]),
              as.data.frame(features_anno[S4Vectors::subjectHits(ov)]))

  # make unique names
  names(lo) <- make.names(names(lo),unique = T)

  # caculate repative positions
  lo <- lo %>%
    dplyr::mutate(p_mid = as.integer((start + end)/2),.after = end) %>%
    dplyr::filter(p_mid >= start.1 & p_mid <= end.1) %>%
    dplyr::mutate(rel_pos = dplyr::case_when(strand.1 == "+" ~ (p_mid - end.1 + tx_len)/f_len,
                                             strand.1 == "-" ~ (start.1 - p_mid + tx_len)/f_len)) %>%
    dplyr::mutate(rel_pos = dplyr::case_when(type == "CDS" ~ rel_pos + 1,
                                             type %in% c("UTR3") ~ rel_pos + 2,
                                             .default = rel_pos))

  # ==============================================================================
  # calculate scale factor for 5UTR and 3UTR
  # ==============================================================================
  # whether scale features to its length
  if(scale_region == TRUE){
    df_lo <- lo %>% dplyr::select(type,f_len) %>% unique() %>%
      dplyr::group_by(type) %>%
      dplyr::summarise(median_len = stats::median(f_len))

    utr5.sf <- df_lo[which(df_lo$type == "UTR5"),]$median_len / df_lo[which(df_lo$type == "CDS"),]$median_len
    utr3.sf <- df_lo[which(df_lo$type == "UTR3"),]$median_len / df_lo[which(df_lo$type == "CDS"),]$median_len

    # whether use custom defined ratios
    if(!is.null(cut_ratio)){
      utr5.sf <- cut_ratio[1]/cut_ratio[2]
      utr3.sf <- cut_ratio[3]/cut_ratio[2]
    }

    plot_df <- lo %>% dplyr::select(type,rel_pos) %>%
      dplyr::mutate(rel_pos = dplyr::case_when(type == "UTR5" ~ scales::rescale(rel_pos,to = c(1 - utr5.sf,1),from = c(0,1)),
                                               type == "UTR3" ~ scales::rescale(rel_pos,to = c(2,2 + utr3.sf),from = c(2,3)),
                                               .default = rel_pos))

    attr(plot_df,"utr5.sf") <- utr5.sf
    attr(plot_df,"utr3.sf") <- utr3.sf
  }else{
    plot_df <- lo %>% dplyr::select(type,rel_pos)
  }

  return(plot_df)
}
