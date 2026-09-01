#' calculate relative position for each features
#'
#' @param bed_file bed format of peaks files, or a `GRanges`/`data.frame` of
#' peaks already loaded into R (must contain `seqnames`, `start`, `end`
#' columns).
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
  # ==============================================================================
  # validate inputs
  # ==============================================================================
  if(is.null(bed_file)){
    stop("`bed_file` must be supplied.",call. = FALSE)
  }

  if(is.null(features_anno) || !methods::is(features_anno,"GRanges")){
    stop("`features_anno` must be a GRanges object from `prepare_features()`.",
        call. = FALSE)
  }

  if(!is.null(cut_ratio) && length(cut_ratio) != 3){
    stop("`cut_ratio` must be a numeric vector of length 3, e.g. c(0.1, 0.8, 0.1).",
        call. = FALSE)
  }

  required_feat_cols <- c("type","f_len","tx_len")
  missing_cols <- setdiff(required_feat_cols,names(S4Vectors::mcols(features_anno)))
  if(length(missing_cols) > 0){
    stop("`features_anno` is missing required column(s): ",
        paste(missing_cols,collapse = ", "),call. = FALSE)
  }

  # load peaks (file path / data.frame / GRanges all supported)
  peaks <- .load_peaks(bed_file)

  # ==============================================================================
  # overlap
  # ==============================================================================
  ov <- IRanges::findOverlaps(query = peaks,subject = features_anno)

  # explicitly select and rename the columns we need, rather than cbind-ing
  # full data frames and disambiguating duplicate names by position, so
  # renamed/reordered upstream columns can't silently misalign the join
  peaks_df <- as.data.frame(peaks)[S4Vectors::queryHits(ov),c("start","end"),drop = FALSE]
  names(peaks_df) <- c("peak_start","peak_end")

  feat_df <- as.data.frame(features_anno)[S4Vectors::subjectHits(ov),
                                          c("start","end","strand","type","f_len","tx_len"),
                                          drop = FALSE]
  names(feat_df) <- c("feat_start","feat_end","feat_strand","type","f_len","tx_len")

  lo <- cbind(peaks_df,feat_df)

  # caculate repative positions
  lo <- lo %>%
    dplyr::mutate(p_mid = as.integer((peak_start + peak_end)/2)) %>%
    dplyr::filter(p_mid >= feat_start & p_mid <= feat_end) %>%
    dplyr::mutate(rel_pos = dplyr::case_when(feat_strand == "+" ~ (p_mid - feat_end + tx_len)/f_len,
                                             feat_strand == "-" ~ (feat_start - p_mid + tx_len)/f_len,
                                             .default = NA_real_))

  if(any(is.na(lo$rel_pos))){
    warning("Some peaks overlap features with unstranded (\"*\") strand; ",
           "their relative position could not be computed and were dropped.",
           call. = FALSE)
    lo <- lo %>% dplyr::filter(!is.na(rel_pos))
  }

  lo <- lo %>%
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

#' Load peaks as a GRanges object
#'
#' Accepts a BED/narrowPeak file path, a `data.frame`, or a `GRanges` object
#' and normalizes it to `GRanges`.
#'
#' @param bed_file file path, `data.frame`, or `GRanges`.
#'
#' @return GRanges object
#' @keywords internal
#' @noRd
.load_peaks <- function(bed_file){
  if(methods::is(bed_file,"GRanges")){
    return(bed_file)
  }

  if(is.data.frame(bed_file)){
    return(GenomicRanges::GRanges(bed_file))
  }

  if(!is.character(bed_file) || length(bed_file) != 1){
    stop("`bed_file` must be a file path, data.frame, or GRanges object.",
        call. = FALSE)
  }

  # load peaks data
  result <- try(rtracklayer::import.bed(bed_file),silent = TRUE)

  # whether is bed file
  if(inherits(result,"try-error")){
    peaks <- utils::read.delim(bed_file,header = FALSE)

    if(ncol(peaks) < 6){
      stop("BED file must have at least 6 columns ",
          "(seqnames, start, end, name, score, strand).",call. = FALSE)
    }

    colnames(peaks)[1:6] <- c("seqnames","start","end","name","score","strand")
    peaks <- GenomicRanges::GRanges(peaks)
  }else{
    peaks <- result
  }

  peaks
}
