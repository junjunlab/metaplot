#' calculate relative position for each features
#'
#' @param bed_file bed format of peaks files, or a `GRanges`/`data.frame` of
#' peaks already loaded into R (must contain `seqnames`, `start`, `end`
#' columns).
#' @param features_anno features_anno object from prepare_features functions.
#' @param scale_region whether do scale for each feature length, default FALSE.
#' @param cut_ratio whether use custom ratio to plot feature regions,
#' for eaxample: c(0.1,0.8,0.1), default NULL.
#' @param scale_factors optional named numeric vector with `utr5.sf` and
#' `utr3.sf` entries. When `NULL` (the default) the factors are estimated from
#' `features_anno`. Pass a precomputed value when several peak files must share
#' exactly the same x-axis scaling (this is what `metaPlot()` does).
#'
#' @return data frame
#' @export
calculate_relative_position <- function(bed_file = NULL,
                                        features_anno = NULL,
                                        scale_region = FALSE,
                                        cut_ratio = NULL,
                                        scale_factors = NULL){
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

  if(!is.null(cut_ratio)){
    if(!is.numeric(cut_ratio) || length(cut_ratio) != 3){
      stop("`cut_ratio` must be a numeric vector of length 3, e.g. c(0.1, 0.8, 0.1).",
          call. = FALSE)
    }
    if(anyNA(cut_ratio) || any(cut_ratio <= 0)){
      stop("`cut_ratio` must contain three positive, non-missing values.",
          call. = FALSE)
    }
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

  if(length(ov) == 0){
    warning("No peak overlaps the annotated features; returning an empty result.",
            call. = FALSE)
  }

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
           "their relative position could not be computed and was dropped.",
           call. = FALSE)
    lo <- lo %>% dplyr::filter(!is.na(rel_pos))
  }

  lo <- lo %>%
    dplyr::mutate(rel_pos = dplyr::case_when(type == "CDS" ~ rel_pos + 1,
                                             type == "UTR3" ~ rel_pos + 2,
                                             .default = rel_pos))

  # ==============================================================================
  # calculate scale factor for 5UTR and 3UTR
  # ==============================================================================
  # whether scale features to its length
  if(isTRUE(scale_region)){
    # always derive the factors from the full annotation, never from the subset
    # of features that happens to overlap these peaks, otherwise every peak
    # file would end up on its own x-axis
    if(is.null(scale_factors)){
      scale_factors <- .feature_scale_factors(features_anno = features_anno,
                                              cut_ratio = cut_ratio)
    }

    utr5.sf <- unname(scale_factors[["utr5.sf"]])
    utr3.sf <- unname(scale_factors[["utr3.sf"]])

    plot_df <- lo %>% dplyr::select(type,rel_pos) %>%
      dplyr::mutate(rel_pos = dplyr::case_when(type == "UTR5" ~ scales::rescale(rel_pos,to = c(1 - utr5.sf,1),from = c(0,1)),
                                               type == "UTR3" ~ scales::rescale(rel_pos,to = c(2,2 + utr3.sf),from = c(2,3)),
                                               .default = rel_pos))

    attr(plot_df,"utr5.sf") <- utr5.sf
    attr(plot_df,"utr3.sf") <- utr3.sf
  }else{
    if(!is.null(cut_ratio)){
      warning("`cut_ratio` is ignored when `scale_region = FALSE`. ",
              "Set `scale_region = TRUE` to use it.",call. = FALSE)
    }

    plot_df <- lo %>% dplyr::select(type,rel_pos)
  }

  # the row names inherited from the feature table are meaningless here
  rownames(plot_df) <- NULL

  return(plot_df)
}

#' Estimate the 5'UTR/CDS and 3'UTR/CDS width ratios
#'
#' The ratios are derived from the **whole** annotation (one median length per
#' feature type) so that every peak file analysed against the same annotation
#' is placed on an identical x-axis.
#'
#' @param features_anno A `GRanges` object from `prepare_features()`.
#' @param cut_ratio Optional numeric vector of length 3
#' (`5'UTR`, `CDS`, `3'UTR`) overriding the estimated ratios.
#'
#' @return A named numeric vector with elements `utr5.sf` and `utr3.sf`.
#' @keywords internal
#' @noRd
.feature_scale_factors <- function(features_anno,cut_ratio = NULL){
  if(!is.null(cut_ratio)){
    return(c(utr5.sf = cut_ratio[1]/cut_ratio[2],
             utr3.sf = cut_ratio[3]/cut_ratio[2]))
  }

  feat <- S4Vectors::mcols(features_anno)
  feat_len <- unique(data.frame(type = as.character(feat$type),
                                f_len = as.numeric(feat$f_len)))
  med <- stats::aggregate(f_len ~ type,data = feat_len,FUN = stats::median)
  med_len <- stats::setNames(med$f_len,med$type)

  if(!"CDS" %in% names(med_len)){
    stop("`features_anno` contains no CDS feature, so `scale_region = TRUE` ",
         "cannot be used. Either supply `cut_ratio` or set `scale_region = FALSE`.",
         call. = FALSE)
  }

  cds_len <- med_len[["CDS"]]
  if(!is.finite(cds_len) || cds_len <= 0){
    stop("The median CDS length in `features_anno` is not a positive number; ",
         "`scale_region = TRUE` cannot be used.",call. = FALSE)
  }

  # a missing UTR collapses that region to zero width instead of failing
  utr5 <- if("UTR5" %in% names(med_len)) unname(med_len[["UTR5"]])/cds_len else 0
  utr3 <- if("UTR3" %in% names(med_len)) unname(med_len[["UTR3"]])/cds_len else 0

  c(utr5.sf = utr5,utr3.sf = utr3)
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
    return(.as_peak_granges(bed_file,source = "`bed_file`"))
  }

  if(!is.character(bed_file) || length(bed_file) != 1 || is.na(bed_file)){
    stop("`bed_file` must be a single file path, a data.frame, or a GRanges object.",
        call. = FALSE)
  }

  if(!file.exists(bed_file)){
    stop("Peak file does not exist: ",bed_file,call. = FALSE)
  }

  # load peaks data
  result <- try(rtracklayer::import.bed(bed_file),silent = TRUE)

  # whether is bed file
  if(inherits(result,"try-error")){
    peaks <- try(utils::read.delim(bed_file,header = FALSE),silent = TRUE)

    if(inherits(peaks,"try-error")){
      stop("Could not read ",bed_file," as a BED file: ",
           conditionMessage(attr(peaks,"condition")),call. = FALSE)
    }

    if(ncol(peaks) < 6){
      stop("BED file must have at least 6 columns ",
           "(seqnames, start, end, name, score, strand).",call. = FALSE)
    }

    colnames(peaks)[1:6] <- c("seqnames","start","end","name","score","strand")
    peaks <- .as_peak_granges(peaks,source = bed_file)
  }else{
    peaks <- result
  }

  peaks
}

#' Coerce a peak data frame to GRanges with a readable error message
#'
#' @param df A data frame of peaks.
#' @param source Description of where `df` came from, used in error messages.
#'
#' @return GRanges object
#' @keywords internal
#' @noRd
.as_peak_granges <- function(df,source = "peaks"){
  required <- c("seqnames","start","end")
  missing_cols <- setdiff(required,names(df))

  if(length(missing_cols) > 0){
    stop(source," must contain the column(s) ",paste(missing_cols,collapse = ", "),
         "; got: ",paste(names(df),collapse = ", "),call. = FALSE)
  }

  # GRanges() would silently leave peaks unstranded when the column is absent
  if(!"strand" %in% names(df)){
    df$strand <- "*"
  }

  # BED files commonly use "." for unstranded, which GRanges() rejects
  strand <- as.character(df$strand)
  strand[is.na(strand) | !strand %in% c("+","-","*")] <- "*"
  df$strand <- strand

  GenomicRanges::GRanges(df)
}
