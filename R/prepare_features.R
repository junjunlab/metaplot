#' prepare annotation features from GTF file
#'
#' @param gtf_file The gtf file, gzip file is accepted.
#' @param features_name The feature name in gtf file,
#' default for ensembl database gtf file
#' features_name= c("five_prime_utr","CDS","three_prime_utr"), if you use
#' ucsc database gtf file, you should change it to c("5UTR","CDS","3UTR").
#' Must be a character vector of length 3, given in the order 5'UTR, CDS, 3'UTR.
#'
#' @return GRanges object
#' @export
prepare_features <- function(gtf_file = NULL,
                             features_name = c("five_prime_utr","CDS","three_prime_utr")){
  # ==============================================================================
  # validate inputs
  # ==============================================================================
  if(is.null(gtf_file)){
    stop("`gtf_file` must be supplied.",call. = FALSE)
  }

  if(!is.character(gtf_file) || length(gtf_file) != 1 || is.na(gtf_file)){
    stop("`gtf_file` must be a single file path.",call. = FALSE)
  }

  if(!file.exists(gtf_file)){
    stop("GTF file does not exist: ",gtf_file,call. = FALSE)
  }

  if(!is.character(features_name) || length(features_name) != 3 || anyNA(features_name)){
    stop("`features_name` must be a character vector of length 3, in the order ",
         "5'UTR, CDS, 3'UTR (e.g. c(\"5UTR\", \"CDS\", \"3UTR\")).",call. = FALSE)
  }

  # assign names for features
  names(features_name) <- c("UTR5","CDS","UTR3")

  # load gtf file
  gtf <- rtracklayer::import.gff(gtf_file,format = "gtf") %>%
    data.frame()

  if(!"transcript_id" %in% names(gtf)){
    stop("`gtf_file` has no `transcript_id` attribute; ",
         "transcript-level features cannot be built.",call. = FALSE)
  }

  # keep the longest-transcript selection meaningful even when the annotation
  # omits gene_id, and avoid a cryptic dplyr error for a missing gene_name
  if(!"gene_id" %in% names(gtf)){
    gtf$gene_id <- gtf$transcript_id
  }

  if(!"gene_name" %in% names(gtf)){
    gtf$gene_name <- gtf$gene_id
  }

  gtf <- gtf %>%
    dplyr::select(seqnames,start,end,width,strand,type,gene_id,gene_name,transcript_id)

  # filter longest transcript
  longest_trans <- .select_longest_transcripts(gtf)

  # filter features
  features <- gtf %>% dplyr::filter(transcript_id %in% longest_trans$transcript_id) %>%
    dplyr::filter(type %in% features_name)

  if(nrow(features) == 0){
    stop("No feature matched `features_name` (",
         paste(features_name,collapse = ", "),") in ",gtf_file,
         ". Check the `type` labels used by this annotation.",call. = FALSE)
  }

  # add new types
  lapply(seq_along(features_name),function(x){
    tmp <- features %>% dplyr::filter(type %in% features_name[x]) %>%
      dplyr::mutate(type = names(features_name)[x])
  }) %>% do.call("rbind",.) -> features

  # get feature total length
  features_len <- features %>% dplyr::group_by(transcript_id,type) %>%
    dplyr::summarise(f_len = sum(width),.groups = "drop")

  # add transcript length
  features_len <- features_len %>%
    dplyr::left_join(y = dplyr::select(longest_trans,transcript_id,trans_len),
                     by = "transcript_id")

  features <- features %>%
    dplyr::left_join(y = features_len,by = c("transcript_id","type"),
                     relationship = "many-to-one") %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::arrange(seqnames,start,end) %>%
    dplyr::relocate(f_len,.after = type)

  # caculate for positive and negtive strand
  pos_starnd <- features %>% dplyr::filter(strand == "+") %>%
    dplyr::arrange(seqnames,start,end) %>%
    dplyr::group_by(transcript_id,type) %>%
    dplyr::mutate(tx_len = cumsum(width),.after = "width")

  neg_starnd <- features %>% dplyr::filter(strand == "-") %>%
    dplyr::arrange(seqnames,-start,-end) %>%
    dplyr::group_by(transcript_id,type) %>%
    dplyr::mutate(tx_len = cumsum(width),.after = "width")

  # combine pos and neg
  mer <- rbind(pos_starnd,neg_starnd) %>%
    dplyr::ungroup() %>%
    GenomicRanges::GRanges()

  return(mer)
}

#' Select the longest transcript per gene by total exon length
#'
#' Sums exon widths per transcript, then keeps a single transcript per gene
#' with the largest total exon length (ties broken by the first row
#' encountered).
#'
#' @param gtf_df A data frame with at least `gene_id`, `transcript_id`,
#' `type` and `width` columns (as produced from an imported GTF).
#'
#' @return A data frame with one row per gene: `gene_id`, `transcript_id`,
#' `trans_len`.
#' @keywords internal
#' @noRd
.select_longest_transcripts <- function(gtf_df){
  gtf_df %>%
    dplyr::filter(type %in% "exon") %>%
    dplyr::group_by(gene_id,transcript_id) %>%
    dplyr::summarise(trans_len = sum(width),.groups = "drop") %>%
    dplyr::group_by(gene_id) %>%
    dplyr::slice_max(trans_len,n = 1,with_ties = FALSE) %>%
    dplyr::ungroup()
}
