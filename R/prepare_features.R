#' prepare annotation features from GTF file
#'
#' @param gtf_file The gtf file, gzip file is accepted.
#' @param features_name The feature name in gtf file,
#' default for ensembl database gtf file
#' features_name= c("five_prime_utr","CDS","three_prime_utr"), if you use
#' ucsc database gtf file, you should change it to c("5UTR","CDS","3UTR").
#'
#' @return GRanges object
#' @export
prepare_features <- function(gtf_file = NULL,
                             features_name = c("five_prime_utr","CDS","three_prime_utr")){
  # assign names for features
  names(features_name) <- c("UTR5","CDS","UTR3")

  # load gtf file
  gtf <- rtracklayer::import.gff(gtf_file,format = "gtf") %>%
    data.frame() %>%
    dplyr::select(seqnames,start,end,width,strand,type,gene_id,gene_name,transcript_id)

  # filter longest transcript
  longest_trans <- gtf %>%
    dplyr::filter(type %in% "exon") %>%
    dplyr::group_by(gene_id,transcript_id) %>%
    dplyr::summarise(trans_len = sum(width)) %>%
    dplyr::slice_head(n = 1)

  # filter features
  features <- gtf %>% dplyr::filter(transcript_id %in% longest_trans$transcript_id) %>%
    dplyr::filter(type %in% features_name)

  # add new types
  # x = 1
  lapply(seq_along(features_name),function(x){
    tmp <- features %>% dplyr::filter(type %in% features_name[x]) %>%
      dplyr::mutate(type = names(features_name)[x])
  }) %>% do.call("rbind",.) -> features

  # get feature total length
  features_len <- features %>% dplyr::group_by(transcript_id,type) %>%
    dplyr::summarise(f_len = sum(width))

  # add transcript length
  features_len <- features_len %>%
    dplyr::left_join(y = longest_trans,by = "transcript_id")

  features <- features %>%
    dplyr::left_join(y = features_len,by = c("transcript_id","type"),
                     relationship = "many-to-many") %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::arrange(seqnames,start,end) %>%
    dplyr::relocate(f_len,.after = type)

  # caculate for positive and negtive strand
  pos_starnd <- features %>% dplyr::filter(strand == "+") %>%
    dplyr::group_by(transcript_id,type) %>%
    dplyr::mutate(tx_len = cumsum(width),.after = "width")

  neg_starnd <- features %>% dplyr::filter(strand == "-") %>%
    dplyr::arrange(seqnames,-start,-end) %>%
    dplyr::group_by(transcript_id,type) %>%
    dplyr::mutate(tx_len = cumsum(width),.after = "width")

  # combine pos and neg
  mer <- rbind(pos_starnd,neg_starnd) %>% GenomicRanges::GRanges()

  return(mer)
}
