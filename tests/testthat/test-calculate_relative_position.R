make_features <- function(strand, start, end, f_len, tx_len, type = "CDS"){
  GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = start, end = end),
    strand = strand,
    type = type,
    f_len = f_len,
    tx_len = tx_len
  )
}

test_that("relative position is computed correctly on the + strand", {
  features_anno <- make_features("+", start = 100, end = 199, f_len = 100, tx_len = 100)
  peaks <- data.frame(seqnames = "chr1", start = 149, end = 150, strand = "+")

  res <- calculate_relative_position(bed_file = peaks, features_anno = features_anno)

  expect_equal(nrow(res), 1)
  # p_mid = 149; (p_mid - feat_end + tx_len) / f_len = (149-199+100)/100 = 0.5, +1 for CDS
  expect_equal(res$rel_pos, 1.5)
})

test_that("relative position is computed correctly on the - strand", {
  features_anno <- make_features("-", start = 500, end = 599, f_len = 100, tx_len = 100)
  peaks <- data.frame(seqnames = "chr1", start = 549, end = 550, strand = "-")

  res <- calculate_relative_position(bed_file = peaks, features_anno = features_anno)

  expect_equal(nrow(res), 1)
  # p_mid = 549; (feat_start - p_mid + tx_len) / f_len = (500-549+100)/100 = 0.51, +1 for CDS
  expect_equal(res$rel_pos, 1.51)
})

test_that("calculate_relative_position accepts file paths, data.frames, and GRanges", {
  features_anno <- make_features("+", start = 100, end = 199, f_len = 100, tx_len = 100)

  peaks_df <- data.frame(seqnames = "chr1", start = 149, end = 150, strand = "+")
  peaks_gr <- GenomicRanges::GRanges(peaks_df)

  res_df <- calculate_relative_position(bed_file = peaks_df, features_anno = features_anno)
  res_gr <- calculate_relative_position(bed_file = peaks_gr, features_anno = features_anno)

  expect_equal(res_df$rel_pos, res_gr$rel_pos)
})

test_that("peaks overlapping unstranded features are dropped with a warning", {
  features_anno <- make_features("*", start = 100, end = 199, f_len = 100, tx_len = 100)
  peaks <- data.frame(seqnames = "chr1", start = 149, end = 150, strand = "*")

  expect_warning(
    res <- calculate_relative_position(bed_file = peaks, features_anno = features_anno),
    "unstranded"
  )
  expect_equal(nrow(res), 0)
})

test_that("calculate_relative_position validates its inputs", {
  features_anno <- make_features("+", start = 100, end = 199, f_len = 100, tx_len = 100)
  peaks <- data.frame(seqnames = "chr1", start = 149, end = 150, strand = "+")

  expect_error(calculate_relative_position(bed_file = NULL, features_anno = features_anno), "bed_file")
  expect_error(calculate_relative_position(bed_file = peaks, features_anno = NULL), "features_anno")
  expect_error(calculate_relative_position(bed_file = peaks, features_anno = GenomicRanges::GRanges()), "required")
  expect_error(
    calculate_relative_position(bed_file = peaks, features_anno = features_anno, cut_ratio = c(0.1, 0.9)),
    "length 3"
  )
})

test_that("column identity is preserved even with extra/reordered peak metadata", {
  # regression test: an earlier implementation relied on cbind() + make.names()
  # to disambiguate duplicate column names by position, which could silently
  # misalign columns if peak metadata columns changed
  features_anno <- make_features("+", start = 100, end = 199, f_len = 100, tx_len = 100)
  peaks <- data.frame(
    seqnames = "chr1", start = 149, end = 150, name = "peak1",
    score = 0, strand = "+", extra_col = "unexpected"
  )

  res <- calculate_relative_position(bed_file = peaks, features_anno = features_anno)

  expect_equal(res$rel_pos, 1.5)
})

make_anno <- function(utr5 = TRUE, utr3 = TRUE){
  start <- 51; end <- 150
  type <- "CDS"; f_len <- 100; tx_len <- 100

  if(utr5){
    start <- c(1, start); end <- c(50, end)
    type <- c("UTR5", type); f_len <- c(50, f_len); tx_len <- c(50, tx_len)
  }
  if(utr3){
    start <- c(start, 151); end <- c(end, 250)
    type <- c(type, "UTR3"); f_len <- c(f_len, 100); tx_len <- c(tx_len, 100)
  }

  GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = start, end = end),
    strand = "+",
    type = type,
    f_len = f_len,
    tx_len = tx_len
  )
}

test_that("scale factors are derived from the whole annotation, not from the peaks", {
  features_anno <- make_anno()

  # a peak that only touches the CDS, and one that only touches the 3'UTR
  peak_cds  <- data.frame(seqnames = "chr1", start = 99, end = 100, strand = "+")
  peak_utr3 <- data.frame(seqnames = "chr1", start = 199, end = 200, strand = "+")

  res_cds <- calculate_relative_position(peak_cds, features_anno, scale_region = TRUE)
  res_utr3 <- calculate_relative_position(peak_utr3, features_anno, scale_region = TRUE)

  # both peak sets must land on the same x-axis
  expect_equal(attr(res_cds, "utr5.sf"), 0.5)
  expect_equal(attr(res_cds, "utr3.sf"), 1)
  expect_equal(attr(res_cds, "utr5.sf"), attr(res_utr3, "utr5.sf"))
  expect_equal(attr(res_cds, "utr3.sf"), attr(res_utr3, "utr3.sf"))
})

test_that("a missing UTR collapses to zero width instead of failing", {
  features_anno <- make_anno(utr5 = FALSE, utr3 = FALSE)
  peaks <- data.frame(seqnames = "chr1", start = 99, end = 100, strand = "+")

  res <- calculate_relative_position(peaks, features_anno, scale_region = TRUE)

  expect_equal(attr(res, "utr5.sf"), 0)
  expect_equal(attr(res, "utr3.sf"), 0)
  expect_equal(nrow(res), 1)
})

test_that("precomputed scale factors are honoured", {
  features_anno <- make_anno()
  peaks <- data.frame(seqnames = "chr1", start = 99, end = 100, strand = "+")

  res <- calculate_relative_position(peaks, features_anno, scale_region = TRUE,
                                     scale_factors = c(utr5.sf = 0.2, utr3.sf = 0.4))

  expect_equal(attr(res, "utr5.sf"), 0.2)
  expect_equal(attr(res, "utr3.sf"), 0.4)
})

test_that("cut_ratio overrides the estimated scale factors", {
  features_anno <- make_anno()
  peaks <- data.frame(seqnames = "chr1", start = 99, end = 100, strand = "+")

  res <- calculate_relative_position(peaks, features_anno, scale_region = TRUE,
                                     cut_ratio = c(0.1, 0.8, 0.1))

  expect_equal(attr(res, "utr5.sf"), 0.1/0.8)
  expect_equal(attr(res, "utr3.sf"), 0.1/0.8)
})

test_that("cut_ratio is reported as ignored when scale_region is FALSE", {
  features_anno <- make_anno()
  peaks <- data.frame(seqnames = "chr1", start = 99, end = 100, strand = "+")

  expect_warning(
    calculate_relative_position(peaks, features_anno, scale_region = FALSE,
                                cut_ratio = c(0.1, 0.8, 0.1)),
    "ignored"
  )
})

test_that("cut_ratio must contain three positive values", {
  features_anno <- make_anno()
  peaks <- data.frame(seqnames = "chr1", start = 99, end = 100, strand = "+")

  expect_error(
    calculate_relative_position(peaks, features_anno, cut_ratio = c(0.1, 0, 0.1)),
    "positive"
  )
  expect_error(
    calculate_relative_position(peaks, features_anno, cut_ratio = "a"),
    "numeric"
  )
})

test_that("a missing peak file is reported clearly", {
  features_anno <- make_anno()

  expect_error(
    calculate_relative_position("no_such_file.bed", features_anno),
    "does not exist"
  )
})

test_that("a peak data frame without a strand column is treated as unstranded peaks", {
  features_anno <- make_anno()
  peaks <- data.frame(seqnames = "chr1", start = 100, end = 101)

  res <- calculate_relative_position(peaks, features_anno)

  expect_equal(nrow(res), 1)
  expect_equal(res$rel_pos, 1.5)
})

test_that("a peak data frame without coordinates is reported clearly", {
  features_anno <- make_anno()
  peaks <- data.frame(chr = "chr1", from = 99, to = 100)

  expect_error(
    calculate_relative_position(peaks, features_anno),
    "seqnames"
  )
})
