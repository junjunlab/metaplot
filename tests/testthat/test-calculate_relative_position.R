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
