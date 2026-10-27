make_anno <- function(){
  GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(1, 51, 151), end = c(50, 150, 250)),
    strand = "+",
    type = c("UTR5", "CDS", "UTR3"),
    f_len = c(50, 100, 100),
    tx_len = c(50, 100, 100)
  )
}

make_peaks <- function(offset = 0){
  data.frame(
    seqnames = "chr1",
    start = c(10, 60, 70, 80, 90, 110, 120, 160, 170, 180) + offset,
    end = c(11, 61, 71, 81, 91, 111, 121, 161, 171, 181) + offset,
    strand = "+"
  )
}

test_that("metaPlot requires bed_file", {
  expect_error(metaPlot(features_anno = make_anno()), "bed_file")
})

test_that("metaPlot accepts a single data.frame instead of a list", {
  res <- metaPlot(bed_file = make_peaks(), features_anno = make_anno())

  expect_type(res, "list")
  expect_true(all(c("data", "plot") %in% names(res)))
  expect_true(all(res$data$group_names == "sample1"))
})

test_that("metaPlot rejects a group_names vector of the wrong length", {
  expect_error(
    metaPlot(bed_file = list(make_peaks(), make_peaks()),
             group_names = "only_one",
             features_anno = make_anno()),
    "must match"
  )
})

test_that("metaPlot warns when cut_ratio is given without scale_region", {
  expect_warning(
    metaPlot(bed_file = list(make_peaks()), features_anno = make_anno(),
             cut_ratio = c(0.1, 0.8, 0.1)),
    "ignored"
  )
})

test_that("metaPlot fails clearly when no peak overlaps the annotation", {
  far_away <- data.frame(seqnames = "chr2", start = 1:5, end = 2:6, strand = "+")

  expect_error(
    suppressWarnings(
      metaPlot(bed_file = list(far_away), features_anno = make_anno())
    ),
    "No peak falls within"
  )
})

test_that("samples touching different feature types share one x-axis", {
  # regression test: scale factors used to be estimated from the features that
  # happened to overlap each peak file, so a CDS-only and a 3'UTR-only sample
  # were rescaled differently and the x-axis labels no longer lined up
  cds_only  <- data.frame(seqnames = "chr1", start = c(60, 70, 80), end = c(61, 71, 81), strand = "+")
  utr3_only <- data.frame(seqnames = "chr1", start = c(160, 170, 180), end = c(161, 171, 181), strand = "+")

  res <- metaPlot(bed_file = list(cds_only, utr3_only),
                  group_names = c("cds", "utr3"),
                  features_anno = make_anno(),
                  scale_region = TRUE)

  expect_false(anyNA(res$data$rel_pos))
  expect_true(all(is.finite(res$data$rel_pos)))
  # CDS positions live in [1, 2], 3'UTR positions in [2, 2 + utr3.sf] = [2, 3]
  expect_true(all(res$data$rel_pos[res$data$group_names == "cds"] >= 1))
  expect_true(all(res$data$rel_pos[res$data$group_names == "cds"] <= 2))
  expect_true(all(res$data$rel_pos[res$data$group_names == "utr3"] >= 2))
  expect_true(all(res$data$rel_pos[res$data$group_names == "utr3"] <= 3))
})

test_that("metaPlot supports faceting", {
  res <- metaPlot(bed_file = list(make_peaks(), make_peaks(offset = 5)),
                  group_names = c("a", "b"),
                  features_anno = make_anno(),
                  facet = TRUE)

  expect_s3_class(res$plot$facet, "FacetWrap")
})
