test_that(".select_longest_transcripts picks the transcript with the largest total exon width per gene", {
  gtf <- data.frame(
    gene_id = c("g1", "g1", "g2"),
    transcript_id = c("t1", "t2", "t3"),
    type = "exon",
    width = c(100, 300, 50),
    stringsAsFactors = FALSE
  )

  res <- metaplot:::.select_longest_transcripts(gtf)

  expect_equal(nrow(res), 2)
  expect_equal(res$transcript_id[res$gene_id == "g1"], "t2")
  expect_equal(res$transcript_id[res$gene_id == "g2"], "t3")
})

test_that(".select_longest_transcripts sums multiple exons per transcript before comparing", {
  gtf <- data.frame(
    gene_id = "g1",
    transcript_id = c("t1", "t1", "t2", "t2"),
    type = "exon",
    width = c(40, 40, 30, 30),
    stringsAsFactors = FALSE
  )

  res <- metaplot:::.select_longest_transcripts(gtf)

  expect_equal(nrow(res), 1)
  expect_equal(res$transcript_id, "t1")
  expect_equal(res$trans_len, 80)
})

test_that(".select_longest_transcripts is not sensitive to input row order", {
  # regression test: an earlier implementation picked the first-encountered
  # transcript per gene rather than the longest one
  gtf <- data.frame(
    gene_id = c("g1", "g1"),
    transcript_id = c("short", "long"),
    type = "exon",
    width = c(10, 500),
    stringsAsFactors = FALSE
  )

  res_forward <- metaplot:::.select_longest_transcripts(gtf)
  res_reversed <- metaplot:::.select_longest_transcripts(gtf[rev(seq_len(nrow(gtf))), ])

  expect_equal(res_forward$transcript_id, "long")
  expect_equal(res_reversed$transcript_id, "long")
})

test_that(".select_longest_transcripts ignores non-exon rows", {
  gtf <- data.frame(
    gene_id = c("g1", "g1", "g1"),
    transcript_id = c("t1", "t1", "t1"),
    type = c("exon", "CDS", "exon"),
    width = c(50, 1000, 50),
    stringsAsFactors = FALSE
  )

  res <- metaplot:::.select_longest_transcripts(gtf)

  expect_equal(res$trans_len, 100)
})
