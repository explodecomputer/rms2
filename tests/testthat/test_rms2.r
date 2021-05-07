context("rms2")
library(rms2)

test_that("rms2", {

	x <- rms2$new("ukb-b-19953")
	x$extract_gwashits()
	x$scan_rsid(x$gwashits$rsid[1])
	# x$coloc_scan(x$gwashits$rsid[1])
	expect_warning(x$colocalise_rsid(x$gwashits$rsid[1], x$rsid_scan[[x$gwashits$rsid[1]]]$id[19]))
	expect_true(is.list(x$coloc_data))
	# x$plot_coloc(x$gwashits$rsid[1], x$rsid_scan[[x$gwashits$rsid[1]]]$id[19])
})
