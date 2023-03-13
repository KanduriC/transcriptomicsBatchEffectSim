test_that("generate_baseline_exprs_matrix works", {
  expect_equal(nrow(generate_baseline_exprs_matrix(n_genes=4, n_examples=5)), 4)
  expect_equal(ncol(generate_baseline_exprs_matrix(n_genes=4, n_examples=5)), 5)
  expect_equal(nrow(generate_baseline_exprs_matrix(n_genes=1500, n_examples=5)), 1500)
})
