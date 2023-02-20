set.seed(123)
sim_params_without_batch <-
  list(
    approx_n_genes = 5000,
    n_examples = 60,
    group1_indices = 1:30,
    group2_indices = 31:60,
    n_true_diff_genes = 500,
    n_true_upreg_genes = 250,
    max_baseline_log2FC_between_groups = 0.25,
    avg_log2FC_between_groups = 1,
    train_split_prop = 0.70,
    batch_effects_exist = FALSE
  )

res_no_batch <- run_single_analysis(sim_params=sim_params_without_batch)

test_that("single analysis without batch effects works", {
  expect_equal(as.integer(res_no_batch$deseq_results$values[1]), 470)
})
