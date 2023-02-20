set.seed(123)
sim_params <-
  list(
    approx_n_genes = 5000,
    n_examples = 60,
    batch1_indices = 1:30,
    batch2_indices = 31:60,
    batch1_group1_indices = 1:20,
    batch1_group2_indices = 21:30,
    batch2_group1_indices = 31:40,
    batch2_group2_indices = 41:60,
    group1_indices = c(1:20, 31:40),
    group2_indices = c(21:30, 41:60),
    batch_difference_threshold = 1,
    n_true_diff_genes = 500,
    n_true_upreg_genes = 250,
    max_baseline_log2FC_between_groups = 0.25,
    avg_log2FC_between_groups = 1,
    train_split_prop = 0.70,
    batch_effects_exist = TRUE
  )

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

res_1 <- rep_simulations(sim_params_list = sim_params, n_times=2)
res_2 <- rep_simulations(sim_params_list = sim_params_without_batch, n_times=2)

deseq_results <- join_results(res_1, res_2, "deseq_results")
logistic_results <- join_results(res_1, res_2, "logistic_results")
# box_plot(deseq_results, title = "Recovering true signals through univariate statistical testing")
# box_plot(logistic_results, title = "Performance of predictive model trained on different types of datasets")

test_that("generating results of replicated simulations works", {
  expect_equal(nrow(deseq_results), 18)
  expect_equal(nrow(logistic_results), 6)
})
