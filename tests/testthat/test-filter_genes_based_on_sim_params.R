sim_params_without_batch <-
  list(
    approx_n_genes = 500,
    n_examples = 40,
    batch1_indices = 1:40,
    batch2_indices = 1:40,
    batch1_group1_indices = 1:20,
    batch1_group2_indices = 21:40,
    batch2_group1_indices = 1:20,
    batch2_group2_indices = 21:40,
    group1_indices = 1:20,
    group2_indices = 21:40,
    batch_difference_threshold = 1,
    n_true_diff_genes = 50,
    n_true_upreg_genes = 25,
    max_baseline_log2FC_between_groups = 0.25,
    avg_log2FC_between_groups = 1
  )

baseline_mat <- generate_baseline_exprs_matrix(n_genes = 500, n_examples = 40)

group_means_filt <- filter_genes_based_on_sim_params(simulated_datamat = baseline_mat,
                                                     sim_params_without_batch,
                                                     batch_effects_exist = FALSE)

test_that("filtering genes based on baseline group differences works", {
  expect_lte(max(abs(group_means_filt$group_fc)), 0.25)
})
