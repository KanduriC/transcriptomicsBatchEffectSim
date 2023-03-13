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

sim_params_with_batch <-
  list(
    approx_n_genes = 500,
    n_examples = 40,
    batch1_indices = 1:20,
    batch2_indices = 21:40,
    batch1_group1_indices = 1:10,
    batch1_group2_indices = 11:20,
    batch2_group1_indices = 21:30,
    batch2_group2_indices = 31:40,
    batch_difference_threshold = 1,
    n_true_diff_genes = 50,
    n_true_upreg_genes = 25,
    max_baseline_log2FC_between_groups = 0.25,
    avg_log2FC_between_groups = 1
  )

baseline_mat <- generate_baseline_exprs_matrix(n_genes = 500, n_examples = 40)

group_means_filt <- get_descriptive_stats(simulated_datamat = baseline_mat,
                                          sim_params=sim_params_without_batch)

sim_results <- induce_group_differences(simulated_datamat = baseline_mat,
                                        group_means_filt = group_means_filt,
                                        sim_params = sim_params_without_batch)
data_objs_list <- prepare_metadata_and_exprsmat(sim_params = sim_params_without_batch,
                              simulated_datamat = sim_results$simulated_datamat_filt)


test_that("preparing metadata works", {
  expect_equal(data_objs_list$batch, rep("batch_1", length(sim_params_without_batch$batch1_indices)))
})
