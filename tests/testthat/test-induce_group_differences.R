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

group_means_filt <- get_descriptive_stats(simulated_datamat = baseline_mat,
                                          sim_params=sim_params_without_batch)

sim_results <- induce_group_differences(simulated_datamat=baseline_mat,
                                        group_means_filt=group_means_filt,
                                        sim_params=sim_params_without_batch)

group_stats <- compute_groupwise_stats(simulated_datamat=sim_results$simulated_datamat_filt,
                                       sim_params=sim_params_without_batch)
pre_fc <- group_means_filt %>% dplyr::select(group_fc) %>% dplyr::mutate(Gene=rownames(group_means_filt))
post_fc <- group_stats %>% dplyr::select(group_fc) %>% dplyr::mutate(Gene=rownames(group_stats))
fc <- dplyr::left_join(pre_fc, post_fc, by="Gene")
n_genes_not_diff <- length(which(fc$group_fc.x==fc$group_fc.y))
n_genes_diff <- nrow(sim_results$simulated_datamat_filt)-n_genes_not_diff

true_diff_genes <- fc[sim_results$diff_exp_indices,]

test_that("inducing group differences work", {
  expect_equal(n_genes_diff, sim_params_without_batch$n_true_diff_genes)
  expect_equal(round(mean(abs(true_diff_genes$group_fc.y))), sim_params_without_batch$avg_log2FC_between_groups)
})
