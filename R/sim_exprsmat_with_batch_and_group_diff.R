simulate_exprsmat_with_batch_and_group_differences <- function(sim_params){
  baseline_exprs_mat <- generate_baseline_exprs_matrix(
    n_genes=sim_params$approx_n_genes,
    n_examples=sim_params$n_examples)
  batch_induced_exprs_mat <- induce_batch_differences(baseline_exprs_mat=baseline_exprs_mat,
                           sim_params=sim_params)
  group_means_filt <- get_descriptive_stats(simulated_datamat = batch_induced_exprs_mat,
                                            sim_params=sim_params)
  sim_results <- induce_group_differences(simulated_datamat=batch_induced_exprs_mat,
                                          group_means_filt=group_means_filt,
                                          sim_params=sim_params,
                                          batch_effects_exist=TRUE)
  data_objs_list <- prepare_metadata_and_exprsmat(
    sim_params = sim_params,
    simulated_datamat = sim_results$batch_simulated_datamat_filt,
    batch_effects_exist=TRUE)
  data_objs_list$group_means_filt <- group_means_filt
  data_objs_list$post_sim_group_means <- compute_groupwise_stats(
    simulated_datamat = data_objs_list$exprs_mat, sim_params = sim_params)
  data_objs_list$diff_exp_indices <- sim_results$diff_exp_indices
  return(data_objs_list)
}
