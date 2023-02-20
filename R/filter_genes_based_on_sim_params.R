fit_loess <- function(mean_and_disperison) {
  mean_disp_loess <-
    loess(
      mean_and_disperison$dispersion ~ mean_and_disperison$gene_means,
      mean_and_disperison
    )
  return(mean_disp_loess)
}

get_group_means_and_diffs <-
  function(count_matrix,
           batch1_group1_indices,
           batch1_group2_indices,
           batch2_group1_indices,
           batch2_group2_indices,
           batch1_indices,
           batch2_indices,
           group1_indices,
           group2_indices) {
    count_matrix <- log2(count_matrix + 0.01)
    group_means <-
      data.frame(
        rowMeans(count_matrix[, batch1_group1_indices]),
        rowMeans(count_matrix[, batch1_group2_indices]),
        rowMeans(count_matrix[, batch2_group1_indices]),
        rowMeans(count_matrix[, batch2_group2_indices]),
        rowMeans(count_matrix[, batch1_indices]),
        rowMeans(count_matrix[, batch2_indices]),
        rowMeans(count_matrix[, group1_indices]),
        rowMeans(count_matrix[, group2_indices])
      )
    colnames(group_means) <-
      c("b1_g1_mean",
        "b1_g2_mean",
        "b2_g1_mean",
        "b2_g2_mean",
        "b1_mean",
        "b2_mean",
        "g1_mean",
        "g2_mean")
    group_means <-
      group_means %>% dplyr::mutate(
        b1_groups_fc = b1_g2_mean - b1_g1_mean,
        b2_groups_fc = b2_g2_mean - b2_g1_mean,
        batch_fc = b2_mean - b1_mean,
        group_fc = g2_mean - g1_mean
      )
    return(group_means)

  }

compute_groupwise_stats <- function(simulated_datamat, sim_params) {
  group_means <-
    get_group_means_and_diffs(
      count_matrix = simulated_datamat,
      batch1_group1_indices = sim_params$batch1_group1_indices,
      batch1_group2_indices = sim_params$batch1_group2_indices,
      batch2_group1_indices = sim_params$batch2_group1_indices,
      batch2_group2_indices = sim_params$batch2_group2_indices,
      batch1_indices = sim_params$batch1_indices,
      batch2_indices = sim_params$batch2_indices,
      group1_indices = sim_params$group1_indices,
      group2_indices = sim_params$group2_indices
    )
  return(group_means)
}

filter_genes_based_on_sim_params_with_batches <- function(group_means, sim_params, mean_disp_loess) {
  group_means_filt <-
    group_means %>% dplyr::filter(
      b1_groups_fc < sim_params$max_baseline_log2FC_between_groups,
      b1_groups_fc > (-(
        sim_params$max_baseline_log2FC_between_groups
      )),
      b2_groups_fc < sim_params$max_baseline_log2FC_between_groups,
      b2_groups_fc > (-(
        sim_params$max_baseline_log2FC_between_groups
      )),
      abs(batch_fc) > sim_params$batch_difference_threshold,
      abs(batch_fc) < sim_params$batch_difference_threshold + 1
    )

  group_means_filt <- group_means_filt %>%
    dplyr::mutate(
      b1_dispersion = predict(mean_disp_loess, 2 ^ b1_mean),
      b2_dispersion = predict(mean_disp_loess, 2 ^ b2_mean)
    )

  group_means_filt$b1_dispersion[which(!is.finite(group_means_filt$b1_dispersion))] <-
    max(group_means_filt$b1_dispersion, na.rm = T)

  group_means_filt$b2_dispersion[which(!is.finite(group_means_filt$b2_dispersion))] <-
    max(group_means_filt$b2_dispersion, na.rm = T)
  return(group_means_filt)
}

filter_genes_based_on_sim_params_without_batches <- function(group_means, sim_params, mean_disp_loess) {
  group_means_filt <-
    group_means %>% dplyr::filter(
      group_fc < sim_params$max_baseline_log2FC_between_groups,
      group_fc > (-(
        sim_params$max_baseline_log2FC_between_groups
      ))
    )

  group_means_filt <- group_means_filt %>%
    dplyr::mutate(
      g1_dispersion = predict(mean_disp_loess, 2 ^ g1_mean),
      g2_dispersion = predict(mean_disp_loess, 2 ^ g2_mean)
    )

  group_means_filt$g1_dispersion[which(!is.finite(group_means_filt$g1_dispersion))] <-
    max(group_means_filt$g1_dispersion, na.rm = T)

  group_means_filt$g2_dispersion[which(!is.finite(group_means_filt$g2_dispersion))] <-
    max(group_means_filt$g2_dispersion, na.rm = T)
  return(group_means_filt)
}

filter_genes_based_on_sim_params <- function(simulated_datamat, sim_params,
                                             batch_effects_exist=FALSE){
  group_stats <- compute_groupwise_stats(simulated_datamat=simulated_datamat, sim_params)
  mean_disp_loess <- fit_loess(mean_and_disperison=mean_and_disperison)
  if(batch_effects_exist){
    group_means_filt <- filter_genes_based_on_sim_params_with_batches(group_means=group_stats, sim_params, mean_disp_loess=mean_disp_loess)
  } else {
    group_means_filt <- filter_genes_based_on_sim_params_without_batches(group_means=group_stats, sim_params, mean_disp_loess=mean_disp_loess)
  }
  if(nrow(group_means_filt) > sim_params$approx_n_genes) {
    gene_indices <- sort(sample(1:nrow(group_means_filt), sim_params$approx_n_genes))
    group_means_filt <- group_means_filt[gene_indices, ]
  }
  return(group_means_filt)
}
