induce_group_differences_without_batches <- function(simulated_datamat, group_means_filt, sim_params) {
  diff_exp_indices <-
    sample(1:nrow(simulated_datamat), sim_params$n_true_diff_genes, replace = F)
  g1_up_indices <- diff_exp_indices[1:sim_params$n_true_upreg_genes]
  g2_up_indices <-
    diff_exp_indices[(sim_params$n_true_upreg_genes + 1):sim_params$n_true_diff_genes]
  fold_change <- sim_params$avg_log2FC_between_groups + 0.25

  group_means_filt[g1_up_indices, ] <- group_means_filt[g1_up_indices, ] %>% dplyr::mutate(g1_mean = g1_mean + fold_change)
  group_means_filt[g2_up_indices, ] <- group_means_filt[g2_up_indices, ] %>% dplyr::mutate(g2_mean = g2_mean + fold_change)

  mean_disp_loess <- fit_loess(mean_and_disperison=mean_and_disperison)
  group_means_filt <- group_means_filt %>%
    dplyr::mutate(
      g1_dispersion = predict(mean_disp_loess, 2 ^ g1_mean),
      g2_dispersion = predict(mean_disp_loess, 2 ^ g2_mean)
    )
  group_means_filt$g1_dispersion[which(!is.finite(group_means_filt$g1_dispersion))] <-
    mean(group_means_filt$g1_dispersion, na.rm = T)
  group_means_filt$g2_dispersion[which(!is.finite(group_means_filt$g2_dispersion))] <-
    mean(group_means_filt$g2_dispersion, na.rm = T)

  group1_mean_dispersion <-
    group_means_filt[g1_up_indices, ] %>% dplyr::select(g1_mean, g1_dispersion)
  group2_mean_dispersion <-
    group_means_filt[g2_up_indices, ] %>% dplyr::select(g2_mean, g2_dispersion)
  group1_mean_dispersion_list <-
    split(group1_mean_dispersion, 1:nrow(group1_mean_dispersion))
  group2_mean_dispersion_list <-
    split(group2_mean_dispersion, 1:nrow(group2_mean_dispersion))

  g1_upreg_counts  <-
    t(sapply(group1_mean_dispersion_list, function(gene) {
      rnbinom(
        length(sim_params$group1_indices),
        mu = 2 ^ (gene$g1_mean),
        size = 1 / gene$g1_dispersion
      )
    }))
  g2_upreg_counts  <-
    t(sapply(group2_mean_dispersion_list, function(gene) {
      rnbinom(
        length(sim_params$group2_indices),
        mu = 2 ^ (gene$g2_mean),
        size = 1 / gene$g2_dispersion
      )
    }))

  simulated_datamat[g1_up_indices, sim_params$group1_indices] <- g1_upreg_counts

  simulated_datamat[g2_up_indices, sim_params$group2_indices] <- g2_upreg_counts

  return(list(simulated_datamat_filt = simulated_datamat,
              g1_up_indices = g1_up_indices, g2_up_indices = g2_up_indices, diff_exp_indices = diff_exp_indices))
}


induce_group_differences_with_batches <- function(batch_simulated_datamat, group_means_filt, sim_params) {
  diff_exp_indices <-
    sample(1:nrow(batch_simulated_datamat),
           sim_params$n_true_diff_genes,
           replace = F)
  g1_up_indices <- diff_exp_indices[1:sim_params$n_true_upreg_genes]
  g2_up_indices <-
    diff_exp_indices[(sim_params$n_true_upreg_genes + 1):sim_params$n_true_diff_genes]
  fold_change <- sim_params$avg_log2FC_between_groups + 0.25

  group_means_filt[g1_up_indices, ] <- group_means_filt[g1_up_indices, ] %>% dplyr::mutate(b1_g1_mean = b1_g1_mean + fold_change, b2_g1_mean = b2_g1_mean + fold_change)
  group_means_filt[g2_up_indices, ] <- group_means_filt[g2_up_indices, ] %>% dplyr::mutate(b1_g2_mean = b1_g2_mean + fold_change, b2_g2_mean = b2_g2_mean + fold_change)

  mean_disp_loess <- fit_loess(mean_and_disperison=mean_and_disperison)
  group_means_filt <- group_means_filt %>%
    dplyr::mutate(
      b1_g1_dispersion = predict(mean_disp_loess, 2 ^ b1_g1_mean),
      b2_g1_dispersion = predict(mean_disp_loess, 2 ^ b2_g1_mean),
      b1_g2_dispersion = predict(mean_disp_loess, 2 ^ b1_g2_mean),
      b2_g2_dispersion = predict(mean_disp_loess, 2 ^ b2_g2_mean)
    )
  group_means_filt$b1_g1_dispersion[which(!is.finite(group_means_filt$b1_g1_dispersion))] <-
    mean(group_means_filt$b1_g1_dispersion, na.rm = T)
  group_means_filt$b2_g1_dispersion[which(!is.finite(group_means_filt$b2_g1_dispersion))] <-
    mean(group_means_filt$b2_g1_dispersion, na.rm = T)
  group_means_filt$b1_g2_dispersion[which(!is.finite(group_means_filt$b1_g2_dispersion))] <-
    mean(group_means_filt$b1_g2_dispersion, na.rm = T)
  group_means_filt$b2_g2_dispersion[which(!is.finite(group_means_filt$b2_g2_dispersion))] <-
    mean(group_means_filt$b2_g2_dispersion, na.rm = T)

  group1_mean_dispersion <-
    group_means_filt[g1_up_indices, ] %>% dplyr::select(b1_g1_mean, b2_g1_mean, b1_g1_dispersion, b2_g1_dispersion)
  group2_mean_dispersion <-
    group_means_filt[g2_up_indices, ] %>% dplyr::select(b1_g2_mean, b2_g2_mean, b1_g2_dispersion, b2_g2_dispersion)
  group1_mean_dispersion_list <-
    split(group1_mean_dispersion, 1:nrow(group1_mean_dispersion))
  group2_mean_dispersion_list <-
    split(group2_mean_dispersion, 1:nrow(group2_mean_dispersion))

  b1_g1_upreg_counts  <-
    t(sapply(group1_mean_dispersion_list, function(gene) {
      rnbinom(
        length(sim_params$batch1_group1_indices),
        mu = 2 ^ (gene$b1_g1_mean),
        size = 1 / gene$b1_g1_dispersion
      )
    }))
  b2_g1_upreg_counts  <-
    t(sapply(group1_mean_dispersion_list, function(gene) {
      rnbinom(
        length(sim_params$batch2_group1_indices),
        mu = 2 ^ (gene$b2_g1_mean),
        size = 1 / gene$b2_g1_dispersion
      )
    }))
  b1_g2_upreg_counts  <-
    t(sapply(group2_mean_dispersion_list, function(gene) {
      rnbinom(
        length(sim_params$batch1_group2_indices),
        mu = 2 ^ (gene$b1_g2_mean),
        size = 1 / gene$b1_g2_dispersion
      )
    }))
  b2_g2_upreg_counts  <-
    t(sapply(group2_mean_dispersion_list, function(gene) {
      rnbinom(
        length(sim_params$batch2_group2_indices),
        mu = 2 ^ (gene$b2_g2_mean),
        size = 1 / gene$b2_g2_dispersion
      )
    }))

  batch_simulated_datamat[g1_up_indices, sim_params$batch1_group1_indices] <- b1_g1_upreg_counts
  batch_simulated_datamat[g1_up_indices, sim_params$batch2_group1_indices] <- b2_g1_upreg_counts
  batch_simulated_datamat[g2_up_indices, sim_params$batch1_group2_indices] <- b1_g2_upreg_counts
  batch_simulated_datamat[g2_up_indices, sim_params$batch2_group2_indices] <- b2_g2_upreg_counts

  return(list(batch_simulated_datamat_filt = batch_simulated_datamat,
              g1_up_indices = g1_up_indices, g2_up_indices = g2_up_indices, diff_exp_indices = diff_exp_indices))
}

induce_group_differences <- function(simulated_datamat, group_means_filt,
                                     sim_params, batch_effects_exist=FALSE){
  if(batch_effects_exist){
    sim_res <- induce_group_differences_with_batches(batch_simulated_datamat=simulated_datamat, group_means_filt, sim_params)
  } else {
    sim_res <- induce_group_differences_without_batches(simulated_datamat=simulated_datamat, group_means_filt, sim_params)
  }
  return(sim_res)
}
