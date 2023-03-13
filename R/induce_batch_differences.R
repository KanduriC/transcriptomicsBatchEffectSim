# induce_batch_differences <- function(baseline_exprs_mat, sim_params) {
#   simulated_row_means <- rowMeans(baseline_exprs_mat)
#   descending_order <- order(simulated_row_means, decreasing = TRUE)
#   ordered_sim_matrix <- baseline_exprs_mat[descending_order, ]
#   batch_1 <- ordered_sim_matrix[, sim_params$batch1_indices]
#   batch_2 <- ordered_sim_matrix[, sim_params$batch2_indices]
#   rev_batch_1 <- apply(batch_1, 2, rev)
#   batch_simulated_datamat <- data.frame(rev_batch_1, batch_2)
#   return(batch_simulated_datamat)
# }

induce_batch_differences <- function(baseline_exprs_mat, sim_params) {
  b1_up_indices <- 1:(nrow(baseline_exprs_mat)/2)
  b2_up_indices <- ((nrow(baseline_exprs_mat)/2)+1):nrow(baseline_exprs_mat)
  mean_disp_loess <- fit_loess(mean_and_disperison=mean_and_disperison)
  group_stats <- compute_groupwise_stats(simulated_datamat=baseline_exprs_mat, sim_params) %>% dplyr::select(b1_mean, b2_mean)
  group_stats[b1_up_indices, ] <- group_stats[b1_up_indices, ] %>% dplyr::mutate(b1_mean = b1_mean + sim_params$batch_difference_threshold + 0.25)
  group_stats[b2_up_indices, ] <- group_stats[b2_up_indices, ] %>% dplyr::mutate(b2_mean = b2_mean + sim_params$batch_difference_threshold + 0.25)
  group_stats <- group_stats %>%
    dplyr::mutate(
      b1_dispersion = predict(mean_disp_loess, 2 ^ b1_mean),
      b2_dispersion = predict(mean_disp_loess, 2 ^ b2_mean)
    )
  group_stats$b1_dispersion[which(!is.finite(group_stats$b1_dispersion))] <-
    mean(group_stats$b1_dispersion, na.rm = T)
  group_stats$b2_dispersion[which(!is.finite(group_stats$b2_dispersion))] <-
    mean(group_stats$b2_dispersion, na.rm = T)
  b1_mean_dispersion <-
    group_stats[b1_up_indices, ] %>% dplyr::select(b1_mean, b1_dispersion)
  b2_mean_dispersion <-
    group_stats[b2_up_indices, ] %>% dplyr::select(b2_mean, b2_dispersion)
  b1_mean_dispersion_list <-
    split(b1_mean_dispersion, 1:nrow(b1_mean_dispersion))
  b2_mean_dispersion_list <-
    split(b2_mean_dispersion, 1:nrow(b2_mean_dispersion))

  b1_upreg_counts  <-
    t(sapply(b1_mean_dispersion_list, function(gene) {
      rnbinom(
        length(sim_params$batch1_indices),
        mu = 2 ^ (gene$b1_mean),
        size = 1 / gene$b1_dispersion
      )
    }))
  b2_upreg_counts  <-
    t(sapply(b2_mean_dispersion_list, function(gene) {
      rnbinom(
        length(sim_params$batch2_indices),
        mu = 2 ^ (gene$b2_mean),
        size = 1 / gene$b2_dispersion
      )
    }))

  baseline_exprs_mat[b1_up_indices, sim_params$batch1_indices] <- b1_upreg_counts
  baseline_exprs_mat[b2_up_indices, sim_params$batch2_indices] <- b2_upreg_counts
  return(baseline_exprs_mat)
}

