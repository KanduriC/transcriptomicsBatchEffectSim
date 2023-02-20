induce_batch_differences <- function(baseline_exprs_mat, sim_params) {
  simulated_row_means <- rowMeans(baseline_exprs_mat)
  descending_order <- order(simulated_row_means, decreasing = TRUE)
  ordered_sim_matrix <- baseline_exprs_mat[descending_order, ]
  batch_1 <- ordered_sim_matrix[, sim_params$batch1_indices]
  batch_2 <- ordered_sim_matrix[, sim_params$batch2_indices]
  rev_batch_1 <- apply(batch_1, 2, rev)
  batch_simulated_datamat <- data.frame(rev_batch_1, batch_2)
  return(batch_simulated_datamat)
}