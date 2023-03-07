deseq_feature_selection <- function(simulated_datamat, coldata, group_means_filt, diff_exp_indices, sim_params, results_list) {
  dds <-
    DESeq2::DESeqDataSetFromMatrix(countData = simulated_datamat,
                           colData = coldata,
                           design = ~ group)
  dds <- DESeq2::DESeq(dds)
  res <-
    DESeq2::results(
      dds,
      contrast = c("group", "group_1", "group_2"),
      independentFiltering = F
    )
  resLFC <-
    DESeq2::lfcShrink(
      dds,
      contrast = c("group", "group_1", "group_2"),
      res = res,
      type = "ashr"
    )

  res.sig_without_batch <- data.frame(resLFC)
  res.sig_without_batch <-
    res.sig_without_batch[res.sig_without_batch$padj < 0.10,]
  res.sig_without_batch <- res.sig_without_batch[complete.cases(res.sig_without_batch), ]
  n_sig_without_batch <- nrow(res.sig_without_batch)
  n_true_sig_without_batch <- 0
  if (n_sig_without_batch > 0) {
    n_true_sig_without_batch <-
      length(intersect(
        rownames(res.sig_without_batch),
        rownames(group_means_filt)[diff_exp_indices]
      ))
  }
  n_false_pos_without_batch <-
    n_sig_without_batch - n_true_sig_without_batch
  specificity_without_batch <-
    (1 - (
      n_false_pos_without_batch / (nrow(resLFC) - sim_params$n_true_diff_genes)
    )) * 100
  true_neg <- nrow(resLFC) - sim_params$n_true_diff_genes - n_false_pos_without_batch
  accuracy <- ((n_true_sig_without_batch + true_neg)/nrow(resLFC))*100
  precision <- (n_true_sig_without_batch/n_sig_without_batch) * 100
  sensitivity_without_batch <-
    (n_true_sig_without_batch / sim_params$n_true_diff_genes) * 100
  f1_score <- 2 * (precision * sensitivity_without_batch) / (precision + sensitivity_without_batch)
  metric <- c("n_significant","n_true_significant","sensitivity", "specificity", "accuracy", "f1_score")
  values <- c(n_sig_without_batch, n_true_sig_without_batch, sensitivity_without_batch, specificity_without_batch, accuracy, f1_score)
  metrics <- data.frame(metric, values)
  results_list$metrics <- metrics
  results_list$signif_feature_mat <- simulated_datamat[rownames(res.sig_without_batch),]
  return(results_list)
}
