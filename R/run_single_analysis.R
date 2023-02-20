#' Run a single iteration of analysis with the parameter specifications
#'
#' @param sim_params A list object in R containing parameter specifications to run the analysis, where the legal parameters are described in the documentation
#'
#' @return A list containing results
#' @export
#'
run_single_analysis <- function(sim_params){
  library(tidyr)
  if(sim_params$batch_effects_exist){
    results <- single_analysis_with_batch(sim_params=sim_params)
  } else {
    results <- single_analysis_no_batch(sim_params=sim_params)
  }
  return(results)
}

single_analysis_with_batch <- function(sim_params){
  results_list <- list()
  data_objs_batch <- simulate_exprsmat_with_batch_and_group_differences(sim_params = sim_params)
  data_objs_batch <- split_train_test(sim_params = sim_params, data_objs = data_objs_batch, batch_effects_exist=TRUE)
  data_objs_batch$train_data_corrected <- correct_batch_effects(
    batch_simulated_datamat = data_objs_batch$train_data,
    batch = data_objs_batch$train_metadata$batch,
    group = data_objs_batch$train_metadata$group)
  logistic_res_corrected <- logistic_regression(x_train=data_objs_batch$train_data_corrected, y_train = data_objs_batch$train_metadata$group, x_test = data_objs_batch$test_data, y_test = data_objs_batch$test_metadata$group, perform_lasso = TRUE)
  results_list <- deseq_feature_selection(data_objs_batch$train_data_corrected, data_objs_batch$train_metadata, data_objs_batch$group_means_filt, data_objs_batch$diff_exp_indices, sim_params, results_list)
  uncorrected_deseq_res <- deseq_feature_selection(data_objs_batch$train_data, data_objs_batch$train_metadata, data_objs_batch$group_means_filt, data_objs_batch$diff_exp_indices, sim_params, results_list=list())
  uncorrected_deseq_res <- uncorrected_deseq_res$metrics %>% dplyr::mutate(type="uncorrected")
  results_list$metrics <- results_list$metrics %>% dplyr::mutate(type="corrected")
  data_objs_batch$train_data_corrected_signif <- results_list$signif_feature_mat
  signif_features <- rownames(data_objs_batch$train_data_corrected_signif)
  data_objs_batch$test_data_signif <- data_objs_batch$test_data[signif_features,]
  logistic_res_corrected_signif <- logistic_regression(x_train=data_objs_batch$train_data_corrected_signif, y_train = data_objs_batch$train_metadata$group, x_test = data_objs_batch$test_data_signif, y_test = data_objs_batch$test_metadata$group, perform_lasso = TRUE)
  logistic_res_corrected_signif_no_reg <- logistic_regression(x_train=data_objs_batch$train_data_corrected_signif, y_train = data_objs_batch$train_metadata$group, x_test = data_objs_batch$test_data_signif, y_test = data_objs_batch$test_metadata$group, perform_lasso = FALSE)
  logistic_res_uncorrected <- logistic_regression(x_train=data_objs_batch$train_data, y_train = data_objs_batch$train_metadata$group, x_test = data_objs_batch$test_data, y_test = data_objs_batch$test_metadata$group, perform_lasso = TRUE)
  res <- list(uncorrected=logistic_res_uncorrected, corrected=logistic_res_corrected, corrected_signif=logistic_res_corrected_signif, corrected_signif_no_reg= logistic_res_corrected_signif_no_reg)
  res <- extract_metric_coef_lists(res)
  logistic_reg_metrics <- do.call(rbind, Map(data.frame, res$metric_lists))
  logistic_reg_metrics <- reshape2::melt(t(logistic_reg_metrics))
  colnames(logistic_reg_metrics) <- c("metric", "type", "values")
  batch_group_corr <- round(length(sim_params$batch1_group1_indices)/length(sim_params$batch1_indices) * 100)
  logistic_reg_metrics <- logistic_reg_metrics %>% dplyr::mutate(group_diff = sim_params$avg_log2FC_between_groups, batch_diff = sim_params$batch_difference_threshold, batch_group_corr = batch_group_corr, seed=sim_params$seed)
  deseq_metrics <- dplyr::bind_rows(results_list$metrics, uncorrected_deseq_res) %>% dplyr::mutate(group_diff = sim_params$avg_log2FC_between_groups, batch_diff = sim_params$batch_difference_threshold, batch_group_corr = batch_group_corr, seed=sim_params$seed)
  results <- list(deseq_results = deseq_metrics, logistic_reg_metrics = logistic_reg_metrics, logistic_coefs = res$coef_lists,simulated_data_objs=data_objs_batch, sim_params = sim_params)
  return(results)
}


 single_analysis_no_batch <- function(sim_params){
   results_list <- list()
   data_objs <- simulate_exprsmat_with_group_differences(sim_params = sim_params)
   data_objs <- split_train_test(sim_params = sim_params, data_objs = data_objs, batch_effects_exist=FALSE)
   results_list <- deseq_feature_selection(data_objs$train_data, data_objs$train_metadata, data_objs$group_means_filt, data_objs$diff_exp_indices, sim_params, results_list)
   data_objs$train_data_signif <- results_list$signif_feature_mat
   signif_features <- rownames(data_objs$train_data_signif)
   data_objs$test_data_signif <- data_objs$test_data[signif_features,]
   logistic_res_signif <- logistic_regression(x_train=data_objs$train_data_signif, y_train = data_objs$train_metadata$group, x_test = data_objs$test_data_signif, y_test = data_objs$test_metadata$group, perform_lasso = TRUE)
   logistic_res_signif_no_reg <- logistic_regression(x_train=data_objs$train_data_signif, y_train = data_objs$train_metadata$group, x_test = data_objs$test_data_signif, y_test = data_objs$test_metadata$group, perform_lasso = FALSE)
   logistic_res_uncorrected <- logistic_regression(x_train=data_objs$train_data, y_train = data_objs$train_metadata$group, x_test = data_objs$test_data, y_test = data_objs$test_metadata$group, perform_lasso = TRUE)
   res <- list(no_batch_effects=logistic_res_uncorrected, no_batch_signif=logistic_res_signif, no_batch_signif_no_reg= logistic_res_signif_no_reg)
   res <- extract_metric_coef_lists(res)
   logistic_reg_metrics <- do.call(rbind, Map(data.frame, res$metric_lists))
   logistic_reg_metrics <- reshape2::melt(t(logistic_reg_metrics))
   colnames(logistic_reg_metrics) <- c("metric", "type", "values")
   logistic_reg_metrics <- logistic_reg_metrics %>% dplyr::mutate(group_diff = sim_params$avg_log2FC_between_groups, seed=sim_params$seed)
   results_list$metrics <- results_list$metrics %>% dplyr::mutate(group_diff = sim_params$avg_log2FC_between_groups, seed=sim_params$seed, type = "no batch effect")
   results <- list(deseq_results = results_list$metrics, logistic_reg_metrics = logistic_reg_metrics, logistic_coefs = res$coef_lists, simulated_data_objs=data_objs, sim_params = sim_params)
   return(results)
 }

 extract_metric_coef_lists <- function(list_of_logistic_res_lists){
   metric_lists <- lapply(list_of_logistic_res_lists, function(l) {
     l[!names(l) %in% "logistic_coef"]
   })
   coef_lists <- lapply(list_of_logistic_res_lists, function(l) {
     l[names(l) %in% "logistic_coef"]
   })
   res <- list(metric_lists = metric_lists, coef_lists = coef_lists)
   return(res)
 }
