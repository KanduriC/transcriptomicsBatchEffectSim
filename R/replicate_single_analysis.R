#' Repeats the single analysis (run_single_analysis) desired number of times for quantifying uncertainty
#'
#' @param sim_params_list A list object in R containing parameter specifications to run the analysis, where the legal parameters are described in the documentation
#' @param n_times An integer specifying the number of replications of the analysis
#'
#' @return A list containing results
#' @export
#'

rep_simulations <- function(sim_params_list, n_times){
  library(tidyr)
  results_list <- list()
  seed_list <- sample(1:1000, n_times)
  for (seed in seed_list) {
    set.seed(seed = seed)
    sim_params_list$seed <- seed
    results_list[[as.character(seed)]] <- run_single_analysis(sim_params=sim_params_list)
  }
  deseq_res <- extract_deseq_results(results_list)
  logistic_res <- extract_logistic_results(results_list)
  log_coef_overlap_stats <- extract_logistic_coef_overlap_metrics(results_list)
  return(list(deseq_results=deseq_res, logistic_results=logistic_res, log_coef_overlap_stats = log_coef_overlap_stats, detailed_results=results_list))
}

extract_deseq_results <- function(rep_sim_res_list){
  res_obj <- lapply(rep_sim_res_list, function(x){x$deseq_results})
  res_obj <- res_obj %>% dplyr::bind_rows() %>% dplyr::group_by(across(c(-values, -seed))) %>% dplyr::filter(grepl('sensitivity|specificity|accuracy', metric))
  return(res_obj)
}

extract_logistic_results <- function(rep_sim_res_list){
  res_obj <- lapply(rep_sim_res_list, function(x){x$logistic_reg_metrics})
  # res_obj <- res_obj %>% dplyr::bind_rows() %>% dplyr::group_by(across(c(-values, -seed))) %>% dplyr::filter(grepl('balanced_accuracy', metric))
  res_obj <- res_obj %>% dplyr::bind_rows()
  return(res_obj)
}

extract_logistic_coef_overlap_metrics <- function(rep_sim_res_list){
  res_obj <- lapply(rep_sim_res_list, function(x){x$log_coef_overlap_metrics})
  res_obj <- res_obj %>% dplyr::bind_rows()
  return(res_obj)
}

multi_param_rep_sim <- function(params_grid, sim_params, n_times) {
  grid_res <- list()
  for(each_row in 1:nrow(params_grid)){
    sim_params$batch_difference_threshold <- params_grid$Var1[each_row]
    sim_params$avg_log2FC_between_groups <- params_grid$Var2[each_row]
    grid_res[[each_row]] <- rep_simulations(sim_params_list = sim_params, n_times=n_times)
  }
  return(grid_res)
}
