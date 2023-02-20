prepare_metadata_and_exprsmat_with_batches <- function(sim_params, batch_simulated_datamat_filt) {
  sample_names <-
    c(
      paste0("b1_g1_s", 1:length(sim_params$batch1_group1_indices)),
      paste0("b1_g2_s", 1:length(sim_params$batch1_group2_indices)),
      paste0("b2_g1_s", 1:length(sim_params$batch2_group1_indices)),
      paste0("b2_g2_s", 1:length(sim_params$batch2_group2_indices))
    )

  batch <-
    c(rep("batch_1", length(sim_params$batch1_indices)), rep("batch_2", length(sim_params$batch2_indices)))
  group <-
    c(
      rep("group_1", length(sim_params$batch1_group1_indices)),
      rep("group_2", length(sim_params$batch1_group2_indices)),
      rep("group_1", length(sim_params$batch2_group1_indices)),
      rep("group_2", length(sim_params$batch2_group2_indices))
    )
  coldata <- data.frame(batch, group, row.names = sample_names)
  colnames(batch_simulated_datamat_filt) <- sample_names
  return(list(exprs_mat=batch_simulated_datamat_filt, metadata=coldata, batch=batch, group=group))
}

prepare_metadata_and_exprsmat_without_batches <- function(sim_params, simulated_datamat_filt) {
  sample_names <-
    c(
      paste0("g1_s", 1:length(sim_params$group1_indices)),
      paste0("g2_s", 1:length(sim_params$group2_indices))
    )

  batch <- c(rep("batch_1", sim_params$n_examples))
  group <-
    c(
      rep("group_1", length(sim_params$group1_indices)),
      rep("group_2", length(sim_params$group2_indices))
    )
  coldata <- data.frame(batch, group, row.names = sample_names)
  colnames(simulated_datamat_filt) <- sample_names
  return(list(exprs_mat=simulated_datamat_filt, metadata=coldata, batch=batch, group=group))
}

prepare_metadata_and_exprsmat <- function(sim_params, simulated_datamat,
                                          batch_effects_exist=FALSE){
  if(batch_effects_exist){
    data_objs_list <- prepare_metadata_and_exprsmat_with_batches(sim_params, batch_simulated_datamat_filt=simulated_datamat)
  } else {
    data_objs_list <- prepare_metadata_and_exprsmat_without_batches(sim_params, simulated_datamat_filt=simulated_datamat)
  }

  return(data_objs_list)
}
