determine_test_indices_with_batches <- function(sim_params){
  train_split_prop <- sim_params$train_split_prop
  n_examples_to_sample <- floor((sim_params$n_examples-sim_params$n_examples*train_split_prop)/4)
  b1_g1_indices_test <- sort(sample(sim_params$batch1_group1_indices, size=n_examples_to_sample, replace = F))
  b1_g2_indices_test <- sort(sample(sim_params$batch1_group2_indices, size=n_examples_to_sample, replace = F))
  b1_test_indices <- c(b1_g1_indices_test, b1_g2_indices_test)
  b2_g1_indices_test <- sort(sample(sim_params$batch2_group1_indices, size=n_examples_to_sample, replace = F))
  b2_g2_indices_test <- sort(sample(sim_params$batch2_group2_indices, size=n_examples_to_sample, replace = F))
  b2_test_indices <- c(b2_g1_indices_test, b2_g2_indices_test)
  test_indices <- c(b1_test_indices, b2_test_indices)
  return(test_indices)
}

determine_test_indices_without_batches <- function(sim_params){
  train_split_prop <- sim_params$train_split_prop
  n_examples_to_sample <- floor((sim_params$n_examples-sim_params$n_examples*train_split_prop)/2)
  g1_indices_test <- sort(sample(sim_params$group1_indices, size=n_examples_to_sample, replace = F))
  g2_indices_test <- sort(sample(sim_params$group2_indices, size=n_examples_to_sample, replace = F))
  test_indices <- c(g1_indices_test, g2_indices_test)
  return(test_indices)
}

determine_test_indices <- function(sim_params, batch_effects_exist=FALSE){
  if(batch_effects_exist){
    test_indices <- determine_test_indices_with_batches(sim_params)
  } else {
    test_indices <- determine_test_indices_without_batches(sim_params)
  }
  return(test_indices)
}

split_train_test <- function(sim_params, data_objs, batch_effects_exist=FALSE){
  test_index <- determine_test_indices(sim_params, batch_effects_exist=batch_effects_exist)
  data_objs$train_metadata <- data_objs$metadata[-test_index, ]
  data_objs$test_metadata <- data_objs$metadata[test_index, ]
  data_objs$train_data <- data_objs$exprs_mat[,rownames(data_objs$train_metadata)]
  data_objs$test_data <- data_objs$exprs_mat[,rownames(data_objs$test_metadata)]
  return(data_objs)
}

