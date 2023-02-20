set.seed(123)
sim_params <-
  list(
    approx_n_genes = 5000,
    n_examples = 60,
    batch1_indices = 1:30,
    batch2_indices = 31:60,
    batch1_group1_indices = 1:20,
    batch1_group2_indices = 21:30,
    batch2_group1_indices = 31:40,
    batch2_group2_indices = 41:60,
    group1_indices = c(1:20, 31:40),
    group2_indices = c(21:30, 41:60),
    batch_difference_threshold = 1,
    n_true_diff_genes = 500,
    n_true_upreg_genes = 250,
    max_baseline_log2FC_between_groups = 0.25,
    avg_log2FC_between_groups = 1,
    train_split_prop = 0.70
  )

data_objs <- simulate_exprsmat_with_batch_and_group_differences(sim_params = sim_params)
data_objs <- split_train_test(sim_params=sim_params, data_objs = data_objs, batch_effects_exist=TRUE)

test_that("determining test indices work", {
  expect_equal(nrow(data_objs$test_metadata), 16)
})
