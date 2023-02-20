correct_batch_effects <- function(batch_simulated_datamat, batch, group){
  batch_adjusted <-
    sva::ComBat_seq(as.matrix(batch_simulated_datamat),
               batch = batch,
               group = group)
  return(batch_adjusted)
}
