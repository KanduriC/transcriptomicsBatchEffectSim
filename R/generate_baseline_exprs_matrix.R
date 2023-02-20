upsample_gene_parameters <- function(mean_and_disperison, approx_n_genes) {
  gene_indices <-
    sample(1:nrow(mean_and_disperison),
           approx_n_genes,
           replace = T)
  upsampled_mean_and_dispersion <- mean_and_disperison[gene_indices, ]
  mean_and_dispersion_list <-
    split(upsampled_mean_and_dispersion,
          1:nrow(upsampled_mean_and_dispersion))
  return(mean_and_dispersion_list)
}

sample_gene_counts <- function(gene, n_examples) {
  return(rnbinom(
    n_examples,
    mu = gene$gene_means,
    size = 1 / gene$dispersion
  ))
}

sample_n_gene_counts <- function(n_examples){
    sample_gene_counts <- function(gene) {
      return(rnbinom(
        n_examples,
        mu = gene$gene_means,
        size = 1 / gene$dispersion
      ))
    }
    return(sample_gene_counts)
}

generate_exprs_matrix <- function(mean_and_dispersion_list, n_examples) {
  my_sampling_func <- sample_n_gene_counts(n_examples)
  new_dat <- t(sapply(mean_and_dispersion_list, my_sampling_func))
  rownames(new_dat) <- paste0("gene_", 1:nrow(new_dat))
  return(new_dat)
}

generate_baseline_exprs_matrix <- function(n_genes, n_examples){
  mean_and_dispersion_list <- upsample_gene_parameters(mean_and_disperison = mean_and_disperison,
                                                       approx_n_genes = n_genes*12)
  baseline_mat <- generate_exprs_matrix(mean_and_dispersion_list=mean_and_dispersion_list, n_examples=n_examples)
  return(baseline_mat)
}
