## transcriptomicsBatchEffectSim

Scripts to simulate batch and biological effects in transcriptomics data to explore the problem terrain of handling batch effects.

### Installation

One common way to install R packages hosted on github is to use the popular `devtools` package.

```
devtools::install_github("https://github.com/KanduriC/transcriptomicsBatchEffectSim.git")
```

Alternatively, the package can be downloaded and installed in the following way:

```
devtools::install("path/to/transcriptomicsBatchEffectSim")
```

### Paramters specification

The `transcriptomicsBatchEffectSim` needs a list of parameters defined for the simulation of biological effects and/or batch effects. These parameters will define how many number of observations/samples need to be simulated, how many number of genes, the indices of two different batches, and the indices of biological groups within the batches, how much should the batches differ in log2 scale, what is the average fold-change in log2 scale that one expects for the biological effects, number of differentially expressed genes and the numbers that are up-regulated or down-regulated and so on. The code chunk below shows examples of parameter specifications for scenarios with and without batch effects.

```
## example parameters when including batch effects

sim_params_with_batch <-
  list(
    approx_n_genes = 5000, # number of genes
    n_examples = 60, # number of samples
    batch1_indices = 1:30, # the indices of batch-1
    batch2_indices = 31:60, # the indices of batch-2
    batch1_group1_indices = 1:24, # the indices of biological group-1 in batch-1 
    batch1_group2_indices = 25:30, # the indices of biological group-2 in batch-1
    batch2_group1_indices = 31:36, # the indices of biological group-1 in batch-2
    batch2_group2_indices = 37:60, # the indices of biological group-2 in batch-2
    group1_indices = c(1:24, 31:36), # the indices of biological group-1
    group2_indices = c(25:30, 37:60), # the indices of biological group-2
    batch_difference_threshold = 1, # the log2 fold-difference between batches +/- this parameter. Here, the gene expression between batches will vary between 0 and 2 (1-1, 1+1)
    n_true_diff_genes = 500, # number of true differentially expressed genes
    n_true_upreg_genes = 250, # number of true up-regulated genes
    max_baseline_log2FC_between_groups = 0.25, # how much the gene expression between groups will vary before introducing any biological effects
    avg_log2FC_between_groups = 1, average log2 fold-change of biological differences
    train_split_prop = 0.70, # when training a predictive model on the simulated data, what proportion of data will be included for training
    batch_effects_exist = TRUE # whether batch effects to exist in the data
  )

sim_params_without_batch <-
  list(
    approx_n_genes = 5000,
    n_examples = 60,
    group1_indices = 1:30,
    group2_indices = 31:60,
    n_true_diff_genes = 500,
    n_true_upreg_genes = 250,
    max_baseline_log2FC_between_groups = 0.25,
    avg_log2FC_between_groups = 1,
    train_split_prop = 0.70,
    batch_effects_exist = FALSE
  )
```

### A single analysis workflow

`transcriptomicsBatchEffectSim` does different operations depending upon whether batch effects are to be simulated or not. When batch effects are to be simulated, the following sequence of steps happen:

- simulates a baseline gene expression matrix, where read counts are simulated directly by sampling from negative bionomial distribution. The parameters for the negative binomial distribution are calibrated based on a real-world experimental RNAseq dataset to match the specific properties of RNAseq experiments in terms of the relation between mean and dispersion that varies depending on the expression level of genes.
- Induces batch differences to the desired extent
- Induces biological group differences to the desired extent
- corrects for batch effects using Combat-Seq implemented in sva package
- performs differential expression testing using DESeq2 before and after batch effect correction
- trains a lasso logistic regression model (where the regularization parameter is chosen through cross-validation) on a training split (as specified) and predicts on test data. This is performed both on data that is corrected and not corrected for batch effects and tested on data that contains batch effects

If the parameter specifications specify a simulated dataset without batch effects, the same steps above are performed except for introducing batch effects or their correction.

The following specification performs one single iteration of the analysis:

```
res_no_batch <- run_single_analysis(sim_params=sim_params_without_batch)
```

The results object `res_no_batch` contains all the results and all the intermediate results including simulation matrices. These can be written to disk as desired. 

### Multiple repetitions of the single analysis for uncertainty estimation

To quantify the uncertainty, the single analysis can be repeated multiple times, where a new independent dataset will be simulated each time and the same analyses as described above will be performed. The following code chunk shows how to accomplish this:

```
multi_res_no_batch <- rep_simulations(sim_params_list = sim_params_without_batch, n_times=10) #here n_time specifies the number of repetitions desired
```

## Correspondence

In case of questions/issues, contact Chakravarthi Kanduri (skanduri (at) uio (dot) no)
