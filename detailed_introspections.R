```{r, message=FALSE, cache=FALSE, echo=TRUE,warning=FALSE}
#return sum of nonzero coefficients other than intercept
# look how true signal distribution is in uncorrected and corrected data in train and test data
# return predicted probabilities and classes

signif_features <- rownames(res_1_80$detailed_results$`443`$simulated_data_objs$group_means_filt)[res_1_80$detailed_results$`443`$simulated_data_objs$diff_exp_indices]

train_dat <- as.data.frame(scale(t(log2(res_1_80$detailed_results$`443`$simulated_data_objs$train_data+0.01))))
train_dat_corr <- as.data.frame(scale(t(log2(res_1_80$detailed_results$`443`$simulated_data_objs$train_data_corrected+0.01))))
test_dat <- as.data.frame(scale(t(log2(res_1_80$detailed_results$`443`$simulated_data_objs$test_data+0.01))))
train_metadat <- res_1_80$detailed_results$`443`$simulated_data_objs$train_metadata
train_metadat <- tibble::rownames_to_column(train_metadat, "sample_name")

test_metadat <- res_1_80$detailed_results$`443`$simulated_data_objs$test_metadata
test_metadat <- tibble::rownames_to_column(test_metadat, "sample_name")

select_gene <- signif_features[1]
select_gene <- "gene_14884"

train_gene_dat <- train_dat %>% dplyr::select(select_gene) %>% dplyr::mutate(type="train")
train_gene_dat <- tibble::rownames_to_column(train_gene_dat, "sample_name")
train_gene_dat <- train_gene_dat %>% dplyr::left_join(train_metadat)

train_corr_gene_dat <- train_dat_corr %>% dplyr::select(select_gene) %>% dplyr::mutate(type="train_corr")
train_corr_gene_dat <- tibble::rownames_to_column(train_corr_gene_dat, "sample_name")
train_corr_gene_dat <- train_corr_gene_dat %>% dplyr::left_join(train_metadat)

test_gene_dat <- test_dat %>% dplyr::select(select_gene) %>% dplyr::mutate(type="test")
test_gene_dat <- tibble::rownames_to_column(test_gene_dat, "sample_name")
test_gene_dat <- test_gene_dat %>% dplyr::left_join(test_metadat)

gene_dat <- dplyr::bind_rows(train_gene_dat, test_gene_dat)
colnames(gene_dat)[2] <- "value"

gene_dat <- dplyr::bind_rows(train_corr_gene_dat, test_gene_dat)
colnames(gene_dat)[2] <- "value"

p18 <- ggplot(gene_dat, aes(x=group, y=value, color=batch))+geom_boxplot()+geom_point(position = "dodge")+facet_wrap(~type)
ggsave("~/Documents/Projects/causalAIRR/true_signif_example2.pdf", plot = p18)
```