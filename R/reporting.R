#' An ad hoc utility function to write the results of analyses to disk. Writes nested lists using folder structures.
#'
#' @param mylist a list of analyses results
#' @param filepath path to write the results
#'
#' @return writes r objects to disk, but does not return anything
#' @export
#'
write_list_to_disk <- function(mylist, filepath) {
  if (is.data.frame(mylist)) {
    write.table(mylist, file = paste0(filepath, ".txt"), sep = "\t", quote=F, row.names = T)
  } else if (!is.list(mylist)) {
    write.table(mylist, file = paste0(filepath, ".txt"), sep="\t", quote=F, row.names = T)
    }
  else {
    if (!dir.exists(filepath)) {
      dir.create(filepath, showWarnings = FALSE)
    }
    if ("sim_params" %in% names(mylist)) {
      yaml::write_yaml(mylist[["sim_params"]], file.path(filepath, "sim_params.yaml"))
    }
    for (i in names(mylist)[!names(mylist) %in% "sim_params"]) {
      subpath <- file.path(filepath, i)
      write_list_to_disk(mylist[[i]], subpath)
    }
  }
}

#' An ad hoc utility function to combine the results of univariate statistical testing between two different parameter settings to make it easy for plotting
#'
#' @param results_list_1 A list of results obtained through `rep_simulations` for one particular parameter setting
#' @param results_list_2 A list of results obtained through `rep_simulations` for one particular parameter setting
#' @param table_name A string among "deseq_results" or "logistic_results"
#'
#' @return A dataframe containing results binded together
#' @export
#'
join_results <- function(results_list_1, results_list_2, table_name){
  if(table_name=="log_coef_overlap_stats"){
    res <- results_list_1[[table_name]] %>% dplyr::bind_rows(results_list_2[[table_name]]) %>% dplyr::filter(id %in% c("no_batch_effects", "uncorrected", "corrected")) %>% dplyr::mutate(metric = "true signal in non-zero coefs proportion")
  } else {
    res_1 <- results_list_1[[table_name]] %>% dplyr::ungroup() %>% dplyr::select(metric, values, type)
    res_2 <- results_list_2[[table_name]] %>% dplyr::ungroup() %>% dplyr::select(metric, values, type)
    res <- res_1 %>% dplyr::bind_rows(res_2)
  }
  if(table_name=="logistic_results"){
    # res <- res %>% dplyr::filter(type %in% c("no_batch_effects", "uncorrected", "corrected")) %>% dplyr::mutate(metric="Balanced accuracy")
    res <- res %>% dplyr::filter(type %in% c("no_batch_effects", "uncorrected", "corrected"))
  }
  return(res)
}

#' An ad hoc box plot function to visualize results obtained through `rep_simulations` and joined through `join_results`
#'
#' @param results_table a dataframe containing results
#' @param title a string representing title for the plot
#'
#' @return a ggplot
#' @export
#'
box_plot <- function(results_table, title, y_label){
  ggplot(results_table, aes(x=type, y=values)) + geom_boxplot(outlier.shape = NA) + geom_point(width = 0.25, size=0.75)+theme_bw() + facet_wrap(~metric, scales = "free") + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))+ theme(plot.margin = unit(c(2,1,1,1), "cm"))+xlab("Type of dataset") + ylab(y_label)+theme(axis.text.x = element_text(face="bold", size=12, angle = 45, vjust = 0.5, hjust=0.5),axis.text.y = element_text(face="bold", size=12),plot.title = element_text(size=12,face="bold",hjust = 0.5),axis.title.x = element_text(size=12, face="bold"),axis.title.y = element_text(size=12, face="bold"),legend.title = element_text(size=12, face="bold"),legend.text = element_text(size=12, face="bold"))+theme(strip.background =element_rect(colour="black",fill="white"),strip.text.x = element_text(size = 12, colour = "darkblue", face="bold"))+theme(plot.margin=unit(c(1,1,1,1),"cm")) + ggtitle(title)
}
