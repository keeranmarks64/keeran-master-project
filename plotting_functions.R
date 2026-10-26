#PLotting Functions for results

library(rhdf5)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(reshape2)
library(matrixStats)
library(stringr)
library(viridis)
library(dplyr)
library(tidyr)

plot_3metrics <- function(sim_file,results_file,test){
  
  ##Reading in data and parameters
  
  results=h5read(results_file,name='results')
  sim_params = h5read(sim_file,"simulation_parameters")
  temp = h5ls(results_file,recursive=2)
  sim_name = substr(sim_file,start=10,stop=nchar(sim_file)-3)
  
  #extracting all group names
  groups = setdiff(temp$name, c('original_data', 'simulation_parameters','results')) 
  
  #extracting abundance and prevalence scaling values
  ab.scale = sim_params[['sim_params']][['ab.scale']]
  prev.scale = sim_params[['sim_params']][['prev.scale']]
  
  #extracting group, subset, repeat, and rep values
  no.groups = length(groups)
  no.subsets = length(sim_params[['test_idx']][['subsets']])
  subsets = names(results[['ab1_prev1_rep1']])
  subset_sizes = sort(sim_params[['test_idx']][['subsets']])
  subsets = paste0("subset_",as.character(subset_sizes))
  no.repeats = sim_params[['sim_params']][['repeats']]
  no.reps = length(sim_params[['test_idx']][['names']])-1 
  
  ##Extracting and Aggregating Data
  
  #Matrices to store data
  auroc_matrix = data.frame(row.names=as.character(seq(1,no.subsets,1)))
  recall_matrix = data.frame(row.names=as.character(seq(1,no.subsets,1)))
  fdr_matrix = data.frame(row.names=as.character(seq(1,no.subsets,1)))
  
  #iterating through groups to collect AUROC, Recall, and FDR metrics
  for(group in groups){
    group_aurocs = c()
    group_recalls = c()
    group_fdrs = c()
    
    #Iterating through each subset of group
    for(subset in subsets){
      
      #Aggregating metrics across repetitions of each subset
      subset_data = results[[group]][[subset]]
      group_aurocs=c(group_aurocs,mean(subset_data[['auroc']]))
      
      subset_recalls = subset_data[['TP']]/(subset_data[['TP']] + subset_data[['FN']])
      group_recalls = c(group_recalls,mean(subset_recalls))
      
      subset_fdrs = as.numeric(subset_data[['FP']]/(subset_data[['TP']] + subset_data[['FP']]))
      subset_fdrs = replace(subset_fdrs,which(is.nan(subset_fdrs)),0)
      group_fdrs = c(group_fdrs,mean(subset_fdrs))
      
    }
    
    auroc_matrix[[group]] = group_aurocs
    recall_matrix[[group]] = group_recalls
    fdr_matrix[[group]] = group_fdrs
    
  }
  
  #Aggregating metrics across group reps
  auroc_matrix_agg = data.frame(row.names=subsets)
  recall_matrix_agg = data.frame(row.names=subsets)
  fdr_matrix_agg = data.frame(row.names=subsets)
  groups_agg = c()
  
  #Effect combinations
  ab.scales = rep(ab.scale,rep(length(prev.scale),length(ab.scale)))
  prev.scales = rep(prev.scale,length(ab.scale))
  effect_combs=data.frame(ab.scales=ab.scales,prev.scale=prev.scales)
  
  #Iterating across group repeats to aggregate
  for(i in 1:(length(groups)/no.repeats)){
    
    #Determining name of group
    starting_idx = 1+no.repeats*(i-1)
    # name = colnames(auroc_matrix)[starting_idx]
    # name = substr(name,1,nchar(name)-5)
    name = paste0("Ab",effect_combs[i,1],"Prev",effect_combs[i,2])
    
    #Aggregating repeats belonging to particular group
    auroc_matrix_agg[[name]] = rowMeans(auroc_matrix[,starting_idx:(starting_idx+no.repeats-1)])
    recall_matrix_agg[[name]] = rowMeans(recall_matrix[,starting_idx:(starting_idx+no.repeats-1)])
    fdr_matrix_agg[[name]] = rowMeans(fdr_matrix[,starting_idx:(starting_idx+no.repeats-1)])
    groups_agg[i]=name
    
  }
  
  rownames(effect_combs) = groups_agg
  
  metrics = list(auroc=auroc_matrix_agg,
               fdr= fdr_matrix_agg,
               recall = recall_matrix_agg)
  
  ##PLotting Aggregated results for each group 
  for(group_name in groups_agg){
    
    fdr = ggplot(NULL,aes(x=subset_sizes,1-fdr_matrix_agg[[group_name]])) +
      geom_line() + 
      geom_point() +
      ylim(0,1) +
      ggtitle("1-FDR") +
      geom_hline(yintercept=0.95,linetype='dashed',colour='red') +
      scale_x_log10(breaks=subset_sizes) +
      xlab("Sample Size") +
      ylab('1 - Mean(FDR)') +
      theme_linedraw()
    
    recall = ggplot(NULL,aes(x=subset_sizes,recall_matrix_agg[[group_name]])) +
      geom_line() + 
      geom_point() +
      ylim(0,1) +
      ggtitle("Recall") +
      scale_x_log10(breaks=subset_sizes) +
      xlab("Sample Size") +
      ylab('Mean(Recall)') +
      theme_linedraw()

    auroc = ggplot(NULL,aes(x=subset_sizes,auroc_matrix_agg[[group_name]])) +
      geom_line() + 
      geom_point() +
      ylim(0,1) +
      ggtitle("AUROC") +
      scale_x_log10(breaks=subset_sizes) +
      xlab("Sample Size") +
      ylab('Mean(AUROC)') +
      theme_linedraw()
    
    plots = grid.arrange(fdr,recall,auroc,ncol=3, top=paste0(sim_name," ",test," Accuracy metrics for Ab. Scale = ",effect_combs[group_name,1]," and Prev. Shift = ",effect_combs[group_name,2]) )
    #fdr + recall + auroc + plot_annotation(title=paste0("Accuracy metrics for ",test," with Ab. Scale = ",ab.scale," Prev. Shift = ",prev.scale))
    
    filename = paste0("result_plots/",sim_name,"/",test,"/","Ab",effect_combs[group_name,1],"Prev",effect_combs[group_name,2],".png")
    ggsave(filename,plots,create.dir = TRUE)
  }
  return(metrics)
  
}


get_3metrics <- function(sim_file,results_file){
  
  ##Reading in data and parameters
  results=h5read(results_file,name='results')
  sim_params = h5read(sim_file,"simulation_parameters")
  temp = h5ls(results_file,recursive=2)
  sim_name = substr(sim_file,start=10,stop=nchar(sim_file)-3)
  
  #extracting all group names
  groups = setdiff(temp$name, c('original_data', 'simulation_parameters','results')) 
  
  #extracting abundance and prevalence scaling values
  ab.scale = sim_params[['sim_params']][['ab.scale']]
  prev.scale = sim_params[['sim_params']][['prev.scale']]
  
  #extracting group, subset, repeat, and rep values
  no.groups = length(groups)
  no.subsets = length(sim_params[['test_idx']][['subsets']])
  subsets = names(results[['ab1_prev1_rep1']])
  subset_sizes = sort(sim_params[['test_idx']][['subsets']])
  subsets = paste0("subset_",as.character(subset_sizes))
  no.repeats = sim_params[['sim_params']][['repeats']]
  no.reps = length(sim_params[['test_idx']][['names']])-1 
  
  ##Extracting and Aggregating Data
  
  auroc_vals = list()
  recall_vals = list()
  fdr_vals = list()
  
  #iterating through groups to collect AUROC, Recall, and FDR metrics
  for(group in groups){
    
    group_aurocs_vals = list()
    group_recalls_vals = list()
    group_fdrs_vals = list()
    
    #Iterating through each subset of group
    for(subset in subsets){
      
      #Aggregating metrics across repetitions of each subset
      subset_data = results[[group]][[subset]]
      group_aurocs_vals[[subset]]=subset_data[['auroc']]
      
      subset_recalls = subset_data[['TP']]/(subset_data[['TP']] + subset_data[['FN']])
      group_recalls_vals[[subset]]=subset_recalls
      
      subset_fdrs = as.numeric(subset_data[['FP']]/(subset_data[['TP']] + subset_data[['FP']]))
      subset_fdrs = replace(subset_fdrs,which(is.nan(subset_fdrs)),0)
      group_fdrs_vals[[subset]] = subset_fdrs
      
    }
    
    auroc_vals[[group]] = group_aurocs_vals
    recall_vals[[group]] = group_recalls_vals
    fdr_vals[[group]] = group_fdrs_vals
    
  }
  
  #Effect combinations
  ab.scales = rep(ab.scale,rep(length(prev.scale),length(ab.scale)))
  ab.scales.nums = rep(1:length(ab.scale),rep(length(prev.scale),length(ab.scale)))
  
  prev.scales = rep(prev.scale,length(ab.scale))
  prev.scales.nums = rep(1:length(prev.scale),length(ab.scale))
  
  effect_combs=data.frame(ab.scales=ab.scales,
                          prev.scale=prev.scales,
                          ab.scales.nums=ab.scales.nums,
                          prev.scales.nums=prev.scales.nums)
  
  comb_names=c()
  
  for(i in 1:length(rownames(effect_combs))){
    comb_names = c(comb_names,paste0("Ab",effect_combs[i,1],"Prev",effect_combs[i,2]))
  }
  
  rownames(effect_combs) = comb_names
  
  #Extracting mean and variance for each effect/subset combinations
  
  #Setting up dataframes
  melts = list(auroc=melt(auroc_vals)[,2:4], recall = melt(recall_vals)[,2:4], fdr= melt(fdr_vals))
  metrics_means = list()
  metrics_SE = list()
  
  for(metric in c("auroc","recall","fdr")){
    names(melts[[metric]]) = c("value","subset","group")
    
    metric_means = metric_SE = data.frame(matrix(NA, nrow = length(subsets), ncol = length(comb_names)))
    rownames(metric_means) = rownames(metric_SE) = subsets
    colnames(metric_means) = colnames(metric_SE) = comb_names
    
    for(comb_name in comb_names){
      
      #Finding effect combination name to string match on group field
      string_match = paste0("ab",effect_combs[comb_name,3],"_prev",effect_combs[comb_name,4])
      
      #Iterating through subsets
      for(subset_name in subsets){
        
        #Extracting only the data matching the subset and appropriate group name
        effect_data <- melts[[metric]] %>%
          filter(str_detect(group,string_match)) %>%
          filter(str_detect(subset,subset_name))
        
        #Taking Mean and Standard of all these values
        metric_means[subset_name,comb_name] = mean(effect_data$value)
        metric_SE[subset_name,comb_name] = sqrt(var(effect_data$value)/length(effect_data$value))
        
      }
    }
    
    metrics_means[[metric]] = metric_means
    metrics_SE[[metric]] = metric_SE
    
    
  }
  
  return(list(metrics_means=metrics_means,metrics_SE = metrics_SE))
  
}


get_sim_tests <- function(sim_name,tests){
  
  sim_file = paste0("datasets/",sim_name,".h5")
  sim_metrics = list()
  
  for(test in tests){
    
    results_file = paste0("datasets/",sim_name,"_",test,"_results.h5")
    
    test_metrics = get_3metrics(sim_file,results_file)
    
    sim_metrics[[test]] = test_metrics
    
  }

  return(sim_metrics)
  
}

# Helper: shift limma-original lines by a small constant and draw them first (behind)
jitter_original <- function(df, amount = 0.01) {
  is_orig <- grepl("original", df$method, ignore.case = TRUE)
  df$value <- ifelse(is_orig, df$value + amount, df$value)
  df[order(!is_orig), ]   # original rows first => drawn underneath
}


generate_linegraphs_by_effect <- function(sim_name, tests, subsets,grouping){
  sim_metrics = get_sim_tests(sim_name, tests)
  
  for(effect in colnames(sim_metrics[[1]][[1]][[1]])){
    
    effect_df_fdr = data.frame(subsets = subsets)
    effect_df_recall = data.frame(subsets = subsets)
    effect_df_auroc = data.frame(subsets = subsets)
    
    effect_df_fdr_SE = data.frame(subsets = subsets)
    effect_df_recall_SE = data.frame(subsets = subsets)
    effect_df_auroc_SE = data.frame(subsets = subsets)
    
    
    for(test in tests){
      
      means = sim_metrics[[test]][["metrics_means"]]
      SE = sim_metrics[[test]][["metrics_SE"]]
      
      effect_df_fdr[[test]] = means[["fdr"]][[effect]]
      effect_df_recall[[test]] = means[["recall"]][[effect]]
      effect_df_auroc[[test]] = means[["auroc"]][[effect]]
      
      effect_df_fdr_SE[[test]] = SE[["fdr"]][[effect]]
      effect_df_recall_SE[[test]] = SE[["recall"]][[effect]]
      effect_df_auroc_SE[[test]] = SE[["auroc"]][[effect]]
      
    }
    
    effect_df_fdr <- melt(effect_df_fdr,  id.vars = 'subsets')
    effect_df_fdr_SE <- melt(effect_df_fdr_SE,  id.vars = 'subsets')
    effect_df_fdr[['SE']] = effect_df_fdr_SE[['value']]
    
    effect_df_recall <- melt(effect_df_recall,  id.vars = 'subsets')
    effect_df_recall_SE <- melt(effect_df_recall_SE,  id.vars = 'subsets')
    effect_df_recall[['SE']] = effect_df_recall_SE[['value']]
    
    effect_df_auroc <- melt(effect_df_auroc,  id.vars = 'subsets')
    effect_df_auroc_SE <- melt(effect_df_auroc_SE,  id.vars = 'subsets')
    effect_df_auroc[['SE']] = effect_df_auroc_SE[['value']]
    
    #Distinguishing Methods and Normalisations
    
    effect_df_fdr = effect_df_fdr %>%
      mutate(norm = str_match(variable, "^[^_]*_[^_]*_(.*)$")[,2]) %>%
      replace_na(list(norm = "TMM")) %>% 
      mutate(normalisation = factor(norm)) %>% 
      mutate(method = factor(str_extract(variable, "^[^_]*_[^_]*"))) %>% 
      mutate(method = str_replace_all(method, "_"," "))
    
    effect_df_recall = effect_df_recall %>%
      mutate(norm = str_match(variable, "^[^_]*_[^_]*_(.*)$")[,2]) %>%
      replace_na(list(norm = "TMM")) %>% 
      mutate(normalisation = factor(norm)) %>% 
      mutate(method = factor(str_extract(variable, "^[^_]*_[^_]*"))) %>% 
      mutate(method = str_replace_all(method, "_"," "))
    
    effect_df_auroc = effect_df_auroc %>%
      mutate(norm = str_match(variable, "^[^_]*_[^_]*_(.*)$")[,2]) %>%
      replace_na(list(norm = "TMM")) %>% 
      mutate(normalisation = factor(norm)) %>% 
      mutate(method = factor(str_extract(variable, "^[^_]*_[^_]*"))) %>% 
      mutate(method = str_replace_all(method, "_"," "))
    
    #Jitter limma-original so it shows behind limma-trend
    jit = -0.01
    
    # fdr is plotted as 1 - value (precision), so flip the sign to nudge precision upward
    effect_df_fdr    = jitter_original(effect_df_fdr, -jit)
    effect_df_recall = jitter_original(effect_df_recall, jit)
    effect_df_auroc  = jitter_original(effect_df_auroc, jit)
    
    #PLotting
    
    fdr = ggplot(effect_df_fdr, aes(subsets,1-value,col=method,linetype=normalisation,shape=normalisation,group=variable)) +
      geom_line(show.legend=FALSE) + 
      geom_point(size=1.2,show.legend=FALSE) +
      geom_errorbar(aes(ymin=1-(value-1.96*SE),ymax=1-(value+1.96*SE)),show.legend = FALSE,width=0.1) +
      ylim(0,1) +
      ggtitle("Precision") +
      geom_hline(yintercept=0.95,linetype='dashed',colour='red') +
      scale_x_log10(breaks=subsets) +
      xlab("Sample Size") +
      ylab('Mean(Precision)') +
      theme_linedraw()
    
    recall = ggplot(effect_df_recall,aes(subsets,value,col=method,linetype=normalisation,shape=normalisation,group=variable)) +
      geom_line(show.legend=FALSE) + 
      geom_point(size=1.2,show.legend=FALSE) +
      geom_errorbar(aes(ymin=value-1.96*SE,ymax=value+1.96*SE),show.legend=FALSE,width=0.1) +
      ylim(0,1) +
      ggtitle("Recall") +
      scale_x_log10(breaks=subsets) +
      xlab("Sample Size") +
      ylab('Mean(Recall)') +
      theme_linedraw()
    
    auroc = ggplot(effect_df_auroc,aes(x=subsets,y=value,col=method,linetype=normalisation,shape=normalisation,group=variable)) +
      geom_line() + 
      geom_point(size=1.2) +
      geom_errorbar(aes(ymin=value-1.96*SE,ymax=value+1.96*SE),width=0.1) +
      ylim(0,1) +
      ggtitle("AUROC") +
      scale_x_log10(breaks=subsets) +
      xlab("Sample Size") +
      ylab('Mean(AUROC)') +
      theme_linedraw() +
      theme(legend.position = "inside",
            legend.position.inside = c(0.6,0.3))
    
    ab.scale = as.numeric(str_extract(effect, "(?<=Ab)\\d+(?:\\.\\d+)?"))
    prev.scale = as.numeric(str_extract(effect, "(?<=Prev)\\d+(?:\\.\\d+)?"))
    niche = str_extract(sim_name, "(?<=_)[^_]+")
    
    plots = grid.arrange(fdr,recall,auroc,ncol=3, top=paste0(" Accuracy Metrics for Niche : ",niche," with ab.scale = ", ab.scale, " and prev.shift = ",prev.scale) )
  
    filename = paste0("metric_plots/line_by_effect/",sim_name,"_",grouping,"/",effect,".png")
    ggsave(filename,plots,create.dir = TRUE)
    
  }
}


generate_heatmaps_by_test <- function(sim_name, tests, subsets){
  sim_metrics = get_sim_tests(sim_name, tests)
  for(test in tests){
    metrics = c("auroc", "recall", "fdr")
    
    metric_titles = c(auroc="AUROC", recall="Recall", fdr="1-FDR")
    plots_grid = list()
    
    text_scale = 2.5
    base_size = 11 * text_scale  # ggplot's default base_size is 11
    
    for (metric in metrics) {
      
      df = sim_metrics[[test]][["metrics_means"]][[metric]]
      
      df_long = melt(as.matrix(df))
      names(df_long) = c("subset_size", "effect", "value")
      
      df_long = df_long %>% 
        mutate(ab.scale = factor(as.numeric(str_extract(effect, "(?<=Ab)\\d+(?:\\.\\d+)?")))) %>% 
        mutate(prev.scale = factor(as.numeric(str_extract(effect, "(?<=Prev)\\d+(?:\\.\\d+)?")))) %>%
        mutate(subset_size = str_extract(subset_size, "\\d+"))
      
      if(metric =="fdr"){
        df_long = df_long %>% mutate(value = 1-value)
      }
      
      plots_metric = list()
      
      for (i in 1:length(subsets)) {
        subset = subsets[i]
        filtered = df_long %>% filter(subset_size == subset)
        
        p = ggplot(filtered, aes(x = ab.scale, y = prev.scale, fill = value)) +
          geom_tile(show.legend = TRUE) +
          geom_text(aes(label = round(value, 2)), color = "white", size = 4 * text_scale) +
          scale_fill_viridis_c(option = "plasma", limits = c(0, 1.2), name = "Value") +
          ylab("prev.scale") +
          theme_gray(base_size = base_size) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                axis.title.y = element_text(face = "bold"))
        
        # Column title: metric name, only on the top row
        if (i == 1) {
          p = p + ggtitle(metric_titles[[metric]])
        }
        
        plots_metric[[i]] = p
      }
      
      plots_grid[[metric]] = plots_metric
    }
    
    # Build row label plots (horizontal text, blank canvas)
    row_labels = list()
    for (i in 1:length(subsets)) {
      subset = subsets[i]
      row_labels[[i]] = ggplot() +
        annotate("text", x = 0, y = 0, label = paste("Subset size:", subset),
                 angle = 0, hjust = 0.5, size = 4 * text_scale, fontface = "bold") +
        theme_void()
    }
    
    # Interleave: for each row, put label + one plot per metric
    final_plots = list()
    for (i in 1:length(subsets)) {
      final_plots = append(final_plots, list(row_labels[[i]]))
      for (metric in metrics) {
        final_plots = append(final_plots, list(plots_grid[[metric]][[i]]))
      }
    }
    
    combined_plot = wrap_plots(final_plots,
                               nrow = length(subsets),
                               ncol = length(metrics) + 1,
                               byrow = TRUE,
                               axes = 'collect',
                               axis_titles = 'collect',
                               guides = 'collect',
                               widths = c(0.4, rep(1, length(metrics)))) +
      plot_annotation(title = paste0(test," : Accuracy Metric by Ab and Prev Scale Across Subset sizes"),
                      theme = theme(plot.title = element_text(hjust = 0.5, size = 16 * text_scale, face = "bold")))
    
    filename = paste0("metric_plots/heatmaps/",sim_name,"/",test,".png")
    ggsave(filename,combined_plot, width = 30, height = 30, dpi = 300, create.dir = TRUE)
  }
}

generate_linegraphs_by_subset <- function(sim_name, tests,grouping){
  sim_metrics = get_sim_tests(sim_name, tests)
  
  df <- data.frame(subset_size = character(),
                   effect = factor(),
                   AUROC = numeric(),
                   Recall = numeric(),
                   `1-FDR` = numeric(),
                   ab.scale = factor(),
                   prev.scale=factor(),
                   test=factor()
  )
  
  df_se <- data.frame(subset_size = character(),
                      effect = factor(),
                      AUROC = numeric(),
                      Recall = numeric(),
                      `1-FDR` = numeric(),
                      ab.scale = factor(),
                      prev.scale=factor(),
                      test=factor()
  )
  
  colnames(df)[5] = "1-FDR"
  
  for(test in tests){
    metrics = c("auroc", "recall", "fdr")
    
    metric_titles = c(auroc="AUROC", recall="Recall", fdr="1-FDR")
    
    data = sim_metrics[[test]][["metrics_means"]][[metrics[1]]]
    se_data = sim_metrics[[test]][["metrics_SE"]][[metrics[1]]]
    
    df_long = melt(as.matrix(data))
    df_long_se = melt(as.matrix(se_data))
    names(df_long) = c("subset_size", "effect", metric_titles[[metrics[1]]])
    names(df_long_se) = c("subset_size", "effect", metric_titles[[metrics[1]]])
    
    for (metric in metrics[2:length(metrics)]) {
      
      df_temp = sim_metrics[[test]][["metrics_means"]][[metric]]
      df_temp_se = sim_metrics[[test]][["metrics_SE"]][[metric]]
      
      df_temp_long = melt(as.matrix(df_temp))
      names(df_temp_long) = c("subset_size", "effect", "value")
      df_temp_long_se = melt(as.matrix(df_temp_se))
      names(df_temp_long_se) = c("subset_size", "effect", "value")
      
      if(metric == "fdr"){
        df_temp_long = df_temp_long %>% 
          mutate(value = 1 - value)
      }
      
      df_long = df_long %>% 
        mutate(!!metric_titles[[metric]] := df_temp_long[["value"]])
      
      df_long_se = df_long_se %>% 
        mutate(!!metric_titles[[metric]] := df_temp_long_se[["value"]])
    }
    
    df_long = df_long %>% 
      mutate(ab.scale = factor(as.numeric(str_extract(effect, "(?<=Ab)\\d+(?:\\.\\d+)?")))) %>% 
      mutate(prev.scale = factor(as.numeric(str_extract(effect, "(?<=Prev)\\d+(?:\\.\\d+)?")))) %>%
      mutate(subset_size = str_extract(subset_size, "\\d+")) %>% 
      mutate(test = factor(rep(test,dim(df_long)[1])))
    
    df_long_se = df_long_se %>% 
      mutate(ab.scale = factor(as.numeric(str_extract(effect, "(?<=Ab)\\d+(?:\\.\\d+)?")))) %>% 
      mutate(prev.scale = factor(as.numeric(str_extract(effect, "(?<=Prev)\\d+(?:\\.\\d+)?")))) %>%
      mutate(subset_size = str_extract(subset_size, "\\d+")) %>% 
      mutate(test = factor(rep(test,dim(df_long_se)[1])))
    
    df = bind_rows(df,df_long)
    df_se = bind_rows(df_se,df_long_se)
    
  }
  
  df <- df %>%
    pivot_longer(
      cols = c(`1-FDR`,Recall, AUROC),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(metric = factor(metric, levels = c( "1-FDR","Recall","AUROC")),
           prev.scale = factor(prev.scale, levels = rev(sort(unique(prev.scale)))))
  
  df_se <- df_se %>%
    pivot_longer(
      cols = c(`1-FDR`,Recall, AUROC),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(metric = factor(metric, levels = c( "1-FDR","Recall","AUROC")),
           prev.scale = factor(prev.scale, levels = rev(sort(unique(prev.scale)))))
  
  df[['SE']] = df_se[['value']]
  
  df = df %>%
    mutate(norm = str_match(test, "^[^_]*_[^_]*_(.*)$")[,2]) %>%
    replace_na(list(norm = "TMM")) %>% 
    mutate(normalisation = factor(norm)) %>% 
    mutate(method = factor(str_extract(test, "^[^_]*_[^_]*"))) %>% 
    mutate(method = str_replace_all(method, "_"," "))
  
    
  
  # Nudge limma-original up slightly so it shows behind limma-trend
  jit = -0.02
  df = df %>%
    mutate(is_orig = grepl("original", method, ignore.case = TRUE),
           value   = ifelse(is_orig, value + jit, value))
  
  # One set of layers (line + points + error bars) for a given subset of the data
  metric_layers <- function(d) {
    list(
      geom_line(data = d),
      geom_point(data = d, size = 1.5),
      geom_errorbar(data = d, aes(ymin = value - 2*SE, ymax = value + 2*SE), width = 0.3)
    )
  }
  
  
  for (ss in unique(df$subset_size)) {
    
    df_ss <- df %>% filter(subset_size == ss)
    
    p <- ggplot(mapping = aes(x = ab.scale, y = value,
                              col = method, group = test, linetype = normalisation, shape = normalisation)) +
      metric_layers(df_ss %>% filter(is_orig)) +    # drawn first (behind)
      metric_layers(df_ss %>% filter(!is_orig)) +   # drawn on top
      facet_grid(prev.scale ~ metric, scales = "free_y", switch = "y") +
      labs(
        title = paste0("Subset size: ", ss),
        x = "Abundance Scale",
        y = "Prevalence Shift"
      ) +
      theme_linedraw() +
      theme(
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold",color='black'),
        strip.text.y.left = element_text(face = "bold", color = 'black',angle=0),
        strip.placement = "outside",
        axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 7,angle=0)
      )
    
    print(p)  # outputs directly, no saving
    filename = paste0("metric_plots/line_by_subset/",sim_name,"_",grouping,"/subset",ss,".png")
    ggsave(filename,p,create.dir = TRUE)
  }
  
}


