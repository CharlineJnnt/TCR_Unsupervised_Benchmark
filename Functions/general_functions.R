#general functions

#not included in
'%!in%' <- function(x,y)!('%in%'(x,y))

#Self-explanatory
remove_first_column <- function(x) {
  rownames(x) <- x[,1]
  x <- x[,-1]
}

#Self-explanatory
rownames_as_first_column <- function(x, name) {
  y <- as.data.frame(rownames(x))
  colnames(y) <- name
  x <- cbind(y, x)
}

plot_distribution_size <- function(cluster_size, color, method){
  n_cluster <- n_distinct(cluster_size$clusters)
  
  cluster_size <- cluster_size %>% 
    group_by(size) %>% 
    dplyr::mutate(count = n()) %>% 
    select(-clusters) %>% unique()
  
  #pdf(paste0(directory, method, "/beta_bis/cluster_distribution_bis.pdf"), width = 10, height = 5)#barplot + curve count
  plot <- ggplot(data=cluster_size, aes(x=size, y=count)) +
    geom_bar(stat="identity", fill = color, color = "black")+
    geom_line(color ="black")+
    annotate("text", label = paste0("n= ", n_cluster), x = max(cluster_size$size), y = max(cluster_size$count), size = 5, colour = "black", vjust="inward",  hjust = "inward")+
    theme_bw()+
    ggtitle(method) +
    theme(plot.title = element_text(hjust=0.5, face = "bold", size = 20),
          legend.title=element_text(size=15),
          text=element_text(size=15))+
    scale_x_continuous(breaks = seq(0,max(cluster_size$size),10))+ylab("Number of clusters")
  #dev.off()
  return(plot)
}

#add sd and mean on graph
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


UpSet_bis = function(m, 
                     comb_col = "black",
                     pt_size = unit(3, "mm"), lwd = 2, 
                     bg_col = "#F0F0F0", bg_pt_col = "white",
                     set_order = order(set_size(m), decreasing = TRUE), 
                     comb_order = if(attr(m, "param")$set_on_rows) {
                       order.comb_mat(m[set_order, ], decreasing = TRUE)
                     } else {
                       order.comb_mat(m[, set_order], decreasing = TRUE)
                     },
                     top_annotation = upset_top_annotation(m),
                     right_annotation = upset_right_annotation(m),
                     row_names_side = "left", lgd_color,
                     ...) {
  
  param = attr(m, "param")
  set_on_rows = param$set_on_rows
  mode = param$mode
  
  m2 = m
  
  class(m2) = "matrix"
  
  pt_size = pt_size
  lwd = lwd
  
  ###
  lgd = Legend(labels = as.vector(lgd_color), 
               title = "Methods", 
               legend_gp = gpar(fill = names(lgd_color), 
                                col = "black"), 
               labels_gp = gpar(font = 8),
               nrow = 1, 
               title_position = "leftcenter")
  ###
  
  if(set_on_rows) {
    n_comb = ncol(m)
    if(length(comb_col == 1)) comb_col = rep(comb_col, n_comb)
    
    layer_fun = function(j, i, x, y, w, h, fill) {
      nr = round(1/as.numeric(h[1]))
      nc = round(1/as.numeric(w[1]))
      subm = matrix(pindex(m2, i, j), nrow = nr, byrow = FALSE)
      
      for(k in seq_len(nr)) {
        if(k %% 2) {
          grid.rect(y = k/nr, height = 1/nr, just = "top", gp = gpar(fill = bg_col[1], col = NA))
        } else {
          if(length(bg_col) > 1) {
            grid.rect(y = k/nr, height = 1/nr, just = "top", gp = gpar(fill = bg_col[2], col = NA))
          }
        }
      }
      
      
      #lines
      for(k in seq_len(nc)) {
        if(sum(subm[, k]) >= 2) {
          i_min = min(which(subm[, k] > 0))
          i_max = max(which(subm[, k] > 0))
          #grid.lines(c(k - 0.5, k - 0.5)/nc, (nr - c(i_min, i_max) + 0.5)/nr, gp = gpar(col = comb_col[jj[k]], lwd = lwd))#intersect lines between points
          grid.lines(c(k - 0.5, k - 0.5)/nc, (nr - c(i_min, i_max) + 0.5)/nr, gp = gpar(col = "black", lwd = lwd))
        }
      }
      
      
      
      grid.points(x, y, size = pt_size, 
                  pch = 21, gp = gpar(fill = ifelse(pindex(m2, i, j), comb_col[i], bg_pt_col)))#comb_col[j]
      
      #jj = unique(j)
      
    }
    
    # check top annotation
    # if it is specified by upset_top_annotation and gp(col) is not set
    ra = top_annotation
    if(length(ra) == 1) {
      ta_call = substitute(top_annotation)
      ta_call = as.list(ta_call)
      if(as.character(ta_call[[1]]) == "upset_top_annotation") {
        if(!"gp" %in% names(as.list(ta_call))) {
          ra@anno_list[[1]]@fun@var_env$gp$fill = comb_col
          ra@anno_list[[1]]@fun@var_env$gp$col = comb_col
        }
      }
    }
    
    ht = Heatmap(m2, cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(type = "none"),
                 layer_fun = layer_fun, show_heatmap_legend = FALSE,
                 top_annotation = ra,
                 right_annotation = right_annotation,
                 row_names_side = row_names_side,
                 row_order = set_order, column_order = comb_order, show_row_names = FALSE,...)
  } else {
    n_comb = nrow(m)
    if(length(comb_col == 1)) comb_col = rep(comb_col, n_comb)
    
    layer_fun = function(j, i, x, y, w, h, fill) {
      nr = round(1/as.numeric(h[1]))
      nc = round(1/as.numeric(w[1]))
      subm = matrix(pindex(m2, i, j), nrow = nr, byrow = FALSE)
      grid.points(x, y, size = pt_size, pch = 21, gp = gpar(fill = ifelse(pindex(m2, i, j), comb_col[j], "#CCCCCC")))
      for(k in seq_len(nr)) {
        if(sum(subm[k, ]) >= 2) {
          i_min = min(which(subm[k, ] > 0))
          i_max = max(which(subm[k, ] > 0))
          grid.lines((c(i_min, i_max) - 0.5)/nc, (nr - c(k ,k) + 0.5)/nr, gp = gpar(col = "black", lwd = lwd))#col = "black
        }
      }
    }
    
    ra = right_annotation
    if(length(ra) == 1) {
      ta_call = substitute(top_annotation)
      ta_call = as.list(ta_call)
      if(as.character(ta_call[[1]]) == "upset_right_annotation") {
        if(!"gp" %in% names(as.list(ta_call))) {
          ra@anno_list[[1]]@fun@var_env$gp$fill = comb_col
          ra@anno_list[[1]]@fun@var_env$gp$col = comb_col
        }
      }
    }
    
    ht = Heatmap(m2, cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(type = "none"),
                 layer_fun = layer_fun, show_heatmap_legend = FALSE,
                 top_annotation = top_annotation,
                 right_annotation = ra,
                 row_order = comb_order, column_order = set_order, ...)
  }
  draw(ht, annotation_legend_list = list(lgd), annotation_legend_side = "bottom")
}

count_spe_by_cluster_one_chain <- function(output){
  output_spe <- merge(output, dataframe_alpha_beta[,c("CDR3", "Epitope", "chain")], by = "CDR3") %>% unique()
  output_spe <- output_spe %>% group_by(clusters) %>% mutate(nb_spe_cluster = n_distinct(Epitope))
  output_spe <- output_spe %>% select(clusters, nb_spe_cluster) %>% unique()
  df_sum <- output_spe %>% 
    group_by(nb_spe_cluster) %>% 
    summarise(count_cluster = n(), .groups = 'drop') 
  df_sum <- df_sum %>% mutate(perc_cluster = count_cluster/sum(df_sum$count_cluster)*100)
  return(df_sum)
}

count_spe_by_cluster_pairing <- function(output, colonne_pair){
  output_spe <- merge(output[,c(colonne_pair, "clusters", "size")], data_human[,c("CDR3b_CDR3a", "Epitope")], by.x = colonne_pair, by.y = "CDR3b_CDR3a") %>% unique()
  output_spe <- output_spe %>% group_by(clusters) %>% mutate(nb_spe_cluster = n_distinct(Epitope))
  output_spe <- output_spe %>% select(clusters, nb_spe_cluster) %>% unique()
  df_sum <- output_spe %>% 
    group_by(nb_spe_cluster) %>% 
    summarise(count_cluster = n(), .groups = 'drop') 
  df_sum <- df_sum %>% mutate(perc_cluster = count_cluster/sum(df_sum$count_cluster)*100)
  return(df_sum)
}

select_one_spe_by_seq <- function(output){
  output_spe <- merge(output, dataframe_alpha_beta[,c("CDR3", "Epitope", "chain")], by = "CDR3") %>% unique()
  output_spe <- output_spe %>% group_by(clusters) %>% mutate(nb_spe_cluster = n_distinct(Epitope))
  output_spe <- output_spe %>% group_by(clusters)  %>% mutate(epi_majoritaire = Epitope[which.max(nb_spe_cluster)])
  
  output_spe_test <- split(output_spe, output_spe$clusters)
  output_spe_new <- list()
  for(cluster in names(output_spe_test)){
    output_spe_seq <- split(output_spe_test[[cluster]], output_spe_test[[cluster]]$CDR3)
    row_toKeep_list <- list()
    for(i in names(output_spe_seq)){
      #print(i)
      if(nrow(output_spe_seq[[i]])>1){#if multiple specificities
        if(output_spe_seq[[i]]$Epitope==output_spe_seq[[i]]$epi_majoritaire[1]){
          row_toKeep_df <- output_spe_seq[[i]] %>% filter(Epitope == output_spe_seq[[i]]$epi_majoritaire[1])#keep specific epitope
        }else{row_toKeep_df <- output_spe_seq[[i]][1,]}#otherwise keep the first row
        
      }else{row_toKeep_df <- output_spe_seq[[i]][1,]}#otherwise keep the first row
      
      row_toKeep_list[[i]] <- row_toKeep_df 
    }
    
    output_spe_keep <- rbindlist(row_toKeep_list) %>% unique() %>% drop_na() %>% as.data.frame()
    output_spe_new[[cluster]] <- output_spe_keep
  }
  
  output_spe_bis <- rbindlist(output_spe_new)
  
  output_spe_bis <- output_spe_bis %>% group_by(clusters) %>% mutate(nb_spe_cluster_new = n_distinct(Epitope))
  output_spe_new <- output_spe_bis %>% select(clusters, nb_spe_cluster_new) %>% unique()
  df_sum <- output_spe_new  %>% 
    group_by(nb_spe_cluster_new ) %>% 
    summarise(count_cluster = n(), .groups = 'drop') 
  df_sum <- df_sum %>% mutate(perc_cluster = count_cluster/sum(df_sum$count_cluster)*100)
  return(df_sum)
}

select_one_spe_by_seq_bis <- function(output, column_pair){
  output_spe <- merge(output[,c(column_pair, "clusters", "size")], data_human[,c("CDR3b_CDR3a", "Epitope")], by.x = column_pair, by.y = "CDR3b_CDR3a") %>% unique()
  colnames(output_spe)[1] <- "CDR3b_CDR3a"
  output_spe <- output_spe %>% group_by(clusters) %>% mutate(nb_spe_cluster = n_distinct(Epitope))
  output_spe <- output_spe %>% group_by(clusters)  %>% mutate(epi_majoritaire = Epitope[which.max(nb_spe_cluster)])
  
  output_spe_test <- split(output_spe, output_spe$clusters)
  output_spe_new <- list()
  for(cluster in names(output_spe_test)){
    output_spe_seq <- split(output_spe_test[[cluster]], output_spe_test[[cluster]]$CDR3b_CDR3a)
    row_toKeep_list <- list()
    for(i in names(output_spe_seq)){
      #print(i)
      if(nrow(output_spe_seq[[i]])>1){#if multiple specificities
        if(output_spe_seq[[i]]$Epitope==output_spe_seq[[i]]$epi_majoritaire[1]){
          row_toKeep_df <- output_spe_seq[[i]] %>% filter(Epitope == output_spe_seq[[i]]$epi_majoritaire[1])#keep specific epitope
        }else{row_toKeep_df <- output_spe_seq[[i]][1,]}#otherwise keep the first row
        
      }else{row_toKeep_df <- output_spe_seq[[i]][1,]}#otherwise keep the first row
      
      row_toKeep_list[[i]] <- row_toKeep_df 
    }
    
    output_spe_keep <- rbindlist(row_toKeep_list) %>% unique() %>% drop_na() %>% as.data.frame()
    output_spe_new[[cluster]] <- output_spe_keep
  }
  
  output_spe_bis <- rbindlist(output_spe_new)
  
  output_spe_bis <- output_spe_bis %>% group_by(clusters) %>% mutate(nb_spe_cluster_new = n_distinct(Epitope))
  output_spe_new <- output_spe_bis %>% select(clusters, nb_spe_cluster_new) %>% unique()
  df_sum <- output_spe_new  %>% 
    group_by(nb_spe_cluster_new ) %>% 
    summarise(count_cluster = n(), .groups = 'drop') 
  df_sum <- df_sum %>% mutate(perc_cluster = count_cluster/sum(df_sum$count_cluster)*100)
  return(df_sum)
}
