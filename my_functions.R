## function pathway_errorbar3


pathway_errorbar3<-function (abundance, daa_results_df, Group, ko_to_kegg = FALSE, 
                             p_values_threshold = 0.05, order = "group", select = NULL, 
                             p_value_bar = TRUE, colors = NULL, x_lab = NULL) 
{
  # browser()
  print("pathway_errorbar3 is running!\n")
  if (is.null(x_lab)) {
    if (ko_to_kegg == TRUE) {
      x_lab <- "pathway_name"
    }
    else {
      x_lab <- "description"
    }
  }
  #browser()
  Group<-as.vector(Group) %>% unlist()
  if (is.null(colors)) {
    mycolors <- c("#d93c3e", "#3685bc", "#6faa3e", "#e8a825",  "#c973e6", "#ee6b3d", "#2db0a7", "#f25292")[1:nlevels(as.factor(Group))]
  }
  errorbar_abundance_mat <- as.matrix(abundance)
  daa_results_filtered_df <- daa_results_df[daa_results_df$p_adjust <  p_values_threshold, ]
  if (!is.null(select)) {
    daa_results_filtered_sub_df <- daa_results_filtered_df[daa_results_filtered_df$feature %in% 
                                                             select, ]
  }else {
    daa_results_filtered_sub_df <- daa_results_filtered_df
  }
  if (nrow(daa_results_filtered_sub_df) > 30) {
    stop(paste0("The feature with statistically significance are more than 30, the visualization will be terrible.\n Please use select to reduce the number.\n Now you have ", 
                paste(paste0("\"", daa_results_filtered_sub_df$feature, 
                             "\""), collapse = ", ")))
  }
  errorbar_sub_abundance_mat <- errorbar_abundance_mat[rownames(errorbar_abundance_mat) %in% 
                                                         daa_results_filtered_sub_df$feature, ]
  # browser()
  errorbar_sub_relative_abundance_mat <- as.matrix(as.data.frame(transform_sample_counts(phyloseq::otu_table(errorbar_sub_abundance_mat %>% as.matrix(), taxa_are_rows = TRUE), function(x) x/sum(x))))
  error_bar_matrix <- cbind(sample = colnames(errorbar_sub_relative_abundance_mat), 
                            group = Group, t(errorbar_sub_relative_abundance_mat))
  error_bar_df <- as.data.frame(error_bar_matrix)
  colnames(error_bar_df)[2]<-"group"
  #browser()
  error_bar_pivot_longer_df <- tidyr::pivot_longer(error_bar_df,   -c(sample, group))
  error_bar_pivot_longer_tibble <- mutate(error_bar_pivot_longer_df, group = as.factor(group))
  error_bar_pivot_longer_tibble$sample <- factor(error_bar_pivot_longer_tibble$sample)
  error_bar_pivot_longer_tibble$name <- factor(error_bar_pivot_longer_tibble$name)
  error_bar_pivot_longer_tibble$value <- as.numeric(error_bar_pivot_longer_tibble$value)
  error_bar_pivot_longer_tibble_summarised <- error_bar_pivot_longer_tibble %>% 
    group_by(name, group) %>% summarise(mean = mean(value), 
                                        sd = stats::sd(value))
  error_bar_pivot_longer_tibble_summarised <- error_bar_pivot_longer_tibble_summarised %>% 
    mutate(group2 = "nonsense")
  # browser()
  switch(order, p_values = {
    order <- order(daa_results_filtered_sub_df$p_adjust)
  }, name = {
    order <- order(daa_results_filtered_sub_df$feature)
  }, group = {
    daa_results_filtered_sub_df$pro <- 1
    for (i in levels(error_bar_pivot_longer_tibble_summarised$name)) {
      error_bar_pivot_longer_tibble_summarised_sub <- error_bar_pivot_longer_tibble_summarised[error_bar_pivot_longer_tibble_summarised$name ==  i, ]
      pro_group <- error_bar_pivot_longer_tibble_summarised_sub[error_bar_pivot_longer_tibble_summarised_sub$mean == 
                                                                  max(error_bar_pivot_longer_tibble_summarised_sub$mean), ]$group
      pro_group <- as.vector(pro_group)
      daa_results_filtered_sub_df[daa_results_filtered_sub_df$feature == 
                                    i, ]$pro <- pro_group
    }
    order <- order(daa_results_filtered_sub_df$pro, daa_results_filtered_sub_df$p_adjust)
  }, pathway_class = {
    if (!"pathway_class" %in% colnames(daa_results_filtered_sub_df)) {
      stop("Please use pathway_annotation function to annotate the pathway_daa results")
    }
    order <- order(daa_results_filtered_sub_df$pathway_class, 
                   daa_results_filtered_sub_df$p_adjust)
  }, {
    order <- order
  })
  
  daa_results_filtered_sub_df <- daa_results_filtered_sub_df[order, ]
  error_bar_pivot_longer_tibble_summarised_ordered <- data.frame(name = NULL, 
                                                                 group = NULL, mean = NULL, sd = NULL)
  for (i in daa_results_filtered_sub_df$feature) {
    error_bar_pivot_longer_tibble_summarised_ordered <- rbind(error_bar_pivot_longer_tibble_summarised_ordered, error_bar_pivot_longer_tibble_summarised[error_bar_pivot_longer_tibble_summarised$name ==  i, ])
  }
  if (ko_to_kegg == FALSE) {
    error_bar_pivot_longer_tibble_summarised_ordered[, x_lab] <- rep(daa_results_filtered_sub_df[,  x_lab], each = length(levels(factor(error_bar_pivot_longer_tibble_summarised_ordered$group))))
  }
  if (ko_to_kegg == TRUE) {
    error_bar_pivot_longer_tibble_summarised_ordered$pathway_class <- rep(daa_results_filtered_sub_df$pathway_class, 
                                                                          each = length(levels(factor(error_bar_pivot_longer_tibble_summarised_ordered$group))))
  }
  error_bar_pivot_longer_tibble_summarised_ordered$name <- factor(error_bar_pivot_longer_tibble_summarised_ordered$name, 
                                                                  levels = rev(daa_results_filtered_sub_df$feature))
  bar_errorbar <- ggplot2::ggplot(error_bar_pivot_longer_tibble_summarised_ordered, ggplot2::aes(mean, name, fill = group)) + 
    ggplot2::geom_errorbar(ggplot2::aes(xmax = mean + sd, xmin = 0), position = ggplot2::position_dodge(width = 0.8), width = 0.5, size = 0.5, color = "black") + 
    ggplot2::geom_bar(stat = "identity",  position = ggplot2::position_dodge(width = 0.8), width = 0.8) + 
    GGally::geom_stripped_cols() + ggplot2::scale_fill_manual(values = mycolors) + 
    ggplot2::scale_color_manual(values = mycolors) + ggprism::theme_prism() + 
    ggplot2::scale_x_continuous(expand = c(0, 0), guide = "prism_offset_minor", 
    ) + ggplot2::scale_y_discrete(labels = rev(daa_results_filtered_sub_df[, x_lab])) +
    ggplot2::labs(x = "Relative Abundance(%)",  y = NULL) + 
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank(), 
                   axis.line.y = ggplot2::element_blank(), axis.line.x = ggplot2::element_line(size = 0.5), 
                   axis.ticks.x = ggplot2::element_line(size = 0.5), panel.grid.major.y = ggplot2::element_blank(), 
                   panel.grid.major.x = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 10, 
                                                                                                    color = "black"), axis.text.x = ggplot2::element_text(margin = ggplot2::margin(r = 0)), 
                   axis.text.y = ggplot2::element_text(size = 10, color = "black", 
                                                       margin = ggplot2::margin(b = 6)), axis.title.x = ggplot2::element_text(size = 10, color = "black", hjust = 0.5), legend.position = "top", legend.key.size = ggplot2::unit(0.1, "cm"), legend.direction = "vertical",  legend.justification = "left", legend.text = ggplot2::element_text(size = 8, face = "bold"), legend.box.just = "right", plot.margin = ggplot2::margin(0, 0.5, 0.5, 0, unit = "cm")) +
    ggplot2::coord_cartesian(clip = "off")
  if (ko_to_kegg == TRUE) {
    pathway_class_group <- daa_results_filtered_sub_df$pathway_class %>% 
      table() %>% data.frame()
    start <- c(1, rev(pathway_class_group$Freq)[1:(length(pathway_class_group$Freq) -  1)]) %>% cumsum()
    end <- cumsum(rev(pathway_class_group$Freq))
    ymin <- start - 1/2
    ymax <- end + 1/2
    nPoints <- length(start)
    pCol <- c("#D51F26", "#272E6A", "#208A42", "#89288F", 
              "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC", 
              "#90D5E4", "#89C75F", "#F37B7D", "#9983BD", "#D24B27", 
              "#3BBCA8", "#6E4B9E", "#0C727C", "#7E1416", "#D8A767", 
              "#3D3D3D")[1:nPoints]
    pFill <- pCol
    for (i in 1:nPoints) {
      bar_errorbar <- bar_errorbar + ggplot2::annotation_custom(grob = grid::rectGrob(gp = grid::gpar(col = pCol[i], 
                                                                                                      fill = pFill[i], lty = NULL, lwd = NULL, alpha = 0.2)), xmin = ggplot2::unit(-2, "native"), xmax = ggplot2::unit(0,  "native"), ymin = ggplot2::unit(ymin[i], "native"), ymax = ggplot2::unit(ymax[i], "native"))
    }
  }
  # browser()
  daa_results_filtered_sub_df <- cbind(daa_results_filtered_sub_df, 
                                       negative_log10_p = -log10(daa_results_filtered_sub_df$p_adjust), 
                                       group_nonsense = "nonsense", log_2_fold_change = NA)
  for (i in daa_results_filtered_sub_df$feature) {
    mean <- error_bar_pivot_longer_tibble_summarised_ordered[error_bar_pivot_longer_tibble_summarised_ordered$name %in% i, ]$mean
    daa_results_filtered_sub_df[daa_results_filtered_sub_df$feature ==  i, ]$log_2_fold_change <- log2(mean[1]/mean[2])
  }
  daa_results_filtered_sub_df$feature <- factor(daa_results_filtered_sub_df$feature, 
                                                levels = rev(daa_results_filtered_sub_df$feature))
  
  p_values_bar <- daa_results_filtered_sub_df %>% ggplot2::ggplot(ggplot2::aes(feature,  log_2_fold_change, fill = group_nonsense)) + 
    ggplot2::geom_bar(stat = "identity",  position = ggplot2::position_dodge(width = 0.8), width = 0.8) + 
    ggplot2::labs(y = "log2 fold change", x = NULL) + GGally::geom_stripped_cols() + 
    ggplot2::scale_fill_manual(values = "#87ceeb") + ggplot2::scale_color_manual(values = "#87ceeb") + 
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = "dashed", color = "black") + ggprism::theme_prism() +
    # ggplot2::scale_y_continuous(expand = c(0,  0), guide = "prism_offset_minor") + 
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),  axis.line.y = ggplot2::element_blank(), axis.line.x = ggplot2::element_line(size = 0.5),    axis.ticks.x = ggplot2::element_line(size = 0.5), panel.grid.major.y = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 10, color = "black"), axis.text.y = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 10, color = "black",  margin = ggplot2::margin(b = 6)), axis.title.x = ggplot2::element_text(size = 11,  color = "black", hjust = 0.5), legend.position = "non") + 
    ggplot2::coord_flip()
  #browser()
  if (ko_to_kegg == TRUE) {
    pathway_class_y <- (ymax + ymin)/2 - 0.5
    pathway_class_plot_df <- as.data.frame(cbind(nonsense = "nonsense", 
                                                 pathway_class_y = pathway_class_y, pathway_class = rev(unique(daa_results_filtered_sub_df$pathway_class))))
    pathway_class_plot_df$pathway_class_y <- as.numeric(pathway_class_plot_df$pathway_class_y)
    pathway_class_annotation <- pathway_class_plot_df %>% 
      ggplot2::ggplot(ggplot2::aes(nonsense, pathway_class_y)) + 
      ggplot2::geom_text(ggplot2::aes(nonsense, pathway_class_y, 
                                      label = pathway_class), size = 3.5, color = "black", 
                         fontface = "bold", family = "sans") + ggplot2::scale_y_discrete(position = "right") + 
      ggprism::theme_prism() + ggplot2::theme(axis.ticks = ggplot2::element_blank(), 
                                              axis.line = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank(), 
                                              panel.grid.major.x = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                                              axis.text = ggplot2::element_blank(), plot.margin = ggplot2::unit(c(0, 0.2, 0, 0), "cm"), axis.title.y = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(), legend.position = "non")
  }
  daa_results_filtered_sub_df$p_adjust <- as.character(daa_results_filtered_sub_df$p_adjust)
  daa_results_filtered_sub_df$unique <- nrow(daa_results_filtered_sub_df) - seq_len(nrow(daa_results_filtered_sub_df)) + 1
  daa_results_filtered_sub_df$p_adjust <- substr(daa_results_filtered_sub_df$p_adjust, 1, 5)
  p_annotation <- daa_results_filtered_sub_df %>% ggplot2::ggplot(ggplot2::aes(group_nonsense, p_adjust)) + 
    ggplot2::geom_text(ggplot2::aes(group_nonsense,  unique, label = p_adjust), size = 3.5, color = "black", fontface = "bold", family = "sans") + ggplot2::labs(y = "p-value (adjusted)") + 
    ggplot2::scale_y_discrete(position = "right") + ggprism::theme_prism() + 
    ggplot2::theme(axis.ticks = ggplot2::element_blank(), 
                   axis.line = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank(), 
                   panel.grid.major.x = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                   axis.text = ggplot2::element_blank(), plot.margin = ggplot2::unit(c(0,  0.2, 0, 0), "cm"), axis.title.y = ggplot2::element_text(size = 11,   color = "black", vjust = 0), axis.title.x = ggplot2::element_blank(), legend.position = "non")
  #browser()
  if (p_value_bar == TRUE) {
    if (ko_to_kegg == TRUE) {
      combination_bar_plot <- pathway_class_annotation +
        bar_errorbar + p_values_bar + p_annotation +
        patchwork::plot_layout(ncol = 4, widths = c(1,1.2, 0.5, 0.1))
      
    }else {
      combination_bar_plot <- bar_errorbar + p_values_bar + 
        p_annotation + patchwork::plot_layout(ncol = 3, 
                                              widths = c(2.3, 0.7, 0.3))
    }
  }else {
    print("no foldchange plot \n")
    #combination_bar_plot <- bar_errorbar + p_annotation + patchwork::plot_layout(ncol = 2, widths = c(2.5,  0.2))
    
    ##no foldchange plot
    combination_bar_plot <- pathway_class_annotation + 
      bar_errorbar  + p_annotation + 
      patchwork::plot_layout(ncol = 3, widths = c(1,1.2, 0.1))#+ p_values_bar
  }
  return(combination_bar_plot)
}



## function ggpicrust3

ggpicrust3<-function (file, metadata, group, pathway, daa_method = "ALDEx2", 
                      ko_to_kegg = FALSE, p.adjust = "BH", order = "group", p_value_bar = TRUE, 
                      x_lab = "pathway_name", select = NULL, reference = NULL, outname=outname,
                      colors = NULL) 
{
  plot_result_list <- list()
  # browser()
  switch(ko_to_kegg, `TRUE` = {
    plot_result_list <- list()
    abundance <- ko2kegg_abundance(file)
    print("ko2kegg_abundance finished!\n")
    # browser()
    daa_results_df <- pathway_daa(abundance = abundance, 
                                  metadata = metadata, group = group, daa_method = daa_method, 
                                  select = select, p.adjust = p.adjust, reference = reference)
    
    write.csv(daa_results_df,paste(outname,"_daa_results_df",group,daa_method,".csv",sep = "_"),row.names = FALSE  )
    
    daa_results_df_count=dplyr::filter(daa_results_df,p_adjust <= 0.05)
    print(paste0("length of daa_results_df(P<=0.05): ",nrow(daa_results_df_count),"\n") )
    print("daa_results_df finished!\n")
    
    if (sum(as.numeric(daa_results_df$p_adjust <= 0.05)) == 0) {
      stop("There are no statistically significant biomarkers")
    }
    if (x_lab == "pathway_name") {
      daa_results_df <- pathway_annotation(daa_results_df = daa_results_df, ko_to_kegg = TRUE)
      print("pathway_annotation finished!\n")
    }
    j <- 1
    
    for (i in unique(daa_results_df$method)  ) {
      print(paste0("use method:",i,"\n"))
      # if(i=="ALDEx2_Wilcoxon rank test"){browser()}
      daa_sub_method_results_df <- daa_results_df[daa_results_df[, "method"] == i, ]
      daa_sub_method_results_df_sorted <- data.frame()
      
      if (is.null(select)) {
        daa_sub_method_results_df_sorted <- daa_sub_method_results_df
      } else if (select == "Top 10") {
        daa_sub_method_results_df_sorted <- daa_sub_method_results_df[order(daa_sub_method_results_df$p_adjust),  ]
        daa_sub_method_results_df_sorted <- daa_sub_method_results_df_sorted%>% head(10)
      } else if (select == "Top 20") {
        daa_sub_method_results_df_sorted <- daa_sub_method_results_df[order(daa_sub_method_results_df$p_adjust), ]
        daa_sub_method_results_df_sorted <- daa_sub_method_results_df_sorted%>% head(20)
      } else if (select == "Top 30") {
        daa_sub_method_results_df_sorted <- daa_sub_method_results_df[order(daa_sub_method_results_df$p_adjust),  ]
        daa_sub_method_results_df_sorted <- daa_sub_method_results_df_sorted %>% head(30)
      }
      # browser()
      write.csv(daa_sub_method_results_df_sorted,paste(outname,"_daa_results_df_sorted_padj_",group,daa_method,".csv",sep = "_"),row.names = FALSE  )
      
      combination_bar_plot <- pathway_errorbar3(abundance = abundance, 
                                                daa_results_df = daa_sub_method_results_df_sorted, 
                                                Group = metadata[, group], ko_to_kegg = ko_to_kegg, 
                                                p_value_bar = p_value_bar,
                                                order = "pathway_class", colors = colors, x_lab = x_lab)
      sub_list <- list(plot = combination_bar_plot, results = daa_sub_method_results_df)
      plot_result_list[[j]] <- sub_list
      j <- j + 1
    }
    return(plot_result_list)
  }, `FALSE` = {
    # browser()
    plot_result_list <- list()
    abundance <- readr::read_delim(file, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    abundance <- tibble::column_to_rownames(abundance, var = "NAME")
    daa_results_df <- pathway_daa(abundance = abundance, 
                                  metadata = metadata, group = group, daa_method = daa_method, 
                                  select = select, p.adjust = p.adjust, reference = reference)
    daa_results_df_count=dplyr::filter(daa_results_df,p_adjust <= 0.05)
    print(paste0("length of daa_results_df(P<=0.05): ",nrow(daa_results_df_count),"\n") )
    
    daa_results_df2 <- pathway_annotation(pathway = pathway, 
                                         ko_to_kegg = FALSE, daa_results_df = daa_results_df)
    j <- 1
    for (i in unique(daa_results_df$method)) {
      daa_sub_method_results_df <- daa_results_df[daa_results_df[, "method"] == i, ]
      combination_bar_plot <- pathway_errorbar3(abundance, 
                                                daa_sub_method_results_df, metadata[, group], 
                                                p_value_bar = p_value_bar,
                                                ko_to_kegg = ko_to_kegg, order = order, colors = colors, 
                                                x_lab = x_lab)
      sub_list <- list(plot = combination_bar_plot, results = daa_sub_method_results_df)
      plot_result_list[[j]] <- sub_list
      j <- j + 1
    }
    return(plot_result_list)
  })
}

ko2kegg_abundance<-function(file){
  
  file_format <- substr(file, nchar(file) - 3, nchar(file))
  switch(file_format, .txt = abundance <- readr::read_delim(file, 
                                                            delim = "\t", escape_double = FALSE, trim_ws = TRUE), 
         .tsv = abundance <- readr::read_delim(file, delim = "\t", 
                                               escape_double = FALSE, trim_ws = TRUE), .csv = abundance <- readr::read_delim(file, 
                                                                                                                             delim = "\t", escape_double = FALSE, trim_ws = TRUE), 
         stop("Error: Please input file as .tsv, .txt or .csv\nThe best input file is what you get from picrust2 output file 'pred_metagenome_unstrat.tsv'"))
  message("Calculation may take a long time, please be patient.")
  load(system.file("extdata", "kegg_reference.RData", package = "ggpicrust2"))
  sample_names <- colnames(abundance)[-1]
  kegg_names <- ko_to_kegg_reference[, 1]
  print(ko_to_kegg_reference)
  kegg_abundance <- matrix(NA, nrow = nrow(kegg_names), ncol = length(sample_names))
  colnames(kegg_abundance) <- sample_names
  rownames(kegg_abundance) <- as.matrix(kegg_names)
  for (i in seq_len(nrow(kegg_abundance))) {
    for (j in seq_len(ncol(kegg_abundance))) {
      print(i)
      print(j)
      kegg_name <- rownames(kegg_abundance)[i]
      sample_name <- colnames(kegg_abundance)[j]
      ko_to_kegg <- ko_to_kegg_reference[ko_to_kegg_reference[,  1] == kegg_name, -1]
      ko_to_kegg <- ko_to_kegg[!is.na(ko_to_kegg)]
      count_one_clumn<-abundance[ as.matrix(abundance[,  1]) %in% ko_to_kegg, sample_name]
      kegg_abundance[i, j] <- sum(count_one_clumn[,1], na.rm = TRUE)
      print(j)
    }
  }
  print(kegg_abundance)
  # browser()
  kegg_abundance <- kegg_abundance[rowSums(kegg_abundance) !=  0, ]
  message("The kegg pathway with zero abundance in all the different samples has been removed.")
  kegg_abundance <- as.data.frame(kegg_abundance)
  return(kegg_abundance)
  
  
  
}
