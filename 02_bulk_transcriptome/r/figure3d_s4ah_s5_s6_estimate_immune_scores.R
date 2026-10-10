# -----------------------------------------------------------------------------
# ESTIMATE purity and immune-infiltration scores
#
# Purpose   : Merge sample groups with ESTIMATE tumour-purity scores (rna_estimate) and immune-cell infiltration scores (rna_immune), then compare tumour vs adjacent normal, tumour only, and normal only across the three progression groups.
#
# Inputs    : objects from the upstream R session: cdata (sample metadata), rna_estimate (ESTIMATE scores), rna_immune (immune-cell infiltration scores). No files are read from disk.
#
# Outputs   : box-plot PDFs written into the sub-directories 2tumor+normal/, 2tumor/ and 2normal/ (created relative to the working directory)
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript figure3d_s4ah_s5_s6_estimate_immune_scores.R
# -----------------------------------------------------------------------------

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")

setwd(file.path(PROJECT_ROOT, "bulk_rnaseq/immune_infiltration"))
conflicts_prefer(dplyr::select, dplyr::filter, dplyr::rename)# declare which function wins on conflict
clin_T <- cdata |> 
  filter(RNA_Tumor != "NG") |> 
  filter(group!="NG") |>
  select(RNA_Tumor, group) |> 
  rename(ID = RNA_Tumor,GP = group)
clin_N <- cdata |> 
  filter(RNA_Normal != "NG") |> 
  filter(group!="NG") |>
  select(RNA_Normal, group) |> 
  rename(ID = RNA_Normal,GP = group)
clin <- rbind(clin_T, clin_N)
clin <-  clin |> 
  column_to_rownames("ID")
score <- merge(clin,rna_estimate,all.x = T,by ="row.names")
idata <- rna_immune
idata <- as.data.frame(t(idata))
data <- merge(score,idata,by.x = "Row.names",by.y="row.names", all = TRUE)
rownames(data) <- data$Row.names
data <- subset(data, select = -Row.names)
data <- data %>% 
  mutate_if(is.character, as.factor)
data <- data |> filter(!is.na(GP))

# 1.1 First compare immune infiltration between tumour and adjacent normal, ignoring progression status
# Tumor+Normal
dir.create("./2tumor+normal")
# continuous variables
data$T_N <- ifelse(grepl("_N", row.names(data)), "N", ifelse(grepl("_T", row.names(data)), "T", NA))
vars <- colnames(data)[!grepl("Pattern|GP|T_N", colnames(data))]
plot_list <- list()
for (var in vars) 
{
  data2 = data[,c(var,"T_N")]
  data2 <- na.omit(data2)
  colnames(data2)[1] <- "value"
  p <- ggplot(data2, aes(x = T_N, y = value, fill = T_N)) + # when mapping, do not pass a variable name, otherwise the loop would keep the last variable of the for loop
    geom_boxplot(alpha = 1) +  # add the box plot
    geom_jitter(aes(group = T_N),alpha = 0.5,position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.75)) +  # two-level grouped points: the x axis is split by one discrete variable, and within each group another discrete variable splits it again (built into ggplot2 as position_jitterdodge)
    labs(x = "", y =  var, fill = NULL) +
    scale_fill_manual(values = c("#2874C5","#C6524A")) +
    #scale_y_continuous(expand = c(0, 0.08)) +
    theme_bw() +
    theme(axis.text=element_text(size=14,color="black"),
          axis.title.y=element_text(size=14,color ="black"),
          axis.text.x = element_text(),
          strip.text = element_text(size = 14),
          strip.background = element_blank(),
          axis.title.x = element_blank(),
          #legend.key.size = unit(0.3, "cm"),
          legend.position = "none" ) + 
    stat_compare_means(aes(group = T_N),
                       label = "p.signif",
                       hide.ns = FALSE,
                       method = "wilcox",
                       label.x = c(1.5)
    )
  plot_list[[var]] <- p
}
plots <- do.call(grid.arrange, c(plot_list, ncol = 6))
filename <- paste("./2tumor+normal/TvsN_without_GP_BoxPlot.pdf", sep = "")
ggsave(filename, plots, width = 12, height = 18)

# Tumor+Normal
dir.create("./2tumor+normal")
# continuous variables
data$T_N <- ifelse(grepl("_N", row.names(data)), "N", ifelse(grepl("_T", row.names(data)), "T", NA))
vars <- colnames(data)[!grepl("Pattern|GP|T_N", colnames(data))]
groups <- c("GP")
plot_list <- list()
for (group in groups){
  for (var in vars) 
  {
    data2 = data[,c(group,var,"T_N")]
    data2 <- na.omit(data2)
    colnames(data2)[1] <- "group"
    colnames(data2)[2] <- "value"
    p <- ggplot(data2, aes(x = group, y = value, fill = T_N)) + # when mapping, do not pass a variable name, otherwise the loop would keep the last variable of the for loop
      geom_boxplot(alpha = 1) +  # add the box plot
      geom_jitter(aes(group = T_N),alpha = 0.5,position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.75)) +  # two-level grouped points: the x axis is split by one discrete variable, and within each group another discrete variable splits it again (built into ggplot2 as position_jitterdodge)
      labs(x = "", y =  var, fill = NULL) +
      scale_fill_manual(values = c("#2874C5","#C6524A")) +
      theme_bw() +
      theme(axis.text=element_text(face="bold", size=14,color="black"),
            axis.title.y=element_text(face="bold", size=14,color ="black"),
            axis.text.x = element_text(),
            strip.text = element_text(size = 14,face = "bold"),
            strip.background = element_blank(),
            axis.title.x = element_blank(),
            #legend.key.size = unit(0.3, "cm"),
            legend.position = "top" ) + 
      stat_compare_means(aes(group = T_N),
                         label = "p.signif",
                         hide.ns = FALSE,
                         method = "wilcox"
                         #label.y = c(0.2)
      )
    plot_list[[var]] <- p
  }
  plots <- do.call(grid.arrange, c(plot_list, ncol = 4))
  filename <- paste("./2tumor+normal/", group, "_BoxPlot.pdf", sep = "")
  ggsave(filename, plots, width = 12, height = 28)
}

# Tumor
dir.create("./2tumor")
dataT <- data[grepl("_T", rownames(data)),]
# continuous variables
vars <- colnames(data)[!grepl("Pattern|GP|T_N", colnames(data))]
groups <- c("GP")
comparisons <- list(c("S","P"),c("S","MP"),c("P","MP"))
plot_list <- list()
for (group in groups){
  for (var in vars) 
  {
    data2 = dataT[,c(group,var)]
    data2 <- na.omit(data2)
    colnames(data2)[1] <- "group"
    colnames(data2)[2] <- "value"
    
    p <- ggplot(data2, aes(x = group, y = value, fill = group)) + # when mapping, do not pass a variable name, otherwise the loop would keep the last variable of the for loop
      geom_boxplot(alpha = 1) +  # add the box plot
      geom_jitter(width = 0.1, height = 0, alpha = 0.4) +  # add the scatter points
      labs(x = "", y =  var) +
      scale_fill_npg() +
      theme_bw() +
      theme(axis.text=element_text(face="bold", size=14,color="black"),
            axis.title.y=element_text(face="bold", size=14,color ="black"),
            axis.text.x = element_text(),
            strip.text = element_text(size = 14,face = "bold"),
            strip.background = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "none" ) + 
      stat_compare_means(label =  "p.signif",comparisons = comparisons,method = "wilcox", p.adjust.method = "BH")
    plot_list[[var]] <- p
  }
  plots <- do.call(grid.arrange, c(plot_list, ncol = 4))
  filename <- paste("./2tumor/", group, "_BoxPlot.pdf", sep = "")
  ggsave(filename, plots, width = 12, height = 28)
}

#Normal
dir.create("./2normal")
dataN <- data[grepl("_N", rownames(data)),]
# continuous variables
vars <- colnames(data)[!grepl("Pattern|GP|T_N", colnames(data))]
groups <- c("GP")
comparisons <- list(c("S","P"),c("S","MP"),c("P","MP"))
plot_list <- list()
for (group in groups){
  for (var in vars) 
  {
    data2 = dataN[,c(group,var)]
    data2 <- na.omit(data2)
    colnames(data2)[1] <- "group"
    colnames(data2)[2] <- "value"
    
    p <- ggplot(data2, aes(x = group, y = value, fill = group)) + # when mapping, do not pass a variable name, otherwise the loop would keep the last variable of the for loop
      geom_boxplot(alpha = 1) +  # add the box plot
      geom_jitter(width = 0.1, height = 0, alpha = 0.4) +  # add the scatter points
      labs(x = "", y =  var) +
      scale_fill_npg() +
      theme_bw() +
      theme(axis.text=element_text(face="bold", size=14,color="black"),
            axis.title.y=element_text(face="bold", size=14,color ="black"),
            axis.text.x = element_text(),
            strip.text = element_text(size = 14,face = "bold"),
            strip.background = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "none" ) + 
      stat_compare_means(label = "p.signif", comparisons = comparisons,
                         method = "wilcox", p.adjust.method = "BH")
    plot_list[[var]] <- p
  }
  plots <- do.call(grid.arrange, c(plot_list, ncol = 4))
  filename <- paste("./2normal/", group, "_BoxPlot.pdf", sep = "")
  ggsave(filename, plots, width = 12, height = 32)}
