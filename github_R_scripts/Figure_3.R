library(ggplot2)
library(dplyr)
library(tidyverse)
library(vegan)
library(rcartocolor)
library(RColorBrewer)
library(patchwork)
library(ggrepel)

#PLOTTING 8 PLANTS AUG-NOV
data<-read.table("AMR_counts_all_samples_aug20_aug21.txt", sep = '\t', quote = "", 
                 header = TRUE, dec = ".", row.names = 1)

var<-read.table("AMR_metadata_all_samples.txt", header = TRUE, row.names = 1, sep ='\t')

data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)

data_norm<-data_norm*1000000
data_norm<-as.data.frame(data_norm)

data_norm<-data_norm %>%
  filter(var$Month == "AUG_20" | var$Month == "SEP_20" | var$Month == "OCT_20" | 
           var$Month == "NOV_20")
data_norm<-as.matrix(data_norm)

var<-var %>%
  filter(var$Month == "AUG_20" | var$Month == "SEP_20" | var$Month == "OCT_20" | 
           var$Month == "NOV_20")


data_nmds <- metaMDS(data_norm, distance = 'bray', autotransform = F, 
                     k = 2, noshare = F, trymax = 100, parallel = 6)
nmds_scores<-scores(data_nmds)
nmds_df<-as.data.frame(nmds_scores)

attach(var)
for_plotting<-merge(nmds_df,var, by = "row.names")

p1<-ggplot(for_plotting, aes(NMDS1, NMDS2, color = Plant))+ 
  scale_color_brewer(palette = "Dark2")+
  geom_label_repel(aes(label = Month, fontface = "bold"), max.overlaps = Inf)


p2<-p1 + theme_bw()
p3<-p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p4<-p3
p5_mds<-p4 + theme(axis.title.x = element_text(face = "bold", size = 14),
                   legend.text = element_text(face = "bold", size = 12),
                   axis.title.y = element_text(face = "bold", size = 14), 
                   legend.title  = element_text(face = "bold", size = 14))

p5_mds_8plants<-p5_mds

#PLOTTING 3 PLANTS AUG20-AUG21
data<-read.table("AMR_counts_all_samples_aug20_aug21.txt", sep = '\t', quote = "", 
                 header = TRUE, row.names = 1)

var<-read.table("AMR_metadata_all_samples.txt", header = TRUE, row.names = 1, sep ='\t')

data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)

data_norm<-data_norm*1000000

data_norm<-as.data.frame(data_norm)
data_norm<-data_norm %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")
data_norm<-as.matrix(data_norm)

var<-var %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")

data_nmds <- metaMDS(data_norm, distance = 'bray', autotransform = F, 
                     k = 2, noshare = F, trymax = 100, parallel = 6)
nmds_scores<-scores(data_nmds)
nmds_df<-as.data.frame(nmds_scores)

attach(var)
for_plotting<-merge(nmds_df,var, by = "row.names")

p1<-ggplot(for_plotting, aes(NMDS1, NMDS2, color = Plant)) + 
  scale_color_brewer(palette = "Dark2") +
  geom_label_repel(aes(label = Month, fontface = "bold"), max.overlaps = Inf)

p2<-p1 + theme_bw()
p3<-p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p4<-p3
p5_mds_long<-p4 + theme(axis.title.x = element_text(face = "bold", size = 14),
                        legend.text = element_text(face = "bold", size = 12),
                        axis.title.y = element_text(face = "bold", size = 14), 
                        legend.title  = element_text(face = "bold", size = 14))
p5_mds_long

#AMR HEATMAPS
source("ancom_v2.1.R")
library(pheatmap)
library(ggplotify)

data<-read.table("AMR_counts_ARO-name_no-rrna.txt", 
                 header = TRUE, sep = '\t',row.names = 1, quote = "", dec = ".")

var<-read.table("AMR_metadata_all_samples.txt", header = TRUE,
                sep = '\t')
rownames(data)<-as.character(rownames(data))
data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)

data_norm<-as.data.frame(data_norm)
data_norm<-data_norm %>%
  filter(var$Month == "AUG_20" | var$Month == "SEP_20" | var$Month == "OCT_20" | 
           var$Month == "NOV_20")
data_norm<-as.matrix(data_norm)


data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)
data_norm$mean
data_filtered<-data %>%
  filter(data_norm$mean > 0.0001) 

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

data_filtered<-data_filtered %>%
  filter(var$Month == "AUG_20" | var$Month == "SEP_20" | var$Month == "OCT_20" | 
           var$Month == "NOV_20")
data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

var<-var %>%
  filter(var$Month == "AUG_20" | var$Month == "SEP_20" | var$Month == "OCT_20" | 
           var$Month == "NOV_20")

prepro<-feature_table_pre_process(feature_table = data_filtered,meta_data = var,sample_var = "Sample", group_var = NULL,
                                  out_cut = 0.05, zero_cut = 0.50, lib_cut = 1, neg_lb = FALSE)

feature_table = prepro$feature_table 
meta_data = prepro$meta_data 
struc_zero = prepro$structure_zeros 

ancom_results<-read.table("ancom_results_AMR_no_rRNA.txt", header = TRUE, sep = '\t',row.names = 1, quote = "", dec = ".")

ancom_results
rownames(ancom_results)<-ancom_results$taxa_id
heatmap_data<-feature_table %>%
  filter(ancom_results$detected_0.8 == "TRUE")

heatmap_data<-as.data.frame(t(as.matrix(heatmap_data)))

heatmap_data<- heatmap_data %>%
  arrange(var$Plant,var$Days_from_AUG_3_2020)

heatmap_data<-t(as.matrix(heatmap_data))
heatmap_data<-heatmap_data+1

newnames_row <- lapply(
  rownames(heatmap_data),
  function(x) bquote(bold(.(x))))
newnames_col <- lapply(
  colnames(heatmap_data),
  function(x) bquote(bold(.(x))))

no_rrna_heatmap<-pheatmap(log10(heatmap_data), cluster_rows = TRUE, cluster_cols = FALSE, 
                          gaps_col = c(16,31,44,46,60,68,75),
                          border_color = "black", 
                          cellheight = 12, cellwidth = 6,
                          fontsize_col = 6,
                          fontsize_row = 10,
                          drop_levels = TRUE,
                          legend = TRUE,
                          main = "Log10 of Differentially Abundant AMR Genes",
                          labels_row = as.expression(newnames_row),
                          labels_col = as.expression(newnames_col))
no_rrna_heatmap<-as.ggplot(no_rrna_heatmap)

(p5_mds_8plants+p5_mds_long)/no_rrna_heatmap

source("ancom_v2.1.R")
library(pheatmap)
data<-read.table("AMR_counts_ARO-name.txt", 
                 header = TRUE, sep = '\t',row.names = 1, quote = "", dec = ".")

var<-read.table("AMR_metadata_all_samples.txt", header = TRUE,
                sep = '\t')
rownames(data)<-as.character(rownames(data))
data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)

data_norm<-as.data.frame(data_norm)
data_norm<-data_norm %>%
  filter(var$Month == "AUG_20" | var$Month == "SEP_20" | var$Month == "OCT_20" | 
           var$Month == "NOV_20")
data_norm<-as.matrix(data_norm)


data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)
data_norm$mean
data_filtered<-data %>%
  filter(data_norm$mean > 0.0001) #%>%

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

data_filtered<-data_filtered %>%
  filter(var$Month == "AUG_20" | var$Month == "SEP_20" | var$Month == "OCT_20" | 
           var$Month == "NOV_20")
data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

var<-var %>%
  filter(var$Month == "AUG_20" | var$Month == "SEPT_20" | var$Month == "OCT_20" | 
           var$Month == "NOV_20")

prepro<-feature_table_pre_process(feature_table = data_filtered,meta_data = var,sample_var = "Sample", group_var = NULL,
                                  out_cut = 0.05, zero_cut = 0.50, lib_cut = 1, neg_lb = FALSE)

feature_table = prepro$feature_table 
meta_data = prepro$meta_data 
struc_zero = prepro$structure_zeros 

ancom_results<-read.table("ancom_results_AMR.txt", header = TRUE, sep = '\t',row.names = 1, quote = "", dec = ".")

ancom_results
rownames(ancom_results)<-ancom_results$taxa_id
heatmap_data<-feature_table %>%
  filter(ancom_results$detected_0.8 == "TRUE")

heatmap_data<-as.data.frame(t(as.matrix(heatmap_data)))

heatmap_data<- heatmap_data %>%
  arrange(var$Plant,var$Days_from_AUG_3_2020)

heatmap_data<-t(as.matrix(heatmap_data))
heatmap_data<-heatmap_data+1

newnames_row <- lapply(
  rownames(heatmap_data),
  function(x) bquote(bold(.(x))))
newnames_col <- lapply(
  colnames(heatmap_data),
  function(x) bquote(bold(.(x))))

rrna_ancom_heatmap<-pheatmap(log10(heatmap_data), cluster_rows = TRUE, cluster_cols = FALSE, 
                             gaps_col = c(13,23,33,34,43,49,54),
                             border_color = "black", 
                             cellheight = 12, cellwidth = 6,
                             fontsize_col = 6,
                             fontsize_row = 10,
                             drop_levels = TRUE,
                             legend = TRUE,
                             main = "Log10 of Differentially Abundant rRNA Genes",
                             labels_row = as.expression(newnames_row),
                             labels_col = as.expression(newnames_col))
rrna_ancom_heatmap<-as.ggplot(rrna_ancom_heatmap)

(p5_mds_8plants+p5_mds_long)/(no_rrna_heatmap+rrna_ancom_heatmap)

#AMR Maaslin plots
data<-read.table("rrna_results.txt", header = TRUE, sep = '\t',
                 check.names = FALSE)
data<-data %>%
  filter(qval < 0.05)

p1<-ggplot(data, aes(reorder(feature,coef), coef,
                     fill = coef > 0)) +
  geom_bar(stat = "identity", color = "Black", size = 1) + 
  coord_flip(ylim = c(-0.6,1.3)) + 
  scale_fill_manual(values = c("Darkred","Darkgreen"))
p1
p2<-p1 + theme_bw()
p3<-p2 + theme(panel.grid.major.x = element_line(size = 0.25,
                                                 color = "black"),
               panel.grid.minor = element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               axis.ticks.y = element_blank()) +
  geom_hline(yintercept = 0:0.3, size = 1)
p3
p4<- p3 + theme(axis.text = element_text(size = 10, face = "bold",
                                         color = "black"),
                legend.position = "none",
                axis.title = element_text(size = 16, face = "bold",
                                          color = "black"))+
  labs(x = "AMR rRNA Mutation", y = "Linear Model Coefficient")
p4
rrna_maslin<-p4
amr_rrna_maslin<-as.ggplot(rrna_maslin)

data<-read.table("significant_results.tsv", header = TRUE, sep = '\t',
                 check.names = FALSE)
data<-data %>%
  filter(qval < 0.05)

p1<-ggplot(data, aes(reorder(feature,coef), coef,
                     fill = coef > 0)) +
  geom_bar(stat = "identity", color = "Black", size = 1) + 
  coord_flip(ylim = c(-0.3,1.85)) + 
  scale_fill_manual(values = c("Darkred","Darkgreen"))
p1
p2<-p1 + theme_bw()
p3<-p2 + theme(panel.grid.major.x = element_line(size = 0.25,
                                                 color = "black"),
               panel.grid.minor = element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               axis.ticks.y = element_blank()) +
  geom_hline(yintercept = 0:0.3, size = 1)
p3
p4<- p3 + theme(axis.text = element_text(size = 10, face = "bold",
                                         color = "black"),
                legend.position = "none",
                axis.title = element_text(size = 16, face = "bold",
                                          color = "black"))+
  labs(x = "AMR Gene", y = "Linear Model Coefficient")
p4
amr_nonrrna_maslin<-p4
amr_nonrrna_maslin<-as.ggplot(amr_nonrrna_maslin)

pcombined<-(p5_mds_8plants+p5_mds_long)/(no_rrna_heatmap+rrna_ancom_heatmap)/(amr_nonrrna_maslin+amr_rrna_maslin)
pcombined + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold", color = "black"))
