library(ggplot2)
library(dplyr)
library(tidyverse)
library(vegan)
library(rcartocolor)
library(RColorBrewer)
library(lmerTest)
library(patchwork)
library(ggrepel)
library(ggplotify)
library(pheatmap)

#PLOTTING 8 PLANTS AUG-NOV
data<-read.table("all_samples_humann3_pathabundance.txt", sep = '\t', quote = "", 
                 header = TRUE, row.names = 1)

var<-read.table("all_samples_humann3_pathabundance_metadata.txt", header = TRUE, row.names = 1, sep ='\t')

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
data_norm<-as.data.frame(data_norm)
data_norm<-as.matrix(data_norm)

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
p5_mds_8plants

#PLOTTING 3 PLANTS AUG20-AUG21
data<-read.table("all_samples_humann3_pathabundance.txt", sep = '\t', quote = "", 
                 header = TRUE, row.names = 1)

var<-read.table("all_samples_humann3_pathabundance_metadata.txt", header = TRUE, row.names = 1, sep ='\t')

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


#PATHWAY Maaslin plots
data<-read.table("significant_results.tsv", header = TRUE, sep = '\t')
data<-data %>%
  filter(qval < 0.05)

p1<-ggplot(data, aes(reorder(feature,coef), coef,
                     fill = coef > 0)) +
  geom_bar(stat = "identity", color = "Black", size = 1) + 
  coord_flip(ylim = c(-1.1,1.3)) + 
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

p4<- p3 + theme(axis.text = element_text(size = 10, face = "bold",
                                         color = "black"),
                legend.position = "none",
                axis.title = element_text(size = 16, face = "bold",
                                          color = "black"))+
  labs(x = "Metabolic Pathway", y = "Linear Model Coefficient")
p4
pathway_maslin<-p4
pathway_maslin<-as.ggplot(pathway_maslin)

pcombined<-(p5_mds_8plants+p5_mds_long)/(pathway_maslin)
pcombined + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold", color = "black"))
