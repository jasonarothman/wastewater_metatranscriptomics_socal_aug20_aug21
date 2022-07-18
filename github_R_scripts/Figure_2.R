library(ggplot2)
library(dplyr)
library(tidyverse)
library(rcartocolor)
library(RColorBrewer)
library(patchwork)
library(reshape2)

#DOMAINS AND UNKNOWNS PLOTTING

data<-read.table("domain_counts_plot_percents.txt", sep = '\t', quote = "", 
                 header = TRUE, row.names = 1)

var<-read.table("metadata_all_sewage_aug20_aug21.txt", header = TRUE, row.names = 1, sep ='\t')

melted<-melt(data)
for_plotting<-cbind(melted,var)

safe_pal<-carto_pal(11,"Safe")
p1<-ggplot(for_plotting, aes(x = reorder(as.factor(Date.1),Days_from_AUG_3_2020), value, fill = variable)) + 
  geom_bar(position = "fill",stat = "identity") + 
  scale_fill_manual(values = safe_pal)
p2<-p1 + theme_bw()
p3<-p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p4<-p3 + labs(x = "Sampling Date", y = "Relative Abundance", fill = "Taxonomic Group") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.text.x = element_text(angle = 90, h = 1, v = .5, color = "Black", size=6),
        axis.title = element_text(face = "bold", size = 12))
p5<-p4 + facet_grid(~Plant, 
                    #ncol = 4, 
                    scales = "free_x", space = "free")
p6<-p5 + theme(panel.spacing.x = unit(0.3, "lines"), panel.border = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold", size = 12),
               legend.title = element_text(face = "bold", size = 12),
               legend.text = element_text(face = "bold", size = 10),
               legend.position = "bottom")+
  scale_y_continuous(expand = c(0.001,0.001))
p6_domains<-p6


#AMR 
data<-read.table("top10_AMR_drug_classes.txt", header = TRUE, row.names = 1,
                 sep= '\t', check.names = FALSE)

var<-read.table("AMR_metadata_all_samples.txt", header = TRUE, row.names = 1, sep ='\t')


melted<-melt(data)
for_plotting<-cbind(melted,var)

safe_pal<-carto_pal(11,"Safe")
p1<-ggplot(for_plotting, aes(x = reorder(as.factor(Date.1),Days_from_AUG_3_2020), value, fill = variable)) + 
  geom_bar(position = "fill",stat = "identity") + 
  scale_fill_manual(values = safe_pal)
p2<-p1 + theme_bw()
p3<-p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p4<-p3 + labs(x = "Sampling Date", y = "Relative Abundance", fill = "Antibiotic Class") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.text.x = element_text(angle = 90, h = 1, v = .5, color = "Black", size=6),
        axis.title = element_text(face = "bold", size = 12))
p5<-p4 + facet_grid(~Plant, 
                    scales = "free_x", space = "free")
p6<-p5 + theme(panel.spacing.x = unit(0.3, "lines"), panel.border = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold", size = 12),
               legend.title = element_text(face = "bold", size = 12),
               legend.text = element_text(face = "bold", size = 10),
               legend.position = "bottom")+
  scale_y_continuous(expand = c(0.001,0.001))
p6_amr<-p6

#BACTERIA
data<-read.table("top_10_bacteria_family_relative_abundance.txt", header = TRUE, row.names = 1,
                 sep= '\t', check.names = FALSE)

melted<-melt(data)
for_plotting<-cbind(melted,var)

safe_pal<-carto_pal(11,"Safe")
p1<-ggplot(for_plotting, aes(x = reorder(as.factor(Date.1),Days_from_AUG_3_2020), value, fill = variable)) + 
  geom_bar(position = "fill",stat = "identity") + 
  scale_fill_manual(values = safe_pal)
p2<-p1 + theme_bw()
p3<-p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p4<-p3 + labs(x = "Sampling Date", y = "Relative Abundance", fill = "Bacterial Family") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.text.x = element_text(angle = 90, h = 1, v = .5, color = "Black", size=6),
        axis.title = element_text(face = "bold", size = 12))
p5<-p4 + facet_grid(~Plant, 
                    scales = "free_x", space = "free")
p6<-p5 + theme(panel.spacing.x = unit(0.3, "lines"), panel.border = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold", size = 12),
               legend.title = element_text(face = "bold", size = 12),
               legend.text = element_text(face = "bold", size = 10),
               legend.position = "bottom")+
  scale_y_continuous(expand = c(0.001,0.001))
p6_bacteria<-p6
p6_bacteria

#VIRUSES
data<-read.table("top_10_viruses_relative_abundance.txt", header = TRUE, row.names = 1,
                 sep= '\t', check.names = FALSE)

var<-read.table("metadata_all_sewage_aug20_aug21.txt", header = TRUE, row.names = 1, sep ='\t')

melted<-melt(data)
for_plotting<-cbind(melted,var)

safe_pal<-carto_pal(11,"Safe")
p1<-ggplot(for_plotting, aes(x = reorder(as.factor(Date.1),Days_from_AUG_3_2020), value, fill = variable)) + 
  geom_bar(position = "fill",stat = "identity") + 
  scale_fill_manual(values = safe_pal)
p2<-p1 + theme_bw()
p3<-p2 + theme(panel.grid.minor = element_blank())
p4<-p3 + labs(x = "Virus", y = "Relative Abundance", fill = "Virus") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.text.x = element_text(angle = 90, h = 1, v = .5, color = "Black", size=6),
        axis.title = element_text(face = "bold", size = 12))
p5<-p4 + facet_grid(~Plant, 
                    scales = "free_x", space = "free")
p6<-p5 + theme(panel.spacing.x = unit(0.3, "lines"), panel.border = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold", size = 12),
               legend.title = element_text(face = "bold", size = 12),
               legend.text = element_text(face = "bold", size = 10),
               legend.position = "bottom") +
  scale_y_continuous(expand = c(0.001,0.001))
p6_virus<-p6

pcombined<-p6_domains/p6_amr/p6_bacteria/p6_virus
pcombined + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold", color = "black"))
