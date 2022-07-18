library(ggplot2)
library(dplyr)
library(tidyverse)
library(rcartocolor)
library(RColorBrewer)
library(pals)
library(ggsci)
library(lubridate)

var<-read.table("metadata_all_sewage_aug20_aug21.txt", header = TRUE, row.names = 1, sep ='\t')

var$Date_ymd<-ymd(var$Date_ymd)

p1<-ggplot(data = var,
           aes( x = var$Date_ymd, y = Plant, color = Plant, group = Plant)) +
  geom_point(color= "black", size = 6, shape = "|")+
  geom_point(aes(color=Plant), size=6, shape = "|") + 
  scale_color_brewer(palette = "Dark2")

p1
p2<-p1 + theme_bw(base_size = 18)
p2
p3<-p2 + labs(x = "Sampling Timepoint", 
              y = "Plant",
              color = "Plant") +
  theme(axis.text.x = element_text(angle = 0, v = 0.5, color = "Black", 
                                   face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", color = "Black", size = 16),
        axis.title = element_text(face = "bold", color = "Black", size = 20),
        legend.position = "none",
        panel.grid.minor.x = element_blank())
p3
p4<-p3 + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", color = "Black", size = 12),
        axis.text.x = element_text(angle = 0, size = 16))
p4

p5<-p4 + coord_fixed(ratio = 10) + theme(plot.margin = unit(c(0,0.5,0,0.5), "cm"))+
  scale_y_discrete(limits = rev) +
  scale_x_date(date_breaks = "month", date_labels = "%b%g",expand =c(0.02,0.02))

p5