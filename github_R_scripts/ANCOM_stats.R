source("ancom_v2.1.R")
library(dplyr)
library(tidyverse)

###AMR### 
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
  filter(data_norm$mean > 0.0001) #%>%

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

ancom_results<-ANCOM(feature_table = feature_table, meta_data = meta_data,
                    struc_zero = struc_zero, main_var = "Plant",
                   p_adj_method = "BH", alpha = 0.05, adj_formula = "Month",
                  rand_formula = "~1 | Batch")
ancom_results

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
  filter(var$Month == "AUG_20" | var$Month == "SEP_20" | var$Month == "OCT_20" | 
           var$Month == "NOV_20")

prepro<-feature_table_pre_process(feature_table = data_filtered,meta_data = var,sample_var = "Sample", group_var = NULL,
                                  out_cut = 0.05, zero_cut = 0.50, lib_cut = 1, neg_lb = FALSE)

feature_table = prepro$feature_table 
meta_data = prepro$meta_data 
struc_zero = prepro$structure_zeros 

ancom_results<-ANCOM(feature_table = feature_table, meta_data = meta_data,
                    struc_zero = struc_zero, main_var = "Plant",
                   p_adj_method = "BH", alpha = 0.05, adj_formula = "Month",
                  rand_formula = "~1 | Batch")
ancom_results

###BACTERIA###
data<-read.table("all_samples_combined_bracken_bacteria_genus_counts_aug20_aug21.txt", 
                 header = TRUE, sep = '\t',row.names = 1, quote = "", dec = ".")

var<-read.table("metadata_all_sewage_aug20_aug21.txt", header = TRUE,
                sep = '\t')

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

data_filtered<-data %>%
  filter(data_norm$mean > 0.0001) #%>%

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

ancom_results<-ANCOM(feature_table = feature_table, meta_data = meta_data,
                    struc_zero = struc_zero, main_var = "Plant",
                   p_adj_method = "BH", alpha = 0.05, adj_formula = "Month",
                  rand_formula = "~1 | Batch")
ancom_results

###VIRUSES###
data<-read.table("all_viruses_counts_table_aug20_aug21.txt", 
                 header = TRUE, sep = '\t',row.names = 1, quote = "", dec = ".")

var<-read.table("metadata_all_sewage_aug20_aug21.txt", header = TRUE,
                sep = '\t')

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

data_filtered<-data %>%
  filter(data_norm$mean > 0.0001) #%>%

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

ancom_results<-ANCOM(feature_table = feature_table, meta_data = meta_data,
                    struc_zero = struc_zero, main_var = "Plant",
                   p_adj_method = "BH", alpha = 0.05, adj_formula = "Month",
                  rand_formula = "~1 | Batch")
ancom_results