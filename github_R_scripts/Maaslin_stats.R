library(Maaslin2)
library(dplyr)

###VIRUSES###
data<-read.table("all_viruses_counts_table_aug20_aug21.txt", 
                 header = TRUE, sep = '\t',row.names = 1, quote = "", dec = ".",
                 stringsAsFactors = FALSE)

var<-read.table("metadata_all_sewage_aug20_aug21.txt", header = TRUE,
                sep = '\t', row.names = 1, stringsAsFactors = FALSE)

data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)

data_norm<-as.data.frame(data_norm)
data_norm<-data_norm %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")

data_norm<-as.matrix(data_norm)


data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)

data_filtered<-data %>%
  filter(data_norm$mean > 0.0001) #%>%

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

data_filtered<-data_filtered %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")
data_filtered<-as.data.frame(t(as.matrix(data_filtered)))
var<-var %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")

fit_data <- Maaslin2(input_data = data_filtered,
                     input_metadata = var,
                     output = "ESC_HTP_PL_longitudinal_viruses",
                     min_prevalence = 0.50,
                     random_effects = c("Plant","Batch"),
                     fixed_effects = c("Days_from_AUG_3_2020"),
                     max_significance = 0.05)

###BACTERIAL###
data<-read.table("all_samples_combined_bracken_bacteria_genus_counts_aug20_aug21.txt", 
                 header = TRUE, sep = '\t',row.names = 1, quote = "", dec = ".",
                 stringsAsFactors = FALSE)

var<-read.table("metadata_all_sewage_aug20_aug21.txt", header = TRUE,
                sep = '\t', row.names = 1, stringsAsFactors = FALSE)

data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)

data_norm<-as.data.frame(data_norm)
data_norm<-data_norm %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")

data_norm<-as.matrix(data_norm)

data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)

data_filtered<-data %>%
  filter(data_norm$mean > 0.0001) #%>%

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

data_filtered<-data_filtered %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

var<-var %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")

fit_data <- Maaslin2(input_data = data_filtered,
                     input_metadata = var,
                     output = "ESC_HTP_PL_longitudinal_bacterial_genus",
                     min_prevalence = 0.50,
                     random_effects = c("Plant","Batch"),
                     fixed_effects = c("Days_from_AUG_3_2020"),
                     max_significance = 0.05)

###HUMANN3 PATHWAY ABUNDANCES###
data<-read.table("all_samples_humann3_pathabundance.txt", 
                 header = TRUE, sep = '\t',row.names = 1, quote = "", dec = ".",
                 stringsAsFactors = FALSE)

var<-read.table("all_samples_humann3_pathabundance_metadata.txt", header = TRUE,
                sep = '\t', row.names = 1, stringsAsFactors = FALSE)

rownames(data)<-as.character(rownames(data))

data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)

data_norm<-as.data.frame(data_norm)
data_norm<-data_norm %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")

data_norm<-as.matrix(data_norm)

data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)

data_filtered<-data %>%
  filter(data_norm$mean > 0.0001) #%>%

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

data_filtered<-data_filtered %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

var<-var %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")

fit_data <- Maaslin2(input_data = data_filtered,
                     input_metadata = var,
                     output = "ESC_HTP_PL_longitudinal_humann3_pathabundance",
                     min_prevalence = 0.50,
                     random_effects = c("Plant","Batch"),
                     fixed_effects = c("Days_from_AUG_3_2020"),
                     max_significance = 0.05)

###AMR###

data<-read.table("AMR_counts_ARO-name.txt", 
                 header = TRUE, sep = '\t',row.names = 1, quote = "", dec = ".",
                 stringsAsFactors = FALSE)

var<-read.table("AMR_metadata_all_samples.txt", header = TRUE,
                sep = '\t', row.names = 1, stringsAsFactors = FALSE)

rownames(data)<-as.character(rownames(data))

data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)

data_norm<-as.data.frame(data_norm)
data_norm<-data_norm %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")

data_norm<-as.matrix(data_norm)

data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)

data_filtered<-data %>%
  filter(data_norm$mean > 0.0001) #%>%

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

data_filtered<-data_filtered %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

var<-var %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")

fit_data <- Maaslin2(input_data = data_filtered,
                     input_metadata = var,
                     output = "ESC_HTP_PL_longitudinal_AMR",
                     min_prevalence = 0.50,
                     random_effects = c("Plant","Batch"),
                     fixed_effects = c("Days_from_AUG_3_2020"),
                     max_significance = 0.05)

data<-read.table("AMR_counts_ARO-name_no-rrna.txt", 
                 header = TRUE, sep = '\t',row.names = 1, quote = "", dec = ".",
                 stringsAsFactors = FALSE)

var<-read.table("AMR_metadata_all_samples.txt", header = TRUE,
                sep = '\t', row.names = 1, stringsAsFactors = FALSE)

rownames(data)<-as.character(rownames(data))

data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)

data_norm<-as.data.frame(data_norm)
data_norm<-data_norm %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")

data_norm<-as.matrix(data_norm)

data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)

data_filtered<-data %>%
  filter(data_norm$mean > 0.0001) #%>%

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

data_filtered<-data_filtered %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

var<-var %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")

fit_data <- Maaslin2(input_data = data_filtered,
                     input_metadata = var,
                     output = "ESC_HTP_PL_longitudinal_AMR_no_rRNA",
                     min_prevalence = 0.50,
                     random_effects = c("Plant","Batch"),
                     fixed_effects = c("Days_from_AUG_3_2020"),
                     max_significance = 0.05)
