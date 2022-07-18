library(dplyr)
library(tidyverse)
library(vegan)
library(dunn.test)
library(lmerTest)
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
bray<- vegdist(data_norm, method = 'bray')

ad<-adonis(data_norm ~ Plant * Month + Batch,  
           data = var, parallel = 4, 
           method = "bray")

shannon<-diversity(data_norm, index = "shannon")

for_plotting<-merge(shannon,var, by = "row.names")

shapiro.test(shannon)
hist(shannon)
qqnorm(shannon)

shannon_krusk<-kruskal.test(shannon~Plant, data = for_plotting)
shannon_krusk
dunn.test(for_plotting$shannon, for_plotting$Plant, method = "bh")

shannon_lmer<-lmer(shannon ~ as.numeric(as.factor(Days_from_AUG_3_2020)) + 
                     (1|Plant) + (1|Batch), data = for_plotting)
summary(shannon_lmer)  

bray.df<-as.data.frame(as.matrix(bray))
bray_means<-as.data.frame(((rowMeans(bray.df))))
names(bray_means)[1]<-"Dissimilarity"
bray_means$Similarity<-(1-bray_means$Dissimilarity)

for_plotting<-merge(var,bray_means, by = "row.names")

bray_lmer<-lmer(Similarity ~ as.numeric(as.factor(Days_from_AUG_3_2020)) + 
                  (1|Plant) + (1|Batch), data = for_plotting)
summary(bray_lmer)

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

bray<- vegdist(data_norm, method = 'bray')

ad<-adonis(data_norm ~ Plant * Month + Batch,  
           data = var, parallel = 4, 
           method = "bray")
ad

shannon<-diversity(data_norm, index = "shannon")

for_plotting<-merge(var,shannon, by = "row.names")

shapiro.test(shannon)
hist(shannon)
qqnorm(shannon)

shannon_krusk<-kruskal.test(shannon~Plant, data = for_plotting)
shannon_krusk
dunn.test(for_plotting$shannon, for_plotting$Plant, method = "bh")
shannon_lmer<-lmer(shannon ~ as.numeric(as.factor(Days_from_AUG_3_2020)) + 
                     (1|Plant) + (1|Batch), data = for_plotting)
summary(shannon_lmer)  

bray.df<-as.data.frame(as.matrix(bray))
bray_means<-as.data.frame(((rowMeans(bray.df))))
names(bray_means)[1]<-"Dissimilarity"
bray_means$Similarity<-(1-bray_means$Dissimilarity)

for_plotting<-cbind(var,bray_means, by = "row.names")


bray_lmer<-lmer(Similarity ~ as.numeric(as.factor(Days_from_AUG_3_2020)) + 
                  (1|Plant) + (1|Batch), data = for_plotting)
summary(bray_lmer)

###BACTERIA###
data<-read.table("all_samples_combined_bracken_bacteria_species_counts_aug20_aug21.txt", sep = '\t', quote = "", 
                 header = TRUE, row.names = 1)

var<-read.table("metadata_all_sewage_aug20_aug21.txt", header = TRUE, row.names = 1, sep ='\t')

data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)

data_norm<-as.data.frame(data_norm)
data_norm<-data_norm %>%
  filter(var$Month == "AUG_20" | var$Month == "SEP_20" | var$Month == "OCT_20" | 
           var$Month == "NOV_20")
data_norm<-as.matrix(data_norm)

var<-var %>%
  filter(var$Month == "AUG_20" | var$Month == "SEP_20" | var$Month == "OCT_20" | 
           var$Month == "NOV_20")

data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)
data_norm <-data_norm %>%
  rownames_to_column("Sample") %>%
  filter(data_norm$mean > 0.0001) %>%
  subset(., select = -c(mean))

data_norm = setNames(data.frame(t(data_norm[,-1])),data_norm[,1])

data_nmds <- metaMDS(data_norm, distance = 'bray', autotransform = F, 
                     k = 2, noshare = F, trymax = 100, parallel = 6)
nmds_scores<-scores(data_nmds)
nmds_df<-as.data.frame(nmds_scores)

attach(var)
for_plotting<-merge(nmds_df, var, by = "row.names")

bray<- vegdist(data_norm, method = 'bray')

ad<-adonis(data_norm ~ Plant * Month + Batch,
           data = var, parallel = 4, 
           method = "bray")
ad

shannon<-diversity(data_norm, index = "shannon")

for_plotting<-merge(var,shannon, by = "row.names")

shapiro.test(shannon)
hist(shannon)
qqnorm(shannon)

shannon_krusk<-kruskal.test(shannon~Plant, data = for_plotting)
shannon_krusk
dunn.test(for_plotting$shannon, for_plotting$Plant, method = "bh")

shannon_lmer<-lmer(shannon ~ as.numeric(as.factor(Days_from_AUG_3_2020)) + 
                     (1|Plant) + (1|Batch), data = for_plotting)
summary(shannon_lmer) 

bray.df<-as.data.frame(as.matrix(bray))
bray_means<-as.data.frame(((rowMeans(bray.df))))
names(bray_means)[1]<-"Dissimilarity"
bray_means$Similarity<-(1-bray_means$Dissimilarity)

for_plotting<-merge(var,bray_means, by = "row.names")

 
bray_lmer<-lmer(Similarity ~ as.numeric(as.factor(Days_from_AUG_3_2020)) + 
                  (1|Plant) + (1|Batch), data = for_plotting)
summary(bray_lmer)


data<-read.table("all_samples_combined_bracken_bacteria_species_counts_aug20_aug21.txt", sep = '\t', quote = "", 
                 header = TRUE, row.names = 1)

var<-read.table("metadata_all_sewage_aug20_aug21.txt", header = TRUE, row.names = 1, sep ='\t')

data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)

data_norm<-as.data.frame(data_norm)
data_norm<-data_norm %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")
data_norm<-as.matrix(data_norm)

var<-var %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")

data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)
data_norm <-data_norm %>%
  rownames_to_column("Sample") %>%
  filter(data_norm$mean > 0.0001) %>%
  subset(., select = -c(mean))

data_norm = setNames(data.frame(t(data_norm[,-1])),data_norm[,1])

data_nmds <- metaMDS(data_norm, distance = 'bray', autotransform = F, 
                     k = 2, noshare = F, trymax = 100, parallel = 6)
nmds_scores<-scores(data_nmds)
nmds_df<-as.data.frame(nmds_scores)

attach(var)
for_plotting<-merge(nmds_df, var, by = "row.names")

bray<- vegdist(data_norm, method = 'bray')

ad<-adonis(data_norm ~ Plant * Month + Batch,  
           data = var, parallel = 4, 
           method = "bray")
ad

shannon<-diversity(data_norm, index = "shannon")

for_plotting<-merge(var,shannon, by = "row.names")

shapiro.test(shannon)
hist(shannon)
qqnorm(shannon)

shannon_krusk<-kruskal.test(shannon~Plant, data = for_plotting)
shannon_krusk
dunn.test(for_plotting$shannon, for_plotting$Plant, method = "bh")

shannon_lmer<-lmer(shannon ~ as.numeric(as.factor(Days_from_AUG_3_2020)) + 
                     (1|Plant) + (1|Batch), data = for_plotting)
summary(shannon_lmer) 

bray.df<-as.data.frame(as.matrix(bray))
bray_means<-as.data.frame(((rowMeans(bray.df))))
names(bray_means)[1]<-"Dissimilarity"
bray_means$Similarity<-(1-bray_means$Dissimilarity)

for_plotting<-merge(var,bray_means, by = "row.names")

bray_lmer<-lmer(Similarity ~ as.numeric(as.factor(Days_from_AUG_3_2020)) + 
                  (1|Plant) + (1|Batch), data = for_plotting)
summary(bray_lmer)

###VIRUSES###
data<-read.table("all_viruses_counts_table_aug20_aug21.txt", sep = '\t', quote = "", 
                 header = TRUE, row.names = 1)

var<-read.table("metadata_all_sewage_aug20_aug21.txt", header = TRUE, row.names = 1, sep ='\t')

data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)

data_norm<-as.data.frame(data_norm)
data_norm<-data_norm %>%
  filter(var$Month == "AUG_20" | var$Month == "SEP_20" | var$Month == "OCT_20" | 
           var$Month == "NOV_20")
data_norm<-as.matrix(data_norm)

var<-var %>%
  filter(var$Month == "AUG_20" | var$Month == "SEP_20" | var$Month == "OCT_20" | 
           var$Month == "NOV_20")

data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)
data_norm <-data_norm %>%
  rownames_to_column("Sample") %>%
  filter(data_norm$mean > 0.0001) %>%
  subset(., select = -c(mean))

data_norm = setNames(data.frame(t(data_norm[,-1])),data_norm[,1])

data_nmds <- metaMDS(data_norm, distance = 'bray', autotransform = F, 
                     k = 2, noshare = F, trymax = 100, parallel = 6)
nmds_scores<-scores(data_nmds)
nmds_df<-as.data.frame(nmds_scores)

attach(var)
for_plotting<-merge(nmds_df, var, by = "row.names")

bray<- vegdist(data_norm, method = 'bray')

ad<-adonis(data_norm ~ Plant * Month + Batch,
           data = var, parallel = 4, 
           method = "bray")
ad

shannon<-diversity(data_norm, index = "shannon")

for_plotting<-merge(var,shannon, by = "row.names")

shapiro.test(shannon)
hist(shannon)
qqnorm(shannon)

shannon_krusk<-kruskal.test(shannon~Plant, data = for_plotting)
shannon_krusk
dunn.test(for_plotting$shannon, for_plotting$Plant, method = "bh")

shannon_lmer<-lmer(shannon ~ as.numeric(as.factor(Days_from_AUG_3_2020)) + 
                     (1|Plant) + (1|Batch), data = for_plotting)
summary(shannon_lmer) 

bray.df<-as.data.frame(as.matrix(bray))
bray_means<-as.data.frame(((rowMeans(bray.df))))
names(bray_means)[1]<-"Dissimilarity"
bray_means$Similarity<-(1-bray_means$Dissimilarity)

for_plotting<-merge(var,bray_means, by = "row.names")


bray_lmer<-lmer(Similarity ~ as.numeric(as.factor(Days_from_AUG_3_2020)) + 
                  (1|Plant) + (1|Batch), data = for_plotting)
summary(bray_lmer)


data<-read.table("all_viruses_counts_table_aug20_aug21.txt", sep = '\t', quote = "", 
                 header = TRUE, row.names = 1)

var<-read.table("metadata_all_sewage_aug20_aug21.txt", header = TRUE, row.names = 1, sep ='\t')

data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)

data_norm<-as.data.frame(data_norm)
data_norm<-data_norm %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")
data_norm<-as.matrix(data_norm)

var<-var %>%
  filter(var$Plant == "ESC" | var$Plant == "HTP" | var$Plant == "PL")

data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)
data_norm <-data_norm %>%
  rownames_to_column("Sample") %>%
  filter(data_norm$mean > 0.0001) %>%
  subset(., select = -c(mean))

data_norm = setNames(data.frame(t(data_norm[,-1])),data_norm[,1])

data_nmds <- metaMDS(data_norm, distance = 'bray', autotransform = F, 
                     k = 2, noshare = F, trymax = 100, parallel = 6)
nmds_scores<-scores(data_nmds)
nmds_df<-as.data.frame(nmds_scores)

attach(var)
for_plotting<-merge(nmds_df, var, by = "row.names")

bray<- vegdist(data_norm, method = 'bray')

ad<-adonis(data_norm ~ Plant * Month + Batch,  
           data = var, parallel = 4, 
           method = "bray")
ad

shannon<-diversity(data_norm, index = "shannon")

for_plotting<-merge(var,shannon, by = "row.names")

shapiro.test(shannon)
hist(shannon)
qqnorm(shannon)

shannon_krusk<-kruskal.test(shannon~Plant, data = for_plotting)
shannon_krusk
dunn.test(for_plotting$shannon, for_plotting$Plant, method = "bh")

shannon_lmer<-lmer(shannon ~ as.numeric(as.factor(Days_from_AUG_3_2020)) + 
                     (1|Plant) + (1|Batch), data = for_plotting)
summary(shannon_lmer) 

bray.df<-as.data.frame(as.matrix(bray))
bray_means<-as.data.frame(((rowMeans(bray.df))))
names(bray_means)[1]<-"Dissimilarity"
bray_means$Similarity<-(1-bray_means$Dissimilarity)

for_plotting<-merge(var,bray_means, by = "row.names")

bray_lmer<-lmer(Similarity ~ as.numeric(as.factor(Days_from_AUG_3_2020)) + 
                  (1|Plant) + (1|Batch), data = for_plotting)
summary(bray_lmer)

###HUMANN3 PATHWAYS###
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


data_nmds <- metaMDS(data_norm, distance = 'bray', autotransform = F, 
                     k = 2, noshare = F, trymax = 100, parallel = 6)
nmds_scores<-scores(data_nmds)
nmds_df<-as.data.frame(nmds_scores)

attach(var)
bray<- vegdist(data_norm, method = 'bray')

ad<-adonis(data_norm ~ Plant * Month + Batch,  
           data = var, parallel = 4, 
           method = "bray")

shannon<-diversity(data_norm, index = "shannon")

for_plotting<-merge(shannon,var, by = "row.names")

shapiro.test(shannon)
hist(shannon)
qqnorm(shannon)

shannon_krusk<-kruskal.test(shannon~Plant, data = for_plotting)
shannon_krusk
dunn.test(for_plotting$shannon, for_plotting$Plant, method = "bh")

shannon_lmer<-lmer(shannon ~ as.numeric(as.factor(Days_from_AUG_3_2020)) + 
                     (1|Plant) + (1|Batch), data = for_plotting)
summary(shannon_lmer)  

bray.df<-as.data.frame(as.matrix(bray))
bray_means<-as.data.frame(((rowMeans(bray.df))))
names(bray_means)[1]<-"Dissimilarity"
bray_means$Similarity<-(1-bray_means$Dissimilarity)

for_plotting<-merge(var,bray_means, by = "row.names")

bray_lmer<-lmer(Similarity ~ as.numeric(as.factor(Days_from_AUG_3_2020)) + 
                  (1|Plant) + (1|Batch), data = for_plotting)
summary(bray_lmer)

data<-read.table("all_samples_humann3_pathabundance.txt", sep = '\t', quote = "", 
                 header = TRUE, row.names = 1)

var<-read.table("all_samples_humann3_pathabundance_metadata.txt", header = TRUE, row.names = 1, sep ='\t')

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

bray<- vegdist(data_norm, method = 'bray')

ad<-adonis(data_norm ~ Plant * Month + Batch,  
           data = var, parallel = 4, 
           method = "bray")
ad

shannon<-diversity(data_norm, index = "shannon")

for_plotting<-merge(var,shannon, by = "row.names")

shapiro.test(shannon)
hist(shannon)
qqnorm(shannon)

shannon_krusk<-kruskal.test(shannon~Plant, data = for_plotting)
shannon_krusk
dunn.test(for_plotting$shannon, for_plotting$Plant, method = "bh")
shannon_lmer<-lmer(shannon ~ as.numeric(as.factor(Days_from_AUG_3_2020)) + 
                     (1|Plant) + (1|Batch), data = for_plotting)
summary(shannon_lmer)  

bray.df<-as.data.frame(as.matrix(bray))
bray_means<-as.data.frame(((rowMeans(bray.df))))
names(bray_means)[1]<-"Dissimilarity"
bray_means$Similarity<-(1-bray_means$Dissimilarity)

for_plotting<-cbind(var,bray_means, by = "row.names")


bray_lmer<-lmer(Similarity ~ as.numeric(as.factor(Days_from_AUG_3_2020)) + 
                  (1|Plant) + (1|Batch), data = for_plotting)
summary(bray_lmer)