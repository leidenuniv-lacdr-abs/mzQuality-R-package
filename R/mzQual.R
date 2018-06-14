rm(list = ls())
setwd("~/Documents/Work/abs/mzQuality-R-package/R/")
library(reshape2)
library(tcltk2)
library(gsubfn)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggfortify)

source("pcaplot.R")


com = read.delim("../data/combined.tsv")



com$compound = gsubfn(".", list(" " = "_", "&" = "_", "*" = "_star", "," = "_", "%" = "_perc", "-" = "_"), as.character(com$compound))
com$ratio = com$area / com$area_is
com = com[order(com$batch, com$order),]
com$position = 1:nrow(com)

com <- com %>% 
  group_by(compound) %>% 
  mutate(area = ifelse(area < 0, min(area[which(area >0)])*0.1,area))

if("injection" %in% colnames(com)){
  com$injection2 = paste(com$sample, com$batch, com$injection)}else{
    com$injection2 = paste(com$sample, com$batch)
  }

com$rep = 0
if("replicate" %in% colnames(df) == T){
  com$repl_sample = ifelse(com$replicate %in% c('', '-', '_', 'a'), 0, 1)
  for (un in unique(com$injection2[which(com$type == "sample")])) {
    ab = as.character(unique(com$replicate[which(com$injection2 == un)]))
    if(length(ab) == 2){
      com$rep[which(com$injection2 == un)] = 1
    }
  }
}else{print("no replicates")}


com <- com %>% 
  group_by(compound) %>% 
  mutate(compound_qc_ratio_median = median(ratio[which(com$type == "qc" & ratio > 0)], na.rm = T))## getting median over all QCs in all batches per compound

### calculating blank effects per median of samples ####
com <- com %>%
  group_by(compound, batch) %>%  ### calculating per batch
  mutate(be = mean(area[which(type == "blank")], na.rm = T) / median(area[which(type == "sample")], na.rm= T), ## calculating blank effect
         be_perc = 100*mean(area[which(type == "blank")], na.rm = T) / median(area[which(type == "sample")], na.rm= T), ## calculating blank effect in percent
         qc_correct_factor =  compound_qc_ratio_median / median(area[which(type == "qc" & area > 0)]) ## calculating ratio of median over all QCs in all batches over median over all QCs in batches
  )

blank_effect <- unique(com[,c("compound", "batch", "be", "be_perc")])


### qc calculations performed ####

com <- com %>% 
  group_by(batch, compound) %>% 
  mutate(med_ratio = median(ratio[which(com$type == "qc")], na.rm = T), 
         rt_mean = mean(rt, na.rm = T), 
         rt_stdev = sd(rt, na.rm = T))

com$qc_correct_factor = com$compound_qc_ratio_median / com$med_ratio
com$inter_median_qc_corrected = com$ratio * com$qc_correct_factor

#### rsdrep calculations ####
if(length(which(com$rep == 1)) > 1){
sa <- com[which(com$rep == 1),] %>%
  group_by(compound, injection2) %>%  ### calculating per compound, batch, sample and injection
  mutate(rsdrep_nc = 100 * (sd(area, na.rm = T) / mean(area, na.rm = T)),
         rsdrep_is_corrected = 100 * (sd(ratio, na.rm = T) / mean(ratio, na.rm = T)),
         rsdrep_inter_median_qc_corrected = 100 * (sd(inter_median_qc_corrected, na.rm =T) / mean(inter_median_qc_corrected, na.rm =T))
  )

rsdrep_df <- sa %>%
  group_by(compound, batch) %>%  ### calculating per batch, sample and compound
  summarise_at(vars(c(rsdrep_nc, rsdrep_is_corrected, rsdrep_inter_median_qc_corrected)), funs(mean(., na.rm=TRUE)))
rsdrep_df_all <-  rsdrep_df %>% 
  group_by(compound) %>% 
  summarise_at(vars(c(rsdrep_nc, rsdrep_is_corrected, rsdrep_inter_median_qc_corrected)), funs(mean(., na.rm = T)))
}

#### rsdqc calculations ####
qc <- com[which(com$type == "qc"),] %>% 
  group_by(compound, batch) %>% 
  mutate(
    rsdqc_nc = 100 * (sd(area, na.rm = T)/mean(area, na.rm = T)),
    rsdqc_is_corrected = 100 * (sd(ratio, na.rm = T)/mean(ratio, na.rm = T)),
    rsdqc_inter_median_qc_corrected = 100 * (sd(inter_median_qc_corrected, na.rm = T)/mean(inter_median_qc_corrected, na.rm = T))
  )

rsdqc_df <- qc %>%
  group_by(compound, batch) %>%  ### calculating per batch, sample and compound
  summarise_at(vars(c(rsdqc_nc, rsdqc_is_corrected, rsdqc_inter_median_qc_corrected)), funs(mean(., na.rm=TRUE)))
rsdqc_df_all <-  rsdqc_df %>% 
  group_by(compound) %>% 
  summarise_at(vars(c(rsdqc_nc, rsdqc_is_corrected, rsdqc_inter_median_qc_corrected)), funs(mean(., na.rm = T)))

#### calculating rt shifts ####
com$rt_shift = com$rt - com$rt_mean

rt_shifts = unique(com[,c("batch", "compound", "rt_mean", "rt_shift","rt_stdev",  "sample")])


# write.table(rt_shifts, file = "rt_shifts.txt", quote = F, row.names = F, sep = "\t")
# write.table(rsdrep_df, file = "batch_reprsd.txt", quote = F, row.names = F, sep = "\t")
# write.table(rsdqc_df, file = "batch_repqc.txt", quote = F, row.names = F, sep = "\t")
# write.table(blank_effect, file = "blank_effect.txt", quote = F, row.names = F, sep = "\t")
# write.table(rsdrep_df_all, file = "reprsd.txt", quote = F, row.names = F, sep = "\t")
# write.table(rsdqc_df_all, file = "repqc.txt", quote = F, row.names = F, sep = "\t")


######### return sample by measurements dataframe #######


com2 <- ddply(com[which(com$type == "sample"),], .(compound), function(x){ x$id2 = 1:nrow(x); x})
com2 <- com2[order(com2$batch),]
com2$sample2 = make.unique(as.character(com2$sample))

com2$batch = factor(com2$batch, levels = sort(unique(com2$batch), decreasing = F))
com2$sample2 = factor(com2$sample2, levels = com2$sample2[order(com2$batch)])


ggplot(com2, aes(x = sample2, y = ratio, fill = batch)) + geom_boxplot() + 
  theme_bw() + ggtitle("sample raw ratio") 
ggplot(com2, aes(x = sample2, y = log2(ratio), fill = batch)) + geom_boxplot() + 
  theme_bw() + ggtitle("sample log2 raw ratio")
ggplot(com2, aes(x = sample2, y = inter_median_qc_corrected, fill = batch)) + geom_boxplot() + 
  theme_bw() + ggtitle("sample inter_median_qc_corrected")
ggplot(com2, aes(x = sample2, y = log2(inter_median_qc_corrected), fill = batch)) + geom_boxplot() + 
  theme_bw() + ggtitle("sample log2 inter_median_qc_corrected")


# ggarrange(p1,p2,p3,p4, nrow = 2, ncol = 2, align = "hv")


ggplot(com2, aes(x = sample, y = ratio, fill = batch)) + geom_boxplot() +
  theme_bw() + ggtitle("sample raw ratio")
ggplot(com2, aes(x = sample, y = log2(ratio), fill = batch)) + geom_boxplot() +
  theme_bw() + ggtitle("sample log2 raw ratio")
ggplot(com2, aes(x = sample, y = inter_median_qc_corrected, fill = batch)) + geom_boxplot() +
  theme_bw() + ggtitle("sample inter_median_qc_corrected")
ggplot(com2, aes(x = sample, y = log2(inter_median_qc_corrected), fill = batch)) + geom_boxplot() +
  theme_bw() + ggtitle("sample log2 inter_median_qc_corrected")

##### diagnostic plots ####
com <- ddply(com, .(compound), function(x){ x$id2 = 1:nrow(x); x})

df1 = as.data.frame(acast(com, id2+sample ~compound, value.var = "area"))

df1$sample = sapply(rownames(df1), function(x){a = strsplit(x, split = "_");
b = paste(a[[1]][2:length(a[[1]])], collapse = "_");
b})
df1$type = "sample"
df1$type[grep("BLANK",df1$sample)] = "BLANK"
df1$type[grep("QC",df1$sample)] = "QC"
df1$type[grep("CAL",df1$sample)] = "CAL"

df1$batch = 0
for (i in 1:nrow(df1)){
  df1$batch[i] = unique(com$batch[which(com$sample == df1$sample[i])])
}

df2 =as.data.frame(apply(df1[,-which(colnames(df1) %in% c("sample", "type", "batch"))], 2, function(x) {log2(as.numeric(x))}))
batch = df1$batch
type =as.factor(df1$type)

pcaPlot(df2, grps = as.factor(type), title = "All types - Area, type")
pcaPlot(df2, grps = as.factor(batch), title = "All types - Area, batch")


df2 =scale(as.data.frame(apply(df1[which(df1$type == "sample"),-which(colnames(df1) %in% c("sample", "type", "batch"))], 2, function(x) {log2(as.numeric(x))})))
batch = df1$batch[which(df1$type == "sample")]
type =as.factor(df1$type[which(df1$type == "sample")])

pcaPlot(df2, grps = as.factor(batch), title = "Only samples - Area, batch")


##### Ratio uncorrected 
df1 = as.data.frame(acast(com, id2+sample ~compound, value.var = "ratio"))

df1$sample = sapply(rownames(df1), function(x){a = strsplit(x, split = "_");
b = paste(a[[1]][2:length(a[[1]])], collapse = "_");
b})
df1$type = "sample"
df1$type[grep("BLANK",df1$sample)] = "BLANK"
df1$type[grep("QC",df1$sample)] = "QC"
df1$type[grep("CAL",df1$sample)] = "CAL"

df1$batch = 0
for (i in 1:nrow(df1)){
  df1$batch[i] = unique(com$batch[which(com$sample == df1$sample[i])])
}

df2 =as.data.frame(apply(df1[,-which(colnames(df1) %in% c("sample", "type", "batch"))], 2, function(x) {log2(as.numeric(x))}))
batch = df1$batch
type =as.factor(df1$type)

pcaPlot(df2, grps = as.factor(type), title = "All types - Ratio, type")
pcaPlot(df2, grps = as.factor(batch), title = "All types - Ratio, batch")


df2 =scale(as.data.frame(apply(df1[which(df1$type == "sample"),-which(colnames(df1) %in% c("sample", "type", "batch"))], 2, function(x) {log2(as.numeric(x))})))
batch = df1$batch[which(df1$type == "sample")]
type =as.factor(df1$type[which(df1$type == "sample")])

pcaPlot(df2, grps = as.factor(batch), title = "Only samples - Ratio, batch")




#### QC corrected metabolites ####
df1 = as.data.frame(acast(com, id2+sample ~compound, value.var = "inter_median_qc_corrected"))

df1$sample = sapply(rownames(df1), function(x){a = strsplit(x, split = "_");
                                  b = paste(a[[1]][2:length(a[[1]])], collapse = "_");
                                  b})
df1$type = "sample"
df1$type[grep("BLANK",df1$sample)] = "BLANK"
df1$type[grep("QC",df1$sample)] = "QC"
df1$type[grep("CAL",df1$sample)] = "CAL"

df1$batch = 0
for (i in 1:nrow(df1)){
  df1$batch[i] = unique(com$batch[which(com$sample == df1$sample[i])])
}

df2 =scale(as.data.frame(apply(df1[,-which(colnames(df1) %in% c("sample", "type", "batch"))], 2, function(x) {log2(as.numeric(x))})))
batch = df1$batch
type =as.factor(df1$type)

pcaPlot(df2, grps = as.factor(type), title = "All types - QC corrected, type")
pcaPlot(df2, grps = as.factor(batch), title = "All types - QC corrected, batch")


df2 =as.data.frame(apply(df1[which(df1$type == "sample"),-which(colnames(df1) %in% c("sample", "type", "batch"))], 2, function(x) {log2(as.numeric(x))}))
batch = df1$batch[which(df1$type == "sample")]
type =as.factor(df1$type[which(df1$type == "sample")])

pcaPlot(df2, grps = as.factor(batch), title = "Only samples - QC corrected, batch")


