rm(list = ls())
setwd("~/Documents/Work/abs/mzQuality/R")
com = read.delim("../data/combined.tsv")
com$ratio = com$area / com$area_is

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

com$injection2 = paste(com$sample, com$batch, com$injection)

com$rep = 0
for (un in unique(com$injection2[which(com$type == "sample")])) {
  ab = as.character(unique(com$replicate[which(com$injection2 == un)]))
  if(length(ab) == 2){
    com$rep[which(com$injection2 == un)] = 1
  }
}

### qc calculations performed ####

com <- com %>% 
  group_by(batch, compound) %>% 
  mutate(med_ratio = median(ratio[which(com$type == "qc")], na.rm = T), 
         rt_mean = mean(rt, na.rm = T), 
         rt_stdev = sd(rt, na.rm = T))

com$qc_correct_factor = com$compound_qc_ratio_median / com$med_ratio
com$inter_median_qc_corrected = com$ratio * com$qc_correct_factor

#### rsdrep calculations ####
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


write.table(rt_shifts, file = "rt_shifts.txt", quote = F, row.names = F, sep = "\t")
write.table(rsdrep_df, file = "batch_reprsd.txt", quote = F, row.names = F, sep = "\t")
write.table(rsdqc_df, file = "batch_repqc.txt", quote = F, row.names = F, sep = "\t")
write.table(blank_effect, file = "blank_effect.txt", quote = F, row.names = F, sep = "\t")
write.table(rsdrep_df_all, file = "reprsd.txt", quote = F, row.names = F, sep = "\t")
write.table(rsdqc_df_all, file = "repqc.txt", quote = F, row.names = F, sep = "\t")

