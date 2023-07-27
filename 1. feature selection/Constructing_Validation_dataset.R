#### TESTE DO DADO DO ICGC RECA-EU ####
## instalar VROOM
#install.packages("vroom")
library("vroom")
library("tidyverse")
library("limma")
library("edgeR")

## In case you don't have the data from ICGC
url <- "https://dcc.icgc.org/api/v1/download?fn=/current/Projects/RECA-EU/exp_seq.RECA-EU.tsv.gz"
destfile <- "exp_seq.RECA-EU.tsv.gz"
download.file(url, destfile)

url <- "https://dcc.icgc.org/api/v1/download?fn=/current/Projects/RECA-EU/donor.RECA-EU.tsv.gz"
destfile <- "donor.RECA-EU.tsv.gz"
download.file(url, destfile)

donor_data <- vroom("donor.RECA-EU.tsv.gz")
count_data <- vroom("exp_seq.RECA-EU.tsv")

# hsa-miR-130a-3p = MIR130A (GENCARDS)
# hsa-miR-381-3p = MIR381 (GENCARDS)
transcriptional_signature<- c("RASD1", "BTBD11", "HMMR", "INSR", "HECW2", "RFLNB", "PTTG1",
                              "MIR130A", "MIR381", "SNHG15", "AF117829.1")

#### Donor Data ####
donor_data <- donor_data %>%
  dplyr::select(c("icgc_donor_id", "submitted_donor_id", "donor_survival_time", "donor_vital_status", "donor_age_at_diagnosis", "donor_sex",  "donor_tumour_stage_at_diagnosis" )) %>% 
  dplyr::rename(code = 'icgc_donor_id', 
                obs.time = 'donor_survival_time',
                status = 'donor_vital_status',
                age = 'donor_age_at_diagnosis', 
                gender = 'donor_sex')

donor_data$metastasis <- substr(donor_data$donor_tumour_stage_at_diagnosis, 5, 6)
donor_data$neoplasm <- substr(donor_data$donor_tumour_stage_at_diagnosis, 3, 4)
donor_data$ajcc.stage <- substr(donor_data$donor_tumour_stage_at_diagnosis, 1, 2)
donor_data$donor_tumour_stage_at_diagnosis <- NULL

#### Expression Data ####

# icgc.count$submitted_sample_id 
# RT Primary tumour - solid tissue
# RA Normal - tissue adjacent to primary
# N Normal - blood derived
# > sum(grepl("A$", unique(icgc.count$submitted_sample_id)))
# [1] 45
# > sum(!grepl("A$", unique(icgc.count$submitted_sample_id)))
# [1] 91
# Selecting only RT Primary tumour
# GRCh37 

count_data2 <- count_data %>% 
  dplyr::filter(!grepl("A$", submitted_sample_id)) %>%  
  dplyr::select(c("submitted_sample_id", "gene_id", "raw_read_count")) %>% 
  pivot_wider(names_from = "submitted_sample_id", values_from = "raw_read_count") %>% as.data.frame()

##Adjusting the RNA names

#loading annotables
#devtools::install_github("stephenturner/annotables")
library("annotables")

ensemblAnnot <- grch37 %>% dplyr::select(ensgene, symbol, description)

count_data3 <- merge(count_data2, ensemblAnnot, by.x = "gene_id", by.y = "ensgene", all.x = TRUE)

count_data3 <- count_data3[,c(1,93:94,2:92)]
count_data3 <- count_data3[,-3]

sum(duplicated(count_data3$gene_id))

count_data4 <- count_data3[!duplicated(count_data3$gene_id),]

duplicateds <- count_data4[duplicated(count_data4$symbol),]

for(i in 1:length(count_data4$symbol)){
  if(count_data4$gene_id[i] %in% duplicateds$gene_id || is.na(count_data4$symbol[i])){
    count_data4$symbol[i] <- count_data4$gene_id[i]
  }
}

sum(duplicated(count_data4$symbol))

# AF117829.1 = ENSG00000251136 (ENSEMBL)
# RFLNB = FAM101B = ENSG00000183688 (ENSEMBL)

count_data4[count_data4$gene_id == "ENSG00000251136", "symbol"] <- "AF117829.1"
count_data4[count_data4$gene_id == "ENSG00000183688", "symbol"] <- "RFLNB"

## Final Dataset
rownames(count_data4) <- count_data4$symbol

count_data4$symbol = NULL
count_data4$gene_id = NULL

count_final_ICGC <- count_data4

## Cleaning extra data
rm(count_data, count_data2,count_data3,count_data4,duplicateds,ensemblAnnot)

#### Normalizating the dataset with Voom ####
expr = edgeR::DGEList(counts = count_final_ICGC)
expr = edgeR::calcNormFactors(expr)

exp_norm_ICGC <- limma::voom(expr, design=NULL, plot = FALSE)$E

rm(expr, count_final_ICGC)

#### Constructing final dataset #### 
## Clinical data
donor_data$status <- ifelse(donor_data$status == "alive", 0, 1)
donor_data$gender <- ifelse(donor_data$gender == "female", 0, 1)

## Expression Data
exp_final_ICGC <- as.data.frame(t(exp_norm_ICGC))
exp_final_ICGC <- rownames_to_column(exp_final_ICGC, "patients_id")
exp_final_ICGC$patients_id <-stringr::str_remove(exp_final_ICGC$patients_id, "RT$")

## Merging datasets
exp_ICGC <- merge(exp_final_ICGC, donor_data, by.x = "patients_id", by.y = "submitted_donor_id", all.x = TRUE)
exp_ICGC <- exp_ICGC[,c(1,53600:53607,2:53599)]
exp_ICGC <- exp_ICGC[,-2]
exp_ICGC$dataset <- "ICGC"
exp_ICGC <- exp_ICGC[,c(1:8,53607,9:53606)]

## Final
rownames(exp_ICGC) <- exp_ICGC$patients_id
exp_ICGC$patients_id <- NULL

rm(donor_data, exp_final_ICGC, exp_norm_ICGC)

save(exp_ICGC, file = "~/Teste_ICGC/data/ICGC_Exp.RData")

#### SIGNATURE VALIDATION DATASET ####
#library(devtools)
#devtools::install_github("dalpozz/unbalanced")
suppressMessages(library(unbalanced))

library("tidyverse")

load("~/transcriptonal_sig_ceRNA_KIRC/0. data/Expression_OmicSelector.RData")

data_ICGC <- exp_ICGC
data_TCGA <- Exp

rm(Exp, exp_ICGC)

#### ICGC Dataset ####
data_ICGC <- data_ICGC %>%  
  dplyr::filter(metastasis %in% c("M0", "M1")) %>%
  droplevels(.)

data_ICGC$metastasis <- as.factor(data_ICGC$metastasis)

## Balancing ICGC dataset
predictor_variables <- data_ICGC[,-5] # Select everything except response
response_variable <- data_ICGC$metastasis

levels(response_variable) <- c('0', '1') 

undersampled_data <- ubBalance(predictor_variables, 
                               response_variable, 
                               type='ubUnder',         # Option for undersampling
                               verbose = TRUE)

data_ICGC_balanced <- cbind(undersampled_data$X,    # combine output
                            undersampled_data$Y)

names(data_ICGC_balanced)[names(data_ICGC_balanced) == "undersampled_data$Y"] <- "metastasis" # change name to class
levels(data_ICGC_balanced$metastasis) <- c('M0', 'M1')

rm(predictor_variables, response_variable, undersampled_data)

#### TCGA Dataset ####
data_TCGA <- data_TCGA[,c(5,8,13:16,38:40,56:277)]

data_TCGA <- data_TCGA %>% dplyr::rename(ajcc.stage = 'pathologic_T', 
                                         obs.time = 'days_to_last_follow_up',
                                         status = 'vital_status.x',
                                         age = 'age_at_initial_pathologic_diagnosis', 
                                         metastasis = 'pathologic_M',
                                         neoplasm = 'pathologic_N',
                                         gender = "gender.x",
                                         dataset = "project_id")

data_TCGA <- column_to_rownames(data_TCGA, 'submitter_id') 

data_TCGA <- data_TCGA %>%  
  dplyr::filter(metastasis %in% c("M0", "M1") ) %>%
  droplevels(.)

data_TCGA$status <- ifelse(data_TCGA$status == "Alive", 0, 1)
data_TCGA$gender <- ifelse(data_TCGA$gender == "female", 0, 1)

data_TCGA <- data_TCGA %>% rename("MIR130A" = "hsa-miR-130a-3p")
data_TCGA <- data_TCGA %>% rename("MIR381" = "hsa-miR-381-3p")

#### Final dataset ####
data_ICGC_balanced2 <- data_ICGC_balanced[,c("metastasis", "dataset",intersect(colnames(data_ICGC), unlist(transcriptional_signature$sign)))] 
data_TCGA2 <- data_TCGA[,c("metastasis", "dataset",intersect(colnames(data_ICGC), unlist(transcriptional_signature$sign)))]

data <- rbind(data_ICGC_balanced2,data_TCGA2)

rm(data_ICGC, data_ICGC_balanced, data_ICGC_balanced2, data_TCGA, data_TCGA2)

save(data, file = "Validation_dataset.RData")
