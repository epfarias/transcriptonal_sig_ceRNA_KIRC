############# KIRC_maftools ##################
# MAF files from Kidney Renal Clear Cell Carcinoma

# maftools: Summarize, Analyze and Visualize Mutation Anotated Files (MAF) Files
# URL: https://www.bioconductor.org/packages/release/bioc/html/maftools.html
# version 2.4.05

#### Installing packages ####
packages_bioconductor <- c("TCGAbiolinks","maftools","BSgenome.Hsapiens.UCSC.hg38","SummarizedExperiment")
packages_cran <- c("DT", "tidyverse", "stringr", "data.table", "pheatmap","NMF")

#use this function to check if each package is on the local machine
#if a package is installed, it will be loaded
#if any are not, the missing package(s) will be installed from Bioconductor and loaded
package.check <- lapply(packages_bioconductor, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
package.check <- lapply(packages_cran, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

rm(packages_cran, packages_bioconductor, package.check)

setwd("~/transcriptonal_sig_ceRNA_KIRC/3. genomic alterations")

## Reading KIRC Maf files ---------------------------

# download MAF aligned against hg38
# it saves data inside a GDCdata and project name directory 
query <- GDCquery(
  project = "TCGA-KIRC", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  legacy = FALSE, 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query)

#Remover duplicatas dentro do query 
query_results <- as.data.frame(query[[1]][[1]])
query_results <- query_results[!duplicated(query_results$cases),]
query[[1]][[1]] <- query_results

# GDC
maf <- GDCprepare(query)
sort(colnames(maf))

# MAF object contains main maf file, summarized data and any associated sample annotations
kirc.maf <- read.maf(maf = maf, useAll = T) 

# checking
getSampleSummary(kirc.maf) #  samples
getGeneSummary(kirc.maf) #  genes (hugo)
getClinicalData(kirc.maf) #  samples, no clinical data  
getFields(kirc.maf) #  variables 

# writes an output file
write.mafSummary(maf = kirc.maf, basename = 'kirc.maf')


## Reading clinical indexed data ------------------

# same as in data portal
clinical <- GDCquery_clinic(project = "TCGA-KIRC", type = "clinical", save.csv = FALSE)
sort(colnames(clinical))

colnames(clinical)[2] <- "Tumor_Sample_Barcode"
# clinical$Overall_Survival_Status <- 1
# clinical$Overall_Survival_Status[which(clinical$vital_status == "Dead")] <- 0
clinical$time <- clinical$days_to_death
clinical$time[is.na(clinical$days_to_death)] <- clinical$days_to_last_follow_up[is.na(clinical$days_to_death)]

# create object for survival analysis 
kirc.mafclin <- read.maf(maf = maf, clinicalData = clinical, isTCGA = T)


## Vizualizing -----------------------------------

# displays variants in each sample and variant types summarized by Variant_Classification
plotmafSummary(maf = kirc.maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = TRUE)

# oncoplot for top ten mutated genes (costumize oncoplots!)
gene.signature = c("HMMR", "RFLNB", "PTTG1", "BTBD11", "HECW2", "INSR") 
oncoplot(maf = kirc.maf, genes = gene.signature)

## Saving maf file -------------------------------
save(kirc.maf, kirc.mafclin, file = "~/transcriptonal_sig_ceRNA_KIRC/0. data/kirc_maf_mafclin.RData")

## Reading gistic or CNV -------------------------

all.lesions_95_01 <- read_file("~/transcriptonal_sig_ceRNA_KIRC/0. data/all_lesions.conf_95_01.txt")
amp.genes_95_01 <- read_file("~/transcriptonal_sig_ceRNA_KIRC/0. data/amp_genes.conf_95_01.txt")
del.genes_95_01 <- read_file("~/transcriptonal_sig_ceRNA_KIRC/0. data/del_genes.conf_95_01.txt")
scores.gis_95_01 <- read_file("~/transcriptonal_sig_ceRNA_KIRC/0. data/scores_95_01.gistic")

KIRC_95_01.gistic = readGistic(gisticAllLesionsFile = all.lesions_95_01, gisticAmpGenesFile = amp.genes_95_01,
                               gisticDelGenesFile = del.genes_95_01, gisticScoresFile = scores.gis_95_01, isTCGA = TRUE)

# GISTIC object
KIRC_95_01.gistic

# checking 
getSampleSummary(KIRC_95_01.gistic) # 534 samples (Tumor_Sample_Barcode)
getGeneSummary(KIRC_95_01.gistic) # 1062 genes (hugo)
getCytobandSummary(KIRC_95_01.gistic) # 365 cytobands
write.GisticSummary(gistic = KIRC_95_01.gistic, basename = 'KIRC_95_01.gistic')

# Save gistic data
save(KIRC_95_01.gistic,file = "~/transcriptonal_sig_ceRNA_KIRC/0. data/KIRC_gistic.RData")

## Vizualizing Gistic ---------------------------------

# genome plot
gisticChromPlot(gistic = KIRC_95_01.gistic, markBands = "all")

# bubble plot
gisticBubblePlot(gistic = KIRC_95_01.gistic)

# oncoplot 
dev.new()
gisticOncoPlot(gistic = KIRC_95_01.gistic, clinicalData = kirc.mafclin, 
               clinicalFeatures = 'ajcc_pathologic_stage', sortByAnnotation = TRUE, 
               top = 10, annotationFontSize = 1,legendFontSize = 1.2, fontSize = 0.5)  

gisticOncoPlot(gistic = KIRC_95_01.gistic, top = 10)
