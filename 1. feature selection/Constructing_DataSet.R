## Constructing the data set for OmicSelector

#### Loading Packages ####
library("tidyverse")
library("OmicSelector")

#### Loading the data ####
load("~/transcriptonal_sig_ceRNA_KIRC/0. data/Counts_Exp_KIRC.RData")
load("~/transcriptonal_sig_ceRNA_KIRC/0. data/GDCRNATools_KIRC.RData")

#### Correcting the anottation, selecting the Exps and transposing the matrix's ####

## Constructing the rnaExp Data
rna_Exp <- as.data.frame(rnaExpr)
rna_Exp <- rownames_to_column(rna_Exp, "ENSEMBLID")
rna_Exp <- merge(rna_Exp,nodes, by.y = "gene", by.x = "ENSEMBLID", all.x = TRUE)
rna_Exp <- rna_Exp[rna_Exp$ENSEMBLID %in% nodes$gene,]
rownames(rna_Exp) <- rna_Exp$symbol
rna_Exp <- rna_Exp[,-c(1,604:606)]
rna_Exp <- t(rna_Exp)
rna_Exp <- as.data.frame(rna_Exp)
rna_Exp <- rownames_to_column(rna_Exp,"sample")
#rna_Exp2 <- ## Retirar os sufixos -11

## Constructing the mirExp Data
mir_Exp <-  as.data.frame(mirExpr)
mir_Exp <- rownames_to_column(mir_Exp,"miRNAs_id")
mir_Exp <- mir_Exp[mir_Exp$miRNAs_id %in% nodes$symbol,] 
rownames(mir_Exp) <- mir_Exp$miRNAs_id
mir_Exp <- mir_Exp[,-1]
mir_Exp <- as.data.frame(t(mir_Exp))
mir_Exp <- rownames_to_column(mir_Exp,"sample")
#mir_Exp$bcr_patient_barcode <- substr(mir_Exp$bcr_patient_barcode,1,12)

#### Removing the extra Exps
rm(ceOutput, ceOutput2, deALL, deALL_MIR, DEGAll_DESeq2, DEGAll_DESeq2_MIR, deLNC,
   dePC, edges, rnaCounts, mirCounts,mirExpr,rnaExpr,nodes,metaMatrix.MIR,metaMatrix.RNA)

## Constructing the data frame
data("orginal_TCGA_data")
OmicSelector_table(table(orginal_TCGA_data$primary_site, orginal_TCGA_data$sample_type))

cancer_cases = filter(orginal_TCGA_data, primary_site == "Kidney" & sample_type == "PrimaryTumor")
cancer_cases = cancer_cases[cancer_cases$project_id == "TCGA-KIRC",]
control_cases = filter(orginal_TCGA_data, sample_type == "SolidTissueNormal")
control_cases = control_cases[control_cases$project_id == "TCGA-KIRC",]

cancer_cases$Class = "Case"
control_cases$Class = "Control"

dataset_Exp= rbind(cancer_cases, control_cases)
dataset_Exp = dataset_Exp[,c(1:55,2640)]

Exp <- left_join(mir_Exp,rna_Exp,by = "sample")

Exp <- left_join(dataset_Exp,Exp,by = "sample")

Exp <- Exp[,c(1:55,57:277,56)]

rm(cancer_cases,control_cases,mir_Exp, rna_Exp,orginal_TCGA_data,dataset_Exp)

#### Filtering only the samples with metastatic stage defined and lines NA's ####
Exp <- Exp[!(Exp$pathologic_M == "MX"),] # 30 pacientes sem informação
Exp$Class <- as.factor(Exp$pathologic_M)
lines_NA <- which(is.na(Exp))
sapply(Exp, function(x) sum(is.na(x)))
Exp <- Exp[-c(53,60,238,239,349,440),]

OmicSelector_table(table(Exp$Class), col.names = c("Class", "Number of cases"))

#### Balancear o data set ####
boxplot(Exp$age_at_diagnosis ~ Exp$Class)
barplot(prop.table(table(Exp$Class)),
        col = rainbow(2),
        main = "Metastasis Distribution")

t.test(Exp$age_at_diagnosis ~ Exp$Class)
OmicSelector_table(table(Exp$gender.x, Exp$Class))
chisq.test(Exp$gender.x, Exp$Class)

set.seed(101)
old_Exp = Exp  # backup
Exp = Exp[grepl("Adenocarcinomas", Exp$disease_type), ]
match_by = c("age_at_diagnosis", "gender.x")
tempdane = dplyr::select(Exp, all_of(match_by))
tempdane$Class = ifelse(Exp$Class == "M1", TRUE, FALSE)
suppressMessages(library(mice))
suppressMessages(library(MatchIt))
temp1 = mice(tempdane, m = 1)
temp2 = temp1$data
temp3 = mice::complete(temp1)
temp3 = temp3[complete.cases(temp3), ]
tempform = OmicSelector_create_formula(match_by)
temp3 <- temp3 %>% mutate(gender.x = factor(gender.x))
mod_match <- matchit(tempform, data = temp3)
newdata = match.data(mod_match)
Exp = Exp[as.numeric(rownames(newdata)), ]

barplot(prop.table(table(Exp$Class)),
        col = rainbow(2),
        main = "Metastasis Distribution")

boxplot(Exp$age_at_diagnosis ~ Exp$Class)
t.test(Exp$age_at_diagnosis ~ Exp$Class)
OmicSelector_table(table(Exp$gender.x, Exp$Class))
chisq.test(Exp$gender.x, Exp$Class)

#### Saving the Exp ####
save(Exp,file = "~/transcriptonal_sig_ceRNA_KIRC/0. data/Expression_OmicSelector.RData")
