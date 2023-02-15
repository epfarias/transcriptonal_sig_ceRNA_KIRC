#### Loading Packages ####
library("tidyverse")
library("OmicSelector")
if(!require("ggstatsplot")){remotes::install_github("IndrajeetPatil/ggstatsplot")}
if(!require("tidyverse")){install.packages("tidyverse")}
if(!require("ggfortify")){install.packages("ggfortify")}
if(!require("survival")){install.packages("survival")}
if(!require("cowplot")){install.packages("cowplot")}
if(!require("pals")){install.packages("pals")}
if(!require("colorspace")){install.packages("colorspace")}
if(!require("ggsci")){install.packages("ggsci")}
if(!require('pheatmap')) {install.packages('pheatmap')}

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

## Constructing the mirExp Data
mir_Exp <-  as.data.frame(mirExpr)
mir_Exp <- rownames_to_column(mir_Exp,"miRNAs_id")
mir_Exp <- mir_Exp[mir_Exp$miRNAs_id %in% nodes$symbol,] 
rownames(mir_Exp) <- mir_Exp$miRNAs_id
mir_Exp <- mir_Exp[,-1]
mir_Exp <- as.data.frame(t(mir_Exp))
mir_Exp <- rownames_to_column(mir_Exp,"sample")

#### Removing the extra Exps
rm(ceOutput, ceOutput2, deALL, deALL_MIR, DEGAll_DESeq2, DEGAll_DESeq2_MIR, deLNC,
   dePC, edges, rnaCounts, mirCounts,mirExpr,rnaExpr,nodes,metaMatrix.MIR,metaMatrix.RNA)

#### Constructing the Clinical Data ####
## Downloading the datasets from Xenabrowser
url <- "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.KIRC.sampleMap%2FKIRC_clinicalMatrix"
destfile <- "kirc_phenotype.tsv.gz"
download.file(url, destfile)

url <- "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2FKIRC_survival.txt"
destfile <- "kirc_survival.tsv.gz"
download.file(url, destfile)


## Extracting the data
kirc_phenotype <- read_tsv(gzfile("kirc_phenotype.tsv.gz"))
kirc_survival <- read_tsv(gzfile("kirc_survival.tsv.gz"))

## Merging the dataframes
kirc_clindata <- merge(x = kirc_phenotype, y = kirc_survival, by.x = "sampleID", by.y = "sample", all.y = TRUE)
kirc_clindata <- kirc_clindata[kirc_clindata$sample_type_id == "01",]
kirc_clindata <- kirc_clindata[, c(1,19,38,61:63,114:115)]

data_aarlens <- left_join(mir_Exp,rna_Exp,by = "sample")

data_aarlens <- merge(x = kirc_clindata, y = data_aarlens, by.x = "sampleID", by.y = "sample", all.y = TRUE)

rm(mir_Exp,rna_Exp, kirc_clindata, kirc_phenotype, kirc_survival, destfile, url)

##### Putting the gene signature symbols ####
sign <-  c("HMMR", "PTTG1", "HECW2","RFLNB", "INSR", "BTBD11", "AF117829.1", 
           "SNHG15", "hsa.miR.381.3p", "hsa.miR.130a.3p")

#### Constructing the dataset for Aarlen's Analysis ####
data_aarlens <- data_aarlens[!is.na(data_aarlens$sampleID), ]
data_aarlens <- data_aarlens[!is.na(data_aarlens$OS.time), ]

rownames(data_aarlens) <-  gsub("-", ".", data_aarlens$sampleID, fixed = T)
colnames(data_aarlens) <-  gsub("-", ".", colnames(data_aarlens), fixed = T)

table(data_aarlens$pathologic_T)

data_aarlens <- data_aarlens %>%
  mutate(pathologic_T = fct_collapse(pathologic_T, T1=c('T1', 'T1a', 'T1b'), T2=c('T2','T2a', 'T2b'), 
                                          T3=c('T3','T3a', 'T3b', 'T3c'), ))

data_aarlens <- data_aarlens[!is.na(data_aarlens$pathologic_T), ]

table(data_aarlens$gender)
data_aarlens$gender <-ifelse(data_aarlens$gender == "MALE", 1, 0)

sign <- intersect(sign, colnames(data_aarlens))

data_aarlens <- data_aarlens[,-1]

data_aarlens <- na.omit(data_aarlens)

#### Aalens additive regression model for censored data ####

aa_fit <-aareg(as.formula(paste0("Surv(OS.time, OS) ~ " , paste(sign, collapse= " + "))),
               data = data_aarlens)
aa_fit

summary(aa_fit)  # provides a more complete summary of results

labs=c('alphabet','alphabet2', 'glasbey','kelly','polychrome', 'stepped', 'stepped2', 'stepped3', 'tol', 'watlington')
op=par(mar=c(0,5,3,1))
pal.bands(alphabet(), alphabet2(), glasbey(), kelly(),
          polychrome(), stepped(), stepped2(), stepped3(), 
          tol(), watlington(), labels=labs, show.names=FALSE)

#vars <- c("Intercept", "age", "genderMale",  "metastasisM1", "metastasisMX",  "AL353637.1", "AL138830.2", "DPP6", "FOXJ1", "HHLA2", "MICG", "LIMCH1", "VSX1", "AR", "IL4", "CASR", "CRH", "GNB3", "SAA1")

vars <- c("HMMR", "PTTG1", "HECW2","RFLNB", "INSR", "BTBD11", "AF117829.1", 
          "SNHG15", "hsa.miR.381.3p", "hsa.miR.130a.3p")

variables <- rev(factor(vars, levels=vars))
#-- color decrease of lightness
cols1 <- as.vector(glasbey(15))
cols1 <- readhex(file = textConnection(paste(cols1, collapse = "\n")),
                 class = "RGB")
#transform to hue/lightness/saturation colorspace
cols1 <- as(cols1, "HLS")
#additive decrease of lightness
cols1@coords[, "L"] <- pmax(0, cols1@coords[, "L"] + 0.05)
cols1 <- as(cols1, "RGB")
cols1 <- hex(cols1)
p1 <- ggcoefstats(
  x = aa_fit,
  title = "Aalen's additive regression model",
  subtitle = "(for censored data)",
  only.significant = F,
  #point.args = list(color = "green", shape = 9),
  package = "pals",
  stats.label.args = list(
    max.time = 3,
    direction = "y",
    point.padding = 0.2,
    nudge_x = .01,
    nudge_y = .5,
    segment.curvature = -0.1,
    segment.angle = 10,
    segment.size  = 0.2,
    segment.linetype = 2
    # nudge_x = .15,
    # box.padding = 0.5,
    # nudge_y = 1,
    # segment.curvature = -0.1,
    # segment.ncp = 3,
    # segment.angle = 20
  ),
  palette = "glasbey",
  sort = "none", 
  k = 3
) + ggplot2::theme(text = element_text(size=18),
                   axis.text = element_text(size=16)
)
#+  ggplot2::scale_y_discrete(labels = vars) 

p1[["layers"]][[4]][["data"]][["expression"]][[1]] <- "list(~widehat(italic(beta))=='-5.37 \U00D7 10'^'-4', ~italic(z)=='-0.9230', ~italic(p)=='3.56 \U00D7 10'^'-1')"               
p1[["layers"]][[4]][["data"]][["expression"]][[2]] <- "list(~widehat(italic(beta))=='-1.71 \U00D7 10'^'-5', ~italic(z)=='-0.7650', ~italic(p)=='4.44 \U00D7 10'^'-1')" 
p1[["layers"]][[4]][["data"]][["expression"]][[3]] <- "list(~widehat(italic(beta))=='1.49 \U00D7 10'^'-4', ~italic(z)==' 2.4500', ~italic(p)=='1.43 \U00D7 10'^'-2')"  
p1[["layers"]][[4]][["data"]][["expression"]][[4]] <- "list(~widehat(italic(beta))=='-1.48 \U00D7 10'^'-4', ~italic(z)=='-2.2900', ~italic(p)=='2.18 \U00D7 10'^'-2')" 
p1[["layers"]][[4]][["data"]][["expression"]][[5]] <- "list(~widehat(italic(beta))=='1.04 \U00D7 10'^'-4', ~italic(z)=='1.4100', ~italic(p)=='1.59 \U00D7 10'^'-1')"  
p1[["layers"]][[4]][["data"]][["expression"]][[6]] <- "list(~widehat(italic(beta))=='-1.02 \U00D7 10'^'-4', ~italic(z)=='-0.6700', ~italic(p)=='5.03 \U00D7 10'^'-1')"  
p1[["layers"]][[4]][["data"]][["expression"]][[7]] <- "list(~widehat(italic(beta))=='3.64 \U00D7 10'^'-5', ~italic(z)=='1.4400', ~italic(p)==' 1.50 \U00D7 10'^'-1')"  
p1[["layers"]][[4]][["data"]][["expression"]][[8]] <- "list(~widehat(italic(beta))=='3.25 \U00D7 10'^'-4', ~italic(z)=='4.4800', ~italic(p)=='7.57 \U00D7 10'^'-6')"
p1[["layers"]][[4]][["data"]][["expression"]][[9]] <- "list(~widehat(italic(beta))=='-2.25 \U00D7 10'^'-5', ~italic(z)=='-0.0251', ~italic(p)=='9.80 \U00D7 10'^'-1')"  
p1[["layers"]][[4]][["data"]][["expression"]][[10]] <- "list(~widehat(italic(beta))=='4.23 \U00D7 10'^'-5', ~italic(z)=='1.8600', ~italic(p)=='6.27 \U00D7 10'^'-2')"
p1[["layers"]][[4]][["data"]][["expression"]][[11]] <- "list(~widehat(italic(beta))=='1.99 \U00D7 10'^'-4', ~italic(z)=='2.6100', ~italic(p)=='9.18 \U00D7 10'^'-3')"  



p2 <- ggplot2::autoplot(aa_fit) +
  theme(legend.position="none", 
        text = element_text(size=16),
        axis.text = element_text(size=10))
p2$layers[[1]]$data$variable <- factor(p2$layers[[1]]$data$variable,
                                       levels= variables)
p2$layers[[2]]$data$variable <- factor(p2$layers[[2]]$data$variable,
                                       levels= variables)
p2$data$variable <- factor(p2$data$variable,
                           levels= variables)
p2 <- p2 +
  scale_fill_manual(values=rev(cols1))

pdf("../figs/fig5_mRMR_aareg.pdf", width = 14, height = 8)
cowplot::plot_grid(p1, p2, labels = "auto", nrow=1, rel_widths = c(0.7,1))

plot(p1)

#### Odds Ratio ####
library("finalfit")

kirc_clindata <- kirc_clindata[kirc_clindata$pathologic_M != "MX", ]

data_aarlens <- data_aarlens[data_aarlens$pathologic_M != "MX", ]
data_aarlens$pathologic_M <- ifelse(data_aarlens$pathologic_M == "M0", 0, 1)

explanatory = sign

dependent = 'pathologic_M'

data_aarlens %>%
  or_plot(dependent, explanatory)

