#### Enrichment Analysis ####
library("clusterProfiler")
library("enrichplot")
library("tidyverse")
library("DOSE")
library("rtracklayer")


load("~/transcriptonal_sig_ceRNA_KIRC/0. data/GDCRNATools_KIRC.RData")


genes_sign <- c("INSR","HMMR","PTTG1","hsa.miR.381.3p","HECW2","AF117829.1",
                "RASD1","RFLNB","SNHG15","hsa.miR.130a.3p","BTBD11")

genes_sign <- genes_sign %>% as.data.frame()
colnames(genes_sign) <- "Symbol"

## Anottation -------------------------
gff <- import.gff("~/transcriptonal_sig_ceRNA_KIRC/0. data/gencode.v41.annotation.gff3")
gff<- as.data.frame(gff@elementMetadata) 
#Filtrar tres colunas do gff
gff <- gff[,c("gene_id" ,"gene_type", "gene_name")]
#Tirar linhas duplicadas
gff <- dplyr::distinct(gff)

gff$gene_id <- sub("\\..*", "", gff$gene_id)

genes_enrich <- merge(genes_sign , gff, by.x = "Symbol", by.y = "gene_name", all.x =  TRUE)

genes_enrich <- genes_enrich %>% na.omit()

genes_enrich <- genes_enrich[!(genes_enrich$gene_type == "lncRNA"),]

gene <- bitr(genes_enrich$Symbol, fromType = "SYMBOL", 
             toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

pc_genes <- merge(gene,dePC,by.x = "SYMBOL", by.y = "symbol",all.x = TRUE)

pc_geneList = pc_genes[,5]
names(pc_geneList) = as.character(pc_genes[,2])
pc_geneList = sort(pc_geneList, decreasing = TRUE)


## Lista para o KEGG
genes_kegg <- bitr_kegg(pc_genes$ENTREZID, fromType = "ncbi-geneid", 
                        toType = "kegg", organism = "hsa")

pc_genes_kegg <- merge(genes_kegg,pc_genes,by.x = "ncbi-geneid", by.y = "ENTREZID",all.x = TRUE)

pc_geneslist_kegg = pc_genes_kegg[,6]
names(pc_geneslist_kegg) = as.character(pc_genes_kegg[,2])
pc_geneslist_kegg = sort(pc_geneslist_kegg,decreasing = TRUE)

#### Removing the extra data ####
rm(ceOutput, ceOutput2, deALL, deALL_MIR, DEGAll_DESeq2, DEGAll_DESeq2_MIR, dePC,
   deLNC, edges, gff, nodes)

#### Enrichment Analysis ####

## KEGG 
gene_kegg <- names(pc_geneslist_kegg)

kk <- enrichKEGG(gene = gene_kegg,
                 organism = "hsa",
                 keyType = "ncbi-geneid",
                 pAdjustMethod = "BH")

p1 <- dotplot(kk)

p2 <- barplot(kk)

cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

## Gene Ontology (BP and MF)
go_BP <- enrichGO(gene = genes_enrich$gene_id,
                  OrgDb = "org.Hs.eg.db",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  keyType = 'ENSEMBL')

go_MF <- enrichGO(gene = genes_enrich$gene_id,
                  OrgDb = "org.Hs.eg.db",
                  ont = "MF",
                  pAdjustMethod = "BH",
                  keyType = 'ENSEMBL')

p3<- barplot(go_BP,
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways")

plot(p3)

p4 <- barplot(go_MF,
             drop = TRUE, 
             showCategory = 10, 
             title = "GO Molecular Functions")
plot(p4)

cowplot::plot_grid(p3, p4, ncol=1, labels=LETTERS[1:2])

## Disease Ontology
DO <- enrichDO(gene = names(pc_geneList), 
               ont = "DO",
               pAdjustMethod = "BH")
               
p5 <- barplot(DO, showCategory = 5)
plot(p5)

#### Using gprofiler 2 ####
library("gprofiler2")

gostres <- gost(query = genes_enrich$Symbol , organism = "hsapiens", significant =  FALSE)

gostplot(gostres, capped = TRUE, interactive = TRUE)


#### GDCRNATools ####
library("GDCRNATools")

source("~/transcriptonal_sig_ceRNA_KIRC/5. enrichment analysis/gdcEnrichPlot.R")

enrichment <- gdcEnrichAnalysis(gene = genes_enrich$gene_id, simplify = TRUE)

p7 <- gdcEnrichPlot(enrichment = enrichment, type = "bubble", category = "KEGG")
plot(p7)

p8 <- gdcEnrichPlot(enrichment = enrichment, type = "bubble", category = "GO_BP")
plot(p8)

p9 <- gdcEnrichPlot(enrichment = enrichment, type = "bubble", category = "GO_MF")
plot(p9)

p10 <- gdcEnrichPlot(enrichment = enrichment, type = "bubble", category = "DO")
plot(p10)


#### Enrichment Analisys by miRNA function Using DianaTools Information####

library("ggplot2")
library("tidyverse")
library("forcats")

KEGG_DT <- read_csv("~/transcriptonal_sig_ceRNA_KIRC/0. data/KEGG_DianaTools.csv")
KEGG_DT$p.adjust <- p.adjust(KEGG_DT$`p-value`, method = "BH", n = length(KEGG_DT$`p-value`))
KEGG_DT_mod <- KEGG_DT[-c(1,2,38,40,41,18,14,20,19,34,23,15,21,28,37,24,8,9,13,22,29,33,39,42),] ## removendo informações não importantes para o trabalho

GO_BP_DT <- read_csv("~/transcriptonal_sig_ceRNA_KIRC/0. data/GO_BP_DianaTools.csv")
GO_BP_DT$p.adjust <- p.adjust(GO_BP_DT$`p-value`, method = "BH", n = length(GO_BP_DT$`p-value`))
GO_BP_DT_mod <- GO_BP_DT[-c(1:6,9:16,18,20:24,28,33,37,38,41:46,48,55:57,59,60,62,64,65,68),] ## removendo informações não importantes para o trabalho

GO_MF_DT <- read_csv("~/transcriptonal_sig_ceRNA_KIRC/0. data/GO_MF_DianaTools.csv")
GO_MF_DT$p.adjust <- p.adjust(GO_MF_DT$`p-value`, method = "BH", n = length(GO_MF_DT$`p-value`))


KEGG_DT_plot <- ggplot(KEGG_DT) + geom_col(aes(x = fct_reorder(`KEGG pathway`,-log10(`p-value`)) , y = -log10(`p-value`), fill = -log10(`p-value`))) +
  coord_flip() + 
  labs(x = "KEGG Pathways", y = "-log10(P-value)", title = "KEGG Pathways altered by miRNAs in Gene Signature") + 
  theme(axis.text.y = element_text(size = 15))
plot(KEGG_DT_plot)

KEGG_DT_mod_plot <- ggplot(KEGG_DT_mod) + geom_dotplot(aes(x = fct_reorder(`KEGG pathway`,`#genes`) , y = `#genes`, fill = `p-value`), binaxis = "y", stackdir = "center", binwidth = 0.9) +
  coord_flip() + 
  labs(x = "KEGG Pathways", y = "Genes Regulated in Pathway")+ 
  theme(axis.text.y = element_text(size = 12)) +
  scale_fill_gradient(low = "red", high = "green") + 
  scale_size_area(max_size = 1) + 
  theme_light()
plot(KEGG_DT_mod_plot)

go_bp_plot <- ggplot(GO_BP_DT) + geom_col(aes(x = fct_reorder(`GO Category`,-log10(`p-value`)) , y = -log10(`p-value`), fill = -log10(`p-value`))) + 
  coord_flip() + 
  labs(x = "GO Biological Process", y = "P-value", title = "GO BP altered by miRNAs in Gene Signature") + 
  theme(axis.text.y = element_text(size = 15))
plot(go_bp_plot)

go_bp_mod_plot <- ggplot(GO_BP_DT_mod) + geom_dotplot(aes(x = fct_reorder(`GO Category`,`#genes`) , y = `#genes`, fill = `p-value`), binaxis = "y", stackdir = "center", binwidth = 10) +
  coord_flip() + 
  labs(x = "GO Biological Process", y = "Genes Regulated in Pathway", title = "GO BP altered by miRNAs in Gene Signature")+ 
  theme(axis.text.y = element_text(size = 12)) +
  scale_fill_gradient(low = "red",high = "blue") + 
  scale_size_area(max_size = 1) + 
  theme_light()
plot(go_bp_mod_plot)

go_mf_plot <- ggplot(GO_MF_DT) + geom_col(aes(x = fct_reorder(`GO Category`,-log10(`p-value`)) , y = -log10(`p-value`), fill = -log10(`p-value`))) + 
  coord_flip() + 
  labs(x = "GO Molecular Functions", y = "-log10(P-value)", title = "GO MF altered by miRNAs in Gene Signature")+ 
  theme(axis.title.y = element_text(size = 20))
plot(go_mf_plot)

go_mf_plot <- ggplot(GO_MF_DT) + geom_dotplot(aes(x = fct_reorder(`GO Category`,`#genes`) , y = `#genes`, fill = `p-value`), binaxis = "y", stackdir = "center", binwidth = 20) +
  coord_flip() + 
  labs(x = "GO Biological Process", y = "Genes Regulated in Pathway", title = "GO MF altered by miRNAs in Gene Signature")+ 
  theme(axis.text.y = element_text(size = 50)) +
  scale_fill_gradient(low = "red", high = "blue") + 
  scale_size_area(max_size = 1) + 
  theme_light()
plot(go_mf_plot)

## Figures to paper
# KEGG
cowplot::plot_grid(p7, KEGG_DT_mod_plot, ncol=1, labels=LETTERS[1:2])

#go_bp
cowplot::plot_grid(p3, go_bp_mod_plot, ncol=1, labels=LETTERS[1:2])

#go_mf
cowplot::plot_grid(p4, go_mf_plot, ncol=1, labels=LETTERS[1:2])
