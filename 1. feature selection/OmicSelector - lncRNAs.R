#### Feature Selection with OmicSelector ####
# Autor: Epitácio Farias 
# 19/12/2022

#### Defining workspace ####
setwd("~/transcriptonal_sig_ceRNA_KIRC/1. feature selection/")

#### Loading the packages ####
library("OmicSelector")

#### Loading the data ####
load("~/0. data/Expression_OmicSelector.RData")

#### Spliting the Dataset ####
library(caret)
# Step 1: Get row numbers for the training data
trainRowNumbers <- createDataPartition(Exp$pathologic_M, p=0.6, list=FALSE)

# Step 2: Create the training  dataset
trainData <- Exp[trainRowNumbers,]
trainData$mix = "train"
write.csv(trainData, "mixed_train.csv",row.names = F)

# Step 3: Get row numbers for the test data
temp <- Exp[-trainRowNumbers,]
testRowNumbers <- createDataPartition(temp$pathologic_M, p = 0.5, list = FALSE)

# Step 4: Create the test dataset
testData <- temp[testRowNumbers,]
testData$mix = "test"
write.csv(testData, "mixed_test.csv",row.names = F)

# Step 5: Create the Validation DataSet
valiData <- temp[-testRowNumbers,]
valiData$mix = "valid"
write.csv(valiData, "mixed_valid.csv",row.names = F)

# Step 6: Create the mixed dataset
mixed <- rbind(trainData,testData,valiData)
write.csv(mixed, "mixed.csv",row.names = F)

# See a slipt summary
OmicSelector_table(table(mixed$Class, mixed$mix))

# Cleaning the workspace
rm(temp, testRowNumbers, trainRowNumbers, mixed, testData, trainData, valiData)

#### Basic Analisys Exploratory ####
# First load the function with the adjusts to use all the dataset not only the miRNAs information
source("OmicSelector_Load_DataMix.R") 
dane = OmicSelector_load_datamix() # load mixed_*.csv files
train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]] # get the objects from list to make the code more readable.

#### Feature Selection ####
# m = c(33,45,51,53,55,57,58,59)
# Run the code "FeaturesSelection_Techniques.R" to execute the m techniques.

selected_sets_of_miRNAs = OmicSelector_merge_formulas(max_miRNAs = 100)

all_sets = readRDS("featureselection_formulas_all.RDS")
length(all_sets) # How many feature selection methods completed in time?
final_sets = readRDS("featureselection_formulas_final.RDS")
length(final_sets) # How many feature selection methods completed in time and fulfilled max_miRNA criteria? (remember about fcsig and cfs_sig)
featureselection_formulas_final = fread("featureselection_formulas_final.csv")
OmicSelector_table(featureselection_formulas_final) # show information about selected formulas

hist(featureselection_formulas_final$ile_miRNA, 
     breaks = ncol(train),
     main = "Number of selected microRNAs distribution",
     xlab = "Number of selected microRNAs"
) # Histogram showing how many miRNAs were selected in final set.
psych::describe(featureselection_formulas_final$ile_miRNA) # Descriptive statistics of how many features where selected in the final set.

#### Benchmarking ####
algorithms = c("rf", "xgbTree", "svmRadial")
source("OmicSelector_benchmark.R")
OmicSelector_tutorial_balanced_benchmark = OmicSelector_benchmark(input_formulas = readRDS("featureselection_formulas_final.RDS"),
                                                                  search_iters = 10, # 5 random hyperparameter sets will be checked; 5 is set here for speed purposes.. for real projects use more, like 5000...
                                                                  algorithms = algorithms, # just add ctree, note that logistic regression (glm) is always included
                                                                  output_file = paste0("benchmark.csv")) # the main output

metody = OmicSelector_get_benchmark_methods("benchmark.csv") # gets the methods used in benchmark
par(mfrow = c(2,2))
for(i in 1:length(metody)){
  temp = OmicSelector_get_benchmark("benchmark.csv") # loads benchmark
  temp2 = dplyr::select(temp, starts_with(paste0(metody[i],"_")))
  boxplot(temp[,paste0(metody[i],"_train_Accuracy")], temp[,paste0(metody[i],"_test_Accuracy")], temp[,paste0(metody[i],"_valid_Accuracy")],
          main = paste0("Method: ", metody[i]), names = c("Training","Testing","Validation"), ylab = "Accuracy", ylim = c(0,1))
  tempids = c(match(paste0(metody[i],"_train_Accuracy"), colnames(temp)), match(paste0(metody[i],"_test_Accuracy"), colnames(temp)), match(paste0(metody[i],"_valid_Accuracy"), colnames(temp)))
}
par(mfrow = c(2,2))

par(mfrow = c(2,2))
for(i in 1:length(metody)){
  temp3 = OmicSelector_get_benchmark("benchmark.csv") # loads benchmark
  temp4 = dplyr::select(temp, starts_with(paste0(metody[i],"_")))
  boxplot(temp[,paste0(metody[i],"_train_Sensitivity")], temp[,paste0(metody[i],"_test_Sensitivity")], temp[,paste0(metody[i],"_valid_Sensitivity")],
          main = paste0("Method: ", metody[i]), names = c("Training","Testing","Validation"), ylab = "Sensitivity", ylim = c(0,1))
  tempids = c(match(paste0(metody[i],"_train_Sensitivity"), colnames(temp)), match(paste0(metody[i],"_test_Sensitivity"), colnames(temp)), match(paste0(metody[i],"_valid_Sensitivity"), colnames(temp)))
}
par(mfrow = c(2,2))

par(mfrow = c(2,2))
for(i in 1:length(metody)){
  temp5 = OmicSelector_get_benchmark("benchmark.csv") # loads benchmark
  temp6 = dplyr::select(temp, starts_with(paste0(metody[i],"_")))
  boxplot(temp[,paste0(metody[i],"_train_Specificity")], temp[,paste0(metody[i],"_test_Specificity")], temp[,paste0(metody[i],"_valid_Specificity")],
          main = paste0("Method: ", metody[i]), names = c("Training","Testing","Validation"), ylab = "Specificity", ylim = c(0,1))
  tempids = c(match(paste0(metody[i],"_train_Specificity"), colnames(temp)), match(paste0(metody[i],"_test_Specificity"), colnames(temp)), match(paste0(metody[i],"_valid_Specificity"), colnames(temp)))
}
par(mfrow = c(2,2))


## Assinatura que teve a melhor acurácia no treino, teste e validação:
acc1 = OmicSelector_best_signature_proposals(benchmark_csv = "benchmark.csv", without_train = F)  # generates the benchmark sorted by metaindex
best_signatures = acc1[1:3, ]  # get top 3 methods
OmicSelector_table(best_signatures[, c("metaindex", "method", "miRy")])

## Assinatura com melhor acurácia no teste e validação:
acc1 = OmicSelector_best_signature_proposals(benchmark_csv = "benchmark.csv", without_train = T)  # generates the benchmark sorted by metaindex
best_signatures = acc1[1:3, ]  # get top 3 methods
OmicSelector_table(best_signatures[, c("metaindex", "method", "miRy")])

## Assinatura com melhor sensitividade e especificidade na validação:
acc = OmicSelector_best_signature_proposals_meta11(benchmark_csv = "benchmark.csv")  # generates the benchmark sorted by metaindex
best_signatures = acc[1:4, ]  # get top 3 methods
OmicSelector_table(best_signatures[, c("metaindex", "method", "miRy")])

## Visualizar o over/underfitting baseado no score de acurácia:
for (i in 1:length(metody)) {
  suppressMessages(library(PairedData))
  suppressMessages(library(profileR))
  pd = paired(as.numeric(acc[1:4, paste0(metody[i], "_train_Accuracy")]), as.numeric(acc[1:4, paste0(metody[i], "_test_Accuracy")]))
  pd = cbind(pd,as.numeric(acc[1:4, paste0(metody[i], "_valid_Accuracy")]))
  colnames(pd) = c("Train Accuracy", "Test Accuracy","Valid Accuracy")
  plot1 = OmicSelector_profileplot(pd, Method.id = acc$method[1:4], standardize = F)

  pd2 = paired(as.numeric(acc[1:4, paste0(metody[i], "_train_Sensitivity")]), as.numeric(acc[1:4, paste0(metody[i], "_test_Sensitivity")]))
  pd2 = cbind(pd2,as.numeric(acc[1:4, paste0(metody[i], "_valid_Sensitivity")]))
  colnames(pd2) = c("Train Sensitivity", "Test Sensitivity","Valid Sensitivity")
  plot2 = OmicSelector_profileplot(pd2, Method.id = acc$method[1:4], standardize = F)
  
  pd3 = paired(as.numeric(acc[1:4, paste0(metody[i], "_train_Specificity")]), as.numeric(acc[1:4, paste0(metody[i], "_test_Specificity")]))
  pd3 = cbind(pd3,as.numeric(acc[1:4, paste0(metody[i], "_valid_Specificity")]))
  colnames(pd3) = c("Train Specificity", "Test Specificity","Valid Specificity")
  plot3 = OmicSelector_profileplot(pd3, Method.id = acc$method[1:4], standardize = F)
  
  require(gridExtra)
  grid.arrange(arrangeGrob(plot1,plot2,plot3, ncol = 2, nrow = 2, top = metody[i]))
}

## Visualizar relação entre acurácia de teste e validação:
acc2 = acc[1:4, ]  # get top 3 methods
accmelt = melt(acc2, id.vars = "method") %>%
  filter(variable != "metaindex") %>%
  filter(variable != "miRy")
accmelt = cbind(accmelt, strsplit2(accmelt$variable, "_"))
acctest = accmelt$value[accmelt$`2` == "test"]
accvalid = accmelt$value[accmelt$`2` == "valid"]
accmeth = accmelt$method[accmelt$`2` == "test"]
unique(accmeth)
plot5 = ggplot(,aes(x = as.numeric(acctest), y = as.numeric(accvalid), shape = accmeth)) + 
  geom_point() + scale_x_continuous(name = "Accuracy on test set",limits = c(0.5, 1)) + scale_y_continuous(name = "Accuracy on validation set", limits = c(0.5, 1)) + theme_bw()
grid.arrange(arrangeGrob(plot5, ncol = 1, nrow = 1))

#### Best Signature ####
for (i in 1:4) {
  cat(paste0("\n\n## ", acc$method[i], "\n\n"))
  par(mfrow = c(1, 2))
  acc = OmicSelector_best_signature_proposals_meta11("benchmark.csv")
  metody = OmicSelector_get_benchmark_methods("benchmark.csv")
  ktory_set = match(acc$method[i], OmicSelector_get_benchmark("benchmark.csv")$method)
  # do_ktorej_kolumny = which(colnames(acc) == 'metaindex') barplot(as.numeric(acc[i,1:do_ktorej_kolumny]))
  for (ii in 1:length(metody)) {
    
    temp = OmicSelector_get_benchmark("benchmark.csv") %>%
      dplyr::select(starts_with(paste0(metody[ii], "_t")), starts_with(paste0(metody[ii], "_v")))
    
    ROCtext = paste0("Training AUC ROC: ", round(temp[ktory_set, 1], 2), " (95%CI: ", round(temp[ktory_set, 2], 2), "-", round(temp[ktory_set, 3],
                                                                                                                               2), ")")
    
    temp = temp[, -c(1:4)]
    temp2 = as.numeric(temp[ktory_set, ])
    temp3 = matrix(temp2, nrow = 3, byrow = T)
    colnames(temp3) = c("Accuracy", "Sensitivity", "Specificity")
    rownames(temp3) = c("Training", "Testing", "Validation")
    temp3 = t(temp3)
    
    plot1 = barplot(temp3, beside = T, ylim = c(0, 1), xlab = paste0(ROCtext, "\nBlack - accuracy, blue - sensitivity, green - specificity"), width = 0.85,
                    col = c("black", "blue", "green"), legend = F, args.legend = list(x = "topright", bty = "n", inset = c(0, -0.25)), cex.lab = 0.7, main = paste0(acc$method[i],
                                                                                                                                                                    " - ", metody[ii]), font.lab = 2)
    ## Add text at top of bars
    text(x = plot1, y = as.numeric(temp3), label = paste0(round(as.numeric(temp[ktory_set, ]) * 100, 1), "%"), pos = 3, cex = 0.6, col = "red")
  }
  par(mfrow = c(1, 1))
  
}

## Assessing the features common in the top3 feature selection models
overlap = OmicSelector_signature_overlap(acc$method[1:3], "benchmark.csv")

attr(overlap, "intersections")

#### Constructin an Heatmap with the Gene Signatures and Models ####

method = OmicSelector_get_benchmark_methods("benchmark.csv")

benchmark = read_csv("benchmark.csv")

roc_auc_bench = data.frame()
roc_auc_bench <- as.data.frame(benchmark$method)

for(i in 1:length(method)){
  roc_auc_bench[,i+1] = benchmark[,paste0(method[i],"_train_ROCAUC")]  
}

row.names(roc_auc_bench) <- roc_auc_bench$`benchmark$method`
roc_auc_bench <- roc_auc_bench[,-1]

roc_auc_bench <- data.matrix(roc_auc_bench)

roc_auc_bench <- melt(roc_auc_bench)

p1 <- ggplot(roc_auc_bench, aes(X2, X1)) +
  geom_tile(aes(fill = value)) + 
  geom_text(aes(label = round(value, 3)), colour = "white") +
  scale_fill_gradient(low = "yellow", high = "darkred") +
  ggtitle("Mean AUC on TCGA-KIRC dataset with 10 CV")+
  ylab("Feature Selection Methods") +
  xlab("Models") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_x_discrete(expand=c(0,0))
p1

#### Looking for the unique genes in the signatures ####
signatures <- read_csv("featureselection_formulas_final.csv")
genes <- strsplit(signatures$formula, split =" +")
genes <- lapply(1:length(genes), function(n) setdiff(genes[[n]], c('+','~','Class')))
names(genes) <- signatures$name

df <- data.frame(fromList(genes)) 
rownames(df) <- unique(unlist(genes))
df$gene <- unique(unlist(genes))

df.tidy <- df %>% tidyr::gather(signature, belong, 1:9)
df.tidy$Belongs <- factor(df.tidy$belong, levels = c(0, 1))
df.tidy$Gene <- factor(df.tidy$gene, levels = unique(unlist(genes)) ) 
df.tidy$Signature<- factor(df.tidy$signature, levels = signatures$name)

gg <- ggplot(df.tidy, aes(x=Gene, y=Signature, fill=Belongs)) + 
  coord_equal() +
  geom_tile(color="white", linewidth = 0.5) + 
  #scale_fill_brewer(palette = "Set1", direction = -1) +
  scale_fill_manual(values = c("gray", "red")) +
  theme(axis.ticks=element_blank()) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

plot(gg)

pdf("fig_heatplot.pdf", width = 12, height = 5)
gg
dev.off()


png("fig_heatplot.png", width = 3000, height = 1200, units = "px", res = 300)
gg
dev.off()

