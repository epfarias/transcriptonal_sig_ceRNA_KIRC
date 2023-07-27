#### Validation Train with TCGA Test ICGC ####

#install.packages("mlr3proba", repos = "https://mlr-org.r-universe.dev")
packages_cran = c("mlr3", "mlr3verse", "mlr3learners", "mlr3extralearners", "tidyverse", "kknn", "apcluster", "xgboost", "e1071", "caret")

package.check <- lapply(packages_cran, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

rm(packages_cran, package.check)


#### Signature information and validation dataset ####

# hsa-miR-130a-3p = MIR130A (GENCARDS)
# hsa-miR-381-3p = MIR381 (GENCARDS)
transcriptional_signature<- data.frame(matrix(0, nrow=1, ncol=2))

colnames(transcriptional_signature) <- c("sign", "ID")
transcriptional_signature$sign <- c("RASD1,BTBD11,HMMR,INSR,HECW2,RFLNB,PTTG1,MIR130A,MIR381,SNHG15,AF117829.1")
transcriptional_signature$sign <- stringr::str_split(transcriptional_signature$sign,",")

transcriptional_signature$ID <- "sign.1"

load("~/Teste_ICGC/data/Validation_dataset.RData")
data <- data %>% na.omit()

## Creating Task ##
id <- transcriptional_signature$ID
cols <- colnames(data) %in% unlist(c("metastasis", intersect(colnames(data), unlist(transcriptional_signature$sign))))
assign(id, TaskClassif$new(data[,cols],id = id, target = "metastasis", positive = "M1"))

list.clas.task <- list(sign.1)

rm(list=ls(pattern = "(^filt.)|(^fs.)|(^sign.)"), cols, id)

# Execute benchmark ----

bmr_table <- data.frame(matrix(0, nrow=0, ncol=12))
colnames(bmr_table) <- c("task_id", "learner_id", "resampling_id", "iters", "acc_train",  "acc_test", "bacc_train", "bacc_test", "auc_train", "auc_test", "bbrier_train", "bbrier_test")

measure = msr("classif.acc")

measures = list( 
  msr("classif.acc", id = "acc_train", predict_sets = "train"), 
  msr("classif.acc", id = "acc_test", predict_sets = "test"), 
  msr("classif.bacc", id = "bacc_train", predict_sets =  "train"), 
  msr("classif.bacc", id = "bacc_test", predict_sets =  "test"),
  msr("classif.auc", id = "auc_train", predict_sets =  "train"), 
  msr("classif.auc", id = "auc_test", predict_sets =  "test"),
  msr("classif.bbrier", id = "bbrier_train", predict_sets =  "train"),
  msr("classif.bbrier", id = "bbrier_test", predict_sets =  "test")
)

#  Load Tunning ---
for(i in c(1:length(list.clas.task))){
  list.clas.task[[i]]$print() 
  
  terminator = trm("evals", n_evals = 30)
  tuner = tnr("random_search")
  
  # #lrn("classif.xgboost")$param_set
  # xgb <- lrn("classif.xgboost")
  # xgb$param_set
  # # classif.xgboost ----
  # search_space = ps(
  #   s = p_dbl(lower = 0.001, upper = 0.1),
  #   alpha = p_dbl(lower = 0, upper = 1)
  # )
  # glmnet <- AutoTuner$new(
  #   learner = lrn("surv.glmnet", id="glmnet"),
  #   resampling = rsmp("holdout"),
  #   measure = measure,
  #   search_space = search_space,
  #   terminator = terminator,
  #   tuner = tuner
  # )
  
  
  learners = c("classif.xgboost", "classif.rfsrc", "classif.naive_bayes", "classif.svm", "classif.kknn")
  learners = lapply(learners, lrn,  predict_type = "prob", predict_sets = c("train", "test"))

  # resampling = rsmp("repeated_cv", folds = 3, repeats = 10)  
   resampling = rsmp("custom")
   resampling$instantiate(list.clas.task[[i]],
                          train = list(c(1:nrow(data))[data$dataset == "TCGA-KIRC"]),
                          test = list(c(1:nrow(data))[data$dataset == "ICGC"])
   )  
  
  # create a benchmark design object
  set.seed(1)
  design = benchmark_grid(list.clas.task[i], 
                          learners,
                          resampling = resampling
                          #resampling = rsmp("holdout")
                          #rsmp("repeated_cv", folds = 5, repeats = 20)
  )
  try({
    bmr = benchmark(design) 
    bmr = bmr$aggregate(measures) %>%
      dplyr::select(colnames(bmr_table)) %>%
      as.data.frame()
    
    print(bmr)
    
    bmr_table <- rbind(bmr_table, bmr)
  }, silent = T)
  bmr <- NULL
}


bmr_table <- bmr_table %>%
  dplyr::arrange(desc(auc_test))

## Confusion Matrix
set.seed(1)

for (i in c(1:length(learners))){
  learners[[i]]$print()
  
  learners[[i]]$train(list.clas.task[[1]], c(1:nrow(data))[data$dataset == "TCGA-KIRC"])
}

pred <- list()
for (i in c(1:length(learners))){
  a <- learners[[i]]$predict(list.clas.task[[1]], c(1:nrow(data))[data$dataset == "ICGC"])
  
  pred[[learners[[i]][["id"]]]] <- a
}

for (i in c(1:length(pred))){
  print(paste("Confusion Matrix of", learners[[i]][["id"]]))
  print(pred[[i]][["confusion"]])
}
