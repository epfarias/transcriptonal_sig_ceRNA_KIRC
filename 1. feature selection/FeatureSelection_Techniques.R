set.seed(100)

run_id = "testes"
formulas = list()
stamp = as.numeric(Sys.time())

#### RFE ####
n = 33;
OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
start_time <- Sys.time()
OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting RFE RF")

ctrl <- rfeControl(functions =rfFuncs,
                   method = "cv", number = 10,
                   saveDetails = TRUE,
                   allowParallel = TRUE,
                   returnResamp = "all",
                   verbose = T)

rfProfile <- rfe(trainx[-81,], train[-81,"Class"],
                 sizes = c(1:50), ## Testar com <<50
                 rfeControl = ctrl,
                 metric = "Accuracy") ## Testar com Kappa

plot(rfProfile, type=c("g", "o"))
formulas[["RandomForestRFE"]] = OmicSelector_create_formula(predictors(rfProfile))

end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
if(debug) {save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas)}

#### Boruta ####
n = 45;
OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
start_time <- Sys.time()
OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting Boruta")

library(Boruta)
bor = Boruta(trainx[-86,], train[-86,"Class"])
formulas[["Boruta"]] = OmicSelector_create_formula(names(bor$finalDecision)[as.character(bor$finalDecision) == "Confirmed"])

end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }

#### My StepWise ####
n = 51
OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting My.stepwise")
start_time <- Sys.time()

suppressMessages(library(My.stepwise))
temp = capture.output(My.stepwise.glm(Y = "Class", colnames(trainx[-86,]), data = train[-86,], 
                                      sle = 0.05, sls = 0.05, myfamily = "binomial"))
# temp2 = temp[length(temp)-1]
# temp3 = temp[length(temp)-3]
# temp4 = temp[length(temp)-5]
temp5 = paste(temp[(length(temp)-12):length(temp)], collapse = " ")
wybrane = FALSE
for(i in 1:length(colnames(trainx)))
{
  temp6 = colnames(trainx)[i]
  wybrane[i] = grepl(temp6, temp5)
}
formulas[["Mystepwise_glm_binomial"]] = OmicSelector_create_formula(colnames(trainx)[wybrane])

end_time <- Sys.time()
saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }

#### stepAIC ####
n=53
OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting stepAIC")
start_time <- Sys.time()

temp = glm(Class ~ ., data = train[-79,], family = "binomial")
temp2 = stepAIC(temp)

formulas[["stepAIC_formulas"]] = temp2$formula

end_time <- Sys.time()
saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }

#### IteratedRFECV ####
n = 55
prefer_no_features = 10

OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting MK method (iterated RFE)")
start_time <- Sys.time()

selectedMirsCV <- OmicSelector_iteratedRFE(trainSet = train[-86,], useCV = T, classLab = 'Class', checkNFeatures = prefer_no_features)$topFeaturesPerN[[prefer_no_features]]
selectedMirsTest <- OmicSelector_iteratedRFE(trainSet = train[-86,], testSet = test, classLab = 'Class', checkNFeatures = prefer_no_features)$topFeaturesPerN[[prefer_no_features]]

formulas[["iteratedRFECV"]] = OmicSelector_create_formula(selectedMirsCV)
formulas[["iteratedRFETest"]] = OmicSelector_create_formula(selectedMirsTest)

end_time <- Sys.time()
saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }

#### LASSO ####
n = 57
OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting LASSO without SMOTE")
start_time <- Sys.time()

library("glmnet")
lasso_fit <- cv.glmnet(as.matrix(trainx[-86,]), train[-86,"Class"], family = "binomial", alpha = 1)
plot(lasso_fit)
coef <- predict(lasso_fit, s = "lambda.min", type = "nonzero")
result <- data.frame(GENE = names(as.matrix(coef(lasso_fit, s = "lambda.min"))
                                  [as.matrix(coef(lasso_fit, s = "lambda.min"))[,1]!=0, 1])[-1], 
                     SCORE = as.numeric(as.matrix(coef(lasso_fit, s = "lambda.min"))
                                        [as.matrix(coef(lasso_fit, 
                                                        s = "lambda.min"))[,1]!=0, 1])[-1])
result <- result[order(-abs(result$SCORE)),]
formulas[["LASSO"]] = OmicSelector_create_formula(as.character(result$GENE))

end_time <- Sys.time()
saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }

#### ElasticNet ####
n = 58
OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting ElasticNet without SMOTE")
start_time <- Sys.time()

library(foreach)
library(glmnet)

a <- seq(0.1, 0.9, 0.05)
search <- foreach(i = a, .combine = rbind) %dopar% {
  library(glmnet)
  cv <- cv.glmnet(as.matrix(trainx[-86,]), train[-86,"Class"], family = "binomial", nfold = 10, parallel = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(as.matrix(trainx[-86,]), train[-86,"Class"], family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
result <- data.frame(GENE = names(as.matrix(coef(md3, s = "lambda.min"))
                                  [as.matrix(coef(md3, s = "lambda.min"))[,1]!=0, 1])[-1], 
                     SCORE = as.numeric(as.matrix(coef(md3, s = "lambda.min"))
                                        [as.matrix(coef(md3, 
                                                        s = "lambda.min"))[,1]!=0, 1])[-1])
result <- result[order(-abs(result$SCORE)),]
formulas[["ElasticNet"]] = OmicSelector_create_formula(as.character(result$GENE))

end_time <- Sys.time()
saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }

#### StepLDA ####
n = 59
OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting stepLDA")
start_time <- Sys.time()

tempdb = cbind(`Class` = train[-86,"Class"], trainx[-86,])

library(caret)
slda <- train(Class ~ ., data = tempdb,
              method = "stepLDA",
              trControl = trainControl(method = "loocv"))

formulas[["stepLDA"]] = OmicSelector_create_formula(predictors(slda))

end_time <- Sys.time()
saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }

#### End ####
saveRDS(formulas, paste0("temp/featureselection_formulas",stamp,".RDS"))
