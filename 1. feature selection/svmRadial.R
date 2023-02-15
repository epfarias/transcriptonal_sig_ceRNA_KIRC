set.seed(100)

source("OmicSelector_Load_DataMix.R")

dane = OmicSelector_load_datamix(); 
train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

library(readr)
featureselection_formulas_final <- read_csv("featureselection_formulas_final.csv")

#### tabela de resultados ####
resul_tab = data.frame()
resul_tab <- featureselection_formulas_final$name
resul_tab <- as.data.frame(resul_tab)

temptrain = train[-86,]
temptrainold = temptrain
fit_on = list(rs1 = 1:nrow(temptrain))
pred_on = list(rs1 = (nrow(temptrain)+1):((nrow(temptrain))+nrow(test)))
temptrain = rbind(temptrain,test)

train_control <- trainControl(method="cv", index= fit_on, indexOut = pred_on, indexFinal = fit_on[[1]], verboseIter = TRUE,
                              classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE) 

model1 = caret::train(as.formula(as.matrix(featureselection_formulas_final[1,2])), data=temptrain, trControl=train_control, method="svmRadial", tuneLength = 10)
print(model1$finalModel)

t1_roc = roc(temptrainold$Class ~ predict(model1, newdata = temptrainold , type = "prob")[,2])
resul_tab[1,"svmRadial_train_ROCAUC"] = t1_roc$auc
resul_tab[1,"svmRadial_train_ROCAUC_lower95CI"] = as.character(ci(t1_roc))[1]
resul_tab[1,"svmRadial_train_ROCAUC_upper95CI"] = as.character(ci(t1_roc))[3]

t1 = caret::confusionMatrix(predict(model1, newdata = temptrainold), as.factor(temptrainold$Class), positive = "M1")
resul_tab[1,"svmRadial_train_Accuracy"] = t1$overall["Accuracy"]
resul_tab[1,"svmRadial_train_Sensitivity"] = t1$byClass["Sensitivity"]
resul_tab[1,"svmRadial_train_Specificity"] = t1$byClass["Specificity"]

v1 = caret::confusionMatrix(predict(model1, newdata = test), as.factor(test$Class), positive = "M1")
resul_tab[1,"svmRadial_test_Accuracy"]  = v1$overall["Accuracy"]
resul_tab[1,"svmRadial_test_Sensitivity"]  = v1$byClass["Sensitivity"]
resul_tab[1,"svmRadial_test_Specificity"]  = v1$byClass["Specificity"]

v1 = caret::confusionMatrix(predict(model1, newdata = valid), as.factor(valid$Class), positive = "M1")
resul_tab[1,"svmRadial_valid_Accuracy"]  = v1$overall["Accuracy"]
resul_tab[1,"svmRadial_valid_Sensitivity"]= v1$byClass["Sensitivity"]
resul_tab[1,"svmRadial_valid_Specificity"] = v1$byClass["Specificity"]

model2 = caret::train(as.formula(as.matrix(featureselection_formulas_final[2,2])), data=temptrain, trControl=train_control, method="svmRadial", tuneLength = 10)
print(model2$finalModel)

t2_roc = roc(temptrainold$Class ~ predict(model2, newdata = temptrainold , type = "prob")[,2])
resul_tab[2,"svmRadial_train_ROCAUC"] = t2_roc$auc
resul_tab[2,"svmRadial_train_ROCAUC_lower95CI"] = as.character(ci(t2_roc))[1]
resul_tab[2,"svmRadial_train_ROCAUC_upper95CI"] = as.character(ci(t2_roc))[3]

t2 = caret::confusionMatrix(predict(model2, newdata = temptrainold), as.factor(temptrainold$Class), positive = "M1")
resul_tab[2,"svmRadial_train_Accuracy"] = t2$overall["Accuracy"]
resul_tab[2,"svmRadial_train_Sensitivity"] = t2$byClass["Sensitivity"]
resul_tab[2,"svmRadial_train_Specificity"] = t2$byClass["Specificity"]

v2 = caret::confusionMatrix(predict(model2, newdata = test), as.factor(test$Class), positive = "M1")
resul_tab[2,"svmRadial_test_Accuracy"]  = v2$overall["Accuracy"]
resul_tab[2,"svmRadial_test_Sensitivity"]  = v2$byClass["Sensitivity"]
resul_tab[2,"svmRadial_test_Specificity"]  = v2$byClass["Specificity"]

v2 = caret::confusionMatrix(predict(model2, newdata = valid), as.factor(valid$Class), positive = "M1")
resul_tab[2,"svmRadial_valid_Accuracy"]  = v2$overall["Accuracy"]
resul_tab[2,"svmRadial_valid_Sensitivity"]= v2$byClass["Sensitivity"]
resul_tab[2,"svmRadial_valid_Specificity"] = v2$byClass["Specificity"]

model3 = caret::train(as.formula(as.matrix(featureselection_formulas_final[3,2])), data=temptrain, trControl=train_control, method="svmRadial", tuneLength = 10)
print(model3$finalModel)

t3_roc = roc(temptrainold$Class ~ predict(model3, newdata = temptrainold , type = "prob")[,2])
resul_tab[3,"svmRadial_train_ROCAUC"] = t3_roc$auc
resul_tab[3,"svmRadial_train_ROCAUC_lower95CI"] = as.character(ci(t3_roc))[1]
resul_tab[3,"svmRadial_train_ROCAUC_upper95CI"] = as.character(ci(t3_roc))[3]

t3 = caret::confusionMatrix(predict(model3, newdata = temptrainold), as.factor(temptrainold$Class), positive = "M1")
resul_tab[3,"svmRadial_train_Accuracy"] = t3$overall["Accuracy"]
resul_tab[3,"svmRadial_train_Sensitivity"] = t3$byClass["Sensitivity"]
resul_tab[3,"svmRadial_train_Specificity"] = t3$byClass["Specificity"]

v3 = caret::confusionMatrix(predict(model3, newdata = test), as.factor(test$Class), positive = "M1")
resul_tab[3,"svmRadial_test_Accuracy"]  = v3$overall["Accuracy"]
resul_tab[3,"svmRadial_test_Sensitivity"]  = v3$byClass["Sensitivity"]
resul_tab[3,"svmRadial_test_Specificity"]  = v3$byClass["Specificity"]

v3 = caret::confusionMatrix(predict(model3, newdata = valid), as.factor(valid$Class), positive = "M1")
resul_tab[3,"svmRadial_valid_Accuracy"]  = v3$overall["Accuracy"]
resul_tab[3,"svmRadial_valid_Sensitivity"]= v3$byClass["Sensitivity"]
resul_tab[3,"svmRadial_valid_Specificity"] = v3$byClass["Specificity"]

model4 = caret::train(as.formula(as.matrix(featureselection_formulas_final[4,2])), data=temptrain, trControl=train_control, method="svmRadial", tuneLength = 10)
print(model4$finalModel)

t4_roc = roc(temptrainold$Class ~ predict(model4, newdata = temptrainold , type = "prob")[,2])
resul_tab[4,"svmRadial_train_ROCAUC"] = t4_roc$auc
resul_tab[4,"svmRadial_train_ROCAUC_lower95CI"] = as.character(ci(t4_roc))[1]
resul_tab[4,"svmRadial_train_ROCAUC_upper95CI"] = as.character(ci(t4_roc))[3]

t4 = caret::confusionMatrix(predict(model4, newdata = temptrainold), as.factor(temptrainold$Class), positive = "M1")
resul_tab[4,"svmRadial_train_Accuracy"] = t4$overall["Accuracy"]
resul_tab[4,"svmRadial_train_Sensitivity"] = t4$byClass["Sensitivity"]
resul_tab[4,"svmRadial_train_Specificity"] = t4$byClass["Specificity"]

v4 = caret::confusionMatrix(predict(model4, newdata = test), as.factor(test$Class), positive = "M1")
resul_tab[4,"svmRadial_test_Accuracy"]  = v4$overall["Accuracy"]
resul_tab[4,"svmRadial_test_Sensitivity"]  = v4$byClass["Sensitivity"]
resul_tab[4,"svmRadial_test_Specificity"]  = v4$byClass["Specificity"]

v4 = caret::confusionMatrix(predict(model4, newdata = valid), as.factor(valid$Class), positive = "M1")
resul_tab[4,"svmRadial_valid_Accuracy"]  = v4$overall["Accuracy"]
resul_tab[4,"svmRadial_valid_Sensitivity"]= v4$byClass["Sensitivity"]
resul_tab[4,"svmRadial_valid_Specificity"] = v4$byClass["Specificity"]

model5 = caret::train(as.formula(as.matrix(featureselection_formulas_final[5,2])), data=temptrain, trControl=train_control, method="svmRadial", tuneLength = 10)
print(model5$finalModel)

t5_roc = roc(temptrainold$Class ~ predict(model5, newdata = temptrainold , type = "prob")[,2])
resul_tab[5,"svmRadial_train_ROCAUC"] = t5_roc$auc
resul_tab[5,"svmRadial_train_ROCAUC_lower95CI"] = as.character(ci(t5_roc))[1]
resul_tab[5,"svmRadial_train_ROCAUC_upper95CI"] = as.character(ci(t5_roc))[3]

t5 = caret::confusionMatrix(predict(model5, newdata = temptrainold), as.factor(temptrainold$Class), positive = "M1")
resul_tab[5,"svmRadial_train_Accuracy"] = t5$overall["Accuracy"]
resul_tab[5,"svmRadial_train_Sensitivity"] = t5$byClass["Sensitivity"]
resul_tab[5,"svmRadial_train_Specificity"] = t5$byClass["Specificity"]

v5 = caret::confusionMatrix(predict(model5, newdata = test), as.factor(test$Class), positive = "M1")
resul_tab[5,"svmRadial_test_Accuracy"]  = v5$overall["Accuracy"]
resul_tab[5,"svmRadial_test_Sensitivity"]  = v5$byClass["Sensitivity"]
resul_tab[5,"svmRadial_test_Specificity"]  = v5$byClass["Specificity"]

v5 = caret::confusionMatrix(predict(model5, newdata = valid), as.factor(valid$Class), positive = "M1")
resul_tab[5,"svmRadial_valid_Accuracy"]  = v5$overall["Accuracy"]
resul_tab[5,"svmRadial_valid_Sensitivity"]= v5$byClass["Sensitivity"]
resul_tab[5,"svmRadial_valid_Specificity"] = v5$byClass["Specificity"]

model6 = caret::train(as.formula(as.matrix(featureselection_formulas_final[6,2])), data=temptrain, trControl=train_control, method="svmRadial", tuneLength = 10)
print(model6$finalModel)

t6_roc = roc(temptrainold$Class ~ predict(model6, newdata = temptrainold , type = "prob")[,2])
resul_tab[6,"svmRadial_train_ROCAUC"] = t6_roc$auc
resul_tab[6,"svmRadial_train_ROCAUC_lower95CI"] = as.character(ci(t6_roc))[1]
resul_tab[6,"svmRadial_train_ROCAUC_upper95CI"] = as.character(ci(t6_roc))[3]

t6 = caret::confusionMatrix(predict(model6, newdata = temptrainold), as.factor(temptrainold$Class), positive = "M1")
resul_tab[6,"svmRadial_train_Accuracy"] = t6$overall["Accuracy"]
resul_tab[6,"svmRadial_train_Sensitivity"] = t6$byClass["Sensitivity"]
resul_tab[6,"svmRadial_train_Specificity"] = t6$byClass["Specificity"]

v6 = caret::confusionMatrix(predict(model6, newdata = test), as.factor(test$Class), positive = "M1")
resul_tab[6,"svmRadial_test_Accuracy"]  = v6$overall["Accuracy"]
resul_tab[6,"svmRadial_test_Sensitivity"]  = v6$byClass["Sensitivity"]
resul_tab[6,"svmRadial_test_Specificity"]  = v6$byClass["Specificity"]

v6 = caret::confusionMatrix(predict(model6, newdata = valid), as.factor(valid$Class), positive = "M1")
resul_tab[6,"svmRadial_valid_Accuracy"]  = v6$overall["Accuracy"]
resul_tab[6,"svmRadial_valid_Sensitivity"]= v6$byClass["Sensitivity"]
resul_tab[6,"svmRadial_valid_Specificity"] = v6$byClass["Specificity"]

model7= caret::train(as.formula(as.matrix(featureselection_formulas_final[7,2])), data=temptrain, trControl=train_control, method="svmRadial", tuneLength = 10)
print(model7$finalModel)

t7_roc = roc(temptrainold$Class ~ predict(model7, newdata = temptrainold , type = "prob")[,2])
resul_tab[7,"svmRadial_train_ROCAUC"] = t7_roc$auc
resul_tab[7,"svmRadial_train_ROCAUC_lower95CI"] = as.character(ci(t7_roc))[1]
resul_tab[7,"svmRadial_train_ROCAUC_upper95CI"] = as.character(ci(t7_roc))[3]

t7 = caret::confusionMatrix(predict(model7, newdata = temptrainold), as.factor(temptrainold$Class), positive = "M1")
resul_tab[7,"svmRadial_train_Accuracy"] = t7$overall["Accuracy"]
resul_tab[7,"svmRadial_train_Sensitivity"] = t7$byClass["Sensitivity"]
resul_tab[7,"svmRadial_train_Specificity"] = t7$byClass["Specificity"]

v7 = caret::confusionMatrix(predict(model7, newdata = test), as.factor(test$Class), positive = "M1")
resul_tab[7,"svmRadial_test_Accuracy"]  = v7$overall["Accuracy"]
resul_tab[7,"svmRadial_test_Sensitivity"]  = v7$byClass["Sensitivity"]
resul_tab[7,"svmRadial_test_Specificity"]  = v7$byClass["Specificity"]

v7 = caret::confusionMatrix(predict(model7, newdata = valid), as.factor(valid$Class), positive = "M1")
resul_tab[7,"svmRadial_valid_Accuracy"]  = v7$overall["Accuracy"]
resul_tab[7,"svmRadial_valid_Sensitivity"]= v7$byClass["Sensitivity"]
resul_tab[7,"svmRadial_valid_Specificity"] = v7$byClass["Specificity"]

model8 = caret::train(as.formula(as.matrix(featureselection_formulas_final[8,2])), data=temptrain, trControl=train_control, method="svmRadial", tuneLength = 10)
print(model8$finalModel)

t8_roc = roc(temptrainold$Class ~ predict(model8, newdata = temptrainold , type = "prob")[,2])
resul_tab[8,"svmRadial_train_ROCAUC"] = t8_roc$auc
resul_tab[8,"svmRadial_train_ROCAUC_lower95CI"] = as.character(ci(t8_roc))[1]
resul_tab[8,"svmRadial_train_ROCAUC_upper95CI"] = as.character(ci(t8_roc))[3]

t8 = caret::confusionMatrix(predict(model8, newdata = temptrainold), as.factor(temptrainold$Class), positive = "M1")
resul_tab[8,"svmRadial_train_Accuracy"] = t8$overall["Accuracy"]
resul_tab[8,"svmRadial_train_Sensitivity"] = t8$byClass["Sensitivity"]
resul_tab[8,"svmRadial_train_Specificity"] = t8$byClass["Specificity"]

v8 = caret::confusionMatrix(predict(model8, newdata = test), as.factor(test$Class), positive = "M1")
resul_tab[8,"svmRadial_test_Accuracy"]  = v8$overall["Accuracy"]
resul_tab[8,"svmRadial_test_Sensitivity"]  = v8$byClass["Sensitivity"]
resul_tab[8,"svmRadial_test_Specificity"]  = v8$byClass["Specificity"]

v8 = caret::confusionMatrix(predict(model8, newdata = valid), as.factor(valid$Class), positive = "M1")
resul_tab[8,"svmRadial_valid_Accuracy"]  = v8$overall["Accuracy"]
resul_tab[8,"svmRadial_valid_Sensitivity"]= v8$byClass["Sensitivity"]
resul_tab[8,"svmRadial_valid_Specificity"] = v8$byClass["Specificity"]

model9 = caret::train(as.formula(as.matrix(featureselection_formulas_final[9,2])), data=temptrain, trControl=train_control, method="svmRadial", tuneLength = 10)
print(model9$finalModel)

t9_roc = roc(temptrainold$Class ~ predict(model9, newdata = temptrainold , type = "prob")[,2])
resul_tab[9,"svmRadial_train_ROCAUC"] = t9_roc$auc
resul_tab[9,"svmRadial_train_ROCAUC_lower95CI"] = as.character(ci(t9_roc))[1]
resul_tab[9,"svmRadial_train_ROCAUC_upper95CI"] = as.character(ci(t9_roc))[3]

t9 = caret::confusionMatrix(predict(model9, newdata = temptrainold), as.factor(temptrainold$Class), positive = "M1")
resul_tab[9,"svmRadial_train_Accuracy"] = t9$overall["Accuracy"]
resul_tab[9,"svmRadial_train_Sensitivity"] = t9$byClass["Sensitivity"]
resul_tab[9,"svmRadial_train_Specificity"] = t9$byClass["Specificity"]

v9 = caret::confusionMatrix(predict(model9, newdata = test), as.factor(test$Class), positive = "M1")
resul_tab[9,"svmRadial_test_Accuracy"]  = v9$overall["Accuracy"]
resul_tab[9,"svmRadial_test_Sensitivity"]  = v9$byClass["Sensitivity"]
resul_tab[9,"svmRadial_test_Specificity"]  = v9$byClass["Specificity"]

v9 = caret::confusionMatrix(predict(model9, newdata = valid), as.factor(valid$Class), positive = "M1")
resul_tab[9,"svmRadial_valid_Accuracy"]  = v9$overall["Accuracy"]
resul_tab[9,"svmRadial_valid_Sensitivity"]= v9$byClass["Sensitivity"]
resul_tab[9,"svmRadial_valid_Specificity"] = v9$byClass["Specificity"]