library(readr)
library(dplyr)
library(tidyr)
library(randomForest)
library(cluster)
library(factoextra)
library(party)
library(caret)
library(here)


############## V Do this always V ##############
FScores <- numeric()
FAverages <- numeric()
ExampleSet <- read_csv(file.choose())#("D:/Documents/Huffman Lab/Fluorescence Spectrometry/Procedures and Data Workup Accessories/Ben's Clustering/AllData Base and Trunc.csv")
ExampleSet<-ExampleSet[!(is.na(ExampleSet$type)), ]
# summary(ExampleSet)
Temp<-0
# Temp<-select(ExampleSet, -PartNum, -type, -Date, -OriginalName)
Temp<-select(ExampleSet, -PartNum, -Variety, -Date, -OriginalName)
Temp<-as.data.frame(Temp)
# Temp$Variety = as.character(Temp$Variety)
# Temp$Variety = factor(Temp$Variety)
Temp$type = as.character(Temp$type)
Temp$type = factor(Temp$type)
# smp_size <- floor(0.75 * nrow(Temp))
# train_ind <- sample(seq_len(nrow(Temp)), size = smp_size) #Gives correct number of rows, but randomly selected.
# train <- Temp[train_ind, ]
# test <- Temp[-train_ind, ]
# set.seed(32432)
############## ^ Do this always ^ ##############

############## V Do this for looping V ##############
# NumFeatures <- (function(p,x) floor(sqrt(p)/3*x))(ncol(Temp),1:12)
NumFeatures <- c(100)#(5, 22, floor(sqrt(ncol(Temp))), 100, 500, 1000, 1524)
iterations = 100
FAverages <- numeric(length(NumFeatures))
# FScores <- matrix(NA,nrow=iterations, ncol=length(NumFeatures))
FScores <- numeric(iterations)
SeedsIndex <- matrix(NA,nrow=626,ncol=length(NumFeatures)*iterations)
SeedsCForest <- matrix(NA,nrow=626,ncol=length(NumFeatures)*iterations)

# NumTrees <- c(50,100,200,500,1000)
# iterations = 100
# FAverages <- numeric(length(NumTrees))
# FScores <- matrix(NA,nrow=iterations, ncol=length(NumTrees))
# SeedsIndex <- matrix(NA,nrow=626,ncol=length(NumFeatures)*iterations)
# SeedsCForest <- matrix(NA,nrow=626,ncol=length(NumFeatures)*iterations)

StartTime <- Sys.time()
# for (SubsetVars in NumFeatures) {
# for (j in seq_along(NumTrees)) {
for (j in seq_along(NumFeatures)) {
  for (i in 1:iterations) {
    SeedsIndex[,i*j] <- .Random.seed
    # index_2 <- createDataPartition(Temp$Variety, p = 0.75, list = FALSE) #should randomly select 75% of data to train and test, but selects ~80%; does get ~75% of each species though.
    index_2 <- createDataPartition(Temp$type, p = 0.75, list = FALSE)
    train <- Temp[index_2, ]
    test  <- Temp[-index_2, ]
    SeedsCForest[,i*j] <- .Random.seed
    # output.cforest<-cforest(Variety ~ ., data = train, controls=cforest_unbiased(ntree=500,mtry=SubsetVars))
    # output.cforest<-cforest(Variety ~ ., data = train, controls=cforest_unbiased(ntree=500,mtry=NumFeatures[j]))
    output.cforest<-cforest(type ~ ., data = train, controls=cforest_unbiased(ntree=100,mtry=NumFeatures[j]))
    # output.cforest<-cforest(Variety ~ ., data = train, controls=cforest_unbiased(ntree=NumTrees[j],mtry=5))
    predictions<-0
    predictions<-predict(output.cforest, newdata = test, type = "response", OOB = TRUE)
    # varieties<-test$Variety
    # results <- confusionMatrix(varieties, predictions)
    types<-test$type
    results <- confusionMatrix(types, predictions)
    resultsDF <- data.frame(cbind(t(results$table)))
    resultsStats <- data.frame(cbind(t(results$overall)))
    resultsClasses <- data.frame(cbind(t(results$byClass)))
    resultsClasses[is.na(resultsClasses)] <- 0
    # FScores[[length(FScores)+1]] <- mean(unlist(resultsClasses[7,]), na.rm = TRUE)
    # FScores[i,j] <- mean(unlist(resultsClasses[7,]), na.rm = TRUE)
    results <- confusionMatrix(types, predictions, positive="organic")
    resultsClasses <- data.frame(cbind(t(results$byClass)))
    F1org <- resultsClasses[7]
    results <- confusionMatrix(types, predictions, positive="inorganic")
    resultsClasses <- data.frame(cbind(t(results$byClass)))
    F1inorg <- resultsClasses[7]
    # FScores[i,j] <- (F1org+F1inorg)/2
    FScores[i] <- (F1org+F1inorg)/2
  }
  # FAverages[[length(FAverages)+1]] <- mean(FScores, na.rm = TRUE)
  # FAverages[j] <- mean(FScores, na.rm = TRUE)
  FAverages <- mean(unlist(FScores), na.rm = TRUE)
}
EndTime <- Sys.time()
difftime(EndTime,StartTime, units = "secs")
MaxCoordinates <- which(FScores == max(unlist(FScores)), arr.ind = TRUE)
print(paste("Max of Fscores:", max(FScores),"is located at"))
MaxCoordinates
############## ^ Do this for looping ^ ##############

# print(paste("The function execution duration was:", as.numeric(EndTime-StartTime,units = "secs"), "seconds."))
AvgF <- mean(FScores, na.rm = TRUE)
StdDevF <- sd(FScores)
MaxF <- max(FScores)
MinF <- min(FScores)
print(paste("Average and std. dev. of Fscores:", AvgF,"+-",StdDevF ))
print(paste("Max. F-score was:",MaxF))
print(paste("Min. F-score was:",MinF))


############## V Do this after looping to save best result V ##############
.Random.seed <- SeedsIndex[,as.numeric(MaxCoordinates)]
index_2 <- createDataPartition(Temp$type, p = 0.75, list = FALSE) #should randomly select 75% of data to train and test, but selects ~80%; does get ~75% of each species though.
train <- Temp[index_2, ]
test  <- Temp[-index_2, ]
.Random.seed <- SeedsCForest[,as.numeric(MaxCoordinates)]
output.cforest<-cforest(type ~ ., data = train, controls=cforest_unbiased(ntree=100,mtry=22))
predictions<-0
predictions<-predict(output.cforest, newdata = test, type = "response", OOB = TRUE)
varieties<-test$type
results <- confusionMatrix(varieties, predictions)
resultsDF <- data.frame(cbind(t(results$table)))
resultsStats <- data.frame(cbind(t(results$overall)))
resultsClasses <- data.frame(cbind(t(results$byClass)))
resultsClasses[is.na(resultsClasses)] <- 0


path_out <- choose.dir()
write.csv(resultsDF, file=file.path(path_out, "ClusterPredictions.csv") )
write.csv(resultsStats, file=file.path(path_out, "ClusterStats.csv") )
write.csv(resultsClasses, file=file.path(path_out, "ClusterClasses.csv") )
############## ^ Do this after looping to save best result ^ ##############

predictions<-as.data.frame(predictions)
varieties<-as.data.frame(varieties)
ClusterTable<-cbind(varieties,predictions)
ClusterTable<-as.data.frame(ClusterTable)
write.csv(ClusterTable, file = file.path(path_out, "TestSet.csv") )
print(ClusterTable)


############## V For 2-level classification V ##############
# Initialize data with two-level categorization included.
ExampleSet <- read_csv(file.choose())
ExampleSet<-ExampleSet[!(is.na(ExampleSet$type)), ]
ExampleSet<-ExampleSet[!(is.na(ExampleSet$Variety)), ]
Temp<-0
Temp<-select(ExampleSet, -PartNum, -Date, -OriginalName)
Temp<-as.data.frame(Temp)
Temp$type = as.character(Temp$type)
Temp$type = factor(Temp$type)
Temp$Variety = as.character(Temp$Variety)
Temp$Variety = factor(Temp$Variety)

####### V For looping 2-stage V #######
iterations <- 100
FScores <- numeric(iterations)
cascade_errors <- vector("list",iterations)
StartTime <- Sys.time()
for (i in 1:iterations) {
####### ^ For looping 2-stage ^ #######

  # Partition data and train on type of particle (high-level categorization).
  index_2 <- createDataPartition(Temp$Variety, p = 0.75, list = FALSE)
  train <- Temp[index_2, ]
  test  <- Temp[-index_2, ]
  output.cforest_type<-cforest(type ~ ., data = train %>% select(-Variety), controls=cforest_unbiased(ntree=100,mtry=100))
  
  # Partition data at "species-level" of categorization and train.
  train_organic <- train %>% filter(type == "organic") %>% mutate(Variety = factor(Variety))
  train_inorganic <- train %>% filter(type == "inorganic") %>% mutate(Variety = factor(Variety))
  output.cforest_organic<-cforest(Variety ~ ., data = train_organic %>% select(-type), controls=cforest_unbiased(ntree=100,mtry=100))
  output.cforest_inorganic<-cforest(Variety ~ ., data = train_inorganic %>% select(-type), controls=cforest_unbiased(ntree=100,mtry=100))
  
  # Predict test set sequentially (dependent on previous predictions)
  test$predictions_type<-predict(output.cforest_type, newdata = test, type = "response", OOB = TRUE)
  test$predictions_Variety<- NA_character_
  organic_idx <- which(test$predictions_type == "organic")
  inorganic_idx <- which(test$predictions_type == "inorganic")
  test$predictions_Variety[organic_idx] <- as.character(predict(output.cforest_organic, newdata = test[organic_idx, ], type = "response", OOB = TRUE))
  test$predictions_Variety[inorganic_idx] <- as.character(predict(output.cforest_inorganic, newdata = test[inorganic_idx, ], type = "response", OOB = TRUE))
  test$predictions_Variety <- as.factor(test$predictions_Variety)
  
  # Determining cascade error (category-level failure vs. species-level failure).
  cascade_errors[[i]] <- test %>% mutate(
    type_correct = (type == predictions_type),
    Variety_correct = (Variety == predictions_Variety),
    Error_Type = factor(case_when(
      type_correct & Variety_correct ~ "No error",
      type_correct & !Variety_correct ~ "Variety-level failure",
      !type_correct ~ "Category-level failure"
    ), levels = c("No error", "Variety-level failure", "Category-level failure"))
  ) %>%
  group_by(Error_Type, .drop = FALSE) %>% summarise(
    Ratio = (n() / nrow(test)),
    .groups = "drop"
  )
  
  
  # total_test_samples <- nrow(cascade_set)
  # cascade_errors[[i]] <- cascade_set %>% group_by(Error_Type) %>% summarise(
  #   #Count = n(),
  #   Ratio = (n() / total_test_samples),
  #   .groups = "drop"
  # )
  
  
  
  # Evaluate performance.
  results_type <- confusionMatrix(test$predictions_type, test$type)
  results_type_DF <- data.frame(cbind(t(results_type$table)))
  results_type_Stats <- data.frame(cbind(t(results_type$overall)))
  results_type_Classes <- data.frame(cbind(t(results_type$byClass)))
  results_type_Classes[is.na(results_type_Classes)] <- 0
  
  results_Variety <- confusionMatrix(test$predictions_Variety, test$Variety)
  results_Variety_DF <- data.frame(cbind(t(results_Variety$table)))
  results_Variety_Stats <- data.frame(cbind(t(results_Variety$overall)))
  results_Variety_Classes <- data.frame(cbind(t(results_Variety$byClass)))
  results_Variety_Classes[is.na(results_Variety_Classes)] <- 0
  
####### V For looping 2-stage V #######
  FScores[i] <- mean(unlist(results_Variety_Classes[7,]), na.rm = TRUE)
}
EndTime <- Sys.time()
difftime(EndTime,StartTime, units = "secs")
Final_FScores <- mean(FScores)
Final_cascade_errors <- bind_rows(cascade_errors) #%>% mean(unlist(cascade_errors$Ratio))
final_averages <- Final_cascade_errors %>% group_by(Error_Type) %>% summarise(
  Average_Ratio = mean(Ratio),
  Std_Dev_Ratio = sd(Ratio),
  .groups = "drop"
)
####### ^ For looping 2-stage ^ #######

# Record performance.
path_out <- choose.dir()

write.csv(results_type_DF, file=file.path(path_out, "ClusterPredictions_type.csv") )
write.csv(results_type_Stats, file=file.path(path_out, "ClusterStats_type.csv") )
write.csv(results_type_Classes, file=file.path(path_out, "ClusterClasses_type.csv") )

write.csv(results_Variety_DF, file=file.path(path_out, "ClusterPredictions_Variety.csv") )
write.csv(results_Variety_Stats, file=file.path(path_out, "ClusterStats_Variety.csv") )
write.csv(results_Variety_Classes, file=file.path(path_out, "ClusterClasses_Variety.csv") )

############## ^ For 2-level classification ^ ##############
