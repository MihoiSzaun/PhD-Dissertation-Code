###### Set of commands that worked for me, rev. 3 (Maxwell) ########
library(readr)
library(dplyr)
library(tidyr)
library(randomForest)
library(cluster)
library(factoextra)
library(party)
library(caret)
library(here)

############ RF for Raman data... ###############
ExampleSet <- read_csv(file.choose())
ExampleSet<-ExampleSet[!(is.na(ExampleSet$Type)), ]
# summary(ExampleSet)
Temp<-0

############ ...(by "Species" column). ###############
Temp<-select(ExampleSet, -PartNum, -Identity, -Type, -Date, -OriginalRow)
Temp<-as.data.frame(Temp)
Temp$Species = as.character(Temp$Species)
Temp$Species = factor(Temp$Species)
index_2 <- createDataPartition(Temp$Species, p = 0.75, list = FALSE)
train <- Temp[index_2, ]
test  <- Temp[-index_2, ]
set.seed(32432)
output.cforest<-cforest(Species ~ ., data = train, controls=cforest_unbiased(ntree=1000,mtry=5))
predictions<-0
predictions<-predict(output.cforest, newdata = test, type = "response", OOB = TRUE)
Species<-test$Species
path_out <- choose.dir()
results <- confusionMatrix(Species, predictions)
resultsDF <- data.frame(cbind(t(results$table)))
resultsStats <- data.frame(cbind(t(results$overall)))
resultsClasses <- data.frame(cbind(t(results$byClass)))
write.csv(resultsDF, file=file.path(path_out, "ClusterPredictions.csv") )
write.csv(resultsStats, file=file.path(path_out, "ClusterStats.csv") )
write.csv(resultsClasses, file=file.path(path_out, "ClusterClasses.csv") )
predictions<-as.data.frame(predictions)
Species<-as.data.frame(Species)
ClusterTable<-cbind(Species,predictions)
ClusterTable<-as.data.frame(ClusterTable)
write.csv(ClusterTable, file = file.path(path_out, "TestSet.csv") )
print(ClusterTable)

############ ...(by "Type" column). ###############
Temp<-select(ExampleSet, -PartNum, -Species, -Identity, -Date, -OriginalRow)
Temp<-as.data.frame(Temp)
Temp$Type = as.character(Temp$Type)
Temp$Type = factor(Temp$Type)

# NumFeatures <- c(5,(sqrt(ncol(Temp)-1))/2,sqrt(ncol(Temp)-1),3*sqrt(ncol(Temp)-1),6*sqrt(ncol(Temp)-1),(ncol(Temp)-1)/3 )
# iterations = 100
# FAverages <- numeric(length(NumFeatures))
# FScores <- matrix(NA,nrow=iterations, ncol=length(NumFeatures))
# SeedsIndex <- matrix(NA,nrow=626,ncol=length(NumFeatures)*iterations)
# SeedsCForest <- matrix(NA,nrow=626,ncol=length(NumFeatures)*iterations)

NumTrees <- c(100)
iterations = 100
FAverages <- numeric(length(NumTrees))
FScores <- matrix(NA,nrow=iterations, ncol=length(NumTrees))
# SeedsIndex <- matrix(NA,nrow=626,ncol=length(NumTrees)*iterations)
# SeedsCForest <- matrix(NA,nrow=626,ncol=length(NumTrees)*iterations)

StartTime <- Sys.time()
# for (SubsetVars in NumFeatures) {
for (j in seq_along(NumTrees)) {
# for (j in seq_along(NumFeatures)) {
  for (i in 1:iterations) {
    # SeedsIndex[,i*j] <- .Random.seed
    index_2 <- createDataPartition(Temp$Type, p = 0.75, list = FALSE)
    train <- Temp[index_2, ]
    test  <- Temp[-index_2, ]
    # set.seed(32432)
    # SeedsCForest[,i*j] <- .Random.seed
    # output.cforest<-cforest(Type ~ ., data = train, controls=cforest_unbiased(ntree=1000,mtry=5))
    # output.cforest<-cforest(Type ~ ., data = train, controls=cforest_unbiased(ntree=500,mtry=NumFeatures[j]))
    output.cforest<-cforest(Type ~ ., data = train, controls=cforest_unbiased(ntree=NumTrees[j],mtry=100))
    predictions<-0
    predictions<-predict(output.cforest, newdata = test, type = "response", OOB = TRUE)
    Type<-test$Type
    results <- confusionMatrix(Type, predictions)
    resultsDF <- data.frame(cbind(t(results$table)))
    resultsStats <- data.frame(cbind(t(results$overall)))
    resultsClasses <- data.frame(cbind(t(results$byClass)))
    FScores[i,j] <- mean(unlist(resultsClasses[7,]), na.rm = TRUE)
  }
  # FAverages[[length(FAverages)+1]] <- mean(FScores, na.rm = TRUE)
  FAverages[j] <- mean(FScores, na.rm = TRUE)
}
EndTime <- Sys.time()
difftime(EndTime,StartTime, units = "secs")
MaxCoordinates <- which(FScores == max(FScores), arr.ind = TRUE)
print(paste("Max of Fscores:", max(FScores),"is located at"))
MaxCoordinates

AvgF <- mean(FScores, na.rm = TRUE)
StdDevF <- sd(FScores)
MaxF <- max(FScores)
MinF <- min(FScores)
print(paste("Average and std. dev. of Fscores:", AvgF,"+-",StdDevF ))
print(paste("Max. F-score was:",MaxF))
print(paste("Min. F-score was:",MinF))

path_out <- choose.dir()
write.csv(resultsDF, file=file.path(path_out, "ClusterPredictions.csv") )
write.csv(resultsStats, file=file.path(path_out, "ClusterStats.csv") )
write.csv(resultsClasses, file=file.path(path_out, "ClusterClasses.csv") )
predictions<-as.data.frame(predictions)
Type<-as.data.frame(Type)
ClusterTable<-cbind(Type,predictions)
ClusterTable<-as.data.frame(ClusterTable)
write.csv(ClusterTable, file = file.path(path_out, "TestSet.csv") )
print(ClusterTable)

############ ...(by "Identity" column). ###############
Temp<-select(ExampleSet, -PartNum, -Species, -Type, -Date, -OriginalRow)
Temp<-as.data.frame(Temp)
Temp$Identity = as.character(Temp$Identity)
Temp$Identity = factor(Temp$Identity)
index_2 <- createDataPartition(Temp$Identity, p = 0.75, list = FALSE)
train <- Temp[index_2, ]
test  <- Temp[-index_2, ]
set.seed(32432)
output.cforest<-cforest(Identity ~ ., data = train, controls=cforest_unbiased(ntree=1000,mtry=5))
predictions<-0
predictions<-predict(output.cforest, newdata = test, type = "response", OOB = TRUE)
Identity<-test$Identity
path_out <- choose.dir()
results <- confusionMatrix(Identity, predictions)
resultsDF <- data.frame(cbind(t(results$table)))
resultsStats <- data.frame(cbind(t(results$overall)))
resultsClasses <- data.frame(cbind(t(results$byClass)))
write.csv(resultsDF, file=file.path(path_out, "ClusterPredictions.csv") )
write.csv(resultsStats, file=file.path(path_out, "ClusterStats.csv") )
write.csv(resultsClasses, file=file.path(path_out, "ClusterClasses.csv") )
predictions<-as.data.frame(predictions)
Identity<-as.data.frame(Identity)
ClusterTable<-cbind(Identity,predictions)
ClusterTable<-as.data.frame(ClusterTable)
write.csv(ClusterTable, file = file.path(path_out, "TestSet.csv") )
print(ClusterTable)