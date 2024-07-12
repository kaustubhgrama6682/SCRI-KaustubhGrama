#load libraries
library(tidyverse)
library(caret)
library(e1071)
library(glmnet)
library(MLmetrics)
library(caretEnsemble)
library(kernlab)
library(dplyr)
library(forcats)

#data
data("GermanCredit")


#Select Variables
GermanCredit <- GermanCredit %>%
    dplyr::select(Class, Duration, Amount, Age, ResidenceDuration, NumberExistingCredits, NumberPeopleMaintenance,
                  Telephone, ForeignWorker, Housing.Rent, Housing.Own, Housing.ForFree, Property.RealEstate,
                  Property.Insurance, Property.CarOther, Property.Unknown) %>%
    dplyr::rename("EmploymentDuration" = "Duration")

#Simulate missing data for the variables Age and Employment Duration
n <- nrow(GermanCredit)
agePct <- 3
durationPct <- 7

#Generate Rows that will hold missing data
set.seed(335)
ageMissingPctRows <- sample(1:n, round(agePct/100 * n, 0))

set.seed(335)
durationMissingPctRows <- sample(1:n, round(durationPct/100 * n, 0))

#Make values NA's
GermanCredit[ageMissingPctRows, "Age"] <- NA
GermanCredit[durationMissingPctRows, "EmploymentDuration"] <- NA

summary(GermanCredit)

featurePlot(x = GermanCredit[, c("EmploymentDuration", "Age")],
            y = GermanCredit$Class,
            plot = "box")

featurePlot(x = GermanCredit[, c("EmploymentDuration", "Age")],
            y = GermanCredit$Class,
            plot = "density")

featurePlot(x = GermanCredit[, 13:16],
            y = GermanCredit$Class,
            plot = "density")


#evaluating near 0 variance
nearZeroVar(GermanCredit,freqCut = 95/5, uniqueCut = 10)
nearZeroVar(GermanCredit, freqCut = 80/20, uniqueCut = 10)

#removing foreign worker factor due to lack of variance
GermanCredit <- dplyr::select(GermanCredit, -ForeignWorker)
#merging levels 2, 3, and 4 of of Numberexistingcredits
GermanCredit$NumberExistingCredits <- fct_collapse(GermanCredit$NumberExistingCredits, "2+" = c("2", "3", "4"))
table(GermanCredit$NumberExistingCredits)


#Partition data into training and test set
set.seed(355)
trainIndex <- createDataPartition(GermanCredit$Class, p = 0.7, list = FALSE)
trainingSet <- GermanCredit[trainIndex,]
testSet <- GermanCredit[-trainIndex,]


#Use the bagging imputation method for filling in missing data
set.seed(355)
bagMissing <- preProcess(trainingSet, method = "bagImpute")
trainingSet <- predict(bagMissing, newdata = trainingSet)


#Transform to a one-hot encoding data structure
dummyModel <- dummyVars(class ~., data = trainingSet)










