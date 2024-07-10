#package
install.packages('caret')
library('caret')
install.packages('AppliedPredictiveModeling')
library(AppliedPredictiveModeling)
install.packages('mlbench')
library(mlbench)
install.packages('earth')
library(earth)

#Visualizations--------------------------
str(iris)


transparentTheme(trans = 0.4)

#Scatter plot matrix
featurePlot(x= iris[, 1:4],
            y = iris$Species,
            plot = "pairs",
            ##add a key at the top
            auto.key = list(columns = 3))



#Scatterplot Matrix with Ellipses
featurePlot(x = iris[, 1:4],
            y = iris$Species,
            plot = "ellipse",
            #Add a key at the top
            auto.key = list(columns = 3))

#Overlayed density plots
transparentTheme(trans = 0.9)
featurePlot(x = iris[, 1:4],
            y = iris$Species,
            plot = "density", 
            #Pass in options to xyplot() to make it prettier
            scales = list(x = list(relation="free"),
                          y = list(relation="free")),
            adjust = 1.5,
            pch = "|",
            layout = c(4,1),
            auto.key = list(columns = 3))

#Box plots
featurePlot(x = iris[, 1:4],
            y = iris$Species,
            plot = "box", 
            ##Pass in options to bwplto()
            scales = list(y = list(relation="free"),
                          x = list(rot = 90)),
            layout = c(4,1),
            auto.key = list(columns = 2))

#Scatter plots
data(BostonHousing)
regVar <- c("age", "lstat", "tax")
str(BostonHousing[, regVar])


theme1 <- trellis.par.get()
theme1$plot.symbol$col = rgb(0.2,0.2,0.2,0.4)
theme1$plot.symbol$pch = 16
theme1$plot.line$col = rgb(1, 0, 0, 0.7)
theme1$plot.line$lwd <- 2
trellis.par.set(theme1)
featurePlot(x = BostonHousing[, regVar],
            y = BostonHousing$medv,
            plot = "scatter",
            layout = c(3,1))


#Pre-Processing------------------------------------
data(etitanic)
head(model.matrix(survived ~., data = etitanic))

dummies <- dummyVars(survived ~., data = etitanic)
head(predict(dummies, newdata = etitanic))

data(mdrr)
data.frame(table(mdrrDescr$nR11))
















