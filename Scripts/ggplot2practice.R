install.packages("ggplot2")
library("ggplot2")

#read in dataset
data(iris)


# Plot
IrisPlot <- ggplot(iris, aes(Sepal.Length, Petal.Length, colour=Species)) + geom_point() + geom_smooth(lm)
print(IrisPlot  + labs(x = "Sepal length (cm)", y = "Petal length (cm)") + ggtitle("Petal and sepal length of iris"))

# add text
IrisPlot + annotate("text", x = 6, y = 5, label = "text")

# add reIrisPloteat
IrisPlot + annotate("text", x = 4:6, y = 5:7, label = "text")

# highlight an area
IrisPlot + annotate("rect", xmin = 5, xmax = 7, ymin = 4, ymax = 6, alpha = .5)

# segment
IrisPlot + annotate("segment", x = 5, xend = 7, y = 4, yend = 5, colour = "black")


bp <- ggplot(PlantGrowth, aes(y = group, x = weight)) + geom_point()
bp
ggplot(iris, aes(Sepal.Length, Petal.Length, colour=Species)) + geom_point(shape=1) + geom_smooth(method=lm)