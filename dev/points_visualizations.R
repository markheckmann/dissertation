#### Visualizations of points ####



#### +-- bagplots ####

library(aplpack)
n <- 100
set.seed(0)
x <- rnorm(n, 3) 
y <- rnorm(n, 3)
plot(x,y)
r <- compute.bagplot(x, y, factor = 3, na.rm = TRUE, dkmethod=2)
polygon(r$hull.loop, density=10, angle=45, col="#FF000050")
polygon(r$hull.bag, col="#FF000050")

#library(rainbow)
#fboxplot(data = ElNino, plot.type = "bivariate", type = "hdr")


#### +-- concentration ellipse ####

require(car)
x=rnorm(100)
y=1+.3*x+.3*rnorm(100)
plot(x,y)
d <- dataEllipse(x,y, levels=0.80, draw=F)
polygon(d, density=10, angle=45, col="#FF000050")


#### +-- ellipsoid hull ####

library(cluster)
xy <- xy.coords(x,y)
xy <- remove_distant_points(xy, .1)
exy <- ellipsoidhull(cbind(xy$x,xy$y))
lines(predict(exy), col="blue")
points(rbind(exy$loc), col = "red", cex = 3, pch = 13)











