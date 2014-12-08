
library(RColorBrewer)
library(MASS)

# example 1
n <- 100
set.seed(0)
x <- rnorm(n, 3) 
y <- rnorm(n, 3)
plot(x,y)
den <- kde2d(x, y, n=100)     # estimate non-parameteric density surface via kernel smoothing
cols = set_alpha_color_value(heat.colors(12), .1)
crit_dens_image(den, col=cols, prob=.1)
crit_contour(den, n=100, col="red", drawlabels=FALSE, prob=.1)


# example 2
n <- 200
set.seed(0)
x <- rnorm(n, c(1,3)) 
y <- rnorm(n, c(1,3))
d <- data.frame(x,y, g=1:2)
plot(x,y, pch=16, cex=.7, col=1:2)
g1 <- subset(d, g==1)
g2 <- subset(d, g==2)
den1 <- kde2d(g1$x, g1$y, n=100)     # estimate non-parameteric density surface via kernel smoothing
den2 <- kde2d(g2$x, g2$y, n=100)     # estimate non-parameteric density surface via kernel smoothing

reds <- set_alpha_color_value("red", seq(.05, 1, len=20))
blues <- set_alpha_color_value("blue", seq(.05, 1, len=20))
# cols <- colorRampPalette(c("white", "darkred"))(10)
# cols <- set_alpha_color_value(cols, alpha=.5)
crit_dens_image(den1, col=blues, prob=.3)
crit_contour(den1, n=100, lty=3, drawlabels=FALSE, prob=.3)
crit_dens_image(den2, col=reds, prob=.3)
crit_contour(den2, n=100, lty=3, drawlabels=FALSE, prob=.3)
