# generate data
set.seed(1)
v <- c(1,1)                                         # vector in direction of axis
X <- t(replicate(4, runif(1, -1, 1) * v))           # generate points on axis

## plot 1: add labels parallel to axis

plot(X, pch=16, asp=1, xlim=c(-1,1), ylim=c(-1,1))  
abline(h=0, v=0, col="grey")
abline(0, 1)    

# add point coords as labels
lbls <- apply(X, 1, function(x) paste(round(x,1), collapse=","))
add_axes_labels(v, X, lbls, cex=.7)


## plot 2: add labels perpendicular to axis / without rotation

plot(X, pch=16, asp=1, xlim=c(-1,1), ylim=c(-1,1))  
abline(h=0, v=0, col="grey")
abline(0, 1)  

# add projected values on vector as labels 
vn <- v / sum(v^2)^.5                               # normalize vector for projection             
lbls <- X %*% v                                     # project points on vector
add_axes_labels(v, X, round(lbls,1), cex=.7, 
                below=FALSE, col="blue", ortho=TRUE)
add_axes_labels(v, X, round(lbls, 1), cex=.7, srt=0, 
                line = .5, col="grey")

