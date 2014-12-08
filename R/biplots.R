# biplots

# I need biplot functions that work well in conjunction with generalized 
# procrustes procedures. 
# set colors for all axes and points
# additional points
# visualization only

# library(OpenRepGrid)
# X <- t(getRatingLayer(boeker))

# X = UDV = U D^g D^h V^T    # SVD
# supp.points: rows treated as supplementary points

biplot_calc <- function (X, g = 1, h = 1 - g, use.rows=1:nrow(X), ...) 
{
  m = nrow(X)
  n = ncol(X)
  k = min(m,n)
  # sanity check supp.points
  
  X.red <- X[use.rows, ]    # reduced X, without supplementary points
  #X <- center(x, ...)
  #X <- normalize(X, ...)
  
  dec <- svd(X.red)
  U <- dec$u
  d <- dec$d
  V <- dec$v
  D <- diag(d, ncol=k, nrow=k)
  D.g <- diag(d^g, ncol=k, nrow=k)
  D.h <- diag(d^h, ncol=k, nrow=k)
  
  # left and right matrix
  L <- U %*% D.g
  R <- D.h %*% t(V)
  
  list(X=X, L=L, R=R, U=U, D=D, V=V, g=g, h=h)
}

#biplot_calc(M)

# library(calibrate)
# 
# x <- rnorm(20,1)
# y <- rnorm(20,1)
# x <- x - mean(x)
# y <- y - mean(y)
# z <- x + y
# b <- c(1,1)
# plot(x,y,asp=1,pch=19)
# tm<-seq(-2,2,by=0.5)
# Calibrate.z <- calibrate(b,z,tm,cbind(x,y), axislab="Z",graphics=TRUE)



#' Visualization workhorse for biplots.
#' 
#' @param L      Matrix with row points.
#' @param R      Matrix with variable vectors.
#' @return      NULL.
#' @export
#' @keywords      internal
#'
biplot_visualize <- function(L, R) {
  
  
}






