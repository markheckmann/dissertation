# generate data
set.seed(1)
v <- c(-1,1)                                    # vector in direction of axis
X <- t(replicate(4, runif(1, -1, 1) * v))         # generate points on axis

## plot 1: add labels parallel to axis

plot(NULL, pch=16, xlim=c(-1,1), ylim=c(-1,1), asp=1)  
abline(h=0, v=0, col="grey")
abline(0, v[2]/v[1], grey(.3))   

#add_tick_marks 

lu = c(1,1)
n <- get_2d_normal_vector(v, upwards=TRUE)
m = scale_vec_length_to_mm(n, mm=1)
U = sweep(X, 2, lu[1]*m, "+")
L = sweep(X, 2, lu[1]*m, "-")
i = 1:4
segments(L[i, 1], L[i, 2], U[i, 1], U[i, 2], col=grey(.3))

# add axes labels
add_axes_labels(v, X, labels=i, cex=.5, line=.5)

# v : axis direction vector