
set.seed(2)
x <- rnorm(100)
y <- rnorm(100)
op <- par(mfrow=c(2,3), mar=c(3,2,3,1))

plot(x,y, main="convex hull")
convex_hull(x,y)

plot(x,y, main="colored hull")
convex_hull(x,y, col="#FF000050", density=10, angle=45, border="red")

plot(x,y, main="20% of most distant\npoints removed")
convex_hull(x,y, alpha=.2, col="#FF000050", density=10, angle=45, border="red")

plot(x,y, main="initial hull points\nremoved (peeled off)")
convex_hull(x,y, peel=T, col="#00FF0050", border="red")

plot(x,y, main="smoothed hull\n(spline interpolated)")
convex_hull(x,y, smooth=2, col="#FF000050", border="red")

plot(x,y, main="another smooting\nmethod")
convex_hull(x,y, peel=T, smooth=1, shape=1, col="#00FF0050", border="red")

par(op)
