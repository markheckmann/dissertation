


#### prepare data ####

# X: n x p
n = 10
p = 4
scale= FALSE
center = TRUE

shift=0 #-mx#*1.2               # shift axes by value in usr coordinates (TODO value for each axis and in original coords)
zero=F                          # draw zero axis value?
i = 1:3  # Points
j = NULL  # axes
dim =c(1,2)

# generate

set.seed(2)
X = t(replicate(n, rnorm(p, mean=c(1,1), sd=c(1,1))))   # n x p      
#colnames(X) <- paste0("var", 1:p, " pole 1")
rownames(X) <- paste0("case", 1:n)
Xtrans = scale(X, scale=scale, center=center)
# AB = cbind(a=rep(100,p), b=rep(.1,p))   # matrix with linear transformation coefficients for each variable


#### +-- init and prepare biplot object ####

#for (a in seq(0, pi/10, len=10)) {
  
x <- biplot_obj_init(X, Xtrans)                   # add X and Xtrans to list
x <- biplot_obj_add_biplot_coords(x, g=1, h=0)    # get biplot coords
x <- biplot_obj_add_labels(x, type=1)             # add labels for rows and columns
x$variable.type = "bipolar"                       # unipolar or bipolar items
x <- prep_biplot_object(x)                        # add AB matrix of lienar transformation, mx value etc.
Q = get_rotation_for_2d(p, 1:2, alpha=a, clockwise=TRUE)
x$Q=Q                                             # rotation matrix for whole plot
x$reflect=c(F,F)                                  # reflect axis?
x$dim = dim

#### +-- set up plot ####
  
set_up_plot(x)
add_circumference(x)

#### +-- add prinicpa. axes ####
#add_principal_axes(rads=c(0,x$mx), col="grey")

#### +-- draw biplot axes ####
add_biplot_axes(x, dim=dim, col=grey(.7), j=j, shift=shift, lty=1)
add_axes_calibration(x, mx=x$mx, dim=dim, col=grey(.7), j=j, cex=.6, below=F, srt=0, 
                     ortho=FALSE, shift=shift, zero=zero)
add_biplot_axes_extensions(x, dim=dim, col=grey(0), j=j, shift=shift, lty=1)

#### +-- add biplot vectors ####
add_biplot_vectors(x, dim=dim, j=j, col=1:2, lwd=3)

#### +-- Add axes projections ####
add_point_axes_projections(x, i=1, j=j, dim=dim, col="blue", lty=2, shift=0)

#### +-- draw row points #####
add_rowpoints(x, i=i, dim=dim, pch=15, col="green", use.depth=F, use.alpha=F)

#### +-- label row points ####
add_rowpoint_labels(x, dim=dim, i=i, method=2, cex=.7, inside=T, use.depth=F, use.alpha=F)

#### +-- labels ####
add_axes_pole_label(x, dim=dim, j=j, shift=shift, cex=.7, col="blue", col2="grey")


#}




#L[, 1:2] %*% t(R[, 1:2])

# settings
# set <- list()
# set$rad.outer <- list(col=grey(.98), 
#                       border=1, lty=1,
#                       lwd=1, nv=100)
# set$plot <- list(type="circle",
#                  todo=NA)