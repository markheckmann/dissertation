# TODO: 
# - show untransformed values on axes (check for correctness)
# - individual shifts (given as mu or original values)
# - calc and viz predicitivity for each axes
# - add 3d depth for labels and points
# - alternative drawing of biplot vectors as points
# - lambda scaling
# - radial color gradient (Martin)
# - simple and fast label placement algorithm (Martin)
# - reflection and rotation


#library(UBbipl)
#library(UBFigs)
library(plotrix)

draw_axis <- function(v, o=c(0,0), col="grey") 
{
  u = v / sum(v^2)^.5   # unit vector in axis direction
  s = u * rad.outer     # stretch vector to end on circumference  
  n = ori + s           # add direction vector to origin
  segments(ori[1],  ori[2], n[1], n[2], col=col)   # segment origin pole 1
  segments(ori[1],  ori[2], -n[1], -n[2], col=col) # segment origin pole 2
}


draw_axis_label_segment <- function(v, origin=c(0,0), rads=c(0,1), col="grey", shift=0, ...) 
{
  n = get_2d_normal_vector(v, upwards=TRUE)   # get vector normal to axis 
  u = v / sum(v^2)^.5     # unit vector in axis direction
  c0 = u * rads[1]        # stretch vector to end on circumference  
  c1 = u * rads[2]        # stretch vector to end on circumference  
  n0 = origin + c0        # add direction vector to origin
  n1 = origin + c1        # add direction vector to origin
  s = shift * n
  n0 = n0 + s             # shift axis
  n1 = n1 + s             # shift axis
  segments(n0[1],  n0[2], n1[1], n1[2], col=col, ...)  # segment origin pole 1
}

# ... passed on to segments
draw_axis_label_segments <- function(v, origin=c(0,0), rads=c(0,1), col="grey", shift=0, ...) 
{
  draw_axis_label_segment(v, origin=origin, rads=rads, col=col, shift=shift, ...)
  draw_axis_label_segment(v, origin=origin, rads=rads*-1, col=col, shift=shift, ...)
}


extend_vec_to_length <- function(v, len=1) 
{
  u = v / sum(v^2)^.5   # unit vector in axis direction
  u * len               # stretch vector to multiple of u
}


# principal axes are not rotated, just a horizontal and 
# vertical line
add_principal_axes <- function(rads=c(0,1), col=1) 
{
  draw_axis_label_segments(c(1,0), rads=rads, col=col)
  draw_axis_label_segments(c(0,1), rads=rads, col=col)
}


# get nice values on biplot axis, with clipping to fit within mx radius
# mx = maximal and minimal value displayed on axis
get_biplot_mus <- function(v, mx, clip=TRUE, n=5, zero=TRUE)
{
  # get pretty values for data range, i.e. the scale values we want to display on the biplot axes 
  # we might use the original data range, but then some regions of the axis might remain empty
  # if no values exist. Hence we use the plot radius as minimal and maxmimal value.
  #
  # lambda = multiple of direction vector to get where the axis value mu lays on the axis.
  # See Gower et al., 2010, p.23. for how to determine lambda from a given mu.
  # Determine maximal displayable mu on axis, i.e. <= a value mx.   
  #       lambda * norm(v) <= mx
  #    =  mu / (v^T v) * norm(v) <= mx
  #    =  mu <= mx * (v^T v) / norm(v)
  
  mu.max = mx * (v %*% v) / norm_vec(v)     # get maximal displayable mu
  mu.max = as.vector(mu.max)                # from 1 x 1-matrix to vector
  mu.max.range = mu.max * c(-1,1)
  
  mus = pretty(mu.max.range, n=n)           # use +/- radius max mu as data range to get pretty mus for range 
  if (clip)                                 # remove mu values outside desired plot range                                        
    mus = mus[mus >= -mu.max & mus <= mu.max]
  if (!zero)
    mus = mus[mus != 0]
  mus
}


# get pretty values for data range, i.e. the scale values we want to display on the biplot axes 
# we might us the original data range, but then some regions of the axis might remain empty
# if no values exist. Hence we use the plot radius as minimal and maxmimal value.
#
# lambda = multiple of direction vector to get where the axis value mu lays on the axis.
# See Gower et al., 2010, p.23. for how to determine lambda from a given mu.
# Determine maximal displayable mu on axis, i.e. <= a value mx.   
#       lambda * norm(v) <= mx
#    =  mu / (v^T v) * norm(v) <= mx
#    =  mu <= mx * (v^T v) / norm(v)
#
get_biplot_lambda_and_values <- function(v, mx, ab=c(0,1), clip=TRUE, n=5, zero=TRUE) 
{
  mu.max = mx * (v %*% v) / norm_vec(v)     # get maximal displayable mu
  mu.max = as.vector(mu.max)                # from 1 x 1-matrix to vector
  
  f = linear_convert_x_y(a=ab[1], b=ab[2])  # build conversion function from x to mu(y) (and back)
                                            # c(0,1) defaults to identity 
  mu.max.range = mu.max * c(-1,1)
  orig.mx.range = f(y=mu.max.range)$x       # convert mus (transformed) to original data range 
  origs = pretty(orig.mx.range, n=n)        # get range for original values
  mus = f(x=origs)$y                        # convert orignal values back to transformed mus
  
  if (clip) {                               # remove mu values outside desired plot range                                        
    i = mus >= -mu.max & mus <= mu.max
    mus = mus[i]
    origs = origs[i]
  }
  if (!zero) {                              # remove zero values (for intersection of axes)
    i <- mus != 0
    mus = mus[i]
    origs = origs[i]
  }
  
  lambdas = mus / t(v) %*% v                # get location of mus on biplot axis as multiple of biplot vector (lambda). 
                                            # See Gower et al., 2010, p.23. for rationale 
  list(lambdas=lambdas,                     # multiples of v
       mus=mus,                             # transformed value at lambda*v
       origs=origs)                         # original value at lambda*v                               
}


# search pretty values in original data range

# add_biplot_axis <- function(v, ){
#   for (i in 1:nrow(R))
#     draw_axis_label_segments(R[i, 1:2], rads=c(0, rad.outer), col=col[i])
#   
# }

# mus : values to be displayed as axis marks
add_calibration_to_biplot_axis <- function(v, mus, labels=mus, i=seq_along(mus), 
                                           clip=TRUE, tick.ext=c(1,1), mm=1,
                                           shift=0) 
{
  lambdas = mus / t(v) %*% v                  # get location of mus on biplot axis as multiple of biplot vector (lambda). 
                                              # See Gower et al., 2010, p.23. for rationale
  P = cbind(lambdas) %*% v                    # get position of mus on biplot axis
  n = get_2d_normal_vector(v, upwards=TRUE)   # get vector normal to axis
  m = scale_vec_length_to_mm(n, mm=mm)        # scale vector to given plot length in mm
  
  P = sweep(P, 2, shift*n, "+")               # parallel shift of axis 
  U = sweep(P, 2, tick.ext[1]*m, "+")         # add length orthogonal to mu position for upper tick mark
  L = sweep(P, 2, tick.ext[2]*m, "-")         # lower tick mark position
  segments(L[i, 1], L[i, 2], U[i, 1], U[i, 2], col=grey(.3))      # draw tick marks
  add_axes_labels(v, P[i, ], labels[i], cex=.6, line=.5, srt = 0)    # add labels
}


add_calibration_to_biplot_axis_2 <- function(v, lambdas, labels, i=seq_along(lambdas), 
                                             clip=TRUE, tick.ext=c(1,1), mm=1,
                                             shift=0, ticks.col=grey(.3), srt=0, line=.5,
                                             ortho=FALSE, ...) 
{
  #lambdas = mus / t(v) %*% v                  # get location of mus on biplot axis as multiple of biplot vector (lambda). 
  # See Gower et al., 2010, p.23. for rationale
  P = cbind(lambdas) %*% v                    # get position of mus on biplot axis
  n = get_2d_normal_vector(v, upwards=TRUE)   # get vector normal to axis
  m = scale_vec_length_to_mm(n, mm=mm)        # scale vector to given plot length in mm
  
  P = sweep(P, 2, shift*n, "+")               # parallel shift of axis 
  U = sweep(P, 2, tick.ext[1]*m, "+")         # add length orthogonal to mu position for upper tick mark
  L = sweep(P, 2, tick.ext[2]*m, "-")         # lower tick mark position
  segments(L[i, 1], L[i, 2], U[i, 1], U[i, 2], col=ticks.col)          # draw tick marks
  add_axes_labels(v, P[i, ], labels[i], line=line, srt = srt, ortho=ortho, ...)      # add labels
}


# p: point
#
add_point_axis_projection <- function(x, v, shift=0, ...) 
{
  x = as.vector(x)                            # make sure its a column vector
  v = as.vector(v)                            # make sure its a column vector
  P = v %*% t(v) / as.vector(t(v) %*% v)      # projection matrix
  p = P %*% x                                 # projected vector
  n = get_2d_normal_vector(v, upwards=TRUE)   # get vector normal to axis 
  p = p + shift * n                           # add shift
  segments(x[1], x[2], p[1], p[2], ...)  # draw projection on vector
}
# v = R[1, 1:2]
# x = L[1, 1:2]
# add_point_axis_projection(x, v, shift=shift)
# add_point_axis_projections(x, t(R[, 1:2]), shift=shift)

# vectors in columns of matrix
# l: row points (vector or matrix)
# a: axes (vector or matrix)
add_point_axis_projections <- function(l, r, shift=0, col=1, ...) 
{
  L = matrix(l, ncol=2)
  R = matrix(r, ncol=2)
  #L = as.matrix(l)
  #R = as.matrix(r)

  n = nrow(L)     # number of points (given as matrix of column vectors)
  p = nrow(R)     # number of vectors to project on (given as matrix of column vectors)
  
  if (length(col) < p)
    col = rep(col, length.out=p)
  if (length(shift) < p)
    shift = rep(shift, p)
  
  for (i in 1:n) 
    for (j in 1:p) {
      add_point_axis_projection(L[i, ], R[j, ], shift=shift[j], col=col[j], ...)   
    }
}


# Linear transformation y = a + bx to convert x to y
#
# Given the data for a and x get the lineare transformation  
# coefficients
# example
# x <- rnorm(10) 
# s <- scale(x)
# ab = determine_linear_trans_from_data(x, s)
# cbind(x*ab[2] + ab[1], s)
#
determine_linear_trans_from_data <- function(x, y) 
{
  # y = a + b*x
  i = c(which.min(x), which.max(x))   # get x values which are not identical
  b = diff(y[i]) / diff(x[i])         # b = (y.1 - y.2) / (x.1 - x.2)
  a = y[1] - b*x[1]                   # a = y.1 - x.1*b
  c(a=a, b=b)
}


# convert between original data and mu value 
# mu is the data after linear transformation
# a:  translation amount
# b:  scaling (stretching / shrinking of data)
# example:
# x <- rnorm(10) 
# y <- scale(x) 
# f = linear_convert_x_y(x=x, y=y)
# cbind(x,y, f(x=x)$y, f(y=y)$x)
#
linear_convert_x_y <- function(a, b, x=NULL, y=NULL) 
{
  if (!is.null(x) & !is.null(y)) {     # if data is given get calc a and b
    ab = determine_linear_trans_from_data(x, y) 
    a = ab["a"]
    b = ab["b"]
  }
  
  function(x=(y-a)/b, y = a + b*x) {    # closure for specific conversion
    list(x=x, y=y)
  }
}


# From the original and transformed data calculate infer coefficients a and b 
# to convert original data (x) into transformed data (y) using y = a + b*x.
#
# orig:  original data 
# trans:  transformed data used to construct biplot
#
get_matrix_of_linear_trans_coefs <- function(orig, trans) 
{
  p = ncol(orig)
  ab = matrix(NA, nrow=p, ncol=2) 
  colnames(ab) = c("a", "b")
  for (j in 1L:p)
    ab[j, 1:2] <- determine_linear_trans_from_data(orig[, j], trans[, j])
  ab
}


#### +-- Drawing functions ####

draw_star <- function(x, ...) 
{
  ctr <- colMeans(x, na.rm = TRUE)    # centroid
  segments(ctr[1], ctr[2], x[, 1], x[, 2], ...)
}


# x : biplot object
add_biplot_vectors <- function(x, dim=1:2, col=1, j=NULL, ...) 
{
  R = x$R
  p = nrow(R)
  Q = x$Q
  if (is.null(Q))             # rotation matrix for whole biplot
    Q = diag(1, p, p)
  R = R %*% Q                 # rotate biplot (see Gower et al. 2011, p. 20)
  
  
  if (is.null(j))
    j <- 1L:nrow(R)
  if (length(col) < p)
    col= rep(col, length.out=p)

  js <- 1L:nrow(R)    # all row indexes
  js <- js[j]         # indexes of vectors to draw
  for (j in js) {
    v <- R[j, dim]
    arrows(0, 0, v[1], v[2], col=col[j], ...)
  }
}


# x : biplot object
add_rowpoints <- function(x, dim=1:2, z=3, i=NULL, col=1, 
                          cex=1, alpha=1,
                          use.depth=FALSE, z.cex.to=c(.7, 1), 
                          use.alpha=FALSE,  z.alpha.to=c(.5, 1), 
                          ...) 
{
  L = x$L
  n = nrow(L)
  R = x$R
  p = nrow(R)
  Q = x$Q
  if (is.null(Q))             # rotation matrix for whole biplot
    Q = diag(1, p, p)
  L = L %*% Q                 # rotate biplot (see Gower et al. 2011, p. 20)
  
  if (is.null(i))
    i <- 1L:n
  if (length(col) < n)
    col <- rep(col, length.out=n)
  if (length(cex) < n)
    cex <- rep(cex, length.out=n)
  if (length(alpha) < n)
    alpha <- rep(alpha, length.out=n)
  
  if (use.depth) {
    zv <- as.vector(L[ , z])
    cex <- scales::rescale(zv, to=z.cex.to)
  }
  if (use.alpha) {
    zv <- as.vector(L[ , z])
    alpha <- scales::rescale(zv, to=z.alpha.to)
  }
  
  col <- set_alpha_color_value(col, alpha)  
  is <- 1L:n          # all row indexes
  is <- is[i]         # indexes of points to draw
  #for (i in 1:n)
  points(L[is, dim, drop=FALSE], col=col, cex=cex, ...)
}


# TODO: rad.outer argument
#
add_biplot_axes <- function(x, dim=1:2, j=NULL, col=1, shift=0, ...) 
{
  R = x$R
  Q = x$Q
  p = nrow(R)
  if (is.null(Q))     # rotation matrix for whole biplot
    Q = diag(1, p, p)
  R = R %*% Q         # rotate biplot (see Gower et al. 2011, p. 20)
  radius = x$mx
  if (is.null(j))
    j <- 1L:p
  if (length(col) < p)
    col= rep(col, length.out=p)
  if (length(shift) < p)
    shift= rep(shift, length.out=p)
  
  js <- 1L:p          # all row indexes
  js <- js[j]         # indexes of vectors to draw
  for (j in js)
    draw_axis_label_segments(R[j, dim], rads=c(0, radius), 
                             col=col[j], shift=shift[j], ...)
}


add_biplot_axes_extensions <- function(x, dim=1:2, j=NULL, col=1, shift=0, rad.extra.mm=5, ...) 
{
  R = x$R
  Q = x$Q
  p = nrow(R)
  if (is.null(Q))     # rotation matrix for whole biplot
    Q = diag(1, p, p)
  R = R %*% Q         # rotate biplot (see Gower et al. 2011, p. 20)
  radius = x$mx
  if (is.null(j))
    j <- 1L:p
  if (length(col) < p)
    col= rep(col, length.out=p)
  if (length(shift) < p)
    shift= rep(shift, length.out=p)
  
  js <- 1L:p          # all row indexes
  js <- js[j]         # indexes of vectors to draw
  rad.extra <- scale_vec_length_to_mm(c(1,0), mm=rad.extra.mm)[1]  # extra radius for axis extensions (in mm)
  for (j in js)
    draw_axis_label_segments(R[j, dim], rads=c(radius, radius + rad.extra), 
                             col=col[j], shift=shift[j], ...)
}


# col: colors for axes
add_axes_calibration <- function(x, mx, dim=1:2, j=NULL, shift=0, n=10, digits=3, col=1, 
                                 tick.ext=c(1,1), mm=1,
                                 ticks.col=col, srt=0, line=.5, ortho=FALSE, 
                                 zero=TRUE, ...) 
{
  R = x$R
  Q = x$Q
  p = nrow(R)
  if (is.null(Q))     # rotation matrix for whole biplot
    Q = diag(1, p, p)
  R = R %*% Q         # rotate biplot (see Gower et al. 2011, p. 20)
  
  AB= x$AB            # get linear trans coefficients for axis (calculate here or at beginning?)
  
  if (is.null(j))
    j <- 1L:p
  if (length(col) < p)
    col= rep(col, length.out=p)
  if (length(ticks.col) < p)
    ticks.col= rep(ticks.col, length.out=p)
  if (length(shift) < p)
    shift= rep(shift, length.out=p)
  if (length(zero) < p)
    zero= rep(zero, length.out=p)
  
  js <- 1L:p          # all row indexes
  js <- js[j]         # indexes of vectors to draw
  
  for (j in js) {
    v = R[j, dim]             # current biplot vector
    ab = AB[j, ]              # get linear trans coefficients for axis (from orig to transformed data)
    l = get_biplot_lambda_and_values(v, mx, ab=ab, n=n, zero=zero[j])   # for 
    lambdas = l$lambdas
    labels = round(l$origs, digits)
    # given the lambdas, i.e. multiples of biplot vectors, and according value labels
    add_calibration_to_biplot_axis_2(v, lambdas=lambdas, labels=labels, mm=mm, shift=shift[j], 
                                     tick.ext=tick.ext, ticks.col=ticks.col[j], srt=srt, line=line, 
                                     ortho=ortho, col=col[j], ...)
  }
}


# ... passed o to segments
# col color for each axis
add_point_axes_projections <- function(x, i=NULL, j=NULL, dim=1:2, shift=0, col=1, ...) 
{
  L = x$L
  R = x$R
  p = nrow(R)
  n = nrow(L)
  
  Q = x$Q
  if (is.null(Q))             # rotation matrix for whole biplot
    Q = diag(1, p, p)
  L = L %*% Q                 # rotate biplot (see Gower et al. 2011, p. 20)
  R = R %*% Q         
  
  if (is.null(i))             # index for which row points to project
    i <- 1L:n
  if (is.null(j))             # index for which biplot vectors to project onto
    j <- 1L:p
  if (length(shift) < p)      # shift for each axis
    shift= rep(shift, length.out=p)
  if (length(col) < p)
    col= rep(col, length.out=p)
  
  is <- 1L:n          # all row indexes
  is <- is[i]         # indexes of vectors to draw
  js <- 1L:p          # all row indexes
  js <- js[j]         # indexes of vectors to draw
  add_point_axis_projections(L[is, dim, drop=FALSE], R[js, dim, drop=FALSE], 
                             shift=shift, col=col, ...)
}


# getter for rownames
get_row_labels <- function(x) 
{
  l <- x$rowlabels
  l$labels
}

# getter for column names
get_column_labels <- function(x) 
{
  x$columnlabels
}


get_label_positions <- function(x, dim=1:2, i=NULL, inside=TRUE, n.helper=200, rad.extra.mm=10) 
{
  L = x$L
  R = x$R
  p = nrow(R)
  n = nrow(L)
  
  Q = x$Q
  if (is.null(Q))             # rotation matrix for whole biplot
    Q = diag(1, p, p)
  L = L %*% Q                 # rotate biplot (see Gower et al. 2011, p. 20)
  
  rad.outer = x$mx
  ori = x$ori
  
  if (is.null(i))             # index for which row points to project
    i <- 1L:n
  
  r1 = rad.outer
  r.extra = (mm_to_usr(rad.extra.mm) - mm_to_usr(0))[1]
  r2 = rad.outer + r.extra
  
  # calc coords for circle points around circumference
  rs = seq(0, 2*pi, length=n.helper)
  x1 = cos(rs) 
  y1 = sin(rs)
  x2 = cos(rs + rs[2]/2) 
  y2 = sin(rs + rs[2]/2)
  cc = cbind(x1 * r1 + ori[1], y1 * r1 + ori[2])        # circle coordinates to matrix
  cc2 = cbind(x2 * r2 + ori[1], y2 * r2 + ori[2])       # circle coordinates to matrix
  
  # add invisible points on the circle so the labels end up inside inside the circle        
  HP = rbind(cc, cc2)                       # helper points
  #points(HP, pch=16)
  helper.labels = rep("  ", 2*n.helper)     # helper ghost labels
  
  # row point labels
  point.labels = get_row_labels(x)
  if (length(point.labels) != n)                      # check, else pointLabels might go astray
    stop("Number of points (", n, ") does not match the number of labels (", length(point.labels), ")", call. = FALSE)
  
  # calc label positions
  if (!inside) {
    helper.labels=NULL
    HP = NULL
  }
  
  # select which points to draw
  is <- 1L:n          # all row indexes
  is <- is[i]         # indexes of points to draw

  labels = c(point.labels[is], helper.labels) 
  AP = rbind(L[is, dim, drop=FALSE], HP)        
  P <- dissertation::pointLabel(AP, labels=labels, doPlot=FALSE)
  xy.coords(P$x[is], P$y[is])
}


# just add some vertical offset to labels
#
get_label_positions_standard <- function(x, dim=1:2, i = NULL) 
{
  n = nrow(x$L)
  if (is.null(i))             # index for which row points to project
    i <- 1L:n
  L = x$L
  Q = x$Q

  L = L %*% Q                 # rotate biplot (see Gower et al. 2011, p. 20)
  y.off <- strheight("m")
  P = L - y.off
  P[i, dim]
}


# n.helper:   number of helper points on circle
# TODO: fix arguments passing
# build simpler function that is faster
# inside:   position labels so they do not trespass the circumference
#
add_rowpoint_labels <- function(x, dim=1:2, z=3, i=NULL, col=1, 
                                method=1, 
                                  cex=1, alpha=1,
                                  use.depth=FALSE, z.cex.to=c(.7, 1), 
                                  use.alpha=FALSE,  z.alpha.to=c(.5, 1),
                                  n.helper = 100, rad.extra.mm=10, inside=TRUE, ...) 
{
  L = x$L
  R = x$R
  p = nrow(R)
  n = nrow(L)
  
  Q = x$Q
  if (is.null(Q))             # rotation matrix for whole biplot
    Q = diag(1, p, p)
  L = L %*% Q                 # rotate biplot (see Gower et al. 2011, p. 20)

  rad.outer = x$mx
  ori = x$ori
  
  r1 = rad.outer
  r.extra = (mm_to_usr(rad.extra.mm) - mm_to_usr(0))[1]
  r2 = rad.outer + r.extra
  
  # other args
  if (is.null(i))
    i <- 1L:n
  if (length(col) < n)
    col <- rep(col, length.out=n)
  if (length(cex) < n)
    cex <- rep(cex, length.out=n)
  if (length(alpha) < n)
    alpha <- rep(alpha, length.out=n)
  
  if (use.depth) {
    zv <- as.vector(L[ , z])
    cex <- scales::rescale(zv, to=z.cex.to)
  }
  if (use.alpha) {
    zv <- as.vector(L[ , z])
    alpha <- scales::rescale(zv, to=z.alpha.to)
  }
  
  col <- set_alpha_color_value(col, alpha)  
  
  # select which points to draw
  is <- 1L:n          # all row indexes
  is <- is[i]         # indexes of points to draw
  
  if (method == 1)
    P = get_label_positions(x, dim=dim, i = is, inside = inside, n.helper = n.helper, 
                            rad.extra.mm = rad.extra.mm)
  if (method ==2)
    P = get_label_positions_standard(x, dim=dim, i=is)
  
  point.labels = get_row_labels(x)
  text(P, labels = point.labels[is], col=col[is], cex=cex[is], ...)
}
# Label placement
# maptools::pointLabel
# pointLabelTeachingDemos::spread.labels
# plotrix::spread.labels


# spread.labs( tmp.y[ tmp ], 1.2*strheight('A'),
#              maxiter=1000, min=0 )



# position labels by quadrant
pos_by_quadrant <- function(x, y=NULL, ori=c(0,0)) 
{
  xy = xy.coords(rbind(x), y)
  x = xy$x
  y = xy$y
  pos = rep(NA, length(x))
  p[x <= ori[1]]  <- 2   # left quadrants
  p[x > ori[1]] <- 4     # right quadrants
  p
}


# we can either label one or both axis poles (i.e. unipolar or bipolar variables)
# TODO : place offset in biplot object?
#
add_axes_pole_label <- function(x, dim=dim, j=NULL, j2=j, rad.extra.mm=5, rad.offset.mm=2, 
                                shift=0, col=1, col2=1, ...) 
{
  R = x$R
  p = nrow(R)
  
  Q = x$Q
  if (is.null(Q))             # rotation matrix for whole biplot
    Q = diag(1, p, p)
  R = R %*% Q                 # rotate biplot (see Gower et al. 2011, p. 20)
  
  rad.outer = x$mx
  rad.extra <- scale_vec_length_to_mm(c(1,0), mm=rad.extra.mm)[1]           # extra radius for axis extensions (in mm)
  rad.offset.lbls <- scale_vec_length_to_mm(c(1,0), mm=rad.offset.mm)[1]    # offset or labels (in mm) in direction  of axes
  rad.labels = rad.outer + rad.extra + rad.offset.lbls                      # complete radius for labels
  
  if (is.null(j))
    j <- 1L:p
  if (is.null(j2))
    j2 <- 1L:p
  if (length(col) < p)
    col= rep(col, length.out=p)
  if (length(col2) < p)
    col2= rep(col2, length.out=p)
  if (length(shift) < p)
    shift= rep(shift, length.out=p)
  
  columnlabels <- get_column_labels(x)
  pole1 = columnlabels$labels
  pole2 = columnlabels$labels2
  
  # unipolar
  js <- 1L:p          # all row indexes
  js <- js[j]         # indexes of vectors to draw
  for (j in js) {
    v = R[j, dim]                                   # get biplot axis vector
    v = extend_vec_to_length(v, len = rad.labels)   # extend to length of axis plus overhead
    n = get_2d_normal_vector(v, upwards=TRUE)       # get vector normal to axis 
    vt = v + shift[j] * n                           # add normal shift (shift is unit vector)
    pos = pos_by_quadrant(vt)
    text(t(vt), labels=pole1[j], col=col[j], pos=pos, offset=0, ...)   # draw axis labels
  }  
  
  # add on for bipolar variables (second pole)
  if (!is.null(pole2)) {
    js2 <- 1L:p          # all row indexes
    js2 <- js2[j2]         # indexes of vectors to draw
    for (j in js2) {
      v = R[j, dim]                                   # get biplot axis vector
      v = extend_vec_to_length(v, len = rad.labels)   # extend to length of axis plus overhead
      n = get_2d_normal_vector(v, upwards=TRUE)       # get vector normal to axis 
      vt = -v + shift[j] * n                           # add normal shift (shift is unit vector)
      pos = pos_by_quadrant(vt)
      text(t(vt), labels=pole2[j], col=col2[j], pos=pos, offset=0, ...)   # draw axis labels
    }    
  }
}
# calc all data for spreadlabels
# V = t(apply(R[, 1:2], 1, extend_vec_to_length, len = rad.labels))
# N = t(apply(V, 1, get_2d_normal_vector, upwards=TRUE))
# VT = V + diag(shift) %*% N
# pos = apply(VT, 1, pos_by_quadrant)

# maptools::pointLabel
# TeachingDemos::spread.labs
# plotrix::spread.labels


#### List objects ####

# pass labels as row and column names in matrix 
new_labels_obj <- function(labels, ...) 
{
  l <- list(labels=labels, ...)
  class(l) <- "labels.obj"
  l
}


#### setup ####

prep_biplot_object <- function(x, dim=1:2, ori=c(0,0)) 
{
  x$AB = get_matrix_of_linear_trans_coefs(x$X, x$Xtrans)   # linear transformation coefficients to convert X.orig into X.trans
  vns <- norm_vec(t(L[, dim]))    # find plot limits by most distant point from origin 
  x$mx = max(vns)                 # maximum length of biplot vector for given dimenions
  x$ori = ori
  x
}


set_up_plot <- function(x, dim=1:2, ext=1.2) 
{ 
  mxplus = x$mx * ext              # add some extra space
  lim = c(-mxplus, mxplus)        # usr limits of plotting region
  par(mar=rep(0,4))
  plot.new()
  plot.window(xlim=lim, ylim=lim, asp=1)
}


add_circumference <- function(x) 
{
  ori <- x$ori                   # origin of biplot
  rad.outer <- x$mx             # outer radius
  circle.coords <- plotrix::draw.circle(ori[1], ori[2],radius=rad.outer,
                                        nv=100, border=1, col=grey(.98), lty=1, lwd=1)  
  #add_gradient_circle(ori[1], ori[2],radius=rad.outer, col=grey(c(.98, .5, .2)), values=c(0, .7, 1))
}




##############

library(plotrix)
library(scales)

# see draw.circle, gradient_n_pal for arguments
# n: number of circles for gradient
# examples
#
# lim = c(-2,2)
# plot(NULL, xlim=lim, ylim=lim)
# add_gradient_circle()
#
# 
# set.seed(0)
# n = 20
# lim = c(-4,4)
# plot(NULL, xlim=lim, ylim=lim)
# x = runif(n, -2, 2) 
# y = runif(n, -2, 2) 
# z = runif(n, -2, 2) 
# add_gradient_circle(r=3, col=grey(c(.98, .5, .2)), values=c(0, .7, 1))
# cex = rescale(z, to=c(1, 1.3))
# alpha = rescale(z, to=c(.2, .7))
# col = set_alpha_color_value("blue", alpha)  
# points(x, y, cex=cex, pch=16, col=col)
# labels = paste("case", 1:n)
# cex = rescale(z, to=c(.7, 1))
# text(x, y, labels=labels, cex=cex, pch=16, col=col, pos=1)
#
add_gradient_circle <- function(x=0, y=0, radius=1, col=c("white", "black"), values=NULL, n = 50)
{  
  colfun <- scales::gradient_n_pal(col, values=values)
  is <- seq(1, 0, length=n)
  dc <- Vectorize(draw.circle)
  dummy <- dc(x, y, radius = is*radius, col=colfun(is), border=NA)
}



#### +----------- Generate biplot -------------- ####

# generate rotation matrix that will rotate biplot in p dimenions
# In order to just target the displayed dimenions, a 2D rotation matrix 
# is introduced in the according positions
#
# round(get_rotation_for_2d(4, c(1,3), alpha=pi/4), 2)
#
get_rotation_for_2d <- function(p, dim=1:2, alpha=0, clockwise=TRUE) 
{
  Q <- diag(1, p, p)
  Q[dim, dim] <-  rotation_matrix_2d(alpha=alpha, clockwise=clockwise)
  Q
}


biplot_obj_init <- function(X, Xtrans) 
{
  l <- list(X=X,
            Xtrans=Xtrans)
  class(l) <- "biplot.data.obj"
  l
}



biplot_obj_add_biplot_coords <- function (x, g = 0, h = 1 - g)
{
  # retrieve
  X <- x$X
  Xtrans <- x$Xtrans
  
  # biplot decomposition
  dec <- svd(Xtrans)
  V = dec$v
  U = dec$u
  r = length(dec$d) 
  D = diag(dec$d, nrow = r, ncol = r)
  DL = diag(dec$d^g, nrow = r, ncol = r)   # to circumentvent D^0 problem
  DR = diag(dec$d^h, nrow = r, ncol = r)
  L <- U %*% DL
  R <- V %*% DR
  
  # return (more code but cleaner)
  x$U = U
  x$V = V
  x$D = D
  x$L = L
  x$R = R
  x$g = g
  x$h = h
  x
}


biplot_obj_add_labels <- function(x, rlabs=rownames(x$X), 
                                  clabs=colnames(x$X), clabs2=NULL, type=2) 
{
  # defaults for row and column labels
  if (type == 1) 
  {
    if (is.null(rlabs))
      rlabs <- rep("", n)
    if (is.null(clabs))
      clabs <- rep("", p)
    if (is.null(clabs2))
      clabs2 <- rep("", p) 
  } 
  if (type == 2) 
  {
    if (is.null(rlabs))
      rlabs <- paste("case", 1:n)
    if (is.null(clabs))
      clabs <- paste("var", 1:p)
    if (is.null(clabs2))
      clabs2 <- paste("var", 1:p) 
  }
  
  # add to biplot object
  x$rowlabels <- new_labels_obj(labels=rlabs)
  x$columnlabels  <- new_labels_obj(labels=clabs,
                                    labels2=clabs2)
  x
}





