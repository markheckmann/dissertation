#### Visualization function ####

#### +-- convex hulls ####

# xy    xy.coords object
#
convex_hull_standard <- function(xy, ...)
{
  hpts <- chull(xy$x, xy$y)
  xy.coords(xy$x[hpts], xy$y[hpts])
}


# "peeled" convex hull: Removing those points that create the hull 
#
# xy    xy.coords object
#
convex_hull_peeled <- function(xy, ...)
{
  hpts <- chull(xy$x, xy$y)               # create intitial hull
  i <- seq_along(xy$x)                    # index for al points
  i.pld <- i[-hpts]                       # index for non-hull points
  j <- chull(xy$x[i.pld], xy$y[i.pld])    # another hull points without initial hull points
  hpts.pld <- i.pld[j]                    # get indexes of new hull points
  xy.coords(xy$x[hpts.pld], xy$y[hpts.pld])
}


# Splining a polygon
#
#   The rows of 'xy' give coordinates of the boundary vertices, in order.
#   'vertices' is the number of spline vertices to create.
#              (Not all are used: some are clipped from the ends.)
#   'k' is the number of points to wrap around the ends to obtain
#       a smooth periodic spline.
#   Returns an array of points. 
# source: 
# http://gis.stackexchange.com/questions/24827/how-to-smooth-the-polygons-in-a-contour-map
#
spline_poly <- function(x, y=NULL, vertices=100, k=3, ...) 
{
  xyc <- xy.coords(x, y, recycle = TRUE)
  xy <- cbind(xyc$x, xyc$y)    # xy is an n by 2 matrix with n >= k.
  
  # Wrap k vertices around each end.
  n <- dim(xy)[1]
  if (k >= 1) {
    data <- rbind(xy[(n-k+1):n,], xy, xy[1:k, ])
  } else {
    data <- xy
  }
  
  # Spline the x and y coordinates.
  data.spline <- spline(1:(n+2*k), data[,1], n=vertices, ...)
  x <- data.spline$x
  x1 <- data.spline$y
  x2 <- spline(1:(n+2*k), data[,2], n=vertices, ...)$y
  
  i <- k < x & x <= n+k     # Retain only the middle part.
  xy.coords(x1[i], x2[i])   # return as xy.coords object
}


# another approach to generating s spline polygon
spline2 <- function(xy, n=100) 
{
  xy.coords(x= spline(xy$x, n=n)$y, 
            y= spline(xy$y, n=n)$y)
}

# Greenacre (2006, p.179). "Tying up the loose ends ..." suggest to create 
# convex hulls with alpha-% of the points by removing the most outlying 
# points. Function removes proportion of points that is most distant
# from the centroid.
#
# xy :    xy.coords object
#
remove_distant_points <- function(xy, alpha=.05) 
{
  x = xy$x
  y = xy$y
  d <- sqrt((x - mean(x))^2 + (y - mean(y))^2)  # Euclidean distance from centroid
  i <- order(d, decreasing=TRUE)
  n.rm <- ceiling(length(d) * alpha)            # get number of points to remove
  if (n.rm > 0)                                 # more than 0 points to be removed?
    i <- i[-(1:n.rm)]
  xy.coords(x[i], y[i])  
}


#' Draw different types of convex hulls around a set of points
#' 
#' @param x,y  The x and y coordinates of a set of points. 
#'  Alternatively, a single argument x can be provided.
#' @param alpha Proportion of points most distant from centroid to 
#'  be removed before creating the (peeled) convex hull.
#' @param peel Remove hull points? (useful to remove outliers, 
#'  default \code{FALSE})
#' @param smooth Integer (1 to 3) indicating the method for creating splines 
#'  around the hull.
#' @param shape Shape parameter used for \code{smooth=1} 
#'  (see \code{xspline} for details).
#' @param draw Draw the hull polygon (default \code{TRUE}).
#' @param ... additional parameters passed to \code{polygon} to draw the hull.
#' @return A list as given by \code{\link{xy.coords}} containing the 
#'  hull points is returned invisibly. 
#' @export
#' @example examples/example-convex-hull.R
#'   
convex_hull <- function(x, y=NULL, alpha=0, peel=FALSE, 
                        smooth=0, shape=1, draw=TRUE, ...)
{
  xy <- xy.coords(x, y, recycle=TRUE)
  if (alpha > 0)            # proportion of most distanc points to be removed
    xy <- remove_distant_points(xy, alpha)
  if (peel)                 # peel off convex hull points
    ch <- convex_hull_peeled(xy)
  else 
    ch <- convex_hull_standard(xy)
  # spline interpolations of convex hull points
  if (smooth == 1) 
    ch <- xspline(ch, shape=shape, open = FALSE, draw = FALSE)
  if (smooth == 2) 
    ch <- spline_poly(ch)
  if (smooth == 3)
    ch <- spline2(ch)
  polygon(ch, ...)
  invisible(ch)
}  


# identify outliers and remove from hull?

# library(mvoutlier)
# ?pcout
# data(bsstop)
# x=bsstop[,5:14]
# # identify multivariate outliers
# x.out = pcout(x, makeplot=FALSE)
# which(x.out$wfinal01 == 0)




#### +-- density plots ####

#' Calculate a critical denisty value from the 2D density distribution.
#' The critical value can be used.
#' @param den Object returned by \code{kde2d}.
#' @param n A bit uncelar needs a check.
#' @param prob Probability used to determine critical value.
#' @export
#' @keywords internal
#' @section TODO:
#' double check meaning of n
#' http://grokbase.com/t/r/r-help/1233hxv1a9/r-contour-for-plotting-confidence-interval-on-scatter-plot-of-bivariate-normal-distribution
#' n might have to equal nobs
#' http://quantitative-ecology.blogspot.de/2008/03/plotting-contours.html
#'
get_crit_density_value <- function(den, n=100, prob=.05)
{
  den.z <- den$z
  z <- vector("numeric")
  for (i in 1:n) {
    suppressWarnings({     
      z.x <- max(which(den$x < x[i]));
      z.y <- max(which(den$y < y[i]));
      z[i] <- den$z[z.x, z.y]
    })
  }
  # store class/level borders of confidence interval in variables
  crit.dens.value <- quantile(z, probs=prob, na.rm = TRUE) # +-1sd
  crit.dens.value
}


#' Draw density contour plot using a cut-off value.
#' 
#' @param den Object returned by \code{kde2d}.
#' @param n A bit uncelar needs a check.
#' @param prob Probability used to determine critical cut-off value for density.
#' @param ... Passes to image or contour.
#' @export
#' @rdname crit-density
#' @example examples/example-crit-density.R
#' 
crit_contour <- function(den, n=100, prob=.05, ...) 
{
  crit <- get_crit_density_value(den, n=n, prob=prob)
  contour(den, levels=crit, add = TRUE, ...)
}

#' @rdname crit-density
#' 
crit_dens_image <- function(den, prob=.05, ...) 
{
  crit <- get_crit_density_value(den, prob=prob)
  den$z[den$z < crit] <- NA
  image(den, add = T, ...)
}



