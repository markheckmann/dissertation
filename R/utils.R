# utility functions

xy_coords <- function(x, y = NULL, xlab = NULL, ylab = NULL, 
                      log = NULL, recycle = FALSE) 
{
  xy <- xy.coords(x, y, xlab, ylab, log, recycle)
  class(xy) <- "xy.coords"
  xy
}


"[.xy.coords" <- function(x, i) 
{
  x$x <- x$x[i]
  x$y <- x$y[i]
  x
}

as.data.frame.xy.coords <- function(x) 
{
  data.frame(x=x$x, y=x$y)
}

# xy = xy_coords(1:5, 5:1)
# xy[1:3]
# as.data.frame(xy)


#' Set alpha value to color
#' 
#' Given a color in any standard format, the color is 
#' converted to hex and an alpha value is added.
#' 
#' @param col A color in a standard format 
#'  (e.g. \code{"#000000"}, \code{"black"}, \code{1}).
#' @param alpha Alpha value to add to color 
#' (default \code{=1}, i.e. opaque).
#' @return hex color value.
#' @export
#' @examples
#' 
#'  set_alpha_color_value("black", .5)
#'  set_alpha_color_value(1:3, .5)
#'  set_alpha_color_value("#000000", .5)
#'  
set_alpha_color_value <- function(col, alpha=1) 
{
  k <- col2rgb(col) / 255
  rgb(red = k["red", ], 
      green = k["green", ],
      blue = k["blue", ],
      alpha = alpha, maxColorValue=1)  
}


#' Euclidean vector norm
#' @param x A vector or a matrix. If a matrix, each column is 
#'    treated as a vector and the Euclidean norm is calculated for each column.
#' @export
#' 
norm_vec <- function(x) 
{
  norm_vec_ <- function(x) sqrt(sum(x^2))
  if (is.vector(x))
    x <- cbind(x)
  apply(x, 2, norm_vec_)
}


#' Convert usr coordinates to mm
#' @param v A 2D vector 
#' @export
#' @keywords      internal
#' 
usr_to_mm <- function(v)
{
  n = c(grconvertX(v[1], from = "user", to = "inches"),
        grconvertY(v[2], from = "user", to = "inches"))
  n * 25.4        # convert inch to mm
}


#' Convert mm to usr coordinates
#' @param v A 2D vector 
#' @export
#' @keywords      internal
#' 
mm_to_usr <- function(v)
{
  n = v / 25.4    # convert mm to inch
  c(grconvertX(n[1], from = "inches", to = "user"),
    grconvertY(n[2], from = "inches", to = "user"))
}


#' Scale length of vector to given mm
#' 
#' This function will scale a vector in usr coordinates to mm
#' so the plotted length will be as desired. This will be particularly
#' useful in situtations where the plot aspect ratio is unequal to 1.
#' 
#' @param v   A 2D vector 
#' @param mm  Numeric. Length to standardize vector to in mm.
#' @export
#' @example examples/example-scale-vec-length-to-mm.R
#' 
scale_vec_length_to_mm <- function(v, mm=10)
{
  origin = usr_to_mm(c(0,0))
  n = usr_to_mm(v) - origin
  nn = n / sum(n^2)^.5
  mm * mm_to_usr(nn + origin) 
}


#' Get plot aspect ratio and 1mm in usr coordinates
#' 
#' @return A list with the following elements:
#'    \item{\code{asp}}{Plot aspect ratio}
#'    \item{\code{usr.1mm}}{width and height of 1 mm in usr coordinates}
#' @export
#' @keywords      internal
#' @examples
#' 
#' plot(1:10, asp=2)
#' get_plot_asp_ratio()
#' 
get_plot_asp_ratio <- function() 
{
  w <- par("pin")[1] / diff(par("usr")[1:2]) * 25.4   # mm per x usr-unit
  h <- par("pin")[2] / diff(par("usr")[3:4]) * 25.4   # mm per y usr-unit
  w.mm = 1/w                    # 1mm width in usr units 
  h.mm = 1/h                    # 1mm height in usr units 
  list(asp=h/w,                 # aspect ratio
       usr.1mm=c(w.mm, h.mm) )  # width and height of 1 mm in usr units
}



#' Calculate the trace of a matrix.
#' @param x       A matrix.
#' @return        The trace of matrix \code{x}
#' @export
#' @keywords      internal
#'
tr <- function(x){
  if (!is.matrix(x))
    stop("x must be a matrix")
  sum(diag(x)) 
} 


#' Conversion from degrees to radians
#' 
#' @param x     Degrees.
#' @export
#' @examples
#' 
#' degree_to_rad(seq(0, 360, length=45))
#' 
degree_to_rad <- function(x){
  x *  pi/180
}


#' Conversion from radians to degrees
#' 
#' @param x     Radians.
#' @export
#' @examples
#' 
#' rad_to_degree(seq(0, 2*pi, length=pi/2))
#'
rad_to_degree <- function(x){
  x * 180/pi  
}



polar_to_cartesian <- function(rad, r) {
  r * c(cos(rad), sin(rad))
}


calc_rss <- function(A, B) {
  sum(diag(t(A - B) %*%  (A - B)))
}

#' Generate a 2D rotation matrix
#' @param alpha Angle in radians, clockwise
#' @param clockwise Rotate clockwise?
#' @return A 2D rotation matrix
#' @export
#' @examples
#' v <- c(0,1)                                  # point
#' as <- seq(0, pi/2, length=10)                # radians to rotate
#' plot(NULL, xlim=c(-1,1), ylim=c(-1,1), asp=1)
#' for (a in as)
#'   points(v %*% rotation_matrix_2d(a))
#' for (a in as)
#'   points(v %*% rotation_matrix_2d(a, FALSE), col="red")
#'   
rotation_matrix_2d <- function(alpha=0, clockwise=TRUE) 
{
  R = matrix(c(cos(alpha), sin(alpha), -sin(alpha), cos(alpha)), 2)
  if (! clockwise)
    R = R %*% diag(c(1,-1))
  R
}




#' Generate a all rotation planes for 3D rotation matrix
#' @param a Angle in radians
#' @return A list with rotation upon all of the three axes.
#' @export
#' 
#' 
rotation_matrix_3d <- function(a) 
{
  R.x <- matrix( c(1, 0, 0, 0, cos(a), sin(a), 0, -sin(a), cos(a) ), 3)
  R.y <- matrix( c(cos(a), 0, -sin(a), 0, 1, 0, sin(a), 0, cos(a) ), 3)
  R.z <- matrix( c(cos(a), sin(a), 0, -sin(a), cos(a), 0, 0, 0, 1 ), 3)
  list(R.x=R.x, R.y=R.y, R.z=R.z)
}


interpolate_rotation_matrix <- function(R, x=.5) 
{
  e <- eigen(R)
  D <- diag(e$values)
  U <- e$vectors
  Rx <- Re(U %*% D^x %*% Conj(t(U)))  
  Rx
}



# TODO Martin
ArrowPlus <- function(x,y, mode="pp",...)
{  
  if (mode == "pv") {   
    #Punkt mit Vektor 
    Arrows(x[1],x[2],x[1]+y[1],x[2]+y[2], ...)
  } else if (mode == "pp") { 
    #Punkt zu Punkt
    Arrows(x[1],x[2],y[1],y[2], ...)
  } 
}



#' Angle / radian between two vectors.
#' 
#' Function comes in two flavors. Either it calculates the smallest angle 
#' [0, 180) or the anti-clockwise angle with respect to the second 
#' vector [0, 360).
#' 
#' @param a,b  Vector
#' @param rad  Angular measure is radian?  (default \code{TRUE}).
#' @param mode 1= smallest angle in [0, 180); 2=anti-clockwise angle with respect 
#' #' vector b in [0, 360).
#' @return Angle
#' @export
#' @example examples/example-vecangle.R
#' 
vec_angle <- function(a, b, rad=TRUE, mode=1) 
{
  # minimal radian between vectors
  if (mode == 1) {
    r <- acos(sum(a * b) /     
                ( sum(a^2)^.5 * sum(b^2)^.5) )    
  } 
  # anti-clockwise radian with respect to vector b
  if (mode == 2) {
    #r = atan2(a[1]*b[2] - b[1]*a[2], a[1]*b[1] + a[2]*b[2])
    r = atan2(b[1]*a[2] - a[1]*b[2], b[1]*a[1] + b[2]*a[2])
    
    #     r = (-180/pi * ang) %% 360     # restrict to 0 to 260 range     
    #     r = r * pi /180                # convert to radians
    r = r %% (2*pi)   # restrict to 2pi range     
  }

  if (!rad)   # return degress if requested
    r <- r * 180 / pi  # convert radians to angle
  r
}



#' Adjust text srt (rotation) to get readable directions
#' 
#' @details
#' srt values are mapped as follows:
#' 0 to 90 = keep, 
#' 271 to 260 = keep, 
#' 91 to 270 = add 180 degrees.
#' 
#' @param a  Vector of angles (anti-clockwise).
#' @param ortho Rotate so text is orthogonal to the axis?
#' @export
#' @example examples/example-harmonize-text.R
#'    
harmonize_text_srt <- function(a, ortho=FALSE)
{
  a.new <- ifelse(a > 90 & a <=270, a + 180, a)
  if (ortho)
    a.new <- ifelse( (a>90 & a<=180) | 
                    (a>=270 & a<360), a.new + 90, a.new -90)
  a.new
}


#' Generate offset for axis labels in usr coords.
#' 
#' @param v Vector of axis direction.
#' @param line Offset in line heights
#' @param below Offset in direction below thes axis?
#' @details
#' Generates an offset to place labels below or above an axis.
#' The offest is given in multiple of the lineheight.
#' @return A vector containing width and height for single character.
#' @examples
#'  axes_label_offset(c(1,1))
#'  
axes_label_offset <- function(v, line=1, below=TRUE) 
{
  v <- v / sum(v^2)^.5            # to unit vec
  R <- rotation_matrix_2d(pi/2)   # 90 degree rotation
  n <- rbind(v %*% R)             # get normal vector
  if (n[2] > 0)                   # reflect if not pointing upwards
    n <- n * -1
  if (!below)                     # reflect normal vector if label 
    n <- n * -1                   # should go above axis
  ch <- par()$cxy                 # char width and height in usr coords 
  n * ch * line
}


#' Get normal vector for given 2D vector.
#' 
#' @param v A 2D vector.
#' @param upwards Vector direction is selected so it points in 
#'  positive y-direction by default.
#' @return A 2D vector of length one.
#' @export
#' @examples
#' 
#'  v <- c(1,1)
#'  u <- get_2d_normal_vector(v)
#'  d <- get_2d_normal_vector(v, upwards=FALSE)
#'  plot(NULL, xlim=c(-1,1), ylim=c(-1,1), asp=1)
#'  arrows(0, 0, v[1], v[2])
#'  arrows(0, 0, u[1], u[2], col=2)
#'  arrows(0, 0, d[1], d[2], col=3)
#'  
get_2d_normal_vector <- function(v, upwards=TRUE)
{
  v <- as.vector(v)               # make sure its a column vector, if e.g. a 1x2 matrix is passed
  v <- v / sum(v^2)^.5            # to unit vec
  R <- rotation_matrix_2d(pi/2)   # 90 degree rotation
  n <- rbind(v %*% R)             # get normal vector
  if (n[2] < 0)                   # reflect if not pointing upwards
    n <- n * -1
  if (!upwards)                     # reflect normal vector if label 
    n <- n * -1                   # should go above axis  
  as.vector(n)
}



#' Add labels to an axis at given positions
#' 
#' @details
#' Places value labels at given points and in direction given by
#' vector. The labels are rotated automatically to be readable.
#' 
#' @param v Direction vector of axis.
#' @param U Points on axis in user coordinates
#' @param labels Labels to place on axis.
#' @param i Indexes of which labels to plot.
#' @param ortho Labels orthogonal to axis? (default \code{FALSE}).
#' @param below Should the label appear below the axis?
#' @param line Number of lines label is placed below axis.
#' @param srt A numerical value specifying (in degrees) how the text 
#'  should be rotated.
#' @param ... Passed on to \code{\link{text}} to plot labels.
#' @section TODO Respect cex argument, correct char width (SO question)
#' @export
#' @example examples/example-add-axis-labels.R
#' 
add_axes_labels <- function(v, U, labels, i=1:nrow(U), 
                            ortho=FALSE, below=TRUE, 
                            line=1, srt=NULL, ...)
{
  if (is.null(srt)) {
    srt <- vec_angle(v, c(1,0), rad=FALSE, mode=2)        # anti-clock angle between axis vector and x-axis 
    srt <- harmonize_text_srt(srt, ortho=ortho)           # make text direction readable
  } 
  d <- axes_label_offset(v, line = line, below=below)   # get offset to place text under or over axis
  Uoff <- sweep(U, 2, d, "+")                           # add x and y offest to all axis points
  text(Uoff[i, ,drop=F], labels=labels[i], srt=srt, ...)
}


#vec_angle(c(1,1), c(1,0), rad=T, mode=2)
