## ----setup, echo=FALSE, warning=FALSE, message=FALSE---------------------
library(mat2tex)

# house shapes
draw_configuration <- function(m, ...) {
  x <- rbind(m, head(m, 1))
  polygon(x, ...)
}


## ----eval=FALSE----------------------------------------------------------
#  library(devtools)
#  install_github("procrustes", "markheckmann")

## ----eval=FALSE----------------------------------------------------------
#  library(procrustes)

## ----echo=FALSE, results='asis'------------------------------------------
opt <- mat2tex_options(digits=0)  
A <- matrix(c(0,1,-1,-1,1,-0), by=TRUE, 3)
xx("\\B{A} = ", A, e=1)

