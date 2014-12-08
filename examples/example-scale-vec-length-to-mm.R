par(mfrow=c(1,2))

plot(NULL, xlim=c(-1,1), ylim=c(-1,1), asp=2,
     main="vectors plot length\n varies by direction")

rs <- seq(0, 2*pi, length=41)
v = c(1,0)
for (r in rs) {
  R = rotation_matrix_2d(r, clockwise = F)  
  n = v %*% R
  segments(0,0, n[1], n[2]) 
}

plot(NULL, xlim=c(-1,1), ylim=c(-1,1), asp=2,
     main="vectors have 10mm plot\nlength in all directions")
for (r in rs) {
  R = rotation_matrix_2d(r, clockwise = F)  
  n = v %*% R
  b = scale_vec_length_to_mm(n, mm=10)
  segments(0,0, b[1], b[2]) 
}
