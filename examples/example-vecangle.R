v <- c(1,0)                                  # point vector
as <- seq(0, 2*pi, length=10)                # radians to rotate
# mode 1: between 0 and 180 degrees, minimal angle / radian
for (a in as)         
  print(vec_angle(v %*% rotation_matrix_2d(a, F), v, mode=1))
# mode 2: between 0 and 360 degrees, anti-clockwise 
for (a in as)
  print(vec_angle(v %*% rotation_matrix_2d(a, F), v, mode=2))
