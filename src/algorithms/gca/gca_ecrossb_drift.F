! this "function" takes
! ... the field quantities: `bx0`, `by0`, `bz0`, `ex0`, `ey0`, `ez0` ...
! ... their norms `e0_SQR` and `b0_SQR` ...
! ... and returns the `vE` components and the `gammaE`

dummy_ = 1.0 / (b0_SQR + TINYFLD)

vE_x = (bz0 * ey0 - by0 * ez0) * dummy_
vE_y = (-bz0 * ex0 + bx0 * ez0) * dummy_
vE_z = (by0 * ex0 - bx0 * ey0) * dummy_

gammaE = 1.0 / sqrt(abs(1.0 - vE_x**2 - vE_y**2 - vE_z**2) + TINYFLD)
