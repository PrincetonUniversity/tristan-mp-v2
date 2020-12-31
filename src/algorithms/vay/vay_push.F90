! this "function" takes
! ... the field quantities: `bx0`, `by0`, `bz0`, `ex0`, `ey0`, `ez0` ...
! ... and the velocities: `u0`, `v0`, `w0` ...
! ... and returns the updated velocities `u0`, `v0`, `w0`
dummy_ = 0.5 * q_over_m * B_norm
ex0 = ex0 * dummy_; ey0 = ey0 * dummy_; ez0 = ez0 * dummy_
dummy_ = dummy_ * CCINV
bx0 = bx0 * dummy_; by0 = by0 * dummy_; bz0 = bz0 * dummy_

dummy_ = 1.0 / sqrt(1.0 + u0**2 + v0**2 + w0**2)

u1 = CC * u0 + 2.0 * ex0 + (CC * v0 * dummy_) * bz0 - (CC * w0 * dummy_) * by0
v1 = CC * v0 + 2.0 * ey0 + (CC * w0 * dummy_) * bx0 - (CC * u0 * dummy_) * bz0
w1 = CC * w0 + 2.0 * ez0 + (CC * u0 * dummy_) * by0 - (CC * v0 * dummy_) * bx0

dummy_ = CCINV * ( u1 * bx0 + v1 * by0 + w1 * bz0)
dummy2_ = CCINV * CCINV * (CC**2 + u1**2 + v1**2 + w1**2) - (bx0**2 + by0**2 + bz0**2)

dummy_ = 1.0 / sqrt(0.5 * (dummy2_ + sqrt(dummy2_**2 + 4.0 * (bx0**2 + by0**2 + bz0**2 + dummy_**2))))
dummy2_ = 1.0 / (1.0 + (bx0 * dummy_)**2 + (by0 * dummy_)**2 + (bz0 * dummy_)**2)

u0 = dummy2_ * (u1 + (u1 * (bx0 * dummy_) + v1 * (by0 * dummy_) + w1 * (bz0 * dummy_)) * (bx0 * dummy_) + v1 * (bz0 * dummy_) - w1 * (by0 * dummy_))
v0 = dummy2_ * (v1 + (u1 * (bx0 * dummy_) + v1 * (by0 * dummy_) + w1 * (bz0 * dummy_)) * (by0 * dummy_) + w1 * (bx0 * dummy_) - u1 * (bz0 * dummy_))
w0 = dummy2_ * (w1 + (u1 * (bx0 * dummy_) + v1 * (by0 * dummy_) + w1 * (bz0 * dummy_)) * (bz0 * dummy_) + u1 * (by0 * dummy_) - v1 * (bx0 * dummy_))

u0 = u0 * CCINV
v0 = v0 * CCINV
w0 = w0 * CCINV
