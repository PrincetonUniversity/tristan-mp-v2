! this "function" takes
! ... the field quantities: `bx0`, `by0`, `bz0`, `ex0`, `ey0`, `ez0` ...
! ... and the velocities: `u0`, `v0`, `w0` ...
! ... and returns the updated velocities `u0`, `v0`, `w0`
dummy_ = 0.5 * q_over_m * B_norm
ex0 = ex0 * dummy_; ey0 = ey0 * dummy_; ez0 = ez0 * dummy_
dummy_ = dummy_ * CCINV
bx0 = bx0 * dummy_; by0 = by0 * dummy_; bz0 = bz0 * dummy_

! half acceleration:
u0 = CC * u0 + ex0
v0 = CC * v0 + ey0
w0 = CC * w0 + ez0

! first half magnetic rotation:
g_temp = CC / sqrt(CC**2 + u0**2 + v0**2 + w0**2)

bx0 = g_temp * bx0
by0 = g_temp * by0
bz0 = g_temp * bz0
dummy_ = 2.0 / (1.0 + bx0 * bx0 + by0 * by0 + bz0 * bz0)
u1 = (u0 + v0 * bz0 - w0 * by0) * dummy_
v1 = (v0 + w0 * bx0 - u0 * bz0) * dummy_
w1 = (w0 + u0 * by0 - v0 * bx0) * dummy_
! second half magnetic rotation + half acceleration:

u0 = (u0 + v1 * bz0 - w1 * by0 + ex0) * CCINV
v0 = (v0 + w1 * bx0 - u1 * bz0 + ey0) * CCINV
w0 = (w0 + u1 * by0 - v1 * bx0 + ez0) * CCINV
