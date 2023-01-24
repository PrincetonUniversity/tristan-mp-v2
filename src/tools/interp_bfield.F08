! this "function" takes `dx`, `dy`, `dz` and `lind` ...
! ... and "returns" `bx0`, `by0`, `bz0`

! linear indices used instead of 3d
#ifdef oneD
  !b_x
  c0 = bx(lind, 0, 0)
  c1 = bx(lind + 1, 0, 0)
  bx0 = c0 * (1 - dx) + c1 * dx

  !b_y
  c0 = 0.5 * (by(lind - 1, 0, 0) + by(lind, 0, 0))
  c1 = 0.5 * (by(lind, 0, 0) + by(lind + 1, 0, 0))
  by0 = c0 * (1 - dx) + c1 * dx

  !b_z
  c0 = 0.5 * (bz(lind - 1, 0, 0) + bz(lind, 0, 0))
  c1 = 0.5 * (bz(lind, 0, 0) + bz(lind + 1, 0, 0))
  bz0 = c0 * (1 - dx) + c1 * dx
#elif defined(twoD)
  !b_x
  c000 = 0.5 * (bx(lind, -NGHOST, 0) + bx(lind - iy, -NGHOST, 0))
  c100 = 0.5 * (bx(lind + 1, -NGHOST, 0) + bx(lind + 1 - iy, -NGHOST, 0))
  c010 = 0.5 * (bx(lind, -NGHOST, 0) + bx(lind + iy, -NGHOST, 0))
  c110 = 0.5 * (bx(lind + 1, -NGHOST, 0) + bx(lind + iy + 1, -NGHOST, 0))
  c00 = c000 * (1 - dx) + c100 * dx
  c10 = c010 * (1 - dx) + c110 * dx
  bx0 = c00 * (1 - dy) + c10 * dy

  !b_y
  c000 = 0.5 * (by(lind - 1, -NGHOST, 0) + by(lind, -NGHOST, 0))
  c100 = 0.5 * (by(lind, -NGHOST, 0) + by(lind + 1, -NGHOST, 0))
  c010 = 0.5 * (by(lind + iy - 1, -NGHOST, 0) + by(lind + iy, -NGHOST, 0))
  c110 = 0.5 * (by(lind + iy, -NGHOST, 0) + by(lind + iy + 1, -NGHOST, 0))
  c00 = c000 * (1 - dx) + c100 * dx
  c10 = c010 * (1 - dx) + c110 * dx
  by0 = c00 * (1 - dy) + c10 * dy

  !b_z
  c000 = 0.25 * (bz(lind - 1 - iy, -NGHOST, 0) + bz(lind - 1, -NGHOST, 0) + &
                 bz(lind - iy, -NGHOST, 0) + bz(lind, -NGHOST, 0))
  c100 = 0.25 * (bz(lind - iy, -NGHOST, 0) + bz(lind, -NGHOST, 0) + &
                 bz(lind + 1 - iy, -NGHOST, 0) + bz(lind + 1, -NGHOST, 0))
  c010 = 0.25 * (bz(lind - 1, -NGHOST, 0) + bz(lind + iy - 1, -NGHOST, 0) + &
                 bz(lind, -NGHOST, 0) + bz(lind + iy, -NGHOST, 0))
  c110 = 0.25 * (bz(lind, -NGHOST, 0) + bz(lind + iy, -NGHOST, 0) + &
                 bz(lind + 1, -NGHOST, 0) + bz(lind + iy + 1, -NGHOST, 0))
  c00 = c000 * (1 - dx) + c100 * dx
  c10 = c010 * (1 - dx) + c110 * dx
  bz0 = c00 * (1 - dy) + c10 * dy
#elif defined(threeD)
  !b_x
  c000 = 0.25 * (bx(lind, -NGHOST, -NGHOST) + bx(lind - iy, -NGHOST, -NGHOST) + &
                 bx(lind - iz, -NGHOST, -NGHOST) + bx(lind - iy - iz, -NGHOST, -NGHOST))
  c100 = 0.25 * (bx(lind + 1, -NGHOST, -NGHOST) + bx(lind + 1 - iy, -NGHOST, -NGHOST) + &
                 bx(lind + 1 - iz, -NGHOST, -NGHOST) + bx(lind + 1 - iy - iz, -NGHOST, -NGHOST))
  c001 = 0.25 * (bx(lind, -NGHOST, -NGHOST) + bx(lind + iz, -NGHOST, -NGHOST) + &
                 bx(lind - iy, -NGHOST, -NGHOST) + bx(lind - iy + iz, -NGHOST, -NGHOST))
  c101 = 0.25 * (bx(lind + 1, -NGHOST, -NGHOST) + bx(lind + 1 + iz, -NGHOST, -NGHOST) + &
                 bx(lind + 1 - iy, -NGHOST, -NGHOST) + bx(lind + 1 - iy + iz, -NGHOST, -NGHOST))
  c010 = 0.25 * (bx(lind, -NGHOST, -NGHOST) + bx(lind + iy, -NGHOST, -NGHOST) + &
                 bx(lind - iz, -NGHOST, -NGHOST) + bx(lind + iy - iz, -NGHOST, -NGHOST))
  c110 = 0.25 * (bx(lind + 1, -NGHOST, -NGHOST) + bx(lind + 1 - iz, -NGHOST, -NGHOST) + &
                 bx(lind + 1 + iy - iz, -NGHOST, -NGHOST) + bx(lind + 1 + iy, -NGHOST, -NGHOST))
  c011 = 0.25 * (bx(lind, -NGHOST, -NGHOST) + bx(lind + iy, -NGHOST, -NGHOST) + &
                 bx(lind + iy + iz, -NGHOST, -NGHOST) + bx(lind + iz, -NGHOST, -NGHOST))
  c111 = 0.25 * (bx(lind + 1, -NGHOST, -NGHOST) + bx(lind + 1 + iy, -NGHOST, -NGHOST) + &
                 bx(lind + 1 + iy + iz, -NGHOST, -NGHOST) + bx(lind + 1 + iz, -NGHOST, -NGHOST))
  c00 = c000 * (1.0 - dx) + c100 * dx
  c01 = c001 * (1.0 - dx) + c101 * dx
  c10 = c010 * (1.0 - dx) + c110 * dx
  c11 = c011 * (1.0 - dx) + c111 * dx
  c0 = c00 * (1.0 - dy) + c10 * dy
  c1 = c01 * (1.0 - dy) + c11 * dy
  bx0 = c0 * (1.0 - dz) + c1 * dz

  !b_y
  c000 = 0.25 * (by(lind - 1 - iz, -NGHOST, -NGHOST) + by(lind - 1, -NGHOST, -NGHOST) + &
                 by(lind - iz, -NGHOST, -NGHOST) + by(lind, -NGHOST, -NGHOST))
  c100 = 0.25 * (by(lind - iz, -NGHOST, -NGHOST) + by(lind, -NGHOST, -NGHOST) + &
                 by(lind + 1 - iz, -NGHOST, -NGHOST) + by(lind + 1, -NGHOST, -NGHOST))
  c001 = 0.25 * (by(lind - 1, -NGHOST, -NGHOST) + by(lind - 1 + iz, -NGHOST, -NGHOST) + &
                 by(lind, -NGHOST, -NGHOST) + by(lind + iz, -NGHOST, -NGHOST))
  c101 = 0.25 * (by(lind, -NGHOST, -NGHOST) + by(lind + iz, -NGHOST, -NGHOST) + &
                 by(lind + 1, -NGHOST, -NGHOST) + by(lind + 1 + iz, -NGHOST, -NGHOST))
  c010 = 0.25 * (by(lind - 1 + iy - iz, -NGHOST, -NGHOST) + by(lind - 1 + iy, -NGHOST, -NGHOST) + &
                 by(lind + iy - iz, -NGHOST, -NGHOST) + by(lind + iy, -NGHOST, -NGHOST))
  c110 = 0.25 * (by(lind + iy - iz, -NGHOST, -NGHOST) + by(lind + iy, -NGHOST, -NGHOST) + &
                 by(lind + 1 + iy - iz, -NGHOST, -NGHOST) + by(lind + 1 + iy, -NGHOST, -NGHOST))
  c011 = 0.25 * (by(lind - 1 + iy, -NGHOST, -NGHOST) + by(lind - 1 + iy + iz, -NGHOST, -NGHOST) + &
                 by(lind + iy, -NGHOST, -NGHOST) + by(lind + iy + iz, -NGHOST, -NGHOST))
  c111 = 0.25 * (by(lind + iy, -NGHOST, -NGHOST) + by(lind + iy + iz, -NGHOST, -NGHOST) + &
                 by(lind + 1 + iy, -NGHOST, -NGHOST) + by(lind + 1 + iy + iz, -NGHOST, -NGHOST))
  c00 = c000 * (1.0 - dx) + c100 * dx
  c01 = c001 * (1.0 - dx) + c101 * dx
  c10 = c010 * (1.0 - dx) + c110 * dx
  c11 = c011 * (1.0 - dx) + c111 * dx
  c0 = c00 * (1.0 - dy) + c10 * dy
  c1 = c01 * (1.0 - dy) + c11 * dy
  by0 = c0 * (1.0 - dz) + c1 * dz

  !b_z
  c000 = 0.25 * (bz(lind - 1 - iy, -NGHOST, -NGHOST) + bz(lind - 1, -NGHOST, -NGHOST) + &
                 bz(lind - iy, -NGHOST, -NGHOST) + bz(lind, -NGHOST, -NGHOST))
  c100 = 0.25 * (bz(lind - iy, -NGHOST, -NGHOST) + bz(lind, -NGHOST, -NGHOST) + &
                 bz(lind + 1 - iy, -NGHOST, -NGHOST) + bz(lind + 1, -NGHOST, -NGHOST))
  c001 = 0.25 * (bz(lind - 1 - iy + iz, -NGHOST, -NGHOST) + bz(lind - 1 + iz, -NGHOST, -NGHOST) + &
                 bz(lind - iy + iz, -NGHOST, -NGHOST) + bz(lind + iz, -NGHOST, -NGHOST))
  c101 = 0.25 * (bz(lind - iy + iz, -NGHOST, -NGHOST) + bz(lind + iz, -NGHOST, -NGHOST) + &
                 bz(lind + 1 - iy + iz, -NGHOST, -NGHOST) + bz(lind + 1 + iz, -NGHOST, -NGHOST))
  c010 = 0.25 * (bz(lind - 1, -NGHOST, -NGHOST) + bz(lind - 1 + iy, -NGHOST, -NGHOST) + &
                 bz(lind, -NGHOST, -NGHOST) + bz(lind + iy, -NGHOST, -NGHOST))
  c110 = 0.25 * (bz(lind, -NGHOST, -NGHOST) + bz(lind + iy, -NGHOST, -NGHOST) + &
                 bz(lind + 1, -NGHOST, -NGHOST) + bz(lind + 1 + iy, -NGHOST, -NGHOST))
  c011 = 0.25 * (bz(lind - 1 + iz, -NGHOST, -NGHOST) + bz(lind - 1 + iy + iz, -NGHOST, -NGHOST) + &
                 bz(lind + iz, -NGHOST, -NGHOST) + bz(lind + iy + iz, -NGHOST, -NGHOST))
  c111 = 0.25 * (bz(lind + iz, -NGHOST, -NGHOST) + bz(lind + iy + iz, -NGHOST, -NGHOST) + &
                 bz(lind + 1 + iz, -NGHOST, -NGHOST) + bz(lind + 1 + iy + iz, -NGHOST, -NGHOST))
  c00 = c000 * (1.0 - dx) + c100 * dx
  c01 = c001 * (1.0 - dx) + c101 * dx
  c10 = c010 * (1.0 - dx) + c110 * dx
  c11 = c011 * (1.0 - dx) + c111 * dx
  c0 = c00 * (1.0 - dy) + c10 * dy
  c1 = c01 * (1.0 - dy) + c11 * dy
  bz0 = c0 * (1.0 - dz) + c1 * dz
#endif
