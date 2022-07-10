! this "function" takes `dx`, `dy`, `dz` and `lind` ...
! ... and "returns" `ex0`, `ey0`, `ez0`

! linear indices used instead of 3d
#ifdef oneD
  !e_x
  c0 = 0.5 * (ex(lind, 0, 0) + ex(lind - 1, 0, 0))
  c1 = 0.5 * (ex(lind, 0, 0) + ex(lind + 1, 0, 0))
  ex0 = c0 * (1 - dx) + c1 * dx

  !e_y
  c0 = ey(lind, 0, 0)
  c1 = ey(lind + 1, 0, 0)
  ey0 = c0 * (1 - dx) + c1 * dx

  !e_z
  c0 = ez(lind, 0, 0)
  c1 = ez(lind + 1, 0, 0)
  ez0 = c0 * (1 - dx) + c1 * dx
#elif twoD
  !e_x
  c000 = 0.5 * (ex(lind, -NGHOST, 0) + ex(lind - 1, -NGHOST, 0))
  c100 = 0.5 * (ex(lind, -NGHOST, 0) + ex(lind + 1, -NGHOST, 0))
  c010 = 0.5 * (ex(lind + iy, -NGHOST, 0) + ex(lind + iy - 1, -NGHOST, 0))
  c110 = 0.5 * (ex(lind + iy, -NGHOST, 0) + ex(lind + iy + 1, -NGHOST, 0))
  c00 = c000 * (1 - dx) + c100 * dx
  c10 = c010 * (1 - dx) + c110 * dx
  ex0 = c00 * (1 - dy) + c10 * dy

  !e_y
  c000 = 0.5 * (ey(lind, -NGHOST, 0) + ey(lind - iy, -NGHOST, 0))
  c100 = 0.5 * (ey(lind + 1, -NGHOST, 0) + ey(lind + 1 - iy, -NGHOST, 0))
  c010 = 0.5 * (ey(lind, -NGHOST, 0) + ey(lind + iy, -NGHOST, 0))
  c110 = 0.5 * (ey(lind + 1, -NGHOST, 0) + ey(lind + iy + 1, -NGHOST, 0))
  c00 = c000 * (1 - dx) + c100 * dx
  c10 = c010 * (1 - dx) + c110 * dx
  ey0 = c00 * (1 - dy) + c10 * dy

  !e_z
  c000 = ez(lind, -NGHOST, 0)
  c100 = ez(lind + 1, -NGHOST, 0)
  c010 = ez(lind + iy, -NGHOST, 0)
  c110 = ez(lind + iy + 1, -NGHOST, 0)
  c00 = c000 * (1 - dx) + c100 * dx
  c10 = c010 * (1 - dx) + c110 * dx
  ez0 = c00 * (1 - dy) + c10 * dy
#elif threeD
  !e_x
  c000 = 0.5 * (ex(lind, -NGHOST, -NGHOST) + ex(lind - 1, -NGHOST, -NGHOST))
  c100 = 0.5 * (ex(lind, -NGHOST, -NGHOST) + ex(lind + 1, -NGHOST, -NGHOST))
  c010 = 0.5 * (ex(lind + iy, -NGHOST, -NGHOST) + ex(lind - 1 + iy, -NGHOST, -NGHOST))
  c110 = 0.5 * (ex(lind + iy, -NGHOST, -NGHOST) + ex(lind + 1 + iy, -NGHOST, -NGHOST))
  c00 = c000 * (1.0 - dx) + c100 * dx
  c10 = c010 * (1.0 - dx) + c110 * dx
  c0 = c00 * (1.0 - dy) + c10 * dy
  c001 = 0.5 * (ex(lind + iz, -NGHOST, -NGHOST) + ex(lind - 1 + iz, -NGHOST, -NGHOST))
  c101 = 0.5 * (ex(lind + iz, -NGHOST, -NGHOST) + ex(lind + 1 + iz, -NGHOST, -NGHOST))
  c011 = 0.5 * (ex(lind + iy + iz, -NGHOST, -NGHOST) + ex(lind - 1 + iy + iz, -NGHOST, -NGHOST))
  c111 = 0.5 * (ex(lind + iy + iz, -NGHOST, -NGHOST) + ex(lind + 1 + iy + iz, -NGHOST, -NGHOST))
  c01 = c001 * (1.0 - dx) + c101 * dx
  c11 = c011 * (1.0 - dx) + c111 * dx
  c1 = c01 * (1.0 - dy) + c11 * dy
  ex0 = c0 * (1.0 - dz) + c1 * dz

  !e_y
  c000 = 0.5 * (ey(lind, -NGHOST, -NGHOST) + ey(lind - iy, -NGHOST, -NGHOST))
  c100 = 0.5 * (ey(lind + 1, -NGHOST, -NGHOST) + ey(lind + 1 - iy, -NGHOST, -NGHOST))
  c010 = 0.5 * (ey(lind, -NGHOST, -NGHOST) + ey(lind + iy, -NGHOST, -NGHOST))
  c110 = 0.5 * (ey(lind + 1, -NGHOST, -NGHOST) + ey(lind + 1 + iy, -NGHOST, -NGHOST))
  c00 = c000 * (1.0 - dx) + c100 * dx
  c10 = c010 * (1.0 - dx) + c110 * dx
  c0 = c00 * (1.0 - dy) + c10 * dy
  c001 = 0.5 * (ey(lind + iz, -NGHOST, -NGHOST) + ey(lind - iy + iz, -NGHOST, -NGHOST))
  c101 = 0.5 * (ey(lind + 1 + iz, -NGHOST, -NGHOST) + ey(lind + 1 - iy + iz, -NGHOST, -NGHOST))
  c011 = 0.5 * (ey(lind + iz, -NGHOST, -NGHOST) + ey(lind + iy + iz, -NGHOST, -NGHOST))
  c111 = 0.5 * (ey(lind + 1 + iz, -NGHOST, -NGHOST) + ey(lind + 1 + iy + iz, -NGHOST, -NGHOST))
  c01 = c001 * (1.0 - dx) + c101 * dx
  c11 = c011 * (1.0 - dx) + c111 * dx
  c1 = c01 * (1.0 - dy) + c11 * dy
  ey0 = c0 * (1.0 - dz) + c1 * dz

  !e_z
  c000 = 0.5 * (ez(lind, -NGHOST, -NGHOST) + ez(lind - iz, -NGHOST, -NGHOST))
  c100 = 0.5 * (ez(lind + 1, -NGHOST, -NGHOST) + ez(lind + 1 - iz, -NGHOST, -NGHOST))
  c010 = 0.5 * (ez(lind + iy, -NGHOST, -NGHOST) + ez(lind + iy - iz, -NGHOST, -NGHOST))
  c110 = 0.5 * (ez(lind + 1 + iy, -NGHOST, -NGHOST) + ez(lind + 1 + iy - iz, -NGHOST, -NGHOST))
  c001 = 0.5 * (ez(lind, -NGHOST, -NGHOST) + ez(lind + iz, -NGHOST, -NGHOST))
  c101 = 0.5 * (ez(lind + 1, -NGHOST, -NGHOST) + ez(lind + 1 + iz, -NGHOST, -NGHOST))
  c011 = 0.5 * (ez(lind + iy, -NGHOST, -NGHOST) + ez(lind + iy + iz, -NGHOST, -NGHOST))
  c111 = 0.5 * (ez(lind + 1 + iy, -NGHOST, -NGHOST) + ez(lind + 1 + iy + iz, -NGHOST, -NGHOST))
  c00 = c000 * (1.0 - dx) + c100 * dx
  c01 = c001 * (1.0 - dx) + c101 * dx
  c10 = c010 * (1.0 - dx) + c110 * dx
  c11 = c011 * (1.0 - dx) + c111 * dx
  c0 = c00 * (1.0 - dy) + c10 * dy
  c1 = c01 * (1.0 - dy) + c11 * dy
  ez0 = c0 * (1.0 - dz) + c1 * dz
#endif
