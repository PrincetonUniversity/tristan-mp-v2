! this "function" takes
! ... the start and end coordinates: `x1`, `x2`, `y1`, `y2`, `z1`, `z2` ...
! ... the start and end cells: `i1`, `i2`, `j1`, `j2`, `k1`, `k2` ...
! ... the start and end cells + 1: `i1p1`, `i2p1` etc ...
! ... the weighted_chargeed charge: `weighted_charge = weight * charge_sp * unit_charge / Bnorm`
! ... and deposits proper currents to corresponding components

xr = min(REAL(min(i1, i2) + 1), max(REAL(max(i1, i2)), 0.5 * (x1 + x2)))
yr = min(REAL(min(j1, j2) + 1), max(REAL(max(j1, j2)), 0.5 * (y1 + y2)))
zr = min(REAL(min(k1, k2) + 1), max(REAL(max(k1, k2)), 0.5 * (z1 + z2)))

#ifdef oneD
  Wx1 = 0.5 * (x1 + xr) - i1
  Wx2 = 0.5 * (x2 + xr) - i2
  onemWx1 = 1 - Wx1
  onemWx2 = 1 - Wx2

  ! deposit with a "-" sign
  Fx1 = -(xr - x1) * weighted_charge
  Fy1 = -(yr - y1) * weighted_charge
  Fz1 = -(zr - z1) * weighted_charge
  Fx2 = -(x2 - xr) * weighted_charge
  Fy2 = -(y2 - yr) * weighted_charge
  Fz2 = -(z2 - zr) * weighted_charge

  jx(i1  , j1  , k1) = jx(i1  , j1  , k1) + Fx1
  jx(i2  , j2  , k2) = jx(i2  , j2  , k2) + Fx2

  jy(i1  , j1  , k1) = jy(i1  , j1  , k1) + Fy1 * onemWx1
  jy(i1p1, j1  , k1) = jy(i1p1, j1  , k1) + Fy1 * Wx1
  jy(i2  , j2  , k2) = jy(i2  , j2  , k2) + Fy2 * onemWx2
  jy(i2p1, j2  , k2) = jy(i2p1, j2  , k2) + Fy2 * Wx2

  jz(i1  , j1  , k1) = jz(i1  , j1  , k1) + Fz1 * onemWx1
  jz(i1p1, j1  , k1) = jz(i1p1, j1  , k1) + Fz1 * Wx1
  jz(i2  , j2  , k2) = jz(i2  , j2  , k2) + Fz2 * onemWx2
  jz(i2p1, j2  , k2) = jz(i2p1, j2  , k2) + Fz2 * Wx2
#elif twoD
  Wx1 = 0.5 * (x1 + xr) - i1;   Wy1 = 0.5 * (y1 + yr) - j1
  Wx2 = 0.5 * (x2 + xr) - i2;   Wy2 = 0.5 * (y2 + yr) - j2
  onemWx1 = 1 - Wx1;            onemWy1 = 1 - Wy1
  onemWx2 = 1 - Wx2;            onemWy2 = 1 - Wy2

  ! deposit with a "-" sign
  Fx1 = -(xr - x1) * weighted_charge
  Fy1 = -(yr - y1) * weighted_charge
  Fz1 = -(zr - z1) * weighted_charge
  Fx2 = -(x2 - xr) * weighted_charge
  Fy2 = -(y2 - yr) * weighted_charge
  Fz2 = -(z2 - zr) * weighted_charge

  jx(i1  , j1  , k1) = jx(i1  , j1  , k1) + Fx1 * onemWy1
  jx(i1  , j1p1, k1) = jx(i1  , j1p1, k1) + Fx1 * Wy1

  jy(i1  , j1  , k1) = jy(i1  , j1  , k1) + Fy1 * onemWx1
  jy(i1p1, j1  , k1) = jy(i1p1, j1  , k1) + Fy1 * Wx1

  jx(i2  , j2  , k2) = jx(i2  , j2  , k2) + Fx2 * onemWy2
  jx(i2  , j2p1, k2) = jx(i2  , j2p1, k2) + Fx2 * Wy2

  jy(i2  , j2  , k2) = jy(i2  , j2  , k2) + Fy2 * onemWx2
  jy(i2p1, j2  , k2) = jy(i2p1, j2  , k2) + Fy2 * Wx2

  jz(i1  , j1  , k1) = jz(i1  , j1  , k1) + Fz1 * onemWx1 * onemWy1
  jz(i1p1, j1  , k1) = jz(i1p1, j1  , k1) + Fz1 * Wx1 * onemWy1
  jz(i1  , j1p1, k1) = jz(i1  , j1p1, k1) + Fz1 * onemWx1 * Wy1
  jz(i1p1, j1p1, k1) = jz(i1p1, j1p1, k1) + Fz1 * Wx1 * Wy1

  jz(i2  , j2  , k2) = jz(i2  , j2  , k2) + Fz2 * onemWx2 * onemWy2
  jz(i2p1, j2  , k2) = jz(i2p1, j2  , k2) + Fz2 * Wx2 * onemWy2
  jz(i2  , j2p1, k2) = jz(i2  , j2p1, k2) + Fz2 * onemWx2 * Wy2
  jz(i2p1, j2p1, k2) = jz(i2p1, j2p1, k2) + Fz2 * Wx2 * Wy2
#elif threeD
  Wx1 = 0.5 * (x1 + xr) - i1;   Wy1 = 0.5 * (y1 + yr) - j1;   Wz1 = 0.5 * (z1 + zr) - k1
  Wx2 = 0.5 * (x2 + xr) - i2;   Wy2 = 0.5 * (y2 + yr) - j2;   Wz2 = 0.5 * (z2 + zr) - k2
  onemWx1 = 1 - Wx1;            onemWy1 = 1 - Wy1;            onemWz1 = 1 - Wz1
  onemWx2 = 1 - Wx2;            onemWy2 = 1 - Wy2;            onemWz2 = 1 - Wz2

  ! deposit with a "-" sign
  Fx1 = -(xr - x1) * weighted_charge
  Fy1 = -(yr - y1) * weighted_charge
  Fz1 = -(zr - z1) * weighted_charge
  Fx2 = -(x2 - xr) * weighted_charge
  Fy2 = -(y2 - yr) * weighted_charge
  Fz2 = -(z2 - zr) * weighted_charge

  jx(i1  , j1  , k1  ) = jx(i1  , j1  , k1  ) + Fx1 * onemWy1 * onemWz1
  jx(i1  , j1p1, k1  ) = jx(i1  , j1p1, k1  ) + Fx1 * Wy1 * onemWz1
  jx(i1  , j1  , k1p1) = jx(i1  , j1  , k1p1) + Fx1 * onemWy1 * Wz1
  jx(i1  , j1p1, k1p1) = jx(i1  , j1p1, k1p1) + Fx1 * Wy1 * Wz1

  jy(i1  , j1  , k1  ) = jy(i1  , j1  , k1  ) + Fy1 * onemWx1 * onemWz1
  jy(i1p1, j1  , k1  ) = jy(i1p1, j1  , k1  ) + Fy1 * Wx1 * onemWz1
  jy(i1  , j1  , k1p1) = jy(i1  , j1  , k1p1) + Fy1 * onemWx1 * Wz1
  jy(i1p1, j1  , k1p1) = jy(i1p1, j1  , k1p1) + Fy1 * Wx1 * Wz1

  jz(i1  , j1  , k1  ) = jz(i1  , j1  , k1  ) + Fz1 * onemWx1 * onemWy1
  jz(i1p1, j1  , k1  ) = jz(i1p1, j1  , k1  ) + Fz1 * Wx1 * onemWy1
  jz(i1  , j1p1, k1  ) = jz(i1  , j1p1, k1  ) + Fz1 * onemWx1 * Wy1
  jz(i1p1, j1p1, k1  ) = jz(i1p1, j1p1, k1  ) + Fz1 * Wx1 * Wy1

  jx(i2  , j2  , k2  ) = jx(i2  , j2  , k2  ) + Fx2 * onemWy2 * onemWz2
  jx(i2  , j2p1, k2  ) = jx(i2  , j2p1, k2  ) + Fx2 * Wy2 * onemWz2
  jx(i2  , j2  , k2p1) = jx(i2  , j2  , k2p1) + Fx2 * onemWy2 * Wz2
  jx(i2  , j2p1, k2p1) = jx(i2  , j2p1, k2p1) + Fx2 * Wy2 * Wz2

  jy(i2  , j2  , k2  ) = jy(i2  , j2  , k2  ) + Fy2 * onemWx2 * onemWz2
  jy(i2p1, j2  , k2  ) = jy(i2p1, j2  , k2  ) + Fy2 * Wx2 * onemWz2
  jy(i2  , j2  , k2p1) = jy(i2  , j2  , k2p1) + Fy2 * onemWx2 * Wz2
  jy(i2p1, j2  , k2p1) = jy(i2p1, j2  , k2p1) + Fy2 * Wx2 * Wz2

  jz(i2  , j2  , k2  ) = jz(i2  , j2  , k2  ) + Fz2 * onemWx2 * onemWy2
  jz(i2p1, j2  , k2  ) = jz(i2p1, j2  , k2  ) + Fz2 * Wx2 * onemWy2
  jz(i2  , j2p1, k2  ) = jz(i2  , j2p1, k2  ) + Fz2 * onemWx2 * Wy2
  jz(i2p1, j2p1, k2  ) = jz(i2p1, j2p1, k2  ) + Fz2 * Wx2 * Wy2
#endif
