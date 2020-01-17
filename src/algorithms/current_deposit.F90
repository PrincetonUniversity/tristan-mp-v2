#include "../defs.F90"

module m_currentdeposit
  use m_globalnamespace
  use m_aux
  use m_domain
  use m_fields
  use m_particles
  implicit none
contains
  subroutine depositCurrents()
    implicit none
    integer :: s, p, ti, tj, tk
    integer(kind=2), pointer, contiguous  :: pt_xi(:), pt_yi(:), pt_zi(:)
    real, pointer, contiguous             :: pt_dx(:), pt_dy(:), pt_dz(:),&
                                           & pt_u(:), pt_v(:), pt_w(:)
    real                                  :: xr, yr, zr, x1, y1, z1, x2, y2, z2
    real                                  :: gamma_inv, temp_charge
    integer(kind=2)                       :: i1, i2, j1, j2, k1, k2
    integer(kind=2)                       :: i1p1, i2p1, j1p1, j2p1, k1p1, k2p1
    real                                  :: Wx1, Wy1, Wz1, Wx2, Wy2, Wz2
    real                                  :: onemWx1, onemWy1, onemWz1, onemWx2, onemWy2, onemWz2
    real                                  :: Fx1, Fy1, Fz1, Fx2, Fy2, Fz2

    jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0

    do s = 1, nspec ! loop over species
      if (species(s)%ch_sp .eq. 0) cycle
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
            pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
            pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi

            pt_dx => species(s)%prtl_tile(ti, tj, tk)%dx
            pt_dy => species(s)%prtl_tile(ti, tj, tk)%dy
            pt_dz => species(s)%prtl_tile(ti, tj, tk)%dz

            pt_u => species(s)%prtl_tile(ti, tj, tk)%u
            pt_v => species(s)%prtl_tile(ti, tj, tk)%v
            pt_w => species(s)%prtl_tile(ti, tj, tk)%w

            temp_charge = species(s)%ch_sp * unit_ch / B_norm
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              ! push the particle back
              gamma_inv = 1.0 / sqrt(1.0 + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
              #ifndef threeD
                x2 = REAL(pt_xi(p)) + pt_dx(p);       y2 = REAL(pt_yi(p)) + pt_dy(p);       z2 = REAL(pt_zi(p)) + pt_dz(p)
                x1 = x2 - pt_u(p) * CC * gamma_inv;   y1 = y2 - pt_v(p) * CC * gamma_inv;   z1 = z2 - pt_w(p) * CC * gamma_inv

                i1 = INT(x1, 2);  i2 = pt_xi(p)
                j1 = INT(y1, 2);  j2 = pt_yi(p)
                k1 = 0;           k2 = 0
                i1p1 = i1 + 1_2;  i2p1 = i2 + 1_2
                j1p1 = j1 + 1_2;  j2p1 = j2 + 1_2

                xr = min(REAL(min(i1, i2) + 1), max(REAL(max(i1, i2)), 0.5 * (x1 + x2)))
                yr = min(REAL(min(j1, j2) + 1), max(REAL(max(j1, j2)), 0.5 * (y1 + y2)))
                zr = min(REAL(min(k1, k2) + 1), max(REAL(max(k1, k2)), 0.5 * (z1 + z2)))

                Wx1 = 0.5 * (x1 + xr) - i1;   Wy1 = 0.5 * (y1 + yr) - j1
                Wx2 = 0.5 * (x2 + xr) - i2;   Wy2 = 0.5 * (y2 + yr) - j2
                onemWx1 = 1 - Wx1;            onemWy1 = 1 - Wy1
                onemWx2 = 1 - Wx2;            onemWy2 = 1 - Wy2

                ! deposit with a "-" sign
                Fx1 = -temp_charge * (xr - x1); Fy1 = -temp_charge * (yr - y1);  Fz1 = -temp_charge * (zr - z1)
                Fx2 = -temp_charge * (x2 - xr); Fy2 = -temp_charge * (y2 - yr);  Fz2 = -temp_charge * (z2 - zr)

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
              #else
                x2 = REAL(pt_xi(p)) + pt_dx(p);       y2 = REAL(pt_yi(p)) + pt_dy(p);       z2 = REAL(pt_zi(p)) + pt_dz(p)
                x1 = x2 - pt_u(p) * CC * gamma_inv;   y1 = y2 - pt_v(p) * CC * gamma_inv;   z1 = z2 - pt_w(p) * CC * gamma_inv

                i1 = INT(x1, 2);  i2 = pt_xi(p)
                j1 = INT(y1, 2);  j2 = pt_yi(p)
                k1 = INT(z1, 2);  k2 = pt_zi(p)
                i1p1 = i1 + 1_2;  i2p1 = i2 + 1_2
                j1p1 = j1 + 1_2;  j2p1 = j2 + 1_2
                k1p1 = k1 + 1_2;  k2p1 = k2 + 1_2

                xr = min(REAL(min(i1, i2) + 1), max(REAL(max(i1, i2)), 0.5 * (x1 + x2)))
                yr = min(REAL(min(j1, j2) + 1), max(REAL(max(j1, j2)), 0.5 * (y1 + y2)))
                zr = min(REAL(min(k1, k2) + 1), max(REAL(max(k1, k2)), 0.5 * (z1 + z2)))

                Wx1 = 0.5 * (x1 + xr) - i1; Wy1 = 0.5 * (y1 + yr) - j1; Wz1 = 0.5 * (z1 + zr) - k1
                Wx2 = 0.5 * (x2 + xr) - i2; Wy2 = 0.5 * (y2 + yr) - j2; Wz2 = 0.5 * (z2 + zr) - k2
                onemWx1 = 1 - Wx1; onemWy1 = 1 - Wy1; onemWz1 = 1 - Wz1
                onemWx2 = 1 - Wx2; onemWy2 = 1 - Wy2; onemWz2 = 1 - Wz2

                ! deposit with a "-" sign
                Fx1 = -temp_charge * (xr - x1); Fy1 = -temp_charge * (yr - y1); Fz1 = -temp_charge * (zr - z1)
                Fx2 = -temp_charge * (x2 - xr); Fy2 = -temp_charge * (y2 - yr); Fz2 = -temp_charge * (z2 - zr)

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
            end do
            pt_xi => null(); pt_yi => null(); pt_zi => null()
            pt_dx => null(); pt_dy => null(); pt_dz => null()
            pt_u => null(); pt_v => null(); pt_w => null()
          end do
        end do
      end do
    end do ! species loop
  end subroutine depositCurrents
end module m_currentdeposit
