#include "../defs.F90"

module m_mover
  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_domain
  use m_particles
  use m_fields
  use m_userfile

  ! extra physics
  #ifdef RADIATION
    use m_radiation
  #endif

  implicit none
contains
  subroutine moveParticles()
    implicit none
    integer                               :: s, p, temp_i, ti, tj, tk
    real                                  :: g_temp, over_g_temp, over_e_temp, temp_r
    integer(kind=2), pointer, contiguous  :: pt_xi(:), pt_yi(:), pt_zi(:)
    real, pointer, contiguous             :: pt_dx(:), pt_dy(:), pt_dz(:),&
                                           & pt_u(:), pt_v(:), pt_w(:)
    real                                  :: ex0, ey0, ez0, bx0, by0, bz0, q_over_m
    real                                  :: u0, v0, w0, u1, v1, w1, dummy_
    real                                  :: ex_rad, ey_rad, ez_rad, bx_rad, by_rad, bz_rad
    real                                  :: u_init, v_init, w_init, du_rad, dv_rad, dw_rad
    logical                               :: dummy_flag
    real                                  :: ex_ext, ey_ext, ez_ext
    real                                  :: bx_ext, by_ext, bz_ext

    #ifdef RADIATION
      dummy_flag = .true.
      do s = 1, nspec
        if (species(s)%m_sp .ne. 0) then
          call computeDensity(s, reset=dummy_flag)
          dummy_flag = .false.
        end if
      end do
      call exchangeArray()
    #endif

    do s = 1, nspec
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

            if (species(s)%m_sp .eq. 0) then
              if (species(s)%ch_sp .ne. 0) then
                ! if massless but charged -> don't do anything
                !   can be used for, e.g., stationary ions
                cycle
              end if
              ! routine for massless particles
              ! !$omp simd
              !dir$ vector aligned
              do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
                ! move particle
                over_e_temp = 1.0 / sqrt(pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)

                pt_dx(p) = pt_dx(p) + CC * pt_u(p) * over_e_temp
                temp_i = INT(pt_dx(p))
                temp_r = MAX(SIGN(1., pt_dx(p)) + temp_i, REAL(temp_i)) - 1
                temp_i = INT(temp_r)
                pt_xi(p) = pt_xi(p) + temp_i
                pt_dx(p) = pt_dx(p) - temp_r

                pt_dy(p) = pt_dy(p) + CC * pt_v(p) * over_e_temp
                temp_i = INT(pt_dy(p))
                temp_r = MAX(SIGN(1., pt_dy(p)) + temp_i, REAL(temp_i)) - 1
                temp_i = INT(temp_r)
                pt_yi(p) = pt_yi(p) + temp_i
                pt_dy(p) = pt_dy(p) - temp_r

                #ifdef threeD
                  pt_dz(p) = pt_dz(p) + CC * pt_w(p) * over_e_temp
                  temp_i = INT(pt_dz(p))
                  temp_r = MAX(SIGN(1., pt_dz(p)) + temp_i, REAL(temp_i)) - 1
                  temp_i = INT(temp_r)
                  pt_zi(p) = pt_zi(p) + temp_i
                  pt_dz(p) = pt_dz(p) - temp_r
                #endif
              end do
            else
              ! routine for massive particles
              q_over_m = species(s)%ch_sp / species(s)%m_sp
              ! !$omp simd
              !dir$ vector aligned
              do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp

                call interpFromEdges(pt_dx(p), pt_dy(p), pt_dz(p),&
                                   & pt_xi(p), pt_yi(p), pt_zi(p),&
                                   & ex, ey, ez, ex0, ey0, ez0)
                call interpFromFaces(pt_dx(p), pt_dy(p), pt_dz(p),&
                                   & pt_xi(p), pt_yi(p), pt_zi(p),&
                                   & bx, by, bz, bx0, by0, bz0)

                #ifdef EXTERNALFIELDS
                  call userExternalFields(REAL(pt_xi(p)) + pt_dx(p),&
                                        & REAL(pt_yi(p)) + pt_dy(p),&
                                        & REAL(pt_zi(p)) + pt_dz(p),&
                                        & ex_ext, ey_ext, ez_ext,&
                                        & bx_ext, by_ext, bz_ext)
                  ex0 = ex0 + ex_ext; ey0 = ey0 + ey_ext; ez0 = ez0 + ez_ext
                  bx0 = bx0 + bx_ext; by0 = by0 + by_ext; bz0 = bz0 + bz_ext
                #endif

                #ifdef RADIATION
                  ex_rad = ex0; ey_rad = ey0; ez_rad = ez0
                  bx_rad = bx0; by_rad = by0; bz_rad = bz0

                  u_init = pt_u(p)
                  v_init = pt_v(p)
                  w_init = pt_w(p)
                #endif

                dummy_ = 0.5 * q_over_m * B_norm
                ex0 = ex0 * dummy_; ey0 = ey0 * dummy_; ez0 = ez0 * dummy_
                dummy_ = dummy_ * CCINV
                bx0 = bx0 * dummy_; by0 = by0 * dummy_; bz0 = bz0 * dummy_

                ! BORIS PUSHER >
                ! half acceleration:
                u0 = CC * pt_u(p) + ex0
                v0 = CC * pt_v(p) + ey0
                w0 = CC * pt_w(p) + ez0

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

                u0 = u0 + v1 * bz0 - w1 * by0 + ex0
                v0 = v0 + w1 * bx0 - u1 * bz0 + ey0
                w0 = w0 + u1 * by0 - v1 * bx0 + ez0
                ! </ BORIS PUSHER

                pt_u(p) = u0 * CCINV
                pt_v(p) = v0 * CCINV
                pt_w(p) = w0 * CCINV

                ! RADIATION >
                #ifdef RADIATION
                  #ifdef SYNCHROTRON
                    if (species(s)%cool_sp) then
                     call particleRadiateSync(s,&
                                            & pt_u(p), pt_v(p), pt_w(p), u_init, v_init, w_init,&
                                            & pt_dx(p), pt_dy(p), pt_dz(p), pt_xi(p), pt_yi(p), pt_zi(p),&
                                            & bx_rad, by_rad, bz_rad, ex_rad, ey_rad, ez_rad)
                    end if
                  #endif
                  #ifdef INVERSECOMPTON
                    if (species(s)%cool_sp) then
                     call particleRadiateIC(s,&
                                          & pt_u(p), pt_v(p), pt_w(p), u_init, v_init, w_init,&
                                          & pt_dx(p), pt_dy(p), pt_dz(p), pt_xi(p), pt_yi(p), pt_zi(p),&
                                          & bx_rad, by_rad, bz_rad, ex_rad, ey_rad, ez_rad)
                    end if
                  #endif
                #endif
                ! </ RADIATION

                ! move particle
                g_temp = sqrt(1.0 + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
                over_g_temp = 1.0 / g_temp

                pt_dx(p) = pt_dx(p) + CC * pt_u(p) * over_g_temp
                temp_i = INT(pt_dx(p))
                temp_r = MAX(SIGN(1.0, pt_dx(p)) + temp_i, REAL(temp_i)) - 1
                temp_i = INT(temp_r)
                pt_xi(p) = pt_xi(p) + temp_i
                pt_dx(p) = pt_dx(p) - temp_r

                pt_dy(p) = pt_dy(p) + CC * pt_v(p) * over_g_temp
                temp_i = INT(pt_dy(p))
                temp_r = MAX(SIGN(1.0, pt_dy(p)) + temp_i, REAL(temp_i)) - 1
                temp_i = INT(temp_r)
                pt_yi(p) = pt_yi(p) + temp_i
                pt_dy(p) = pt_dy(p) - temp_r

                #ifdef threeD
                  pt_dz(p) = pt_dz(p) + CC * pt_w(p) * over_g_temp
                  temp_i = INT(pt_dz(p))
                  temp_r = MAX(SIGN(1.0, pt_dz(p)) + temp_i, REAL(temp_i)) - 1
                  temp_i = INT(temp_r)
                  pt_zi(p) = pt_zi(p) + temp_i
                  pt_dz(p) = pt_dz(p) - temp_r
                #endif
              end do
            end if
            pt_xi => null(); pt_yi => null(); pt_zi => null()
            pt_dx => null(); pt_dy => null(); pt_dz => null()
            pt_u => null(); pt_v => null(); pt_w => null()
          end do ! tk
        end do ! tj
      end do ! ti
    end do ! species
    call printDiag((mpi_rank .eq. 0), "moveParticles()", .true.)
  end subroutine moveParticles
end module m_mover
