module m_mover
  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_domain
  use m_particles
  use m_fields
  use m_userfile
  use m_exchangearray, only: exchangeArray

  ! extra physics
#ifdef RADIATION
  use m_radiation
#endif

  implicit none
contains
  subroutine moveParticles(timestep)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in) :: timestep
    integer :: s, p, ti, tj, tk
    integer(kind=2) :: temp_i
    real :: g_temp, over_e_temp, temp_r
    integer(kind=2), pointer, contiguous :: pt_xi(:), pt_yi(:), pt_zi(:)
    real, pointer, contiguous :: pt_dx(:), pt_dy(:), pt_dz(:), &
                                 pt_u(:), pt_v(:), pt_w(:), pt_wei(:)
    real :: ex0, ey0, ez0, bx0, by0, bz0, q_over_m
    real :: u0, v0, w0, u1, v1, w1, dummy_, dx, dy, dz
#if defined (EXTERNALFIELDS)
    real :: ex_ext, ey_ext, ez_ext
    real :: bx_ext, by_ext, bz_ext
#endif
    real :: c000, c100, c001, c101, c010, c110, c011, c111, &
            c00, c01, c10, c11, c0, c1
    integer :: iy, iz, lind

#ifdef PRTLPAYLOADS
    real, pointer, contiguous :: pt_pld1(:), pt_pld2(:), pt_pld3(:)
    real :: incr_pld1, incr_pld2, incr_pld3
#endif

#if defined (RADIATION) || defined (GCA)
    integer, pointer, contiguous :: pt_proc(:)
#endif

#ifdef RADIATION
    logical :: dummy_flag
    real :: ex_rad, ey_rad, ez_rad, bx_rad, by_rad, bz_rad
    real :: u_init, v_init, w_init
    integer, pointer, contiguous :: pt_ind(:)
#endif

    if (.false.) print *, timestep

    iy = this_meshblock % ptr % i2 - this_meshblock % ptr % i1 + 1
    iz = iy * (this_meshblock % ptr % j2 - this_meshblock % ptr % j1 + 1)

#ifdef RADIATION
    if (rad_dens_lim .gt. 0) then
      ! if density limit is enabled
      dummy_flag = .true.
      do s = 1, nspec
        if (species(s) % m_sp .ne. 0) then
          call computeDensity(s, reset=dummy_flag)
          dummy_flag = .false.
        end if
      end do
      call exchangeArray()
    end if
#endif

    do s = 1, nspec
      if (.not. species(s) % move_sp) cycle
      if (species(s) % m_sp .eq. 0) then
        do tk = 1, species(s) % tile_nz
          do tj = 1, species(s) % tile_ny
            do ti = 1, species(s) % tile_nx
              pt_xi => species(s) % prtl_tile(ti, tj, tk) % xi
              pt_yi => species(s) % prtl_tile(ti, tj, tk) % yi
              pt_zi => species(s) % prtl_tile(ti, tj, tk) % zi

              pt_dx => species(s) % prtl_tile(ti, tj, tk) % dx
              pt_dy => species(s) % prtl_tile(ti, tj, tk) % dy
              pt_dz => species(s) % prtl_tile(ti, tj, tk) % dz

              pt_u => species(s) % prtl_tile(ti, tj, tk) % u
              pt_v => species(s) % prtl_tile(ti, tj, tk) % v
              pt_w => species(s) % prtl_tile(ti, tj, tk) % w

#ifdef PRTLPAYLOADS
              pt_pld1 => species(s) % prtl_tile(ti, tj, tk) % payload1
              pt_pld2 => species(s) % prtl_tile(ti, tj, tk) % payload2
              pt_pld3 => species(s) % prtl_tile(ti, tj, tk) % payload3
#endif

              ! routine for massless particles
#ifdef PRTLPAYLOADS
              !$omp simd private(over_e_temp, temp_r, temp_i, incr_pld1, incr_pld2, incr_pld3, u0, v0, w0)
#else
              !$omp simd private(over_e_temp, temp_r, temp_i)
#endif
              !dir$ vector aligned
              do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
                ! move particle
                over_e_temp = 1.0 / sqrt(pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
#ifdef PRTLPAYLOADS
                u0 = pt_u(p); v0 = pt_v(p); w0 = pt_w(p)
                call usrSetPhPld(u0, v0, w0, over_e_temp, incr_pld1, incr_pld2, incr_pld3)
                pt_pld1(p) = pt_pld1(p) + incr_pld1
                pt_pld2(p) = pt_pld2(p) + incr_pld2
                pt_pld3(p) = pt_pld3(p) + incr_pld3
#endif
                ! this "function" takes
                ! ... inverse energy: `over_e_temp` ...
                ! ... reads the velocities from: `pt_*(p)` ...
                ! ... and updates the particle position `pt_*(p)`
#include "position_update.F08"
              end do ! p
              pt_xi => null(); pt_yi => null(); pt_zi => null()
              pt_dx => null(); pt_dy => null(); pt_dz => null()
              pt_u => null(); pt_v => null(); pt_w => null()
#ifdef PRTLPAYLOADS
              pt_pld1 => null(); pt_pld2 => null(); pt_pld3 => null()
#endif
            end do ! ti
          end do ! tj
        end do ! tk
      else ! massive particles
        do tk = 1, species(s) % tile_nz
          do tj = 1, species(s) % tile_ny
            do ti = 1, species(s) % tile_nx
              pt_xi => species(s) % prtl_tile(ti, tj, tk) % xi
              pt_yi => species(s) % prtl_tile(ti, tj, tk) % yi
              pt_zi => species(s) % prtl_tile(ti, tj, tk) % zi

              pt_dx => species(s) % prtl_tile(ti, tj, tk) % dx
              pt_dy => species(s) % prtl_tile(ti, tj, tk) % dy
              pt_dz => species(s) % prtl_tile(ti, tj, tk) % dz

              pt_u => species(s) % prtl_tile(ti, tj, tk) % u
              pt_v => species(s) % prtl_tile(ti, tj, tk) % v
              pt_w => species(s) % prtl_tile(ti, tj, tk) % w

              pt_wei => species(s) % prtl_tile(ti, tj, tk) % weight

#ifdef PRTLPAYLOADS
              pt_pld1 => species(s) % prtl_tile(ti, tj, tk) % payload1
              pt_pld2 => species(s) % prtl_tile(ti, tj, tk) % payload2
              pt_pld3 => species(s) % prtl_tile(ti, tj, tk) % payload3
#endif

#if defined(RADIATION) || defined(GCA)
              pt_proc => species(s) % prtl_tile(ti, tj, tk) % proc
#endif

#if defined(RADIATION)
              pt_ind => species(s) % prtl_tile(ti, tj, tk) % ind
#endif

              ! routine for massive particles
              q_over_m = species(s) % ch_sp / species(s) % m_sp
#if !defined(RADIATION) && !defined(EXTERNALFIELDS) && !defined(GCA)
              !$omp simd private(lind, dummy_, g_temp, over_e_temp,&
              !$omp  temp_r, temp_i, u0, v0, w0, u1, v1, w1,&
              !$omp  ex0, ey0, ez0, bx0, by0, bz0,&
#ifdef PRTLPAYLOADS
              !$omp  incr_pld1, incr_pld2, incr_pld3,&
#endif
              !$omp  c000, c100, c001, c101, c010, c110, c011, c111,&
              !$omp  c00, c01, c10, c11, c0, c1)
              !dir$ vector aligned
#endif
              do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
#ifndef DEBUG
                ! use "fast" interpolation
#ifdef oneD
                lind = pt_xi(p)
#elif defined(twoD)
                lind = pt_xi(p) + (NGHOST + pt_yi(p)) * iy
#elif defined(threeD)
                lind = pt_xi(p) + (NGHOST + pt_yi(p)) * iy + (NGHOST + pt_zi(p)) * iz
#endif
                dx = pt_dx(p); dy = pt_dy(p); dz = pt_dz(p)
                ! these "functions" take
                ! ... coordinates and linear index: `dx`, `dy`, `dz` and `lind` ...
                ! ... and "return" `bx0`, `by0`, `bz0`, `ex0`, `ey0`, `ez0`
#include "interp_efield.F08"
#include "interp_bfield.F08"
#else
                dx = pt_dx(p); dy = pt_dy(p); dz = pt_dz(p)
                call interpFromEdges(dx, dy, dz, pt_xi(p), pt_yi(p), pt_zi(p), &
                                     ex, ey, ez, ex0, ey0, ez0)
                call interpFromFaces(dx, dy, dz, pt_xi(p), pt_yi(p), pt_zi(p), &
                                     bx, by, bz, bx0, by0, bz0)
#endif

#if defined (EXTERNALFIELDS) && !defined (GCA)
                call userExternalFields(REAL(pt_xi(p)) + pt_dx(p), &
                                        REAL(pt_yi(p)) + pt_dy(p), &
                                        REAL(pt_zi(p)) + pt_dz(p), &
                                        ex_ext, ey_ext, ez_ext, &
                                        bx_ext, by_ext, bz_ext)
                ex0 = ex0 + ex_ext; ey0 = ey0 + ey_ext; ez0 = ez0 + ez_ext
                bx0 = bx0 + bx_ext; by0 = by0 + by_ext; bz0 = bz0 + bz_ext
#endif

#ifdef RADIATION
                ! save fields at time `t = n`
                ex_rad = ex0; ey_rad = ey0; ez_rad = ez0
                bx_rad = bx0; by_rad = by0; bz_rad = bz0

                ! save velocities before the push
                u_init = pt_u(p); v_init = pt_v(p); w_init = pt_w(p)
#endif

                ! . . . . simple Boris/Vay pusher . . . .
! IN:
!     pt_u(p), pt_v(p), pt_w(p):      velocities
!     ex0, ey0, ez0:                  electric field
!     bx0, by0, bz0:                  magnetic field
! IN (RADIATION):
!     pt_xi(p), pt_yi(p), pt_zi(p):   coordinates
!     pt_dx(p), pt_dy(p), pt_dz(p):   cell offsets
!     u_init, v_init, w_init:         velocities before the push
!     ex_rad, ey_rad, ez_rad:         electric field (backup)
!     bx_rad, by_rad, bz_rad:         magnetic field (backup)
! OUT:
!     pt_u(p), pt_v(p), pt_w(p):      velocities
! NOTE:
!     fields ex0, ey0, ez0, bx0, by0, bz0 get modified
#include "momentum_update.F08"
                over_e_temp = 1.0 / sqrt(1.0 + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
#ifdef PRTLPAYLOADS
                ! !TODO: using ex0, ey0, ez0, bx0, by0, bz0 is not correct here
                call usrSetElPld(q_over_m, pt_u(p), pt_v(p), pt_w(p), over_e_temp, ex0, ey0, ez0, bx0, by0, bz0, incr_pld1, incr_pld2, incr_pld3)
                pt_pld1(p) = pt_pld1(p) + incr_pld1
                pt_pld2(p) = pt_pld2(p) + incr_pld2
                pt_pld3(p) = pt_pld3(p) + incr_pld3
#endif
                ! IN:
                !     over_e_temp:                    inverse of energy
                !     pt_u(p), pt_v(p), pt_w(p):      velocities
                ! OUT:
                !     pt_xi(p), pt_yi(p), pt_zi(p):   coordinates
                !     pt_dx(p), pt_dy(p), pt_dz(p):   cell offsets
#include "position_update.F08"

              end do
              pt_xi => null(); pt_yi => null(); pt_zi => null()
              pt_dx => null(); pt_dy => null(); pt_dz => null()
              pt_u => null(); pt_v => null(); pt_w => null()

              pt_wei => null()

#ifdef PRTLPAYLOADS
              pt_pld1 => null(); pt_pld2 => null(); pt_pld3 => null()
#endif

#if defined (RADIATION) || defined (GCA)
              pt_proc => null(); 
#endif

#if defined(RADIATION)
              pt_ind => null()
#endif

            end do ! ti
          end do ! tj
        end do ! tk
      end if
    end do ! species
    call printDiag("moveParticles()", 2)

  end subroutine moveParticles
end module m_mover

