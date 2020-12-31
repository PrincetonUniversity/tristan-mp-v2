#include "../defs.F90"

module m_mover
  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_domain
  use m_particles
  use m_fields
  use m_userfile
  use m_exchangearray, only: exchangeArray

  implicit none
contains
  subroutine moveParticles(timestep)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in)                   :: timestep
    integer                               :: s, p, temp_i, ti, tj, tk
    real                                  :: g_temp, over_e_temp, temp_r
    integer(kind=2), pointer, contiguous  :: pt_xi(:), pt_yi(:), pt_zi(:)
    real, pointer, contiguous             :: pt_dx(:), pt_dy(:), pt_dz(:),&
                                           & pt_u(:), pt_v(:), pt_w(:), pt_wei(:)
    real                                  :: ex0, ey0, ez0, bx0, by0, bz0, q_over_m
    real                                  :: u0, v0, w0, u1, v1, w1, dummy_, dummy2_, dx, dy, dz
    logical                               :: dummy_flag
    real                                  :: ex_ext, ey_ext, ez_ext
    real                                  :: bx_ext, by_ext, bz_ext
    real                                  :: c000, c100, c001, c101, c010, c110, c011, c111,&
                                           & c00, c01, c10, c11, c0, c1
    integer                               :: iy, iz, lind

    iy = this_meshblock%ptr%sx + 2 * NGHOST
    iz = iy * (this_meshblock%ptr%sy + 2 * NGHOST)

    do s = 1, nspec
      if (.not. species(s)%move_sp) cycle
      if (species(s)%m_sp .eq. 0) then
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

              ! routine for massless particles
              !$omp simd private(over_e_temp, temp_r, temp_i)
              !dir$ vector aligned
              do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
                ! move particle
                over_e_temp = 1.0 / sqrt(pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
                ! this "function" takes
                ! ... inverse energy: `over_e_temp` ...
                ! ... reads the velocities from: `pt_*(p)` ...
                ! ... and updates the particle position `pt_*(p)`
                include "position_update.F"
              end do ! p
              pt_xi => null();  pt_yi => null();  pt_zi => null()
              pt_dx => null();  pt_dy => null();  pt_dz => null()
              pt_u => null();   pt_v => null();   pt_w => null()
            end do ! tk
          end do ! tj
        end do ! ti
      else ! massive particles
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

              pt_wei => species(s)%prtl_tile(ti, tj, tk)%weight

              ! routine for massive particles
              q_over_m = species(s)%ch_sp / species(s)%m_sp
              #if !defined(EXTERNALFIELDS)
              !$omp simd private(lind, dummy_, dummy2_, g_temp, over_e_temp,&
              !$omp  temp_r, temp_i, u0, v0, w0, u1, v1, w1,&
              !$omp  ex0, ey0, ez0, bx0, by0, bz0,&
              !$omp  c000, c100, c001, c101, c010, c110, c011, c111,&
              !$omp  c00, c01, c10, c11, c0, c1)
              !dir$ vector aligned
              #endif
              do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
                #ifdef oneD
                  lind = pt_xi(p)
                #elif twoD
                  lind = pt_xi(p) + (NGHOST + pt_yi(p)) * iy
                #elif threeD
                  lind = pt_xi(p) + (NGHOST + pt_yi(p)) * iy + (NGHOST + pt_zi(p)) * iz
                #endif
                dx = pt_dx(p); dy = pt_dy(p); dz = pt_dz(p)
                ! these "functions" take
                ! ... coordinates and linear index: `dx`, `dy`, `dz` and `lind` ...
                ! ... and "return" `bx0`, `by0`, `bz0`, `ex0`, `ey0`, `ez0`
                include "interp_efield.F"
                include "interp_bfield.F"

                #if defined (EXTERNALFIELDS)
                  call userExternalFields(REAL(pt_xi(p)) + pt_dx(p),&
                                        & REAL(pt_yi(p)) + pt_dy(p),&
                                        & REAL(pt_zi(p)) + pt_dz(p),&
                                        & ex_ext, ey_ext, ez_ext,&
                                        & bx_ext, by_ext, bz_ext)
                  ex0 = ex0 + ex_ext; ey0 = ey0 + ey_ext; ez0 = ez0 + ez_ext
                  bx0 = bx0 + bx_ext; by0 = by0 + by_ext; bz0 = bz0 + bz_ext
                #endif

                ! . . . . simple Boris pusher . . . .
                u0 = pt_u(p); v0 = pt_v(p); w0 = pt_w(p)
                ! this "function" takes
                ! ... the field quantities: `bx0`, `by0`, `bz0`, `ex0`, `ey0`, `ez0` ...
                ! ... and the velocities: `u0`, `v0`, `w0` ...
                ! ... and returns the updated velocities `u0`, `v0`, `w0`
                #ifdef VAY
                  include "vay_push.F"
                #else
                  include "boris_push.F"
                #endif
                pt_u(p) = u0; pt_v(p) = v0; pt_w(p) = w0
                over_e_temp = 1.0 / sqrt(1.0 + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
                ! this "function" takes
                ! ... inverse energy: `over_e_temp` ...
                ! ... reads the velocities from: `pt_*(p)` ...
                ! ... and updates the particle position `pt_*(p)`
                include "position_update.F"

              end do
              pt_xi => null(); pt_yi => null(); pt_zi => null()
              pt_dx => null(); pt_dy => null(); pt_dz => null()
              pt_u => null(); pt_v => null(); pt_w => null()

              pt_wei => null()
              
            end do ! tk
          end do ! tj
        end do ! ti
      end if
    end do ! species
    call printDiag((mpi_rank .eq. 0), "moveParticles()", .true.)

  end subroutine moveParticles
end module m_mover
