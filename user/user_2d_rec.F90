#include "../src/defs.F90"

module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real    :: nCS_over_nUP, current_width, upstream_T, cs_x
  real    :: injector_x1, injector_x2, injector_sx, injector_betax
  integer :: injector_reset_interval

  private :: nCS_over_nUP, current_width, upstream_T, cs_x
  private :: injector_x1, injector_x2, injector_sx, injector_betax
  private :: injector_reset_interval
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userInitParticles, userInitFields, userReadInput,&
           & userSpatialDistribution
  !...............................................................!
contains
  subroutine userInitialize()
    implicit none
    call userReadInput()
    call userInitParticles()
    call userInitFields()
  end subroutine userInitialize

  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'upstream_T', upstream_T)
    call getInput('problem', 'nCS_nUP', nCS_over_nUP)
    call getInput('problem', 'current_width', current_width)
    call getInput('problem', 'injector_sx', injector_sx)
    call getInput('problem', 'injector_betax', injector_betax)
    cs_x = 0.5
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob,&
                                 & dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional  :: x_glob, y_glob, z_glob
    real, intent(in), optional  :: dummy1, dummy2, dummy3
    if (present(x_glob) .and. present(dummy1) .and. present(dummy2)) then
      userSpatialDistribution = 1.0 / (cosh((x_glob - dummy1) / dummy2))**2
    else
      call throwError("ERROR: variable not present in `userSpatialDistribution()`")
    end if
    return
  end function userSpatialDistribution

  subroutine userInitParticles()
    implicit none
    real                :: nUP, nCS
    type(region)        :: back_region
    real                :: sx_glob, sy_glob, shift_gamma, shift_beta, current_sheet_T
    procedure (spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    nUP = 0.5 * ppc0
    nCS = nUP * nCS_over_nUP

    sx_glob = REAL(global_mesh%sx)
    sy_glob = REAL(global_mesh%sy)
    
    back_region%x_min = 0.0
    back_region%y_min = 0.0
    back_region%x_max = sx_glob * cs_x
    back_region%y_max = sy_glob
    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, nUP, upstream_T)
    back_region%x_min = sx_glob * cs_x
    back_region%y_min = 0.0
    back_region%x_max = sx_glob
    back_region%y_max = sy_glob
    call fillRegionWithThermalPlasma(back_region, (/3, 4/), 2, nUP, upstream_T)

    shift_beta = sqrt(sigma) * c_omp / (current_width * nCS_over_nUP)
    if (shift_beta .ge. 1) then
      call throwError('ERROR: `shift_beta` >= 1 in `userInitParticles()`')
    end if
    shift_gamma = 1.0 / sqrt(1.0 - shift_beta**2)
    current_sheet_T = 0.5 * sigma / nCS_over_nUP

    back_region%x_min = sx_glob * cs_x - 10 * current_width
    back_region%x_max = sx_glob * cs_x + 10 * current_width
    back_region%y_min = 0.0
    back_region%y_max = sy_glob
    call fillRegionWithThermalPlasma(back_region, (/5, 6/), 2, nCS, current_sheet_T,&
                                   & shift_gamma = shift_gamma, shift_dir = 3,&
                                   & spat_distr_ptr = spat_distr_ptr,&
                                   & dummy1 = cs_x * sx_glob, dummy2 = current_width)
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob
    real    :: x_glob, sx_glob
    ex(:,:,:) = 0; ey(:,:,:) = 0; ez(:,:,:) = 0
    bx(:,:,:) = 0; by(:,:,:) = 0; bz(:,:,:) = 0
    jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0

    injector_x1 = injector_sx - 1.0e-5
    injector_x2 = REAL(global_mesh%sx) - injector_sx + 1.0e-5
    injector_reset_interval = INT(injector_sx / (injector_betax * CC))

    k = 0
    sx_glob = REAL(global_mesh%sx)
    do i = -NGHOST, this_meshblock%ptr%sx - 1 + NGHOST
      i_glob = i + this_meshblock%ptr%x0
      x_glob = REAL(i_glob)
      !do j = -NGHOST, this_meshblock%ptr%sy - 1 + NGHOST
      by(i,:,:) = tanh((x_glob - cs_x * sx_glob) / current_width)
      !end do
    end do
  end subroutine userInitFields
  !............................................................!

  !--- driving ------------------------------------------------!
  subroutine userDriveParticles(step)
    implicit none
    integer, optional, intent(in) :: step
    ! ... dummy loop ...
    ! integer :: s, ti, tj, tk, p
    ! do s = 1, nspec
    !   do ti = 1, species(s)%tile_nx
    !     do tj = 1, species(s)%tile_ny
    !       do tk = 1, species(s)%tile_nz
    !         do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
    !           ...
    !         end do
    !       end do
    !     end do
    !   end do
    ! end do
  end subroutine userDriveParticles
  
  subroutine userExternalFields(xp, yp, zp,&
                              & ex_ext, ey_ext, ez_ext,&
                              & bx_ext, by_ext, bz_ext)
    implicit none
    real, intent(in)  :: xp, yp, zp
    real, intent(out) :: ex_ext, ey_ext, ez_ext
    real, intent(out) :: bx_ext, by_ext, bz_ext
    ! some functions of xp, yp, zp
    ex_ext = 0.0; ey_ext = 0.0; ez_ext = 0.0
    bx_ext = 0.0; by_ext = 0.0; bz_ext = 0.0
  end subroutine userExternalFields
  !............................................................!

  !--- boundaries ---------------------------------------------!
  subroutine userParticleBoundaryConditions(step)
    implicit none
    real                            :: nUP, old_x1, old_x2, x_glob
    integer                         :: s, ti, tj, tk, p, nUP_tot
    integer                         :: injector_i1_glob, injector_i2_glob
    type(region)                    :: back_region
    integer, optional, intent(in)             :: step
    procedure (spatialDistribution), pointer  :: spat_distr_ptr => null()
    
    ! reset the injector position every once in a while
    if ((modulo(step, injector_reset_interval) .eq. 0) .and. (step .gt. 0)) then
      injector_x1 = injector_x1 +&
                        & REAL(injector_reset_interval) * CC * injector_betax
      injector_x2 = injector_x2 -&
                        & REAL(injector_reset_interval) * CC * injector_betax
    end if
    
    ! move the injectors
    old_x1 = injector_x1; old_x2 = injector_x2
    injector_x1 = injector_x1 - injector_betax * CC
    injector_x2 = injector_x2 + injector_betax * CC

    injector_i1_glob = INT(injector_x1)
    injector_i2_glob = INT(injector_x2)

    if (((injector_i1_glob .ge. this_meshblock%ptr%x0) .and.&
       & (injector_i1_glob .lt. this_meshblock%ptr%x0 + this_meshblock%ptr%sx)) .or.&
      & ((injector_i2_glob .ge. this_meshblock%ptr%x0) .and.&
       & (injector_i2_glob .lt. this_meshblock%ptr%x0 + this_meshblock%ptr%sx))) then
      
      if (modulo(step, injector_reset_interval) .eq. 0) then
        ! remove particles left and right from the injectors every once in a while
        do s = 1, nspec
          do ti = 1, species(s)%tile_nx
            do tj = 1, species(s)%tile_ny
              do tk = 1, species(s)%tile_nz
                do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp 
                  x_glob = REAL(species(s)%prtl_tile(ti, tj, tk)%xi(p) + this_meshblock%ptr%x0)&
                         & + species(s)%prtl_tile(ti, tj, tk)%dx(p)
                  if ((x_glob .le. old_x1) .or. (x_glob .gt. old_x2)) then
                    species(s)%prtl_tile(ti, tj, tk)%proc(p) = -1
                  end if
                end do
              end do
            end do
          end do
        end do
      end if
    end if

    ! inject background particles at the injectors' positions
    nUP = 0.5 * ppc0
    
    ! left injector
    back_region%x_min = injector_x1
    back_region%x_max = old_x1
    back_region%y_min = 0.0
    back_region%y_max = REAL(global_mesh%sy)

    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, nUP, upstream_T)

    ! right injector
    back_region%x_min = old_x2
    back_region%x_max = injector_x2
    back_region%y_min = 0.0
    back_region%y_max = REAL(global_mesh%sy)

    call fillRegionWithThermalPlasma(back_region, (/3, 4/), 2, nUP, upstream_T)
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step)
    implicit none
    real                          :: sx_glob, x_glob
    integer                       :: i, j, k
    integer                       :: i_glob, injector_i1_glob, injector_i2_glob
    integer, optional, intent(in) :: step
    
    injector_i1_glob = INT(injector_x1)
    injector_i2_glob = INT(injector_x2)

    if ((injector_i1_glob .le. this_meshblock%ptr%x0 + this_meshblock%ptr%sx) .or.&
      & (injector_i2_glob .ge. this_meshblock%ptr%x0)) then
      ! reset fields left and right from the injectors
      sx_glob = REAL(global_mesh%sx)
      do i = -NGHOST, this_meshblock%ptr%sx - 1 + NGHOST
        i_glob = i + this_meshblock%ptr%x0
        x_glob = REAL(i_glob)
        if ((i_glob .lt. injector_i1_glob) .or. (i_glob .gt. injector_i2_glob)) then
          ex(i, :, :) = 0.0; ey(i, :, :) = 0.0; ez(i, :, :) = 0.0
          bx(i, :, :) = 0.0; bz(i, :, :) = 0.0
          by(i, :, :) = tanh((x_glob - cs_x * sx_glob) / current_width)
        end if
      end do
    end if
  end subroutine userFieldBoundaryConditions
  !............................................................!
end module m_userfile
