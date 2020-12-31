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
  real :: upstream_T, current_width, cs_x1, cs_x2, nCS_over_nUP
  private :: upstream_T, current_width, cs_x1, cs_x2, nCS_over_nUP
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'upstream_T', upstream_T)      ! <- upstream temperature
    call getInput('problem', 'nCS_nUP', nCS_over_nUP)       ! <- density of the CS / density upstream
    call getInput('problem', 'current_width', current_width)
    cs_x1 = 0.25; cs_x2 = 0.75
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob,&
                                 & dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional  :: x_glob, y_glob, z_glob
    real, intent(in), optional  :: dummy1, dummy2, dummy3
    ! dummy1 -> position of the current sheet
    ! dummy2 -> width of the current sheet
    userSpatialDistribution = 1.0 / (cosh((x_glob - dummy1) / dummy2))**2
    return
  end function

  subroutine userInitParticles()
    implicit none
    real            :: nUP, nCS
    type(region)    :: fill_region
    real            :: sx_glob, sy_glob, shift_gamma, shift_beta, current_sheet_T
    procedure (spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    nUP = 0.5 * ppc0
    nCS = nUP * nCS_over_nUP

    ! initial drift velocity of the current sheet
    !     to ensure `curlB = j`
    shift_beta = sqrt(sigma) * c_omp / (current_width * nCS_over_nUP)
    if (shift_beta .gt. 1) then
      call throwError('ERROR: `shift_beta` > 1 in `userInitParticles()`')
    end if
    shift_gamma = 1.0 / sqrt(1.0 - shift_beta**2)
    !     to ensure `B^2/2 = sigma T n me c^2`
    current_sheet_T = 0.5 * sigma / nCS_over_nUP

    sx_glob = REAL(global_mesh%sx)
    sy_glob = REAL(global_mesh%sy)

    fill_region%x_min = 0.0
    fill_region%y_min = 0.0
    fill_region%x_max = sx_glob
    fill_region%y_max = sy_glob

    ! now filling region with plasma
    !     upstream:
    call fillRegionWithThermalPlasma(fill_region, (/1, 2/), 2, nUP, upstream_T)
    !     current sheet #1:
    call fillRegionWithThermalPlasma(fill_region, (/1, 2/), 2, nCS, current_sheet_T,&
                                   & shift_gamma = shift_gamma, shift_dir = 3,&
                                   & spat_distr_ptr = spat_distr_ptr,&
                                   & dummy1 = cs_x1 * sx_glob, dummy2 = current_width)
    !     current sheet #2:
    call fillRegionWithThermalPlasma(fill_region, (/1, 2/), 2, nCS, current_sheet_T,&
                                   & shift_gamma = shift_gamma, shift_dir = -3,&
                                   & spat_distr_ptr = spat_distr_ptr,&
                                   & dummy1 = cs_x2 * sx_glob, dummy2 = current_width)

  end subroutine userInitParticles

subroutine userInitFields()
  implicit none
  integer :: i, j, k
  real    :: sx_glob, x_glob
  ex(:,:,:) = 0; ey(:,:,:) = 0; ez(:,:,:) = 0
  bx(:,:,:) = 0; by(:,:,:) = 0; bz(:,:,:) = 0
  jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0

  sx_glob = REAL(global_mesh%sx)
  do i = 0, this_meshblock%ptr%sx - 1
    x_glob = i + this_meshblock%ptr%x0 + 0.5
    do j = 0, this_meshblock%ptr%sy - 1
      do k = 0, this_meshblock%ptr%sz - 1
        by(i, j, k) = tanh((x_glob - cs_x1 * sx_glob) / current_width) -&
                    & tanh((x_glob - cs_x2 * sx_glob) / current_width) - 1
      end do
    end do
  end do
end subroutine userInitFields
  !............................................................!

  !--- driving ------------------------------------------------!
  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
    ! called after particles move and deposit ...
    ! ... and before the currents are added to the electric field
  end subroutine userCurrentDeposit

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
  !............................................................!

  !--- boundaries ---------------------------------------------!
  subroutine userParticleBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
    logical                       :: updateE_, updateB_

    if (present(updateE)) then
      updateE_ = updateE
    else
      updateE_ = .true.
    end if

    if (present(updateB)) then
      updateB_ = updateB
    else
      updateB_ = .true.
    end if
  end subroutine userFieldBoundaryConditions
  !............................................................!
end module m_userfile
