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
  real    :: upstream_T, amplitude, inject_rate, n_inject, min_inj
  private :: upstream_T, amplitude, inject_rate, n_inject, min_inj
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'upstream_T', upstream_T)
    call getInput('problem', 'amplitude', amplitude)
    call getInput('problem', 'inject_rate', inject_rate)
    call getInput('problem', 'n_inject', n_inject)
    call getInput('problem', 'min_inj', min_inj)
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob,&
                                 & dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional  :: x_glob, y_glob, z_glob
    real, intent(in), optional  :: dummy1, dummy2, dummy3

    return
  end function

  function userSLBload(x_glob, y_glob, z_glob,&
                     & dummy1, dummy2, dummy3)
    real :: userSLBload
    ! global coordinates
    real, intent(in), optional  :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional  :: dummy1, dummy2, dummy3
    return
  end function

  subroutine userInitParticles()
    implicit none
    real                :: nUP
    integer             :: s, ti, tj, tk, p
    real                :: xp, kx, u_
    type(region)        :: back_region
    real                :: sx_glob, shift_gamma, shift_beta, current_sheet_T

    procedure (spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    !nUP = 0.5 * ppc0
    ! back_region%x_min = 0.0
    ! back_region%x_max = global_mesh%sx
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    ex(:,:,:) = amplitude; ey(:,:,:) = 0; ez(:,:,:) = 0
    bx(:,:,:) = 0; by(:,:,:) = 0; bz(:,:,:) = 0
    jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0
  end subroutine userInitFields
  !............................................................!

  !--- driving ------------------------------------------------!
  subroutine userDriveParticles(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userDriveParticles
  !............................................................!

  !--- boundaries ---------------------------------------------!
  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
    ! called after particles move and deposit ...
    ! ... and before the currents are added to the electric field
  end subroutine userCurrentDeposit

  subroutine userParticleBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
    type(region)                  :: back_region

    n_inject = n_inject + 0.5 * inject_rate

    if (n_inject .gt. min_inj) then
      back_region%x_min = 0.0
      back_region%x_max = REAL(global_mesh%sx)
      call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, n_inject, upstream_T)
      n_inject = 0.0
    end if
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

  !--- user-specific output -----------------------------------!
  subroutine userOutput(step)
    implicit none
    integer, optional, intent(in) :: step
    ! ...
  end subroutine userOutput
  !............................................................!
end module m_userfile
