! Configuration for this userfile:
! ```
!   $ python configure.py --nghosts=5 --user=user_two_bulbs -qed -bwpp
! ```

module m_userfile
  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real :: ph_energy_1, ph_energy_2

  private :: ph_energy_1, ph_energy_2
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'ph_e1', ph_energy_1)
    call getInput('problem', 'ph_e2', ph_energy_2)
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob, &
                                   dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3

    return
  end function

  function userSLBload(x_glob, y_glob, z_glob, &
                       dummy1, dummy2, dummy3)
    real :: userSLBload
    ! global coordinates
    real, intent(in), optional :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional :: dummy1, dummy2, dummy3
    return
  end function

  subroutine userInitParticles()
    implicit none
    real :: xg, yg, zg, kx, ky, kz, U_, TH_
    integer :: ntot, n
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    ntot = (global_mesh % sx * global_mesh % sy * global_mesh % sz * ppc0) / 2

    do n = 1, ntot
      xg = random(dseed) * (global_mesh % sx)
      yg = random(dseed) * (global_mesh % sy)
      zg = random(dseed) * (global_mesh % sz)

      U_ = 2 * (random(dseed) - 0.5)
      TH_ = 2 * M_PI * random(dseed)
      kx = ph_energy_1 * sqrt(1 - U_**2) * cos(TH_)
      ky = ph_energy_1 * sqrt(1 - U_**2) * sin(TH_)
      kz = ph_energy_1 * U_
      call injectParticleGlobally(1, xg, yg, zg, kx, ky, kz)
    end do

    do n = 1, ntot
      xg = random(dseed) * (global_mesh % sx)
      yg = random(dseed) * (global_mesh % sy)
      zg = random(dseed) * (global_mesh % sz)

      U_ = 2 * (random(dseed) - 0.5)
      TH_ = 2 * M_PI * random(dseed)
      kx = ph_energy_2 * sqrt(1 - U_**2) * cos(TH_)
      ky = ph_energy_2 * sqrt(1 - U_**2) * sin(TH_)
      kz = ph_energy_2 * U_
      call injectParticleGlobally(2, xg, yg, zg, kx, ky, kz)
    end do
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = 0; by(:, :, :) = 0; bz(:, :, :) = 0
    jx(:, :, :) = 0; jy(:, :, :) = 0; jz(:, :, :) = 0
    ! ... dummy loop ...
    ! do i = 0, this_meshblock%ptr%sx - 1
    !   i_glob = i + this_meshblock%ptr%x0
    !   do j = 0, this_meshblock%ptr%sy - 1
    !     j_glob = j + this_meshblock%ptr%y0
    !     do k = 0, this_meshblock%ptr%sz - 1
    !       k_glob = k + this_meshblock%ptr%z0
    !       ...
    !     end do
    !   end do
    ! end do
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

  subroutine userExternalFields(xp, yp, zp, &
                                ex_ext, ey_ext, ez_ext, &
                                bx_ext, by_ext, bz_ext)
    implicit none
    real, intent(in) :: xp, yp, zp
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
    integer, optional, intent(in) :: step
    integer :: i
    real :: thet, rnd, u_, v_, w_
    real :: x1_g, y1_g, x2_g, y2_g
    real :: x1_l, y1_l, x2_l, y2_l
    real :: dx_, dy_, dz_, x_, y_
    integer(kind=2) :: xi_, yi_, zi_

  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
    logical :: updateE_, updateB_

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

#include "optional.F"

end module m_userfile
