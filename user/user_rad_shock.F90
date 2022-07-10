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
  real :: gamma_up, dgamma_up
  private :: gamma_up, dgamma_up
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'gamma', gamma_up)
    call getInput('problem', 'dgamma', dgamma_up)
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
    real :: dens
    type(region) :: back_region
    real :: eph, xg, yg, zg, kx, ky, kz, U_, TH_
    integer :: ntot, n
    real :: p_dot_k, px
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    dens = 0.5 * ppc0
#if defined(twoD) || defined (threeD)
    back_region % y_min = 0.0
    back_region % y_max = REAL(global_mesh % sy)
#endif
    back_region % x_min = 0.0
    back_region % x_max = 0.5 * REAL(global_mesh % sx)
    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, dens, dgamma_up, &
                                     shift_gamma=gamma_up, shift_dir=1, zero_current=.true.)
    back_region % x_min = 0.5 * REAL(global_mesh % sx)
    back_region % x_max = REAL(global_mesh % sx)
    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, dens, dgamma_up, &
                                     shift_gamma=gamma_up, shift_dir=-1, zero_current=.true.)

    ntot = global_mesh % sx * global_mesh % sy * global_mesh % sz * ppc0 * 0.5
    do n = 1, ntot
      xg = random(dseed) * (global_mesh % sx) * 0.5
      yg = random(dseed) * (global_mesh % sy)
      zg = 0.5
      U_ = 2 * (random(dseed) - 0.5)
      TH_ = 2 * M_PI * random(dseed)
      eph = dgamma_up
      kx = eph * sqrt(1 - U_**2) * cos(TH_)
      ky = eph * sqrt(1 - U_**2) * sin(TH_)
      kz = eph * U_
      px = -sqrt(gamma_up**2 - 1.0)
      p_dot_k = kx * px
      kx = kx + (p_dot_k / (gamma_up + 1.0) - eph) * px
      call injectParticleGlobally(3, xg, yg, zg, kx, ky, kz)
    end do
    do n = 1, ntot
      xg = (random(dseed) + 1.0) * (global_mesh % sx) * 0.5
      yg = random(dseed) * (global_mesh % sy)
      zg = 0.5
      U_ = 2 * (random(dseed) - 0.5)
      TH_ = 2 * M_PI * random(dseed)
      eph = dgamma_up
      kx = eph * sqrt(1 - U_**2) * cos(TH_)
      ky = eph * sqrt(1 - U_**2) * sin(TH_)
      kz = eph * U_
      px = sqrt(gamma_up**2 - 1.0)
      p_dot_k = kx * px
      kx = kx + (p_dot_k / (gamma_up + 1.0) - eph) * px
      call injectParticleGlobally(3, xg, yg, zg, kx, ky, kz)
    end do
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = 0; by(:, :, :) = 0; bz(:, :, :) = 0
    jx(:, :, :) = 0; jy(:, :, :) = 0; jz(:, :, :) = 0
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
end module m_userfile
