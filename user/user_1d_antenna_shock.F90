module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
  use m_exchangearray
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real :: antenna_n
  real :: antenna_A
  real :: background_T
  real :: antenna_omega
  real :: flow_gamma

  private :: antenna_n
  private :: antenna_A
  private :: background_T
  private :: antenna_omega
  private :: flow_gamma
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'antenna_n', antenna_n)
    call getInput('problem', 'antenna_A', antenna_A)
    call getInput('problem', 'antenna_omega', antenna_omega)
    call getInput('problem', 'background_T', background_T)
    call getInput('problem', 'flow_gamma', flow_gamma)
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob,&
                                 & dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3

    return
  end function userSpatialDistribution

  function userSLBload(x_glob, y_glob, z_glob,&
                     & dummy1, dummy2, dummy3)
    real :: userSLBload
    ! global coordinates
    real, intent(in), optional :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional :: dummy1, dummy2, dummy3

    return
  end function

  subroutine userInitParticles()
    implicit none
    integer :: s, i, j, k, ti, tj, tk, p, x_iglob, ierr, density_int, n, x_int, y_int, z_int
    real :: x_glob, y_glob, sy_glob, u_, v_, w_, antenna_k, dv, dgamma, du, beta
    real, allocatable :: densavg_loc(:), densavg_glob(:)
    real :: local_sigma, local_EB, b_sq, weight_inject, nUP, pulsator, density_frac, dx, dy, dz, npairs, nelectron, npositron
    logical :: dummy_flag
    type(region) :: back_region
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()

    type(maxwellian) :: mwl
    mwl % dimension = 1
    mwl % generated = .false.
    mwl % temperature = background_T
    mwl % shift_gamma = flow_gamma
    mwl % shift_flag = .true.
    mwl % shift_dir = -1

    antenna_k = 8.0 * antenna_n * atan(1.0) / REAL(global_mesh % sx)

    call computeDensity(1, reset=.true., ds=0)
    call computeDensity(2, reset=.false., ds=0)
    call exchangeArray()

    ! inject particles in the domain continuously when below threshold
    do i = 0, this_meshblock % ptr % sx - 1
      do j = 0, this_meshblock % ptr % sy - 1
        do k = 0, this_meshblock % ptr % sz - 1
          x_glob = REAL(i + this_meshblock % ptr % x0)
          y_glob = REAL(j + this_meshblock % ptr % y0)

          npairs = 0.5 * ppc0
          nelectron = 0.0
          npositron = 0.0

          density_int = INT(npairs - 0.5 * lg_arr(i, j, k))
          density_frac = (npairs - 0.5 * lg_arr(i, j, k)) - REAL(density_int)

          ! inject integer amount of particles
          do n = 1, density_int
            dx = random(dseed); dy = 0.5
            dz = 0.5
            u_ = 0.0; v_ = 0.0; w_ = 0.0; 
            call generateFromMaxwellian(mwl, u_, v_, w_)
            call createParticle(1, INT(i, 2), INT(j, 2), INT(k, 2), dx, dy, dz, u_, v_, w_)
            call createParticle(2, INT(i, 2), INT(j, 2), INT(k, 2), dx, dy, dz, u_, v_, w_)
          end do

          ! inject float amount of particles (statistical)
          if (random(dseed) .lt. density_frac) then
            dx = random(dseed); dy = 0.5
            dz = 0.5
            u_ = 0.0; v_ = 0.0; w_ = 0.0; 
            call generateFromMaxwellian(mwl, u_, v_, w_)
            call createParticle(1, INT(i, 2), INT(j, 2), INT(k, 2), dx, dy, dz, u_, v_, w_)
            call createParticle(2, INT(i, 2), INT(j, 2), INT(k, 2), dx, dy, dz, u_, v_, w_)
          end if

        end do
      end do
    end do

  end subroutine userInitParticles

  !--- driving ------------------------------------------------!
  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
    ! integer :: i, j, i_glob, j_glob
    ! real :: x_glob, y_glob
    ! ! called after particles move and deposit ...
    ! ! ... and before the currents are added to the electric field

  end subroutine userCurrentDeposit

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob
    real :: x_glob, y_glob, antenna_k, beta_drift

    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = 0; by(:, :, :) = 0; bz(:, :, :) = 1
    jx(:, :, :) = 0; jy(:, :, :) = 0; jz(:, :, :) = 0

    antenna_k = 8.0 * antenna_n * atan(1.0) / REAL(global_mesh % sx)
    beta_drift = sqrt(1.0 - 1.0 / flow_gamma**2)

    do i = -NGHOST, this_meshblock % ptr % sx - 1 + NGHOST
      i_glob = i + this_meshblock % ptr % x0
      x_glob = REAL(i_glob)

      ey(i, :, :) = (-beta_drift) * bz(i, :, :)
      ez(i, :, :) = -(-beta_drift) * by(i, :, :)

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
                              & bx_ext, by_ext, bz_ext, step)
    implicit none
    integer, intent(in) :: step
    real, intent(in) :: xp, yp, zp
    real, intent(out) :: ex_ext, ey_ext, ez_ext
    real, intent(out) :: bx_ext, by_ext, bz_ext
    real :: antenna_k, x_glob, y_glob
    ! some functions of xp, yp, zp
    ex_ext = 0.0; ey_ext = 0.0; ez_ext = 0.0
    bx_ext = 0.0; by_ext = 0.0; bz_ext = 0.0

  end subroutine userExternalFields
  !............................................................!

  function userBoundX(step, i)
    integer, intent(in) :: step
    integer, intent(in) :: i
    integer :: i_glob
    real :: userBoundX
    userBoundX = 0.0
    return
  end function userBoundX

  function userBoundY(step, i)
    integer, intent(in) :: step
    integer, intent(in) :: i
    integer :: i_glob
    real :: userBoundY
    userBoundY = 0.0
    return
  end function userBoundY

  !--- boundaries ---------------------------------------------!
  subroutine userParticleBoundaryConditions(step)
    implicit none
    integer :: s, i, j, k, ti, tj, tk, p, x_iglob, ierr, density_int, n
    real :: x_glob, y_glob, sy_glob, u_, v_, w_, t0, tm
    real, allocatable :: densavg_loc(:), densavg_glob(:)
    real :: local_sigma, local_EB, b_sq, weight_inject, nUP, pulsator, density_frac, dx, dy, dz
    logical :: dummy_flag
    real, dimension(0:this_meshblock % ptr % sx - 1, 0:this_meshblock % ptr % sy - 1) :: densitymax
    type(region) :: back_region
    integer, optional, intent(in) :: step
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()

#ifdef MPI08
    type(MPI_REQUEST) :: mpi_req
#endif

#ifdef MPI
    integer :: mpi_req
#endif

  end subroutine userParticleBoundaryConditions

  subroutine userGridManagement(step)
    implicit none
    integer, optional, intent(in) :: step
    return
  end subroutine userGridManagement

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
#include "optional.F"
  !............................................................!
end module m_userfile
