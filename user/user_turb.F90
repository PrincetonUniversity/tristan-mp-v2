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
  real, private     :: background_T, sigma_z, background_w0
  real, private     :: background_bz, delta_b_rms
  integer, private  :: N_modes
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    real :: gamma_mean
    call getInput('problem', 'background_T', background_T)
    call getInput('problem', 'sigma_z', sigma_z)
    call getInput('problem', 'N_modes', N_modes)
    if (background_T .lt. 0.7) then
      gamma_mean = 1.0 + 1.5 * background_T + 1.875 * background_T**2 -&
                 & 1.875 * background_T**3 + 1.05469 * background_T**4 + 1.40625 * background_T**5
    else
      gamma_mean = 0.5 / background_T + 3.0 * background_T +&
                 & (0.0625 * (-1.23186 - 2.0 * log(background_T))) / background_T**5
    end if
    background_w0 = background_T + gamma_mean

    ! in units of `B_norm`:
    ! ... (delta_Brms / B_norm) = sqrt(w0)
    ! ... (Bz / B_norm) = sqrt(w0 * sigma_z / sigma_0)
    delta_b_rms = sqrt(background_w0)
    background_bz = sqrt(background_w0 * sigma_z / sigma)
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
    real                :: background_n
    type(region)        :: back_region
    procedure (spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    background_n = REAL(ppc0) * 0.5
    back_region%x_min = 0.0
    back_region%y_min = 0.0
    back_region%x_max = REAL(global_mesh%sx)
    back_region%y_max = REAL(global_mesh%sy)
    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, background_n, background_T)
  end subroutine userInitParticles

  real function BxMode(x, y, z, phases_1, phases_2)
    implicit none
    real, intent(in)                :: x, y, z
    real, allocatable, intent(in)   :: phases_1(:,:), phases_2(:,:)
    integer                         :: n, m
    real                            :: dummy, beta_mn, kn, km
    dummy = 0.0
    do n = 1, N_modes
      kn = 2.0 * M_PI * n / REAL(global_mesh%sy)
      do m = 1, N_modes
        km = 2.0 * M_PI * m / REAL(global_mesh%sx)
        beta_mn = 2.0 * delta_b_rms / (N_modes * sqrt(REAL(n**2 + m**2)))
        dummy = dummy + beta_mn * n * sin(km * x + phases_1(m, n)) * cos(kn * y + phases_2(m, n))
      end do
    end do
    BxMode = dummy
    return
  end function

  real function ByMode(x, y, z, phases_1, phases_2)
    implicit none
    real, intent(in)                :: x, y, z
    real, allocatable, intent(in)   :: phases_1(:,:), phases_2(:,:)
    integer                         :: n, m
    real                            :: dummy, beta_mn, kn, km
    dummy = 0.0
    do n = 1, N_modes
      kn = 2.0 * M_PI * n / REAL(global_mesh%sy)
      do m = 1, N_modes
        km = 2.0 * M_PI * m / REAL(global_mesh%sx)
        beta_mn = 2.0 * delta_b_rms / (N_modes * sqrt(REAL(n**2 + m**2)))
        dummy = dummy - beta_mn * n * cos(km * x + phases_1(m, n)) * sin(kn * y + phases_2(m, n))
      end do
    end do
    ByMode = dummy
    return
  end function

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    integer :: root_rnk = 0, n, m, ierr
    real    :: x_, y_, z_
    real, allocatable :: rand_phases_1(:,:), rand_phases_2(:,:)

    allocate(rand_phases_1(N_modes, N_modes), rand_phases_2(N_modes, N_modes))
    if (mpi_rank .eq. root_rnk) then
      do n = 1, N_modes
        do m = 1, N_modes
          rand_phases_1(n, m) = 2 * M_PI * random(dseed)
          rand_phases_2(n, m) = 2 * M_PI * random(dseed)
        end do
      end do
    end if
    call MPI_Bcast(rand_phases_1, N_modes**2, MPI_REAL, root_rnk, MPI_COMM_WORLD, i)
    call MPI_Bcast(rand_phases_2, N_modes**2, MPI_REAL, root_rnk, MPI_COMM_WORLD, i)

    ex(:,:,:) = 0; ey(:,:,:) = 0; ez(:,:,:) = 0
    bx(:,:,:) = 0; by(:,:,:) = 0; bz(:,:,:) = 0
    jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0

    bz(:,:,:) = background_bz

    ! ... dummy loop ...
    do i = 0, this_meshblock%ptr%sx - 1
      i_glob = i + this_meshblock%ptr%x0
      do j = 0, this_meshblock%ptr%sy - 1
        j_glob = j + this_meshblock%ptr%y0
        do k = 0, this_meshblock%ptr%sz - 1
          k_glob = k + this_meshblock%ptr%z0

          x_ = REAL(i_glob);  y_ = REAL(j_glob) + 0.5;  z_ = REAL(k_glob) + 0.5
          bx(i, j, k) = BxMode(x_, y_, z_, rand_phases_1, rand_phases_2)

          x_ = REAL(i_glob) + 0.5;  y_ = REAL(j_glob);  z_ = REAL(k_glob) + 0.5
          by(i, j, k) = ByMode(x_, y_, z_, rand_phases_1, rand_phases_2)
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
