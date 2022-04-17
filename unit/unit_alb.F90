#include "../src/defs.F90"

module m_userfile
  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_loadbalancing
  use m_particlelogistics
  use m_exchangeparts, only: redistributeParticlesBetweenMeshblocks
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real, private :: sph_radius
  real, private :: fly_x, fly_y, fly_z, velocity
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput("problem", "radius", sph_radius, 5.0)
    velocity = 2.0
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob,&
                                 & dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional  :: x_glob, y_glob, z_glob
    real, intent(in), optional  :: dummy1, dummy2, dummy3
    real :: r
    if (present(x_glob) .and. present(y_glob) .and. present(z_glob) .and.&
      & present(dummy1) .and. present(dummy2) .and. present(dummy3)) then
      r = sqrt((x_glob - dummy1)**2 + (y_glob - dummy2)**2 + (z_glob - dummy3)**2)
      if (r .lt. sph_radius) then
        userSpatialDistribution = 1.0
      else
        userSpatialDistribution = 0.0
      end if
    else
      call throwError('ERROR: smth wrong in `userSpatialDistribution()`.')
    end if
    return
  end function

  function userSLBload(x_glob, y_glob, z_glob,&
                     & dummy1, dummy2, dummy3)
    real :: userSLBload
    ! global coordinates
    real, intent(in), optional  :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional  :: dummy1, dummy2, dummy3
    userSLBload = 100.0 / (10.0 + (x_glob)**2)
    return
  end function

  subroutine userInitParticles()
    implicit none
    type(region)    :: back_region
    integer         :: i, inds(2)
    procedure (spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    back_region%x_min = 0.0
    back_region%x_max = REAL(global_mesh%sx)

    #ifdef twoD
      back_region%y_min = 0.0
      back_region%y_max = REAL(global_mesh%sy)
    #endif

    #ifdef threeD
      back_region%y_min = 0.0
      back_region%y_max = REAL(global_mesh%sy)
      back_region%z_min = 0.0
      back_region%z_max = REAL(global_mesh%sz)
    #endif

    inds(1) = 1; inds(2) = 2
    call fillRegionWithThermalPlasma(back_region, inds, 2, ppc0, 0.0, &
                                   & spat_distr_ptr = spat_distr_ptr,&
                                   & dummy1 = 0.5 * REAL(global_mesh%sx),&
                                   & dummy2 = 0.5 * REAL(global_mesh%sy),&
                                   & dummy3 = 0.5 * REAL(global_mesh%sz))
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    ex(:,:,:) = 0; ey(:,:,:) = 0; ez(:,:,:) = 0
    bx(:,:,:) = 0; by(:,:,:) = 0; bz(:,:,:) = 0
    jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0
    ! ... dummy loop ...
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
    integer :: s, ti, tj, tk, p
    if (step .eq. 0) then
      fly_x = 0.5; fly_y = 0.0; fly_z = 0.87
    else if (step .eq. 100) then
      fly_x = 0.87; fly_y = 0.5; fly_z = 0.0
    else if (step .eq. 200) then
      fly_x = 0.25; fly_y = 0.76; fly_z = 0.6
    else if (step .ge. 300) then
      fly_x = 0.0; fly_y = 0.0; fly_z = 0.0
    end if

    if ((step .le. 400) .and. modulo(step, 100) .eq. 0) then
      do s = 1, nspec
        do ti = 1, species(s)%tile_nx
          do tj = 1, species(s)%tile_ny
            do tk = 1, species(s)%tile_nz
              do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
                species(s)%prtl_tile(ti, tj, tk)%u(p) = fly_x * velocity
                species(s)%prtl_tile(ti, tj, tk)%v(p) = fly_y * velocity
                species(s)%prtl_tile(ti, tj, tk)%w(p) = fly_z * velocity
              end do
            end do
          end do
        end do
      end do
    end if
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

  #ifdef GCA
    logical function userEnforceGCA(xi, yi, zi, dx, dy, dz, u, v, w, weight)
      implicit none
      integer(kind=2), intent(in), optional   :: xi, yi, zi
      real, intent(in), optional              :: dx, dy, dz, u, v, w
      real, intent(in), optional              :: weight
      userEnforceGCA = .false.
    end function userEnforceGCA
  #endif
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
    integer, allocatable          :: left_group(:), right_group(:)
    integer                       :: ierr, rnk, r

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
