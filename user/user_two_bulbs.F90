#include "../src/defs.F90"

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

  procedure (spatialDistribution), pointer :: user_slb_load_ptr => null()

  !--- PRIVATE variables -----------------------------------------!
  integer   :: ph_ndot, inject_interval
  real      :: ph_energy, del_x1, del_x2

  private   :: ph_ndot, ph_energy, del_x1, del_x2
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
    user_slb_load_ptr => userSLBload
  end subroutine userInitialize

  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'ndot', ph_ndot)
    call getInput('problem', 'energy', ph_energy)
    call getInput('problem', 'dx1', del_x1)
    call getInput('problem', 'dx2', del_x2)
    call getInput('problem', 'inj_interval', inject_interval, 100000)
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
    procedure (spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    ex(:,:,:) = 0; ey(:,:,:) = 0; ez(:,:,:) = 0
    bx(:,:,:) = 0; by(:,:,:) = 0; bz(:,:,:) = 0
    jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0
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
    integer                       :: i
    real                          :: thet, rnd, u_, v_, w_
    real                          :: x1_g, y1_g, x2_g, y2_g
    real                          :: x1_l, y1_l, x2_l, y2_l
    real                          :: dx_, dy_, dz_, x_, y_
    integer(kind=2)               :: xi_, yi_, zi_
    
    if (step .lt. inject_interval) then
      x1_g = global_mesh%sx * del_x1; y1_g = global_mesh%sy * 0.5
      x2_g = global_mesh%sx * del_x2; y2_g = global_mesh%sy * 0.5

      dz_ = 0.5; zi_ = 0

      do i = 1, ph_ndot
        rnd = random(dseed)
        thet = M_PI * (rnd - 0.5)
        u_ = cos(thet) * ph_energy
        v_ = sin(thet) * ph_energy
        w_ = 0.0

        x_ = x1_g + 10 * random(dseed) * cos(thet)
        y_ = y1_g + 10 * random(dseed) * sin(thet)

        call injectParticleGlobally(1, x_, y_, 0.5, u_, v_, w_)

        rnd = random(dseed)
        thet = M_PI * (rnd - 0.5)
        u_ = -cos(thet) * ph_energy
        v_ = sin(thet) * ph_energy
        w_ = 0.0

        x_ = x2_g - 10 * random(dseed) * cos(thet)
        y_ = y2_g + 10 * random(dseed) * sin(thet)

        call injectParticleGlobally(2, x_, y_, 0.5, u_, v_, w_)
      end do
    end if
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userFieldBoundaryConditions
  !............................................................!
end module m_userfile
