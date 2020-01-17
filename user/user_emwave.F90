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
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob,&
                                 & dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional  :: x_glob, y_glob, z_glob
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
    real :: kx, ky, ex_norm, ey_norm, exy_norm

    ex(:,:,:) = -1; ey(:,:,:) = -1; ez(:,:,:) = -1
    bx(:,:,:) = -1; by(:,:,:) = -1; bz(:,:,:) = -1

    kx = 5; ky = 2
    kx = kx * 2 * M_PI / global_mesh%sx
    ky = ky * 2 * M_PI / global_mesh%sy
    if (ky .ne. 0) then
      ex_norm = 1; ey_norm = (-kx / ky)
    else
      ey_norm = 1; ex_norm = (-ky / kx)
    end if
    exy_norm = sqrt(ex_norm**2 + ey_norm**2)
    ex_norm = ex_norm / exy_norm
    ey_norm = ey_norm / exy_norm

    do i = 0, this_meshblock%ptr%sx - 1
      i_glob = i + this_meshblock%ptr%x0
      do j = 0, this_meshblock%ptr%sy - 1
        j_glob = j + this_meshblock%ptr%y0
        do k = 0, this_meshblock%ptr%sz - 1
          k_glob = k + this_meshblock%ptr%z0
          ex(i, j, k) = ex_norm * sin((i_glob) * kx + (j_glob - 0.5) * ky)
          ey(i, j, k) = ey_norm * sin((i_glob - 0.5) * kx + (j_glob) * ky)
          ez(i, j, k) = 0
          bx(i, j, k) = 0
          by(i, j, k) = 0
          bz(i, j, k) = sin((i_glob) * kx + (j_glob) * ky)
        end do
      end do
    end do
  end subroutine userInitFields
  !............................................................!

  !--- driving ------------------------------------------------!
  subroutine userDriveParticles()
    implicit none
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
  subroutine userParticleBoundaryConditions()
    implicit none
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions()
    implicit none
  end subroutine userFieldBoundaryConditions
  !............................................................!
end module m_userfile
