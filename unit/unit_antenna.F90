#include "../src/defs.F90"

module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_particlelogistics
  implicit none

  !--- PRIVATE functions -----------------------------------------!
  private :: userInitParticles, userInitFields
  !...............................................................!
contains
  subroutine userInitialize()
    implicit none
    integer        :: npart
    npart = INT(this_meshblock%ptr%sx * this_meshblock%ptr%sy * this_meshblock%ptr%sz * ppc0)
    call userInitParticles(npart)
    call userInitFields()
  end subroutine userInitialize

  subroutine userDriveParticles()
    implicit none
    ! integer :: s, ti, tj, tk, p
    ! do s = 1, nspec
    !   do ti = 1, species(s)%tile_nx
    !     do tj = 1, species(s)%tile_ny
    !       do tk = 1, species(s)%tile_nz
    !         do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
    !           species(s)%prtl_tile(ti, tj, tk)%u(p) = 1.0
    !         end do
    !       end do
    !     end do
    !   end do
    ! end do
  end subroutine userDriveParticles

  subroutine userInitParticles(npart)
    implicit none
    integer, intent(in) :: npart
    integer(kind=2)     :: xi_, yi_, zi_
    integer(kind=2)     :: i, j
    integer             :: p, s
    real                :: dx_, dy_, dz_, u_, v_

    #ifndef threeD
      dz_ = 0.5; zi_ = 0
    #else
      dz_ = 10; zi_ = 0.43
    #endif

    if (mpi_rank .eq. 0) then
      i = 10; j = 10
      call createParticle(1, i, j, zi_, 0.5, 0.5, dz_, 0.1, 0.0, 0.0)
      call createParticle(2, i, j, zi_, 0.5, 0.5, dz_, -0.1, 0.0, 0.0)
    endif
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: ind1, ind2, ind3
    integer :: i_glob, j_glob, k_glob
    ex(:,:,:) = 0; ey(:,:,:) = 0; ez(:,:,:) = 0
    bx(:,:,:) = 0; by(:,:,:) = 0; bz(:,:,:) = 0
  end subroutine userInitFields
end module m_userfile
