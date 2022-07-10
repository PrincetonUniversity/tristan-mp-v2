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
  integer :: nwaves
  real :: upstream_T, amplitude
  private :: nwaves, upstream_T, amplitude
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'upstream_T', upstream_T)
    call getInput('problem', 'nwaves', nwaves)
    call getInput('problem', 'amplitude', amplitude)
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
    real :: nUP
    integer :: s, ti, tj, tk, p
    real :: xp, yp, zp, kx, u_, v_, w_
    type(region) :: back_region
    real :: sx_glob, sy_glob, shift_gamma, shift_beta, current_sheet_T
    type(maxwellian) :: maxw

    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    back_region % x_min = 0.0
    back_region % x_max = global_mesh % sx
#ifdef twoD
    back_region % y_min = 0.0
    back_region % y_max = global_mesh % sy
#endif

    nUP = ppc0
    kx = 2.0 * M_PI * REAL(nwaves) / REAL(global_mesh % sx)
    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, nUP, upstream_T)
    s = 1
    do ti = 1, species(s) % tile_nx
      do tj = 1, species(s) % tile_ny
        do tk = 1, species(s) % tile_nz
          do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
            xp = REAL(species(s) % prtl_tile(ti, tj, tk) % xi(p)) + species(s) % prtl_tile(ti, tj, tk) % dx(p)
            xp = xp + REAL(this_meshblock % ptr % x0)
            u_ = amplitude * sin(kx * xp)
            species(s) % prtl_tile(ti, tj, tk) % u(p) = species(s) % prtl_tile(ti, tj, tk) % u(p) + u_
          end do
        end do
      end do
    end do
    ! "removing" ions
    s = 2
    do ti = 1, species(s) % tile_nx
      do tj = 1, species(s) % tile_ny
        do tk = 1, species(s) % tile_nz
          species(s) % prtl_tile(ti, tj, tk) % npart_sp = 0
        end do
      end do
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
