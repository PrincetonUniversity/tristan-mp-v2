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
  real, private :: backgr_T, prtl_beta, pitch_angle, backgr_N, B_Z
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    integer :: h, SEED
    real :: temp

    call getInput('problem', 'SEED', SEED)
    do h = 1, SEED
      temp = random(dseed)
    end do

    call getInput('problem', 'backgr_T', backgr_T)
    call getInput('problem', 'backgr_N', backgr_N)
    call getInput('problem', 'PITCH', pitch_angle, 0.0)
    ! convert pitch angle to radians
    pitch_angle = pitch_angle * M_PI / 180.0
    call getInput('problem', 'prtl_beta', prtl_beta)
    call getInput('problem', 'B_Z', B_Z, 0.0)
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob,&
                                 & dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3
    if (.false.) print *, x_glob, y_glob, z_glob, dummy1, dummy2, dummy3
    userSpatialDistribution = 0.0
    return
  end function

  function userSLBload(x_glob, y_glob, z_glob,&
                     & dummy1, dummy2, dummy3)
    real :: userSLBload
    ! global coordinates
    real, intent(in), optional :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional :: dummy1, dummy2, dummy3
    if (.false.) print *, x_glob, y_glob, z_glob, dummy1, dummy2, dummy3
    userSLBload = 0.0
    return
  end function

  subroutine userInitParticles()
    implicit none
    real :: vx, vy, vz, gamma
    real :: xg, yg, zg
    integer :: npart
    real :: nUP, sx_glob, sy_glob, sz_glob, rho_larmor
    type(region) :: back_region
    integer :: s, ti, tj, tk
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    nUP = 0.5 * ppc0 * backgr_N

    sx_glob = REAL(global_mesh % sx)
    sy_glob = REAL(global_mesh % sy)
    sz_glob = REAL(global_mesh % sz)

    back_region % x_min = 0.0
    back_region % x_max = sx_glob

#if defined(TWO_D) || defined(THREE_D)
    back_region % y_min = 0.0
    back_region % y_max = sy_glob
#endif

#if defined(THREE_D)
    back_region % z_min = 0.0
    back_region % z_max = sz_glob
#endif

    ! initialize all particle velocities to 0
    if (backgr_N .ne. 0) then
      if (backgr_T .eq. 0) then
        call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, nUP, 1e-5)
        do s = 1, 2
          do ti = 1, species(s) % tile_nx
            do tj = 1, species(s) % tile_ny
              do tk = 1, species(s) % tile_nz
                npart = species(s) % prtl_tile(ti, tj, tk) % npart_sp
                species(s) % prtl_tile(ti, tj, tk) % u(1:npart) = 0.0
                species(s) % prtl_tile(ti, tj, tk) % v(1:npart) = 0.0
                species(s) % prtl_tile(ti, tj, tk) % w(1:npart) = 0.0
              end do
            end do
          end do
        end do
      else
        call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, nUP, backgr_T)
      end if
    end if

    vx = prtl_beta * cos(pitch_angle)
    vy = prtl_beta * sin(pitch_angle)
    vz = 0.0
    gamma = 1.0 / sqrt(1.0 - prtl_beta**2)

    !rho_larmor = c_omp * gamma * prtl_beta / sqrt(sigma)

    xg = 0.1 * sx_glob
    !yg = (0.5 * sy_glob) - rho_larmor
    yg = 0.5 * sy_glob
    zg = 0.5 * sz_glob

    vx = gamma * vx
    vy = gamma * vy
    vz = gamma * vz
    call injectParticleGlobally(3, xg, yg, zg, vx, vy, vz)
    call injectParticleGlobally(4, xg, yg, zg, 0.0, 0.0, 0.0)
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    ! integer :: i, j, k
    ! integer :: i_glob, j_glob, k_glob

    ! uncomment when adding external B-field
    ! real    :: pitch_deg, pitch_rad
    ! real, parameter :: PI = 3.1415927
    ! call getInput('problem', 'PITCH', pitch_deg) ! PITCH is in degrees
    ! ! convert to radians
    ! pitch_rad = pitch_deg * PI / 180
    ! ex(:,:,:) = 0; ey(:,:,:) = 0; ez(:,:,:) = 0
    ! bx(:,:,:) = cos(pitch_rad); by(:,:,:) = 0; bz(:,:,:) = sin(pitch_rad)
    ! jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0

    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = 0; by(:, :, :) = 0; bz(:, :, :) = B_Z
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
    if (.false.) print *, step
    ! called after particles move and deposit ...
    ! ... and before the currents are added to the electric field
  end subroutine userCurrentDeposit

  subroutine userDriveParticles(step)
    implicit none
    integer, optional, intent(in) :: step
    if (.false.) print *, step
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
    real, intent(in) :: xp, yp, zp
    real, intent(out) :: ex_ext, ey_ext, ez_ext
    real, intent(out) :: bx_ext, by_ext, bz_ext
    if (.false.) print *, xp, yp, zp
    ! some functions of xp, yp, zp
    ex_ext = 0.0; ey_ext = 0.0; ez_ext = 0.0
    bx_ext = 0.0; by_ext = 0.0; bz_ext = 0.0
  end subroutine userExternalFields
  !............................................................!

  !--- boundaries ---------------------------------------------!
  subroutine userParticleBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
    if (.false.) print *, step
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
    logical :: updateE_, updateB_
    if (.false.) print *, step, updateE, updateB

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
