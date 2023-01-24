module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
  use m_helpers
  use m_exchangearray, only: exchangeArray
#ifdef USROUTPUT
  use m_writeusroutput
#endif
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real, private :: ph_energy, ph_fraction, gammabeta
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'ph_energy', ph_energy, 0.01)
    call getInput('problem', 'ph_fraction', ph_fraction, -1.0)
    call getInput('problem', 'gammabeta', gammabeta, 10.0)
  end subroutine userReadInput

  function userSLBload(x_glob, y_glob, z_glob, &
                       dummy1, dummy2, dummy3)
    real :: userSLBload
    ! global coordinates
    real, intent(in), optional :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional :: dummy1, dummy2, dummy3
    return
  end function

  function userSpatialDistribution(x_glob, y_glob, z_glob, &
                                   dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3
    return
  end function userSpatialDistribution

  subroutine userInitParticles()
    implicit none
    real :: nUP_elec, nUP_pos, nCS_elec, nCS_pos
    type(region) :: back_region
    real :: sx_glob, sy_glob, shift_gamma, shift_beta, current_sheet_T
    integer :: s, ti, tj, tk, p
    real :: ux, uy, uz, gamma
    integer :: n, ncells, nphotons
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    if (ph_fraction .gt. 0.0) then
      ncells = global_mesh % sx * global_mesh % sy * global_mesh % sz
      nphotons = INT(REAL(ph_fraction, 8) * REAL(ppc0, 8) * REAL(ncells, 8))
      do n = 1, nphotons
        call injectPhoton(0)
      end do
    end if
    species(3) % move_sp = .false.
    species(3) % compton_sp = .true.

    call injectParticleGlobally(1, REAL(global_mesh % sx) * 0.1, &
                                   REAL(global_mesh % sy) * 0.25, &
                                   REAL(global_mesh % sz) * 0.5, &
                                   gammabeta, 0.0, 0.0)
    species(1) % compton_sp = .false.
    species(1) % cool_sp = .true.
    species(1) % deposit_sp = .false.
    call injectParticleGlobally(2, REAL(global_mesh % sx) * 0.1, &
                                   REAL(global_mesh % sy) * 0.75, &
                                   REAL(global_mesh % sz) * 0.5, &
                                   gammabeta, 0.0, 0.0)
    species(2) % compton_sp = .true.
    species(2) % cool_sp = .false.
    species(2) % deposit_sp = .false.
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
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
    real(kind=8) :: injector_x1, injector_x2
    real(kind=8) :: old_x1, old_x2
    real :: x_glob, y_glob
    real :: vpart
    real :: nUP_elec, nUP_pos, nUP_ions
    real :: ux, uy, uz, gamma
    integer :: s, ti, tj, tk, p, nUP_tot
    integer :: injector_i1_glob, injector_i2_glob, old_x1_glob, old_x2_glob
    real :: average_electron
    integer :: ncells, lackofparticles, addedelectron
    type(region) :: back_region
    integer :: i, j, k, i_glob, density_int, n, j_glob
    integer :: injector_i1_up, injector_i2_up
    integer :: ncells_up, lackofparticles_up, average_electron_up
    real :: density_frac, dx, dy, dz
    integer, optional, intent(in) :: step
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()

    integer :: nphotons
    real :: xg, yg, zg, age, kx, ky, kz

    s = 3
    do ti = 1, species(s) % tile_nx
      do tj = 1, species(s) % tile_ny
        do tk = 1, species(s) % tile_nz
          do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
            ! reset scattered photon momentum
            if (species(s) % prtl_tile(ti, tj, tk) % payload2(p) .ge. 0.0) then
              call generateRandomVelocity(kx, ky, kz)
              species(s) % prtl_tile(ti, tj, tk) % payload2(p) = -1.0
              species(s) % prtl_tile(ti, tj, tk) % u(p) = kx
              species(s) % prtl_tile(ti, tj, tk) % v(p) = ky
              species(s) % prtl_tile(ti, tj, tk) % w(p) = kz
            end if
          end do
        end do
      end do
    end do
  end subroutine userParticleBoundaryConditions

  subroutine generateRandomVelocity(kx, ky, kz)
    implicit none
    real, intent(out) :: kx, ky, kz
    real :: rand_costh, rand_phi
    rand_costh = 2.0 * random(dseed) - 1.0
    rand_phi = 2.0 * M_PI * random(dseed)
    kx = ph_energy * sqrt(1.0 - rand_costh**2) * cos(rand_phi)
    ky = ph_energy * sqrt(1.0 - rand_costh**2) * sin(rand_phi)
    kz = ph_energy * rand_costh
  end subroutine generateRandomVelocity

  subroutine injectPhoton(step)
    implicit none
    integer, optional, intent(in) :: step
    real :: xg, yg, zg, kx, ky, kz

    xg = random(dseed) * REAL(global_mesh % sx)
    yg = random(dseed) * REAL(global_mesh % sy)
    zg = 0.5
    call generateRandomVelocity(kx, ky, kz)
    call injectParticleGlobally(3, xg, yg, zg, &
                                kx, ky, kz, &
                                1.0, REAL(step), -1.0, 0.0)
    !                                ----------   ---  ---
    !                                   ^          ^    ^
    !                                   |          |    |
    !                        creation time         |    number of scatterings
    !                                              |
    !                                           last scattering timestep
  end subroutine injectPhoton

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
  end subroutine userFieldBoundaryConditions
  !............................................................!

  subroutine writeUsrRestart(rst_file)
    implicit none
    integer, intent(in) :: rst_file
  end subroutine writeUsrRestart

  subroutine readUsrRestart(rst_file)
    implicit none
    integer, intent(in) :: rst_file
  end subroutine readUsrRestart

  subroutine userDeallocate()
    implicit none
  end subroutine userDeallocate

  elemental subroutine usrSetPhPld(u0, v0, w0, over_e_temp, &
                                   incr_pld1, incr_pld2, incr_pld3)
    !$omp declare simd(usrSetPhPld)
    real, intent(in) :: u0, v0, w0, over_e_temp
    real, intent(out) :: incr_pld1, incr_pld2, incr_pld3
    incr_pld1 = 0.0; incr_pld2 = 0.0; incr_pld3 = 0.0
  end subroutine

  elemental subroutine usrSetElPld(q_over_m, u0, v0, w0, over_e_temp, &
                                   ex0, ey0, ez0, bx0, by0, bz0, &
                                   incr_pld1, incr_pld2, incr_pld3)
    !$omp declare simd(usrSetElPld)
    real, intent(in) :: q_over_m, u0, v0, w0, over_e_temp, &
                        ex0, ey0, ez0, &
                        bx0, by0, bz0
    real, intent(out) :: incr_pld1, incr_pld2, incr_pld3
    incr_pld1 = 0.0; incr_pld2 = 0.0; incr_pld3 = 0.0
  end subroutine

end module m_userfile
