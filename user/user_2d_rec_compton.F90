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
  real, private :: nCS_over_nUP, current_width, upstream_T, cs_x
  real, private :: boost_y_Gamma, boost_y_beta
  real, private :: injector_sx, measure_x
  real(kind=8), private :: injector_x1_fld, injector_x2_fld
  real, private :: fraction_ions
  integer, private :: injector_reset_interval, open_boundaries, no_cooling
  integer, private :: cs_lecs, cs_ions, cs_heavy, up_lecs, up_ions, up_heavy, ph_index
  logical, private :: perturb, simple_bc, fake_photons, use_planckian
  real, private :: ph_energy, ph_fraction, ph_maxage, ph_temperature
  real, private :: bguide
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution, generateRandomDirection, generateRandomMomentum
  private :: injectPhoton, planckSample
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'upstream_T', upstream_T)
    call getInput('problem', 'current_width', current_width)
    call getInput('problem', 'nCS_nUP', nCS_over_nUP, 0.0)
    if (boundary_x .ne. 1) then
      call getInput('problem', 'injector_sx', injector_sx)
    end if
    call getInput('problem', 'cs_lecs', cs_lecs, 1)
    call getInput('problem', 'cs_ions', cs_ions, 2)
    call getInput('problem', 'up_lecs', up_lecs, 3)
    call getInput('problem', 'up_ions', up_ions, 4)
    call getInput('problem', 'measure_x', measure_x, 0.2)
    call getInput('problem', 'open_boundaries', open_boundaries, -1)
    call getInput('problem', 'perturb', perturb, .true.)
    call getInput('problem', 'simple_bc', simple_bc, .false.)
    call getInput('problem', 'no_cooling', no_cooling, -1)

    call getInput('problem', 'bguide', bguide, 0.0)
    call getInput('problem', 'fake_photons', fake_photons, .false.)

    call getInput('problem', 'ph_index', ph_index, 5)
    ! energy for the monoenergetic photon distribution (delta function)
    call getInput('problem', 'ph_energy', ph_energy, -1.0)
    ! temperature for the Planckian photon distribution
    call getInput('problem', 'ph_temperature', ph_temperature, -1.0)
    call getInput('problem', 'ph_fraction', ph_fraction, -1.0)

    if ((ph_energy .lt. 0.0) .and. (ph_temperature .lt. 0.0)) then
      call throwError("ERROR: either `ph_energy` or `ph_temperature` have to be > 0.")
    end if
    if ((ph_energy .gt. 0.0) .and. (ph_temperature .gt. 0.0)) then
      if (mpi_rank .eq. 0) then
        print *, "WARNING: both `ph_energy` and `ph_temperature` are > 0. Picking Planckian distribution and ignoring `ph_energy`."
      end if
    end if
    use_planckian = (ph_temperature .gt. 0.0)

    call getInput('problem', 'ph_maxage', ph_maxage, 10000.0)
    if (current_width .lt. 0.0) then
      call throwError("ERROR: `current_width` has to be > 0.")
    end if
    cs_x = 0.5
  end subroutine userReadInput

  function userSLBload(x_glob, y_glob, z_glob, &
                       dummy1, dummy2, dummy3)
    real :: userSLBload
    ! global coordinates
    real, intent(in), optional :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional :: dummy1, dummy2, dummy3
    real :: distance, sx_glob

    sx_glob = REAL(global_mesh % sx)
    distance = (x_glob - 0.5 * sx_glob)
    userSLBload = 10 * exp(-distance**2 / (20 * current_width)**2) + 1
    return
  end function

  function userSpatialDistribution(x_glob, y_glob, z_glob, &
                                   dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3
    real :: rad2
    if (present(x_glob) .and. present(y_glob) .and. &
        present(dummy1) .and. present(dummy2)) then
      if (present(dummy3) .and. (dummy3 .ne. 0.0)) then
        ! underdensity at the center
        rad2 = (x_glob - dummy1)**2 + (y_glob - dummy3)**2
        userSpatialDistribution = 1.0 / (cosh((x_glob - dummy1) / dummy2))**2 * &
                                  (1.0 - exp(-rad2 / (5.0 * dummy2)**2))
      else
        ! uniform current layer
        userSpatialDistribution = 1.0 / (cosh((x_glob - dummy1) / dummy2))**2
      end if
    else
      call throwError("ERROR: variable not present in `userSpatialDistribution()`")
    end if
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

    nUP_elec = 0.5 * ppc0
    nCS_elec = nUP_elec * nCS_over_nUP

    nUP_pos = 0.5 * ppc0
    nCS_pos = nUP_pos * nCS_over_nUP

    sx_glob = REAL(global_mesh % sx)
    sy_glob = REAL(global_mesh % sy)

    back_region % x_min = 0.0
    back_region % y_min = 0.0
    back_region % x_max = sx_glob
    back_region % y_max = sy_glob
    call fillRegionWithThermalPlasma(back_region, (/up_lecs, up_ions/), 2, nUP_pos, upstream_T)

    if (nCS_over_nUP .ne. 0) then
      shift_beta = sqrt(sigma) * c_omp / (current_width * nCS_over_nUP)
      if (shift_beta .ge. 1) then
        call throwError('ERROR: `shift_beta` >= 1 in `userInitParticles()`')
      end if
      shift_gamma = 1.0 / sqrt(1.0 - shift_beta**2)
      current_sheet_T = 0.5 * sigma * (1.0 + bguide**2) / nCS_over_nUP

      back_region % x_min = sx_glob * cs_x - 10 * current_width
      back_region % x_max = sx_glob * cs_x + 10 * current_width
      back_region % y_min = 0
      back_region % y_max = sy_glob
      if (perturb) then
        call fillRegionWithThermalPlasma(back_region, (/cs_lecs, cs_ions/), 2, nCS_pos, current_sheet_T, &
                                         shift_gamma=shift_gamma, shift_dir=3, &
                                         spat_distr_ptr=spat_distr_ptr, &
                                         dummy1=cs_x * sx_glob, dummy2=current_width, dummy3=cs_x * sy_glob)
      else
        call fillRegionWithThermalPlasma(back_region, (/cs_lecs, cs_ions/), 2, nCS_pos, current_sheet_T, &
                                         shift_gamma=shift_gamma, shift_dir=3, &
                                         spat_distr_ptr=spat_distr_ptr, &
                                         dummy1=cs_x * sx_glob, dummy2=current_width)
      end if
    end if

    if (fake_photons) then
      if (ph_fraction .gt. 0.0) then
        ncells = global_mesh % sx * global_mesh % sy * global_mesh % sz
        nphotons = INT(REAL(ph_fraction, 8) * REAL(ppc0, 8) * REAL(ncells, 8))
        do n = 1, nphotons
          call injectPhoton(0)
        end do
        species(ph_index) % move_sp = .false.
      end if
    end if

  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob
    real :: x_glob, sx_glob
    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = 0; by(:, :, :) = 0; bz(:, :, :) = 0
    jx(:, :, :) = 0; jy(:, :, :) = 0; jz(:, :, :) = 0

    k = 0
    sx_glob = REAL(global_mesh % sx)
    do i = -NGHOST, this_meshblock % ptr % sx - 1 + NGHOST
      i_glob = i + this_meshblock % ptr % x0
      x_glob = REAL(i_glob) + 0.5
      by(i, :, :) = tanh((x_glob - cs_x * sx_glob) / current_width)
    end do
    bz(:, :, :) = bguide
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

#ifdef RADIATION
    species(cs_lecs) % cool_sp = .false.
    species(cs_ions) % cool_sp = .false.
    if (step .lt. no_cooling) then
      species(up_lecs) % cool_sp = .false.
      species(up_ions) % cool_sp = .false.
    else
      species(up_lecs) % cool_sp = .true.
      species(up_ions) % cool_sp = .true.
    end if
#endif

    !
    ! Moving injector boundary conditions
    !
    if (boundary_x .ne. 1) then
      old_x1 = REAL(injector_sx, 8) - 1.0d-7
      old_x2 = REAL(global_mesh % sx, 8) - REAL(injector_sx, 8) + 1.0d-7

      ! inject background particles at the injectors' positions
      nUP_elec = 0.5 * ppc0
      nUP_pos = 0.5 * ppc0

      call computeDensity(up_lecs, reset=.true., ds=0, charge=.false.)

      old_x1_glob = INT(old_x1) + 1
      old_x2_glob = INT(old_x2)
      injector_i1_glob = old_x1_glob + 1
      injector_i2_glob = old_x2_glob - 1
      injector_i1_up = old_x1_glob + 10
      injector_i2_up = old_x2_glob - 10

      ncells = 0
      average_electron = 0

      average_electron_up = 0
      ncells_up = 0

      do i = 0, this_meshblock % ptr % sx - 1
        i_glob = i + this_meshblock % ptr % x0
        do j = 0, this_meshblock % ptr % sy - 1
          j_glob = j + this_meshblock % ptr % y0
          do k = 0, this_meshblock % ptr % sz - 1
            if (((i_glob .le. injector_i1_glob) .and. (i_glob .gt. old_x1_glob)) .or. ((i_glob .ge. injector_i2_glob) .and. (i_glob .lt. old_x2_glob))) then
              average_electron = average_electron + lg_arr(i, j, k)
              ncells = ncells + 1
            end if
            if (((i_glob .gt. injector_i1_glob) .and. (i_glob .le. injector_i1_up)) .or. ((i_glob .lt. injector_i2_glob) .and. (i_glob .ge. injector_i2_up))) then
              average_electron_up = average_electron_up + lg_arr(i, j, k)
              ncells_up = ncells_up + 1
            end if
          end do
        end do
      end do

      lackofparticles = ncells * nUP_elec - average_electron
      addedelectron = 0
      lackofparticles_up = ncells_up * nUP_elec - average_electron_up

      do i = 0, this_meshblock % ptr % sx - 1
        i_glob = i + this_meshblock % ptr % x0
        do j = 0, this_meshblock % ptr % sy - 1
          j_glob = j + this_meshblock % ptr % y0
          do k = 0, this_meshblock % ptr % sz - 1
            if ((((i_glob .le. injector_i1_glob) .and. (i_glob .gt. old_x1_glob)) .or. &
                 ((i_glob .ge. injector_i2_glob) .and. (i_glob .lt. old_x2_glob))) .and. &
                (lackofparticles .gt. addedelectron) .and. (lackofparticles_up .ge. 0)) then
              if (lg_arr(i, j, k) .lt. nUP_elec) then
                ! underdensity
                density_int = INT(nUP_elec - lg_arr(i, j, k))
                density_frac = (nUP_elec - lg_arr(i, j, k)) - REAL(density_int)
                do n = 1, density_int
                  ! inject integer amount of particles
                  dx = random(dseed); dy = random(dseed)
                  dz = 0.5
                  ! inject electron
                  call createParticle(up_lecs, INT(i, 2), INT(j, 2), INT(k, 2), dx, dy, dz, 0.0, 0.0, 0.0)
                  ! inject ions/positron
                  call createParticle(up_ions, INT(i, 2), INT(j, 2), INT(k, 2), dx, dy, dz, 0.0, 0.0, 0.0)
                  addedelectron = addedelectron + 1
                end do
                if (random(dseed) .lt. density_frac) then
                  ! inject float amount of particles
                  ! if we need to add 0.3 particles per cell, we add a particle
                  ! in each cell with probability 0.3
                  dx = random(dseed); dy = random(dseed)
                  dz = 0.5
                  call createParticle(up_lecs, INT(i, 2), INT(j, 2), INT(k, 2), dx, dy, dz, 0.0, 0.0, 0.0)
                  call createParticle(up_ions, INT(i, 2), INT(j, 2), INT(k, 2), dx, dy, dz, 0.0, 0.0, 0.0)
                  addedelectron = addedelectron + 1
                end if
              end if
            end if
          end do
        end do
      end do

      !
      ! Clearing and resetting particles
      !
      do s = 1, nspec
        do ti = 1, species(s) % tile_nx
          do tj = 1, species(s) % tile_ny
            do tk = 1, species(s) % tile_nz
              do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
                x_glob = REAL(species(s) % prtl_tile(ti, tj, tk) % xi(p) + this_meshblock % ptr % x0) &
                         + species(s) % prtl_tile(ti, tj, tk) % dx(p)
                y_glob = REAL(species(s) % prtl_tile(ti, tj, tk) % yi(p) + this_meshblock % ptr % y0) &
                         + species(s) % prtl_tile(ti, tj, tk) % dy(p)
                ! clear particles beyond injectors just ONCE in the beginning
                if (step .eq. 1) then
                  if ((x_glob .le. injector_sx * 0.9) .or. (x_glob .gt. global_mesh % sx - injector_sx * 0.9)) then
                    species(s) % prtl_tile(ti, tj, tk) % proc(p) = -1
                  end if
                end if

#ifndef ABSORB
                ! remove particles within the absorbing layer manually (with a given probability)
                if ((y_glob .lt. ds_abs) .or. (y_glob .ge. global_mesh % sy - ds_abs)) then
                  if (random(dseed) .lt. CC / ds_abs) then
                    species(s) % prtl_tile(ti, tj, tk) % proc(p) = -1
                  end if
                end if
#endif

                if (s .eq. ph_index) then
                  if (.not. fake_photons) then ! * * * * dynamic photons
                    ! remove old photons
                    age = REAL(step) - species(s) % prtl_tile(ti, tj, tk) % payload1(p)
                    if (age .gt. ph_maxage) then
                      species(s) % prtl_tile(ti, tj, tk) % proc(p) = -1
                    end if
                  else ! * * * * photons as static background
                    ! reset scattered photon momenta
                    if (species(s) % prtl_tile(ti, tj, tk) % payload2(p) .ge. 0.0) then
                      call generateRandomMomentum(kx, ky, kz)
                      species(s) % prtl_tile(ti, tj, tk) % payload2(p) = -1.0
                      species(s) % prtl_tile(ti, tj, tk) % u(p) = kx
                      species(s) % prtl_tile(ti, tj, tk) % v(p) = ky
                      species(s) % prtl_tile(ti, tj, tk) % w(p) = kz
                    end if
                  end if
                end if
              end do
            end do
          end do
        end do
      end do
    end if

    ! inject photons
    ! (only if photons modeled dynamically)
    if (.not. fake_photons) then
      if (ph_fraction .gt. 0) then
        ncells = global_mesh % sx * global_mesh % sy * global_mesh % sz
        nphotons = INT(REAL(ph_fraction, 8) * REAL(ppc0, 8) * REAL(ncells, 8))
        do n = 1, nphotons
          call injectPhoton(step)
        end do
      end if
    end if

  end subroutine userParticleBoundaryConditions

  real function planckSample()
    implicit none
    real :: prob, n
    real :: rnd
    rnd = random(dseed)
    n = 0.0
    prob = 0.0
    do while ((prob .lt. rnd) .and. (n .lt. 40.0))
      n = n + 1.0
      prob = prob + 1.0 / (1.20206 * n**3)
    end do
    planckSample = -log(random(dseed) * random(dseed) * random(dseed) + 1e-16) / n
    return
  end function planckSample

  subroutine generateRandomDirection(kx, ky, kz)
    implicit none
    real, intent(out) :: kx, ky, kz
    real :: rand_costh, rand_phi
    rand_costh = 2.0 * random(dseed) - 1.0
    rand_phi = 2.0 * M_PI * random(dseed)
    kx = sqrt(1.0 - rand_costh**2) * cos(rand_phi)
    ky = sqrt(1.0 - rand_costh**2) * sin(rand_phi)
    kz = rand_costh
  end subroutine generateRandomDirection

  subroutine generateRandomMomentum(kx, ky, kz)
    implicit none
    real, intent(out) :: kx, ky, kz
    real :: energy
    call generateRandomDirection(kx, ky, kz)
    if (use_planckian) then
      energy = ph_temperature * planckSample()
    else
      energy = ph_energy
    end if
    kx = kx * energy
    ky = ky * energy
    kz = kz * energy
  end subroutine generateRandomMomentum

  subroutine injectPhoton(step)
    implicit none
    integer, optional, intent(in) :: step
    real :: xg, yg, zg, kx, ky, kz

    xg = random(dseed) * REAL(global_mesh % sx)
    yg = random(dseed) * REAL(global_mesh % sy)
    zg = 0.5
    call generateRandomMomentum(kx, ky, kz)
    call injectParticleGlobally(ph_index, xg, yg, zg, &
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
    real :: sx_glob, x_glob, delta_x, y_glob
    real(kind=8) :: injector_x1, injector_x2
    integer :: i, j, k
    integer :: i_glob, injector_i1_glob, injector_i2_glob, j_glob
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
    logical :: updateE_, updateB_
    real :: lambdaIJ, lambdaIpJ, lambdaIJp, lambdaIpJp
    real :: bx_target, by_target, bz_target, ex_target, ey_target, ez_target

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

    if ((step .ge. open_boundaries) .and. (open_boundaries .ge. 0)) then
      absorb_y = 1
      boundary_y = 0
      call reassignNeighborsForAll(meshblocks)
    end if

    if (boundary_x .ne. 1) then
      ! --------------------------------------------------------------------------
      !                        boundaries near the injector
      ! --------------------------------------------------------------------------

      injector_x1 = REAL(injector_sx, 8) - 1.0d-7
      injector_x2 = REAL(global_mesh % sx, 8) - REAL(injector_sx, 8) + 1.0d-7

      injector_i1_glob = INT(injector_x1) + 1
      injector_i2_glob = INT(injector_x2)

      if ((injector_i1_glob .lt. this_meshblock % ptr % x0 + this_meshblock % ptr % sx) .or. &
          (injector_i2_glob .ge. this_meshblock % ptr % x0)) then
        ! reset fields left and right from the injectors
        if (updateB_) then
          sx_glob = REAL(global_mesh % sx)
          do i = -NGHOST, this_meshblock % ptr % sx - 1 + NGHOST
            i_glob = i + this_meshblock % ptr % x0
            x_glob = REAL(i_glob) + 0.5
            if (i_glob .le. injector_i1_glob) then
              if (simple_bc) then
                bx(i, :, :) = 0.0
                by(i, :, :) = tanh((x_glob - cs_x * sx_glob) / current_width)
                bz(i, :, :) = bguide
              else
                delta_x = 4.0 * REAL(i_glob) / MAX(REAL(injector_i1_glob), 0.1)
                bx(i, :, :) = tanh(delta_x) * bx(i, :, :)
                by(i, :, :) = (1.0 - tanh(delta_x)) * tanh((x_glob - cs_x * sx_glob) / current_width) + tanh(delta_x) * by(i, :, :)
                bz(i, :, :) = tanh(delta_x) * bz(i, :, :) + (1.0 - tanh(delta_x)) * bguide
              end if
            else if (i_glob .ge. injector_i2_glob) then
              if (simple_bc) then
                bx(i, :, :) = 0.0
                by(i, :, :) = tanh((x_glob - cs_x * sx_glob) / current_width)
                bz(i, :, :) = bguide
              else
                delta_x = 4.0 * REAL(global_mesh % sx - 1 - i_glob) / MAX(REAL(global_mesh % sx - 1 - injector_i2_glob), 0.1)
                bx(i, :, :) = tanh(delta_x) * bx(i, :, :)
                by(i, :, :) = (1.0 - tanh(delta_x)) * tanh((x_glob - cs_x * sx_glob) / current_width) + tanh(delta_x) * by(i, :, :)
                bz(i, :, :) = tanh(delta_x) * bz(i, :, :) + (1.0 - tanh(delta_x)) * bguide
              end if
            end if
          end do
        end if
        if (updateE_) then
          sx_glob = REAL(global_mesh % sx)
          do i = -NGHOST, this_meshblock % ptr % sx - 1 + NGHOST
            i_glob = i + this_meshblock % ptr % x0
            if (i_glob .lt. injector_i1_glob) then
              if (simple_bc) then
                ex(i, :, :) = 0.0
                ey(i, :, :) = 0.0
                ez(i, :, :) = 0.0
              else
                delta_x = 4.0 * REAL(i_glob) / MAX(REAL(injector_i1_glob), 0.1)
                ex(i, :, :) = tanh(delta_x) * ex(i, :, :)
                ey(i, :, :) = tanh(delta_x) * ey(i, :, :)
                ez(i, :, :) = tanh(delta_x) * ez(i, :, :)
              end if
            else if (i_glob .gt. injector_i2_glob) then
              if (simple_bc) then
                ex(i, :, :) = 0.0
                ey(i, :, :) = 0.0
                ez(i, :, :) = 0.0
              else
                delta_x = 4.0 * REAL(global_mesh % sx - 1 - i_glob) / MAX(REAL(global_mesh % sx - 1 - injector_i2_glob), 0.1)
                ex(i, :, :) = tanh(delta_x) * ex(i, :, :)
                ey(i, :, :) = tanh(delta_x) * ey(i, :, :)
                ez(i, :, :) = tanh(delta_x) * ez(i, :, :)
              end if
            end if
          end do
        end if
      end if
    end if

#ifndef ABSORB
    ! --------------------------------------------------------------------------
    !                 user-defined absorbing boundary conditions
    ! --------------------------------------------------------------------------
    if (boundary_y .ne. 1) then
      do i = -NGHOST, this_meshblock % ptr % sx - 1 + NGHOST
        do j = -NGHOST, this_meshblock % ptr % sy - 1 + NGHOST
          i_glob = i + this_meshblock % ptr % x0
          j_glob = j + this_meshblock % ptr % y0
          x_glob = REAL(i_glob)
          y_glob = REAL(j_glob)

          ! i, j + 1/2
          bx_target = 0.1 * tanh((y_glob + 0.5 - 0.5 * REAL(global_mesh % sy)) / current_width)
          ! i + 1/2, j
          by_target = tanh((x_glob + 0.5 - 0.5 * REAL(global_mesh % sx)) / current_width)
          ! i + 1/2, j + 1/2
          bz_target = bguide
          ! i + 1/2, j
          ex_target = 0.0
          ! i, j + 1/2
          ey_target = -0.1 * bguide * tanh((x_glob - 0.5 * REAL(global_mesh % sx)) / current_width)
          ! i, j
          ez_target = 0.1

          if (y_glob .lt. ds_abs) then
            lambdaIJ = (2.0 * CC / ds_abs) * abs((ds_abs - y_glob) / ds_abs)**3.0
            lambdaIpJ = lambdaIJ
            lambdaIJp = (2.0 * CC / ds_abs) * abs((ds_abs - y_glob - 0.5) / ds_abs)**3.0
            lambdaIpJp = lambdaIJp
          else if (y_glob .gt. REAL(global_mesh % sy) - ds_abs) then
            lambdaIJ = (2.0 * CC / ds_abs) * abs((REAL(global_mesh % sy) - ds_abs - y_glob) / ds_abs)**3.0
            lambdaIpJ = lambdaIJ
            lambdaIJp = (2.0 * CC / ds_abs) * abs((REAL(global_mesh % sy) - ds_abs - y_glob - 0.5) / ds_abs)**3.0
            lambdaIpJp = lambdaIJp
          end if

          if ((y_glob .lt. ds_abs) .or. (y_glob .gt. REAL(global_mesh % sy) - ds_abs)) then
            if (updateB_) then
              bx(i, j, :) = exp(-lambdaIJp) * bx(i, j, :) + (1.0 - exp(-lambdaIJp)) * bx_target
              by(i, j, :) = exp(-lambdaIpJ) * by(i, j, :) + (1.0 - exp(-lambdaIpJ)) * by_target
              bz(i, j, :) = exp(-lambdaIpJp) * bz(i, j, :) + (1.0 - exp(-lambdaIpJp)) * bz_target
            end if
            if (updateE_) then
              ex(i, j, :) = exp(-lambdaIpJ) * ex(i, j, :) + (1.0 - exp(-lambdaIpJ)) * ex_target
              ey(i, j, :) = exp(-lambdaIJp) * ey(i, j, :) + (1.0 - exp(-lambdaIJp)) * ey_target
              ez(i, j, :) = exp(-lambdaIJ) * ez(i, j, :) + (1.0 - exp(-lambdaIJ)) * ez_target
            end if
          end if
        end do
      end do
    end if
#endif
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
