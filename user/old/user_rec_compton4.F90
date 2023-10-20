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
  use m_outputnamespace, only: tot_output_interval
  use m_qednamespace, only: QED_tau0
#ifdef USROUTPUT
  use m_writeusroutput
#endif
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  ! current sheet parameters
  real, private :: cs_overdensity, cs_width, up_temperature, b_guide
  integer, private :: top_lecs, btm_lecs, top_ions, btm_ions, cs_lecs, cs_ions
  ! injector parameters
  real, private :: injector_padding
  real(kind=8), private :: injector_x1_fld, injector_x2_fld
  integer, private :: open_bc_y_step
  ! photon variables
  integer, private :: ph_index, ph_index_esc
  real, private :: ph_temperature, ph_fraction, ph_Lbox
  logical, private :: ph_planck
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution, planckSample
  !...............................................................!
contains
  ! --- Utilities ------------------------------------------------!
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

  subroutine injectPhoton(step)
    implicit none
    integer, optional, intent(in) :: step
    real :: xg, yg, zg, kx, ky, kz, energy

    ! inject in the midplane
    xg = 0.5 * REAL(global_mesh % sx)
    yg = random(dseed) * REAL(global_mesh % sy)
    zg = random(dseed) * ph_Lbox * 2.0

    if (ph_planck) then
      ! planck distribution
      energy = ph_temperature * planckSample()
    else
      ! monoenergetic
      energy = ph_temperature
    end if

    call generateRandomDirection(kx, ky, kz)
    call injectParticleGlobally(ph_index, xg, yg, 0.5, &
                                energy * kx, energy * ky, energy * kz, &
                                1.0, & ! < weight
                                zg, 0.0, 0.0)
    !  payloads:
    !    payload1 -- z position
    !    payload2 -- lifetime
    !    payload3 -- # of scatterings
  end subroutine injectPhoton

  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    ! guide field
    call getInput('problem', 'b_guide', b_guide, 0.0)

    ! current sheet
    call getInput('problem', 'cs_width', cs_width)
    call getInput('problem', 'cs_overdensity', cs_overdensity, 0.0)
    call getInput('problem', 'cs_lecs', cs_lecs, 5)
    call getInput('problem', 'cs_ions', cs_ions, 6)

    ! upstream
    call getInput('problem', 'up_temperature', up_temperature)
    call getInput('problem', 'top_lecs', top_lecs, 1)
    call getInput('problem', 'top_ions', top_ions, 2)
    call getInput('problem', 'btm_lecs', btm_lecs, 3)
    call getInput('problem', 'btm_ions', btm_ions, 4)

    ! replenisher
    if (boundary_x .ne. 1) then
      call getInput('problem', 'injector_padding', injector_padding)
      if (injector_padding .lt. nfilter + 4) then
        print *, 'WARNING: injector_padding < nfilter + 4, setting injector_padding = nfilter + 4'
        injector_padding = nfilter + 4
      end if
    end if

    ! outflow boundaries
    call getInput('problem', 'open_bc_y_step', open_bc_y_step, -1)

    ! photons
    call getInput('problem', 'ph_index', ph_index, 7)
    call getInput('problem', 'ph_index_esc', ph_index_esc, 8)
    call getInput('problem', 'ph_fraction', ph_fraction)
    ! temperature for the Planckian photon distribution
    ! negative temperature means that the monoenergetic distribution is used
    call getInput('problem', 'ph_temperature', ph_temperature)
    ph_planck = (ph_temperature .gt. 0.0)
    ph_temperature = ABS(ph_temperature)
    ph_Lbox = REAL(global_mesh % sy) / 2.0
  end subroutine userReadInput

  ! function userSLBload(x_glob, y_glob, z_glob, &
  !                      dummy1, dummy2, dummy3)
  !   real :: userSLBload
  !   ! global coordinates
  !   real, intent(in), optional :: x_glob, y_glob, z_glob
  !   ! global box dimensions
  !   real, intent(in), optional :: dummy1, dummy2, dummy3
  !   real :: distance, sx_glob
  !   return
  ! end function

  function userSpatialDistribution(x_glob, y_glob, z_glob, &
                                   dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3
    real :: rad2
    if (present(x_glob) .and. present(y_glob) .and. &
        present(dummy1) .and. present(dummy2) .and. &
        present(dummy3) .and. (dummy3 .ne. 0.0)) then
      ! underdensity at the center
      rad2 = (x_glob - dummy1)**2 + (y_glob - dummy3)**2
      userSpatialDistribution = 1.0 / (cosh((x_glob - dummy1) / dummy2))**2 * &
                                (1.0 - exp(-rad2 / (5.0 * dummy2)**2))
    else
      call throwError("ERROR: variable not present in `userSpatialDistribution()`")
    end if
    return
  end function userSpatialDistribution

  subroutine userInitParticles()
    implicit none
    ! real :: nUP_elec, nUP_pos, nCS_elec, nCS_pos
    type(region) :: back_region
    real :: sx_glob, sy_glob, shift_gamma, shift_beta, cs_temperature
    integer :: s, ti, tj, tk, p
    real :: ux, uy, uz, gamma
    integer :: n, ncells, nphotons
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    ! nUP_elec = 0.5 * ppc0
    ! nCS_elec = nUP_elec * nCS_over_nUP

    ! nUP_pos = 0.5 * ppc0
    ! nCS_pos = nUP_pos * nCS_over_nUP

    ! sx_glob = REAL(global_mesh % sx)
    ! sy_glob = REAL(global_mesh % sy)

    back_region % x_min = 0.0
    back_region % y_min = 0.0
    back_region % x_max = 0.5 * REAL(global_mesh % sx)
    back_region % y_max = REAL(global_mesh % sy)
    call fillRegionWithThermalPlasma(back_region, (/btm_lecs, btm_ions/), 2, 0.5 * ppc0, up_temperature)

    back_region % x_min = 0.5 * REAL(global_mesh % sx)
    back_region % x_max = REAL(global_mesh % sx)
    call fillRegionWithThermalPlasma(back_region, (/top_lecs, top_ions/), 2, 0.5 * ppc0, up_temperature)

    if (cs_overdensity .ne. 0) then
      shift_beta = sqrt(sigma) * c_omp / (cs_width * cs_overdensity)
      if (shift_beta .ge. 1) then
        call throwError('ERROR: `shift_beta` >= 1 in `userInitParticles()`')
      end if
      shift_gamma = 1.0 / sqrt(1.0 - shift_beta**2)
      cs_temperature = 0.5 * sigma * (1.0 + b_guide**2) / cs_overdensity

      back_region % x_min = REAL(global_mesh % sx) * 0.5 - 10 * cs_width
      back_region % x_max = REAL(global_mesh % sx) * 0.5 + 10 * cs_width
      back_region % y_min = 0
      back_region % y_max = REAL(global_mesh % sy)
      call fillRegionWithThermalPlasma(back_region, (/cs_lecs, cs_ions/), 2, 0.5 * ppc0, cs_temperature, &
                                       shift_gamma=shift_gamma, shift_dir=3, &
                                       spat_distr_ptr=spat_distr_ptr, &
                                       dummy1=0.5 * REAL(global_mesh % sx), &
                                       dummy2=cs_width, &
                                       dummy3=0.5 * REAL(global_mesh % sy))
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
      x_glob = REAL(i_glob)
      by(i, :, :) = tanh(((x_glob + 0.5) - 0.5 * sx_glob) / cs_width)
    end do
    bz(:, :, :) = b_guide
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
    real(kind=8) :: nphotons_r
    integer :: nphotons, n
    integer :: s, ti, tj, tk

    ! inject new photons
    nphotons_r = 4.0 * REAL(ppc0, 8) * ph_Lbox * CC * REAL(ph_fraction, 8)
    nphotons = INT(nphotons_r)
    do n = 1, nphotons
      call injectPhoton(step)
    end do
    if (random(dseed) .lt. (nphotons_r - REAL(nphotons, 8))) then
      call injectPhoton(step)
    end if

    ! remove escaping photons after output written

    species(ph_index_esc) % move_sp = .false.
    species(ph_index_esc) % compton_sp = .false.
    if (modulo(step, tot_output_interval) .eq. 1) then
      do ti = 1, species(ph_index_esc) % tile_nx
        do tj = 1, species(ph_index_esc) % tile_ny
          do tk = 1, species(ph_index_esc) % tile_nz
            species(ph_index_esc) % prtl_tile(ti, tj, tk) % npart_sp = 0
          end do
        end do
      end do
    end if
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
    integer, optional, intent(in) :: step
    type(region) :: imin_inj_region, imax_inj_region
    integer :: s, ti, tj, tk
    integer :: i, j, k, p
    integer :: i_glob
    integer :: imin_inj, imax_inj, imin_inj_local, imax_inj_local
    integer :: npart_sp, ncells
    real :: x_glob, y_glob, z_glob, injector_sx_flds
    real :: energy, kx, ky, kz
    real :: xg, yg, zg, dummy1, dummy2
    logical :: leaving_in_x, leaving_in_y, leaving_in_z
    real(kind=8) :: dens_imin, dens_imax

    call computeNpart(top_lecs, reset=.true., ds=0)
    call computeNpart(top_ions, reset=.false., ds=0)
    if (top_lecs .ne. btm_lecs) then
      call computeNpart(btm_lecs, reset=.false., ds=0)
      call computeNpart(btm_ions, reset=.false., ds=0)
    end if

    imin_inj = INT(injector_padding)
    imax_inj = global_mesh % sx - INT(injector_padding)

    dens_imin = 0.0
    dens_imax = 0.0
    if ((imin_inj .ge. this_meshblock % ptr % x0) .and. (imin_inj .lt. this_meshblock % ptr % x0 + this_meshblock % ptr % sx)) then
      ! imin contained within the current meshblock
      ncells = 0
      imin_inj_local = imin_inj - this_meshblock % ptr % x0
      do i = imin_inj_local, imin_inj_local + 10
        do j = 0, this_meshblock % ptr % sy - 1
          do k = 0, this_meshblock % ptr % sz - 1
            dens_imin = dens_imin + REAL(lg_arr(i, j, k), 8)
            ncells = ncells + 1
          end do
        end do
      end do
      dens_imin = dens_imin / REAL(ncells, 8)

      if (ppc0 .gt. dens_imin) then
        imin_inj_region % x_min = REAL(imin_inj)
        imin_inj_region % x_max = REAL(imin_inj) + 1.0
        imin_inj_region % y_min = this_meshblock % ptr % y0
        imin_inj_region % y_max = this_meshblock % ptr % y0 + this_meshblock % ptr % sy
        call fillRegionWithThermalPlasma(imin_inj_region, (/btm_lecs, btm_ions/), 2, &
                                         0.5 * REAL(ppc0 - dens_imin), up_temperature)
      end if
    end if

    if ((imax_inj .ge. this_meshblock % ptr % x0) .and. (imax_inj .lt. this_meshblock % ptr % x0 + this_meshblock % ptr % sx)) then
      ! imax contained within the current meshblock
      ncells = 0
      imax_inj_local = imax_inj - this_meshblock % ptr % x0
      do i = imax_inj_local - 11, imax_inj_local - 1
        do j = 0, this_meshblock % ptr % sy - 1
          do k = 0, this_meshblock % ptr % sz - 1
            dens_imax = dens_imax + REAL(lg_arr(i, j, k), 8)
            ncells = ncells + 1
          end do
        end do
      end do
      dens_imax = dens_imax / REAL(ncells, 8)

      if (ppc0 .gt. dens_imax) then
        imax_inj_region % x_min = REAL(imax_inj) - 1.0
        imax_inj_region % x_max = REAL(imax_inj)
        imax_inj_region % y_min = this_meshblock % ptr % y0
        imax_inj_region % y_max = this_meshblock % ptr % y0 + this_meshblock % ptr % sy
        call fillRegionWithThermalPlasma(imax_inj_region, (/top_lecs, top_ions/), 2, &
                                         0.5 * REAL(ppc0 - dens_imax), up_temperature)
      end if
    end if

    ! remove particles
    injector_sx_flds = REAL(injector_padding - nfilter)
    do s = 1, ph_index
      do ti = 1, species(s) % tile_nx
        do tj = 1, species(s) % tile_ny
          do tk = 1, species(s) % tile_nz
            npart_sp = species(s) % prtl_tile(ti, tj, tk) % npart_sp
            do p = 1, npart_sp
              x_glob = REAL(species(s) % prtl_tile(ti, tj, tk) % xi(p) + this_meshblock % ptr % x0) &
                       + species(s) % prtl_tile(ti, tj, tk) % dx(p)
              if (s .lt. ph_index) then
                ! . . . . . .
                ! . . . pairs
                ! . . . . . .
                ! advect pairs outwards at the injector boundaries
                if (boundary_x .eq. 1) then
                  cycle
                end if
                if ((x_glob .lt. injector_padding)) then
                  species(s) % prtl_tile(ti, tj, tk) % u(p) = -0.1
                  species(s) % prtl_tile(ti, tj, tk) % v(p) = 0
                  species(s) % prtl_tile(ti, tj, tk) % w(p) = 0
                else if (x_glob .ge. global_mesh % sx - injector_padding) then
                  species(s) % prtl_tile(ti, tj, tk) % u(p) = 0.1
                  species(s) % prtl_tile(ti, tj, tk) % v(p) = 0
                  species(s) % prtl_tile(ti, tj, tk) % w(p) = 0
                end if
                ! remove pairs once they are close to boundaries (with nfilter buffer)
                if ((x_glob .lt. injector_sx_flds) .or. (x_glob .ge. global_mesh % sx - injector_sx_flds)) then
                  species(s) % prtl_tile(ti, tj, tk) % proc(p) = -1
                end if
              else if (s .eq. ph_index) then
                ! . . . . . . .
                ! . . . photons
                ! . . . . . . .
                y_glob = REAL(species(s) % prtl_tile(ti, tj, tk) % yi(p) + this_meshblock % ptr % y0) &
                         + species(s) % prtl_tile(ti, tj, tk) % dy(p)
                z_glob = species(s) % prtl_tile(ti, tj, tk) % payload1(p)
                leaving_in_x = (x_glob .lt. REAL(global_mesh % sx) * 0.5 - ph_Lbox * 0.5) .or. &
                               (x_glob .ge. REAL(global_mesh % sx) * 0.5 + ph_Lbox * 0.5)
                leaving_in_y = (boundary_y .ne. 1) .and. ((y_glob .lt. 2.0 * ds_abs) .or. &
                                                          (y_glob .ge. REAL(global_mesh % sy) - 2.0 * ds_abs))
                leaving_in_z = (z_glob .lt. 0.0) .or. (z_glob .ge. 2 * ph_Lbox)

                if (leaving_in_x .or. leaving_in_y .or. leaving_in_z) then
                  ! remove photons
                  call createParticle(ph_index_esc, &
                                      species(s) % prtl_tile(ti, tj, tk) % xi(p), &
                                      species(s) % prtl_tile(ti, tj, tk) % yi(p), &
                                      species(s) % prtl_tile(ti, tj, tk) % zi(p), &
                                      species(s) % prtl_tile(ti, tj, tk) % dx(p), &
                                      species(s) % prtl_tile(ti, tj, tk) % dy(p), &
                                      species(s) % prtl_tile(ti, tj, tk) % dz(p), &
                                      species(s) % prtl_tile(ti, tj, tk) % u(p), &
                                      species(s) % prtl_tile(ti, tj, tk) % v(p), &
                                      species(s) % prtl_tile(ti, tj, tk) % w(p), &
                                      ind=species(s) % prtl_tile(ti, tj, tk) % ind(p), &
                                      proc=species(s) % prtl_tile(ti, tj, tk) % proc(p), &
                                      weight=species(s) % prtl_tile(ti, tj, tk) % weight(p), &
                                      payload1=species(s) % prtl_tile(ti, tj, tk) % payload1(p), &
                                      payload2=species(s) % prtl_tile(ti, tj, tk) % payload2(p), &
                                      payload3=species(s) % prtl_tile(ti, tj, tk) % payload3(p))
                  species(s) % prtl_tile(ti, tj, tk) % proc(p) = -1
                end if
              end if
            end do
          end do
        end do
      end do
    end do
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    real :: sx_glob, x_glob, delta_x, y_glob
    integer :: i, j, k
    integer :: i_glob, j_glob
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
    real :: lambdaIJ, lambdaIpJ, lambdaIJp, lambdaIpJp
    real :: bx_target, by_target, bz_target, ex_target, ey_target, ez_target
    real :: kappa, injector_sx_flds
    real :: x1min, x1max, y1min, y1max, y2min, y2max

    if ((step .ge. open_bc_y_step) .and. (open_bc_y_step .ge. 0)) then
      absorb_y = 1
      boundary_y = 0
      call reassignNeighborsForAll(meshblocks)
    end if

    ! --------------------------------------------------------------------------
    !                        boundaries near the injector
    ! --------------------------------------------------------------------------
    injector_sx_flds = REAL(injector_padding - nfilter)
    x1min = injector_sx_flds
    x1max = REAL(global_mesh % sx) - injector_sx_flds

    do i = -NGHOST, this_meshblock % ptr % sx - 1 + NGHOST
      i_glob = i + this_meshblock % ptr % x0
      x_glob = REAL(i_glob)
      kappa = 10.0

      by_target = tanh(((x_glob + 0.5) - 0.5 * REAL(global_mesh % sx)) / cs_width)
      bz_target = b_guide

      if ((x_glob .lt. x1min) .or. (x_glob .ge. x1max)) then
        if (x_glob .lt. x1min) then
          lambdaIJ = kappa * (abs(x1min - x_glob) / (1.5 * injector_sx_flds))**3
          lambdaIpJ = kappa * (abs(x1min - (x_glob + 0.5)) / (1.5 * injector_sx_flds))**3
          lambdaIJp = lambdaIJ
          lambdaIpJp = lambdaIpJ
        else if (x_glob .ge. x1max) then
          lambdaIJ = kappa * (abs(x_glob - x1max) / (1.5 * injector_sx_flds))**3
          lambdaIpJ = kappa * (abs((x_glob + 0.5) - x1max) / (1.5 * injector_sx_flds))**3
          lambdaIJp = lambdaIJ
          lambdaIpJp = lambdaIpJ
        end if
        bx(i, :, :) = exp(-lambdaIJp) * bx(i, :, :)
        by(i, :, :) = (1.0 - exp(-lambdaIpJ)) * by_target + exp(-lambdaIpJ) * by(i, :, :)
        bz(i, :, :) = (1.0 - exp(-lambdaIpJp)) * bz_target + exp(-lambdaIpJp) * bz(i, :, :)
        ex(i, :, :) = exp(-lambdaIpJ) * ex(i, :, :)
        ey(i, :, :) = exp(-lambdaIJp) * ey(i, :, :)
        ez(i, :, :) = exp(-lambdaIJ) * ez(i, :, :)
      end if
    end do

#ifndef ABSORB
    ! --------------------------------------------------------------------------
    !                 user-defined absorbing boundary conditions
    ! --------------------------------------------------------------------------
    if (boundary_y .ne. 1) then
      kappa = 10.0

      do i = -NGHOST, this_meshblock % ptr % sx - 1 + NGHOST
        do j = -NGHOST, this_meshblock % ptr % sy - 1 + NGHOST
          i_glob = i + this_meshblock % ptr % x0
          j_glob = j + this_meshblock % ptr % y0
          x_glob = REAL(i_glob)
          y_glob = REAL(j_glob)

          ! i, j + 1/2
          bx_target = 0.1 * tanh(((y_glob + 0.5) - 0.5 * REAL(global_mesh % sy)) / cs_width)
          ! i + 1/2, j
          by_target = tanh(((x_glob + 0.5) - 0.5 * REAL(global_mesh % sx)) / cs_width)
          ! i + 1/2, j + 1/2
          bz_target = b_guide
          ! i + 1/2, j
          ex_target = 0.0
          ! i, j + 1/2
          ey_target = -0.1 * b_guide * tanh((x_glob - 0.5 * REAL(global_mesh % sx)) / cs_width)
          ! i, j
          ez_target = 0.1
          !
          ! here's how the absorbing boundaries work:
          !
          !         y1min   y2min                                    y2max   y1max
          ! |-------|-------|----------------------------------------|-------|-------|
          !     ^       ^                                                ^       ^
          !     |       |                                            ___/        |
          ! e/b-target  |                                           /        e/b-target
          !             |                                          /
          !        exp(-lambda) * e/b + (1 - exp(-lambda)) * e/b-target

          y1min = ds_abs
          y2min = 2.0 * ds_abs
          y2max = REAL(global_mesh % sy) - 2.0 * ds_abs
          y1max = REAL(global_mesh % sy) - ds_abs

          if ((y_glob .lt. y1min) .or. (y_glob .ge. y1max)) then
            bx(i, j, :) = bx_target
            by(i, j, :) = by_target
            bz(i, j, :) = bz_target
            ex(i, j, :) = ex_target
            ey(i, j, :) = ey_target
            ez(i, j, :) = ez_target
          else if ((y_glob .lt. y2min) .or. (y_glob .ge. y2max)) then
            if (y_glob .lt. y2min) then
              lambdaIJ = kappa * abs((y2min - y_glob) / ds_abs)**3.0
              lambdaIJp = kappa * abs((y2min - (y_glob + 0.5)) / ds_abs)**3.0
            else if (y_glob .ge. y2max) then
              lambdaIJ = kappa * abs((y_glob - y2max) / ds_abs)**3.0
              lambdaIJp = kappa * abs(((y_glob + 0.5) - y2max) / ds_abs)**3.0
            end if
            lambdaIpJ = lambdaIJ
            lambdaIpJp = lambdaIJp

            bx(i, j, :) = exp(-lambdaIJp) * bx(i, j, :) + (1.0 - exp(-lambdaIJp)) * bx_target
            by(i, j, :) = exp(-lambdaIpJ) * by(i, j, :) + (1.0 - exp(-lambdaIpJ)) * by_target
            bz(i, j, :) = exp(-lambdaIpJp) * bz(i, j, :) + (1.0 - exp(-lambdaIpJp)) * bz_target
            ex(i, j, :) = exp(-lambdaIpJ) * ex(i, j, :) + (1.0 - exp(-lambdaIpJ)) * ex_target
            ey(i, j, :) = exp(-lambdaIJp) * ey(i, j, :) + (1.0 - exp(-lambdaIJp)) * ey_target
            ez(i, j, :) = exp(-lambdaIJ) * ez(i, j, :) + (1.0 - exp(-lambdaIJ)) * ez_target
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
    incr_pld1 = CC * w0 * over_e_temp; incr_pld2 = 1.0; incr_pld3 = 0.0
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
