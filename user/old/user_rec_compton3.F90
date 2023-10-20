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
  real, private :: nCS_over_nUP, current_width, upstream_T
  real, private :: injector_sx, measure_x
  real(kind=8), private :: injector_x1_fld, injector_x2_fld
  integer, private :: injector_reset_interval, open_boundaries
  integer, private :: cs_lecs, cs_ions, cs_heavy, up_lecs, up_ions, up_heavy
  real, private :: bguide
  ! photon variables
  integer, private :: ph_index, ph_EscInY_index, ph_EscInZ_index, ph_EscInX_index
  real, private :: ph_temperature, ph_fraction, ph_injdistance
  real, private :: ph_injector_x1, ph_injector_x2
  logical, private :: ph_planck
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  private :: planckSample
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

  subroutine generateRandomDirection(kx, ky, kz, positive_kx)
    implicit none
    real, intent(out) :: kx, ky, kz
    real :: rand_costh, rand_phi
    logical, optional, intent(in) :: positive_kx

    rand_costh = 2.0 * random(dseed) - 1.0

    if (present(positive_kx)) then
      if (positive_kx) then
        rand_phi = M_PI * random(dseed) - M_PI * 0.5
      else
        ! only negative kx
        rand_phi = M_PI * random(dseed) + M_PI * 0.5
      end if
    else
      rand_phi = 2.0 * M_PI * random(dseed)
    end if

    kx = sqrt(1.0 - rand_costh**2) * cos(rand_phi)
    ky = sqrt(1.0 - rand_costh**2) * sin(rand_phi)
    kz = rand_costh
  end subroutine generateRandomDirection

  subroutine injectPhoton(step)
    implicit none
    integer, optional, intent(in) :: step
    real :: xg, yg, zg, kx, ky, kz, energy

    if (random(dseed) .lt. 0.5) then
      xg = ph_injector_x1
    else
      xg = ph_injector_x2
    end if
    yg = random(dseed) * REAL(global_mesh % sy)
    zg = 0.5

    if (ph_planck) then
      energy = ph_temperature * planckSample()
    else
      energy = ph_temperature
    end if

    ! call generateRandomDirection(kx, ky, kz, xg .lt. REAL(global_mesh % sx) * 0.5)
    if (xg .lt. REAL(global_mesh % sx) * 0.5) then
      kx = 1.0
    else
      kx = -1.0
    end if
    ky = 0.0
    kz = 0.0
    call injectParticleGlobally(ph_index, xg, yg, zg, &
                                energy * kx, energy * ky, energy * kz, &
                                1.0, 0.0, 0.0, 0.0)
  end subroutine injectPhoton

  ! subroutine resetPhoton(step, s, ti, tj, tk, p)
  !   implicit none
  !   integer, intent(in) :: step, s, ti, tj, tk, p
  !   real :: energy, kx, ky, kz
  !   real :: xg, yg

  !   xg = REAL(species(s) % prtl_tile(ti, tj, tk) % xi(p) + this_meshblock % ptr % x0) &
  !        + species(s) % prtl_tile(ti, tj, tk) % dx(p)
  !   yg = REAL(species(s) % prtl_tile(ti, tj, tk) % yi(p) + this_meshblock % ptr % y0) &
  !        + species(s) % prtl_tile(ti, tj, tk) % dy(p)

  !   if (ph_planck) then
  !     energy = ph_temperature * planckSample()
  !   else
  !     energy = ph_temperature
  !   end if

  !   call generateRandomDirection(kx, ky, kz)

  !   species(s) % prtl_tile(ti, tj, tk) % payload1(p) = xg
  !   species(s) % prtl_tile(ti, tj, tk) % payload2(p) = yg
  !   species(s) % prtl_tile(ti, tj, tk) % payload3(p) = 0.5
  !   species(s) % prtl_tile(ti, tj, tk) % u(p) = energy * kx
  !   species(s) % prtl_tile(ti, tj, tk) % v(p) = energy * ky
  !   species(s) % prtl_tile(ti, tj, tk) % w(p) = energy * kz
  ! end subroutine resetPhoton

  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    ! guide field
    call getInput('problem', 'bguide', bguide, 0.0)

    ! current sheet
    call getInput('problem', 'current_width', current_width)
    call getInput('problem', 'nCS_nUP', nCS_over_nUP, 0.0)
    call getInput('problem', 'cs_lecs', cs_lecs, 1)
    call getInput('problem', 'cs_ions', cs_ions, 2)

    ! upstream
    call getInput('problem', 'upstream_T', upstream_T)
    call getInput('problem', 'up_lecs', up_lecs, 3)
    call getInput('problem', 'up_ions', up_ions, 4)

    ! replenisher
    if (boundary_x .ne. 1) then
      call getInput('problem', 'injector_sx', injector_sx)
      if (injector_sx .lt. nfilter + 4) then
        print *, 'WARNING: injector_sx < nfilter + 4, setting injector_sx = nfilter + 4'
        injector_sx = nfilter + 4
      end if
    end if

    ! outflow boundaries
    call getInput('problem', 'open_boundaries', open_boundaries, -1)

    ! photons
    call getInput('problem', 'ph_index', ph_index, 5)
    ! temperature for the Planckian photon distribution
    call getInput('problem', 'ph_EscInX_index', ph_EscInX_index, 6)
    call getInput('problem', 'ph_EscInY_index', ph_EscInY_index, 7)
    call getInput('problem', 'ph_EscInZ_index', ph_EscInZ_index, 8)
    call getInput('problem', 'ph_fraction', ph_fraction)
    ! negative temperature means that the monoenergetic distribution (delta function) is used
    call getInput('problem', 'ph_temperature', ph_temperature)
    call getInput('problem', 'ph_injdistance', ph_injdistance)
    ph_planck = (ph_temperature .gt. 0.0)
    ph_temperature = ABS(ph_temperature)

    ph_injector_x1 = REAL(global_mesh % sx) * (0.5 - ph_injdistance)
    ph_injector_x2 = REAL(global_mesh % sx) * (0.5 + ph_injdistance)
  end subroutine userReadInput

  function userSLBload(x_glob, y_glob, z_glob, &
                       dummy1, dummy2, dummy3)
    real :: userSLBload
    ! global coordinates
    real, intent(in), optional :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional :: dummy1, dummy2, dummy3
    real :: distance, sx_glob
    return
  end function

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

      back_region % x_min = sx_glob * 0.5 - 10 * current_width
      back_region % x_max = sx_glob * 0.5 + 10 * current_width
      back_region % y_min = 0
      back_region % y_max = sy_glob
      call fillRegionWithThermalPlasma(back_region, (/cs_lecs, cs_ions/), 2, nCS_pos, current_sheet_T, &
                                       shift_gamma=shift_gamma, shift_dir=3, &
                                       spat_distr_ptr=spat_distr_ptr, &
                                       dummy1=0.5 * sx_glob, dummy2=current_width, dummy3=0.5 * sy_glob)
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
      by(i, :, :) = tanh(((x_glob + 0.5) - 0.5 * sx_glob) / current_width)
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
    real :: t_inject
    real(kind=8) :: nphotons_r
    integer :: nphotons, n, ncells
    integer :: s, ti, tj, tk

    ! inject new photons
    t_inject = global_mesh % sy / CC
    ncells = global_mesh % sx * global_mesh % sy * global_mesh % sz
    nphotons_r = REAL(ph_fraction, 8) * REAL(ppc0, 8) * 0.5 * REAL(ncells, 8) / REAL(t_inject, 8)
    nphotons = INT(nphotons_r)
    do n = 1, nphotons
      call injectPhoton(step)
    end do
    if (random(dseed) .lt. (nphotons_r - REAL(nphotons, 8))) then
      call injectPhoton(step)
    end if

    ! remove escaping photons after output written

    do s = ph_EscInX_index, ph_EscInZ_index
      species(s) % move_sp = .false.
      species(s) % compton_sp = .false.
    end do
    if (modulo(step, tot_output_interval) .eq. 1) then
      do s = ph_EscInX_index, ph_EscInZ_index
        do ti = 1, species(s) % tile_nx
          do tj = 1, species(s) % tile_ny
            do tk = 1, species(s) % tile_nz
              species(s) % prtl_tile(ti, tj, tk) % npart_sp = 0
            end do
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
    integer :: npart_sp, ncells, ph_esc_index
    real :: x_glob, y_glob, injector_sx_flds
    real :: distance_z, energy, kx, ky, kz
    real :: xg, yg, zg, dummy1, dummy2
    logical :: leaving_in_x, leaving_in_y, leaving_in_z
    real(kind=8) :: dens_imin, dens_imax

    call computeNpart(up_lecs, reset=.true., ds=0)
    call computeNpart(up_ions, reset=.false., ds=0)

    imin_inj = INT(injector_sx)
    imax_inj = global_mesh % sx - INT(injector_sx)

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
        call fillRegionWithThermalPlasma(imin_inj_region, (/up_lecs, up_ions/), 2, 0.5 * REAL(ppc0 - dens_imin), upstream_T)
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
        call fillRegionWithThermalPlasma(imax_inj_region, (/up_lecs, up_ions/), 2, 0.5 * REAL(ppc0 - dens_imax), upstream_T)
      end if
    end if

    injector_sx_flds = REAL(injector_sx - nfilter)
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
                if ((x_glob .lt. injector_sx)) then
                  species(s) % prtl_tile(ti, tj, tk) % u(p) = -0.1
                  species(s) % prtl_tile(ti, tj, tk) % v(p) = 0
                  species(s) % prtl_tile(ti, tj, tk) % w(p) = 0
                else if (x_glob .ge. global_mesh % sx - injector_sx) then
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
                distance_z = abs(species(s) % prtl_tile(ti, tj, tk) % payload3(p))
                leaving_in_x = (x_glob .lt. ph_injector_x1 - 2.0) .or. (x_glob .ge. ph_injector_x2 + 2.0)
                leaving_in_y = (boundary_y .ne. 1) .and. ((y_glob .lt. 2.0 * ds_abs) .or. (y_glob .ge. REAL(global_mesh % sy) - 2.0 * ds_abs))
                leaving_in_z = distance_z .gt. global_mesh % sy * 0.5

                if (leaving_in_y) then
                  ph_esc_index = ph_EscInY_index
                else if (leaving_in_z) then
                  ph_esc_index = ph_EscInZ_index
                else if (leaving_in_x) then
                  ph_esc_index = ph_EscInX_index
                else
                  cycle
                end if

                ! remove photons after distance in Z
                call createParticle(ph_esc_index, &
                                    species(s) % prtl_tile(ti, tj, tk) % xi(p), &
                                    species(s) % prtl_tile(ti, tj, tk) % yi(p), &
                                    species(s) % prtl_tile(ti, tj, tk) % zi(p), &
                                    species(s) % prtl_tile(ti, tj, tk) % dx(p), &
                                    species(s) % prtl_tile(ti, tj, tk) % dy(p), &
                                    species(s) % prtl_tile(ti, tj, tk) % dz(p), &
                                    species(s) % prtl_tile(ti, tj, tk) % u(p), &
                                    species(s) % prtl_tile(ti, tj, tk) % v(p), &
                                    species(s) % prtl_tile(ti, tj, tk) % w(p), &
                                    weight=species(s) % prtl_tile(ti, tj, tk) % weight(p), &
                                    payload1=species(s) % prtl_tile(ti, tj, tk) % payload1(p), &
                                    payload2=species(s) % prtl_tile(ti, tj, tk) % payload2(p), &
                                    payload3=species(s) % prtl_tile(ti, tj, tk) % payload3(p))
                species(s) % prtl_tile(ti, tj, tk) % proc(p) = -1
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

    if ((step .ge. open_boundaries) .and. (open_boundaries .ge. 0)) then
      absorb_y = 1
      boundary_y = 0
      call reassignNeighborsForAll(meshblocks)
    end if

    ! --------------------------------------------------------------------------
    !                        boundaries near the injector
    ! --------------------------------------------------------------------------
    injector_sx_flds = REAL(injector_sx - nfilter)
    x1min = injector_sx_flds
    x1max = REAL(global_mesh % sx) - injector_sx_flds

    do i = -NGHOST, this_meshblock % ptr % sx - 1 + NGHOST
      i_glob = i + this_meshblock % ptr % x0
      x_glob = REAL(i_glob)
      kappa = 10.0

      by_target = tanh(((x_glob + 0.5) - 0.5 * REAL(global_mesh % sx)) / current_width)
      bz_target = bguide

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
          bx_target = 0.1 * tanh(((y_glob + 0.5) - 0.5 * REAL(global_mesh % sy)) / current_width)
          ! i + 1/2, j
          by_target = tanh(((x_glob + 0.5) - 0.5 * REAL(global_mesh % sx)) / current_width)
          ! i + 1/2, j + 1/2
          bz_target = bguide
          ! i + 1/2, j
          ex_target = 0.0
          ! i, j + 1/2
          ey_target = -0.1 * bguide * tanh((x_glob - 0.5 * REAL(global_mesh % sx)) / current_width)
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
    incr_pld1 = 0.0; incr_pld2 = 0.0; incr_pld3 = CC * w0 * over_e_temp
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
