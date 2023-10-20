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
  real, private :: cs_overdensity, cs_width, up_temperature
  real, private :: injector_padding, measure_x
  real(kind=8), private :: injector_x1_fld, injector_x2_fld
  integer, private :: injector_reset_interval, open_boundaries
  integer, private :: cs_lecs, cs_ions, cs_heavy, up_lecs, up_ions, up_heavy
  real, private :: b_guide
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    ! guide field
    call getInput('problem', 'b_guide', b_guide, 0.0)

    ! current sheet
    call getInput('problem', 'cs_width', cs_width)
    call getInput('problem', 'cs_overdensity', cs_overdensity, 0.0)
    call getInput('problem', 'cs_lecs', cs_lecs, 1)
    call getInput('problem', 'cs_ions', cs_ions, 2)

    ! upstream
    call getInput('problem', 'up_temperature', up_temperature)
    call getInput('problem', 'up_lecs', up_lecs, 3)
    call getInput('problem', 'up_ions', up_ions, 4)

    ! replenisher
    if (boundary_x .ne. 1) then
      call getInput('problem', 'injector_padding', injector_padding)
      if (injector_padding .lt. nfilter + 4) then
        print *, 'WARNING: injector_padding < nfilter + 4, setting injector_padding = nfilter + 4'
        injector_padding = nfilter + 4
      end if
    end if

    ! outflow boundaries
    call getInput('problem', 'open_boundaries', open_boundaries, -1)
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
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    nUP_elec = 0.5 * ppc0
    nCS_elec = nUP_elec * cs_overdensity

    nUP_pos = 0.5 * ppc0
    nCS_pos = nUP_pos * cs_overdensity

    sx_glob = REAL(global_mesh % sx)
    sy_glob = REAL(global_mesh % sy)

    back_region % x_min = 0.0
    back_region % y_min = 0.0
    back_region % x_max = sx_glob
    back_region % y_max = sy_glob
    call fillRegionWithThermalPlasma(back_region, (/up_lecs, up_ions/), 2, nUP_pos, up_temperature)

    if (cs_overdensity .ne. 0) then
      shift_beta = sqrt(sigma) * c_omp / (cs_width * cs_overdensity)
      if (shift_beta .ge. 1) then
        call throwError('ERROR: `shift_beta` >= 1 in `userInitParticles()`')
      end if
      shift_gamma = 1.0 / sqrt(1.0 - shift_beta**2)
      current_sheet_T = 0.5 * sigma * (1.0 + b_guide**2) / cs_overdensity

      back_region % x_min = sx_glob * 0.5 - 10 * cs_width
      back_region % x_max = sx_glob * 0.5 + 10 * cs_width
      back_region % y_min = 0
      back_region % y_max = sy_glob
      call fillRegionWithThermalPlasma(back_region, (/cs_lecs, cs_ions/), 2, nCS_pos, current_sheet_T, &
                                       shift_gamma=shift_gamma, shift_dir=3, &
                                       spat_distr_ptr=spat_distr_ptr, &
                                       dummy1=0.5 * sx_glob, dummy2=cs_width, dummy3=0.5 * sy_glob)
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
    integer :: ncells
    real :: x_glob, injector_padding_flds
    real(kind=8) :: dens_imin, dens_imax

    call computeNpart(up_lecs, reset=.true., ds=0)
    call computeNpart(up_ions, reset=.false., ds=0)

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
        call fillRegionWithThermalPlasma(imin_inj_region, (/up_lecs, up_ions/), 2, 0.5 * REAL(ppc0 - dens_imin), up_temperature)
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
        call fillRegionWithThermalPlasma(imax_inj_region, (/up_lecs, up_ions/), 2, 0.5 * REAL(ppc0 - dens_imax), up_temperature)
      end if
    end if

    injector_padding_flds = REAL(injector_padding - nfilter)
    do s = 1, nspec
      do ti = 1, species(s) % tile_nx
        do tj = 1, species(s) % tile_ny
          do tk = 1, species(s) % tile_nz
            do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
              x_glob = REAL(species(s) % prtl_tile(ti, tj, tk) % xi(p) + this_meshblock % ptr % x0) &
                       + species(s) % prtl_tile(ti, tj, tk) % dx(p)
              if ((x_glob .lt. injector_padding)) then
                species(s) % prtl_tile(ti, tj, tk) % u(p) = -0.1
                species(s) % prtl_tile(ti, tj, tk) % v(p) = 0
                species(s) % prtl_tile(ti, tj, tk) % w(p) = 0
              else if (x_glob .ge. global_mesh % sx - injector_padding) then
                species(s) % prtl_tile(ti, tj, tk) % u(p) = 0.1
                species(s) % prtl_tile(ti, tj, tk) % v(p) = 0
                species(s) % prtl_tile(ti, tj, tk) % w(p) = 0
              end if
              if ((x_glob .lt. injector_padding_flds) .or. (x_glob .ge. global_mesh % sx - injector_padding_flds)) then
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
    real :: kappa, injector_padding_flds
    real :: x1min, x1max, y1min, y1max, y2min, y2max

    if ((step .ge. open_boundaries) .and. (open_boundaries .ge. 0)) then
      absorb_y = 1
      boundary_y = 0
      call reassignNeighborsForAll(meshblocks)
    end if

    ! --------------------------------------------------------------------------
    !                        boundaries near the injector
    ! --------------------------------------------------------------------------
    injector_padding_flds = REAL(injector_padding - nfilter)
    x1min = injector_padding_flds
    x1max = REAL(global_mesh % sx) - injector_padding_flds

    do i = -NGHOST, this_meshblock % ptr % sx - 1 + NGHOST
      i_glob = i + this_meshblock % ptr % x0
      x_glob = REAL(i_glob)
      kappa = 10.0

      by_target = tanh(((x_glob + 0.5) - 0.5 * REAL(global_mesh % sx)) / cs_width)
      bz_target = b_guide

      if ((x_glob .lt. x1min) .or. (x_glob .ge. x1max)) then
        if (x_glob .lt. x1min) then
          lambdaIJ = kappa * (abs(x1min - x_glob) / (1.5 * injector_padding_flds))**3
          lambdaIpJ = kappa * (abs(x1min - (x_glob + 0.5)) / (1.5 * injector_padding_flds))**3
          lambdaIJp = lambdaIJ
          lambdaIpJp = lambdaIpJ
        else if (x_glob .ge. x1max) then
          lambdaIJ = kappa * (abs(x_glob - x1max) / (1.5 * injector_padding_flds))**3
          lambdaIpJ = kappa * (abs((x_glob + 0.5) - x1max) / (1.5 * injector_padding_flds))**3
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

  !--- user-specific output -----------------------------------!
#ifdef USROUTPUT
  subroutine userOutput(step)
    implicit none
    integer, optional, intent(in) :: step
    integer :: root_rank = 0
    real, allocatable :: y_bins(:), ExB_arr(:), ExB_arr_global(:)
    ! real                          :: dr, x_glob, y_glob, z_glob, r_glob
    real :: dummy_x, dummy_y, dummy_z, dummy
    ! real, allocatable             :: sum_ExBr_f(:), sum_f(:), sum_ExBr_f_global(:), sum_f_global(:)
    ! integer                       :: ri, rnum = 50, i, j, k, ierr
    integer :: x_bin, yi, ynum, i, j, k, ierr

    ynum = INT(global_mesh % sy)

    ! allocate(y_bins(ynum))
    allocate (ExB_arr(ynum))
    allocate (ExB_arr_global(ynum))
    ExB_arr(:) = 0.0

    if (this_meshblock % ptr % x0 .eq. 0) then
      x_bin = INT(measure_x * global_mesh % sx)
      i = x_bin; k = 0
      do j = 0, this_meshblock % ptr % sy - 1
        dummy_x = -(ez(i, j, k) * by(i, j, k)) + ey(i, j, k) * bz(i, j, k)
        dummy = bx(i, j, k)**2 + by(i, j, k)**2 + bz(i, j, k)**2

        yi = j + this_meshblock % ptr % y0
        ExB_arr(yi + 1) = dummy_x / dummy
      end do
    end if

    call MPI_REDUCE(ExB_arr, ExB_arr_global, ynum, default_mpi_real, MPI_SUM, root_rank, MPI_COMM_WORLD, ierr)

    if (mpi_rank .eq. root_rank) then
      call writeUsrOutputTimestep(step)
      ! call writeUsrOutputArray('y', y_bins)
      call writeUsrOutputArray('ExB', ExB_arr_global)
      call writeUsrOutputEnd()
    end if
  end subroutine userOutput

  logical function userExcludeParticles(s, ti, tj, tk, p)
    implicit none
    integer, intent(in) :: s, ti, tj, tk, p
    userExcludeParticles = .true.
  end function userExcludeParticles
#endif
  !............................................................!

end module m_userfile
