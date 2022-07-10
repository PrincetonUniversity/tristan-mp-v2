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
  real, private :: nCS_over_nUP, current_width, upstream_T, cs_x, cs_x1, cs_x2
  real, private :: boost_y_Gamma, boost_y_beta
  real, private :: injector_sx, measure_x
  real(kind=8), private :: injector_x1_fld, injector_x2_fld
  real, private :: fraction_ions
  integer, private :: injector_reset_interval, open_boundaries, no_cooling
  integer, private :: cs_lecs, cs_ions, cs_heavy, up_lecs, up_ions, up_heavy
  logical, private :: perturb, simple_bc
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
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
    call getInput('problem', 'up_lecs', up_lecs, 1)
    call getInput('problem', 'up_ions', up_ions, 2)
    call getInput('problem', 'measure_x', measure_x, 0.2)
    call getInput('problem', 'open_boundaries', open_boundaries, -1)
    call getInput('problem', 'perturb', perturb, .false.)
    call getInput('problem', 'simple_bc', simple_bc, .true.)
    call getInput('problem', 'no_cooling', no_cooling, -1)
    if (current_width .lt. 0.0) then
      call throwError("ERROR: `current_width` has to be > 0.")
    end if
    if (boundary_x .eq. 1) then
      ! double periodic
      cs_x1 = 0.25; cs_x2 = 0.75
    else
      ! single current sheet
      cs_x = 0.5
    end if
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
        rad2 = (x_glob - dummy1)**2 + (y_glob - dummy3)**2
        userSpatialDistribution = 1.0 / (cosh((x_glob - dummy1) / dummy2))**2 * &
                                  (1.0 - exp(-rad2 / (5.0 * dummy2)**2))
      else
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
      current_sheet_T = 0.5 * sigma / nCS_over_nUP

      if (boundary_x .eq. 1) then
        back_region % x_min = sx_glob * cs_x1 - 10 * current_width
        back_region % x_max = sx_glob * cs_x1 + 10 * current_width
        back_region % y_min = 0
        back_region % y_max = sy_glob
        call fillRegionWithThermalPlasma(back_region, (/cs_lecs, cs_ions/), 2, nCS_pos, current_sheet_T, &
                                         shift_gamma=shift_gamma, shift_dir=3, &
                                         spat_distr_ptr=spat_distr_ptr, &
                                         dummy1=cs_x1 * sx_glob, dummy2=current_width)

        back_region % x_min = sx_glob * cs_x2 - 10 * current_width
        back_region % x_max = sx_glob * cs_x2 + 10 * current_width
        call fillRegionWithThermalPlasma(back_region, (/cs_lecs, cs_ions/), 2, nCS_pos, current_sheet_T, &
                                         shift_gamma=shift_gamma, shift_dir=-3, &
                                         spat_distr_ptr=spat_distr_ptr, &
                                         dummy1=cs_x2 * sx_glob, dummy2=current_width)

      else
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
      if (boundary_x .eq. 1) then
        by(i, :, :) = tanh((x_glob - cs_x1 * sx_glob) / current_width) - &
                      tanh((x_glob - cs_x2 * sx_glob) / current_width) - 1.0
      else
        by(i, :, :) = tanh((x_glob - cs_x * sx_glob) / current_width)
      end if
    end do
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
    real :: x_glob
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

    if ((step .ge. open_boundaries) .and. (open_boundaries .gt. 0)) then
      boundary_y = 0
    end if

    if ((step .ge. open_boundaries) .and. (open_boundaries .gt. 0)) then
      call reassignNeighborsForAll(meshblocks)
      open_boundaries = 0
    end if

#ifdef RADIATION
    if (step .lt. no_cooling) then
      species(1) % cool_sp = .false.
      species(2) % cool_sp = .false.
    else
      species(1) % cool_sp = .true.
      species(2) % cool_sp = .true.
    end if
#endif

    if (boundary_x .ne. 1) then
      old_x1 = REAL(injector_sx, 8) - 1.0d-7
      old_x2 = REAL(global_mesh % sx, 8) - REAL(injector_sx, 8) + 1.0d-7

      ! inject background particles at the injectors' positions
      nUP_elec = 0.5 * ppc0
      nUP_pos = 0.5 * ppc0

      call computeDensity(1, reset=.true., ds=0, charge=.false.)

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
                vpart = 0.0
                do n = 1, density_int
                  ! inject integer amount of particles
                  dx = random(dseed); dy = random(dseed)
                  dz = 0.5
                  ! inject electron
                  call createParticle(1, INT(i, 2), INT(j, 2), INT(k, 2), dx, dy, dz, vpart, 0.0, 0.0)
                  ! inject ions/positron
                  call createParticle(2, INT(i, 2), INT(j, 2), INT(k, 2), dx, dy, dz, vpart, 0.0, 0.0)
                  addedelectron = addedelectron + 1
                end do
                if (random(dseed) .lt. density_frac) then
                  ! inject float amount of particles
                  ! if we need to add 0.3 particles per cell, we add a particle
                  ! in each cell with probability 0.3
                  dx = random(dseed); dy = random(dseed)
                  dz = 0.5
                  call createParticle(1, INT(i, 2), INT(j, 2), INT(k, 2), dx, dy, dz, vpart, 0.0, 0.0)
                  call createParticle(2, INT(i, 2), INT(j, 2), INT(k, 2), dx, dy, dz, vpart, 0.0, 0.0)
                  addedelectron = addedelectron + 1
                end if
              end if
            end if
          end do
        end do
      end do

      if (step .eq. 1) then
        ! clear particles beyond injectors just once in the beginning
        do s = 1, nspec
          do ti = 1, species(s) % tile_nx
            do tj = 1, species(s) % tile_ny
              do tk = 1, species(s) % tile_nz
                do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
                  x_glob = REAL(species(s) % prtl_tile(ti, tj, tk) % xi(p) + this_meshblock % ptr % x0) &
                           + species(s) % prtl_tile(ti, tj, tk) % dx(p)
                  if ((x_glob .le. injector_sx * 0.9) .or. (x_glob .gt. global_mesh % sx - injector_sx * 0.9)) then
                    species(s) % prtl_tile(ti, tj, tk) % proc(p) = -1
                  end if
                end do
              end do
            end do
          end do
        end do

      end if
    end if
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    real :: sx_glob, x_glob, delta_x
    real(kind=8) :: injector_x1, injector_x2
    integer :: i, j, k
    integer :: i_glob, injector_i1_glob, injector_i2_glob
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

    if ((step .ge. open_boundaries) .and. (boundary_y .ne. 1)) then
      boundary_y = 0
      call reassignNeighborsForAll(meshblocks)
    end if

    if (boundary_x .ne. 1) then

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
                bz(i, :, :) = 0.0
                by(i, :, :) = tanh((x_glob - cs_x * sx_glob) / current_width)
              else
                delta_x = 4.0 * REAL(i_glob) / MAX(REAL(injector_i1_glob), 0.1)
                bx(i, :, :) = tanh(delta_x) * bx(i, :, :)
                bz(i, :, :) = tanh(delta_x) * bz(i, :, :)
                by(i, :, :) = (1.0 - tanh(delta_x)) * tanh((x_glob - cs_x * sx_glob) / current_width) + tanh(delta_x) * by(i, :, :)
              end if
            else if (i_glob .ge. injector_i2_glob) then
              if (simple_bc) then
                bx(i, :, :) = 0.0
                bz(i, :, :) = 0.0
                by(i, :, :) = tanh((x_glob - cs_x * sx_glob) / current_width)
              else
                delta_x = 4.0 * REAL(global_mesh % sx - 1 - i_glob) / MAX(REAL(global_mesh % sx - 1 - injector_i2_glob), 0.1)
                bx(i, :, :) = tanh(delta_x) * bx(i, :, :)
                bz(i, :, :) = tanh(delta_x) * bz(i, :, :)
                by(i, :, :) = (1.0 - tanh(delta_x)) * tanh((x_glob - cs_x * sx_glob) / current_width) + tanh(delta_x) * by(i, :, :)
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
  end subroutine userFieldBoundaryConditions
  !............................................................!

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
