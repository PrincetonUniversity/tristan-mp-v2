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
#ifdef USROUTPUT
  use m_writeusroutput
#endif
  implicit none

!--- PRIVATE variables -----------------------------------------!
  real, private :: nCS_over_nUP, current_width, upstream_T, cs_x, cs_x1, cs_x2
  real, private :: boost_y_Gamma, boost_y_beta
  real, private :: injector_sz, injector_betax, measure_x, indent
  real, private :: fraction_ions
  integer, private :: injector_reset_interval, open_boundaries, no_cooling
  integer, private :: cs_lecs, cs_ions, up_lecs, up_ions
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
    if (boundary_z .ne. 1) then
      call getInput('problem', 'injector_sz', injector_sz)
    else
      call throwError("ERROR: 3D reconnection allows 1 current sheet only")
    end if
    call getInput('problem', 'fraction_ions', fraction_ions, 0.1)
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
    if (boundary_z .eq. 1) then
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
    real :: distance, sz_glob

    sz_glob = REAL(global_mesh % sz)
    distance = (z_glob - 0.5 * sz_glob)
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
        rad2 = (z_glob - dummy1)**2 + (x_glob - dummy3)**2
        userSpatialDistribution = 1.0 / (cosh((z_glob - dummy1) / dummy2))**2 * &
                                  (1.0 - exp(-rad2 / (5.0 * dummy2)**2))
      else
        userSpatialDistribution = 1.0 / (cosh((z_glob - dummy1) / dummy2))**2
      end if
    else
      call throwError("ERROR: variable not present in `userSpatialDistribution()`")
    end if
    return
  end function userSpatialDistribution

  subroutine userInitParticles()
    implicit none
    real :: nUP_elec, nUP_pos, nUP_ions, nCS_elec, nCS_pos, nCS_ions
    type(region) :: back_region
    real :: sx_glob, sy_glob, sz_glob, shift_gamma, shift_beta, current_sheet_T
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
    sz_glob = REAL(global_mesh % sz)

    back_region % x_min = 0.0
    back_region % y_min = 0.0
    back_region % z_min = injector_sz
    back_region % x_max = sx_glob
    back_region % y_max = sy_glob
    back_region % z_max = sz_glob - injector_sz
    call fillRegionWithThermalPlasma(back_region, (/up_lecs, up_ions/), 2, nUP_pos, upstream_T)

    if (nCS_over_nUP .ne. 0) then
      shift_beta = sqrt(sigma) * c_omp / (current_width * nCS_over_nUP)
      if (shift_beta .ge. 1) then
        call throwError('ERROR: `shift_beta` >= 1 in `userInitParticles()`')
      end if
      shift_gamma = 1.0 / sqrt(1.0 - shift_beta**2)
      current_sheet_T = 0.5 * sigma / (nCS_over_nUP + 1)

      back_region % x_min = 0.0
      back_region % x_max = sx_glob

      back_region % y_min = 0.0
      back_region % y_max = sy_glob

      back_region % z_min = sz_glob * cs_x - 10 * current_width
      back_region % z_max = sz_glob * cs_x + 10 * current_width
      if (perturb) then
        call fillRegionWithThermalPlasma(back_region, (/cs_lecs, cs_ions/), 2, nCS_pos, current_sheet_T, &
                                         shift_gamma=shift_gamma, shift_dir=3, &
                                         spat_distr_ptr=spat_distr_ptr, &
                                         dummy1=cs_x * sz_glob, dummy2=current_width, dummy3=cs_x * sx_glob)
      else
        call fillRegionWithThermalPlasma(back_region, (/cs_lecs, cs_ions/), 2, nCS_pos, current_sheet_T, &
                                         shift_gamma=shift_gamma, shift_dir=3, &
                                         spat_distr_ptr=spat_distr_ptr, &
                                         dummy1=cs_x * sz_glob, dummy2=current_width)
      end if
    end if
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob
    real :: z_glob, sz_glob
    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = 0; by(:, :, :) = 0; bz(:, :, :) = 0
    jx(:, :, :) = 0; jy(:, :, :) = 0; jz(:, :, :) = 0

    k = 0
    sz_glob = REAL(global_mesh % sz)
    do i = -NGHOST, this_meshblock % ptr % sz - 1 + NGHOST
      i_glob = i + this_meshblock % ptr % z0
      z_glob = REAL(i_glob) + 0.5
      bx(:, :, i) = tanh((z_glob - cs_x * sz_glob) / current_width)
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
    real :: z_glob
    real(kind=8) :: old_z1, old_z2
    real :: nUP_elec, nUP_pos, nUP_ions
    real :: ux, uy, uz, gamma, vpart
    integer :: s, ti, tj, tk, p, nUP_tot
    integer :: injector_i1_glob, injector_i2_glob, old_z1_glob, old_z2_glob
    type(region) :: back_region
    real :: average_electron
    integer :: ncells, lackofparticles, addedelectron
    integer :: i, j, k, k_glob, density_int, n, j_glob
    integer :: injector_i1_up, injector_i2_up
    integer :: ncells_up, lackofparticles_up, average_electron_up
    real :: density_frac, dx, dy, dz
    integer, optional, intent(in) :: step
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()

    if ((step .ge. open_boundaries) .and. (open_boundaries .ge. 0)) then
      boundary_x = 0
    end if

    if ((step .eq. open_boundaries) .and. (open_boundaries .ge. 0)) then
      call reassignNeighborsForAll(meshblocks)
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

    if (boundary_z .ne. 1) then
      ! reset the injector position every once in a while
      old_z1 = injector_sz - 1.0d-7
      old_z2 = REAL(global_mesh % sz, 8) - injector_sz + 1.0d-7
      ! inject background particles at the injectors' positions
      nUP_elec = 0.5 * ppc0
      nUP_pos = 0.5 * ppc0

      call computeDensity(1, reset=.true., ds=0, charge=.false.)

      old_z1_glob = INT(old_z1) + 1
      old_z2_glob = INT(old_z2)
      injector_i1_glob = old_z1_glob + 1
      injector_i2_glob = old_z2_glob - 1
      injector_i1_up = old_z1_glob + 10
      injector_i2_up = old_z2_glob - 10

      ncells = 0
      average_electron = 0

      average_electron_up = 0
      ncells_up = 0

      do i = 0, this_meshblock % ptr % sx - 1
        do j = 0, this_meshblock % ptr % sy - 1
          do k = 0, this_meshblock % ptr % sz - 1
            k_glob = k + this_meshblock % ptr % z0
            if (((k_glob .le. injector_i1_glob) .and. (k_glob .gt. old_z1_glob)) .or. ((k_glob .ge. injector_i2_glob) .and. (k_glob .lt. old_z2_glob))) then
              average_electron = average_electron + lg_arr(i, j, k)
              ncells = ncells + 1
            end if
            if (((k_glob .gt. injector_i1_glob) .and. (k_glob .le. injector_i1_up)) .or. ((k_glob .lt. injector_i2_glob) .and. (k_glob .ge. injector_i2_up))) then
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
        do j = 0, this_meshblock % ptr % sy - 1
          do k = 0, this_meshblock % ptr % sz - 1
            k_glob = k + this_meshblock % ptr % z0
            if ((((k_glob .le. injector_i1_glob) .and. (k_glob .gt. old_z1_glob)) .or. &
                 ((k_glob .ge. injector_i2_glob) .and. (k_glob .lt. old_z2_glob))) .and. &
                (lackofparticles .gt. addedelectron) .and. (lackofparticles_up .ge. 0)) then
              if (lg_arr(i, j, k) .lt. nUP_elec) then
                ! underdensity
                density_int = INT(nUP_elec - lg_arr(i, j, k))
                density_frac = (nUP_elec - lg_arr(i, j, k)) - REAL(density_int)
                vpart = 0.0
                do n = 1, density_int
                  ! inject integer amount of particles
                  dx = random(dseed); dy = random(dseed)
                  dz = random(dseed)
                  ! inject electron
                  call createParticle(1, INT(i, 2), INT(j, 2), INT(k, 2), dx, dy, dz, vpart, 0.0, 0.0)
                  addedelectron = addedelectron + 1
                  ! inject ions/positron
                  call createParticle(2, INT(i, 2), INT(j, 2), INT(k, 2), dx, dy, dz, vpart, 0.0, 0.0)
                end do
                if (random(dseed) .lt. density_frac) then
                  ! inject float amount of particles
                  ! if we need to add 0.3 particles per cell, we add a particle
                  ! in each cell with probability 0.3
                  dx = random(dseed); dy = random(dseed)
                  dz = random(dseed)
                  call createParticle(1, INT(i, 2), INT(j, 2), INT(k, 2), dx, dy, dz, vpart, 0.0, 0.0)
                  addedelectron = addedelectron + 1
                  call createParticle(2, INT(i, 2), INT(j, 2), INT(k, 2), dx, dy, dz, vpart, 0.0, 0.0)
                end if
              end if
            end if
          end do
        end do
      end do
    end if
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    real :: sz_glob, z_glob, delta_x
    integer :: i, j, k
    integer :: i_glob, injector_i1_glob, injector_i2_glob
    real(kind=8) :: injector_z1, injector_z2
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

    if ((step .ge. open_boundaries) .and. (open_boundaries .ge. 0)) then
      boundary_x = 0
    end if

    if ((step .eq. open_boundaries) .and. (open_boundaries .ge. 0)) then
      call reassignNeighborsForAll(meshblocks)
    end if

    if (boundary_z .ne. 1) then

      injector_z1 = REAL(injector_sz, 8) - 1.0d-7
      injector_z2 = REAL(global_mesh % sz, 8) - REAL(injector_sz, 8) + 1.0d-7

      injector_i1_glob = INT(injector_z1) + 1
      injector_i2_glob = INT(injector_z2)

      if ((injector_i1_glob .lt. this_meshblock % ptr % z0 + this_meshblock % ptr % sz) .or. &
          (injector_i2_glob .ge. this_meshblock % ptr % z0)) then
        ! reset fields left and right from the injectors
        if (updateB_) then
          sz_glob = REAL(global_mesh % sz)
          do i = -NGHOST, this_meshblock % ptr % sz - 1 + NGHOST
            i_glob = i + this_meshblock % ptr % z0
            z_glob = REAL(i_glob) + 0.5
            if (i_glob .le. injector_i1_glob) then
              if (simple_bc) then
                bx(:, :, i) = tanh((z_glob - cs_x * sz_glob) / current_width)
                bz(:, :, i) = 0.0
                by(:, :, i) = 0.0
              else
                delta_x = 4.0 * REAL(i_glob) / MAX(REAL(injector_i1_glob), 0.1)
                bx(:, :, i) = (1.0 - tanh(delta_x)) * tanh((z_glob - cs_x * sz_glob) / current_width) + tanh(delta_x) * bx(:, :, i)
                by(:, :, i) = tanh(delta_x) * by(:, :, i)
                bz(:, :, i) = tanh(delta_x) * bz(:, :, i)
              end if
            else if (i_glob .ge. injector_i2_glob) then
              if (simple_bc) then
                bx(:, :, i) = tanh((z_glob - cs_x * sz_glob) / current_width)
                bz(:, :, i) = 0.0
                by(:, :, i) = 0.0
              else
                delta_x = 4.0 * REAL(global_mesh % sz - 1 - i_glob) / MAX(REAL(global_mesh % sz - 1 - injector_i2_glob), 0.1)
                bx(:, :, i) = (1.0 - tanh(delta_x)) * tanh((z_glob - cs_x * sz_glob) / current_width) + tanh(delta_x) * bx(:, :, i)
                by(:, :, i) = tanh(delta_x) * by(:, :, i)
                bz(:, :, i) = tanh(delta_x) * bz(:, :, i)
              end if
            end if
          end do
        end if
        if (updateE_) then
          sz_glob = REAL(global_mesh % sz)
          do i = -NGHOST, this_meshblock % ptr % sz - 1 + NGHOST
            i_glob = i + this_meshblock % ptr % z0
            if (i_glob .lt. injector_i1_glob) then
              if (simple_bc) then
                ex(:, :, i) = 0.0
                ey(:, :, i) = 0.0
                ez(:, :, i) = 0.0
              else
                delta_x = 4.0 * REAL(i_glob) / MAX(REAL(injector_i1_glob), 0.1)
                ex(:, :, i) = tanh(delta_x) * ex(:, :, i)
                ey(:, :, i) = tanh(delta_x) * ey(:, :, i)
                ez(:, :, i) = tanh(delta_x) * ez(:, :, i)
              end if
            else if (i_glob .gt. injector_i2_glob) then
              if (simple_bc) then
                ex(:, :, i) = 0.0
                ey(:, :, i) = 0.0
                ez(:, :, i) = 0.0
              else
                delta_x = 4.0 * REAL(global_mesh % sz - 1 - i_glob) / MAX(REAL(global_mesh % sz - 1 - injector_i2_glob), 0.1)
                ex(:, :, i) = tanh(delta_x) * ex(:, :, i)
                ey(:, :, i) = tanh(delta_x) * ey(:, :, i)
                ez(:, :, i) = tanh(delta_x) * ez(:, :, i)
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

    ! do yi = 0, ynum - 1
    !   y_bins(yi + 1) = 0.5 + REAL(yi)
    ! end do

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
