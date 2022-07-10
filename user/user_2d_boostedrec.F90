module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
#ifdef USROUTPUT
  use m_writeusroutput
#endif
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real, private :: nCS_over_nUP, current_width, upstream_T, cs_x, cs_x1, cs_x2
  real, private :: boost_y_Gamma, boost_y_beta, measure_x
  integer, private :: cs_lecs, cs_ions, up_lecs, up_ions
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'current_width', current_width)
    call getInput('problem', 'upstream_T', upstream_T, 0.0)
    call getInput('problem', 'nCS_nUP', nCS_over_nUP, 0.0)
    call getInput('problem', 'cs_lecs', cs_lecs, 1)
    call getInput('problem', 'cs_ions', cs_ions, 2)
    call getInput('problem', 'up_lecs', up_lecs, 1)
    call getInput('problem', 'up_ions', up_ions, 2)
    call getInput('problem', 'measure_x', measure_x, 0.2)
    call getInput('problem', 'boost_Gamma', boost_y_Gamma, 1.0)
    if (boost_y_Gamma .ne. 1.0) then
      boost_y_beta = sqrt(1.0 - 1.0 / boost_y_Gamma**2)
    else
      boost_y_beta = 0.0
    end if

    if (current_width .le. 0.0) then
      call throwError("ERROR: `current_width` has to be > 0.")
    end if
    ! double periodic
    cs_x1 = 0.25; cs_x2 = 0.75
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
    real :: nUP, nCS
    type(region) :: back_region
    real :: sx_glob, sy_glob, shift_gamma, shift_beta, current_sheet_T
    integer :: s, ti, tj, tk, p
    real :: ux, uy, uz, gamma
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    nUP = 0.5 * ppc0
    nCS = 0.5 * ppc0 * nCS_over_nUP

    sx_glob = REAL(global_mesh % sx)
    sy_glob = REAL(global_mesh % sy)

    ! background is NOT boosted
    !back_region%x_min = 0.0
    !back_region%y_min = 0.0
    !back_region%x_max = sx_glob
    !back_region%y_max = sy_glob
    !call fillRegionWithThermalPlasma(back_region, (/up_lecs, up_ions/), 2, nUP, 0.0)

    if (nCS_over_nUP .ne. 0) then
      ! this is already in the lab frame
      shift_beta = sqrt(sigma) * c_omp / (current_width * nCS_over_nUP)
      if (shift_beta .ge. 1) then
        call throwError('ERROR: `shift_beta` >= 1 in `userInitParticles()`')
      end if
      shift_gamma = 1.0 / sqrt(1.0 - shift_beta**2)
      current_sheet_T = 0.5 * sigma / (nCS_over_nUP)

      back_region % x_min = sx_glob * cs_x1 - 10 * current_width
      back_region % x_max = sx_glob * cs_x1 + 10 * current_width
      back_region % y_min = 0.0
      back_region % y_max = sy_glob
      call fillRegionWithThermalPlasma(back_region, (/cs_lecs, cs_ions/), 2, nCS, current_sheet_T, &
                                       shift_gamma=shift_gamma, shift_dir=3, &
                                       spat_distr_ptr=spat_distr_ptr, &
                                       dummy1=cs_x1 * sx_glob, dummy2=current_width)
      back_region % x_min = sx_glob * cs_x2 - 10 * current_width
      back_region % x_max = sx_glob * cs_x2 + 10 * current_width
      call fillRegionWithThermalPlasma(back_region, (/cs_lecs, cs_ions/), 2, nCS, current_sheet_T, &
                                       shift_gamma=shift_gamma, shift_dir=-3, &
                                       spat_distr_ptr=spat_distr_ptr, &
                                       dummy1=cs_x2 * sx_glob, dummy2=current_width)
    end if

    back_region % x_min = 0.0
    back_region % y_min = 0.0
    back_region % x_max = sx_glob
    back_region % y_max = sy_glob
    call fillRegionWithThermalPlasma(back_region, (/up_lecs, up_ions/), 2, nUP, 0.0)

    ! loop over all particles
    do s = 1, nspec
      do ti = 1, species(s) % tile_nx
        do tj = 1, species(s) % tile_ny
          do tk = 1, species(s) % tile_nz
            do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
              ux = species(s) % prtl_tile(ti, tj, tk) % u(p)
              uy = species(s) % prtl_tile(ti, tj, tk) % v(p)
              uz = species(s) % prtl_tile(ti, tj, tk) % w(p)
              gamma = sqrt(1.0 + ux**2 + uy**2 + uz**2)
              uy = (boost_y_Gamma * uy + boost_y_beta * boost_y_Gamma * gamma)
              species(s) % prtl_tile(ti, tj, tk) % v(p) = uy
            end do
          end do
        end do
      end do
    end do

  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, i_glob
    real :: x_glob, sx_glob
    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = 0; by(:, :, :) = 0; bz(:, :, :) = 0
    jx(:, :, :) = 0; jy(:, :, :) = 0; jz(:, :, :) = 0

    sx_glob = REAL(global_mesh % sx)
    do i = -NGHOST, this_meshblock % ptr % sx - 1 + NGHOST
      i_glob = i + this_meshblock % ptr % x0
      x_glob = REAL(i_glob) + 0.5
      by(i, :, :) = tanh((x_glob - cs_x1 * sx_glob) / current_width) - &
                    tanh((x_glob - cs_x2 * sx_glob) / current_width) - 1.0
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
    integer, optional, intent(in) :: step
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
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
