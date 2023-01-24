module m_writehistory
  use m_globalnamespace
  use m_outputnamespace, only: hst_enable, hst_interval
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_readinput, only: getInput
  use m_helpers
  use m_restart, only: rst_simulation
  use m_qednamespace
  implicit none

  ! real, private :: Etot_0
  logical, private :: first_step = .true.

contains
  subroutine initializeHistory()
    implicit none
    call getInput('output', 'hst_enable', hst_enable, .false.)
    call getInput('output', 'hst_interval', hst_interval, 1)
  end subroutine initializeHistory

  subroutine writeHistory(step)
    implicit none
    integer, intent(in) :: step
    real :: e_energy(3), b_energy(3)
    real :: global_e_energy(3), global_b_energy(3)
    real, allocatable :: prtl_energy(:)
    real, allocatable :: global_prtl_energy(:)
    real, allocatable :: prtl_num(:)
    real, allocatable :: global_prtl_num(:)
    real :: Etot
    integer :: ierr, s, column_width
    logical :: photons_present
    character :: vert_div, hor_div
    character(len=STR_MAX) :: FMT, dummy1, dummy2, dummy3, filename
    procedure(getFMT), pointer :: get_fmt_ptr => null()
    real :: volume

    get_fmt_ptr => getFMTForRealScientific

    allocate (global_prtl_energy(nspec))
    allocate (global_prtl_num(nspec))
    global_prtl_energy(:) = 0
    global_prtl_num(:) = 0

    call computeEnergyInBox(e_energy, b_energy, prtl_energy, prtl_num)

    call MPI_REDUCE(e_energy, global_e_energy, 3, default_mpi_real, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(b_energy, global_b_energy, 3, default_mpi_real, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(prtl_energy, global_prtl_energy, nspec, default_mpi_real, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(prtl_num, global_prtl_num, nspec, default_mpi_real, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (mpi_rank .eq. 0) then
      ! mean energy densities in units of n_0 m_e c^2 and particle densities in units of n_0:
      volume = REAL(global_mesh % sx) * REAL(global_mesh % sy) * REAL(global_mesh % sz)
      global_e_energy(:) = global_e_energy(:) * 0.5 * (B_norm * CCINV**2 * c_omp)**2 / volume
      global_b_energy(:) = global_b_energy(:) * 0.5 * (B_norm * CCINV**2 * c_omp)**2 / volume
      do s = 1, nspec
        global_prtl_num(s) = global_prtl_num(s) / ppc0 / volume
        if (species(s) % m_sp .ne. 0) then
          global_prtl_energy(s) = global_prtl_energy(s) * species(s) % m_sp / ppc0 / volume
        else
          global_prtl_energy(s) = global_prtl_energy(s) / ppc0 / volume
        end if
      end do

      filename = trim(output_dir_name)//'/history'
      column_width = 13

      if (first_step .and. (.not. rst_simulation)) then
        first_step = .false.
        open (UNIT_history, file=filename, status="replace", form="formatted")

        FMT = "(A7,A"//trim(STR(column_width))//",A"//trim(STR(column_width))//",A"//trim(STR(column_width))//")"
        write (UNIT_history, FMT, advance='no') '[time]', '[Ex^2]', '[Ey^2]', '[Ez^2]'
        FMT = "(A"//trim(STR(column_width))//",A"//trim(STR(column_width))//",A"//trim(STR(column_width))//")"
        write (UNIT_history, FMT, advance='no') '[Bx^2]', '[By^2]', '[Bz^2]'
        do s = 1, nspec
          FMT = "(A"//trim(STR(column_width))//")"
          write (UNIT_history, FMT, advance='no') '[Esp'//trim(STR(s))//']'
        end do
        do s = 1, nspec
          FMT = "(A"//trim(STR(column_width))//")"
          write (UNIT_history, FMT, advance='no') '[Nsp'//trim(STR(s))//']'
        end do
        FMT = "(A"//trim(STR(column_width))//")"
        write (UNIT_history, FMT) '[Etot]'

        close (UNIT_history)
      end if  ! first_step

      open (UNIT_history, file=filename, status="old", position="append", form="formatted")

      dummy1 = get_fmt_ptr(global_e_energy(1), column_width)
      dummy2 = get_fmt_ptr(global_e_energy(2), column_width)
      dummy3 = get_fmt_ptr(global_e_energy(3), column_width)
      FMT = "(I7,"//trim(dummy1)//","//trim(dummy2)//","//trim(dummy3)//")"
      write (UNIT_history, FMT, advance='no') step, global_e_energy(1), global_e_energy(2), global_e_energy(3)
      dummy1 = get_fmt_ptr(global_b_energy(1), column_width)
      dummy2 = get_fmt_ptr(global_b_energy(2), column_width)
      dummy3 = get_fmt_ptr(global_b_energy(3), column_width)
      FMT = "("//trim(dummy1)//","//trim(dummy2)//","//trim(dummy3)//")"
      write (UNIT_history, FMT, advance='no') global_b_energy(1), global_b_energy(2), global_b_energy(3)
      do s = 1, nspec
        dummy1 = get_fmt_ptr(global_prtl_energy(s), column_width)
        FMT = "("//trim(dummy1)//")"
        write (UNIT_history, FMT, advance='no') global_prtl_energy(s)
      end do
      do s = 1, nspec
        dummy1 = get_fmt_ptr(global_prtl_num(s), column_width)
        FMT = "("//trim(dummy1)//")"
        write (UNIT_history, FMT, advance='no') global_prtl_num(s)
      end do
      Etot = global_e_energy(1) + global_e_energy(2) + global_e_energy(3) 
      Etot = Etot + global_b_energy(1) + global_b_energy(2) + global_b_energy(3) 
      do s = 1, nspec
        Etot = Etot + global_prtl_energy(s)
      end do
      dummy1 = get_fmt_ptr(Etot, column_width)
      FMT = "("//trim(dummy1)//")"
      write (UNIT_history, FMT) Etot

      close (UNIT_history)
    end if  !  mpi_rank = 0

    call printDiag("writeHistory()", 2)
  end subroutine writeHistory

  subroutine computeEnergyInBox(e_energy, b_energy, prtl_energy, prtl_num)
    implicit none
    real, intent(out) :: e_energy(3), b_energy(3)
    real, allocatable, intent(out) :: prtl_energy(:)
    real, allocatable, intent(out) :: prtl_num(:)
    real :: gamma, wei
    integer :: s, ti, tj, tk, p
    integer(kind=2) :: i, j, k

    allocate (prtl_energy(nspec))
    allocate (prtl_num(nspec))

    do s = 1, nspec
      prtl_energy(s) = 0
      prtl_num(s) = 0
      if (.not. species(s) % output_sp_hist) cycle
      if (species(s) % m_sp .eq. 0) then
        do tk = 1, species(s) % tile_nz
          do tj = 1, species(s) % tile_ny
            do ti = 1, species(s) % tile_nx
              do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
                gamma = species(s) % prtl_tile(ti, tj, tk) % u(p)**2 + &
                        species(s) % prtl_tile(ti, tj, tk) % v(p)**2 + &
                        species(s) % prtl_tile(ti, tj, tk) % w(p)**2
                wei = species(s) % prtl_tile(ti, tj, tk) % weight(p)
                prtl_energy(s) = prtl_energy(s) + sqrt(gamma) * wei
                prtl_num(s) = prtl_num(s) + wei
              end do
            end do
          end do
        end do
      else
        do tk = 1, species(s) % tile_nz
          do tj = 1, species(s) % tile_ny
            do ti = 1, species(s) % tile_nx
              do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
                gamma = 1.0 + species(s) % prtl_tile(ti, tj, tk) % u(p)**2 + &
                        species(s) % prtl_tile(ti, tj, tk) % v(p)**2 + &
                        species(s) % prtl_tile(ti, tj, tk) % w(p)**2
                prtl_energy(s) = prtl_energy(s) + sqrt(gamma) * &
                                 species(s) % prtl_tile(ti, tj, tk) % weight(p)
                prtl_num(s) = prtl_num(s) + species(s) % prtl_tile(ti, tj, tk) % weight(p)
              end do
            end do
          end do
        end do
      end if
    end do

    e_energy(:) = 0
    b_energy(:) = 0
    do k = 0, INT(this_meshblock % ptr % sz - 1, 2)
      do j = 0, INT(this_meshblock % ptr % sy - 1, 2)
        do i = 0, INT(this_meshblock % ptr % sx - 1, 2)
          e_energy(1) = e_energy(1) + ex(i, j, k)**2
          e_energy(2) = e_energy(2) + ey(i, j, k)**2
          e_energy(3) = e_energy(3) + ez(i, j, k)**2
          b_energy(1) = b_energy(1) + bx(i, j, k)**2
          b_energy(2) = b_energy(2) + by(i, j, k)**2
          b_energy(3) = b_energy(3) + bz(i, j, k)**2
        end do
      end do
    end do
  end subroutine computeEnergyInBox

end module m_writehistory
