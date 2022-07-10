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
  implicit none

  ! real, private :: Etot_0
  ! logical, private :: first_step = .true.

contains
  subroutine initializeHistory()
    implicit none
    call getInput('output', 'hst_enable', hst_enable, .false.)
    call getInput('output', 'hst_interval', hst_interval, 1)
  end subroutine initializeHistory

  ! FIX2 total E^2, total B^2, total E_kin
  ! time, total mass, momenta, kinetic energies in three directions, total energy, and magnetic energies in three directions
  subroutine writeHistory(step)
    implicit none
    integer, intent(in) :: step
    if (.false.) print *, step
    ! real :: e_energy, b_energy, Eprt1, Eprt2, Etot
    ! real :: lec_e, ion_e, massive_e, massless_e
    ! real :: global_e_energy, global_b_energy
    ! real, allocatable :: prtl_energy(:)
    ! real, allocatable :: global_prtl_energy(:)
    ! integer :: ierr, s, column_width
    ! logical :: photons_present
    ! character :: vert_div, hor_div
    ! character(len=STR_MAX) :: FMT, dummy1, dummy2, dummy3, dummy4, filename
    ! procedure(getFMT), pointer :: get_fmt_ptr => null()

    ! get_fmt_ptr => getFMTForRealScientific

    ! allocate (global_prtl_energy(nspec))

    ! call computeEnergyInBox(e_energy, b_energy, prtl_energy)

    ! call MPI_REDUCE(e_energy, global_e_energy, 1, default_mpi_real, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    ! call MPI_REDUCE(b_energy, global_b_energy, 1, default_mpi_real, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    ! call MPI_REDUCE(prtl_energy, global_prtl_energy, nspec, default_mpi_real, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    ! if (mpi_rank .eq. 0) then
    !   global_e_energy = global_e_energy * B_norm**2 / 2
    !   global_b_energy = global_b_energy * B_norm**2 / 2
    !   do s = 1, nspec
    !     if (species(s) % m_sp .ne. 0) then
    !       global_prtl_energy(s) = global_prtl_energy(s) * unit_ch * CC**2 * species(s) % m_sp
    !     else
    !       global_prtl_energy(s) = global_prtl_energy(s) * unit_ch * CC**2
    !     end if
    !   end do

    !   ! if photons are present -- save massive/massless
    !   ! ... otherwise -- electrons/ions
    !   photons_present = .false.
    !   do s = 1, nspec
    !     if ((species(s) % m_sp .eq. 0) .and. (species(s) % ch_sp .eq. 0)) then
    !       photons_present = .true.
    !     end if
    !   end do

    !   if (.not. photons_present) then
    !     ! only save electron/ion energies
    !     lec_e = 0
    !     ion_e = 0
    !     do s = 1, nspec
    !       if (species(s) % ch_sp .lt. 0) then
    !         lec_e = lec_e + global_prtl_energy(s)
    !       else
    !         ion_e = ion_e + global_prtl_energy(s)
    !       end if
    !     end do
    !     Eprt1 = lec_e
    !     Eprt2 = ion_e
    !   else
    !     ! save massive/massless energies
    !     massive_e = 0
    !     massless_e = 0
    !     do s = 1, nspec
    !       if ((species(s) % m_sp .ne. 0) .or. (species(s) % ch_sp .ne. 0)) then
    !         massive_e = massive_e + global_prtl_energy(s)
    !       else
    !         massless_e = massless_e + global_prtl_energy(s)
    !       end if
    !     end do
    !     Eprt1 = massive_e
    !     Eprt2 = massless_e
    !   end if

    !   Etot = Eprt1 + Eprt2 + global_e_energy + global_b_energy

    !   filename = trim(output_dir_name)//'/history'
    !   column_width = 13
    !   vert_div = '|'
    !   hor_div = '='

    !   if (first_step) then
    !     first_step = .false.
    !     Etot_0 = Etot
    !     open (UNIT_history, file=filename, status="replace", form="formatted")
    !     do s = 1, column_width * 5
    !       FMT(s:s) = hor_div
    !     end do
    !     write (UNIT_history, '(A)') FMT(1:column_width * 5)

    !     write (UNIT_history, '(A10,A'//trim(STR(column_width * 3))//',A3)') vert_div, ' ', vert_div

    !     FMT = "(A7,A3,A"//trim(STR(column_width))//",A"//trim(STR(column_width))//",A"//trim(STR(column_width))//",A3)"
    !     write (UNIT_history, FMT) '[time]', vert_div, '[E^2]', '[B^2]', '[E^2+B^2]', vert_div

    !     FMT = "(A10,A"//trim(STR(column_width))//",A"//trim(STR(column_width))//",A" &
    !           //trim(STR(column_width))//",A3,A"//trim(STR(column_width))//")"
    !     write (UNIT_history, FMT) vert_div, '[% Etot]', '[% Etot]', '[% Etot]', vert_div, '[Etot]'

    !     write (UNIT_history, '(A10,A'//trim(STR(column_width * 3))//',A3)') vert_div, ' ', vert_div

    !     FMT = "(A10,A"//trim(STR(column_width))//",A"//trim(STR(column_width))//",A" &
    !           //trim(STR(column_width))//",A3,A"//trim(STR(column_width))//")"
    !     if (.not. photons_present) then
    !       write (UNIT_history, FMT) vert_div, '[lecs]', '[ions]', '[tot part]', vert_div, '[% dEtot]'
    !     else
    !       write (UNIT_history, FMT) vert_div, '[massive]', '[massless]', '[tot part]', vert_div, '[% dEtot]'
    !     end if

    !     FMT = "(A10,A"//trim(STR(column_width))//",A"//trim(STR(column_width))//",A"//trim(STR(column_width))//",A3)"
    !     write (UNIT_history, FMT) vert_div, '[% Etot]', '[% Etot]', '[% Etot]', vert_div

    !     write (UNIT_history, '(A10,A'//trim(STR(column_width * 3))//',A3)') vert_div, ' ', vert_div

    !     do s = 1, column_width * 5
    !       FMT(s:s) = hor_div
    !     end do
    !     write (UNIT_history, '(A)') FMT(1:column_width * 5)

    !     close (UNIT_history)
    !   end if

    !   open (UNIT_history, file=filename, status="old", position="append", form="formatted")

    !   write (UNIT_history, '(A10,A'//trim(STR(column_width * 3))//',A3)') vert_div, ' ', vert_div

    !   dummy1 = get_fmt_ptr(column_width)
    !   dummy2 = get_fmt_ptr(column_width)
    !   dummy3 = get_fmt_ptr(column_width)
    !   FMT = "(I7,A3,"//trim(dummy1)//","//trim(dummy2)//","//trim(dummy3)//",A3)"
    !   write (UNIT_history, FMT) step, vert_div, &
    !     global_e_energy, global_b_energy, &
    !     global_e_energy + global_b_energy, vert_div

    !   dummy1 = get_fmt_ptr(column_width - 1)
    !   dummy2 = get_fmt_ptr(column_width - 1)
    !   dummy3 = get_fmt_ptr(column_width - 1)
    !   dummy4 = get_fmt_ptr(column_width)
    !   FMT = "(A10,"//trim(dummy1)//",A1,"//trim(dummy2)//",A1,"//trim(dummy3)//",A1,A3,"//trim(dummy4)//")"
    !   write (UNIT_history, FMT) vert_div, global_e_energy * 100 / Etot, '%', &
    !     global_b_energy * 100 / Etot, '%', &
    !     (global_e_energy + global_b_energy) * 100 / Etot, '%', &
    !     vert_div, Etot

    !   write (UNIT_history, '(A10,A'//trim(STR(column_width * 3))//',A3)') vert_div, ' ', vert_div

    !   dummy1 = get_fmt_ptr(column_width)
    !   dummy2 = get_fmt_ptr(column_width)
    !   dummy3 = get_fmt_ptr(column_width)
    !   dummy4 = get_fmt_ptr(column_width - 1)
    !   FMT = "(A10,"//trim(dummy1)//","//trim(dummy2)//","//trim(dummy3)//",A3,"//trim(dummy4)//",A1)"
    !   write (UNIT_history, FMT) vert_div, Eprt1, Eprt2, Eprt1 + Eprt2, vert_div, &
    !     (Etot - Etot_0) * 100 / Etot_0, '%'

    !   dummy1 = get_fmt_ptr(column_width - 1)
    !   dummy2 = get_fmt_ptr(column_width - 1)
    !   dummy3 = get_fmt_ptr(column_width - 1)
    !   FMT = "(A10,"//trim(dummy1)//",A1,"//trim(dummy2)//",A1,"//trim(dummy3)//",A1,A3)"
    !   write (UNIT_history, FMT) vert_div, Eprt1 * 100 / Etot, '%', &
    !     Eprt2 * 100 / Etot, '%', &
    !     (Eprt1 + Eprt2) * 100 / Etot, '%', vert_div

    !   write (UNIT_history, '(A10,A'//trim(STR(column_width * 3))//',A3)') vert_div, ' ', vert_div

    !   do s = 1, column_width * 5
    !     FMT(s:s) = hor_div
    !   end do
    !   write (UNIT_history, '(A)') FMT(1:column_width * 5)

    !   close (UNIT_history)
    ! end if
    call printDiag("writeHistory()", 2)
  end subroutine writeHistory

  subroutine computeEnergyInBox(e_energy, b_energy, prtl_energy)
    implicit none
    real, intent(out) :: e_energy, b_energy
    real, allocatable, intent(out) :: prtl_energy(:)
    real :: gamma, ex0, ey0, ez0, bx0, by0, bz0
    integer :: s, ti, tj, tk, p
    integer(kind=2) :: i, j, k

    allocate (prtl_energy(nspec))

    do s = 1, nspec
      prtl_energy(s) = 0
      if (species(s) % m_sp .eq. 0) then
        do ti = 1, species(s) % tile_nx
          do tj = 1, species(s) % tile_ny
            do tk = 1, species(s) % tile_nz
              do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
                gamma = species(s) % prtl_tile(ti, tj, tk) % u(p)**2 + &
                        species(s) % prtl_tile(ti, tj, tk) % v(p)**2 + &
                        species(s) % prtl_tile(ti, tj, tk) % w(p)**2
                prtl_energy(s) = prtl_energy(s) + sqrt(gamma) * &
                                 species(s) % prtl_tile(ti, tj, tk) % weight(p)
              end do
            end do
          end do
        end do
      else
        do ti = 1, species(s) % tile_nx
          do tj = 1, species(s) % tile_ny
            do tk = 1, species(s) % tile_nz
              do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
                gamma = 1.0 + species(s) % prtl_tile(ti, tj, tk) % u(p)**2 + &
                        species(s) % prtl_tile(ti, tj, tk) % v(p)**2 + &
                        species(s) % prtl_tile(ti, tj, tk) % w(p)**2
                prtl_energy(s) = prtl_energy(s) + sqrt(gamma) * &
                                 species(s) % prtl_tile(ti, tj, tk) % weight(p)
              end do
            end do
          end do
        end do
      end if
    end do

    e_energy = 0
    b_energy = 0
    do i = 0, INT(this_meshblock % ptr % sx - 1, 2)
      do j = 0, INT(this_meshblock % ptr % sy - 1, 2)
        do k = 0, INT(this_meshblock % ptr % sz - 1, 2)
          call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
          call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
          e_energy = e_energy + (ex0**2 + ey0**2 + ez0**2)
          b_energy = b_energy + (bx0**2 + by0**2 + bz0**2)
        end do
      end do
    end do
  end subroutine computeEnergyInBox

end module m_writehistory
