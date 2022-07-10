module m_restart
#ifdef IFPORT
  use ifport, only: makedirqq
#endif
  use m_globalnamespace
  use m_outputnamespace, only: tot_output_index, slice_index
  use m_readinput, only: getInput
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_particlelogistics, only: allocateParticlesOnEmptyTile
  use m_helpers
  implicit none

  ! # of cpu simultaneously accessing filesystem
  integer :: rst_cpu_group
  logical :: rst_simulation = .false.
  logical :: rst_separate, rst_enable = .false.
  integer :: rst_interval, rst_start

  private :: writeFldRestart

contains
  subroutine initializeRestart()
    implicit none
    call getInput('restart', 'do_restart', rst_simulation, .false.)
    call getInput('restart', 'enable', rst_enable, .false.)
    call getInput('restart', 'start', rst_start, 0)
    call getInput('restart', 'interval', rst_interval, 10000)
    call getInput('restart', 'rewrite', rst_separate, .false.)
    call getInput('restart', 'cpu_group', rst_cpu_group, 50)
    rst_separate = (.not. rst_separate)
  end subroutine initializeRestart

  subroutine writeRestart(timestep)
    implicit none
    integer, intent(in) :: timestep
    character(len=STR_MAX) :: stepchar, rst_dir
#ifdef IFPORT
    logical :: result
#endif
    integer :: ierr, rnk, rnk_cnt, recv_count = 0
    integer :: dummy(1)
#ifdef MPI08
    type(MPI_STATUS) :: istat
#endif

#ifdef MPI
    integer :: istat(MPI_STATUS_SIZE)
#endif

    ! determine the directory to write
    if (rst_separate) then
      write (stepchar, "(i5.5)") INT(timestep / rst_interval)
    else
      write (stepchar, "(i5.5)") 0
    end if
    rst_dir = trim(restart_dir_name)//'/step_'//trim(stepchar)

    ! create directories for all cpus from rank #0
    if (mpi_rank .eq. 0) then
#ifdef IFPORT
      result = makedirqq(trim(rst_dir))
#else
      call system('mkdir -p '//trim(rst_dir))
#endif
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (mpi_rank .eq. 0) then
      dummy(1) = 1
      rnk = 1; rnk_cnt = 0
      do while (rnk .lt. mpi_size)
        do rnk_cnt = 0, rst_cpu_group - 1
          if (rnk + rnk_cnt .ge. mpi_size) cycle
          call MPI_SEND(dummy, 1, MPI_INTEGER, rnk + rnk_cnt, 1, MPI_COMM_WORLD, ierr)
        end do
        do rnk_cnt = 0, rst_cpu_group - 1
          if (rnk + rnk_cnt .ge. mpi_size) cycle
          call MPI_RECV(dummy, 1, MPI_INTEGER, rnk + rnk_cnt, 2, MPI_COMM_WORLD, istat, ierr)
          recv_count = recv_count + 1
        end do
        rnk = rnk + rst_cpu_group
      end do
      call writeFldRestart(timestep, rst_dir)
      call writePrtlRestart(timestep, rst_dir)
    else
      call MPI_RECV(dummy, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, istat, ierr)
      call writeFldRestart(timestep, rst_dir)
      call writePrtlRestart(timestep, rst_dir)
      dummy(1) = 2
      call MPI_SEND(dummy, 1, MPI_INTEGER, 0, 2, MPI_COMM_WORLD, ierr)
    end if

#ifdef DEBUG
    print *, 'rank #', mpi_rank, 'wrote restart'
    if (mpi_rank .eq. 0) then
      print *, 'received from:', recv_count
    end if
#endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (mpi_rank .eq. 0) then
      print *, '[DONE] Restart written into ', trim(rst_dir)
    end if
    call printDiag("writeRestart()", 2)
  end subroutine writeRestart

  subroutine writeFldRestart(timestep, rst_dir)
    implicit none
    integer, intent(in) :: timestep
    character(len=STR_MAX), intent(in) :: rst_dir
    character(len=STR_MAX) :: filename, mpichar

    write (mpichar, "(i8.8)") mpi_rank

    filename = trim(rst_dir)//'/flds.rst.'//trim(mpichar)
    open (UNIT_restart_fld, file=filename, status="replace", form="unformatted")
    write (UNIT_restart_fld) timestep, dseed, tot_output_index, slice_index
    write (UNIT_restart_fld) ex, ey, ez, bx, by, bz
    write (UNIT_restart_fld) CC, ppc0, c_omp, sigma
    close (UNIT_restart_fld)
  end subroutine writeFldRestart

  subroutine writePrtlRestart(timestep, rst_dir)
    ! [DEP_PRT]
    implicit none
    integer, intent(in) :: timestep
    character(len=STR_MAX), intent(in) :: rst_dir
    character(len=STR_MAX) :: filename, mpichar
    integer :: s, ti, tj, tk, num

    if (.false.) print *, timestep

    write (mpichar, "(i8.8)") mpi_rank

    filename = trim(rst_dir)//'/prtl.rst.'//trim(mpichar)
    open (UNIT_restart_prtl, file=filename, status="replace", form="unformatted")
    do s = 1, nspec
      write (UNIT_restart_prtl) species(s) % cntr_sp
      write (UNIT_restart_prtl) species(s) % tile_sx, species(s) % tile_sy, species(s) % tile_sz
      write (UNIT_restart_prtl) species(s) % tile_nx, species(s) % tile_ny, species(s) % tile_nz
      do ti = 1, species(s) % tile_nx
        do tj = 1, species(s) % tile_ny
          do tk = 1, species(s) % tile_nz
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % spec
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % maxptl_sp
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % npart_sp
            num = species(s) % prtl_tile(ti, tj, tk) % npart_sp
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % x1
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % x2
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % y1
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % y2
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % z1
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % z2

            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % xi(1:num)
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % yi(1:num)
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % zi(1:num)
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % dx(1:num)
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % dy(1:num)
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % dz(1:num)
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % u(1:num)
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % v(1:num)
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % w(1:num)
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % weight(1:num)
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % ind(1:num)
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % proc(1:num)

#ifdef PRTLPAYLOADS
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % payload1(1:num)
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % payload2(1:num)
            write (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % payload3(1:num)
#endif
          end do
        end do
      end do
    end do
    close (UNIT_restart_prtl)
  end subroutine writePrtlRestart

  subroutine restartSimulation()
    implicit none
    character(len=STR_MAX) :: mpichar, filename
    integer :: s, ti, tj, tk, num, ierr
    integer :: dummy_int1, dummy_int2, dummy_int3
    write (mpichar, "(i8.8)") mpi_rank

    if (mpi_rank .eq. 0) then
      print *, 'Reading restart data...'
    end if

    ! loading fields
    filename = trim(restart_from)//'/flds.rst.'//trim(mpichar)
    open (UNIT_restart_fld, file=filename, form="unformatted")
    rewind (UNIT_restart_fld)
    read (UNIT_restart_fld) start_timestep, dseed, tot_output_index, slice_index
    read (UNIT_restart_fld) ex, ey, ez, bx, by, bz
    read (UNIT_restart_fld) CC, ppc0, c_omp, sigma
    close (UNIT_restart_fld)
    start_timestep = start_timestep + 1

    call renormalizeUnits()

    if (mpi_rank .eq. 0) then
      print *, '`CC`, `ppc0`, `c_omp` & `sigma` are read from restart ...'
      print *, '... values read from input are ignored.'
    end if

    ! loading particles
    filename = trim(restart_from)//'/prtl.rst.'//trim(mpichar)
    open (UNIT_restart_prtl, file=filename, form="unformatted")
    do s = 1, nspec
      read (UNIT_restart_prtl) species(s) % cntr_sp
      ! check that the tile sizes are the same
      read (UNIT_restart_prtl) dummy_int1, dummy_int2, dummy_int3
      if ((dummy_int1 .ne. species(s) % tile_sx) .or. &
          (dummy_int2 .ne. species(s) % tile_sy) .or. &
          (dummy_int3 .ne. species(s) % tile_sz)) then
        call throwError('ERROR. Wrong tile sizes after the restart')
      end if
      ! check that the # of tiles are the same
      read (UNIT_restart_prtl) dummy_int1, dummy_int2, dummy_int3
      if ((dummy_int1 .ne. species(s) % tile_nx) .or. &
          (dummy_int2 .ne. species(s) % tile_ny) .or. &
          (dummy_int3 .ne. species(s) % tile_nz)) then
        call throwError('ERROR. Wrong # of tiles after the restart')
      end if

      ! loop through all the tiles and read
      do ti = 1, species(s) % tile_nx
        do tj = 1, species(s) % tile_ny
          do tk = 1, species(s) % tile_nz
            read (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % spec

            ! reallocate the tile if necessary
            read (UNIT_restart_prtl) dummy_int1
            if (dummy_int1 .ne. species(s) % prtl_tile(ti, tj, tk) % maxptl_sp) then
              call allocateParticlesOnEmptyTile(species(s) % prtl_tile(ti, tj, tk), dummy_int1)
            end if
            read (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % npart_sp
            num = species(s) % prtl_tile(ti, tj, tk) % npart_sp

            ! read out and check tile dimensions
            read (UNIT_restart_prtl) dummy_int1
            if (dummy_int1 .ne. species(s) % prtl_tile(ti, tj, tk) % x1) then
              call throwError('ERROR. Wrong tile dimension after the restart')
            end if
            read (UNIT_restart_prtl) dummy_int1
            if (dummy_int1 .ne. species(s) % prtl_tile(ti, tj, tk) % x2) then
              call throwError('ERROR. Wrong tile dimension after the restart')
            end if
            read (UNIT_restart_prtl) dummy_int1
            if (dummy_int1 .ne. species(s) % prtl_tile(ti, tj, tk) % y1) then
              call throwError('ERROR. Wrong tile dimension after the restart')
            end if
            read (UNIT_restart_prtl) dummy_int1
            if (dummy_int1 .ne. species(s) % prtl_tile(ti, tj, tk) % y2) then
              call throwError('ERROR. Wrong tile dimension after the restart')
            end if
            read (UNIT_restart_prtl) dummy_int1
            if (dummy_int1 .ne. species(s) % prtl_tile(ti, tj, tk) % z1) then
              call throwError('ERROR. Wrong tile dimension after the restart')
            end if
            read (UNIT_restart_prtl) dummy_int1
            if (dummy_int1 .ne. species(s) % prtl_tile(ti, tj, tk) % z2) then
              call throwError('ERROR. Wrong tile dimension after the restart')
            end if

            ! finally read out all the particles
            ! DEP_PRT [particle-dependent]
            read (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % xi(1:num)
            read (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % yi(1:num)
            read (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % zi(1:num)
            read (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % dx(1:num)
            read (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % dy(1:num)
            read (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % dz(1:num)
            read (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % u(1:num)
            read (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % v(1:num)
            read (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % w(1:num)
            read (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % weight(1:num)
            read (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % ind(1:num)
            read (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % proc(1:num)

#ifdef PRTLPAYLOADS
            read (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % payload1(1:num)
            read (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % payload2(1:num)
            read (UNIT_restart_prtl) species(s) % prtl_tile(ti, tj, tk) % payload3(1:num)
#endif
          end do
        end do
      end do
    end do
    close (UNIT_restart_prtl)

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    call printDiag("restartSimulation()", 1)
  end subroutine restartSimulation

end module m_restart
