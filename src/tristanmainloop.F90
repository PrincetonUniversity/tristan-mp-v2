#include "defs.F90"

module m_mainloop
  use m_globalnamespace
  use m_helpers
  use m_aux
  use m_writeoutput
  use m_fldsolver
  use m_mover
  use m_currentdeposit
  use m_exchangeparts
  use m_exchangefields
  use m_exchangecurrents
  use m_particlelogistics
  use m_filtering
  use m_userfile
  use m_errors

  ! extra physics
  #ifdef QED
    use m_qedphysics
  #endif
  implicit none

  integer       :: timestep

  real(kind=8)  :: t_fullstep, t_movestep,&
                 & t_depositstep, t_filterstep,&
                 & t_outputstep, t_fldexchstep,&
                 & t_prtlexchxtep, t_fldslvrstep,&
                 & t_usrfuncs
  #ifdef QED
    real (kind=8) :: t_qedstep
  #endif


  !--- PRIVATE functions -----------------------------------------!
  private :: makeReport
  !...............................................................!

  !--- PRIVATE variables -----------------------------------------!
  private :: t_fullstep, t_movestep, &
           & t_depositstep, t_filterstep, &
           & t_outputstep, t_fldexchstep,&
           & t_prtlexchxtep, t_fldslvrstep,&
           & t_usrfuncs
  #ifdef QED
    private :: t_qedstep
  #endif
  !...............................................................!
contains
  subroutine mainloop()
    implicit none
    integer       :: ierr, i
    integer       :: s, ti, tj, tk, p

    ! ADD needs to be changed for restart
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call printReport((mpi_rank .eq. 0), "Starting mainloop()")

    t_fullstep = 0;       t_movestep = 0
    t_depositstep = 0;    t_filterstep = 0
    t_outputstep = 0;     t_fldexchstep = 0
    t_prtlexchxtep = 0;   t_fldslvrstep = 0
    t_usrfuncs = 0

    #ifdef QED
      t_qedstep = 0
    #endif

    do timestep = 0, final_timestep
        t_fullstep = MPI_WTIME()

      ! MAINLOOP >
      !-------------------------------------------------
      ! User defined boundary conditions for fields
        t_usrfuncs = MPI_WTIME()
      call userFieldBoundaryConditions(timestep)
        t_usrfuncs = MPI_WTIME() - t_usrfuncs
      !.................................................

      !-------------------------------------------------
      ! Exchanging `E` and `B`-fields
        t_fldexchstep = MPI_WTIME()
      call exchangeFields(.true., .true.)
        t_fldexchstep = MPI_WTIME() - t_fldexchstep
      !.................................................

      !-------------------------------------------------
      ! Advancing 1st halfstep of `dB / dt = curl E`
        t_fldslvrstep = MPI_WTIME()
      call advanceBHalfstep()
        t_fldslvrstep = MPI_WTIME() - t_fldslvrstep
      !.................................................

      !-------------------------------------------------
      ! Exchanging `B`-fields
        t_fldexchstep = MPI_WTIME() - t_fldexchstep
      call exchangeFields(.false., .true.)
        t_fldexchstep = MPI_WTIME() - t_fldexchstep
      !.................................................

      !-------------------------------------------------
      ! QED business
      #ifdef QED
          t_qedstep = MPI_WTIME()
        call QEDstep(timestep)
        call clearGhostParticles()
          t_qedstep = MPI_WTIME() - t_qedstep
      #endif
      !.................................................

      !-------------------------------------------------
      ! Pushing particles
        t_movestep = MPI_WTIME()
      call moveParticles()
        t_movestep = MPI_WTIME() - t_movestep
      !.................................................

      !-------------------------------------------------
      ! Advancing 2nd halfstep of `dB / dt = curl E`
        t_fldslvrstep = MPI_WTIME() - t_fldslvrstep
      call advanceBHalfstep()
        t_fldslvrstep = MPI_WTIME() - t_fldslvrstep
      !.................................................

      !-------------------------------------------------
      ! Exchanging `B`-fields
        t_fldexchstep = MPI_WTIME() - t_fldexchstep
      call exchangeFields(.false., .true.)
        t_fldexchstep = MPI_WTIME() - t_fldexchstep
      !.................................................

      !-------------------------------------------------
      ! Advancing fullstep of `dE / dt = -curl B`
        t_fldslvrstep = MPI_WTIME() - t_fldslvrstep
      call advanceEFullstep()
        t_fldslvrstep = MPI_WTIME() - t_fldslvrstep
      !.................................................

      !-------------------------------------------------
      ! Exchanging `E`-fields
        t_fldexchstep = MPI_WTIME() - t_fldexchstep
      call exchangeFields(.true., .false.)
        t_fldexchstep = MPI_WTIME() - t_fldexchstep
      !.................................................

      !-------------------------------------------------
      ! Depositing current: `j_s = rho_s * v_s`
        t_depositstep = MPI_WTIME()
      call depositCurrents()
        t_depositstep = MPI_WTIME() - t_depositstep
      !.................................................

      !-------------------------------------------------
      ! Exchanging currents
        t_fldexchstep = MPI_WTIME() - t_fldexchstep
      call exchangeCurrents()
        t_fldexchstep = MPI_WTIME() - t_fldexchstep
      !.................................................

      !-------------------------------------------------
      ! Filtering currents
        t_filterstep = MPI_WTIME()
      call filterCurrents()
        t_filterstep = MPI_WTIME() - t_filterstep
      !.................................................

      !-------------------------------------------------
      ! Adding currents: `dE / dt += -j`
        t_fldslvrstep = MPI_WTIME() - t_fldslvrstep
      call addCurrents()
        t_fldslvrstep = MPI_WTIME() - t_fldslvrstep
      !.................................................

      !-------------------------------------------------
      ! Exchanging `E`-fields
        t_fldexchstep = MPI_WTIME() - t_fldexchstep
      call exchangeFields(.true., .false.)
        t_fldexchstep = MPI_WTIME() - t_fldexchstep
      !.................................................

      !-------------------------------------------------
      ! Exchanging particles
        t_prtlexchxtep = MPI_WTIME()
      call exchangeParticles()
      call clearGhostParticles()
      call checkTileSizes()
        t_prtlexchxtep = MPI_WTIME() - t_prtlexchxtep
      !.................................................

      !-------------------------------------------------
      ! User defined driving and ...
      !     ... boundary conditions for particles
        t_usrfuncs = MPI_WTIME() - t_usrfuncs
      call userParticleBoundaryConditions(timestep)
      call clearGhostParticles() ! hack
      call userDriveParticles(timestep)
        t_usrfuncs = MPI_WTIME() - t_usrfuncs
      !.................................................

      !-------------------------------------------------
      ! Output
      t_outputstep = 0
      if ((modulo(timestep, output_interval) .eq. 0) .and.&
        & (timestep .ge. output_start)) then
        t_outputstep = MPI_WTIME()
        call writeOutput(timestep)
        t_outputstep = MPI_WTIME() - t_outputstep
      end if
      !.................................................
      ! </ MAINLOOP

        t_fullstep = MPI_WTIME() - t_fullstep

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (ierr .eq. MPI_SUCCESS) then
        call makeReport(timestep)
      end if
    end do
  end subroutine mainloop

  subroutine makeReport(tstep)
    implicit none
    integer, intent(in)           :: tstep
    integer                       :: ierr, s, ti, tj, tk
    real                          :: fullstep
    integer, allocatable          :: nprt_sp(:), nprt_sp_global(:,:)
    real(kind=8), allocatable     :: dt_fullstep(:), dt_movestep(:),&
                                   & dt_depositstep(:), dt_filterstep(:),&
                                   & dt_outputstep(:), dt_fldexchstep(:),&
                                   & dt_prtlexchxtep(:), dt_fldslvrstep(:),&
                                   & dt_usrfuncs(:)

    #ifdef QED
      real(kind=8), allocatable     :: dt_qedstep(:)
    #endif

    ! full # of particles for each species
    allocate(nprt_sp(nspec), nprt_sp_global(nspec, mpi_size))
    nprt_sp(:) = 0
    do s = 1, nspec
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            nprt_sp(s) = nprt_sp(s) + species(s)%prtl_tile(ti, tj, tk)%npart_sp
          end do
        end do
      end do
    end do

    call MPI_GATHER(nprt_sp, nspec, MPI_INTEGER,&
                  & nprt_sp_global, nspec, MPI_INTEGER,&
                  & 0, MPI_COMM_WORLD, ierr)

    allocate(dt_fullstep(mpi_size), dt_movestep(mpi_size))
    allocate(dt_depositstep(mpi_size), dt_filterstep(mpi_size))
    allocate(dt_outputstep(mpi_size), dt_fldexchstep(mpi_size))
    allocate(dt_prtlexchxtep(mpi_size), dt_fldslvrstep(mpi_size))
    allocate(dt_usrfuncs(mpi_size))

    call MPI_GATHER(t_fullstep, 1, MPI_REAL8,&
                  & dt_fullstep, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER(t_movestep, 1, MPI_REAL8,&
                  & dt_movestep, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER(t_depositstep, 1, MPI_REAL8,&
                  & dt_depositstep, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER(t_filterstep, 1, MPI_REAL8,&
                  & dt_filterstep, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER(t_outputstep, 1, MPI_REAL8,&
                  & dt_outputstep, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER(t_fldexchstep, 1, MPI_REAL8,&
                  & dt_fldexchstep, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER(t_prtlexchxtep, 1, MPI_REAL8,&
                  & dt_prtlexchxtep, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER(t_fldslvrstep, 1, MPI_REAL8,&
                  & dt_fldslvrstep, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER(t_usrfuncs, 1, MPI_REAL8,&
                  & dt_usrfuncs, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)

    #ifdef QED
      allocate(dt_qedstep(mpi_size))
      call MPI_GATHER(t_qedstep, 1, MPI_REAL8,&
                    & dt_qedstep, 1, MPI_REAL8,&
                    & 0, MPI_COMM_WORLD, ierr)
    #endif

    if (mpi_rank .eq. 0) then
      fullstep = SUM(dt_fullstep) * 1000 / mpi_size
      call printTimeHeader(tstep)
      ! call printReport(.true., "timestep: " // STR(tstep))
      call printTime(dt_fullstep, "Full_step: ")
      call printTime(dt_movestep, "  move_step: ", fullstep)
      call printTime(dt_depositstep, "  deposit_step: ", fullstep)
      call printTime(dt_filterstep, "  filter_step: ", fullstep)
      call printTime(dt_fldexchstep, "  fld_exchange: ", fullstep)
      call printTime(dt_prtlexchxtep, "  prtl_exchange: ", fullstep)
      call printTime(dt_fldslvrstep, "  fld_solver: ", fullstep)
      call printTime(dt_usrfuncs, "  usr_funcs: ", fullstep)
      call printTime(dt_outputstep, "  output_step: ", fullstep)

      #ifdef QED
        call printTime(dt_qedstep, "  qed_step: ", fullstep)
      #endif

      do s = 1, nspec
        if (s .ne. nspec) then
          call printNpart(nprt_sp_global(s, :),&
                        & "  nprt " // trim(STR(s)) // " [core]: ")
        else
          call printNpart(nprt_sp_global(s, :),&
                        & "  nprt " // trim(STR(s)) // " [core]: ")
        end if
      end do

      call printTimeFooter()

      print *, ""
    end if

    deallocate(dt_fullstep, dt_movestep)
    deallocate(dt_depositstep, dt_filterstep)
    deallocate(dt_outputstep, dt_fldexchstep)
    deallocate(dt_prtlexchxtep, dt_fldslvrstep)
    deallocate(dt_usrfuncs)

    #ifdef QED
      deallocate(dt_qedstep)
    #endif
  end subroutine makeReport

end module m_mainloop
