#include "../defs.F90"

module m_mainloop
  use m_globalnamespace
  use m_outputnamespace, only: slice_output_enable, slice_output_start, slice_output_interval,&
                             & tot_output_enable, tot_output_start, tot_output_interval,&
                             & hst_enable, hst_interval, usrout_enable, usrout_interval
  use m_helpers
  use m_aux
  use m_writeslice, only: writeSlices
  use m_writetot, only: writeTotOutput
  use m_writehistory, only: writeHistory
  use m_restart, only: writeRestart, rst_enable, rst_interval, rst_start
  use m_fldsolver
  use m_mover
  use m_currentdeposit
  use m_exchangeparts
  use m_exchangefields
  use m_exchangecurrents
  use m_particlelogistics
  use m_filtering
  use m_userfile, only: userDriveParticles, userParticleBoundaryConditions,&
                      & userFieldBoundaryConditions, userCurrentDeposit
  #ifdef USROUTPUT
    use m_userfile, only: userOutput
  #endif
  use m_errors

  #ifdef DOWNSAMPLING
    use m_particledownsampling
  #endif

  #ifdef ALB
    use m_adaptivelb, only: redistributeMeshblocksALB
  #endif

  implicit none

  !--- PRIVATE functions -----------------------------------------!
  private :: makeReport, startTimer, flushTimer
  !...............................................................!

  !--- PRIVATE variables -----------------------------------------!
  integer, private       :: timestep
  real(kind=8), private  :: timers(20), d_timers(20)
  !...............................................................!
contains
  subroutine mainloop()
    implicit none
    integer       :: ierr, i
    integer       :: s, ti, tj, tk, p

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call printDiag("Starting mainloop()", 0)

    ! timer numbering ([*] = optional):
    !  1 = full step
    !  2 = move substep
    !  3 = deposit substep
    !  4 = filtering substep 
    !  5 = output
    !  6 = field exchange/communications
    !  7 = particle exchange/communications
    !  8 = field solver
    !  9 = user-specific routines
    ! 10 = ...
    ! 11 = particle downsampling substep [*]
    ! 12 = adaptive load balancing substep [*]

    do timestep = start_timestep, final_timestep
      call printDiag("", 20)
      call printDiag("Starting timestep # "//STR(timestep), 0)
        timers(:) = 0d0
        d_timers(:) = 0d0
        call startTimer(1)

      ! MAINLOOP >

      !-------------------------------------------------
      ! Dynamic balancing of processor loads
      #ifdef ALB
          call startTimer(12)
        call redistributeMeshblocksALB(timestep)
        call exchangeFields(exchangeE=.true., exchangeB=.true.)
          call flushTimer(12)
      #endif

      !-------------------------------------------------
      ! User defined boundary conditions for B-field
        call startTimer(9)
      call userFieldBoundaryConditions(timestep, updateE=.true., updateB=.true.)
        call flushTimer(9)
      !.................................................

      !-------------------------------------------------
      ! Reset currents to zero
        call startTimer(3)
      if (enable_currentdeposit) then
        call resetCurrents()
      end if
        call flushTimer(3)
      !.................................................

      !-------------------------------------------------
      ! Exchanging `E` and `B`-fields
        call startTimer(6)
      call exchangeFields(exchangeE=.true., exchangeB=.true.)
        call flushTimer(6)
      !.................................................

      !-------------------------------------------------
      ! Advancing 1st halfstep of `dB / dt = curl E`
        call startTimer(8)
      if (enable_fieldsolver) call advanceBHalfstep()
        call flushTimer(8)
      !.................................................

      !-------------------------------------------------
      ! User defined boundary conditions for B-field
        call startTimer(9)
      if (enable_fieldsolver) call userFieldBoundaryConditions(timestep, updateE=.false., updateB=.true.)
        call flushTimer(9)
      !.................................................

      !-------------------------------------------------
      ! Exchanging `B`-fields
        call startTimer(6)
      if (enable_fieldsolver) call exchangeFields(exchangeE=.false., exchangeB=.true.)
        call flushTimer(6)
      !.................................................

      !-------------------------------------------------
      ! Pushing particles
        call startTimer(2)
      call moveParticles(timestep)
        call flushTimer(2)
      !.................................................

      !-------------------------------------------------
      ! User defined driving for particles
        call startTimer(9)
      call userDriveParticles(timestep)
        call flushTimer(9)
      !.................................................

      !-------------------------------------------------
      ! Advancing 2nd halfstep of `dB / dt = curl E`
        call startTimer(8)
      if (enable_fieldsolver) call advanceBHalfstep()
        call flushTimer(8)
      !.................................................

      !-------------------------------------------------
      ! User defined boundary conditions for B-field
        call startTimer(9)
      call userFieldBoundaryConditions(timestep, updateE=.false., updateB=.true.)
        call flushTimer(9)
      !.................................................

      !-------------------------------------------------
      ! Exchanging `B`-fields
        call startTimer(6)
      call exchangeFields(exchangeE=.false., exchangeB=.true.)
        call flushTimer(6)
      !.................................................

      !-------------------------------------------------
      ! Advancing fullstep of `dE / dt = -curl B`
        call startTimer(8)
      if (enable_fieldsolver) call advanceEFullstep()
        call flushTimer(8)
      !.................................................

      !-------------------------------------------------
      ! Exchanging `E`-fields
        call startTimer(6)
      if (enable_fieldsolver) call exchangeFields(exchangeE=.true., exchangeB=.false.)
        call flushTimer(6)
      !.................................................

      !-------------------------------------------------
      ! Depositing current: `j_s = rho_s * v_s`
        call startTimer(3)
      if (enable_currentdeposit) then
        call depositCurrents()
      end if
        call flushTimer(3)
      !.................................................

      !-------------------------------------------------
      ! Additional user-specific current deposition routine
        call startTimer(9)
      call userCurrentDeposit(timestep)
        call flushTimer(9)
      !.................................................

      !-------------------------------------------------
      ! Particle downsampling
      #ifdef DOWNSAMPLING
          call startTimer(1)
        call exchangeParticles()
        call clearGhostParticles()
        call checkTileSizes()
        call downsamplingStep(timestep)
        call clearGhostParticles()
          call flushTimer(11)
      #endif
      !.................................................

      !-------------------------------------------------
      ! Exchanging currents
        call startTimer(6)
      if (enable_currentdeposit) call exchangeCurrents()
        call flushTimer(6)
      !.................................................

      !-------------------------------------------------
      ! Filtering currents
        call startTimer(4)
      if (enable_currentdeposit) call filterCurrents()
        call flushTimer(4)
      !.................................................

      !-------------------------------------------------
      ! Adding currents: `dE / dt += -j`
        call startTimer(8)
      if (enable_fieldsolver) call addCurrents()
        call flushTimer(8)
      !.................................................

      !-------------------------------------------------
      ! User defined boundary conditions for E-field
        call startTimer(9)
      call userFieldBoundaryConditions(timestep, updateE=.true., updateB=.false.)
        call flushTimer(9)
      !.................................................

      !-------------------------------------------------
      ! Exchanging `E`-fields
        call startTimer(6)
      call exchangeFields(exchangeE=.true., exchangeB=.false.)
        call flushTimer(6)
      !.................................................

      !-------------------------------------------------
      ! Exchanging particles
        call startTimer(7)
      call exchangeParticles()
      call clearGhostParticles()
      call checkTileSizes()
        call flushTimer(7)
      !.................................................

      !-------------------------------------------------
      ! User defined boundary conditions for particles
        call startTimer(9)
      call userParticleBoundaryConditions(timestep)
      call clearGhostParticles()
        call flushTimer(9)
      !.................................................

      !-------------------------------------------------
      ! Tot output
      if ((tot_output_enable) .and.&
        & (modulo(timestep, tot_output_interval) .eq. 0) .and.&
        & (timestep .ge. tot_output_start)) then
          call startTimer(5)
        call writeTotOutput(timestep)
          call flushTimer(5)
      end if
      !.................................................

      !-------------------------------------------------
      ! History
      if ((hst_enable) .and.&
        & (modulo(timestep, hst_interval) .eq. 0)) then
          call startTimer(5)
        call writeHistory(timestep)
          call flushTimer(5)
      end if
      !.................................................

      !-------------------------------------------------
      ! User-defined output
      #ifdef USROUTPUT
        if ((usrout_enable) .and.&
          & (modulo(timestep, usrout_interval) .eq. 0)) then
            call startTimer(9)
          call userOutput(timestep)
            call flushTimer(9)
        end if
      #endif
      !.................................................

      !-------------------------------------------------
      ! Slices
      if ((slice_output_enable) .and.&
        & (timestep .ge. slice_output_start) .and.&
        & (modulo(timestep, slice_output_interval) .eq. 0)) then
          call startTimer(5)
        call writeSlices(timestep)
          call flushTimer(5)
      end if
      !.................................................

      !-------------------------------------------------
      ! Restart
      if ((rst_enable) .and.&
        & (timestep .ge. rst_start) .and.&
        & (modulo(timestep - rst_start, rst_interval) .eq. 0) .and.&
        & (timestep .gt. 0)) then
          call startTimer(5)
        call writeRestart(timestep)
          call flushTimer(5)
      end if
      !.................................................
      ! </ MAINLOOP

        call flushTimer(1)

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (ierr .eq. MPI_SUCCESS) then
        call makeReport(timestep)
      end if

      call printDiag("Finished timestep # "//STR(timestep), 0)
    end do
  end subroutine mainloop

  subroutine startTimer(timer_id)
    implicit none
    integer, intent(in) :: timer_id
    integer :: ierr
    d_timers(timer_id) = MPI_WTIME()
  end subroutine startTimer

  subroutine flushTimer(timer_id)
    implicit none
    integer, intent(in) :: timer_id
    real(kind=8)        :: dt
    integer :: ierr
    dt = (MPI_WTIME() - d_timers(timer_id))
    timers(timer_id) = timers(timer_id) + dt
  end subroutine flushTimer

  subroutine makeReport(tstep)
    implicit none
    integer, intent(in)           :: tstep
    integer                       :: ierr, s, ti, tj, tk
    real                          :: fullstep
    integer(kind=8), allocatable  :: nprt_sp(:), nprt_sp_global(:,:)

    real(kind=8)  :: t_fullstep, t_movestep,&
                   & t_depositstep, t_filterstep,&
                   & t_outputstep, t_fldexchstep,&
                   & t_prtlexchxtep, t_fldslvrstep,&
                   & t_usrfuncs
    #ifdef DOWNSAMPLING
      real(kind=8) :: t_dwnstep
    #endif

    #ifdef ALB
      real(kind=8) :: t_albstep
    #endif

    real(kind=8), allocatable     :: dt_fullstep(:), dt_movestep(:),&
                                   & dt_depositstep(:), dt_filterstep(:),&
                                   & dt_outputstep(:), dt_fldexchstep(:),&
                                   & dt_prtlexchxtep(:), dt_fldslvrstep(:),&
                                   & dt_usrfuncs(:)

    #ifdef DOWNSAMPLING
      real(kind=8), allocatable     :: dt_dwnstep(:)
    #endif
    
    #ifdef ALB
      real(kind=8), allocatable     :: dt_albstep(:)
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

    call MPI_GATHER(nprt_sp, nspec, MPI_INTEGER8,&
                  & nprt_sp_global, nspec, MPI_INTEGER8,&
                  & 0, MPI_COMM_WORLD, ierr)

    ! accumulate timers
    t_fullstep = timers(1);     t_movestep = timers(2)
    t_depositstep = timers(3);  t_filterstep = timers(4)
    t_outputstep = timers(5);   t_fldexchstep = timers(6)
    t_prtlexchxtep = timers(7); t_fldslvrstep = timers(8)
    t_usrfuncs = timers(9)

    #ifdef DOWNSAMPLING
      t_dwnstep = timers(11)
    #endif

    #ifdef ALB
      t_albstep = timers(12)
    #endif

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

    #ifdef DOWNSAMPLING
      allocate(dt_dwnstep(mpi_size))
      call MPI_GATHER(t_dwnstep, 1, MPI_REAL8,&
                    & dt_dwnstep, 1, MPI_REAL8,&
                    & 0, MPI_COMM_WORLD, ierr)
    #endif
    
    #ifdef ALB
      allocate(dt_albstep(mpi_size))
      call MPI_GATHER(t_albstep, 1, MPI_REAL8,&
                    & dt_albstep, 1, MPI_REAL8,&
                    & 0, MPI_COMM_WORLD, ierr)
    #endif

    if (mpi_rank .eq. 0) then
      fullstep = SUM(dt_fullstep) * 1000 / mpi_size
      call printTimeHeader(tstep)

      call printTime(dt_fullstep, "Full_step: ")
      call printTime(dt_movestep, "  move_step: ", fullstep)
      call printTime(dt_depositstep, "  deposit_step: ", fullstep)
      call printTime(dt_filterstep, "  filter_step: ", fullstep)
      call printTime(dt_fldexchstep, "  fld_exchange: ", fullstep)
      call printTime(dt_prtlexchxtep, "  prtl_exchange: ", fullstep)
      call printTime(dt_fldslvrstep, "  fld_solver: ", fullstep)
      call printTime(dt_usrfuncs, "  usr_funcs: ", fullstep)
      call printTime(dt_outputstep, "  output_step: ", fullstep)

      #ifdef DOWNSAMPLING
        call printTime(dt_dwnstep, "  dwn_step: ", fullstep)
      #endif

      #ifdef ALB
        call printTime(dt_albstep, "  alb_step: ", fullstep)
      #endif

      call printNpartHeader()
      do s = 1, nspec
        call printNpart(nprt_sp_global(s, :),&
                      & "  species # " // trim(STR(s)))
      end do

      call printTimeFooter()
      print *, ""
    end if

    call printWarnings(tstep)

    deallocate(dt_fullstep, dt_movestep)
    deallocate(dt_depositstep, dt_filterstep)
    deallocate(dt_outputstep, dt_fldexchstep)
    deallocate(dt_prtlexchxtep, dt_fldslvrstep)
    deallocate(dt_usrfuncs)

    #ifdef DOWNSAMPLING
      deallocate(dt_dwnstep)
    #endif

    #ifdef ALB
      deallocate(dt_albstep)
    #endif

  end subroutine makeReport

end module m_mainloop
