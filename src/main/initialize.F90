module m_initialize
#ifdef IFPORT
  use ifport, only: makedirqq
#endif
  use m_globalnamespace
  use m_outputnamespace, only: tot_output_index, slice_index, &
                               tot_output_enable, slice_output_enable, hst_enable
  use m_qednamespace

  use m_helpers
  use m_errors

  use m_aux
  use m_readinput, only: getInput, readCommandlineArgs

  use m_domain
  use m_loadbalancing, only: initializeLB
#ifdef SLB
  use m_staticlb, only: redistributeMeshblocksSLB
#endif

  use m_fieldlogistics, only: initializeFields
  use m_particlelogistics, only: initializeParticles
  use m_exchangeparts, only: initializePrtlExchange
  use m_currentdeposit, only: resetCurrents
  use m_exchangefields

  use m_userfile, only: userReadInput, userInitParticles, userInitFields, userFieldBoundaryConditions

#ifdef SLB
  use m_userfile, only: user_slb_load_ptr => userSLBload
#endif

  use m_outputlogistics, only: initializeOutput, initializeSlice
  use m_writehistory, only: initializeHistory

#ifdef USROUTPUT
  use m_writeusroutput, only: initializeUsrOutput
#endif

  use m_restart, only: initializeRestart, restartSimulation, rst_simulation, rst_enable, rst_tlim_enable

  use m_particlebinning

#ifdef DOWNSAMPLING
  use m_particledownsampling
#endif

  ! extra physics
#ifdef RADIATION
  use m_radiation, only: initializeRadiation
  use m_outputnamespace, only: rad_spectra, glob_rad_spectra
#endif

  implicit none

  !--- PRIVATE functions -----------------------------------------!
  private :: initializeCommunications, initializeOutput, &
             initializeDirectories, &
             printParams, initializeSlice, &
             distributeMeshblocks, initializeDomain, &
             initializeSimulation

#ifdef DOWNSAMPLING
  private :: initializeDownsampling
#endif
  !...............................................................!
contains
  ! initialize all the necessary arrays and variables
  subroutine initializeAll()
    implicit none
    call readCommandlineArgs()

    ! initializing the simulation parameters class ...
    ! ... which stores all the input values for the simulation
    call initializeSimulationParameters()
    call initializeWarnings()

    call initializeDomain()

    call initializeCommunications()

    call getInput('particles', 'nspec', nspec, 2)

    call initializeOutput()
    call initializeHistory()
#ifdef USROUTPUT
    call initializeUsrOutput()
#endif
    call initializeSlice()
    call initializeRestart()
    call initializeDirectories()
    call initializeSimulation()

    call distributeMeshblocks()

    call userReadInput()

    call initializeLB()

#ifdef SLB
    call redistributeMeshblocksSLB(user_slb_load_ptr)
#endif

    call initializeFields()

    call initializeParticles()

    call resetCurrents()

#ifdef RADIATION
    call initializeRadiation()
#endif

#ifdef DOWNSAMPLING
    call initializeDownsampling()
#endif

    call initializePrtlExchange()

    call initializeRandomSeed(mpi_rank)

    if (.not. rst_simulation) then
      call userInitFields()
      call userFieldBoundaryConditions(0, updateE=.true., updateB=.true.)
      call exchangeFields(exchangeE=.true., exchangeB=.true.)
      call userInitParticles()
      call printDiag("userInitialize()", 1)
    else
      call restartSimulation()
    end if
    call userReadInput()

    call checkEverything()

    call printParams()

    call printDiag("InitializeAll()", 0)
  end subroutine initializeAll

  subroutine printParams()
    implicit none
    integer :: n
    character(len=STR_MAX) :: FMT
    real :: dummy
    ! printing simulation parameters in the report

    if (mpi_rank .eq. 0) then
      FMT = '== Full simulation parameters =========================================='
      write (*, '(A)') trim(FMT)
      do n = 1, sim_params % count
        if (sim_params % param_type(n) .eq. 1) then
          FMT = '(A30,A1,A20,A1,I10)'
          write (*, FMT) trim(sim_params % param_group(n) % str), ':', &
            trim(sim_params % param_name(n) % str), ':', &
            sim_params % param_value(n) % value_int
        else if (sim_params % param_type(n) .eq. 2) then
          FMT = getFMTForReal(sim_params % param_value(n) % value_real)
          FMT = '(A30,A1,A20,A1,'//trim(FMT)//')'
          write (*, FMT) trim(sim_params % param_group(n) % str), ':', &
            trim(sim_params % param_name(n) % str), ':', &
            sim_params % param_value(n) % value_real
        else if (sim_params % param_type(n) .eq. 3) then
          FMT = '(A30,A1,A20,A1,L10)'
          write (*, FMT) trim(sim_params % param_group(n) % str), ':', &
            trim(sim_params % param_name(n) % str), ':', &
            sim_params % param_value(n) % value_bool
        else
          call throwError('ERROR. Unknown `param_type` in `saveAllParameters`.')
        end if
      end do
      FMT = '........................................................................'
      write (*, '(A)') trim(FMT)

      print *, ""
      FMT = '== Fiducial physical parameters ========================================'
      write (*, '(A)') trim(FMT)

      dummy = c_omp
      FMT = getFMTForReal(dummy)
      FMT = '(A35,'//trim(FMT)//')'
      write (*, FMT) trim('skin depth [dx]:'), dummy

      dummy = 2.0 * M_PI * c_omp / CC
      FMT = getFMTForReal(dummy)
      FMT = '(A35,'//trim(FMT)//')'
      write (*, FMT) trim('plasma oscillation period [dt]:'), dummy

      dummy = c_omp / sqrt(sigma)
      FMT = getFMTForReal(dummy)
      FMT = '(A35,'//trim(FMT)//')'
      write (*, FMT) trim('gyroradius [dx]:'), dummy

      dummy = 2.0 * M_PI * c_omp / (CC * sqrt(sigma))
      FMT = getFMTForReal(dummy)
      FMT = '(A35,'//trim(FMT)//')'
      write (*, FMT) trim('gyration period [dt]:'), dummy

#ifdef RADIATION
#ifdef SYNCHROTRON
      dummy = cool_gamma_syn
      FMT = getFMTForReal(dummy)
      FMT = '(A35,'//trim(FMT)//')'
      write (*, FMT) trim('synchrotron break:'), dummy

      dummy = emit_gamma_syn
      FMT = getFMTForReal(dummy)
      FMT = '(A35,'//trim(FMT)//')'
      write (*, FMT) trim('synchrotron peak is mc^2 for:'), dummy
#endif

#ifdef INVERSECOMPTON
      dummy = cool_gamma_ic
      FMT = getFMTForReal(dummy)
      FMT = '(A35,'//trim(FMT)//')'
      write (*, FMT) trim('inverse Compton break:'), dummy

      dummy = emit_gamma_ic
      FMT = getFMTForReal(dummy)
      FMT = '(A35,'//trim(FMT)//')'
      write (*, FMT) trim('inverse Compton peak is mc^2 for:'), dummy
#endif
#endif

      FMT = '........................................................................'
      write (*, '(A)') trim(FMT)
      print *, ""
    end if
  end subroutine printParams

  subroutine initializeCommunications()
    implicit none
    integer :: ierr
    real :: dummy
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)
    mpi_statsize = MPI_STATUS_SIZE
    if (mpi_size .ne. sizex * sizey * sizez) then
      call throwError('ERROR: # of processors is not equal to the number of processors from input')
    end if
    if (SIZEOF(dummy) .eq. 4) then
      default_mpi_real = MPI_REAL
    else if (SIZEOF(dummy) .eq. 8) then
      default_mpi_real = MPI_REAL8
    else
      call throwError('ERROR: SIZEOF(dummy) is not 4 or 8')
    end if
  end subroutine initializeCommunications

  subroutine initializeDomain()
    implicit none
    sizex = 1; sizey = 1; sizez = 1
    global_mesh % x0 = 0; global_mesh % y0 = 0; global_mesh % z0 = 0
    global_mesh % sx = 1; global_mesh % sy = 1; global_mesh % sz = 1

#if defined(oneD) || defined (twoD) || defined (threeD)
    call getInput('node_configuration', 'sizex', sizex)
    call getInput('grid', 'mx0', global_mesh % sx)
#endif
#if defined(twoD) || defined (threeD)
    call getInput('node_configuration', 'sizey', sizey)
    call getInput('grid', 'my0', global_mesh % sy)
#endif
#if defined(threeD)
    call getInput('node_configuration', 'sizez', sizez)
    call getInput('grid', 'mz0', global_mesh % sz)
#endif

    if ((modulo(global_mesh % sx, sizex) .ne. 0) .or. &
        (modulo(global_mesh % sy, sizey) .ne. 0) .or. &
        (modulo(global_mesh % sz, sizez) .ne. 0)) then
      call throwError('ERROR: grid size is not evenly divisible by the number of cores')
    end if

    call getInput('grid', 'boundary_x', boundary_x, 1)
    call getInput('grid', 'boundary_y', boundary_y, 1)
    call getInput('grid', 'boundary_z', boundary_z, 1)
#ifdef oneD
    boundary_y = 1
    boundary_z = 1
#elif defined(twoD)
    boundary_z = 1
#endif

    call getInput('grid', 'abs_thick', ds_abs, 10.0)
    call getInput('grid', 'absorb_x', absorb_x, 0)
    call getInput('grid', 'absorb_y', absorb_y, 0)
    call getInput('grid', 'absorb_z', absorb_z, 0)
#ifndef ABSORB
    if (absorb_x .ne. 0 .or. absorb_y .ne. 0 .or. absorb_z .ne. 0) then
      call throwError('ERROR: `absorb_x`, `absorb_y`, and `absorb_z` are only valid if ABSORB is defined')
    end if
#endif

    if ((absorb_x .eq. 2) .or. (absorb_y .eq. 2) .or. (absorb_z .eq. 2)) then
      absorb_x = 2
      absorb_y = 2
      absorb_z = 2
    end if
  end subroutine initializeDomain

  subroutine distributeMeshblocks()
    implicit none
    integer, dimension(3) :: ind, m
    integer :: rnk
    m(1) = global_mesh % sx / sizex
    m(2) = global_mesh % sy / sizey
    m(3) = global_mesh % sz / sizez
    allocate (meshblocks(mpi_size))
    if (.not. allocated(new_meshblocks)) allocate (new_meshblocks(mpi_size))
    this_meshblock % ptr => meshblocks(mpi_rank + 1)
    do rnk = 0, mpi_size - 1
      ind = rnkToInd(rnk)
      meshblocks(rnk + 1) % rnk = rnk
      ! find sizes and corner coords
      meshblocks(rnk + 1) % sx = m(1)
      meshblocks(rnk + 1) % sy = m(2)
      meshblocks(rnk + 1) % sz = m(3)
      meshblocks(rnk + 1) % x0 = ind(1) * m(1) + global_mesh % x0
      meshblocks(rnk + 1) % y0 = ind(2) * m(2) + global_mesh % y0
      meshblocks(rnk + 1) % z0 = ind(3) * m(3) + global_mesh % z0
    end do
    ! assign all neighbors
    call reassignNeighborsForAll(meshblocks)

    call printDiag("distributeMeshblocks()", 1)
  end subroutine distributeMeshblocks

  subroutine initializeSimulation()
    implicit none
    call getInput('time', 'last', final_timestep, 1000)
    call getInput('time', 'wall_t_max', wall_t_max, real(0,8))
    call getInput('algorithm', 'nfilter', nfilter, 16)
    call getInput('algorithm', 'c', CC, 0.45)
    call getInput('algorithm', 'corr', CORR, 1.025)
    call getInput('algorithm', 'fieldsolver', enable_fieldsolver, .true.)
    call getInput('algorithm', 'currdeposit', enable_currentdeposit, .true.)
    call getInput('plasma', 'ppc0', ppc0)
    call getInput('plasma', 'sigma', sigma)
    if (sigma .le. 0.0) then
      call throwError('Reference sigma value must be > 0.')
    end if
    call getInput('plasma', 'c_omp', c_omp)
    call renormalizeUnits()

    call getInput('grid', 'resize_tiles', resize_tiles, .false.)
    call getInput('grid', 'shrink_tiles', shrink_tiles, .false.)
    call getInput('grid', 'min_tile_nprt', min_tile_nprt, 100)

    call printDiag("initializeSimulation()", 1)
  end subroutine initializeSimulation

  subroutine initializeDirectories()
    ! create output/restart directories
    !   if does not already exist
    !     note: some compilers may not support IFPORT
#ifdef IFPORT
    logical :: result
#endif
    diag_file_name = trim(output_dir_name)//'/'//trim(diag_file_name)
    warn_file_name = trim(output_dir_name)//'/'//trim(warn_file_name)
    if (mpi_rank .eq. 0) then
#ifdef IFPORT
      result = makedirqq(trim(output_dir_name))
      output_dir_spec = trim(output_dir_name)//'/spec'
      output_dir_flds = trim(output_dir_name)//'/flds'
      output_dir_prtl = trim(output_dir_name)//'/prtl'
      result = makedirqq(trim(output_dir_spec))
      result = makedirqq(trim(output_dir_flds))
      result = makedirqq(trim(output_dir_prtl))
      if (rst_enable .or. rst_tlim_enable) then
        result = makedirqq(trim(restart_dir_name))
      end if
      if (slice_output_enable) then
        result = makedirqq(trim(slice_dir_name))
      end if
#else
      call system('mkdir -p '//trim(output_dir_name))
      output_dir_spec = trim(output_dir_name)//'/spec'
      output_dir_flds = trim(output_dir_name)//'/flds'
      output_dir_prtl = trim(output_dir_name)//'/prtl'
      call system('mkdir -p '//trim(output_dir_spec))
      call system('mkdir -p '//trim(output_dir_flds))
      call system('mkdir -p '//trim(output_dir_prtl))
      if (rst_enable .or. rst_tlim_enable) then
        call system('mkdir -p '//trim(restart_dir_name))
      end if
      if (slice_output_enable) then
        call system('mkdir -p '//trim(slice_dir_name))
      end if
#endif
      ! diagnostics file
      open (UNIT_diag, file=diag_file_name, status="replace", form="formatted")
      close (UNIT_diag)
      ! warnings file
      open (UNIT_warn, file=warn_file_name, status="replace", form="formatted")
      close (UNIT_warn)
    end if
  end subroutine initializeDirectories

#ifdef DOWNSAMPLING
  subroutine initializeDownsampling()
    implicit none
    call getInput('downsampling', 'interval', dwn_interval, 1)
    call getInput('downsampling', 'start', dwn_start, 0)

    call getInput('downsampling', 'max_weight', dwn_maxweight, 1e2)

    call getInput('downsampling', 'cartesian_bins', dwn_cartesian_bins, .false.)

    call getInput('downsampling', 'dynamic_bins', dwn_dynamic_bins, .false.)
    if (dwn_cartesian_bins) then
      call getInput('downsampling', 'mom_bins', dwn_n_mom_bins, 5)
      if (dwn_dynamic_bins) then
        call getInput('downsampling', 'mom_spread', dwn_mom_spread, 0.1)
      end if
    else
      call getInput('downsampling', 'angular_bins', dwn_n_angular_bins, 5)
      call getInput('downsampling', 'energy_bins', dwn_n_energy_bins, 5)
      call getInput('downsampling', 'log_e_bins', dwn_log_e_bins, .true.)
    end if
    call getInput('downsampling', 'energy_max', dwn_energy_max, 1e2)
    call getInput('downsampling', 'energy_min', dwn_energy_min, 1e-2)
    call getInput('downsampling', 'int_weights', dwn_int_weights, .false.)

    call printDiag("initializeDownsampling()", 1)
  end subroutine initializeDownsampling
#endif
end module m_initialize
