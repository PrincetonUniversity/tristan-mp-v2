#include "../defs.F90"

module m_initialize
  #ifdef IFPORT
    use ifport, only: makedirqq
  #endif
  use m_globalnamespace
  use m_outputnamespace, only: tot_output_index, slice_index,&
                             & tot_output_enable, slice_output_enable, hst_enable
  use m_writerestart, only: rst_simulation, rst_enable
  use m_aux
  use m_readinput
  use m_domain
  use m_loadbalancing
  use m_particles
  use m_particlelogistics
  use m_fields
  use m_userfile, only: userReadInput, userInitParticles,&
                      & userInitFields, user_slb_load_ptr => userSLBload
  use m_outputlogistics, only: initializeOutput, initializeSlice
  use m_writehistory, only: initializeHistory
  use m_writerestart, only: initializeRestart
  use m_helpers
  use m_errors
  use m_particlebinning

  #ifdef DOWNSAMPLING
    use m_particledownsampling
  #endif
  implicit none

  !--- PRIVATE functions -----------------------------------------!
  private :: initializeCommunications, initializeOutput,&
           & firstRankInitialize, initializeParticles,&
           & initializeLB, printParams, initializeSlice,&
           & distributeMeshblocks, initializeDomain,&
           & initializePrtlExchange, initializeFields,&
           & initializeSimulation, checkEverything,&
           & restartSimulation

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
    sim_params%count = 0
    allocate(sim_params%param_type(1000))
    allocate(sim_params%param_group(1000))
    allocate(sim_params%param_name(1000))
    allocate(sim_params%param_value(1000))

    call initializeDomain()

    call initializeCommunications()
      call printDiag((mpi_rank .eq. 0), "initializeCommunications()", .true.)

    call distributeMeshblocks()
      call printDiag((mpi_rank .eq. 0), "distributeMeshblocks()", .true.)

    call initializeLB()
      call printDiag((mpi_rank .eq. 0), "initializeLB()", .true.)

    #ifdef SLB
      call redistributeMeshblocksSLB(user_slb_load_ptr)
        call printDiag((mpi_rank .eq. 0), "redistributeMeshblocksSLB()", .true.)
    #endif

    call initializeSimulation()
      call printDiag((mpi_rank .eq. 0), "initializeSimulation()", .true.)

    call initializeOutput()
      call printDiag((mpi_rank .eq. 0), "initializeOutput()", .true.)
    call initializeHistory()
      call printDiag((mpi_rank .eq. 0), "initializeHistory()", .true.)
    call initializeSlice()
      call printDiag((mpi_rank .eq. 0), "initializeSlice()", .true.)
    call initializeRestart()
      call printDiag((mpi_rank .eq. 0), "initializeRestart()", .true.)

    call initializeFields()
      call printDiag((mpi_rank .eq. 0), "initializeFields()", .true.)

    call initializeParticles()
      call printDiag((mpi_rank .eq. 0), "initializeParticles()", .true.)

    #ifdef DOWNSAMPLING
      call initializeDownsampling()
        call printDiag((mpi_rank .eq. 0), "initializeDownsampling()", .true.)
    #endif

    call initializePrtlExchange()
      call printDiag((mpi_rank .eq. 0), "initializePrtlExchange()", .true.)

    call initializeRandomSeed(mpi_rank)
      call printDiag((mpi_rank .eq. 0), "initializeRandomSeed()", .true.)

    if (mpi_rank .eq. 0) then
      call firstRankInitialize()
      call printDiag(.true., "firstRankInitialize()", .true.)
    end if

    if (.not. rst_simulation) then
      call userReadInput()
      call userInitParticles()
      call userInitFields()
        call printDiag((mpi_rank .eq. 0), "userInitialize()", .true.)
    else
      call userReadInput()
      call restartSimulation()
        call printDiag((mpi_rank .eq. 0), "restartSimulation()", .true.)
    end if

    call checkEverything()
      call printDiag((mpi_rank .eq. 0), "checkEverything()", .true.)

    call printParams()

    call printReport((mpi_rank .eq. 0), "InitializeAll()")
  end subroutine initializeAll

  subroutine printParams()
    implicit none
    integer                 :: n
    character(len=STR_MAX)  :: FMT
    real                    :: dummy
    ! printing simulation parameters in the report

    if (mpi_rank .eq. 0) then
      FMT = '== Full simulation parameters =========================================='
      write(*, '(A)') trim(FMT)
      do n = 1, sim_params%count
        if (sim_params%param_type(n) .eq. 1) then
          FMT = '(A30,A1,A20,A1,I10)'
          write (*, FMT) trim(sim_params%param_group(n)%str), ':',&
                       & trim(sim_params%param_name(n)%str), ':',&
                       & sim_params%param_value(n)%value_int
        else if (sim_params%param_type(n) .eq. 2) then
          FMT = getFMTForReal(sim_params%param_value(n)%value_real)
          FMT = '(A30,A1,A20,A1,' // trim(FMT) // ')'
          write (*, FMT) trim(sim_params%param_group(n)%str), ':',&
                       & trim(sim_params%param_name(n)%str), ':',&
                       & sim_params%param_value(n)%value_real
        else if (sim_params%param_type(n) .eq. 3) then
          FMT = '(A30,A1,A20,A1,L10)'
          write (*, FMT) trim(sim_params%param_group(n)%str), ':',&
                       & trim(sim_params%param_name(n)%str), ':',&
                       & sim_params%param_value(n)%value_bool
        else
          call throwError('ERROR. Unknown `param_type` in `saveAllParameters`.')
        end if
      end do
      FMT = '........................................................................'
      write(*, '(A)') trim(FMT)

      print *, ""
      FMT = '== Fiducial physical parameters ========================================'
      write(*, '(A)') trim(FMT)

      dummy = c_omp
      FMT = getFMTForReal(dummy)
      FMT = '(A35,' // trim(FMT) // ')'
      write (*, FMT) trim('skin depth [dx]:'), dummy

      dummy = 2.0 * M_PI * c_omp / CC
      FMT = getFMTForReal(dummy)
      FMT = '(A35,' // trim(FMT) // ')'
      write (*, FMT) trim('plasma oscillation period [dt]:'), dummy

      dummy = c_omp / sqrt(sigma)
      FMT = getFMTForReal(dummy)
      FMT = '(A35,' // trim(FMT) // ')'
      write (*, FMT) trim('gyroradius [dx]:'), dummy

      dummy = 2.0 * M_PI * c_omp / (CC * sqrt(sigma))
      FMT = getFMTForReal(dummy)
      FMT = '(A35,' // trim(FMT) // ')'
      write (*, FMT) trim('gyration period [dt]:'), dummy

      FMT = '........................................................................'
      write(*, '(A)') trim(FMT)
      print *, ""
    end if
  end subroutine printParams

  subroutine initializeCommunications()
    implicit none
    integer :: ierr
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)
    mpi_statsize = MPI_STATUS_SIZE
    if (mpi_size .ne. sizex * sizey * sizez) then
      call throwError('ERROR: # of processors is not equal to the number of processors from input')
    end if
  end subroutine initializeCommunications

  subroutine initializeDomain()
    implicit none
    sizex = 1; sizey = 1; sizez = 1
    global_mesh%x0 = 0; global_mesh%y0 = 0; global_mesh%z0 = 0
    global_mesh%sx = 1; global_mesh%sy = 1; global_mesh%sz = 1

    #if defined(oneD) || defined (twoD) || defined (threeD)
      call getInput('node_configuration', 'sizex', sizex)
      call getInput('grid', 'mx0', global_mesh%sx)
    #endif

    #if defined(twoD) || defined (threeD)
      call getInput('node_configuration', 'sizey', sizey)
      call getInput('grid', 'my0', global_mesh%sy)
    #endif
    
    #if defined(threeD)
      call getInput('node_configuration', 'sizez', sizez)
      call getInput('grid', 'mz0', global_mesh%sz)
    #endif

    if ((modulo(global_mesh%sx, sizex) .ne. 0) .or.&
      & (modulo(global_mesh%sy, sizey) .ne. 0) .or.&
      & (modulo(global_mesh%sz, sizez) .ne. 0)) then
      call throwError('ERROR: grid size is not evenly divisible by the number of cores')
    end if

    call getInput('grid', 'abs_thick', ds_abs, 10.0)
    call getInput('grid', 'boundary_x', boundary_x, 1)
    call getInput('grid', 'boundary_y', boundary_y, 1)
    call getInput('grid', 'boundary_z', boundary_z, 1)
    #ifdef oneD
      boundary_y = 1
      boundary_z = 1
      if (boundary_x .eq. 2) then
        #ifndef ABSORB
          call throwError('ERROR. define `-DABSORB` flag during compilation for absorbing boundaries.')
        #endif
      end if
    #elif twoD
      boundary_z = 1
      if ((boundary_x .eq. 2) .or. (boundary_y .eq. 2)) then
        boundary_x = 2
        boundary_y = 2
        #ifndef ABSORB
          call throwError('ERROR. define `-DABSORB` flag during compilation for absorbing boundaries.')
        #endif
      end if
    #elif threeD
      if ((boundary_x .eq. 2) .or. (boundary_y .eq. 2) .or. (boundary_z .eq. 2)) then
        boundary_x = 2
        boundary_y = 2
        boundary_z = 2
        #ifndef ABSORB
          call throwError('ERROR. define `-DABSORB` flag during compilation for absorbing boundaries.')
        #endif
      end if
    #endif
  end subroutine initializeDomain

  subroutine distributeMeshblocks()
    implicit none
    integer, dimension(3) :: ind, m
    integer               :: rnk
    m(1) = global_mesh%sx / sizex
    m(2) = global_mesh%sy / sizey
    m(3) = global_mesh%sz / sizez
    allocate(meshblocks(mpi_size))
    if (.not. allocated(new_meshblocks)) allocate(new_meshblocks(mpi_size))
    this_meshblock%ptr => meshblocks(mpi_rank + 1)
    do rnk = 0, mpi_size - 1
      ind = rnkToInd(rnk)
      meshblocks(rnk + 1)%rnk = rnk
      ! find sizes and corner coords
      meshblocks(rnk + 1)%sx = m(1)
      meshblocks(rnk + 1)%sy = m(2)
      meshblocks(rnk + 1)%sz = m(3)
      meshblocks(rnk + 1)%x0 = ind(1) * m(1) + global_mesh%x0
      meshblocks(rnk + 1)%y0 = ind(2) * m(2) + global_mesh%y0
      meshblocks(rnk + 1)%z0 = ind(3) * m(3) + global_mesh%z0
    end do
    ! assign all neighbors
    call reassignNeighborsForAll()
  end subroutine distributeMeshblocks

  subroutine initializeLB()
    implicit none
    ! initializing static LB variables
    slb_x = .false.; slb_sxmin = -1
    slb_y = .false.; slb_symin = -1
    slb_z = .false.; slb_szmin = -1

    alb_x = .false.; alb_sxmin = -1; alb_int_x = -1; alb_start_x = -1
    alb_y = .false.; alb_symin = -1; alb_int_y = -1; alb_start_y = -1
    alb_z = .false.; alb_szmin = -1; alb_int_z = -1; alb_start_z = -1
    #if defined(oneD) || defined (twoD) || defined (threeD)
      call getInput('static_load_balancing', 'in_x', slb_x, .false.)
      call getInput('static_load_balancing', 'sx_min', slb_sxmin, 10)

      call getInput('adaptive_load_balancing', 'in_x', alb_x, .false.)
      call getInput('adaptive_load_balancing', 'sx_min', alb_sxmin, 10)
      call getInput('adaptive_load_balancing', 'interval_x', alb_int_x, 1000)
      call getInput('adaptive_load_balancing', 'start_x', alb_start_x, 0)
    #endif

    #if defined(twoD) || defined (threeD)
      call getInput('static_load_balancing', 'in_y', slb_y, .false.)
      call getInput('static_load_balancing', 'sy_min', slb_symin, 10)

      call getInput('adaptive_load_balancing', 'in_y', alb_y, .false.)
      call getInput('adaptive_load_balancing', 'sy_min', alb_symin, 10)
      call getInput('adaptive_load_balancing', 'interval_y', alb_int_y, 1000)
      call getInput('adaptive_load_balancing', 'start_y', alb_start_y, 0)
    #endif

    #if defined(threeD)
      call getInput('static_load_balancing', 'in_z', slb_z, .false.)
      call getInput('static_load_balancing', 'sz_min', slb_szmin, 10)

      call getInput('adaptive_load_balancing', 'in_z', alb_z, .false.)
      call getInput('adaptive_load_balancing', 'sz_min', alb_szmin, 10)
      call getInput('adaptive_load_balancing', 'interval_z', alb_int_z, 1000)
      call getInput('adaptive_load_balancing', 'start_z', alb_start_z, 0)
    #endif
  end subroutine initializeLB

  subroutine initializeSimulation()
    implicit none
    call getInput('time', 'last', final_timestep, 1000)
    call getInput('algorithm', 'nfilter', nfilter, 16)
    call getInput('algorithm', 'c', CC, 0.45)
    call getInput('algorithm', 'corr', CORR, 1.025)
    call getInput('algorithm', 'fieldsolver', enable_fieldsolver, .true.)
    call getInput('algorithm', 'currdeposit', enable_currentdeposit, .true.)
    call getInput('plasma', 'ppc0', ppc0)
    call getInput('plasma', 'sigma', sigma)
    if (sigma .le. 0.0) then
      call throwError('Reference sigma value must be > 0.')
    endif
    call getInput('plasma', 'c_omp', c_omp)
    call renormalizeUnits()

    call getInput('grid', 'resize_tiles', resize_tiles, .false.)
    call getInput('grid', 'min_tile_nprt', min_tile_nprt, 100)
  end subroutine initializeSimulation

  subroutine initializeParticles()
    implicit none
    integer                 :: s, ti, tj, tk
    character(len=STR_MAX)  :: var_name
    integer                 :: maxptl_

    call getInput('particles', 'nspec', nspec, 2)

    allocate(species(nspec))
    do s = 1, nspec
      #ifdef oneD
        call getInput('grid', 'tileX', species(s)%tile_sx)
        species(s)%tile_sy = 1
        species(s)%tile_sz = 1
      #elif twoD
        call getInput('grid', 'tileX', species(s)%tile_sx)
        call getInput('grid', 'tileY', species(s)%tile_sy)
        species(s)%tile_sz = 1
      #elif threeD
        call getInput('grid', 'tileX', species(s)%tile_sx)
        call getInput('grid', 'tileY', species(s)%tile_sy)
        call getInput('grid', 'tileZ', species(s)%tile_sz)
      #endif
      species(s)%tile_nx = ceiling(real(this_meshblock%ptr%sx) / real(species(s)%tile_sx))
      species(s)%tile_ny = ceiling(real(this_meshblock%ptr%sy) / real(species(s)%tile_sy))
      species(s)%tile_nz = ceiling(real(this_meshblock%ptr%sz) / real(species(s)%tile_sz))
      allocate(species(s)%prtl_tile(species(s)%tile_nx,&
                                  & species(s)%tile_ny,&
                                  & species(s)%tile_nz))
    end do

    do s = 1, nspec
      write (var_name, "(A6,I1)") "maxptl", s
      call getInput('particles', var_name, maxptl_)
      write (var_name, "(A1,I1)") "m", s
      call getInput('particles', var_name, species(s)%m_sp)
      write (var_name, "(A2,I1)") "ch", s
      call getInput('particles', var_name, species(s)%ch_sp)

      write (var_name, "(A7,I1)") "deposit", s
      call getInput('particles', var_name, species(s)%deposit_sp, (species(s)%ch_sp .ne. 0))
      write (var_name, "(A4,I1)") "move", s
      call getInput('particles', var_name, species(s)%move_sp, .true.)
      write (var_name, "(A6,I1)") "output", s
      call getInput('particles', var_name, species(s)%output_sp, .true.)

      if ((species(s)%m_sp .eq. 0) .and. (species(s)%ch_sp .ne. 0)) then
        call throwError('ERROR: massless charged particles are not allowed')
      end if
      if ((species(s)%m_sp .ne. 0) .and. (species(s)%ch_sp .eq. 0)) then
        call throwError('ERROR: massive zero-charge particles are not allowed')
      end if
      if ((species(s)%ch_sp .eq. 0) .and. (species(s)%deposit_sp)) then
        call throwError('ERROR: zero-charged particles cannot deposit current')
      end if

      #ifdef DOWNSAMPLING
        write (var_name, "(A3,I1)") "dwn", s
        call getInput('particles', var_name, species(s)%dwn_sp, .false.)
      #endif

      ! extra physics properties

      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            call createEmptyTile(s, ti, tj, tk, maxptl_)
          end do
        end do
      end do
      species(s)%cntr_sp = 0
    end do
  end subroutine initializeParticles

  subroutine initializePrtlExchange()
    implicit none
    integer             :: buffsize, buffsize_x, buffsize_y
    integer             :: buffsize_xy
    integer             :: buffsize_z, buffsize_xz, buffsize_yz
    integer             :: buffsize_xyz
    integer             :: multiplier, ierr, ind1, ind2, ind3, ind
    integer             :: additional_real, additional_int, additional_int2

    #ifdef MPI08
      type(MPI_DATATYPE), dimension(0:2)            :: oldtypes
    #endif

    #ifdef MPI
      integer, dimension(0:2)                       :: oldtypes
    #endif

    integer, dimension(0:2)                         :: blockcounts
    integer(kind=MPI_ADDRESS_KIND), dimension(0:2)  :: offsets
    integer(kind=MPI_ADDRESS_KIND)                  :: extent_int2, extent_real, lb

    call getInput('grid', 'max_buff', max_buffsize, 100)

    multiplier = max(INT(ppc0), 1) * max_buffsize
    ! FIX this might change over time (due to load balancing)
    buffsize_x = 0
    buffsize_y = 0; buffsize_xy = 0
    buffsize_z = 0; buffsize_xz = 0; buffsize_yz = 0; buffsize_xyz = 0
    #if defined (oneD) || defined (twoD) || defined (threeD)
      buffsize_x = this_meshblock%ptr%sy * this_meshblock%ptr%sz * multiplier
    #endif

    #if defined (twoD) || defined (threeD)
      buffsize_y = this_meshblock%ptr%sx * this_meshblock%ptr%sz * multiplier
      buffsize_xy = this_meshblock%ptr%sz * multiplier
    #endif

    #if defined(threeD)
      buffsize_z = this_meshblock%ptr%sx * this_meshblock%ptr%sy * multiplier
      buffsize_xz = this_meshblock%ptr%sz * multiplier
      buffsize_yz = this_meshblock%ptr%sx * multiplier
      buffsize_xyz = multiplier
    #endif

    #ifdef oneD
      buffsize = multiplier
    #elif twoD
      buffsize = MAX0(this_meshblock%ptr%sx, this_meshblock%ptr%sy, this_meshblock%ptr%sz) * multiplier
    #elif threeD
      buffsize = MAX0(this_meshblock%ptr%sx, this_meshblock%ptr%sy, this_meshblock%ptr%sz)**2 * multiplier
    #endif

    allocate(recv_enroute(buffsize))

    do ind1 = -1, 1
      do ind2 = -1, 1
        do ind3 = -1, 1
          if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
          #ifdef oneD
            if ((ind2 .ne. 0) .or. (ind3 .ne. 0)) cycle
          #elif twoD
            if (ind3 .ne. 0) cycle
          #endif
          if ((ind2 .eq. 0) .and. (ind3 .eq. 0)) then
            buffsize = buffsize_x
          else if ((ind1 .eq. 0) .and. (ind3 .eq. 0)) then
            buffsize = buffsize_y
          else if ((ind1 .eq. 0) .and. (ind2 .eq. 0)) then
            buffsize = buffsize_z
          else if (ind3 .eq. 0) then
            buffsize = buffsize_xy
          else if (ind2 .eq. 0) then
            buffsize = buffsize_xz
          else if (ind1 .eq. 0) then
            buffsize = buffsize_yz
          else
            buffsize = buffsize_xyz
          end if
          enroute_bot%get(ind1, ind2, ind3)%max_send = buffsize
          allocate(enroute_bot%get(ind1, ind2, ind3)%send_enroute(buffsize))
        end do
      end do
    end do

    ! DEP_PRT [particle-dependent]
    ! new type for myMPI_ENROUTE
    additional_real = 0; additional_int = 0; additional_int2 = 0

    call MPI_TYPE_GET_EXTENT(MPI_INTEGER2, lb, extent_int2, ierr)
    call MPI_TYPE_GET_EXTENT(MPI_REAL, lb, extent_real, ierr)

    #ifdef PRTLPAYLOADS
      additional_real = additional_real + 3
    #endif

    !     # of blockcounts = 3:
    !       3  x integer2  [xi, yi, zi]
    !       7  x real      [dx, dy, dz, u, v, w, weight]
    !                                                         | + 3 if PRTLPAYLOADS
    !       2  x integer   [ind, proc]
    blockcounts(0) = 3 + additional_int2
    oldtypes(0) = MPI_INTEGER2
    blockcounts(1) = 7 + additional_real
    oldtypes(1) = MPI_REAL
    blockcounts(2) = 2 + additional_int
    oldtypes(2) = MPI_INTEGER

    offsets(0) = 0
    offsets(1) = blockcounts(0) * extent_int2 + offsets(0)
    offsets(2) = blockcounts(1) * extent_real + offsets(1)
    call MPI_TYPE_CREATE_STRUCT(3, blockcounts, offsets, oldtypes, myMPI_ENROUTE, ierr)
    call MPI_TYPE_COMMIT(myMPI_ENROUTE, ierr)
  end subroutine initializePrtlExchange

  subroutine initializeFields()
    implicit none
    if (allocated(ex)) deallocate(ex)
    if (allocated(ey)) deallocate(ey)
    if (allocated(ez)) deallocate(ez)
    if (allocated(bx)) deallocate(bx)
    if (allocated(by)) deallocate(by)
    if (allocated(bz)) deallocate(bz)
    if (allocated(jx)) deallocate(jx)
    if (allocated(jy)) deallocate(jy)
    if (allocated(jz)) deallocate(jz)
    if (allocated(jx_buff)) deallocate(jx_buff)
    if (allocated(jy_buff)) deallocate(jy_buff)
    if (allocated(jz_buff)) deallocate(jz_buff)
    #ifdef oneD
      allocate(ex(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(ey(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(ez(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(bx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(by(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(bz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(jx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(jy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(jz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(jx_buff(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(jy_buff(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(jz_buff(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(lg_arr(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      ! 20 = max # of fields sent in each direction
      sendrecv_offsetsz = NGHOST * 20
      ! 2 (~5) directions to send/recv in 1D
      sendrecv_buffsz = sendrecv_offsetsz * 5
    #elif twoD
      allocate(ex(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(ey(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(ez(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(bx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(by(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(bz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(jx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(jy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(jz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(jx_buff(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                     & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(jy_buff(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                     & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(jz_buff(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                     & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(lg_arr(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                    & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))

     ! 20 = max # of fields sent in each direction
     sendrecv_offsetsz = MAX0(this_meshblock%ptr%sx, this_meshblock%ptr%sy, this_meshblock%ptr%sz) * NGHOST * 20
     ! 8 (~10) directions to send/recv in 2D
     sendrecv_buffsz = sendrecv_offsetsz * 10
    #elif threeD
      allocate(ex(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(ey(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(ez(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(bx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(by(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(bz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(jx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(jy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(jz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(jx_buff(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                     & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                     & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(jy_buff(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                     & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                     & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(jz_buff(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                     & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                     & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(lg_arr(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                    & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                    & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      ! 20 = max # of fields sent in each direction
      sendrecv_offsetsz = MAX0(this_meshblock%ptr%sx, this_meshblock%ptr%sy, this_meshblock%ptr%sz)**2 * NGHOST * 20
      ! 26 (~30) directions to send/recv in 3D
      sendrecv_buffsz = sendrecv_offsetsz * 30
    #endif

    allocate(sm_arr(0:this_meshblock%ptr%sx - 1, 0:this_meshblock%ptr%sy - 1, 0:this_meshblock%ptr%sz - 1))

    if (allocated(send_fld)) deallocate(send_fld)
    if (allocated(recv_fld)) deallocate(recv_fld)
    allocate(send_fld(sendrecv_buffsz))
    allocate(recv_fld(sendrecv_offsetsz))
  end subroutine initializeFields

  subroutine firstRankInitialize()
    ! create output/restart directories
    !   if does not already exist
    !     note: some compilers may not support IFPORT
    #ifdef IFPORT
      logical :: result
      if (tot_output_enable .or. hst_enable) then
        result = makedirqq(trim(output_dir_name))
      end if
      if (rst_enable) then
        result = makedirqq(trim(restart_dir_name))
      end if
      if (slice_output_enable) then
        result = makedirqq(trim(slice_dir_name))
      end if
    #else
      if (tot_output_enable .or. hst_enable) then
        call system('mkdir -p ' // trim(output_dir_name))
      end if
      if (rst_enable) then
        call system('mkdir -p ' // trim(restart_dir_name))
      end if
      if (slice_output_enable) then
        call system('mkdir -p ' // trim(slice_dir_name))
      end if
    #endif
  end subroutine firstRankInitialize

  subroutine restartSimulation()
    implicit none
    character(len=STR_MAX)      :: mpichar, filename
    integer                     :: s, ti, tj, tk, num, pid
    integer                     :: dummy_int1, dummy_int2, dummy_int3
    real                        :: dummy_real
    write(mpichar, "(i8.8)") mpi_rank

    if (mpi_rank .eq. 0) then
      print *, 'Reading restart data...'
    end if

    ! loading fields
    filename = trim(restart_from) // '/flds.rst.' // trim(mpichar)
    open(UNIT_restart_fld, file=filename, form="unformatted")
    rewind(UNIT_restart_fld)
    read(UNIT_restart_fld) start_timestep, dseed, tot_output_index, slice_index
    read(UNIT_restart_fld) ex, ey, ez, bx, by, bz
    read(UNIT_restart_fld) CC, ppc0, c_omp, sigma
    close(UNIT_restart_fld)
    start_timestep = start_timestep + 1

    call renormalizeUnits()

    if (mpi_rank .eq. 0) then
      print *, '`CC`, `ppc0`, `c_omp` & `sigma` are read from restart ...'
      print *, '... values read from input are ignored.'
    end if

    ! loading particles
    filename = trim(restart_from) // '/prtl.rst.' // trim(mpichar)
    open(UNIT_restart_prtl, file=filename, form="unformatted")
    do s = 1, nspec
      read(UNIT_restart_prtl) species(s)%cntr_sp
      ! check that the tile sizes are the same
      read(UNIT_restart_prtl) dummy_int1, dummy_int2, dummy_int3
      if ((dummy_int1 .ne. species(s)%tile_sx) .or.&
        & (dummy_int2 .ne. species(s)%tile_sy) .or.&
        & (dummy_int3 .ne. species(s)%tile_sz)) then
        call throwError('ERROR. Wrong tile sizes after the restart')
      end if
      ! check that the # of tiles are the same
      read(UNIT_restart_prtl) dummy_int1, dummy_int2, dummy_int3
      if ((dummy_int1 .ne. species(s)%tile_nx) .or.&
        & (dummy_int2 .ne. species(s)%tile_ny) .or.&
        & (dummy_int3 .ne. species(s)%tile_nz)) then
        call throwError('ERROR. Wrong # of tiles after the restart')
      end if

      ! loop through all the tiles and read
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            read(UNIT_restart_prtl) species(s)%prtl_tile(ti, tj, tk)%spec

            ! reallocate the tile if necessary
            read(UNIT_restart_prtl) dummy_int1
            if (dummy_int1 .ne. species(s)%prtl_tile(ti, tj, tk)%maxptl_sp) then
              call allocateParticlesOnEmptyTile(s, species(s)%prtl_tile(ti, tj, tk), dummy_int1)
            end if
            read(UNIT_restart_prtl) species(s)%prtl_tile(ti, tj, tk)%npart_sp
            num = species(s)%prtl_tile(ti, tj, tk)%npart_sp

            ! read out and check tile dimensions
            read(UNIT_restart_prtl) dummy_int1
            if (dummy_int1 .ne. species(s)%prtl_tile(ti, tj, tk)%x1) then
              call throwError('ERROR. Wrong tile dimension after the restart')
            end if
            read(UNIT_restart_prtl) dummy_int1
            if (dummy_int1 .ne. species(s)%prtl_tile(ti, tj, tk)%x2) then
              call throwError('ERROR. Wrong tile dimension after the restart')
            end if
            read(UNIT_restart_prtl) dummy_int1
            if (dummy_int1 .ne. species(s)%prtl_tile(ti, tj, tk)%y1) then
              call throwError('ERROR. Wrong tile dimension after the restart')
            end if
            read(UNIT_restart_prtl) dummy_int1
            if (dummy_int1 .ne. species(s)%prtl_tile(ti, tj, tk)%y2) then
              call throwError('ERROR. Wrong tile dimension after the restart')
            end if
            read(UNIT_restart_prtl) dummy_int1
            if (dummy_int1 .ne. species(s)%prtl_tile(ti, tj, tk)%z1) then
              call throwError('ERROR. Wrong tile dimension after the restart')
            end if
            read(UNIT_restart_prtl) dummy_int1
            if (dummy_int1 .ne. species(s)%prtl_tile(ti, tj, tk)%z2) then
              call throwError('ERROR. Wrong tile dimension after the restart')
            end if

            ! finally read out all the particles
            ! DEP_PRT [particle-dependent]
            read(UNIT_restart_prtl) species(s)%prtl_tile(ti, tj, tk)%xi(1:num)
            read(UNIT_restart_prtl) species(s)%prtl_tile(ti, tj, tk)%yi(1:num)
            read(UNIT_restart_prtl) species(s)%prtl_tile(ti, tj, tk)%zi(1:num)
            read(UNIT_restart_prtl) species(s)%prtl_tile(ti, tj, tk)%dx(1:num)
            read(UNIT_restart_prtl) species(s)%prtl_tile(ti, tj, tk)%dy(1:num)
            read(UNIT_restart_prtl) species(s)%prtl_tile(ti, tj, tk)%dz(1:num)
            read(UNIT_restart_prtl) species(s)%prtl_tile(ti, tj, tk)%u(1:num)
            read(UNIT_restart_prtl) species(s)%prtl_tile(ti, tj, tk)%v(1:num)
            read(UNIT_restart_prtl) species(s)%prtl_tile(ti, tj, tk)%w(1:num)
            read(UNIT_restart_prtl) species(s)%prtl_tile(ti, tj, tk)%weight(1:num)
            read(UNIT_restart_prtl) species(s)%prtl_tile(ti, tj, tk)%ind(1:num)
            read(UNIT_restart_prtl) species(s)%prtl_tile(ti, tj, tk)%proc(1:num)

            #ifdef PRTLPAYLOADS
              read(UNIT_restart_prtl) species(s)%prtl_tile(ti, tj, tk)%payload1(1:num)
              read(UNIT_restart_prtl) species(s)%prtl_tile(ti, tj, tk)%payload2(1:num)
              read(UNIT_restart_prtl) species(s)%prtl_tile(ti, tj, tk)%payload3(1:num)
            #endif
          end do
        end do
      end do
    end do
    close(UNIT_restart_prtl)
  end subroutine restartSimulation

  subroutine checkEverything()
    implicit none
    ! check that the domain size is larger than the number of ghost zones
    #ifdef oneD
      if (this_meshblock%ptr%sx .lt. NGHOST) then
        call throwError('ERROR: ghost zones overflow the domain size in ' // trim(STR(mpi_rank)))
      end if
    #elif twoD
      if ((this_meshblock%ptr%sx .lt. NGHOST) .or.&
        & (this_meshblock%ptr%sy .lt. NGHOST)) then
        call throwError('ERROR: ghost zones overflow the domain size in ' // trim(STR(mpi_rank)))
      end if
    #elif threeD
      if ((this_meshblock%ptr%sx .lt. NGHOST) .or.&
        & (this_meshblock%ptr%sy .lt. NGHOST) .or.&
        & (this_meshblock%ptr%sz .lt. NGHOST)) then
        call throwError('ERROR: ghost zones overflow the domain size in ' // trim(STR(mpi_rank)))
      end if
    #endif
  end subroutine checkEverything

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
    end subroutine initializeDownsampling
  #endif
end module m_initialize
