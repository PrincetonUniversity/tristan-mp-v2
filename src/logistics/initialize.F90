#include "../defs.F90"

module m_initialize
  #ifdef IFPORT
    use ifport
  #endif
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_loadbalancing
  use m_particles
  use m_particlelogistics
  use m_fields
  use m_userfile
  use m_helpers
  use m_errors

  ! extra physics
  #ifdef RADIATION
    use m_radiation
  #endif

  #ifdef BWPAIRPRODUCTION
    use m_bwpairproduction
  #endif

  implicit none

  !--- PRIVATE functions -----------------------------------------!
  private :: initializeCommunications, initializeOutput,&
           & firstRankInitialize, initializeParticles,&
           & initializeLB,&
           & distributeMeshblocks, initializeDomain,&
           & initializePrtlExchange, initializeFields,&
           & initializeSimulation, checkEverything
  #ifdef RADIATION
    private :: initializeRadiation
  #endif

  #ifdef BWPAIRPRODUCTION
    private :: initializeBWPairProduction
  #endif
  !...............................................................!
contains
  ! initialize all the necessary arrays and variables
  subroutine initializeAll()
    implicit none
    call readCommandlineArgs()

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

    ! ADD possibility to define output function in userfile
    ! ADD hst file?

    call initializeFields()
      call printDiag((mpi_rank .eq. 0), "initializeFields()", .true.)

    call initializeParticles()
      call printDiag((mpi_rank .eq. 0), "initializeParticles()", .true.)

    #ifdef RADIATION
      call initializeRadiation()
        call printDiag((mpi_rank .eq. 0), "initializeRadiation()", .true.)
    #endif

    #ifdef BWPAIRPRODUCTION
      call initializeBWPairProduction()
        call printDiag((mpi_rank .eq. 0), "initializeBWPairProduction()", .true.)
    #endif

    call initializePrtlExchange()
      call printDiag((mpi_rank .eq. 0), "initializePrtlExchange()", .true.)

    call initializeRandomSeed(mpi_rank)
      call printDiag((mpi_rank .eq. 0), "initializeRandomSeed()", .true.)

    if (mpi_rank .eq. 0) then
      call firstRankInitialize()
      call printDiag(.true., "firstRankInitialize()", .true.)
    end if

    call userInitialize()
      call printDiag((mpi_rank .eq. 0), "userInitialize()", .true.)

    call checkEverything()
      call printDiag((mpi_rank .eq. 0), "checkEverything()", .true.)

    call printReport((mpi_rank .eq. 0), "InitializeAll()")
  end subroutine initializeAll

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
    call getInput('node_configuration', 'sizex', sizex)
    call getInput('node_configuration', 'sizey', sizey)
    #ifdef threeD
      call getInput('node_configuration', 'sizez', sizez)
    #else
      sizez = 1
    #endif
    global_mesh%x0 = 0
    global_mesh%y0 = 0
    global_mesh%z0 = 0
    call getInput('grid', 'mx0', global_mesh%sx)
    call getInput('grid', 'my0', global_mesh%sy)
    #ifdef threeD
      call getInput('grid', 'mz0', global_mesh%sz)
    #else
      global_mesh%sz = 1
    #endif

    if ((modulo(global_mesh%sx, sizex) .ne. 0) .and.&
      & (modulo(global_mesh%sy, sizey) .ne. 0) .and.&
      & (modulo(global_mesh%sz, sizez) .ne. 0)) then
      call throwError('ERROR: grid size is not evenly divisible by the number of cores')
    end if

    call getInput('grid', 'boundary_x', boundary_x, 1)
    call getInput('grid', 'boundary_y', boundary_y, 1)
    #ifdef threeD
      call getInput('grid', 'boundary_z', boundary_z, 1)
    #else
      boundary_z = 0
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
    call getInput('static_load_balancing', 'in_x', slb_x, .false.)
    call getInput('static_load_balancing', 'sx_min', slb_sxmin, 10)
    call getInput('static_load_balancing', 'in_y', slb_y, .false.)
    call getInput('static_load_balancing', 'sy_min', slb_symin, 10)
    #ifdef threeD
      call getInput('static_load_balancing', 'in_z', slb_z, .false.)
      call getInput('static_load_balancing', 'sz_min', slb_szmin, 10)
    #else
      slb_z = .false.
      slb_szmin = -1
    #endif

    ! initializing adaptive LB variables
    call getInput('adaptive_load_balancing', 'in_x', alb_x, .false.)
    call getInput('adaptive_load_balancing', 'in_y', alb_y, .false.)

    call getInput('adaptive_load_balancing', 'sx_min', alb_sxmin, 10)
    call getInput('adaptive_load_balancing', 'sy_min', alb_symin, 10)

    call getInput('adaptive_load_balancing', 'interval_x', alb_int_x, 1000)
    call getInput('adaptive_load_balancing', 'interval_y', alb_int_y, 1000)

    call getInput('adaptive_load_balancing', 'start_x', alb_start_x, 0)
    call getInput('adaptive_load_balancing', 'start_y', alb_start_y, 0)

    #ifdef threeD
      call getInput('adaptive_load_balancing', 'in_z', alb_z, .false.)
      call getInput('adaptive_load_balancing', 'sz_min', alb_szmin, 10)
      call getInput('adaptive_load_balancing', 'interval_z', alb_int_z, 1000)
      call getInput('adaptive_load_balancing', 'start_z', alb_start_z, 0)
    #else
      alb_z = .false.
      alb_szmin = -1
      alb_int_z = -1
      alb_start_z = -1
    #endif
  end subroutine initializeLB

  subroutine initializeOutput()
    implicit none
    call getInput('output', 'start', output_start, 0)
    call getInput('output', 'interval', output_interval, 10)
    call getInput('output', 'stride', output_stride, 10)
    call getInput('output', 'istep', output_istep, 4)

    call getInput('output', 'spec_min', spec_min, 1e-2)
    call getInput('output', 'spec_max', spec_max, 1e2)
    call getInput('output', 'spec_num', spec_num, 100)
    spec_min = log(spec_min)
    spec_max = log(spec_max)

    call getInput('output', 'flds_at_prtl', flds_at_prtl, .false.)

    #ifdef HDF5

    #ifdef MPI08
      h5comm = MPI_COMM_WORLD%MPI_VAL
      h5info = MPI_INFO_NULL%MPI_VAL
    #endif

    #ifdef MPI
      h5comm = MPI_COMM_WORLD
      h5info = MPI_INFO_NULL
    #endif

    #endif
  end subroutine initializeOutput

  subroutine initializeSimulation()
    implicit none
    call getInput('time', 'last', final_timestep, 1000)
    call getInput('algorithm', 'nfilter', nfilter, 16)
    call getInput('grid', 'resize_tiles', resize_tiles, .false.)
    call getInput('grid', 'min_tile_nprt', min_tile_nprt, 100)
  end subroutine initializeSimulation

  subroutine initializeParticles()
    implicit none
    integer                 :: s, ti, tj, tk
    character(len=STR_MAX)  :: var_name
    integer                  :: maxptl_

    call getInput('plasma', 'ppc0', ppc0)
    call getInput('plasma', 'sigma', sigma)
    call getInput('plasma', 'c_omp', c_omp)
    B_norm = CC**2 * sqrt(sigma) / c_omp
    unit_ch = CC**2 / (ppc0 * c_omp**2)

    call getInput('particles', 'nspec', nspec, 2)

    allocate(species(nspec))
    do s = 1, nspec
      call getInput('grid', 'tileX', species(s)%tile_sx)
      call getInput('grid', 'tileY', species(s)%tile_sy)
      call getInput('grid', 'tileZ', species(s)%tile_sz)
      #ifndef threeD
        species(s)%tile_sz = 1
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

      ! extra physics properties
      #ifdef RADIATION
        write (var_name, "(A4,I1)") "cool", s
        call getInput('particles', var_name, species(s)%cool_sp, .false.)
        if ((species(s)%cool_sp) .and. (species(s)%m_sp .eq. 0)) then
          call throwError('Unable to cool `m=0` particles.')
        end if
      #endif

      #ifdef BWPAIRPRODUCTION
        write (var_name, "(A2,I1)") "bw", s
        call getInput('particles', var_name, species(s)%bw_sp, 0)
        if ((species(s)%bw_sp .ne. 0) .and.&
          & ((species(s)%ch_sp .ne. 0) .or. (species(s)%m_sp .ne. 0))) then
          call throwError('`ch != 0` or `m != 0` particles cannot pair-produce via BW.')
        end if
        if ((species(s)%bw_sp .gt. 2)) then
          call throwError('only two BW photon populations are allowed.')
        end if
      #endif

      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            species(s)%prtl_tile(ti, tj, tk)%maxptl_sp = maxptl_ / &
                              & (species(s)%tile_nx * species(s)%tile_ny * species(s)%tile_nz)
            species(s)%prtl_tile(ti, tj, tk)%npart_sp = 0

            species(s)%prtl_tile(ti, tj, tk)%x1 = (ti - 1) * species(s)%tile_sx
            species(s)%prtl_tile(ti, tj, tk)%x2 = min(ti * species(s)%tile_sx, this_meshblock%ptr%sx)
            species(s)%prtl_tile(ti, tj, tk)%y1 = (tj - 1) * species(s)%tile_sy
            species(s)%prtl_tile(ti, tj, tk)%y2 = min(tj * species(s)%tile_sy, this_meshblock%ptr%sy)
            species(s)%prtl_tile(ti, tj, tk)%z1 = (tk - 1) * species(s)%tile_sz
            species(s)%prtl_tile(ti, tj, tk)%z2 = min(tk * species(s)%tile_sz, this_meshblock%ptr%sz)
            #ifdef DEBUG
              if ((species(s)%prtl_tile(ti, tj, tk)%x1 .eq. 0) .and.&
                & (species(s)%prtl_tile(ti, tj, tk)%x2 .eq. 0) .and.&
                & (species(s)%prtl_tile(ti, tj, tk)%y1 .eq. 0) .and.&
                & (species(s)%prtl_tile(ti, tj, tk)%y2 .eq. 0) .and.&
                & (species(s)%prtl_tile(ti, tj, tk)%z1 .eq. 0) .and.&
                & (species(s)%prtl_tile(ti, tj, tk)%z2 .eq. 0)) then
                print *, ti, tj, tk
                print *, species(s)%prtl_tile(ti, tj, tk)%x1,&
                 & species(s)%prtl_tile(ti, tj, tk)%x2,&
                 & species(s)%prtl_tile(ti, tj, tk)%y1,&
                 & species(s)%prtl_tile(ti, tj, tk)%y2,&
                 & species(s)%prtl_tile(ti, tj, tk)%z1,&
                 & species(s)%prtl_tile(ti, tj, tk)%z2
               call throwError('ERROR IN PRTLINIT')
              end if
            #endif

            call allocateParticles(species(s)%prtl_tile(ti, tj, tk),&
                                 & species(s)%prtl_tile(ti, tj, tk)%maxptl_sp)
          end do
        end do
      end do
      species(s)%cntr_sp = 0
    end do
  end subroutine initializeParticles

  subroutine initializePrtlExchange()
    implicit none
    integer            :: buffsize, buffsize_x, buffsize_y
    integer            :: buffsize_xy
    integer            :: buffsize_z, buffsize_xz, buffsize_yz
    integer            :: buffsize_xyz
    integer            :: multiplier, ierr, ind1, ind2, ind3

    #ifdef MPI08
      type(MPI_DATATYPE), dimension(0:2)            :: oldtypes
    #endif

    #ifdef MPI
      integer, dimension(0:2)                       :: oldtypes
    #endif

    integer, dimension(0:2)                         :: blockcounts
    integer(kind=MPI_ADDRESS_KIND), dimension(0:2)  :: offsets
    integer(kind=MPI_ADDRESS_KIND)                  :: extent_int2, extent_real, lb

    multiplier = max(INT(ppc0), 1) * 1000
    ! FIX this might change over time (due to load balancing)
    buffsize_x = this_meshblock%ptr%sy * this_meshblock%ptr%sz * multiplier
    buffsize_y = this_meshblock%ptr%sx * this_meshblock%ptr%sz * multiplier
    buffsize_xy = this_meshblock%ptr%sz * multiplier
    #ifdef threeD
      buffsize = MAX0(this_meshblock%ptr%sx, this_meshblock%ptr%sy, this_meshblock%ptr%sz)**2 * multiplier

      buffsize_z = this_meshblock%ptr%sx * this_meshblock%ptr%sy * multiplier
      buffsize_xz = this_meshblock%ptr%sz * multiplier
      buffsize_yz = this_meshblock%ptr%sx * multiplier
      buffsize_xyz = multiplier
    #else
      buffsize = MAX0(this_meshblock%ptr%sx, this_meshblock%ptr%sy, this_meshblock%ptr%sz) * multiplier

      buffsize_z = 0; buffsize_xz = 0; buffsize_yz = 0; buffsize_xyz = 0
    #endif

    allocate(recv_enroute(buffsize))

    do ind1 = -1, 1
      do ind2 = -1, 1
        do ind3 = -1, 1
          if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
          #ifndef threeD
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
          allocate(enroute_bot%get(ind1, ind2, ind3)%send_enroute(buffsize))
        end do
      end do
    end do

    ! DEP_PRT [particle-dependent]
    ! new type for myMPI_ENROUTE
    !   BY DEFAULT:
    !     # of blockcounts = 3:
    !       3 x integer2  [xi, yi, zi]
    !       6 x real      [dx, dy, dz, u, v, w]
    !       2 x integer   [ind, proc]
    call MPI_TYPE_GET_EXTENT(MPI_INTEGER2, lb, extent_int2, ierr)
    call MPI_TYPE_GET_EXTENT(MPI_REAL, lb, extent_real, ierr)
    blockcounts(0) = 3
    oldtypes(0) = MPI_INTEGER2
    blockcounts(1) = 6
    oldtypes(1) = MPI_REAL
    blockcounts(2) = 2
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
    allocate(ex(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
              & fldBoundZ))
    allocate(ey(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
              & fldBoundZ))
    allocate(ez(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
              & fldBoundZ))
    allocate(bx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
              & fldBoundZ))
    allocate(by(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
              & fldBoundZ))
    allocate(bz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
              & fldBoundZ))
    allocate(jx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
              & fldBoundZ))
    allocate(jy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
              & fldBoundZ))
    allocate(jz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
              & fldBoundZ))

    ! exchange fields
    ! 20 = max # of fields sent in each direction
    #ifndef threeD
      sendrecv_offsetsz = MAX0(this_meshblock%ptr%sx, this_meshblock%ptr%sy, this_meshblock%ptr%sz) * NGHOST * 20
      sendrecv_buffsz = sendrecv_offsetsz * 10
      ! 8 (~10) directions to send/recv in 2D
    #else
      sendrecv_offsetsz = MAX0(this_meshblock%ptr%sx, this_meshblock%ptr%sy, this_meshblock%ptr%sz)**2 * NGHOST * 20
      sendrecv_buffsz = sendrecv_offsetsz * 30
      ! 26 (~30) directions to send/recv in 3D
    #endif

    if (allocated(send_fld)) deallocate(send_fld)
    allocate(send_fld(sendrecv_buffsz))
    if (allocated(recv_fld)) deallocate(recv_fld)
    allocate(recv_fld(sendrecv_offsetsz))

    ! output fields
    if (allocated(lg_arr)) deallocate(lg_arr)
    if (allocated(sm_arr)) deallocate(sm_arr)
    allocate(lg_arr(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                  & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                  & fldBoundZ))
    allocate(sm_arr(0:this_meshblock%ptr%sx - 1,&
                  & 0:this_meshblock%ptr%sy - 1,&
                  & 0:this_meshblock%ptr%sz - 1))
  end subroutine initializeFields

  subroutine firstRankInitialize()
    ! create output/restart directories
    !   if does not already exist
    !     note: some compilers may not support IFPORT
    #ifdef IFPORT
      logical :: result
      result = makedirqq(trim(output_dir_name))
      result = makedirqq(trim(restart_dir_name))
    #else
      call system('mkdir -p ' // trim(output_dir_name))
      call system('mkdir -p ' // trim(restart_dir_name))
    #endif
  end subroutine firstRankInitialize

  subroutine checkEverything()
    implicit none
    ! check that the domain size is larger than the number of ghost zones
    #ifndef threeD
      if ((this_meshblock%ptr%sx .lt. NGHOST) .or.&
        & (this_meshblock%ptr%sy .lt. NGHOST)) then
        call throwError('ERROR: ghost zones overflow the domain size in ' // trim(STR(mpi_rank)))
      end if
    #else
      if ((this_meshblock%ptr%sx .lt. NGHOST) .or.&
        & (this_meshblock%ptr%sy .lt. NGHOST) .or.&
        & (this_meshblock%ptr%sz .lt. NGHOST)) then
        call throwError('ERROR: ghost zones overflow the domain size in ' // trim(STR(mpi_rank)))
      end if
    #endif
  end subroutine checkEverything

  ! extra physics
  #ifdef RADIATION
    subroutine initializeRadiation()
      implicit none
      call getInput('radiation', 'emit_gamma_syn', emit_gamma_syn, 10.0)
      call getInput('radiation', 'emit_gamma_ic', emit_gamma_ic, 10.0)
      call getInput('radiation', 'gamma_syn', cool_gamma_syn, 10.0)
      call getInput('radiation', 'gamma_ic', cool_gamma_ic, 10.0)
      call getInput('radiation', 'beta_rec', rad_beta_rec, 0.1)
      call getInput('radiation', 'dens_limit', rad_dens_lim, 1e8)
      #ifdef EMIT
        call getInput('radiation', 'photon_sp', rad_photon_sp, 3)
        if ((nspec .lt. rad_photon_sp) .or.&
          & (species(rad_photon_sp)%ch_sp .ne. 0) .or.&
          & (species(rad_photon_sp)%m_sp .ne. 0)) then
          call throwError('Wrong choice of `photon_sp`.')
        end if
      #endif

      if (.not. allocated(rad_spectra)) allocate(rad_spectra(nspec, spec_num))
      if (.not. allocated(glob_rad_spectra)) allocate(glob_rad_spectra(nspec, spec_num))
      rad_spectra(:, :) = 0.0
      glob_rad_spectra(:, :) = 0.0
    end subroutine initializeRadiation
  #endif

  #ifdef BWPAIRPRODUCTION
    subroutine initializeBWPairProduction()
      implicit none
      call getInput('bw_pp', 'tau_BW', BW_tau)
      call getInput('bw_pp', 'interval', BW_interval, 1)
      call getInput('bw_pp', 'algorithm', BW_algorithm)
      call getInput('bw_pp', 'electron_sp', BW_electron_sp, 1)
      call getInput('bw_pp', 'positron_sp', BW_positron_sp, 2)
    end subroutine initializeBWPairProduction
  #endif

end module m_initialize
