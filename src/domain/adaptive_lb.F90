#include "../defs.F90"

module m_adaptivelb
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_helpers, only: reassignNeighborsForAll, indToRnk

  use m_fieldlogistics, only: restoreFieldsFromBackups, deallocateFieldBackups,&
                            & deallocateFields, reallocateFields, reallocateFieldBuffers,&
                            & backupEBfields

  use m_particlelogistics, only: shiftParticlesX, shiftParticlesY, shiftParticlesZ, backupParticles,&
                               & reallocateParticles, restoreParticlesFromBackup,&
                               & deallocateParticleBackup, clearGhostParticles
  use m_exchangeparts, only: redistributeParticlesBetweenMeshblocks
  use m_exchangefields, only: exchangeFieldSlabInX, exchangeFieldSlabInY, exchangeFieldSlabInZ

  implicit none
contains
  subroutine redistributeMeshblocksALB(step)
    implicit none
    integer, optional, intent(in) :: step
    integer :: valueRSS

    if (alb_x .and. (step .ge. alb_start_x) .and. (mod(step, alb_int_x) .eq. 0)) call metaRedistInX()
    if (alb_y .and. (step .ge. alb_start_y) .and. (mod(step, alb_int_y) .eq. 0)) call metaRedistInY()
    if (alb_z .and. (step .ge. alb_start_z) .and. (mod(step, alb_int_z) .eq. 0)) call metaRedistInZ()

    call printDiag("redistributeMeshblocksALB()", 1)
  end subroutine redistributeMeshblocksALB

  ! the following subroutine computes the load
  !   (e.g. # of particles) for any given domain
  ! this can be arbitrary function for which the algorithm
  !   balances the domain distribution
  subroutine computeLoadALB(load_ALB)
    implicit none
    integer, intent(out) :: load_ALB
    integer              :: s, ti, tj, tk
    load_ALB = 0
    do s = 1, nspec
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            load_ALB = load_ALB + species(s)%prtl_tile(ti, tj, tk)%npart_sp
          end do
        end do
      end do
    end do
  end subroutine computeLoadALB

  ! accumulate loads from all the sources
  subroutine accumulateLoads()
    implicit none
    integer :: ierr
    if (.not. allocated(lb_load_glob)) allocate(lb_load_glob(mpi_size))
    call MPI_ALLGATHER(lb_load, 1, MPI_INTEGER,&
                     & lb_load_glob, 1, MPI_INTEGER,&
                     & MPI_COMM_WORLD, ierr)
  end subroutine accumulateLoads

  ! Return true if loadA is too large and needs balancing
  function isImbalanced(loadA, loadB)
    logical :: isImbalanced
    ! Two loads for comparison
    integer, intent(in)  :: loadA, loadB
    isImbalanced = (loadA .gt. INT(REAL(loadB) * 1.05))
    return
  end function isImbalanced

  subroutine metaRedistInX()
    implicit none
    integer               :: delta_i, i, j, k, cnt, ierr
    integer               :: load_x0, load_x1, sx0_old, sx1_old, inds(3)

    if (.not. allocated(lb_group_x0)) allocate(lb_group_x0(sizey * sizez))
    if (.not. allocated(lb_group_x1)) allocate(lb_group_x1(sizey * sizez))

    do delta_i = 0, sizex - 2
      ! select the left and right domain slabs
      cnt = 1
      do k = 0, sizez - 1
        do j = 0, sizey - 1
          inds(1) = delta_i; inds(2) = j; inds(3) = k
          lb_group_x0(cnt) = indToRnk(inds)
          inds(1) = delta_i + 1;
          lb_group_x1(cnt) = indToRnk(inds)
          cnt = cnt + 1
        end do
      end do
      ! now all the actions are between these two slabs

      ! old dimensions
      sx0_old = meshblocks(lb_group_x0(1) + 1)%sx
      sx1_old = meshblocks(lb_group_x1(1) + 1)%sx

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call computeLoadALB(lb_load)
      call accumulateLoads()

      ! cumulative loads for each slab
      load_x0 = 0; load_x1 = 0
      do cnt = 1, sizey * sizez
        load_x0 = load_x0 + lb_load_glob(lb_group_x0(cnt) + 1)
        load_x1 = load_x1 + lb_load_glob(lb_group_x1(cnt) + 1)
      end do

      ! Apply balancing driver and checks befor requesting re-allocation
      if (isImbalanced(load_x0, load_x1) .and. &
       & (sx0_old - NGHOST - 1 .gt. alb_slab_x).and. &
       & (sx0_old - alb_slab_x .gt. alb_sxmin)) then
        call reshapeInX(lb_group_x0, lb_group_x1, -alb_slab_x)
      else if (isImbalanced(load_x1, load_x0) .and. &
            & (sx1_old - NGHOST - 1 .gt. alb_slab_x) .and. &
            & (sx1_old - alb_slab_x .gt. alb_sxmin)) then
        call reshapeInX(lb_group_x0, lb_group_x1, alb_slab_x)
      end if
    end do

    if (allocated(lb_group_x0)) deallocate(lb_group_x0)
    if (allocated(lb_group_x1)) deallocate(lb_group_x1)

    call printDiag("metaRedistInX()", 2)
  end subroutine metaRedistInX

  subroutine metaRedistInY()
    implicit none
    integer               :: delta_j, i, j, k, cnt, ierr
    integer               :: load_y0, load_y1, sy0_old, sy1_old, inds(3)

    if (.not. allocated(lb_group_y0)) allocate(lb_group_y0(sizex * sizez))
    if (.not. allocated(lb_group_y1)) allocate(lb_group_y1(sizex * sizez))

    do delta_j = 0, sizey - 2
      ! select the left and right domain slabs
      cnt = 1
      do k = 0, sizez - 1
        do i = 0, sizex - 1
          inds(1) = i; inds(2) = delta_j; inds(3) = k
          lb_group_y0(cnt) = indToRnk(inds)
          inds(2) = delta_j + 1;
          lb_group_y1(cnt) = indToRnk(inds)
          cnt = cnt + 1
        end do
      end do
      ! now all the actions are between these two slabs

      ! old dimensions
      sy0_old = meshblocks(lb_group_y0(1) + 1)%sy
      sy1_old = meshblocks(lb_group_y1(1) + 1)%sy

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call computeLoadALB(lb_load)
      call accumulateLoads()

      ! cumulative loads for each slab
      load_y0 = 0; load_y1 = 0
      do cnt = 1, sizex * sizez
        load_y0 = load_y0 + lb_load_glob(lb_group_y0(cnt) + 1)
        load_y1 = load_y1 + lb_load_glob(lb_group_y1(cnt) + 1)
      end do

      ! Apply balancing driver and checks befor requesting re-allocation
      if (isImbalanced(load_y0, load_y1) .and. &
       & (sy0_old - NGHOST - 1 .gt. alb_slab_y).and. &
       & (sy0_old - alb_slab_y .gt. alb_symin)) then
        call reshapeInY(lb_group_y0, lb_group_y1, -alb_slab_y)
      else if (isImbalanced(load_y1, load_y0) .and. &
            & (sy1_old - NGHOST - 1 .gt. alb_slab_y) .and. &
            & (sy1_old - alb_slab_y .gt. alb_symin)) then
        call reshapeInY(lb_group_y0, lb_group_y1, alb_slab_y)
      end if
    end do

    if (allocated(lb_group_y0)) deallocate(lb_group_y0)
    if (allocated(lb_group_y1)) deallocate(lb_group_y1)

    call printDiag("metaRedistInY()", 2)
  end subroutine metaRedistInY

  subroutine metaRedistInZ()
    implicit none
    integer               :: delta_k, i, j, k, cnt, ierr
    integer               :: load_z0, load_z1, sz0_old, sz1_old, inds(3)

    if (.not. allocated(lb_group_z0)) allocate(lb_group_z0(sizex * sizey))
    if (.not. allocated(lb_group_z1)) allocate(lb_group_z1(sizex * sizey))

    do delta_k = 0, sizez - 2
      ! select the left and right domain slabs
      cnt = 1
      do j = 0, sizey - 1
        do i = 0, sizex - 1
          inds(1) = i; inds(2) = j; inds(3) = delta_k
          lb_group_z0(cnt) = indToRnk(inds)
          inds(3) = delta_k + 1;
          lb_group_z1(cnt) = indToRnk(inds)
          cnt = cnt + 1
        end do
      end do
      ! now all the actions are between these two slabs

      ! old dimensions
      sz0_old = meshblocks(lb_group_z0(1) + 1)%sz
      sz1_old = meshblocks(lb_group_z1(1) + 1)%sz

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call computeLoadALB(lb_load)
      call accumulateLoads()

      ! cumulative loads for each slab
      load_z0 = 0; load_z1 = 0
      do cnt = 1, sizex * sizez
        load_z0 = load_z0 + lb_load_glob(lb_group_z0(cnt) + 1)
        load_z1 = load_z1 + lb_load_glob(lb_group_z1(cnt) + 1)
      end do

      ! Apply balancing driver and checks befor requesting re-allocation
      if (isImbalanced(load_z0, load_z1) .and. &
       & (sz0_old - NGHOST - 1 .gt. alb_slab_z).and. &
       & (sz0_old - alb_slab_z .gt. alb_szmin)) then
        call reshapeInZ(lb_group_z0, lb_group_z1, -alb_slab_z)
      else if (isImbalanced(load_z1, load_z0) .and. &
            & (sz1_old - NGHOST - 1 .gt. alb_slab_z) .and. &
            & (sz1_old - alb_slab_z .gt. alb_szmin)) then
        call reshapeInZ(lb_group_z0, lb_group_z1, alb_slab_z)
      end if
    end do

    if (allocated(lb_group_z0)) deallocate(lb_group_z0)
    if (allocated(lb_group_z1)) deallocate(lb_group_z1)

    call printDiag("metaRedistInZ()", 2)
  end subroutine metaRedistInZ

  subroutine reshapeGlobalInX(proc_group, INJECT)
    integer, allocatable, intent(in)    :: proc_group(:)
    integer, intent(in)                 :: INJECT
    integer                             :: i, j, k, q, nproc_group, rnk, inds(3), ierr
    integer                             :: i1_from, i2_from, j1_from, j2_from, k1_from, k2_from,&
                                            & i1_to, i2_to, j1_to, j2_to, k1_to, k2_to
    logical                             :: leftmost, rightmost

    nproc_group = size(proc_group)

    leftmost = .true.
    rightmost = .true.
    do q = 1, nproc_group
      rnk = proc_group(q)
      if ((global_mesh%sx).ne.(meshblocks(rnk + 1)%x0 + meshblocks(rnk + 1)%sx)) rightmost = .false.
      if ((0).ne.(meshblocks(rnk + 1)%x0)) leftmost = .false.
    end do

    if ((.not. leftmost) .and. (.not. rightmost)) then
      call throwError('ERROR: Wrong ranks in reshapeGlobalInX')
    else if ((leftmost) .and. (rightmost)) then
      leftmost = .false.
    end if

    ! backup the fields with current sizes
    do q = 1, nproc_group
      rnk = proc_group(q)
      if (mpi_rank .eq. rnk) then
        call backupEBfields()
      end if
    end do

    ! get new meshblock dimensions
    new_meshblocks(:) = meshblocks(:)
    call reassignNeighborsForAll(new_meshblocks)
    do q = 1, nproc_group
      rnk = proc_group(q)
      new_meshblocks(rnk + 1)%sx = new_meshblocks(rnk + 1)%sx + INJECT
    end do

    if (leftmost) then
      do k = 0, sizez - 1
        do j = 0, sizey - 1
          do i = 1, sizex - 1
            inds(1) = i; inds(2) = j; inds(3) = k
            rnk = indToRnk(inds)
            new_meshblocks(rnk + 1)%x0 = new_meshblocks(rnk + 1)%x0 + INJECT
          end do
        end do
      end do
    end if

    ! reallocate field arrays given the new meshblock dimensions
    do q = 1, nproc_group
      rnk = proc_group(q)
      if (mpi_rank .eq. rnk) then
        call deallocateFields()
        call reallocateFields(new_meshblocks(mpi_rank + 1))
        call reallocateFieldBuffers(new_meshblocks(mpi_rank + 1))
      end if
    end do

    ! ... and recover fields from backup
    do q = 1, nproc_group
      rnk = proc_group(q)
      if ((mpi_rank .eq. rnk)) then
        ! recover from backup
        if (INJECT .gt. 0) then
          if (rightmost) then
            i1_to = 0; i2_to = this_meshblock%ptr%sx - 1 - abs(INJECT)
            j1_to = 0; j2_to = this_meshblock%ptr%sy - 1
            k1_to = 0; k2_to = this_meshblock%ptr%sz - 1
            i1_from = i1_to; i2_from = i2_to
            j1_from = j1_to; j2_from = j2_to
            k1_from = k1_to; k2_from = k2_to
          else if (leftmost) then
            i1_to = abs(INJECT); i2_to = new_meshblocks(rnk + 1)%sx - 1
            j1_to = 0; j2_to = new_meshblocks(rnk + 1)%sy - 1
            k1_to = 0; k2_to = new_meshblocks(rnk + 1)%sz - 1
            i1_from = 0; i2_from = this_meshblock%ptr%sx - 1
            j1_from = 0; j2_from = this_meshblock%ptr%sy - 1
            k1_from = 0; k2_from = this_meshblock%ptr%sz - 1
          end if
        else
          ! `left_rnk` shrinks
          ! `right_rnk` inflates
          if (rightmost) then
            i1_to = 0; i2_to = this_meshblock%ptr%sx - 1
            j1_to = 0; j2_to = this_meshblock%ptr%sy - 1
            k1_to = 0; k2_to = this_meshblock%ptr%sz - 1
            i1_from = i1_to; i2_from = i2_to
            j1_from = j1_to; j2_from = j2_to
            k1_from = k1_to; k2_from = k2_to
          else if (leftmost) then
            i1_to = 0; i2_to = new_meshblocks(rnk + 1)%sx - 1
            j1_to = 0; j2_to = new_meshblocks(rnk + 1)%sy - 1
            k1_to = 0; k2_to = new_meshblocks(rnk + 1)%sz - 1
            i1_from = abs(INJECT); i2_from = this_meshblock%ptr%sx - 1
            j1_from = 0; j2_from = this_meshblock%ptr%sy - 1
            k1_from = 0; k2_from = this_meshblock%ptr%sz - 1
          end if
        end if

        call restoreFieldsFromBackups(i1_from, i2_from, j1_from, j2_from, k1_from, k2_from,&
                                    & i1_to, i2_to, j1_to, j2_to, k1_to, k2_to)
      end if
    end do

    ! resize the meshblocks
    meshblocks(:) = new_meshblocks(:)

    ! deallocate buffers and redistribute particles
    do q = 1, nproc_group
      rnk = proc_group(q)
      if ((mpi_rank .eq. rnk)) then
        call deallocateFieldBackups()

        if (leftmost) then
          ! INJECT particles
          call shiftParticlesX(INJECT)
        end if

        ! backup particles
        call backupParticles()

        ! reshuffle particle tiles
        call reallocateParticles(this_meshblock%ptr)
        ! restore particles from backup
        call restoreParticlesFromBackup()
        call deallocateParticleBackup()

      end if
    end do

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    ! put particles back on proper meshblocks
    call redistributeParticlesBetweenMeshblocks()
    call clearGhostParticles()

    call printDiag("reshapeGlobalInX()", 3)

    global_mesh%sx = global_mesh%sx + INJECT

  end subroutine reshapeGlobalInX

  subroutine reshapeInX(left_group, right_group, SHIFT)
    integer, allocatable, intent(in)    :: left_group(:), right_group(:)
    integer, intent(in)                 :: SHIFT
    integer                       :: nproc_group, q, left_rnk, right_rnk, ierr
    integer                       :: new_sx, new_sy, new_sz
    integer                       :: i1_from, i2_from, j1_from, j2_from, k1_from, k2_from,&
                                   & i1_to, i2_to, j1_to, j2_to, k1_to, k2_to
    nproc_group = size(left_group)
    #ifdef DEBUG
      if (size(left_group) .ne. size(right_group)) then
        call throwError('ERROR: wrong groups specified for `reshapeInX`.')
      end if
    #endif

    ! backup the fields with current sizes
    do q = 1, nproc_group
      left_rnk = left_group(q)
      right_rnk = right_group(q)
      if ((mpi_rank .eq. left_rnk) .or. (mpi_rank .eq. right_rnk)) then
        call backupEBfields()
      end if
    end do

    ! get new meshblock dimensions
    new_meshblocks(:) = meshblocks(:)
    call reassignNeighborsForAll(new_meshblocks)
    do q = 1, nproc_group
      left_rnk = left_group(q)
      right_rnk = right_group(q)
      new_meshblocks(left_rnk + 1)%sx = new_meshblocks(left_rnk + 1)%sx + SHIFT
      new_meshblocks(right_rnk + 1)%sx = new_meshblocks(right_rnk + 1)%sx - SHIFT
      new_meshblocks(right_rnk + 1)%x0 = new_meshblocks(right_rnk + 1)%x0 + SHIFT
    end do
    ! at this point DO NOT CHANGE `meshblocks` ...
    ! ... as the `exchangeFieldSlabIn*` still assumes old dimensions

    ! reallocate field arrays given the new meshblock dimensions
    do q = 1, nproc_group
      left_rnk = left_group(q)
      right_rnk = right_group(q)
      if ((mpi_rank .eq. left_rnk) .or. (mpi_rank .eq. right_rnk)) then
        call deallocateFields()
        call reallocateFields(new_meshblocks(mpi_rank + 1))
        call reallocateFieldBuffers(new_meshblocks(mpi_rank + 1))
      end if
    end do

    ! send/recv missing fields
    ! ... and recover from backup
    do q = 1, nproc_group
      left_rnk = left_group(q)
      right_rnk = right_group(q)
      if ((mpi_rank .eq. left_rnk) .or. (mpi_rank .eq. right_rnk)) then
        ! exchange slab
        call exchangeFieldSlabInX(left_rnk, right_rnk, SHIFT)

        ! recover from backup
        if (SHIFT .gt. 0) then
          ! `left_rnk` inflates
          ! `right_rnk` shrinks
          if (mpi_rank .eq. left_rnk) then
            i1_to = 0; i2_to = this_meshblock%ptr%sx - 1
            j1_to = 0; j2_to = this_meshblock%ptr%sy - 1
            k1_to = 0; k2_to = this_meshblock%ptr%sz - 1
            i1_from = i1_to; i2_from = i2_to
            j1_from = j1_to; j2_from = j2_to
            k1_from = k1_to; k2_from = k2_to
          else if (mpi_rank .eq. right_rnk) then
            i1_to = 0; i2_to = new_meshblocks(right_rnk + 1)%sx - 1
            j1_to = 0; j2_to = new_meshblocks(right_rnk + 1)%sy - 1
            k1_to = 0; k2_to = new_meshblocks(right_rnk + 1)%sz - 1
            i1_from = SHIFT; i2_from = this_meshblock%ptr%sx - 1
            j1_from = 0; j2_from = this_meshblock%ptr%sy - 1
            k1_from = 0; k2_from = this_meshblock%ptr%sz - 1
          end if
        else
          ! `left_rnk` shrinks
          ! `right_rnk` inflates
          if (mpi_rank .eq. left_rnk) then
            i1_to = 0; i2_to = this_meshblock%ptr%sx - 1 - abs(SHIFT)
            j1_to = 0; j2_to = this_meshblock%ptr%sy - 1
            k1_to = 0; k2_to = this_meshblock%ptr%sz - 1
            i1_from = i1_to; i2_from = i2_to
            j1_from = j1_to; j2_from = j2_to
            k1_from = k1_to; k2_from = k2_to
          else if (mpi_rank .eq. right_rnk) then
            i1_to = abs(SHIFT); i2_to = new_meshblocks(right_rnk + 1)%sx - 1
            j1_to = 0; j2_to = new_meshblocks(right_rnk + 1)%sy - 1
            k1_to = 0; k2_to = new_meshblocks(right_rnk + 1)%sz - 1
            i1_from = 0; i2_from = this_meshblock%ptr%sx - 1
            j1_from = 0; j2_from = this_meshblock%ptr%sy - 1
            k1_from = 0; k2_from = this_meshblock%ptr%sz - 1
          end if
        end if

        call restoreFieldsFromBackups(i1_from, i2_from, j1_from, j2_from, k1_from, k2_from,&
                                    & i1_to, i2_to, j1_to, j2_to, k1_to, k2_to)
      end if
    end do

    ! resize the meshblocks
    meshblocks(:) = new_meshblocks(:)

    ! deallocate buffers and redistribute particles
    do q = 1, nproc_group
      left_rnk = left_group(q)
      right_rnk = right_group(q)
      if ((mpi_rank .eq. left_rnk) .or. (mpi_rank .eq. right_rnk)) then
        call deallocateFieldBackups()

        ! shift particles
        if (mpi_rank .eq. right_rnk) then
          call shiftParticlesX(-SHIFT)
        end if

        ! backup particles
        call backupParticles()

        ! reshuffle particle tiles
        call reallocateParticles(this_meshblock%ptr)
        ! restore particles from backup
        call restoreParticlesFromBackup()
        call deallocateParticleBackup()

      end if
    end do

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    ! put particles back on proper meshblocks
    call redistributeParticlesBetweenMeshblocks()
    call clearGhostParticles()

    call printDiag("reshapeInX()", 3)
  end subroutine reshapeInX

  subroutine reshapeInY(left_group, right_group, SHIFT)
    integer, allocatable, intent(in)    :: left_group(:), right_group(:)
    integer, intent(in)                 :: SHIFT
    integer                       :: nproc_group, q, left_rnk, right_rnk, ierr
    integer                       :: new_sx, new_sy, new_sz
    integer                       :: i1_from, i2_from, j1_from, j2_from, k1_from, k2_from,&
                                   & i1_to, i2_to, j1_to, j2_to, k1_to, k2_to
    nproc_group = size(left_group)
    #ifdef DEBUG
      if (size(left_group) .ne. size(right_group)) then
        call throwError('ERROR: wrong groups specified for `reshapeInY`.')
      end if
    #endif

    ! backup the fields with current sizes
    do q = 1, nproc_group
      left_rnk = left_group(q)
      right_rnk = right_group(q)
      if ((mpi_rank .eq. left_rnk) .or. (mpi_rank .eq. right_rnk)) then
        call backupEBfields()
      end if
    end do

    ! get new meshblock dimensions
    new_meshblocks(:) = meshblocks(:)
    call reassignNeighborsForAll(new_meshblocks)
    do q = 1, nproc_group
      left_rnk = left_group(q)
      right_rnk = right_group(q)
      new_meshblocks(left_rnk + 1)%sy = new_meshblocks(left_rnk + 1)%sy + SHIFT
      new_meshblocks(right_rnk + 1)%sy = new_meshblocks(right_rnk + 1)%sy - SHIFT
      new_meshblocks(right_rnk + 1)%y0 = new_meshblocks(right_rnk + 1)%y0 + SHIFT
    end do
    ! at this point DO NOT CHANGE `meshblocks` ...
    ! ... as the `exchangeFieldSlabIn*` still assumes old dimensions

    ! reallocate field arrays given the new meshblock dimensions
    do q = 1, nproc_group
      left_rnk = left_group(q)
      right_rnk = right_group(q)
      if ((mpi_rank .eq. left_rnk) .or. (mpi_rank .eq. right_rnk)) then
        call deallocateFields()
        call reallocateFields(new_meshblocks(mpi_rank + 1))
        call reallocateFieldBuffers(new_meshblocks(mpi_rank + 1))
      end if
    end do

    ! send/recv missing fields
    ! ... and recover from backup
    do q = 1, nproc_group
      left_rnk = left_group(q)
      right_rnk = right_group(q)
      if ((mpi_rank .eq. left_rnk) .or. (mpi_rank .eq. right_rnk)) then
        ! exchange slab
        call exchangeFieldSlabInY(left_rnk, right_rnk, SHIFT)

        ! recover from backup
        if (SHIFT .gt. 0) then
          ! `left_rnk` inflates
          ! `right_rnk` shrinks
          if (mpi_rank .eq. left_rnk) then
            i1_to = 0; i2_to = this_meshblock%ptr%sx - 1
            j1_to = 0; j2_to = this_meshblock%ptr%sy - 1
            k1_to = 0; k2_to = this_meshblock%ptr%sz - 1
            i1_from = i1_to; i2_from = i2_to
            j1_from = j1_to; j2_from = j2_to
            k1_from = k1_to; k2_from = k2_to
          else if (mpi_rank .eq. right_rnk) then
            i1_to = 0; i2_to = new_meshblocks(right_rnk + 1)%sx - 1
            j1_to = 0; j2_to = new_meshblocks(right_rnk + 1)%sy - 1
            k1_to = 0; k2_to = new_meshblocks(right_rnk + 1)%sz - 1
            i1_from = 0; i2_from = this_meshblock%ptr%sx - 1
            j1_from = SHIFT; j2_from = this_meshblock%ptr%sy - 1
            k1_from = 0; k2_from = this_meshblock%ptr%sz - 1
          end if
        else
          ! `left_rnk` shrinks
          ! `right_rnk` inflates
          if (mpi_rank .eq. left_rnk) then
            i1_to = 0; i2_to = this_meshblock%ptr%sx - 1
            j1_to = 0; j2_to = this_meshblock%ptr%sy - 1 - abs(SHIFT)
            k1_to = 0; k2_to = this_meshblock%ptr%sz - 1
            i1_from = i1_to; i2_from = i2_to
            j1_from = j1_to; j2_from = j2_to
            k1_from = k1_to; k2_from = k2_to
          else if (mpi_rank .eq. right_rnk) then
            i1_to = 0; i2_to = new_meshblocks(right_rnk + 1)%sx - 1
            j1_to = abs(SHIFT); j2_to = new_meshblocks(right_rnk + 1)%sy - 1
            k1_to = 0; k2_to = new_meshblocks(right_rnk + 1)%sz - 1
            i1_from = 0; i2_from = this_meshblock%ptr%sx - 1
            j1_from = 0; j2_from = this_meshblock%ptr%sy - 1
            k1_from = 0; k2_from = this_meshblock%ptr%sz - 1
          end if
        end if

        call restoreFieldsFromBackups(i1_from, i2_from, j1_from, j2_from, k1_from, k2_from,&
                                    & i1_to, i2_to, j1_to, j2_to, k1_to, k2_to)
      end if
    end do

    ! resize the meshblocks
    meshblocks(:) = new_meshblocks(:)

    ! deallocate buffers and redistribute particles
    do q = 1, nproc_group
      left_rnk = left_group(q)
      right_rnk = right_group(q)
      if ((mpi_rank .eq. left_rnk) .or. (mpi_rank .eq. right_rnk)) then
        call deallocateFieldBackups()

        ! shift particles
        if (mpi_rank .eq. right_rnk) then
          call shiftParticlesY(-SHIFT)
        end if

        ! backup particles
        call backupParticles()

        ! reshuffle particle tiles
        call reallocateParticles(this_meshblock%ptr)
        ! restore particles from backup
        call restoreParticlesFromBackup()
        call deallocateParticleBackup()

      end if
    end do

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    ! put particles back on proper meshblocks
    call redistributeParticlesBetweenMeshblocks()
    call clearGhostParticles()

    call printDiag("reshapeInY()", 3)
  end subroutine reshapeInY

  subroutine reshapeInZ(left_group, right_group, SHIFT)
    integer, allocatable, intent(in)    :: left_group(:), right_group(:)
    integer, intent(in)                 :: SHIFT
    integer                       :: nproc_group, q, left_rnk, right_rnk, ierr
    integer                       :: new_sx, new_sy, new_sz
    integer                       :: i1_from, i2_from, j1_from, j2_from, k1_from, k2_from,&
                                   & i1_to, i2_to, j1_to, j2_to, k1_to, k2_to
    nproc_group = size(left_group)
    #ifdef DEBUG
      if (size(left_group) .ne. size(right_group)) then
        call throwError('ERROR: wrong groups specified for `reshapeInZ`.')
      end if
    #endif

    ! backup the fields with current sizes
    do q = 1, nproc_group
      left_rnk = left_group(q)
      right_rnk = right_group(q)
      if ((mpi_rank .eq. left_rnk) .or. (mpi_rank .eq. right_rnk)) then
        call backupEBfields()
      end if
    end do

    ! get new meshblock dimensions
    new_meshblocks(:) = meshblocks(:)
    call reassignNeighborsForAll(new_meshblocks)
    do q = 1, nproc_group
      left_rnk = left_group(q)
      right_rnk = right_group(q)
      new_meshblocks(left_rnk + 1)%sz = new_meshblocks(left_rnk + 1)%sz + SHIFT
      new_meshblocks(right_rnk + 1)%sz = new_meshblocks(right_rnk + 1)%sz - SHIFT
      new_meshblocks(right_rnk + 1)%z0 = new_meshblocks(right_rnk + 1)%z0 + SHIFT
    end do
    ! at this point DO NOT CHANGE `meshblocks` ...
    ! ... as the `exchangeFieldSlabIn*` still assumes old dimensions

    ! reallocate field arrays given the new meshblock dimensions
    do q = 1, nproc_group
      left_rnk = left_group(q)
      right_rnk = right_group(q)
      if ((mpi_rank .eq. left_rnk) .or. (mpi_rank .eq. right_rnk)) then
        call deallocateFields()
        call reallocateFields(new_meshblocks(mpi_rank + 1))
        call reallocateFieldBuffers(new_meshblocks(mpi_rank + 1))
      end if
    end do

    ! send/recv missing fields
    ! ... and recover from backup
    do q = 1, nproc_group
      left_rnk = left_group(q)
      right_rnk = right_group(q)
      if ((mpi_rank .eq. left_rnk) .or. (mpi_rank .eq. right_rnk)) then
        ! exchange slab
        call exchangeFieldSlabInZ(left_rnk, right_rnk, SHIFT)

        ! recover from backup
        if (SHIFT .gt. 0) then
          ! `left_rnk` inflates
          ! `right_rnk` shrinks
          if (mpi_rank .eq. left_rnk) then
            i1_to = 0; i2_to = this_meshblock%ptr%sx - 1
            j1_to = 0; j2_to = this_meshblock%ptr%sy - 1
            k1_to = 0; k2_to = this_meshblock%ptr%sz - 1
            i1_from = i1_to; i2_from = i2_to
            j1_from = j1_to; j2_from = j2_to
            k1_from = k1_to; k2_from = k2_to
          else if (mpi_rank .eq. right_rnk) then
            i1_to = 0; i2_to = new_meshblocks(right_rnk + 1)%sx - 1
            j1_to = 0; j2_to = new_meshblocks(right_rnk + 1)%sy - 1
            k1_to = 0; k2_to = new_meshblocks(right_rnk + 1)%sz - 1
            i1_from = 0; i2_from = this_meshblock%ptr%sx - 1
            j1_from = 0; j2_from = this_meshblock%ptr%sy - 1
            k1_from = SHIFT; k2_from = this_meshblock%ptr%sz - 1
          end if
        else
          ! `left_rnk` shrinks
          ! `right_rnk` inflates
          if (mpi_rank .eq. left_rnk) then
            i1_to = 0; i2_to = this_meshblock%ptr%sx - 1
            j1_to = 0; j2_to = this_meshblock%ptr%sy - 1
            k1_to = 0; k2_to = this_meshblock%ptr%sz - 1 - abs(SHIFT)
            i1_from = i1_to; i2_from = i2_to
            j1_from = j1_to; j2_from = j2_to
            k1_from = k1_to; k2_from = k2_to
          else if (mpi_rank .eq. right_rnk) then
            i1_to = 0; i2_to = new_meshblocks(right_rnk + 1)%sx - 1
            j1_to = 0; j2_to = new_meshblocks(right_rnk + 1)%sy - 1
            k1_to = abs(SHIFT); k2_to = new_meshblocks(right_rnk + 1)%sz - 1
            i1_from = 0; i2_from = this_meshblock%ptr%sx - 1
            j1_from = 0; j2_from = this_meshblock%ptr%sy - 1
            k1_from = 0; k2_from = this_meshblock%ptr%sz - 1
          end if
        end if

        call restoreFieldsFromBackups(i1_from, i2_from, j1_from, j2_from, k1_from, k2_from,&
                                    & i1_to, i2_to, j1_to, j2_to, k1_to, k2_to)
      end if
    end do

    ! resize the meshblocks
    meshblocks(:) = new_meshblocks(:)

    ! deallocate buffers and redistribute particles
    do q = 1, nproc_group
      left_rnk = left_group(q)
      right_rnk = right_group(q)
      if ((mpi_rank .eq. left_rnk) .or. (mpi_rank .eq. right_rnk)) then
        call deallocateFieldBackups()

        ! shift particles
        if (mpi_rank .eq. right_rnk) then
          call shiftParticlesZ(-SHIFT)
        end if

        ! backup particles
        call backupParticles()

        ! reshuffle particle tiles
        call reallocateParticles(this_meshblock%ptr)
        ! restore particles from backup
        call restoreParticlesFromBackup()
        call deallocateParticleBackup()

      end if
    end do

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    ! put particles back on proper meshblocks
    call redistributeParticlesBetweenMeshblocks()
    call clearGhostParticles()

    call printDiag("reshapeInZ()", 3)
  end subroutine reshapeInZ

end module m_adaptivelb
