#include "../defs.F90"

module m_exchangeparts
  use m_globalnamespace
  use m_readinput, only: getInput
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_particlelogistics
  !--- PRIVATE functions -----------------------------------------!
  !...............................................................!
contains
  subroutine initializePrtlExchange()
    implicit none
    integer             :: ierr, additional_real, additional_int, additional_int2

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

    call reallocateEnrouteArray(this_meshblock%ptr)

    ! DEP_PRT [particle-dependent]
    ! new type for myMPI_ENROUTE
    additional_real = 0; additional_int = 0; additional_int2 = 0

    call MPI_TYPE_GET_EXTENT(MPI_INTEGER2, lb, extent_int2, ierr)
    call MPI_TYPE_GET_EXTENT(MPI_REAL, lb, extent_real, ierr)

#ifdef GCA
    additional_int2 = additional_int2 + 3
    additional_real = additional_real + 8
#endif

#ifdef PRTLPAYLOADS
    additional_real = additional_real + 3
#endif

    !     # of blockcounts = 3:
    !       3  x integer2  [xi, yi, zi]                       | + 3 if GCA [xi_past, yi_past, zi_past]
    !       7  x real      [dx, dy, dz, u, v, w, weight]      | + 8 if GCA [dx_past, dy_past, dz_past, u_eff, v_eff, w_eff, u_par, u_perp]
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

    call printDiag("initializePrtlExchange()", 1)
  end subroutine initializePrtlExchange

#ifndef MPINONBLOCK
  ! blocking MPI communication
  subroutine exchangeParticles()
    implicit none
    integer(kind=2), pointer, contiguous  :: pt_xi(:), pt_yi(:), pt_zi(:)
    integer, pointer, contiguous          :: pt_proc(:)
    real, pointer, contiguous             :: pt_dx(:), pt_dy(:), pt_dz(:)
    integer             :: s, p, send_x, send_y, send_z, ti, tj, tk, ti_p, tj_p, tk_p
    integer             :: mpi_sendto, mpi_recvfrom, mpi_sendtag, mpi_recvtag, mpi_tag, mpi_tag2
    integer             :: ierr, ind1, ind2, ind3, cntr, temp_xyz
    integer             :: cnt_recv_enroute, cnt_send_enroute
    logical :: should_send, should_recv
    integer(kind=2)     :: new_xyz

#ifdef MPI08
    type(MPI_REQUEST), allocatable  :: mpi_req(:)
    type(MPI_STATUS)                :: istat
#endif

#ifdef MPI
    integer, allocatable            :: mpi_req(:)
    integer                         :: istat(MPI_STATUS_SIZE)
#endif

    logical, allocatable              :: mpi_sendflags(:), mpi_recvflags(:)
    logical                           :: quit_loop

    allocate(mpi_req(sendrecv_neighbors))
    allocate(mpi_sendflags(sendrecv_neighbors))
    allocate(mpi_recvflags(sendrecv_neighbors))

    do s = 1, nspec ! loop over species

#ifdef LOWMEM
      enroute_bot%get(:,:,:)%cnt = 0
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
            pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
            pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi
            pt_dx => species(s)%prtl_tile(ti, tj, tk)%dx
            pt_dy => species(s)%prtl_tile(ti, tj, tk)%dy
            pt_dz => species(s)%prtl_tile(ti, tj, tk)%dz
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              ! send_* = -1 / 0 / +1
              send_x = 0; send_y = 0; send_z = 0
#if defined(oneD) || defined(twoD) || defined(threeD)
              send_x = (ISIGN(1, pt_xi(p) - this_meshblock%ptr%sx) + 1) / 2 - (ISIGN(1, -pt_xi(p) - 1) + 1) / 2
#endif
#if defined(twoD) || defined(threeD)
              send_y = (ISIGN(1, pt_yi(p) - this_meshblock%ptr%sy) + 1) / 2 - (ISIGN(1, -pt_yi(p) - 1) + 1) / 2
#endif
#if defined(threeD)
              send_z = (ISIGN(1, pt_zi(p) - this_meshblock%ptr%sz) + 1) / 2 - (ISIGN(1, -pt_zi(p) - 1) + 1) / 2
#endif
              if ((send_x .ne. 0) .or. (send_y .ne. 0) .or. (send_z .ne. 0)) then
                if (.not. associated(this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr)) then
                  cycle
                end if
                ! count the # of particles to be sent in each direction
                enroute_bot%get(send_x, send_y, send_z)%cnt = enroute_bot%get(send_x, send_y, send_z)%cnt + 1
              end if
            end do ! particles
            pt_xi => null(); pt_yi => null(); pt_zi => null()
            pt_dx => null(); pt_dy => null(); pt_dz => null()
          end do
        end do
      end do
      ! allocate the corresponding buffer arrays for each direction
      do send_x = -1, 1
        do send_y = -1, 1
          do send_z = -1, 1
            call reallocateEnroute(send_x, send_y, send_z, enroute_bot%get(send_x, send_y, send_z)%cnt + 1)
          end do
        end do
      end do
#endif ! LOWMEM

      enroute_bot%get(:,:,:)%cnt = 0
      ! particle crosses MPI blocks //
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
            pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
            pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi

            pt_dx => species(s)%prtl_tile(ti, tj, tk)%dx
            pt_dy => species(s)%prtl_tile(ti, tj, tk)%dy
            pt_dz => species(s)%prtl_tile(ti, tj, tk)%dz

            pt_proc => species(s)%prtl_tile(ti, tj, tk)%proc
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              ! send_* = -1 / 0 / +1
              send_x = 0; send_y = 0; send_z = 0
#if defined(oneD) || defined(twoD) || defined(threeD)
              send_x = (ISIGN(1, pt_xi(p) - this_meshblock%ptr%sx) + 1) / 2 - (ISIGN(1, -pt_xi(p) - 1) + 1) / 2
#endif
#if defined(twoD) || defined(threeD)
              send_y = (ISIGN(1, pt_yi(p) - this_meshblock%ptr%sy) + 1) / 2 - (ISIGN(1, -pt_yi(p) - 1) + 1) / 2
#endif
#if defined(threeD)
              send_z = (ISIGN(1, pt_zi(p) - this_meshblock%ptr%sz) + 1) / 2 - (ISIGN(1, -pt_zi(p) - 1) + 1) / 2
#endif
              ! FIX1 check for null() boundaries
              if ((send_x .ne. 0) .or. (send_y .ne. 0) .or. (send_z .ne. 0)) then
                if (.not. associated(this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr)) then
                  ! make ghost particle
                  pt_proc(p) = -1
                  cycle
                end if
                ! copy this particle to temporary `enroute_bot` array
                enroute_bot%get(send_x, send_y, send_z)%cnt = enroute_bot%get(send_x, send_y, send_z)%cnt + 1
                if (enroute_bot%get(send_x, send_y, send_z)%cnt .ge.&
                & enroute_bot%get(send_x, send_y, send_z)%max) then
                  print *, "RANK: ", mpi_rank
                  print *, "MBLOCK: ", this_meshblock%ptr%x0, this_meshblock%ptr%y0
                  print *, "SPEC: ", s
                  print *, "DIR: ", send_x, send_y, send_z
                  call throwError('ERROR: particle send buffer array too small: '//&
                  & trim(STR(enroute_bot%get(send_x, send_y, send_z)%max))//'.')
                end if
                cntr = enroute_bot%get(send_x, send_y, send_z)%cnt
                call copyToEnroute(s, ti, tj, tk, p, enroute_bot%get(send_x, send_y, send_z)%enroute(cntr))
                ! make ghost particle
                pt_proc(p) = -1

                ! shift coordinates to fit the new grid
#if defined(oneD) || defined(twoD) || defined(threeD)
                new_xyz = enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi
                temp_xyz = this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sx
                new_xyz = -(send_x - 1) * (2 + send_x) * (new_xyz * (send_x + 1) - (temp_xyz - 1) * send_x) / 2
                enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi = new_xyz
#endif

#if defined(twoD) || defined(threeD)
                new_xyz = enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi
                temp_xyz = this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sy
                new_xyz = -(send_y - 1) * (2 + send_y) * (new_xyz * (send_y + 1) - (temp_xyz - 1) * send_y) / 2
                enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi = new_xyz
#endif

#if defined(threeD)
                new_xyz = enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi
                temp_xyz = this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sz
                new_xyz = -(send_z - 1) * (2 + send_z) * (new_xyz * (send_z + 1) - (temp_xyz - 1) * send_z) / 2
                enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi = new_xyz
#endif

#ifdef GCA
                ! this thing below can probably be done better
                ! shift past coordinates to fit the new grid
#if defined(oneD) || defined(twoD) || defined(threeD)
                if (send_x .eq. 1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi_past =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi_past -&
                  & this_meshblock%ptr%sx
                else if (send_x .eq. -1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi_past =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi_past +&
                  & this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sx
                end if
#endif

#if defined(twoD) || defined(threeD)
                if (send_y .eq. 1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi_past =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi_past -&
                  & this_meshblock%ptr%sy
                else if (send_y .eq. -1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi_past =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi_past +&
                  & this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sy
                end if
#endif

#if defined(threeD)
                if (send_z .eq. 1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi_past =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi_past -&
                  & this_meshblock%ptr%sz
                else if (send_z .eq. -1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi_past =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi_past +&
                  & this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sz
                end if
#endif
#endif
              end if
            end do ! particles
            pt_xi => null(); pt_yi => null(); pt_zi => null()
            pt_dx => null(); pt_dy => null(); pt_dz => null()
            pt_proc => null()
          end do ! tk
        end do ! tj
      end do ! ti
      ! // particle crosses MPI blocks

      cntr = 0
      do ind1 = -1, 1
        do ind2 = -1, 1
          do ind3 = -1, 1
            if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
#ifdef oneD
            if ((ind2 .ne. 0) .or. (ind3 .ne. 0)) cycle
#elif twoD
            if (ind3 .ne. 0) cycle
#endif
            cntr = cntr + 1
            mpi_tag = (ind3 + 2) + 3 * (ind2 + 1) + 9 * (ind1 + 1)
            mpi_tag2 = (ind3 + 2) + 3 * (ind2 + 1) + 9 * (ind1 + 1)+200
            should_send = associated(this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr)
            should_recv = associated(this_meshblock%ptr%neighbor(-ind1,-ind2,-ind3)%ptr)
            cnt_send_enroute = enroute_bot%get(ind1,ind2,ind3)%cnt
            if (should_send .and. should_recv) then
              mpi_sendto = this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr%rnk
              mpi_recvfrom = this_meshblock%ptr%neighbor(-ind1,-ind2,-ind3)%ptr%rnk

              call MPI_SENDRECV(cnt_send_enroute, 1, MPI_INTEGER, mpi_sendto, mpi_tag2,&
              & cnt_recv_enroute, 1, MPI_INTEGER, mpi_recvfrom, mpi_tag2, MPI_COMM_WORLD, istat, ierr)
              if (cnt_recv_enroute .ge. size(recv_enroute%enroute)) then
                call throwError('ERROR: particle had rcv array too small.')
              end if
              call MPI_SENDRECV(enroute_bot%get(ind1,ind2,ind3)%enroute(1:enroute_bot%get(ind1,ind2,ind3)%cnt),&
              & enroute_bot%get(ind1,ind2,ind3)%cnt, myMPI_ENROUTE,&
              & mpi_sendto, mpi_tag,recv_enroute%enroute(1:cnt_recv_enroute), &
              & cnt_recv_enroute, myMPI_ENROUTE, mpi_recvfrom, mpi_tag, MPI_COMM_WORLD, istat, ierr)
              if (cnt_recv_enroute .gt. 0) then
                call extractParticlesFromEnroute(cnt_recv_enroute, s)
              end if ! if > 0 particles received
            else if ((.not. should_send) .and. should_recv) then
              ! proc # 8
              mpi_recvfrom = this_meshblock%ptr%neighbor(-ind1,-ind2,-ind3)%ptr%rnk ! 7
              call MPI_RECV(cnt_recv_enroute, 1, MPI_INTEGER, mpi_recvfrom, mpi_tag2, MPI_COMM_WORLD, istat, ierr)
              call MPI_RECV(recv_enroute%enroute(1:cnt_recv_enroute),cnt_recv_enroute, myMPI_ENROUTE,&
              & mpi_recvfrom, mpi_tag, MPI_COMM_WORLD, istat, ierr)
              if (cnt_recv_enroute .gt. 0) then
                call extractParticlesFromEnroute(cnt_recv_enroute, s)
              end if
            else if ((.not. should_recv) .and. should_send) then
              mpi_sendto = this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr%rnk
              call MPI_SEND(cnt_send_enroute, 1, MPI_INTEGER, mpi_sendto, mpi_tag2, MPI_COMM_WORLD, ierr)
              call MPI_SEND(enroute_bot%get(ind1,ind2,ind3)%enroute(1:enroute_bot%get(ind1,ind2,ind3)%cnt),&
              & enroute_bot%get(ind1,ind2,ind3)%cnt, myMPI_ENROUTE,&
              & mpi_sendto, mpi_tag, MPI_COMM_WORLD, ierr)
            end if
          end do
        end do
      end do

      ! particle moves between tiles within a single MPI block //
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
            pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
            pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi
            pt_proc => species(s)%prtl_tile(ti, tj, tk)%proc
            ! FIX1 make sure this is vectorized
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              if (pt_proc(p) .eq. -1) cycle
              ti_p = FLOOR(REAL(pt_xi(p)) / REAL(species(s)%tile_sx)) + 1
              tj_p = FLOOR(REAL(pt_yi(p)) / REAL(species(s)%tile_sy)) + 1
              tk_p = FLOOR(REAL(pt_zi(p)) / REAL(species(s)%tile_sz)) + 1
              if ((ti_p .ne. ti) .or. (tj_p .ne. tj) .or. (tk_p .ne. tk)) then
                call moveParticleBetweenTiles(s, ti, tj, tk, p)
              end if
            end do ! particles
            pt_xi => null(); pt_yi => null(); pt_zi => null()
            pt_proc => null()
          end do ! tk
        end do ! tj
      end do ! ti
      ! // particle moves between tiles within a single MPI block

#ifdef DEBUG
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
            pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
            pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi
            pt_proc => species(s)%prtl_tile(ti, tj, tk)%proc
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              if (pt_proc(p) .eq. -1) cycle
              if ((pt_xi(p) .lt. species(s)%prtl_tile(ti, tj, tk)%x1) .or. &
              & (pt_xi(p) .ge. species(s)%prtl_tile(ti, tj, tk)%x2) .or. &
              & (pt_yi(p) .lt. species(s)%prtl_tile(ti, tj, tk)%y1) .or. &
              & (pt_yi(p) .ge. species(s)%prtl_tile(ti, tj, tk)%y2) .or. &
              & (pt_zi(p) .lt. species(s)%prtl_tile(ti, tj, tk)%z1) .or. &
              & (pt_zi(p) .ge. species(s)%prtl_tile(ti, tj, tk)%z2)) then
                call throwError('ERROR: particle in wrong tile after exchange')
              end if
            end do
          end do
        end do
      end do
#endif

#ifdef LOWMEM
      do send_x = -1, 1
        do send_y = -1, 1
          do send_z = -1, 1
            if (allocated(enroute_bot%get(send_x, send_y, send_z)%enroute)) then
              deallocate(enroute_bot%get(send_x, send_y, send_z)%enroute)
            end if
          end do
        end do
      end do
#endif

    end do ! loop over species
    call printDiag("exchangeParticles()", 2)
  end subroutine exchangeParticles
#else
  ! non-blocking MPI communication
  subroutine exchangeParticles()
    implicit none
    integer(kind=2), pointer, contiguous  :: pt_xi(:), pt_yi(:), pt_zi(:)
    integer, pointer, contiguous          :: pt_proc(:)
    real, pointer, contiguous             :: pt_dx(:), pt_dy(:), pt_dz(:)
    integer             :: s, p, send_x, send_y, send_z, ti, tj, tk, ti_p, tj_p, tk_p
    integer             :: mpi_sendto, mpi_recvfrom, mpi_sendtag, mpi_recvtag
    integer             :: ierr, ind1, ind2, ind3, cntr, temp_xyz
    integer             :: cnt_recv_enroute
    integer(kind=2)     :: new_xyz

#ifdef MPI08
    type(MPI_REQUEST), allocatable  :: mpi_req(:)
    type(MPI_STATUS)                :: istat
#endif

#ifdef MPI
    integer, allocatable            :: mpi_req(:)
    integer                         :: istat(MPI_STATUS_SIZE)
#endif

    logical, allocatable              :: mpi_sendflags(:), mpi_recvflags(:)
    logical                           :: quit_loop

    allocate(mpi_req(sendrecv_neighbors))
    allocate(mpi_sendflags(sendrecv_neighbors))
    allocate(mpi_recvflags(sendrecv_neighbors))

    do s = 1, nspec ! loop over species

#ifdef LOWMEM
      enroute_bot%get(:,:,:)%cnt = 0
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
            pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
            pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi
            pt_dx => species(s)%prtl_tile(ti, tj, tk)%dx
            pt_dy => species(s)%prtl_tile(ti, tj, tk)%dy
            pt_dz => species(s)%prtl_tile(ti, tj, tk)%dz
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              ! send_* = -1 / 0 / +1
              send_x = 0; send_y = 0; send_z = 0
#if defined(oneD) || defined(twoD) || defined(threeD)
              send_x = (ISIGN(1, pt_xi(p) - this_meshblock%ptr%sx) + 1) / 2 - (ISIGN(1, -pt_xi(p) - 1) + 1) / 2
#endif
#if defined(twoD) || defined(threeD)
              send_y = (ISIGN(1, pt_yi(p) - this_meshblock%ptr%sy) + 1) / 2 - (ISIGN(1, -pt_yi(p) - 1) + 1) / 2
#endif
#if defined(threeD)
              send_z = (ISIGN(1, pt_zi(p) - this_meshblock%ptr%sz) + 1) / 2 - (ISIGN(1, -pt_zi(p) - 1) + 1) / 2
#endif
              if ((send_x .ne. 0) .or. (send_y .ne. 0) .or. (send_z .ne. 0)) then
                if (.not. associated(this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr)) then
                  cycle
                end if
                ! count the # of particles to be sent in each direction
                enroute_bot%get(send_x, send_y, send_z)%cnt = enroute_bot%get(send_x, send_y, send_z)%cnt + 1
              end if
            end do ! particles
            pt_xi => null(); pt_yi => null(); pt_zi => null()
            pt_dx => null(); pt_dy => null(); pt_dz => null()
          end do
        end do
      end do
      ! allocate the corresponding buffer arrays for each direction
      do send_x = -1, 1
        do send_y = -1, 1
          do send_z = -1, 1
            call reallocateEnroute(send_x, send_y, send_z, enroute_bot%get(send_x, send_y, send_z)%cnt + 1)
          end do
        end do
      end do
#endif ! LOWMEM

      enroute_bot%get(:,:,:)%cnt = 0
      ! particle crosses MPI blocks //
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
            pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
            pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi

            pt_dx => species(s)%prtl_tile(ti, tj, tk)%dx
            pt_dy => species(s)%prtl_tile(ti, tj, tk)%dy
            pt_dz => species(s)%prtl_tile(ti, tj, tk)%dz

            pt_proc => species(s)%prtl_tile(ti, tj, tk)%proc
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              ! send_* = -1 / 0 / +1
              send_x = 0; send_y = 0; send_z = 0
#if defined(oneD) || defined(twoD) || defined(threeD)
              send_x = (ISIGN(1, pt_xi(p) - this_meshblock%ptr%sx) + 1) / 2 - (ISIGN(1, -pt_xi(p) - 1) + 1) / 2
#endif
#if defined(twoD) || defined(threeD)
              send_y = (ISIGN(1, pt_yi(p) - this_meshblock%ptr%sy) + 1) / 2 - (ISIGN(1, -pt_yi(p) - 1) + 1) / 2
#endif
#if defined(threeD)
              send_z = (ISIGN(1, pt_zi(p) - this_meshblock%ptr%sz) + 1) / 2 - (ISIGN(1, -pt_zi(p) - 1) + 1) / 2
#endif
              ! FIX1 check for null() boundaries
              if ((send_x .ne. 0) .or. (send_y .ne. 0) .or. (send_z .ne. 0)) then
                if (.not. associated(this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr)) then
                  ! make ghost particle
                  pt_proc(p) = -1
                  cycle
                end if
                ! copy this particle to temporary `enroute_bot` array
                enroute_bot%get(send_x, send_y, send_z)%cnt = enroute_bot%get(send_x, send_y, send_z)%cnt + 1
                if (enroute_bot%get(send_x, send_y, send_z)%cnt .ge.&
                & enroute_bot%get(send_x, send_y, send_z)%max) then
                  print *, "RANK: ", mpi_rank
                  print *, "MBLOCK: ", this_meshblock%ptr%x0, this_meshblock%ptr%y0
                  print *, "SPEC: ", s
                  print *, "DIR: ", send_x, send_y, send_z
                  call throwError('ERROR: particle send buffer array too small: '//&
                  & trim(STR(enroute_bot%get(send_x, send_y, send_z)%max))//'.')
                end if
                cntr = enroute_bot%get(send_x, send_y, send_z)%cnt
                call copyToEnroute(s, ti, tj, tk, p, enroute_bot%get(send_x, send_y, send_z)%enroute(cntr))
                ! make ghost particle
                pt_proc(p) = -1

                ! shift coordinates to fit the new grid
#if defined(oneD) || defined(twoD) || defined(threeD)
                new_xyz = enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi
                temp_xyz = this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sx
                new_xyz = -(send_x - 1) * (2 + send_x) * (new_xyz * (send_x + 1) - (temp_xyz - 1) * send_x) / 2
                enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi = new_xyz
#endif

#if defined(twoD) || defined(threeD)
                new_xyz = enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi
                temp_xyz = this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sy
                new_xyz = -(send_y - 1) * (2 + send_y) * (new_xyz * (send_y + 1) - (temp_xyz - 1) * send_y) / 2
                enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi = new_xyz
#endif

#if defined(threeD)
                new_xyz = enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi
                temp_xyz = this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sz
                new_xyz = -(send_z - 1) * (2 + send_z) * (new_xyz * (send_z + 1) - (temp_xyz - 1) * send_z) / 2
                enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi = new_xyz
#endif

#ifdef GCA
                ! this thing below can probably be done better
                ! shift past coordinates to fit the new grid
#if defined(oneD) || defined(twoD) || defined(threeD)
                if (send_x .eq. 1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi_past =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi_past -&
                  & this_meshblock%ptr%sx
                else if (send_x .eq. -1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi_past =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi_past +&
                  & this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sx
                end if
#endif

#if defined(twoD) || defined(threeD)
                if (send_y .eq. 1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi_past =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi_past -&
                  & this_meshblock%ptr%sy
                else if (send_y .eq. -1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi_past =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi_past +&
                  & this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sy
                end if
#endif

#if defined(threeD)
                if (send_z .eq. 1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi_past =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi_past -&
                  & this_meshblock%ptr%sz
                else if (send_z .eq. -1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi_past =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi_past +&
                  & this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sz
                end if
#endif
#endif
              end if
            end do ! particles
            pt_xi => null(); pt_yi => null(); pt_zi => null()
            pt_dx => null(); pt_dy => null(); pt_dz => null()
            pt_proc => null()
          end do ! tk
        end do ! tj
      end do ! ti
      ! // particle crosses MPI blocks

      ! start sending //
      cntr = 0
      do ind1 = -1, 1
        do ind2 = -1, 1
          do ind3 = -1, 1
            if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
#ifdef oneD
            if ((ind2 .ne. 0) .or. (ind3 .ne. 0)) cycle
#elif twoD
            if (ind3 .ne. 0) cycle
#endif
            if (.not. associated(this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr)) cycle
            cntr = cntr + 1
            mpi_sendtag = (ind3 + 2) + 3 * (ind2 + 1) + 9 * (ind1 + 1)

            ! post non-blocking send requests
            mpi_sendto = this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr%rnk
            call MPI_ISEND(enroute_bot%get(ind1,ind2,ind3)%enroute(1:enroute_bot%get(ind1,ind2,ind3)%cnt),&
            & enroute_bot%get(ind1,ind2,ind3)%cnt, myMPI_ENROUTE,&
            & mpi_sendto, mpi_sendtag, MPI_COMM_WORLD, mpi_req(cntr), ierr)
          end do
        end do
      end do
      ! // start sending

      ! particle moves between tiles within a single MPI block //
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
            pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
            pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi
            pt_proc => species(s)%prtl_tile(ti, tj, tk)%proc
            ! FIX1 make sure this is vectorized
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              if (pt_proc(p) .eq. -1) cycle
              ti_p = FLOOR(REAL(pt_xi(p)) / REAL(species(s)%tile_sx)) + 1
              tj_p = FLOOR(REAL(pt_yi(p)) / REAL(species(s)%tile_sy)) + 1
              tk_p = FLOOR(REAL(pt_zi(p)) / REAL(species(s)%tile_sz)) + 1
              if ((ti_p .ne. ti) .or. (tj_p .ne. tj) .or. (tk_p .ne. tk)) then
                call moveParticleBetweenTiles(s, ti, tj, tk, p)
              end if
            end do ! particles
            pt_xi => null(); pt_yi => null(); pt_zi => null()
            pt_proc => null()
          end do ! tk
        end do ! tj
      end do ! ti
      ! // particle moves between tiles within a single MPI block

#ifdef DEBUG
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
            pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
            pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi
            pt_proc => species(s)%prtl_tile(ti, tj, tk)%proc
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              if (pt_proc(p) .eq. -1) cycle
              if ((pt_xi(p) .lt. species(s)%prtl_tile(ti, tj, tk)%x1) .or. &
              & (pt_xi(p) .ge. species(s)%prtl_tile(ti, tj, tk)%x2) .or. &
              & (pt_yi(p) .lt. species(s)%prtl_tile(ti, tj, tk)%y1) .or. &
              & (pt_yi(p) .ge. species(s)%prtl_tile(ti, tj, tk)%y2) .or. &
              & (pt_zi(p) .lt. species(s)%prtl_tile(ti, tj, tk)%z1) .or. &
              & (pt_zi(p) .ge. species(s)%prtl_tile(ti, tj, tk)%z2)) then
                call throwError('ERROR: particle in wrong tile after exchange')
              end if
            end do
          end do
        end do
      end do
#endif

      ! wait to send & receive all the MPI calls and write data to memory
      quit_loop = .false.
      mpi_sendflags(:) = .false.
      mpi_recvflags(:) = .false.
      do while (.not. quit_loop)
        ! try to receive while not done
        quit_loop = .true.
        cntr = 0
        do ind1 = -1, 1
          do ind2 = -1, 1
            do ind3 = -1, 1
              if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
#ifdef oneD
              if ((ind2 .ne. 0) .or. (ind3 .ne. 0)) cycle
#elif twoD
              if (ind3 .ne. 0) cycle
#endif
              if (.not. associated(this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr)) cycle
              cntr = cntr + 1

              if (.not. mpi_sendflags(cntr)) then
                ! check if the message has been sent
                quit_loop = .false.
                call MPI_TEST(mpi_req(cntr), mpi_sendflags(cntr), istat, ierr)
              end if

              mpi_recvfrom = this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr%rnk
              mpi_recvtag = (-ind3 + 2) + 3 * (-ind2 + 1) + 9 * (-ind1 + 1)

              if (.not. mpi_recvflags(cntr)) then
                quit_loop = .false.
                call MPI_IPROBE(mpi_recvfrom, mpi_recvtag, MPI_COMM_WORLD, mpi_recvflags(cntr), istat, ierr)
                if (mpi_recvflags(cntr)) then
                  call MPI_GET_COUNT(istat, myMPI_ENROUTE, cnt_recv_enroute, ierr)

                  if (cnt_recv_enroute .ge. size(recv_enroute%enroute)) then
                    call throwError('ERROR: particle had rcv array too small.')
                  end if

                  call MPI_RECV(recv_enroute%enroute(1:cnt_recv_enroute), cnt_recv_enroute, myMPI_ENROUTE,&
                  & mpi_recvfrom, mpi_recvtag, MPI_COMM_WORLD, istat, ierr)

                  ! write received data to local memory
                  if (cnt_recv_enroute .gt. 0) then
                    call extractParticlesFromEnroute(cnt_recv_enroute, s)
                  end if ! if > 0 particles received
                end if ! if the message can be received -> get the size & receive it
              end if ! if the message has already been received

            end do ! ind3
          end do ! ind2
        end do ! ind1
      end do ! global loop

#ifdef LOWMEM
      do send_x = -1, 1
        do send_y = -1, 1
          do send_z = -1, 1
            if (allocated(enroute_bot%get(send_x, send_y, send_z)%enroute)) then
              deallocate(enroute_bot%get(send_x, send_y, send_z)%enroute)
            end if
          end do
        end do
      end do
#endif

    end do ! loop over species
    call printDiag("exchangeParticles()", 2)
  end subroutine exchangeParticles

#endif

  subroutine redistributeParticlesBetweenMeshblocks()
    implicit none
    integer(kind=2), pointer, contiguous  :: pt_xi(:), pt_yi(:), pt_zi(:)
    integer, pointer, contiguous          :: pt_proc(:)
    real, pointer, contiguous             :: pt_dx(:), pt_dy(:), pt_dz(:)
    integer             :: s, p, send_x, send_y, send_z, ti, tj, tk, ti_p, tj_p, tk_p
    integer             :: mpi_sendto, mpi_recvfrom, mpi_sendtag, mpi_recvtag
    integer             :: ierr, ind1, ind2, ind3, cntr
    integer             :: cnt_recv_enroute

#ifdef MPI08
    type(MPI_REQUEST), allocatable  :: mpi_req(:)
    type(MPI_STATUS)                :: istat
#endif

#ifdef MPI
    integer, allocatable            :: mpi_req(:)
    integer                         :: istat(MPI_STATUS_SIZE)
#endif

    logical, allocatable              :: mpi_sendflags(:), mpi_recvflags(:)
    logical                           :: quit_loop

    allocate(mpi_req(sendrecv_neighbors))
    allocate(mpi_sendflags(sendrecv_neighbors))
    allocate(mpi_recvflags(sendrecv_neighbors))

    do s = 1, nspec ! loop over species

#ifdef LOWMEM
      enroute_bot%get(:,:,:)%cnt = 0
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
            pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
            pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi
            pt_dx => species(s)%prtl_tile(ti, tj, tk)%dx
            pt_dy => species(s)%prtl_tile(ti, tj, tk)%dy
            pt_dz => species(s)%prtl_tile(ti, tj, tk)%dz
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              ! send_* = -1 / 0 / +1
              send_x = 0; send_y = 0; send_z = 0
#if defined(oneD) || defined(twoD) || defined(threeD)
              send_x = (ISIGN(1, pt_xi(p) - this_meshblock%ptr%sx) + 1) / 2 - (ISIGN(1, -pt_xi(p) - 1) + 1) / 2
#endif
#if defined(twoD) || defined(threeD)
              send_y = (ISIGN(1, pt_yi(p) - this_meshblock%ptr%sy) + 1) / 2 - (ISIGN(1, -pt_yi(p) - 1) + 1) / 2
#endif
#if defined(threeD)
              send_z = (ISIGN(1, pt_zi(p) - this_meshblock%ptr%sz) + 1) / 2 - (ISIGN(1, -pt_zi(p) - 1) + 1) / 2
#endif
              if ((send_x .ne. 0) .or. (send_y .ne. 0) .or. (send_z .ne. 0)) then
                if (.not. associated(this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr)) then
                  cycle
                end if
                ! count the # of particles to be sent in each direction
                enroute_bot%get(send_x, send_y, send_z)%cnt = enroute_bot%get(send_x, send_y, send_z)%cnt + 1
              end if
            end do ! particles
            pt_xi => null(); pt_yi => null(); pt_zi => null()
            pt_dx => null(); pt_dy => null(); pt_dz => null()
          end do
        end do
      end do
      ! allocate the corresponding buffer arrays for each direction
      do send_x = -1, 1
        do send_y = -1, 1
          do send_z = -1, 1
            call reallocateEnroute(send_x, send_y, send_z,&
            & enroute_bot%get(send_x, send_y, send_z)%cnt+1)
          end do
        end do
      end do
#endif ! LOWMEM

      enroute_bot%get(:,:,:)%cnt = 0
      ! particle crosses MPI blocks //
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
            pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
            pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi

            pt_dx => species(s)%prtl_tile(ti, tj, tk)%dx
            pt_dy => species(s)%prtl_tile(ti, tj, tk)%dy
            pt_dz => species(s)%prtl_tile(ti, tj, tk)%dz

            pt_proc => species(s)%prtl_tile(ti, tj, tk)%proc
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              ! send_* = -1 / 0 / +1
              send_x = 0; send_y = 0; send_z = 0
#if defined(oneD) || defined(twoD) || defined(threeD)
              if (pt_xi(p) .lt. 0) then
                send_x = -1
              else if (pt_xi(p) .ge. this_meshblock%ptr%sx) then
                send_x = 1
              end if
#endif
#if defined(twoD) || defined(threeD)
              if (pt_yi(p) .lt. 0) then
                send_y = -1
              else if (pt_yi(p) .ge. this_meshblock%ptr%sy) then
                send_y = 1
              end if
#endif
#if defined(threeD)
              if (pt_zi(p) .lt. 0) then
                send_z = -1
              else if (pt_zi(p) .ge. this_meshblock%ptr%sz) then
                send_z = 1
              end if
#endif
              ! FIX1 check for null() boundaries
              if ((send_x .ne. 0) .or. (send_y .ne. 0) .or. (send_z .ne. 0)) then
                if (.not. associated(this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr)) then
                  ! make ghost particle
                  pt_proc(p) = -1
                  cycle
                end if

#ifdef DEBUG
                ! check that we're not sending across two meshblocks
                if (send_x .eq. -1) then
                  if (-pt_xi(p) .gt. this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sx) then
                    call throwError('ERROR: particle is trying to be sent accross more than one meshblock in '//trim(STR(send_x))//'.')
                  end if
                else if (send_x .eq. 1) then
                  if (pt_xi(p) - this_meshblock%ptr%sx .ge. this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sx) then
                    call throwError('ERROR: particle is trying to be sent accross more than one meshblock in '//trim(STR(send_x))//'.')
                  end if
                end if
                if (send_y .eq. -1) then
                  if (-pt_yi(p) .gt. this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sy) then
                    call throwError('ERROR: particle is trying to be sent accross more than one meshblock in '//trim(STR(send_y))//'.')
                  end if
                else if (send_y .eq. 1) then
                  if (pt_yi(p) - this_meshblock%ptr%sy .ge. this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sy) then
                    call throwError('ERROR: particle is trying to be sent accross more than one meshblock in '//trim(STR(send_y))//'.')
                  end if
                end if
                if (send_z .eq. -1) then
                  if (-pt_zi(p) .gt. this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sz) then
                    call throwError('ERROR: particle is trying to be sent accross more than one meshblock in '//trim(STR(send_z))//'.')
                  end if
                else if (send_z .eq. 1) then
                  if (pt_zi(p) - this_meshblock%ptr%sz .ge. this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sz) then
                    call throwError('ERROR: particle is trying to be sent accross more than one meshblock in '//trim(STR(send_z))//'.')
                  end if
                end if
#endif

                ! copy this particle to temporary `enroute_bot` array
                enroute_bot%get(send_x, send_y, send_z)%cnt = enroute_bot%get(send_x, send_y, send_z)%cnt + 1
                if (enroute_bot%get(send_x, send_y, send_z)%cnt .ge.&
                & enroute_bot%get(send_x, send_y, send_z)%max) then
                  call throwError('ERROR: particle send buffer array too small in `redistributeParticlesBetweenMeshblocks()`: '//&
                  & trim(STR(enroute_bot%get(send_x, send_y, send_z)%max))//'.')
                end if
                cntr = enroute_bot%get(send_x, send_y, send_z)%cnt
                call copyToEnroute(s, ti, tj, tk, p, enroute_bot%get(send_x, send_y, send_z)%enroute(cntr))
                ! make ghost particle
                pt_proc(p) = -1

                ! shift coordinates to fit the new grid
#if defined(oneD) || defined(twoD) || defined(threeD)
                if (send_x .eq. -1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi =&
                  & INT(enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi +&
                  & this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sx, 2)
                else if (send_x .eq. 1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi =&
                  & INT(enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi - this_meshblock%ptr%sx, 2)
                end if
#endif

#if defined(twoD) || defined(threeD)
                if (send_y .eq. -1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi +&
                  & this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sy
                else if (send_y .eq. 1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi - this_meshblock%ptr%sy
                end if
#endif

#if defined(threeD)
                if (send_z .eq. -1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi +&
                  & this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sz
                else if (send_z .eq. 1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi - this_meshblock%ptr%sz
                end if
#endif

#ifdef GCA
#if defined(oneD) || defined(twoD) || defined(threeD)
                if (send_x .eq. 1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi_past =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi_past -&
                  & this_meshblock%ptr%sx
                else if (send_x .eq. -1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi_past =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%xi_past +&
                  & this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sx
                end if
#endif

#if defined(twoD) || defined(threeD)
                if (send_y .eq. 1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi_past =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi_past -&
                  & this_meshblock%ptr%sy
                else if (send_y .eq. -1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi_past =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%yi_past +&
                  & this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sy
                end if
#endif

#if defined(threeD)
                if (send_z .eq. 1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi_past =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi_past -&
                  & this_meshblock%ptr%sz
                else if (send_z .eq. -1) then
                  enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi_past =&
                  & enroute_bot%get(send_x, send_y, send_z)%enroute(cntr)%zi_past +&
                  & this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sz
                end if
#endif
#endif
              end if
            end do ! particles
            pt_xi => null(); pt_yi => null(); pt_zi => null()
            pt_dx => null(); pt_dy => null(); pt_dz => null()
            pt_proc => null()
          end do ! tk
        end do ! tj
      end do ! ti
      ! // particle crosses MPI blocks

      ! start sending //
      cntr = 0
      do ind1 = -1, 1
        do ind2 = -1, 1
          do ind3 = -1, 1
            if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
#ifdef oneD
            if ((ind2 .ne. 0) .or. (ind3 .ne. 0)) cycle
#elif twoD
            if (ind3 .ne. 0) cycle
#endif
            if (.not. associated(this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr)) cycle
            cntr = cntr + 1
            mpi_sendtag = (ind3 + 2) + 3 * (ind2 + 1) + 9 * (ind1 + 1) + 100 * s

            ! post non-blocking send requests
            mpi_sendto = this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr%rnk
            call MPI_ISEND(enroute_bot%get(ind1,ind2,ind3)%enroute(1:enroute_bot%get(ind1,ind2,ind3)%cnt),&
            & enroute_bot%get(ind1,ind2,ind3)%cnt, myMPI_ENROUTE,&
            & mpi_sendto, mpi_sendtag, MPI_COMM_WORLD, mpi_req(cntr), ierr)
          end do
        end do
      end do
      ! // start sending

      ! particle moves between tiles within a single MPI block //
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
            pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
            pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi
            pt_proc => species(s)%prtl_tile(ti, tj, tk)%proc
            ! FIX1 make sure this is vectorized
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              if (pt_proc(p) .eq. -1) cycle
              ti_p = FLOOR(REAL(pt_xi(p)) / REAL(species(s)%tile_sx)) + 1
              tj_p = FLOOR(REAL(pt_yi(p)) / REAL(species(s)%tile_sy)) + 1
              tk_p = FLOOR(REAL(pt_zi(p)) / REAL(species(s)%tile_sz)) + 1
              if ((ti_p .ne. ti) .or. (tj_p .ne. tj) .or. (tk_p .ne. tk)) then
                call moveParticleBetweenTiles(s, ti, tj, tk, p)
              end if
            end do ! particles
            pt_xi => null(); pt_yi => null(); pt_zi => null()
            pt_proc => null()
          end do ! tk
        end do ! tj
      end do ! ti
      ! // particle moves between tiles within a single MPI block

      ! wait to send & receive all the MPI calls and write data to memory
      quit_loop = .false.
      mpi_sendflags(:) = .false.
      mpi_recvflags(:) = .false.
      do while (.not. quit_loop)
        ! try to receive while not done
        quit_loop = .true.
        cntr = 0
        do ind1 = -1, 1
          do ind2 = -1, 1
            do ind3 = -1, 1
              if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
#ifdef oneD
              if ((ind2 .ne. 0) .or. (ind3 .ne. 0)) cycle
#elif twoD
              if (ind3 .ne. 0) cycle
#endif
              if (.not. associated(this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr)) cycle
              cntr = cntr + 1

              if (.not. mpi_sendflags(cntr)) then
                ! check if the message has been sent
                quit_loop = .false.
                call MPI_TEST(mpi_req(cntr), mpi_sendflags(cntr), istat, ierr)
              end if

              mpi_recvfrom = this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr%rnk
              mpi_recvtag = (-ind3 + 2) + 3 * (-ind2 + 1) + 9 * (-ind1 + 1) + 100 * s

              if (.not. mpi_recvflags(cntr)) then
                quit_loop = .false.
                call MPI_IPROBE(mpi_recvfrom, mpi_recvtag, MPI_COMM_WORLD, mpi_recvflags(cntr), istat, ierr)
                if (mpi_recvflags(cntr)) then
                  call MPI_GET_COUNT(istat, myMPI_ENROUTE, recv_enroute%cnt, ierr)

                  if (cnt_recv_enroute .ge. size(recv_enroute%enroute)) then
                    call throwError('ERROR: particle had rcv array too small.')
                  end if

                  call MPI_RECV(recv_enroute%enroute(1:recv_enroute%cnt), recv_enroute%cnt, myMPI_ENROUTE,&
                  & mpi_recvfrom, mpi_recvtag, MPI_COMM_WORLD, istat, ierr)
                  ! write received data to local memory
                  if (recv_enroute%cnt .gt. 0) then
                    call extractParticlesFromEnroute(recv_enroute%cnt, s)
                  end if ! if > 0 particles received
                end if ! if the message can be received -> get the size & receive it
              end if ! if the message has already been received

            end do ! ind3
          end do ! ind2
        end do ! ind1
      end do ! global loop

#ifdef DEBUG
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
            pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
            pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi
            pt_proc => species(s)%prtl_tile(ti, tj, tk)%proc
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              if (pt_proc(p) .eq. -1) cycle
              if ((pt_xi(p) .lt. species(s)%prtl_tile(ti, tj, tk)%x1) .or. &
              & (pt_xi(p) .ge. species(s)%prtl_tile(ti, tj, tk)%x2) .or. &
              & (pt_yi(p) .lt. species(s)%prtl_tile(ti, tj, tk)%y1) .or. &
              & (pt_yi(p) .ge. species(s)%prtl_tile(ti, tj, tk)%y2) .or. &
              & (pt_zi(p) .lt. species(s)%prtl_tile(ti, tj, tk)%z1) .or. &
              & (pt_zi(p) .ge. species(s)%prtl_tile(ti, tj, tk)%z2)) then
                print *, 'xyz', pt_xi(p), pt_yi(p), pt_zi(p)
                print *, 'tx:', species(s)%prtl_tile(ti, tj, tk)%x1, species(s)%prtl_tile(ti, tj, tk)%x2
                print *, 'ty:', species(s)%prtl_tile(ti, tj, tk)%y1, species(s)%prtl_tile(ti, tj, tk)%y2
                print *, 'tz:', species(s)%prtl_tile(ti, tj, tk)%z1, species(s)%prtl_tile(ti, tj, tk)%z2
                print *, 'sxyz:', this_meshblock%ptr%sx, this_meshblock%ptr%sy, this_meshblock%ptr%sz
                call throwError('ERROR: particle in wrong tile after exchange')
              end if
            end do
          end do
        end do
      end do
#endif

#ifdef LOWMEM
      do send_x = -1, 1
        do send_y = -1, 1
          do send_z = -1, 1
            if (allocated(enroute_bot%get(send_x, send_y, send_z)%enroute)) then
              deallocate(enroute_bot%get(send_x, send_y, send_z)%enroute)
            end if
          end do
        end do
      end do
#endif

    end do ! loop over species
    call printDiag("redistributeParticlesBetweenMeshblocks()", 3)
  end subroutine redistributeParticlesBetweenMeshblocks

  subroutine moveParticleBetweenTiles(s, ti, tj, tk, p)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in) :: s, ti, tj, tk, p
    call createParticleFromAttributes(s=s, xi=species(s)%prtl_tile(ti, tj, tk)%xi(p),&
    & yi=species(s)%prtl_tile(ti, tj, tk)%yi(p),&
    & zi=species(s)%prtl_tile(ti, tj, tk)%zi(p),&
    & dx=species(s)%prtl_tile(ti, tj, tk)%dx(p),&
    & dy=species(s)%prtl_tile(ti, tj, tk)%dy(p),&
    & dz=species(s)%prtl_tile(ti, tj, tk)%dz(p),&
#ifdef GCA
    & xi_past=species(s)%prtl_tile(ti, tj, tk)%xi_past(p),&
    & yi_past=species(s)%prtl_tile(ti, tj, tk)%yi_past(p),&
    & zi_past=species(s)%prtl_tile(ti, tj, tk)%zi_past(p),&
    & dx_past=species(s)%prtl_tile(ti, tj, tk)%dx_past(p),&
    & dy_past=species(s)%prtl_tile(ti, tj, tk)%dy_past(p),&
    & dz_past=species(s)%prtl_tile(ti, tj, tk)%dz_past(p),&
#endif
    & u=species(s)%prtl_tile(ti, tj, tk)%u(p),&
    & v=species(s)%prtl_tile(ti, tj, tk)%v(p),&
    & w=species(s)%prtl_tile(ti, tj, tk)%w(p),&
#ifdef GCA
    & u_eff=species(s)%prtl_tile(ti, tj, tk)%u_eff(p),&
    & v_eff=species(s)%prtl_tile(ti, tj, tk)%v_eff(p),&
    & w_eff=species(s)%prtl_tile(ti, tj, tk)%w_eff(p),&
    & u_par=species(s)%prtl_tile(ti, tj, tk)%u_par(p),&
    & u_perp=species(s)%prtl_tile(ti, tj, tk)%u_perp(p),&
#endif
#ifdef PRTLPAYLOADS
    & payload1=species(s)%prtl_tile(ti, tj, tk)%payload1(p),&
    & payload2=species(s)%prtl_tile(ti, tj, tk)%payload2(p),&
    & payload3=species(s)%prtl_tile(ti, tj, tk)%payload3(p),&
#endif
    & ind=species(s)%prtl_tile(ti, tj, tk)%ind(p),&
    & proc=species(s)%prtl_tile(ti, tj, tk)%proc(p),&
    & weight=species(s)%prtl_tile(ti, tj, tk)%weight(p))
    ! schedule particle for deletion
    species(s)%prtl_tile(ti, tj, tk)%proc(p) = -1
  end subroutine moveParticleBetweenTiles

end module m_exchangeparts
