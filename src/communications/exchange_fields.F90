module m_exchangefields
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_fields
  
#ifndef MPINONBLOCK
  private :: findCnt, bufferSendArray
#endif

contains

#ifndef MPINONBLOCK
  ! blocking MPI communication
  subroutine findCnt(ind1, ind2, ind3, exchangeE, exchangeB, send_cnt)
    implicit none
    integer :: imin, imax, jmin, jmax, kmin, kmax
    logical, intent(in) :: exchangeE, exchangeB
    integer, intent(in) :: ind1, ind2, ind3
    integer, intent(out) :: send_cnt
    ! highlight the region to send and save to `send_fld`
    imin = -2 * NGHOST
    imax = -2 * NGHOST
    jmin = -2 * NGHOST
    jmax = -2 * NGHOST
    kmin = -2 * NGHOST
    kmax = -2 * NGHOST

    if (ind1 .eq. 0) then
      imin = 0; imax = this_meshblock % ptr % sx - 1
    else if (ind1 .eq. -1) then
      imin = 0; imax = NGHOST - 1
    else if (ind1 .eq. 1) then
      imin = this_meshblock % ptr % sx - NGHOST; imax = this_meshblock % ptr % sx - 1
    end if
    if (ind2 .eq. 0) then
      jmin = 0; jmax = this_meshblock % ptr % sy - 1
    else if (ind2 .eq. -1) then
      jmin = 0; jmax = NGHOST - 1
    else if (ind2 .eq. 1) then
      jmin = this_meshblock % ptr % sy - NGHOST; jmax = this_meshblock % ptr % sy - 1
    end if
    if (ind3 .eq. 0) then
      kmin = 0; kmax = this_meshblock % ptr % sz - 1
    else if (ind3 .eq. -1) then
      kmin = 0; kmax = NGHOST - 1
    else if (ind3 .eq. 1) then
      kmin = this_meshblock % ptr % sz - NGHOST; kmax = this_meshblock % ptr % sz - 1
    end if
#ifdef oneD
    jmin = 0; jmax = 0
    kmin = 0; kmax = 0
#elif twoD
    kmin = 0; kmax = 0
#endif

    if ((imin .eq. -2 * NGHOST) .or. (imax .eq. -2 * NGHOST) .or. &
        (jmin .eq. -2 * NGHOST) .or. (jmax .eq. -2 * NGHOST) .or. &
        (kmin .eq. -2 * NGHOST) .or. (kmax .eq. -2 * NGHOST)) then
      call throwError("Error: invalid index evaluation in findCnt")
    end if

    send_cnt = 1
    if (exchangeE) then
      send_cnt = send_cnt + 3 * (imax - imin + 1) * (jmax - jmin + 1) * (kmax - kmin + 1)
    end if

    if (exchangeB) then
      send_cnt = send_cnt + 3 * (imax - imin + 1) * (jmax - jmin + 1) * (kmax - kmin + 1)
    end if
    send_cnt = send_cnt - 1
  end subroutine findCnt

  subroutine bufferSendArray(offset, ind1, ind2, ind3, exchangeE, exchangeB, send_cnt)
    implicit none
    integer :: imin, imax, jmin, jmax, kmin, kmax, i, j, k
    logical, intent(in) :: exchangeE, exchangeB
    integer, intent(in) :: ind1, ind2, ind3
    integer, intent(in) :: offset
    integer, intent(out) :: send_cnt
    imin = -2 * NGHOST
    imax = -2 * NGHOST
    jmin = -2 * NGHOST
    jmax = -2 * NGHOST
    kmin = -2 * NGHOST
    kmax = -2 * NGHOST

    ! highlight the region to send and save to `send_fld`
    if (ind1 .eq. 0) then
      imin = 0; imax = this_meshblock % ptr % sx - 1
    else if (ind1 .eq. -1) then
      imin = 0; imax = NGHOST - 1
    else if (ind1 .eq. 1) then
      imin = this_meshblock % ptr % sx - NGHOST; imax = this_meshblock % ptr % sx - 1
    end if
    if (ind2 .eq. 0) then
      jmin = 0; jmax = this_meshblock % ptr % sy - 1
    else if (ind2 .eq. -1) then
      jmin = 0; jmax = NGHOST - 1
    else if (ind2 .eq. 1) then
      jmin = this_meshblock % ptr % sy - NGHOST; jmax = this_meshblock % ptr % sy - 1
    end if
    if (ind3 .eq. 0) then
      kmin = 0; kmax = this_meshblock % ptr % sz - 1
    else if (ind3 .eq. -1) then
      kmin = 0; kmax = NGHOST - 1
    else if (ind3 .eq. 1) then
      kmin = this_meshblock % ptr % sz - NGHOST; kmax = this_meshblock % ptr % sz - 1
    end if
#ifdef oneD
    jmin = 0; jmax = 0
    kmin = 0; kmax = 0
#elif twoD
    kmin = 0; kmax = 0
#endif

    if ((imin .eq. -2 * NGHOST) .or. (imax .eq. -2 * NGHOST) .or. &
        (jmin .eq. -2 * NGHOST) .or. (jmax .eq. -2 * NGHOST) .or. &
        (kmin .eq. -2 * NGHOST) .or. (kmax .eq. -2 * NGHOST)) then
      call throwError("Error: invalid index evaluation in bufferSendArray")
    end if

    send_cnt = 1
    do i = imin, imax
      do j = jmin, jmax
        do k = kmin, kmax
          if (exchangeE) then
            send_EB(offset + send_cnt + 0) = ex(i, j, k)
            send_EB(offset + send_cnt + 1) = ey(i, j, k)
            send_EB(offset + send_cnt + 2) = ez(i, j, k)
            send_cnt = send_cnt + 3
          end if
          if (exchangeB) then
            send_EB(offset + send_cnt + 0) = bx(i, j, k)
            send_EB(offset + send_cnt + 1) = by(i, j, k)
            send_EB(offset + send_cnt + 2) = bz(i, j, k)
            send_cnt = send_cnt + 3
          end if
        end do
      end do
    end do
    send_cnt = send_cnt - 1
  end subroutine bufferSendArray

  subroutine extractRecvArray(ind1, ind2, ind3, exchangeE, exchangeB)
    implicit none
    integer :: imin, imax, jmin, jmax, kmin, kmax, i, j, k
    integer, intent(in) :: ind1, ind2, ind3
    logical, intent(in) :: exchangeE, exchangeB
    integer :: cnt
    imin = -2 * NGHOST
    imax = -2 * NGHOST
    jmin = -2 * NGHOST
    jmax = -2 * NGHOST
    kmin = -2 * NGHOST
    kmax = -2 * NGHOST

    ! highlight the region to extract the `recv_fld`
    if (ind1 .eq. 0) then
      imin = 0; imax = this_meshblock % ptr % sx - 1
    else if (ind1 .eq. -1) then
      imin = -NGHOST; imax = -1
    else if (ind1 .eq. 1) then
      imin = this_meshblock % ptr % sx; imax = this_meshblock % ptr % sx + NGHOST - 1
    end if
    if (ind2 .eq. 0) then
      jmin = 0; jmax = this_meshblock % ptr % sy - 1
    else if (ind2 .eq. -1) then
      jmin = -NGHOST; jmax = -1
    else if (ind2 .eq. 1) then
      jmin = this_meshblock % ptr % sy; jmax = this_meshblock % ptr % sy + NGHOST - 1
    end if
    if (ind3 .eq. 0) then
      kmin = 0; kmax = this_meshblock % ptr % sz - 1
    else if (ind3 .eq. -1) then
      kmin = -NGHOST; kmax = -1
    else if (ind3 .eq. 1) then
      kmin = this_meshblock % ptr % sz; kmax = this_meshblock % ptr % sz + NGHOST - 1
    end if
#ifdef oneD
    jmin = 0; jmax = 0
    kmin = 0; kmax = 0
#elif twoD
    kmin = 0; kmax = 0
#endif

    if ((imin .eq. -2 * NGHOST) .or. (imax .eq. -2 * NGHOST) .or. &
        (jmin .eq. -2 * NGHOST) .or. (jmax .eq. -2 * NGHOST) .or. &
        (kmin .eq. -2 * NGHOST) .or. (kmax .eq. -2 * NGHOST)) then
      call throwError("Error: invalid index evaluation in findCnt")
    end if

    ! copy `recv_fld` to ghost cells
    cnt = 1
    do i = imin, imax
      do j = jmin, jmax
        do k = kmin, kmax
          if (exchangeE) then
            ex(i, j, k) = recv_fld(cnt + 0)
            ey(i, j, k) = recv_fld(cnt + 1)
            ez(i, j, k) = recv_fld(cnt + 2)
            cnt = cnt + 3
          end if
          if (exchangeB) then
            bx(i, j, k) = recv_fld(cnt + 0)
            by(i, j, k) = recv_fld(cnt + 1)
            bz(i, j, k) = recv_fld(cnt + 2)
            cnt = cnt + 3
          end if
        end do
      end do
    end do
  end subroutine extractRecvArray

  subroutine exchangeFields(exchangeE, exchangeB)
    implicit none
    logical, intent(in) :: exchangeE, exchangeB
    integer :: ind1, ind2, ind3
    integer :: cnt, ierr
    integer :: mpi_sendto, mpi_recvfrom, mpi_tag
    logical :: should_send, should_recv

#ifdef MPI08
    type(MPI_STATUS) :: istat
#endif

#ifdef MPI
    integer :: istat(MPI_STATUS_SIZE)
#endif

    if ((.not. exchangeE) .and. (.not. exchangeB)) then
      call throwError('ERROR: `exchangeFields()` called with `.false.` and `.false.`')
    end if

    ! looping through all send directions
    do ind1 = -1, 1
      do ind2 = -1, 1
        do ind3 = -1, 1
          if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
#ifdef oneD
          if ((ind2 .ne. 0) .or. (ind3 .ne. 0)) cycle
#elif twoD
          if (ind3 .ne. 0) cycle
#endif
          mpi_tag = 100 + (ind3 + 2) + 3 * (ind2 + 1) + 9 * (ind1 + 1)
          should_send = associated(this_meshblock % ptr % neighbor(ind1, ind2, ind3) % ptr)
          should_recv = associated(this_meshblock % ptr % neighbor(-ind1, -ind2, -ind3) % ptr)
          if (should_send .and. should_recv) then
            mpi_sendto = this_meshblock % ptr % neighbor(ind1, ind2, ind3) % ptr % rnk
            mpi_recvfrom = this_meshblock % ptr % neighbor(-ind1, -ind2, -ind3) % ptr % rnk

            call bufferSendArray(0, ind1, ind2, ind3, exchangeE, exchangeB, cnt)

            call MPI_SENDRECV(send_EB(1:cnt), cnt, default_mpi_real, mpi_sendto, mpi_tag, &
                              recv_fld(1:cnt), cnt, default_mpi_real, mpi_recvfrom, mpi_tag, &
                              MPI_COMM_WORLD, istat, ierr)
            call extractRecvArray(-ind1, -ind2, -ind3, exchangeE, exchangeB)
          else if ((.not. should_send) .and. should_recv) then
            mpi_recvfrom = this_meshblock % ptr % neighbor(-ind1, -ind2, -ind3) % ptr % rnk
            call findCnt(ind1, ind2, ind3, exchangeE, exchangeB, cnt)
            call MPI_RECV(recv_fld(1:cnt), cnt, default_mpi_real, mpi_recvfrom, mpi_tag, MPI_COMM_WORLD, istat, ierr)
            call extractRecvArray(-ind1, -ind2, -ind3, exchangeE, exchangeB)
          else if ((.not. should_recv) .and. should_send) then
            mpi_sendto = this_meshblock % ptr % neighbor(ind1, ind2, ind3) % ptr % rnk
            call bufferSendArray(0, ind1, ind2, ind3, exchangeE, exchangeB, cnt)
            call MPI_SEND(send_EB(1:cnt), cnt, default_mpi_real, mpi_sendto, mpi_tag, MPI_COMM_WORLD, ierr)
          end if
        end do
      end do
    end do
    call printDiag("exchangeFields()", 2)
  end subroutine exchangeFields

#else
  ! non-blocking MPI communication
  subroutine exchangeFields(exchangeE, exchangeB)
    implicit none
    logical, intent(in) :: exchangeE, exchangeB
    integer :: i, j, k, imin, imax, jmin, jmax, kmin, kmax
    integer :: ind1, ind2, ind3, cntr, n_cntr
    integer :: send_cnt, recv_cnt, ierr
    integer :: mpi_sendto, mpi_recvfrom, mpi_sendtag, mpi_recvtag
    integer :: mpi_offset

#ifdef MPI08
    type(MPI_REQUEST), allocatable :: mpi_req(:)
    type(MPI_STATUS) :: istat
#endif

#ifdef MPI
    integer, allocatable :: mpi_req(:)
    integer :: istat(MPI_STATUS_SIZE)
#endif

    logical, allocatable :: mpi_sendflags(:), mpi_recvflags(:)
    logical :: quit_loop

    if ((.not. exchangeE) .and. (.not. exchangeB)) then
      call throwError('ERROR: `exchangeFields()` called with `.false.` and `.false.`')
    end if

    allocate (mpi_req(sendrecv_neighbors))
    allocate (mpi_sendflags(sendrecv_neighbors))
    allocate (mpi_recvflags(sendrecv_neighbors))

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
          if (.not. associated(this_meshblock % ptr % neighbor(ind1, ind2, ind3) % ptr)) cycle
          cntr = cntr + 1

          mpi_sendto = this_meshblock % ptr % neighbor(ind1, ind2, ind3) % ptr % rnk
          mpi_sendtag = (ind3 + 2) + 3 * (ind2 + 1) + 9 * (ind1 + 1)

          ! highlight the region to send and save to `send_fld`
          if (ind1 .eq. 0) then
            imin = 0; imax = this_meshblock % ptr % sx - 1
          else if (ind1 .eq. -1) then
            imin = 0; imax = NGHOST - 1
          else if (ind1 .eq. 1) then
            imin = this_meshblock % ptr % sx - NGHOST; imax = this_meshblock % ptr % sx - 1
          end if
          if (ind2 .eq. 0) then
            jmin = 0; jmax = this_meshblock % ptr % sy - 1
          else if (ind2 .eq. -1) then
            jmin = 0; jmax = NGHOST - 1
          else if (ind2 .eq. 1) then
            jmin = this_meshblock % ptr % sy - NGHOST; jmax = this_meshblock % ptr % sy - 1
          end if
          if (ind3 .eq. 0) then
            kmin = 0; kmax = this_meshblock % ptr % sz - 1
          else if (ind3 .eq. -1) then
            kmin = 0; kmax = NGHOST - 1
          else if (ind3 .eq. 1) then
            kmin = this_meshblock % ptr % sz - NGHOST; kmax = this_meshblock % ptr % sz - 1
          end if
#ifdef oneD
          jmin = 0; jmax = 0
          kmin = 0; kmax = 0
#elif twoD
          kmin = 0; kmax = 0
#endif

          ! write send/recv arrays in/from a given direction
          !     in 3D: 26 directions, in 2D: 8, in 1D: 2
          mpi_offset = (cntr - 1) * sendrecv_offsetsz
          send_cnt = 1
          do i = imin, imax
            do j = jmin, jmax
              do k = kmin, kmax
                if (exchangeE) then
                  send_fld(mpi_offset + send_cnt + 0) = ex(i, j, k)
                  send_fld(mpi_offset + send_cnt + 1) = ey(i, j, k)
                  send_fld(mpi_offset + send_cnt + 2) = ez(i, j, k)
                  send_cnt = send_cnt + 3
                end if
                if (exchangeB) then
                  send_fld(mpi_offset + send_cnt + 0) = bx(i, j, k)
                  send_fld(mpi_offset + send_cnt + 1) = by(i, j, k)
                  send_fld(mpi_offset + send_cnt + 2) = bz(i, j, k)
                  send_cnt = send_cnt + 3
                end if
              end do
            end do
          end do
          send_cnt = send_cnt - 1

          ! post non-blocking send requests
          call MPI_ISEND(send_fld(mpi_offset + 1:mpi_offset + send_cnt), send_cnt, default_mpi_real, &
                         mpi_sendto, mpi_sendtag, MPI_COMM_WORLD, mpi_req(cntr), ierr)
        end do
      end do
    end do

    ! wait to send & receive all the MPI calls and write data to memory
    quit_loop = .false.
    mpi_sendflags(:) = .false.
    mpi_recvflags(:) = .false.
    do while (.not. quit_loop)
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
            if (.not. associated(this_meshblock % ptr % neighbor(ind1, ind2, ind3) % ptr)) cycle
            cntr = cntr + 1

            if (.not. mpi_sendflags(cntr)) then
              ! check if the message has been sent
              quit_loop = .false.
              call MPI_TEST(mpi_req(cntr), mpi_sendflags(cntr), istat, ierr)
            end if

            mpi_recvfrom = this_meshblock % ptr % neighbor(ind1, ind2, ind3) % ptr % rnk
            mpi_recvtag = (-ind3 + 2) + 3 * (-ind2 + 1) + 9 * (-ind1 + 1)

            if (.not. mpi_recvflags(cntr)) then
              quit_loop = .false.
              call MPI_IPROBE(mpi_recvfrom, mpi_recvtag, MPI_COMM_WORLD, mpi_recvflags(cntr), istat, ierr)
              if (mpi_recvflags(cntr)) then
                ! if the message is ready to be received -> get the size & receive it
                call MPI_GET_COUNT(istat, default_mpi_real, recv_cnt, ierr)
                call MPI_RECV(recv_fld(1:recv_cnt), recv_cnt, default_mpi_real, &
                              mpi_recvfrom, mpi_recvtag, MPI_COMM_WORLD, istat, ierr)

                ! write received data to local memory
                ! highlight the region to extract the `recv_fld`
                if (ind1 .eq. 0) then
                  imin = 0; imax = this_meshblock % ptr % sx - 1
                else if (ind1 .eq. -1) then
                  imin = -NGHOST; imax = -1
                else if (ind1 .eq. 1) then
                  imin = this_meshblock % ptr % sx; imax = this_meshblock % ptr % sx + NGHOST - 1
                end if
                if (ind2 .eq. 0) then
                  jmin = 0; jmax = this_meshblock % ptr % sy - 1
                else if (ind2 .eq. -1) then
                  jmin = -NGHOST; jmax = -1
                else if (ind2 .eq. 1) then
                  jmin = this_meshblock % ptr % sy; jmax = this_meshblock % ptr % sy + NGHOST - 1
                end if
                if (ind3 .eq. 0) then
                  kmin = 0; kmax = this_meshblock % ptr % sz - 1
                else if (ind3 .eq. -1) then
                  kmin = -NGHOST; kmax = -1
                else if (ind3 .eq. 1) then
                  kmin = this_meshblock % ptr % sz; kmax = this_meshblock % ptr % sz + NGHOST - 1
                end if
#ifdef oneD
                jmin = 0; jmax = 0
                kmin = 0; kmax = 0
#elif twoD
                kmin = 0; kmax = 0
#endif

                ! copy `recv_fld` to ghost cells
                send_cnt = 1
                do i = imin, imax
                  do j = jmin, jmax
                    do k = kmin, kmax
                      if (exchangeE) then
                        ex(i, j, k) = recv_fld(send_cnt + 0)
                        ey(i, j, k) = recv_fld(send_cnt + 1)
                        ez(i, j, k) = recv_fld(send_cnt + 2)
                        send_cnt = send_cnt + 3
                      end if
                      if (exchangeB) then
                        bx(i, j, k) = recv_fld(send_cnt + 0)
                        by(i, j, k) = recv_fld(send_cnt + 1)
                        bz(i, j, k) = recv_fld(send_cnt + 2)
                        send_cnt = send_cnt + 3
                      end if
                    end do
                  end do
                end do

              end if ! if receive is ready
            end if ! if receive hasn't been done yet
          end do ! ind3
        end do ! ind2
      end do ! ind1
    end do ! global loop
    call printDiag("exchangeFields()", 2)
  end subroutine exchangeFields

#endif

  subroutine exchangeFieldSlabInX(rnk1, rnk2, slab)
    implicit none
    integer, intent(in) :: rnk1, rnk2, slab

    ! `slab < 0` means `rnk1` is sending
    ! `slab > 0` means `rnk1` is receiving
    ! `rnk1` is assumed closer to origin than `rnk2`

    if ((mpi_rank .eq. rnk1) .or. (mpi_rank .eq. rnk2)) then
#ifdef DEBUG
      ! check that sizes along `X` match
      if ((meshblocks(rnk1 + 1) % sy .ne. meshblocks(rnk2 + 1) % sy) .or. &
          (meshblocks(rnk1 + 1) % sz .ne. meshblocks(rnk2 + 1) % sz)) then
        call throwError('ERROR: sizes in `YZ` not matching: '//trim(STR(rnk1))//' '//trim(STR(rnk2))//'.')
      end if
      ! check that `rnk1` and `rnk2` are neighbors
      if (mpi_rank .eq. rnk1) then
        if (this_meshblock % ptr % neighbor(1, 0, 0) % ptr % rnk .ne. rnk2) then
          call throwError('ERROR: #'//trim(STR(rnk1))//' and #'//trim(STR(rnk2))//' are not neighbors.')
        end if
      else if (mpi_rank .eq. rnk2) then
        if (this_meshblock % ptr % neighbor(-1, 0, 0) % ptr % rnk .ne. rnk1) then
          call throwError('ERROR: #'//trim(STR(rnk2))//' and #'//trim(STR(rnk1))//' are not neighbors.')
        end if
      end if
#endif
      if (slab .gt. 0) then
        ! `rnk1` receiving, `rnk2` sending
        call SendRecvSlabInX(rnk2, rnk1, slab, +1)
      else if (slab .lt. 0) then
        ! `rnk1` sending, `rnk2` receiving
        call SendRecvSlabInX(rnk1, rnk2, abs(slab), -1)
      else
        call throwError('ERROR: `slab` cannot be zero in exchangeFieldSlabIn*().')
      end if
    end if
  end subroutine exchangeFieldSlabInX

  subroutine exchangeFieldSlabInY(rnk1, rnk2, slab)
    implicit none
    integer, intent(in) :: rnk1, rnk2, slab

    ! `slab < 0` means `rnk1` is sending
    ! `slab > 0` means `rnk1` is receiving
    ! `rnk1` is assumed closer to origin than `rnk2`

    if ((mpi_rank .eq. rnk1) .or. (mpi_rank .eq. rnk2)) then
#ifdef DEBUG
      ! check that sizes along `X` match
      if ((meshblocks(rnk1 + 1) % sx .ne. meshblocks(rnk2 + 1) % sx) .or. &
          (meshblocks(rnk1 + 1) % sz .ne. meshblocks(rnk2 + 1) % sz)) then
        call throwError('ERROR: sizes in `XZ` not matching: '//trim(STR(rnk1))//' '//trim(STR(rnk2))//'.')
      end if
      ! check that `rnk1` and `rnk2` are neighbors
      if (mpi_rank .eq. rnk1) then
        if (this_meshblock % ptr % neighbor(0, 1, 0) % ptr % rnk .ne. rnk2) then
          call throwError('ERROR: #'//trim(STR(rnk1))//' and #'//trim(STR(rnk2))//' are not neighbors.')
        end if
      else if (mpi_rank .eq. rnk2) then
        if (this_meshblock % ptr % neighbor(0, -1, 0) % ptr % rnk .ne. rnk1) then
          call throwError('ERROR: #'//trim(STR(rnk2))//' and #'//trim(STR(rnk1))//' are not neighbors.')
        end if
      end if
#endif
      if (slab .gt. 0) then
        ! `rnk1` receiving, `rnk2` sending
        call SendRecvSlabInY(rnk2, rnk1, slab, +1)
      else if (slab .lt. 0) then
        ! `rnk1` sending, `rnk2` receiving
        call SendRecvSlabInY(rnk1, rnk2, abs(slab), -1)
      else
        call throwError('ERROR: `slab` cannot be zero in exchangeFieldSlabIn*().')
      end if
    end if
  end subroutine exchangeFieldSlabInY

  subroutine exchangeFieldSlabInZ(rnk1, rnk2, slab)
    implicit none
    integer, intent(in) :: rnk1, rnk2, slab

    ! `slab < 0` means `rnk1` is sending
    ! `slab > 0` means `rnk1` is receiving
    ! `rnk1` is assumed closer to origin than `rnk2`

    if ((mpi_rank .eq. rnk1) .or. (mpi_rank .eq. rnk2)) then
#ifdef DEBUG
      ! check that sizes along `X` match
      if ((meshblocks(rnk1 + 1) % sx .ne. meshblocks(rnk2 + 1) % sx) .or. &
          (meshblocks(rnk1 + 1) % sy .ne. meshblocks(rnk2 + 1) % sy)) then
        call throwError('ERROR: sizes in `XY` not matching: '//trim(STR(rnk1))//' '//trim(STR(rnk2))//'.')
      end if
      ! check that `rnk1` and `rnk2` are neighbors
      if (mpi_rank .eq. rnk1) then
        if (this_meshblock % ptr % neighbor(0, 0, 1) % ptr % rnk .ne. rnk2) then
          call throwError('ERROR: #'//trim(STR(rnk1))//' and #'//trim(STR(rnk2))//' are not neighbors.')
        end if
      else if (mpi_rank .eq. rnk2) then
        if (this_meshblock % ptr % neighbor(0, 0, -1) % ptr % rnk .ne. rnk1) then
          call throwError('ERROR: #'//trim(STR(rnk2))//' and #'//trim(STR(rnk1))//' are not neighbors.')
        end if
      end if
#endif
      if (slab .gt. 0) then
        ! `rnk1` receiving, `rnk2` sending
        call SendRecvSlabInZ(rnk2, rnk1, slab, +1)
      else if (slab .lt. 0) then
        ! `rnk1` sending, `rnk2` receiving
        call SendRecvSlabInZ(rnk1, rnk2, abs(slab), -1)
      else
        call throwError('ERROR: `slab` cannot be zero in exchangeFieldSlabIn*().')
      end if
    end if
  end subroutine exchangeFieldSlabInZ

  subroutine SendRecvSlabInX(rnk_send, rnk_recv, slab, direction)
    implicit none
    integer, intent(in) :: rnk_send, rnk_recv, slab, direction
    real, allocatable :: buffer(:, :, :, :)
    integer :: i, j, k
    integer :: i1_send, i2_send, j1_send, j2_send, k1_send, k2_send
    integer :: i_, dummy_, size_, ierr, mpi_tag

#ifdef MPI08
    type(MPI_STATUS) :: istat
#endif

#ifdef MPI
    integer :: istat(MPI_STATUS_SIZE)
#endif

    if ((mpi_rank .eq. rnk_send) .or. (mpi_rank .eq. rnk_recv)) then
#ifdef DEBUG
      ! check that we're not sending too much
      if (slab .ge. meshblocks(rnk_send + 1) % sx - NGHOST - 1) then
        call throwError('ERROR: cannot send more than `sx - NGHOST - 1` cells in x.')
      end if
#endif

      i1_send = 0; i2_send = slab - 1
      j1_send = 0; j2_send = meshblocks(rnk_send + 1) % sy - 1
      k1_send = 0; k2_send = meshblocks(rnk_send + 1) % sz - 1

      allocate (buffer(i1_send:i2_send, j1_send:j2_send, k1_send:k2_send, 6))
      size_ = (i2_send - i1_send + 1) * (j2_send - j1_send + 1) * (k2_send - k1_send + 1) * 6

      mpi_tag = 1

      if (mpi_rank .eq. rnk_send) then
        ! fill the buffer
        if (direction .lt. 0) then
          dummy_ = this_meshblock % ptr % sx - slab
        else
          dummy_ = 0
        end if
        do k = k1_send, k2_send
          do j = j1_send, j2_send
            do i = i1_send, i2_send
              i_ = dummy_ + i
              buffer(i, j, k, 1) = ex_back(i_, j, k)
              buffer(i, j, k, 2) = ey_back(i_, j, k)
              buffer(i, j, k, 3) = ez_back(i_, j, k)
              buffer(i, j, k, 4) = bx_back(i_, j, k)
              buffer(i, j, k, 5) = by_back(i_, j, k)
              buffer(i, j, k, 6) = bz_back(i_, j, k)
            end do
          end do
        end do
        call MPI_SEND(buffer, size_, default_mpi_real, rnk_recv, mpi_tag, MPI_COMM_WORLD, ierr)
      else if (mpi_rank .eq. rnk_recv) then
        call MPI_RECV(buffer, size_, default_mpi_real, rnk_send, mpi_tag, MPI_COMM_WORLD, istat, ierr)
        ! extract from the buffer
        if (direction .lt. 0) then
          dummy_ = 0
        else
          dummy_ = this_meshblock % ptr % sx
        end if
        do k = k1_send, k2_send
          do j = j1_send, j2_send
            do i = i1_send, i2_send
              i_ = dummy_ + i
              ex(i_, j, k) = buffer(i, j, k, 1)
              ey(i_, j, k) = buffer(i, j, k, 2)
              ez(i_, j, k) = buffer(i, j, k, 3)
              bx(i_, j, k) = buffer(i, j, k, 4)
              by(i_, j, k) = buffer(i, j, k, 5)
              bz(i_, j, k) = buffer(i, j, k, 6)
            end do
          end do
        end do
      end if
      deallocate (buffer)
    end if
  end subroutine SendRecvSlabInX

  subroutine SendRecvSlabInY(rnk_send, rnk_recv, slab, direction)
    implicit none
    integer, intent(in) :: rnk_send, rnk_recv, slab, direction
    real, allocatable :: buffer(:, :, :, :)
    integer :: i, j, k
    integer :: i1_send, i2_send, j1_send, j2_send, k1_send, k2_send
    integer :: j_, dummy_, size_, ierr, mpi_tag

#ifdef MPI08
    type(MPI_STATUS) :: istat
#endif

#ifdef MPI
    integer :: istat(MPI_STATUS_SIZE)
#endif

    if ((mpi_rank .eq. rnk_send) .or. (mpi_rank .eq. rnk_recv)) then
#ifdef DEBUG
      ! check that we're not sending too much
      if (slab .ge. meshblocks(rnk_send + 1) % sy - NGHOST - 1) then
        call throwError('ERROR: cannot send more than `sy - NGHOST - 1` cells in y.')
      end if
#endif

      i1_send = 0; i2_send = meshblocks(rnk_send + 1) % sx - 1
      j1_send = 0; j2_send = slab - 1
      k1_send = 0; k2_send = meshblocks(rnk_send + 1) % sz - 1

      allocate (buffer(i1_send:i2_send, j1_send:j2_send, k1_send:k2_send, 6))
      size_ = (i2_send - i1_send + 1) * (j2_send - j1_send + 1) * (k2_send - k1_send + 1) * 6

      mpi_tag = 1

      if (mpi_rank .eq. rnk_send) then
        ! fill the buffer
        if (direction .lt. 0) then
          dummy_ = this_meshblock % ptr % sy - slab
        else
          dummy_ = 0
        end if
        do k = k1_send, k2_send
          do j = j1_send, j2_send
            do i = i1_send, i2_send
              j_ = dummy_ + j
              buffer(i, j, k, 1) = ex_back(i, j_, k)
              buffer(i, j, k, 2) = ey_back(i, j_, k)
              buffer(i, j, k, 3) = ez_back(i, j_, k)
              buffer(i, j, k, 4) = bx_back(i, j_, k)
              buffer(i, j, k, 5) = by_back(i, j_, k)
              buffer(i, j, k, 6) = bz_back(i, j_, k)
            end do
          end do
        end do
        call MPI_SEND(buffer, size_, default_mpi_real, rnk_recv, mpi_tag, MPI_COMM_WORLD, ierr)
      else if (mpi_rank .eq. rnk_recv) then
        call MPI_RECV(buffer, size_, default_mpi_real, rnk_send, mpi_tag, MPI_COMM_WORLD, istat, ierr)
        ! extract from the buffer
        if (direction .lt. 0) then
          dummy_ = 0
        else
          dummy_ = this_meshblock % ptr % sy
        end if
        do k = k1_send, k2_send
          do j = j1_send, j2_send
            do i = i1_send, i2_send
              j_ = dummy_ + j
              ex(i, j_, k) = buffer(i, j, k, 1)
              ey(i, j_, k) = buffer(i, j, k, 2)
              ez(i, j_, k) = buffer(i, j, k, 3)
              bx(i, j_, k) = buffer(i, j, k, 4)
              by(i, j_, k) = buffer(i, j, k, 5)
              bz(i, j_, k) = buffer(i, j, k, 6)
            end do
          end do
        end do
      end if
      deallocate (buffer)
    end if
  end subroutine SendRecvSlabInY

  subroutine SendRecvSlabInZ(rnk_send, rnk_recv, slab, direction)
    implicit none
    integer, intent(in) :: rnk_send, rnk_recv, slab, direction
    real, allocatable :: buffer(:, :, :, :)
    integer :: i, j, k
    integer :: i1_send, i2_send, j1_send, j2_send, k1_send, k2_send
    integer :: k_, dummy_, size_, ierr, mpi_tag

#ifdef MPI08
    type(MPI_STATUS) :: istat
#endif

#ifdef MPI
    integer :: istat(MPI_STATUS_SIZE)
#endif

    if ((mpi_rank .eq. rnk_send) .or. (mpi_rank .eq. rnk_recv)) then
#ifdef DEBUG
      ! check that we're not sending too much
      if (slab .ge. meshblocks(rnk_send + 1) % sz - NGHOST - 1) then
        call throwError('ERROR: cannot send more than `sz - NGHOST - 1` cells in z.')
      end if
#endif

      i1_send = 0; i2_send = meshblocks(rnk_send + 1) % sx - 1
      j1_send = 0; j2_send = meshblocks(rnk_send + 1) % sy - 1
      k1_send = 0; k2_send = slab - 1

      allocate (buffer(i1_send:i2_send, j1_send:j2_send, k1_send:k2_send, 6))
      size_ = (i2_send - i1_send + 1) * (j2_send - j1_send + 1) * (k2_send - k1_send + 1) * 6

      mpi_tag = 1

      if (mpi_rank .eq. rnk_send) then
        ! fill the buffer
        if (direction .lt. 0) then
          dummy_ = this_meshblock % ptr % sz - slab
        else
          dummy_ = 0
        end if
        do k = k1_send, k2_send
          do j = j1_send, j2_send
            do i = i1_send, i2_send
              k_ = dummy_ + k
              buffer(i, j, k, 1) = ex_back(i, j, k_)
              buffer(i, j, k, 2) = ey_back(i, j, k_)
              buffer(i, j, k, 3) = ez_back(i, j, k_)
              buffer(i, j, k, 4) = bx_back(i, j, k_)
              buffer(i, j, k, 5) = by_back(i, j, k_)
              buffer(i, j, k, 6) = bz_back(i, j, k_)
            end do
          end do
        end do
        call MPI_SEND(buffer, size_, default_mpi_real, rnk_recv, mpi_tag, MPI_COMM_WORLD, ierr)
      else if (mpi_rank .eq. rnk_recv) then
        call MPI_RECV(buffer, size_, default_mpi_real, rnk_send, mpi_tag, MPI_COMM_WORLD, istat, ierr)
        ! extract from the buffer
        if (direction .lt. 0) then
          dummy_ = 0
        else
          dummy_ = this_meshblock % ptr % sz
        end if
        do k = k1_send, k2_send
          do j = j1_send, j2_send
            do i = i1_send, i2_send
              k_ = dummy_ + k
              ex(i, j, k_) = buffer(i, j, k, 1)
              ey(i, j, k_) = buffer(i, j, k, 2)
              ez(i, j, k_) = buffer(i, j, k, 3)
              bx(i, j, k_) = buffer(i, j, k, 4)
              by(i, j, k_) = buffer(i, j, k, 5)
              bz(i, j, k_) = buffer(i, j, k, 6)
            end do
          end do
        end do
      end if
      deallocate (buffer)
    end if
  end subroutine SendRecvSlabInZ

end module m_exchangefields
