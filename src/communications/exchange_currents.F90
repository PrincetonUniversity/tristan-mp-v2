module m_exchangecurrents
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
  subroutine findCnt(ind1, ind2, ind3, fill_ghosts, send_cnt)
    integer :: imin, imax, jmin, jmax, kmin, kmax
    integer, intent(in) :: ind1, ind2, ind3
    logical, intent(in) :: fill_ghosts
    integer, intent(out) :: send_cnt
    imin = -2 * NGHOST
    imax = -2 * NGHOST
    jmin = -2 * NGHOST
    jmax = -2 * NGHOST
    kmin = -2 * NGHOST
    kmax = -2 * NGHOST

    ! highlight the region to send and save to `send_EB`
    if (.not. fill_ghosts) then
      !   sending ghost zones + normal zones
      if (ind1 .eq. 0) then
        imin = -NGHOST; imax = this_meshblock % ptr % sx + NGHOST - 1
      else if (ind1 .eq. -1) then
        imin = -NGHOST; imax = NGHOST - 1
      else if (ind1 .eq. 1) then
        imin = this_meshblock % ptr % sx - NGHOST; imax = this_meshblock % ptr % sx + NGHOST - 1
      end if

      if (ind2 .eq. 0) then
        jmin = -NGHOST; jmax = this_meshblock % ptr % sy + NGHOST - 1
      else if (ind2 .eq. -1) then
        jmin = -NGHOST; jmax = NGHOST - 1
      else if (ind2 .eq. 1) then
        jmin = this_meshblock % ptr % sy - NGHOST; jmax = this_meshblock % ptr % sy + NGHOST - 1
      end if

      if (ind3 .eq. 0) then
        kmin = -NGHOST; kmax = this_meshblock % ptr % sz + NGHOST - 1
      else if (ind3 .eq. -1) then
        kmin = -NGHOST; kmax = NGHOST - 1
      else if (ind3 .eq. 1) then
        kmin = this_meshblock % ptr % sz - NGHOST; kmax = this_meshblock % ptr % sz + NGHOST - 1
      end if
#ifdef oneD
      jmin = 0; jmax = 0
      kmin = 0; kmax = 0
#elif defined(twoD)
      kmin = 0; kmax = 0
#endif
    else
      !   sending just the normal zones
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
#elif defined(twoD)
      kmin = 0; kmax = 0
#endif
    end if

    if ((imin .eq. -2 * NGHOST) .or. (imax .eq. -2 * NGHOST) .or. &
        (jmin .eq. -2 * NGHOST) .or. (jmax .eq. -2 * NGHOST) .or. &
        (kmin .eq. -2 * NGHOST) .or. (kmax .eq. -2 * NGHOST)) then
      call throwError("Error: invalid index evaluation in findCnt")
    end if

    send_cnt = 3 * (imax - imin + 1) * (jmax - jmin + 1) * (kmax - kmin + 1)
  end subroutine findCnt

  subroutine bufferSendArray(ind1, ind2, ind3, fill_ghosts, send_cnt)
    integer :: imin, imax, jmin, jmax, kmin, kmax, i, j, k, idiff
    integer, intent(in) :: ind1, ind2, ind3
    logical, intent(in) :: fill_ghosts
    integer, intent(out) :: send_cnt
    imin = -2 * NGHOST
    imax = -2 * NGHOST
    jmin = -2 * NGHOST
    jmax = -2 * NGHOST
    kmin = -2 * NGHOST
    kmax = -2 * NGHOST

    ! highlight the region to send and save to `send_EB`
    if (.not. fill_ghosts) then
      !   sending ghost zones + normal zones
      if (ind1 .eq. 0) then
        imin = -NGHOST; imax = this_meshblock % ptr % sx + NGHOST - 1
      else if (ind1 .eq. -1) then
        imin = -NGHOST; imax = NGHOST - 1
      else if (ind1 .eq. 1) then
        imin = this_meshblock % ptr % sx - NGHOST; imax = this_meshblock % ptr % sx + NGHOST - 1
      end if

      if (ind2 .eq. 0) then
        jmin = -NGHOST; jmax = this_meshblock % ptr % sy + NGHOST - 1
      else if (ind2 .eq. -1) then
        jmin = -NGHOST; jmax = NGHOST - 1
      else if (ind2 .eq. 1) then
        jmin = this_meshblock % ptr % sy - NGHOST; jmax = this_meshblock % ptr % sy + NGHOST - 1
      end if

      if (ind3 .eq. 0) then
        kmin = -NGHOST; kmax = this_meshblock % ptr % sz + NGHOST - 1
      else if (ind3 .eq. -1) then
        kmin = -NGHOST; kmax = NGHOST - 1
      else if (ind3 .eq. 1) then
        kmin = this_meshblock % ptr % sz - NGHOST; kmax = this_meshblock % ptr % sz + NGHOST - 1
      end if
#ifdef oneD
      jmin = 0; jmax = 0
      kmin = 0; kmax = 0
#elif defined(twoD)
      kmin = 0; kmax = 0
#endif
    else
      !   sending just the normal zones
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
#elif defined(twoD)
      kmin = 0; kmax = 0
#endif
    end if

    if ((imin .eq. -2 * NGHOST) .or. (imax .eq. -2 * NGHOST) .or. &
        (jmin .eq. -2 * NGHOST) .or. (jmax .eq. -2 * NGHOST) .or. &
        (kmin .eq. -2 * NGHOST) .or. (kmax .eq. -2 * NGHOST)) then
      call throwError("Error: invalid index evaluation in bufferSendArray")
    end if

    send_cnt = 1
    idiff = imax - imin
    do k = kmin, kmax
      do j = jmin, jmax
        send_EB(send_cnt:send_cnt + idiff) = jx(imin:imax, j, k)
        send_cnt = send_cnt + idiff + 1
        send_EB(send_cnt:send_cnt + idiff) = jy(imin:imax, j, k)
        send_cnt = send_cnt + idiff + 1
        send_EB(send_cnt:send_cnt + idiff) = jz(imin:imax, j, k)
        send_cnt = send_cnt + idiff + 1
      end do
    end do
    send_cnt = send_cnt - 1
  end subroutine bufferSendArray

  subroutine extractRecvArray(ind1, ind2, ind3, fill_ghosts)
    implicit none
    integer :: imin, imax, jmin, jmax, kmin, kmax, i, j, k, idiff
    integer :: send_cnt
    integer, intent(in) :: ind1, ind2, ind3
    logical, intent(in) :: fill_ghosts
    imin = -2 * NGHOST
    imax = -2 * NGHOST
    jmin = -2 * NGHOST
    jmax = -2 * NGHOST
    kmin = -2 * NGHOST
    kmax = -2 * NGHOST

    if (.not. fill_ghosts) then
      !   write to ghosts + normal zones
      if (ind1 .eq. 0) then
        imin = -NGHOST; imax = this_meshblock % ptr % sx + NGHOST - 1
      else if (ind1 .eq. -1) then
        imin = -NGHOST; imax = NGHOST - 1
      else if (ind1 .eq. 1) then
        imin = this_meshblock % ptr % sx - NGHOST; imax = this_meshblock % ptr % sx + NGHOST - 1
      end if

      if (ind2 .eq. 0) then
        jmin = -NGHOST; jmax = this_meshblock % ptr % sy + NGHOST - 1
      else if (ind2 .eq. -1) then
        jmin = -NGHOST; jmax = NGHOST - 1
      else if (ind2 .eq. 1) then
        jmin = this_meshblock % ptr % sy - NGHOST; jmax = this_meshblock % ptr % sy + NGHOST - 1
      end if

      if (ind3 .eq. 0) then
        kmin = -NGHOST; kmax = this_meshblock % ptr % sz + NGHOST - 1
      else if (ind3 .eq. -1) then
        kmin = -NGHOST; kmax = NGHOST - 1
      else if (ind3 .eq. 1) then
        kmin = this_meshblock % ptr % sz - NGHOST; kmax = this_meshblock % ptr % sz + NGHOST - 1
      end if
#ifdef oneD
      jmin = 0; jmax = 0
      kmin = 0; kmax = 0
#elif defined(twoD)
      kmin = 0; kmax = 0
#endif
    else
      !   write to ghosts
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
#elif defined(twoD)
      kmin = 0; kmax = 0
#endif
    end if
    if ((imin .eq. -2 * NGHOST) .or. (imax .eq. -2 * NGHOST) .or. &
        (jmin .eq. -2 * NGHOST) .or. (jmax .eq. -2 * NGHOST) .or. &
        (kmin .eq. -2 * NGHOST) .or. (kmax .eq. -2 * NGHOST)) then
      call throwError("Error: invalid index evaluation in extractRecvArray")
    end if

    send_cnt = 1
    idiff = imax - imin
    do k = kmin, kmax
      do j = jmin, jmax
        if (.not. fill_ghosts) then
          ! add to existing values
          jx_buff(imin:imax, j, k) = jx_buff(imin:imax, j, k) + recv_fld(send_cnt:send_cnt + idiff)
          send_cnt = send_cnt + idiff + 1
          jy_buff(imin:imax, j, k) = jy_buff(imin:imax, j, k) + recv_fld(send_cnt:send_cnt + idiff)
          send_cnt = send_cnt + idiff + 1
          jz_buff(imin:imax, j, k) = jz_buff(imin:imax, j, k) + recv_fld(send_cnt:send_cnt + idiff)
          send_cnt = send_cnt + idiff + 1
        else
          ! overwrite the existing values
          jx(imin:imax, j, k) = recv_fld(send_cnt:send_cnt + idiff)
          send_cnt = send_cnt + idiff + 1
          jy(imin:imax, j, k) = recv_fld(send_cnt:send_cnt + idiff)
          send_cnt = send_cnt + idiff + 1
          jz(imin:imax, j, k) = recv_fld(send_cnt:send_cnt + idiff)
          send_cnt = send_cnt + idiff + 1
        end if
      end do
    end do

  end subroutine extractRecvArray

  subroutine exchangeCurrents(fill_ghosts_Q)
    implicit none
    integer :: ind1, ind2, ind3
    integer :: cnt, ierr
    integer :: mpi_sendto, mpi_recvfrom, mpi_tag
    logical :: fill_ghosts
    logical :: should_send, should_recv
    logical, optional, intent(in) :: fill_ghosts_Q

#ifdef MPI08
    type(MPI_STATUS) :: istat
#endif

#ifdef MPI
    integer :: istat(MPI_STATUS_SIZE)
#endif

    if (present(fill_ghosts_Q)) then
      fill_ghosts = fill_ghosts_Q
    else
      fill_ghosts = .false.
    end if

    if (.not. fill_ghosts) then
      jx_buff(:, :, :) = 0.0
      jy_buff(:, :, :) = 0.0
      jz_buff(:, :, :) = 0.0
    end if

    do ind3 = -1, 1
      do ind2 = -1, 1
        do ind1 = -1, 1
          if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
#ifdef oneD
          if ((ind2 .ne. 0) .or. (ind3 .ne. 0)) cycle
#elif defined(twoD)
          if (ind3 .ne. 0) cycle
#endif

          mpi_tag = (ind3 + 2) + 3 * (ind2 + 1) + 9 * (ind1 + 1) + 50

          should_send = associated(this_meshblock % ptr % neighbor(ind1, ind2, ind3) % ptr)
          should_recv = associated(this_meshblock % ptr % neighbor(-ind1, -ind2, -ind3) % ptr)
          if (should_send .and. should_recv) then

            mpi_sendto = this_meshblock % ptr % neighbor(ind1, ind2, ind3) % ptr % rnk
            mpi_recvfrom = this_meshblock % ptr % neighbor(-ind1, -ind2, -ind3) % ptr % rnk

            call bufferSendArray(ind1, ind2, ind3, fill_ghosts, cnt)
            call MPI_SENDRECV(send_EB(1:cnt), cnt, default_mpi_real, mpi_sendto, mpi_tag, &
                              recv_fld(1:cnt), cnt, default_mpi_real, mpi_recvfrom, mpi_tag, &
                              MPI_COMM_WORLD, istat, ierr)
            call extractRecvArray(-ind1, -ind2, -ind3, fill_ghosts)
          else if ((.not. should_send) .and. should_recv) then
            mpi_recvfrom = this_meshblock % ptr % neighbor(-ind1, -ind2, -ind3) % ptr % rnk
            call findCnt(-ind1, -ind2, -ind3, fill_ghosts, cnt)
            call MPI_RECV(recv_fld(1:cnt), cnt, default_mpi_real, mpi_recvfrom, mpi_tag, MPI_COMM_WORLD, istat, ierr)
            call extractRecvArray(-ind1, -ind2, -ind3, fill_ghosts)
          else if ((.not. should_recv) .and. should_send) then
            mpi_sendto = this_meshblock % ptr % neighbor(ind1, ind2, ind3) % ptr % rnk
            call bufferSendArray(ind1, ind2, ind3, fill_ghosts, cnt)
            call MPI_SEND(send_EB(1:cnt), cnt, default_mpi_real, mpi_sendto, mpi_tag, MPI_COMM_WORLD, ierr)
          end if
        end do
      end do
    end do

    if (.not. fill_ghosts) then
      jx(:, :, :) = jx(:, :, :) + jx_buff(:, :, :)
      jy(:, :, :) = jy(:, :, :) + jy_buff(:, :, :)
      jz(:, :, :) = jz(:, :, :) + jz_buff(:, :, :)
    end if

    call printDiag("exchangeCurrents()", 2)
  end subroutine exchangeCurrents

#else

  ! non-blocking MPI communication
  subroutine exchangeCurrents(fill_ghosts_Q)
    implicit none
    integer :: i, j, k, imin, imax, jmin, jmax, kmin, kmax, idiff, loc
    integer :: ind1, ind2, ind3, cntr, n_cntr
    integer :: send_cnt, recv_cnt, ierr
    integer :: mpi_sendto, mpi_recvfrom, mpi_sendtag, mpi_recvtag
    integer :: mpi_offset
    logical :: fill_ghosts
    logical, optional, intent(in) :: fill_ghosts_Q

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

    if (present(fill_ghosts_Q)) then
      fill_ghosts = fill_ghosts_Q
    else
      fill_ghosts = .false.
    end if

    allocate (mpi_req(sendrecv_neighbors))
    allocate (mpi_sendflags(sendrecv_neighbors))
    allocate (mpi_recvflags(sendrecv_neighbors))

    if (.not. fill_ghosts) then
      jx_buff(:, :, :) = 0.0
      jy_buff(:, :, :) = 0.0
      jz_buff(:, :, :) = 0.0
    end if

    ! exchange the ghost cells + the real cells
    cntr = 0
    do ind3 = -1, 1
      do ind2 = -1, 1
        do ind1 = -1, 1
          if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
#ifdef oneD
          if ((ind2 .ne. 0) .or. (ind3 .ne. 0)) cycle
#elif defined(twoD)
          if (ind3 .ne. 0) cycle
#endif
          if (.not. associated(this_meshblock % ptr % neighbor(ind1, ind2, ind3) % ptr)) cycle
          cntr = cntr + 1

          mpi_sendto = this_meshblock % ptr % neighbor(ind1, ind2, ind3) % ptr % rnk
          mpi_sendtag = (ind3 + 2) + 3 * (ind2 + 1) + 9 * (ind1 + 1)

          ! highlight the region to send and save to `send_fld`
          if (.not. fill_ghosts) then
            !   sending ghost zones + normal zones
            if (ind1 .eq. 0) then
              imin = -NGHOST; imax = this_meshblock % ptr % sx + NGHOST - 1
            else if (ind1 .eq. -1) then
              imin = -NGHOST; imax = NGHOST - 1
            else if (ind1 .eq. 1) then
              imin = this_meshblock % ptr % sx - NGHOST; imax = this_meshblock % ptr % sx + NGHOST - 1
            end if
            if (ind2 .eq. 0) then
              jmin = -NGHOST; jmax = this_meshblock % ptr % sy + NGHOST - 1
            else if (ind2 .eq. -1) then
              jmin = -NGHOST; jmax = NGHOST - 1
            else if (ind2 .eq. 1) then
              jmin = this_meshblock % ptr % sy - NGHOST; jmax = this_meshblock % ptr % sy + NGHOST - 1
            end if
            if (ind3 .eq. 0) then
              kmin = -NGHOST; kmax = this_meshblock % ptr % sz + NGHOST - 1
            else if (ind3 .eq. -1) then
              kmin = -NGHOST; kmax = NGHOST - 1
            else if (ind3 .eq. 1) then
              kmin = this_meshblock % ptr % sz - NGHOST; kmax = this_meshblock % ptr % sz + NGHOST - 1
            end if
#ifdef oneD
            jmin = 0; jmax = 0
            kmin = 0; kmax = 0
#elif defined(twoD)
            kmin = 0; kmax = 0
#endif
          else
            !   sending just the normal zones
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
#elif defined(twoD)
            kmin = 0; kmax = 0
#endif
          end if

          ! write send/recv arrays in/from a given direction
          !     in 3D: 26 directions, in 2D: 8, in 1D: 2
          mpi_offset = (cntr - 1) * sendrecv_offsetsz
          send_cnt = 1
          idiff = imax - imin
          loc = mpi_offset + send_cnt
          do k = kmin, kmax
            do j = jmin, jmax
              send_fld(loc:loc + idiff) = jx(imin:imax, j, k)
              loc = loc + idiff + 1
              send_fld(loc:loc + idiff) = jy(imin:imax, j, k)
              loc = loc + idiff + 1
              send_fld(loc:loc + idiff) = jz(imin:imax, j, k)
              loc = loc + idiff + 1
              send_cnt = send_cnt + 3 * (idiff + 1)
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
      do ind3 = -1, 1
        do ind2 = -1, 1
          do ind1 = -1, 1
            if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
#ifdef oneD
            if ((ind2 .ne. 0) .or. (ind3 .ne. 0)) cycle
#elif defined(twoD)
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
                if (.not. fill_ghosts) then
                  !   write to ghosts + normal zones
                  if (ind1 .eq. 0) then
                    imin = -NGHOST; imax = this_meshblock % ptr % sx + NGHOST - 1
                  else if (ind1 .eq. -1) then
                    imin = -NGHOST; imax = NGHOST - 1
                  else if (ind1 .eq. 1) then
                    imin = this_meshblock % ptr % sx - NGHOST; imax = this_meshblock % ptr % sx + NGHOST - 1
                  end if
                  if (ind2 .eq. 0) then
                    jmin = -NGHOST; jmax = this_meshblock % ptr % sy + NGHOST - 1
                  else if (ind2 .eq. -1) then
                    jmin = -NGHOST; jmax = NGHOST - 1
                  else if (ind2 .eq. 1) then
                    jmin = this_meshblock % ptr % sy - NGHOST; jmax = this_meshblock % ptr % sy + NGHOST - 1
                  end if
                  if (ind3 .eq. 0) then
                    kmin = -NGHOST; kmax = this_meshblock % ptr % sz + NGHOST - 1
                  else if (ind3 .eq. -1) then
                    kmin = -NGHOST; kmax = NGHOST - 1
                  else if (ind3 .eq. 1) then
                    kmin = this_meshblock % ptr % sz - NGHOST; kmax = this_meshblock % ptr % sz + NGHOST - 1
                  end if
#ifdef oneD
                  jmin = 0; jmax = 0
                  kmin = 0; kmax = 0
#elif defined(twoD)
                  kmin = 0; kmax = 0
#endif
                else
                  !   write to ghosts
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
#elif defined(twoD)
                  kmin = 0; kmax = 0
#endif
                end if

                ! copy `recv_fld` to ghost cells and active cells
                send_cnt = 1
                idiff = imax - imin
                do k = kmin, kmax
                  do j = jmin, jmax
                    if (.not. fill_ghosts) then
                      ! add to existing values
                      jx_buff(imin:imax, j, k) = jx_buff(imin:imax, j, k) + recv_fld(send_cnt:send_cnt + idiff)
                      send_cnt = send_cnt + idiff + 1
                      jy_buff(imin:imax, j, k) = jy_buff(imin:imax, j, k) + recv_fld(send_cnt:send_cnt + idiff)
                      send_cnt = send_cnt + idiff + 1
                      jz_buff(imin:imax, j, k) = jz_buff(imin:imax, j, k) + recv_fld(send_cnt:send_cnt + idiff)
                      send_cnt = send_cnt + idiff + 1
                    else
                      ! overwrite the existing values
                      jx(imin:imax, j, k) = recv_fld(send_cnt:send_cnt + idiff)
                      send_cnt = send_cnt + idiff + 1
                      jy(imin:imax, j, k) = recv_fld(send_cnt:send_cnt + idiff)
                      send_cnt = send_cnt + idiff + 1
                      jz(imin:imax, j, k) = recv_fld(send_cnt:send_cnt + idiff)
                      send_cnt = send_cnt + idiff + 1
                    end if
                  end do
                end do

              end if ! if receive is ready
            end if ! if receive hasn't been done yet
          end do ! ind3
        end do ! ind2
      end do ! ind1
    end do ! global loop

    if (.not. fill_ghosts) then
      jx(:, :, :) = jx(:, :, :) + jx_buff(:, :, :)
      jy(:, :, :) = jy(:, :, :) + jy_buff(:, :, :)
      jz(:, :, :) = jz(:, :, :) + jz_buff(:, :, :)
    end if

    call printDiag("exchangeCurrents()", 2)
  end subroutine exchangeCurrents
#endif
end module m_exchangecurrents
