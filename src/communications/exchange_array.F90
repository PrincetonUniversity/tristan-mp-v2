#include "../defs.F90"

module m_exchangearray
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_fields
contains
#ifndef MPINONBLOCK
  ! blocking MPI communication
  subroutine findCnt(ind1, ind2, ind3, send_cnt)
    integer                       :: imin, imax, jmin, jmax, kmin, kmax, i, j, k
    integer, intent(in)           :: ind1, ind2, ind3    
    integer, intent(out)          :: send_cnt 

    ! highlight the region to send and save to `send_EB`
    if (ind1 .eq. 0) then
       imin = -NGHOST; imax = this_meshblock%ptr%sx + NGHOST - 1
    else if (ind1 .eq. -1) then
       imin = -NGHOST; imax = NGHOST - 1
    else if (ind1 .eq. 1) then
       imin = this_meshblock%ptr%sx - NGHOST; imax = this_meshblock%ptr%sx + NGHOST - 1
    end if
    
    if (ind2 .eq. 0) then
       jmin = -NGHOST; jmax = this_meshblock%ptr%sy + NGHOST - 1
    else if (ind2 .eq. -1) then
       jmin = -NGHOST; jmax = NGHOST - 1
    else if (ind2 .eq. 1) then
       jmin = this_meshblock%ptr%sy - NGHOST; jmax = this_meshblock%ptr%sy + NGHOST - 1
    end if

    if (ind3 .eq. 0) then
       kmin = -NGHOST; kmax = this_meshblock%ptr%sz + NGHOST - 1
    else if (ind3 .eq. -1) then
       kmin = -NGHOST; kmax = NGHOST - 1
    else if (ind3 .eq. 1) then
       kmin = this_meshblock%ptr%sz - NGHOST; kmax = this_meshblock%ptr%sz + NGHOST - 1
    end if
    
#ifdef oneD
    jmin = 0; jmax = 0
    kmin = 0; kmax = 0
#elif twoD
    kmin = 0; kmax = 0
#endif
            
    send_cnt = (imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)
  end subroutine findCnt

  subroutine bufferSendArray( ind1, ind2, ind3, send_cnt)
    integer                       :: imin, imax, jmin, jmax, kmin, kmax, i, j, k
    integer, intent(in)           :: ind1, ind2, ind3
    integer, intent(out)          :: send_cnt

    ! highlight the region to send and save to `send_fld`
    if (ind1 .eq. 0) then
       imin = -NGHOST; imax = this_meshblock%ptr%sx + NGHOST - 1
    else if (ind1 .eq. -1) then
       imin = -NGHOST; imax = NGHOST - 1
    else if (ind1 .eq. 1) then
       imin = this_meshblock%ptr%sx - NGHOST; imax = this_meshblock%ptr%sx + NGHOST - 1
    end if
    
    if (ind2 .eq. 0) then
       jmin = -NGHOST; jmax = this_meshblock%ptr%sy + NGHOST - 1
    else if (ind2 .eq. -1) then
       jmin = -NGHOST; jmax = NGHOST - 1
    else if (ind2 .eq. 1) then
       jmin = this_meshblock%ptr%sy - NGHOST; jmax = this_meshblock%ptr%sy + NGHOST - 1
    end if

    if (ind3 .eq. 0) then
       kmin = -NGHOST; kmax = this_meshblock%ptr%sz + NGHOST - 1
    else if (ind3 .eq. -1) then
       kmin = -NGHOST; kmax = NGHOST - 1
    else if (ind3 .eq. 1) then
       kmin = this_meshblock%ptr%sz - NGHOST; kmax = this_meshblock%ptr%sz + NGHOST - 1
    end if
    
#ifdef oneD
    jmin = 0; jmax = 0
    kmin = 0; kmax = 0
#elif twoD
    kmin = 0; kmax = 0
#endif

    send_cnt = 1
    do i = imin, imax
       do j = jmin, jmax
          do k = kmin, kmax
             send_EB(send_cnt + 0) = lg_arr(i, j, k)
             send_cnt = send_cnt + 1
          end do
       end do
    end do
    send_cnt = send_cnt - 1
  end subroutine bufferSendArray

  subroutine extractRecvArray(ind1, ind2, ind3)
    implicit none
    integer                       :: imin, imax, jmin, jmax, kmin, kmax, i, j, k
    integer                       :: send_cnt
    integer, intent(in)           :: ind1, ind2, ind3


    ! write received data to local memory
    ! highlight the region to extract the `recv_fld`
    !   write to ghosts + normal zones
    if (ind1 .eq. 0) then
       imin = -NGHOST; imax = this_meshblock%ptr%sx + NGHOST - 1
    else if (ind1 .eq. -1) then
       imin = -NGHOST; imax = NGHOST - 1
    else if (ind1 .eq. 1) then
       imin = this_meshblock%ptr%sx - NGHOST; imax = this_meshblock%ptr%sx + NGHOST - 1
    end if

    if (ind2 .eq. 0) then
       jmin = -NGHOST; jmax = this_meshblock%ptr%sy + NGHOST - 1
    else if (ind2 .eq. -1) then
       jmin = -NGHOST; jmax = NGHOST - 1
    else if (ind2 .eq. 1) then
       jmin = this_meshblock%ptr%sy - NGHOST; jmax = this_meshblock%ptr%sy + NGHOST - 1
    end if
    
    if (ind3 .eq. 0) then
       kmin = -NGHOST; kmax = this_meshblock%ptr%sz + NGHOST - 1
    else if (ind3 .eq. -1) then
       kmin = -NGHOST; kmax = NGHOST - 1
    else if (ind3 .eq. 1) then
       kmin = this_meshblock%ptr%sz - NGHOST; kmax = this_meshblock%ptr%sz + NGHOST - 1
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
             jx_buff(i, j, k) = jx_buff(i, j, k) + recv_fld(send_cnt)
             send_cnt = send_cnt + 1
          end do
       end do
    end do
    
  end subroutine extractRecvArray

    
  subroutine exchangeArray()
    implicit none
    integer           :: i, j, k, imin, imax, jmin, jmax, kmin, kmax
    integer           :: ind1, ind2, ind3, cntr, n_cntr
    integer           :: cnt, ierr
    integer           :: mpi_sendto, mpi_recvfrom, mpi_tag
    integer           :: mpi_offset
    logical           :: should_send, should_recv

    #ifdef MPI08
      type(MPI_REQUEST), allocatable  :: mpi_req(:)
      type(MPI_STATUS)                :: istat
    #endif

    #ifdef MPI
      integer, allocatable            :: mpi_req(:)
      integer                         :: istat(MPI_STATUS_SIZE)
    #endif

    jx_buff(:,:,:) = 0.0

    do ind1 = -1, 1
      do ind2 = -1, 1
        do ind3 = -1, 1
          if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
          #ifdef oneD
            if ((ind2 .ne. 0) .or. (ind3 .ne. 0)) cycle
          #elif twoD
            if (ind3 .ne. 0) cycle
          #endif

          mpi_tag = (ind3 + 2) + 3 * (ind2 + 1) + 9 * (ind1 + 1) + 300

          should_send = associated(this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr)
          should_recv = associated(this_meshblock%ptr%neighbor(-ind1,-ind2,-ind3)%ptr)
          if (should_send .and. should_recv) then

             mpi_sendto = this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr%rnk
             mpi_recvfrom = this_meshblock%ptr%neighbor(-ind1,-ind2,-ind3)%ptr%rnk

             call bufferSendArray(ind1, ind2, ind3, cnt)
             call MPI_SENDRECV(send_EB(1 : cnt), cnt, MPI_REAL, mpi_sendto, mpi_tag,&
                  & recv_fld(1 : cnt), cnt, MPI_REAL, mpi_recvfrom, mpi_tag,&
                  & MPI_COMM_WORLD, istat, ierr)
             call extractRecvArray(-ind1, -ind2, -ind3)
          else if ((.not. should_send) .and. should_recv) then
            mpi_recvfrom = this_meshblock%ptr%neighbor(-ind1,-ind2,-ind3)%ptr%rnk 
            call findCnt(ind1, ind2, ind3, cnt)
            call MPI_RECV(recv_fld(1 : cnt),cnt, MPI_REAL,mpi_recvfrom, mpi_tag, MPI_COMM_WORLD, istat, ierr)
            call extractRecvArray(-ind1, -ind2, -ind3)
         else if ((.not. should_recv) .and. should_send) then
            mpi_sendto = this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr%rnk 
            call bufferSendArray(ind1, ind2, ind3, cnt)
            call MPI_SEND(send_EB(1 : cnt), cnt, MPI_REAL,mpi_sendto, mpi_tag, MPI_COMM_WORLD, ierr)
         end if
        end do
      end do
    end do

    lg_arr(:,:,:) = lg_arr(:,:,:) + jx_buff(:,:,:)

    call printDiag("exchangeArray()", 3)
  end subroutine exchangeArray

#else
  ! non-blocking MPI communication

  subroutine exchangeArray()
    implicit none
    integer            :: i, j, k, imin, imax, jmin, jmax, kmin, kmax
    integer            :: ind1, ind2, ind3, cntr, n_cntr
    integer            :: send_cnt, recv_cnt, ierr
    integer            :: mpi_sendto, mpi_recvfrom, mpi_sendtag, mpi_recvtag
    integer            :: mpi_offset

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

    ! using this temporary array
    jx_buff(:, :, :) = 0.0

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

          mpi_sendto = this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr%rnk
          mpi_sendtag = (ind3 + 2) + 3 * (ind2 + 1) + 9 * (ind1 + 1)

          ! highlight the region to send and save to `send_fld`
          if (ind1 .eq. 0) then
            imin = -NGHOST; imax = this_meshblock%ptr%sx + NGHOST - 1
          else if (ind1 .eq. -1) then
            imin = -NGHOST; imax = NGHOST - 1
          else if (ind1 .eq. 1) then
            imin = this_meshblock%ptr%sx - NGHOST; imax = this_meshblock%ptr%sx + NGHOST - 1
          end if
          if (ind2 .eq. 0) then
            jmin = -NGHOST; jmax = this_meshblock%ptr%sy + NGHOST - 1
          else if (ind2 .eq. -1) then
            jmin = -NGHOST; jmax = NGHOST - 1
          else if (ind2 .eq. 1) then
            jmin = this_meshblock%ptr%sy - NGHOST; jmax = this_meshblock%ptr%sy + NGHOST - 1
          end if
          if (ind3 .eq. 0) then
            kmin = -NGHOST; kmax = this_meshblock%ptr%sz + NGHOST - 1
          else if (ind3 .eq. -1) then
            kmin = -NGHOST; kmax = NGHOST - 1
          else if (ind3 .eq. 1) then
            kmin = this_meshblock%ptr%sz - NGHOST; kmax = this_meshblock%ptr%sz + NGHOST - 1
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
                send_fld(mpi_offset + send_cnt) = lg_arr(i, j, k)
                send_cnt = send_cnt + 1
              end do
            end do
          end do
          send_cnt = send_cnt - 1

          ! post non-blocking send requests
          call MPI_ISEND(send_fld(mpi_offset + 1 : mpi_offset + send_cnt), send_cnt, MPI_REAL,&
          & mpi_sendto, mpi_sendtag, MPI_COMM_WORLD, mpi_req(cntr), ierr)
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
                ! if the message is ready to be received -> get the size & receive it
                call MPI_GET_COUNT(istat, MPI_INTEGER, recv_cnt, ierr)
                if (recv_cnt .lt. 0) then
                  print *, "exchange_array", recv_cnt, mpi_rank
                end if
                call MPI_RECV(recv_fld(1:recv_cnt), recv_cnt, MPI_REAL,&
                & mpi_recvfrom, mpi_recvtag, MPI_COMM_WORLD, istat, ierr)

                ! write received data to local memory
                ! highlight the region to extract the `recv_fld`
                !   write to ghosts + normal zones
                if (ind1 .eq. 0) then
                  imin = -NGHOST; imax = this_meshblock%ptr%sx + NGHOST - 1
                else if (ind1 .eq. -1) then
                  imin = -NGHOST; imax = NGHOST - 1
                else if (ind1 .eq. 1) then
                  imin = this_meshblock%ptr%sx - NGHOST; imax = this_meshblock%ptr%sx + NGHOST - 1
                end if
                if (ind2 .eq. 0) then
                  jmin = -NGHOST; jmax = this_meshblock%ptr%sy + NGHOST - 1
                else if (ind2 .eq. -1) then
                  jmin = -NGHOST; jmax = NGHOST - 1
                else if (ind2 .eq. 1) then
                  jmin = this_meshblock%ptr%sy - NGHOST; jmax = this_meshblock%ptr%sy + NGHOST - 1
                end if
                if (ind3 .eq. 0) then
                  kmin = -NGHOST; kmax = this_meshblock%ptr%sz + NGHOST - 1
                else if (ind3 .eq. -1) then
                  kmin = -NGHOST; kmax = NGHOST - 1
                else if (ind3 .eq. 1) then
                  kmin = this_meshblock%ptr%sz - NGHOST; kmax = this_meshblock%ptr%sz + NGHOST - 1
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
                      jx_buff(i, j, k) = jx_buff(i, j, k) + recv_fld(send_cnt)
                      send_cnt = send_cnt + 1
                    end do
                  end do
                end do

              end if ! if receive is ready
            end if ! if receive hasn't been done yet
          end do ! ind3
        end do ! ind2
      end do ! ind1
    end do ! global loop
    lg_arr(:,:,:) = lg_arr(:,:,:) + jx_buff(:,:,:)
    call printDiag("exchangeArray()", 3)
  end subroutine exchangeArray
#endif
end module m_exchangearray

