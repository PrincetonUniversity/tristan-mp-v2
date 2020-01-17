#include "../defs.F90"

module m_exchangecurrents
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_fields
contains
  subroutine exchangeCurrents(fill_ghosts_Q)
    implicit none
    integer           :: i, j, k, imin, imax, jmin, jmax, kmin, kmax
    integer           :: ind1, ind2, ind3, cntr, n_cntr
    integer           :: send_cnt, recv_cnt, ierr
    integer           :: mpi_sendto, mpi_recvfrom, mpi_sendtag, mpi_recvtag
    integer           :: mpi_offset
    logical           :: fill_ghosts
    logical, optional, intent(in) :: fill_ghosts_Q

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

    if (present(fill_ghosts_Q)) then
      fill_ghosts = fill_ghosts_Q
    else
      fill_ghosts = .false.
    end if

    allocate(mpi_req(sendrecv_neighbors))
    allocate(mpi_sendflags(sendrecv_neighbors))
    allocate(mpi_recvflags(sendrecv_neighbors))

    ! exchange the ghost cells + the real cells
    cntr = 0
    do ind1 = -1, 1
      do ind2 = -1, 1
        do ind3 = -1, 1
          if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
          #ifndef threeD
            if (ind3 .ne. 0) cycle
          #endif
          if (.not. associated(this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr)) cycle
          cntr = cntr + 1

          mpi_sendto = this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr%rnk
          mpi_sendtag = 100 * (mpi_rank + 1) + (ind3 + 2) + 3 * (ind2 + 1) + 9 * (ind1 + 1)

          ! highlight the region to send and save to `send_fld`
          if (.not. fill_ghosts) then
            !   sending ghost zones + normal zones
            if (ind1 .eq. 0) then
              imin = 0; imax = this_meshblock%ptr%sx - 1
            else if (ind1 .eq. -1) then
              imin = -NGHOST; imax = NGHOST - 1
            else if (ind1 .eq. 1) then
              imin = this_meshblock%ptr%sx - NGHOST; imax = this_meshblock%ptr%sx + NGHOST - 1
            end if
            if (ind2 .eq. 0) then
              jmin = 0; jmax = this_meshblock%ptr%sy - 1
            else if (ind2 .eq. -1) then
              jmin = -NGHOST; jmax = NGHOST - 1
            else if (ind2 .eq. 1) then
              jmin = this_meshblock%ptr%sy - NGHOST; jmax = this_meshblock%ptr%sy + NGHOST - 1
            end if
            if (ind3 .eq. 0) then
              kmin = 0; kmax = this_meshblock%ptr%sz - 1
            else if (ind3 .eq. -1) then
              kmin = -NGHOST; kmax = NGHOST - 1
            else if (ind3 .eq. 1) then
              kmin = this_meshblock%ptr%sz - NGHOST; kmax = this_meshblock%ptr%sz + NGHOST - 1
            end if
            #ifndef threeD
              kmin = 0; kmax = 0
            #endif
          else
            !   sending just the normal zones
            if (ind1 .eq. 0) then
              imin = 0; imax = this_meshblock%ptr%sx - 1
            else if (ind1 .eq. -1) then
              imin = 0; imax = NGHOST - 1
            else if (ind1 .eq. 1) then
              imin = this_meshblock%ptr%sx - NGHOST; imax = this_meshblock%ptr%sx - 1
            end if
            if (ind2 .eq. 0) then
              jmin = 0; jmax = this_meshblock%ptr%sy - 1
            else if (ind2 .eq. -1) then
              jmin = 0; jmax = NGHOST - 1
            else if (ind2 .eq. 1) then
              jmin = this_meshblock%ptr%sy - NGHOST; jmax = this_meshblock%ptr%sy - 1
            end if
            if (ind3 .eq. 0) then
              kmin = 0; kmax = this_meshblock%ptr%sz - 1
            else if (ind3 .eq. -1) then
              kmin = 0; kmax = NGHOST - 1
            else if (ind3 .eq. 1) then
              kmin = this_meshblock%ptr%sz - NGHOST; kmax = this_meshblock%ptr%sz - 1
            end if
            #ifndef threeD
              kmin = 0; kmax = 0
            #endif
          end if

          ! write send/recv arrays in/from a given direction
          !     in 3D: 26 directions, in 2D: 8, in 1D: 2
          mpi_offset = (cntr - 1) * sendrecv_offsetsz
          send_cnt = 1
          do i = imin, imax
            do j = jmin, jmax
              do k = kmin, kmax
                send_fld(mpi_offset + send_cnt + 0) = jx(i, j, k)
                send_fld(mpi_offset + send_cnt + 1) = jy(i, j, k)
                send_fld(mpi_offset + send_cnt + 2) = jz(i, j, k)
                send_cnt = send_cnt + 3
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
            #ifndef threeD
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
            mpi_recvtag = 100 * (mpi_recvfrom + 1) + (-ind3 + 2) + 3 * (-ind2 + 1) + 9 * (-ind1 + 1)

            if (.not. mpi_recvflags(cntr)) then
              quit_loop = .false.
              call MPI_IPROBE(mpi_recvfrom, mpi_recvtag, MPI_COMM_WORLD, mpi_recvflags(cntr), istat, ierr)
              if (mpi_recvflags(cntr)) then
                ! if the message is ready to be received -> get the size & receive it
                call MPI_GET_COUNT(istat, MPI_REAL, recv_cnt, ierr)
                call MPI_RECV(recv_fld(1:recv_cnt), recv_cnt, MPI_REAL,&
                            & mpi_recvfrom, mpi_recvtag, MPI_COMM_WORLD, istat, ierr)

                ! write received data to local memory
                ! highlight the region to extract the `recv_fld`
                if (.not. fill_ghosts) then
                  !   write to ghosts + normal zones
                  if (ind1 .eq. 0) then
                    imin = 0; imax = this_meshblock%ptr%sx - 1
                  else if (ind1 .eq. -1) then
                    imin = -NGHOST; imax = NGHOST - 1
                  else if (ind1 .eq. 1) then
                    imin = this_meshblock%ptr%sx - NGHOST; imax = this_meshblock%ptr%sx + NGHOST - 1
                  end if
                  if (ind2 .eq. 0) then
                    jmin = 0; jmax = this_meshblock%ptr%sy - 1
                  else if (ind2 .eq. -1) then
                    jmin = -NGHOST; jmax = NGHOST - 1
                  else if (ind2 .eq. 1) then
                    jmin = this_meshblock%ptr%sy - NGHOST; jmax = this_meshblock%ptr%sy + NGHOST - 1
                  end if
                  if (ind3 .eq. 0) then
                    kmin = 0; kmax = this_meshblock%ptr%sz - 1
                  else if (ind3 .eq. -1) then
                    kmin = -NGHOST; kmax = NGHOST - 1
                  else if (ind3 .eq. 1) then
                    kmin = this_meshblock%ptr%sz - NGHOST; kmax = this_meshblock%ptr%sz + NGHOST - 1
                  end if
                  #ifndef threeD
                    kmin = 0; kmax = 0
                  #endif
                else
                  !   write to ghosts
                  if (ind1 .eq. 0) then
                    imin = 0; imax = this_meshblock%ptr%sx - 1
                  else if (ind1 .eq. -1) then
                    imin = -NGHOST; imax = -1
                  else if (ind1 .eq. 1) then
                    imin = this_meshblock%ptr%sx; imax = this_meshblock%ptr%sx + NGHOST - 1
                  end if
                  if (ind2 .eq. 0) then
                    jmin = 0; jmax = this_meshblock%ptr%sy - 1
                  else if (ind2 .eq. -1) then
                    jmin = -NGHOST; jmax = -1
                  else if (ind2 .eq. 1) then
                    jmin = this_meshblock%ptr%sy; jmax = this_meshblock%ptr%sy + NGHOST - 1
                  end if
                  if (ind3 .eq. 0) then
                    kmin = 0; kmax = this_meshblock%ptr%sz - 1
                  else if (ind3 .eq. -1) then
                    kmin = -NGHOST; kmax = -1
                  else if (ind3 .eq. 1) then
                    kmin = this_meshblock%ptr%sz; kmax = this_meshblock%ptr%sz + NGHOST - 1
                  end if
                  #ifndef threeD
                    kmin = 0; kmax = 0
                  #endif
                end if

                ! copy `recv_fld` to ghost cells and active cells
                send_cnt = 1
                do i = imin, imax
                  do j = jmin, jmax
                    do k = kmin, kmax
                      if (.not. fill_ghosts) then
                        ! add to existing values
                        jx(i, j, k) = jx(i, j, k) + recv_fld(send_cnt + 0)
                        jy(i, j, k) = jy(i, j, k) + recv_fld(send_cnt + 1)
                        jz(i, j, k) = jz(i, j, k) + recv_fld(send_cnt + 2)
                      else
                        ! overwrite the existing values
                        jx(i, j, k) = recv_fld(send_cnt + 0)
                        jy(i, j, k) = recv_fld(send_cnt + 1)
                        jz(i, j, k) = recv_fld(send_cnt + 2)
                      end if
                      send_cnt = send_cnt + 3
                    end do
                  end do
                end do
              end if ! if receive is ready
            end if ! if receive hasn't been done yet
          end do ! ind3
        end do ! ind2
      end do ! ind1
    end do ! global loop
    call printDiag((mpi_rank .eq. 0), "exchangeCurrents()", .true.)
  end subroutine exchangeCurrents
end module m_exchangecurrents
