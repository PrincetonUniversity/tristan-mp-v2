#include "../defs.F90"

module m_exchangeparts
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_particlelogistics
  !--- PRIVATE functions -----------------------------------------!
  private :: copyToEnroute, copyFromEnroute,&
           & extractParticlesFromEnroute
  !...............................................................!
contains
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
      enroute_bot%get(:,:,:)%cnt_send = 0
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
                enroute_bot%get(send_x, send_y, send_z)%cnt_send = enroute_bot%get(send_x, send_y, send_z)%cnt_send + 1
                if (enroute_bot%get(send_x, send_y, send_z)%cnt_send .ge.&
                  & enroute_bot%get(send_x, send_y, send_z)%max_send) then
                  call throwError('ERROR: particle send buffer array too small.')
                end if
                cntr = enroute_bot%get(send_x, send_y, send_z)%cnt_send
                call copyToEnroute(s, ti, tj, tk, p, enroute_bot%get(send_x, send_y, send_z)%send_enroute(cntr))
                ! make ghost particle
                pt_proc(p) = -1

                ! shift coordinates to fit the new grid
                #if defined(oneD) || defined(twoD) || defined(threeD)
                  new_xyz = enroute_bot%get(send_x, send_y, send_z)%send_enroute(cntr)%xi
                  temp_xyz = this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sx
                  new_xyz = -(send_x - 1) * (2 + send_x) * (new_xyz * (send_x + 1) - (temp_xyz - 1) * send_x) / 2
                  enroute_bot%get(send_x, send_y, send_z)%send_enroute(cntr)%xi = new_xyz
                #endif

                #if defined(twoD) || defined(threeD)
                  new_xyz = enroute_bot%get(send_x, send_y, send_z)%send_enroute(cntr)%yi
                  temp_xyz = this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sy
                  new_xyz = -(send_y - 1) * (2 + send_y) * (new_xyz * (send_y + 1) - (temp_xyz - 1) * send_y) / 2
                  enroute_bot%get(send_x, send_y, send_z)%send_enroute(cntr)%yi = new_xyz
                #endif

                #if defined(threeD)
                  new_xyz = enroute_bot%get(send_x, send_y, send_z)%send_enroute(cntr)%zi
                  temp_xyz = this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sz
                  new_xyz = -(send_z - 1) * (2 + send_z) * (new_xyz * (send_z + 1) - (temp_xyz - 1) * send_z) / 2
                  enroute_bot%get(send_x, send_y, send_z)%send_enroute(cntr)%zi = new_xyz
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
            call MPI_ISEND(enroute_bot%get(ind1,ind2,ind3)%send_enroute(1:enroute_bot%get(ind1,ind2,ind3)%cnt_send),&
                         & enroute_bot%get(ind1,ind2,ind3)%cnt_send, myMPI_ENROUTE,&
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
                  call MPI_RECV(recv_enroute(1:cnt_recv_enroute), cnt_recv_enroute, myMPI_ENROUTE,&
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
    end do ! loop over species
    call printDiag((mpi_rank .eq. 0), "exchangeParticles()", .true.)
  end subroutine exchangeParticles

  subroutine copyToEnroute(spec_id, ti, tj, tk, prtl_id, enroute)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in)               :: spec_id, prtl_id, ti, tj, tk
    type(prtl_enroute), intent(inout) :: enroute
    enroute%weight = species(spec_id)%prtl_tile(ti, tj, tk)%weight(prtl_id)
    enroute%xi = species(spec_id)%prtl_tile(ti, tj, tk)%xi(prtl_id)
    enroute%yi = species(spec_id)%prtl_tile(ti, tj, tk)%yi(prtl_id)
    enroute%zi = species(spec_id)%prtl_tile(ti, tj, tk)%zi(prtl_id)
    enroute%dx = species(spec_id)%prtl_tile(ti, tj, tk)%dx(prtl_id)
    enroute%dy = species(spec_id)%prtl_tile(ti, tj, tk)%dy(prtl_id)
    enroute%dz = species(spec_id)%prtl_tile(ti, tj, tk)%dz(prtl_id)
    enroute%u = species(spec_id)%prtl_tile(ti, tj, tk)%u(prtl_id)
    enroute%v = species(spec_id)%prtl_tile(ti, tj, tk)%v(prtl_id)
    enroute%w = species(spec_id)%prtl_tile(ti, tj, tk)%w(prtl_id)
    enroute%ind = species(spec_id)%prtl_tile(ti, tj, tk)%ind(prtl_id)
    enroute%proc = species(spec_id)%prtl_tile(ti, tj, tk)%proc(prtl_id)
    #ifdef PRTLPAYLOADS
      enroute%payload1 = species(spec_id)%prtl_tile(ti, tj, tk)%payload1(prtl_id)
      enroute%payload2 = species(spec_id)%prtl_tile(ti, tj, tk)%payload2(prtl_id)
      enroute%payload3 = species(spec_id)%prtl_tile(ti, tj, tk)%payload3(prtl_id)
    #endif
  end subroutine copyToEnroute

  subroutine copyFromEnroute(enroute, spec_id)
    implicit none
    type(prtl_enroute), intent(in)  :: enroute
    integer, intent(in)             :: spec_id
    ! DEP_PRT [particle-dependent]
    call createParticleFromAttributes(spec_id, enroute%xi, enroute%yi, enroute%zi,&
                                             & enroute%dx, enroute%dy, enroute%dz,&
                                             & enroute%u, enroute%v, enroute%w,&
                                             #ifdef PRTLPAYLOADS
                                              & enroute%payload1, enroute%payload2, enroute%payload3,&
                                             #endif
                                             & enroute%ind, enroute%proc, enroute%weight)
  end subroutine copyFromEnroute

  subroutine moveParticleBetweenTiles(s, ti, tj, tk, p)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in) :: s, ti, tj, tk, p
    call createParticleFromAttributes(s, species(s)%prtl_tile(ti, tj, tk)%xi(p),&
                                       & species(s)%prtl_tile(ti, tj, tk)%yi(p),&
                                       & species(s)%prtl_tile(ti, tj, tk)%zi(p),&
                                       & species(s)%prtl_tile(ti, tj, tk)%dx(p),&
                                       & species(s)%prtl_tile(ti, tj, tk)%dy(p),&
                                       & species(s)%prtl_tile(ti, tj, tk)%dz(p),&
                                       & species(s)%prtl_tile(ti, tj, tk)%u(p),&
                                       & species(s)%prtl_tile(ti, tj, tk)%v(p),&
                                       & species(s)%prtl_tile(ti, tj, tk)%w(p),&
                                       #ifdef PRTLPAYLOADS
                                         & species(s)%prtl_tile(ti, tj, tk)%payload1(p),&
                                         & species(s)%prtl_tile(ti, tj, tk)%payload2(p),&
                                         & species(s)%prtl_tile(ti, tj, tk)%payload3(p),&
                                       #endif
                                       & species(s)%prtl_tile(ti, tj, tk)%ind(p),&
                                       & species(s)%prtl_tile(ti, tj, tk)%proc(p),&
                                       & species(s)%prtl_tile(ti, tj, tk)%weight(p))
    ! schedule particle for deletion
    species(s)%prtl_tile(ti, tj, tk)%proc(p) = -1
  end subroutine moveParticleBetweenTiles

  subroutine extractParticlesFromEnroute(cnt, spec_id)
    implicit none
    integer, intent(in)   :: cnt, spec_id
    integer               :: p
    do p = 1, cnt
      call copyFromEnroute(recv_enroute(p), spec_id)
    end do
  end subroutine extractParticlesFromEnroute
end module m_exchangeparts
