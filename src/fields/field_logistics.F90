module m_fieldlogistics
  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_errors
  use m_domain
  use m_fields

contains

  subroutine initializeFields()
    implicit none
    call deallocateFields()
    call reallocateFields(this_meshblock % ptr)
    call reallocateFieldBuffers(this_meshblock % ptr)
    call printDiag("initializeFields()", 1)
  end subroutine initializeFields

  subroutine reallocateFields(meshblock)
    implicit none
    type(mesh), intent(in) :: meshblock
    integer :: i1, i2, j1, j2, k1, k2
    i1 = -NGHOST; i2 = meshblock % sx - 1 + NGHOST
    j1 = -NGHOST; j2 = meshblock % sy - 1 + NGHOST
    k1 = -NGHOST; k2 = meshblock % sz - 1 + NGHOST

#ifdef oneD
    j1 = 0; j2 = 0
    k1 = 0; k2 = 0
#elif twoD
    k1 = 0; k2 = 0
#endif

    allocate (ex(i1:i2, j1:j2, k1:k2))
    allocate (ey(i1:i2, j1:j2, k1:k2))
    allocate (ez(i1:i2, j1:j2, k1:k2))
    allocate (bx(i1:i2, j1:j2, k1:k2))
    allocate (by(i1:i2, j1:j2, k1:k2))
    allocate (bz(i1:i2, j1:j2, k1:k2))
    allocate (jx(i1:i2, j1:j2, k1:k2))
    allocate (jy(i1:i2, j1:j2, k1:k2))
    allocate (jz(i1:i2, j1:j2, k1:k2))
    allocate (jx_buff(i1:i2, j1:j2, k1:k2))
    allocate (jy_buff(i1:i2, j1:j2, k1:k2))
    allocate (jz_buff(i1:i2, j1:j2, k1:k2))
    allocate (lg_arr(i1:i2, j1:j2, k1:k2))

    allocate (sm_arr(0:meshblock % sx - 1, 0:meshblock % sy - 1, 0:meshblock % sz - 1))
  end subroutine reallocateFields

  subroutine reallocateFieldBuffers(meshblock)
    implicit none
    type(mesh), intent(in) :: meshblock
    integer :: max_n_fld = 10

    ! `max_n_fld` = max # of fields sent/received in each direction
#ifdef oneD
    sendrecv_offsetsz = NGHOST * max_n_fld
    ! 2 (~5) directions to send/recv in 1D
    sendrecv_buffsz = sendrecv_offsetsz * 5
#elif twoD
    sendrecv_offsetsz = MAX0(meshblock % sx, meshblock % sy, meshblock % sz) * NGHOST * max_n_fld
    ! 8 (~10) directions to send/recv in 2D
    sendrecv_buffsz = sendrecv_offsetsz * 10
#elif threeD
    sendrecv_offsetsz = MAX0(meshblock % sx, meshblock % sy, meshblock % sz)**2 * NGHOST * max_n_fld
    ! 26 (~30) directions to send/recv in 3D
    sendrecv_buffsz = sendrecv_offsetsz * 30
#endif

    allocate (send_fld(sendrecv_buffsz))
    allocate (recv_fld(sendrecv_offsetsz))
#ifndef MPINONBLOCK
    allocate (send_EB(sendrecv_offsetsz))
#endif
  end subroutine reallocateFieldBuffers

  subroutine backupEBfields()
    implicit none
    integer :: i1, i2, j1, j2, k1, k2
    i1 = 0; i2 = this_meshblock % ptr % sx - 1
    j1 = 0; j2 = this_meshblock % ptr % sy - 1
    k1 = 0; k2 = this_meshblock % ptr % sz - 1

    call deallocateFieldBackups()

    allocate (ex_back(i1:i2, j1:j2, k1:k2))
    allocate (ey_back(i1:i2, j1:j2, k1:k2))
    allocate (ez_back(i1:i2, j1:j2, k1:k2))
    allocate (bx_back(i1:i2, j1:j2, k1:k2))
    allocate (by_back(i1:i2, j1:j2, k1:k2))
    allocate (bz_back(i1:i2, j1:j2, k1:k2))

    ex_back(i1:i2, j1:j2, k1:k2) = ex(i1:i2, j1:j2, k1:k2)
    ey_back(i1:i2, j1:j2, k1:k2) = ey(i1:i2, j1:j2, k1:k2)
    ez_back(i1:i2, j1:j2, k1:k2) = ez(i1:i2, j1:j2, k1:k2)
    bx_back(i1:i2, j1:j2, k1:k2) = bx(i1:i2, j1:j2, k1:k2)
    by_back(i1:i2, j1:j2, k1:k2) = by(i1:i2, j1:j2, k1:k2)
    bz_back(i1:i2, j1:j2, k1:k2) = bz(i1:i2, j1:j2, k1:k2)
  end subroutine backupEBfields

  subroutine restoreFieldsFromBackups(i1_from, i2_from, j1_from, j2_from, k1_from, k2_from, &
                                      i1_to, i2_to, j1_to, j2_to, k1_to, k2_to)
    implicit none
    integer, intent(in) :: i1_from, i2_from, j1_from, j2_from, k1_from, k2_from
    integer, intent(in) :: i1_to, i2_to, j1_to, j2_to, k1_to, k2_to
#ifdef DEBUG
    if ((i2_from - i1_from .ne. i2_to - i1_to) .or. &
        (j2_from - j1_from .ne. j2_to - j1_to) .or. &
        (k2_from - k1_from .ne. k2_to - k1_to)) then
      call throwError('ERROR: wrong dimensions in `restoreFieldsFromBackups()`.')
    end if
#endif
    ex(i1_to:i2_to, j1_to:j2_to, k1_to:k2_to) = ex_back(i1_from:i2_from, j1_from:j2_from, k1_from:k2_from)
    ey(i1_to:i2_to, j1_to:j2_to, k1_to:k2_to) = ey_back(i1_from:i2_from, j1_from:j2_from, k1_from:k2_from)
    ez(i1_to:i2_to, j1_to:j2_to, k1_to:k2_to) = ez_back(i1_from:i2_from, j1_from:j2_from, k1_from:k2_from)
    bx(i1_to:i2_to, j1_to:j2_to, k1_to:k2_to) = bx_back(i1_from:i2_from, j1_from:j2_from, k1_from:k2_from)
    by(i1_to:i2_to, j1_to:j2_to, k1_to:k2_to) = by_back(i1_from:i2_from, j1_from:j2_from, k1_from:k2_from)
    bz(i1_to:i2_to, j1_to:j2_to, k1_to:k2_to) = bz_back(i1_from:i2_from, j1_from:j2_from, k1_from:k2_from)
    jx(:, :, :) = 0.0
    jy(:, :, :) = 0.0
    jz(:, :, :) = 0.0
  end subroutine restoreFieldsFromBackups

  subroutine deallocateFields()
    implicit none
    ! fields
    if (allocated(ex)) deallocate (ex)
    if (allocated(ey)) deallocate (ey)
    if (allocated(ez)) deallocate (ez)
    if (allocated(bx)) deallocate (bx)
    if (allocated(by)) deallocate (by)
    if (allocated(bz)) deallocate (bz)
    if (allocated(jx)) deallocate (jx)
    if (allocated(jy)) deallocate (jy)
    if (allocated(jz)) deallocate (jz)
    if (allocated(jx_buff)) deallocate (jx_buff)
    if (allocated(jy_buff)) deallocate (jy_buff)
    if (allocated(jz_buff)) deallocate (jz_buff)
    if (allocated(lg_arr)) deallocate (lg_arr)

    if (allocated(sm_arr)) deallocate (sm_arr)

    ! send/recv buffers
    if (allocated(send_fld)) deallocate (send_fld)
    if (allocated(recv_fld)) deallocate (recv_fld)
#ifndef MPINONBLOCK
    if (allocated(send_EB)) deallocate (send_EB)
#endif
  end subroutine deallocateFields

  subroutine deallocateFieldBackups()
    implicit none
    ! backup fields
    if (allocated(ex_back)) deallocate (ex_back)
    if (allocated(ey_back)) deallocate (ey_back)
    if (allocated(ez_back)) deallocate (ez_back)
    if (allocated(bx_back)) deallocate (bx_back)
    if (allocated(by_back)) deallocate (by_back)
    if (allocated(bz_back)) deallocate (bz_back)
  end subroutine deallocateFieldBackups

end module m_fieldlogistics
