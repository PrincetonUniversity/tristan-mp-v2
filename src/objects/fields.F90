#include "../defs.F90"

!--- FIELDS ----------------------------------------------------!
! To store all the field related quantities
!   - field indices go
!       from `-NGHOST`
!       to `sx/sy/sz - 1 + NGHOST` inclusively
!...............................................................!

module m_fields
  use m_globalnamespace
  implicit none

  real, allocatable     :: ex(:,:,:), ey(:,:,:), ez(:,:,:),&
                         & bx(:,:,:), by(:,:,:), bz(:,:,:)
  real, allocatable     :: jx(:,:,:), jy(:,:,:), jz(:,:,:)
  real, allocatable     :: recv_fld(:), send_fld(:)
  integer               :: sendrecv_buffsz, sendrecv_offsetsz
  real, allocatable     :: lg_arr(:,:,:)
  real, allocatable     :: sm_arr(:,:,:)
end module m_fields
