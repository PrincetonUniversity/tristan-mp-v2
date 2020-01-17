#include "../defs.F90"

module m_qedphysics
#ifdef QED

  use m_aux
  use m_globalnamespace
  use m_bwpairproduction
  implicit none

  !--- PRIVATE variables/functions -------------------------------!
  !...............................................................!
contains
  subroutine QEDstep(timestep)
    implicit none
    integer, intent(in) :: timestep
    if (modulo(timestep, BW_interval) .eq. 0) then
      call bwPairProduction()
    end if
    call printDiag((mpi_rank .eq. 0), "QEDstep()", .true.)
  end subroutine QEDstep

#endif
end module m_qedphysics
