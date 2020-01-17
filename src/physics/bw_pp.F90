#include "../defs.F90"

module m_bwpairproduction
#ifdef BWPAIRPRODUCTION

  use m_globalnamespace
  use m_aux
  use m_errors
  use m_bincoupling
  use m_particlelogistics
  implicit none

  ! BW parameters
  real    :: BW_tau
  integer :: BW_interval, BW_electron_sp, BW_positron_sp
  integer :: BW_algorithm

  !--- PRIVATE variables/functions -------------------------------!
  !...............................................................!
contains
  subroutine bwPairProduction()
    implicit none

    ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
    !
    !
    ! THIS ROUTINE IS STILL UNDER CONSTRUCTION
    !   if you need access, please contact developers
    !
    !
    ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  end subroutine bwPairProduction

#endif
end module m_bwpairproduction
