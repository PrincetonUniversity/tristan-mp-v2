#include "../defs.F90"

module m_testsite
  use m_globalnamespace
  use m_outputnamespace
  use m_helpers
  use m_aux
  use m_particles
  use m_fields
  use m_domain
  use m_writeslice
  use m_writetot
  use m_writehistory
  use m_restart
  use m_fldsolver
  use m_mover
  use m_currentdeposit
  use m_exchangeparts
  use m_exchangefields
  use m_exchangecurrents
  use m_particlelogistics
  use m_filtering
  use m_userfile
  use m_errors

  implicit none
  !...............................................................!
contains

  subroutine testcode()
    implicit none
    
    ! YOUR TEST CODE GOES HERE ...

    call checkEverything()
  end subroutine testcode

end module m_testsite
