#include "../defs.F90"

program tristan
  use m_initialize
  use m_mainloop
  use m_testsite
  use m_finalize
  implicit none
  !----- main code --------------------------

  call initializeAll()
  #ifndef TESTMODE
    call mainloop()
  #else
    call testcode()
  #endif
  call finalizeAll()

  !..... main code ..........................
end program tristan
