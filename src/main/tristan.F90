#include "../defs.F90"

program tristan
  use m_initialize, only: initializeAll
  use m_mainloop, only: mainloop
  use m_finalize, only: finalizeAll
  implicit none
  !----- main code --------------------------

  call initializeAll()
  call mainloop()
  call finalizeAll()

  !..... main code ..........................
end program tristan
