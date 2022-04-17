#include "../defs.F90"

module m_errors
  use m_globalnamespace
  implicit none
contains
  subroutine throwError(msg)
    character(len=*), intent(in)  :: msg
    print *, msg
    stop
  end subroutine
end module m_errors
