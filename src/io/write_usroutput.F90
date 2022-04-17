#include "../defs.F90"

module m_writeusroutput
#ifdef USROUTPUT
  use m_globalnamespace
  use m_outputnamespace, only: usrout_enable, usrout_interval
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_readinput, only: getInput
  use m_helpers
  implicit none

  character(len=STR_MAX), private :: usrout_filename
  logical, private                :: first_time

contains
  subroutine initializeUsrOutput()
    implicit none
    call getInput('output', 'usr_enable', usrout_enable, .false.)
    call getInput('output', 'usr_interval', usrout_interval, 1)
    usrout_filename = trim(output_dir_name) // '/usroutput'
    first_time = .true.
  end subroutine initializeUsrOutput

  subroutine writeUsrOutputTimestep(step)
    implicit none
    integer, intent(in) :: step
    if (first_time) then
      open(UNIT_usrout, file=usrout_filename, status="replace", form="formatted")
      first_time = .false.
    else
      open(UNIT_usrout, file=usrout_filename, status="old", position="append", form="formatted")
    end if
    write(UNIT_usrout, '(A, I5)') 't =', step
    close(UNIT_usrout)
  end subroutine writeUsrOutputTimestep

  subroutine writeUsrOutputArray(name, arr)
    implicit none
    character(len=*), intent(in)    :: name
    real, intent(in), allocatable   :: arr(:)
    integer                         :: i
    open(UNIT_usrout, file=usrout_filename, status="old", position="append", form="formatted")
    write(UNIT_usrout, '(A)') trim(name) // ':'
    do i = 1, size(arr)
      write(UNIT_usrout, "(ES23.16)", advance="no") arr(i)
      write(UNIT_usrout, '(A)', advance="no") ','
    end do
    write(UNIT_usrout, '(A)') ''
    close(UNIT_usrout)
  end subroutine writeUsrOutputArray

  subroutine writeUsrOutputReal(name, value)
    implicit none
    character(len=*), intent(in)  :: name
    real, intent(in)              :: value
    open(UNIT_usrout, file=usrout_filename, status="old", position="append", form="formatted")
    write(UNIT_usrout, '(A)') trim(name) // ':'
    write(UNIT_usrout, "(ES23.16)") value
    close(UNIT_usrout)
  end subroutine writeUsrOutputReal

  subroutine writeUsrOutputEnd()
    implicit none
    open(UNIT_usrout, file=usrout_filename, status="old", position="append", form="formatted")
    write(UNIT_usrout, '(A)') '=============================='
    close(UNIT_usrout)
  end subroutine writeUsrOutputEnd
#endif
end module m_writeusroutput
