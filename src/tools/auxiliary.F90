#include "../defs.F90"

module m_aux
  use m_globalnamespace
  implicit none
  real(dprec)                  :: dseed
  ! integer, dimension(0:15)        :: state
  ! integer                         :: rand_ind

  abstract interface
    function spatialDistribution(x_glob, y_glob, z_glob,&
                               & dummy1, dummy2, dummy3)
      real :: spatialDistribution
      real, intent(in), optional  :: x_glob, y_glob, z_glob
      real, intent(in), optional  :: dummy1, dummy2, dummy3
    end function spatialDistribution
  end interface

  interface STR
    module procedure intToStr
    module procedure realToStr
  end interface STR

  !--- PRIVATE functions -----------------------------------------!
  private :: intToStr, realToStr
  !...............................................................!
contains
  subroutine printDiag(bool, msg, prepend)
    implicit none
    character(len=*), intent(in)  :: msg
    logical, intent(in)           :: bool
    logical, optional, intent(in) :: prepend
    character(len=STR_MAX)        :: dummy
    integer                       :: sz, i, ierr
    #ifdef DEBUG
      if (bool) then
        sz = len(trim(msg))
        if (present(prepend)) then
          if (prepend) then
            dummy = '...'
            sz = sz + 3
          end if
        else
          dummy = ''
        end if
        dummy = trim(dummy) // trim(msg)
        do i = 1, (38 - sz)
          dummy = trim(dummy) // '.'
        end do
        dummy = trim(dummy) // '[OK]'
        print *, trim(dummy)
      end if
    #endif
  end subroutine printDiag

  subroutine printReport(bool, msg, prepend)
    implicit none
    character(len=*), intent(in)  :: msg
    logical, intent(in)           :: bool
    logical, optional, intent(in) :: prepend
    character(len=STR_MAX)        :: dummy
    integer                       :: sz, i, ierr
    if (bool) then
      sz = len(trim(msg))
      if (present(prepend)) then
        if (prepend) then
          dummy = '...'
          sz = sz + 3
        end if
      else
        dummy = ''
      end if
      dummy = trim(dummy) // trim(msg)
      do i = 1, (38 - sz)
        dummy = trim(dummy) // '.'
      end do
      dummy = trim(dummy) // '[OK]'
      print *, trim(dummy)
    end if
  end subroutine printReport

  subroutine printTimeHeader(tstep)
    implicit none
    integer, intent(in)    :: tstep
    character(len=STR_MAX) :: dummy
    integer                :: sz, i

    ! printing divider
    do i = 68, 70
      dummy(i : i) = ' '
    end do
    do i = 1, 67
      dummy(i : i) = '-'
    end do
    print *, dummy(1:70)

    ! printing timestep
    sz = len(trim("Timestep: " // STR(tstep)))
    do i = 1, 67
      dummy(i : i) = '.'
    end do
    dummy(1 : sz) = trim("Timestep: " // STR(tstep))
    dummy(62:67) = '[DONE]'
    print *, dummy(1:70)

    ! printing header
    do i = 1, 70
      dummy(i : i) = ' '
    end do
    dummy(1:67) = '[ROUTINE]          [TIME, ms]      [MIN/MAX, %]       [FRACTION, %]'
    print *, dummy(1:70)
  end subroutine printTimeHeader

  subroutine printTimeFooter()
    implicit none
    character(len=STR_MAX) :: dummy
    integer                :: i

    do i = 68, 70
      dummy(i : i) = ' '
    end do
    do i = 1, 67
      dummy(i : i) = '.'
    end do
    print *, dummy(1:70)
  end subroutine printTimeFooter

  subroutine printTime(dt_arr, msg, fullstep)
    implicit none
    character(len=*), intent(in)          :: msg
    character(len=STR_MAX)                :: dummy, dummy1
    real(kind=8), intent(in)              :: dt_arr(:)
    real, optional, intent(in)            :: fullstep
    real                                  :: dt_mean, dt_max, dt_min
    integer                               :: pcent_max, pcent_min, sz, sz1, i
    dt_mean = SUM(dt_arr) * 1000 / mpi_size
    dt_max = MAXVAL(dt_arr) * 1000
    dt_min = MINVAL(dt_arr) * 1000
    pcent_max = dt_max * 100 / dt_mean
    pcent_min = dt_min * 100 / dt_mean
    if (present(fullstep)) then
      if (dt_mean / fullstep .lt. 1e-4) then
        dt_mean = 0; pcent_max = 0; pcent_min = 0
      end if
    end if

    do i = 1, 70
      dummy(i : i) = ' '
    end do

    sz = len(msg)
    dummy(1 : sz) = msg

    dummy1 = trim(STR(dt_mean))
    sz = len(trim(dummy1))
    dummy(20 : 20 + sz - 1) = trim(dummy1)
    dummy(36 : 36) = '-'

    dummy1 = trim(STR(pcent_min))
    sz = len(trim(dummy1))
    dummy(37 : 37 + sz - 1) = trim(dummy1)

    dummy(44 : 44) = '+'
    dummy1 = trim(STR(pcent_max))
    sz = len(trim(dummy1))
    dummy(45 : 45 + sz - 1) = trim(dummy1)
    if (present(fullstep)) then
      dummy1 = trim(STR(dt_mean * 100 / fullstep))
      sz1 = len(trim(dummy1))
      dummy(55 : 55 + sz1 - 1) = trim(dummy1)
    end if

    print *, dummy(1:70)
  end subroutine printTime

  subroutine printNpart(npart_arr, msg)
    implicit none
    character(len=*), intent(in)          :: msg
    character(len=STR_MAX)                :: dummy, dummy1
    integer, intent(in)                   :: npart_arr(:)
    integer                               :: npart_mean, npart_max, npart_min
    integer                               :: pcent_max, pcent_min, sz, sz1, i
    npart_mean = SUM(npart_arr) / mpi_size
    npart_max = MAXVAL(npart_arr)
    npart_min = MINVAL(npart_arr)
    if (npart_mean .ne. 0) then
      pcent_max = npart_max * 100 / npart_mean
    else
      pcent_max = 0
    end if
    if (npart_mean .ne. 0) then
      pcent_min = npart_min * 100 / npart_mean
    else
      pcent_min = 0
    end if

    do i = 1, 70
      dummy(i : i) = ' '
    end do

    sz = len(msg)
    dummy(1 : sz) = msg

    dummy1 = trim(STR(npart_mean))
    sz = len(trim(dummy1))
    dummy(20 : 20 + sz - 1) = trim(dummy1)
    dummy(36 : 36) = '-'

    dummy1 = trim(STR(pcent_min))
    sz = len(trim(dummy1))
    dummy(37 : 37 + sz - 1) = trim(dummy1)

    dummy(44 : 44) = '+'
    dummy1 = trim(STR(pcent_max))
    sz = len(trim(dummy1))
    dummy(45 : 45 + sz - 1) = trim(dummy1)

    print *, dummy(1:70)
  end subroutine printNpart

  function intToStr(my_int) result(string)
    implicit none
    integer, intent(in)       :: my_int
    character(:), allocatable :: string
    character(len=STR_MAX)    :: temp
    write(temp, '(i0)') my_int
    string = trim(temp)
  end function intToStr

  function realToStr(my_real) result(string)
    implicit none
    real, intent(in)          :: my_real
    character(:), allocatable :: string
    character(len=STR_MAX)    :: temp
    write(temp, '(G0.2)') my_real
    string = trim(temp)
  end function realToStr

  function STRtoINT(my_str) result(my_int)
    implicit none
    character(len=*), intent(in)  :: my_str
    integer                       :: my_int
    read (my_str, *) my_int
  end function STRtoINT

  !--- Taken from Zeltron -------------------------------------------------!
  ! Reference: http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
  !........................................................................!

  real(dprec) function randomNum(DSEED)
    ! implicit none
    ! integer, intent(in) :: dseed
    ! call random_number(random)
  	implicit none
  	real(dprec)    :: DSEED
  	integer        :: I
  	real(dprec)    :: S2P31, S2P31M, SEED
  	DATA              S2P31M/2147483647.D0/,S2P31/2147483648.D0/
  	SEED = DSEED
    SEED = DMOD(16807.D0*SEED,S2P31M)
    randomNum = SEED / S2P31
  	DSEED = SEED
  	return
  end function randomNum

  real function random(DSEED)
  	implicit none
  	real(dprec)    :: DSEED
    real           :: rnd
    rnd = 1.0
    do while(rnd .eq. 1.0)
      rnd = randomNum(DSEED)
    end do
    random = rnd
    return
  end function random

  integer function randomInt(DSEED, amin, amax)
    implicit none
    real(dprec)         :: DSEED
    integer, intent(in) :: amin, amax
    randomInt = amin + INT((amax - amin) * random(dseed))
    return
  end function randomInt

  real function poisson(num)
    implicit none
    real, intent(in) :: num
    real(kind=8)     :: Lps, pps
    real             :: kps, ups
    Lps = EXP(-REAL(num, 8))
    kps = 0
    pps = 1
    do while (pps .ge. Lps)
      kps = kps + 1
      ups = random(dseed)
      pps = pps * ups
    end do
    poisson = kps - 1
    return
  end function poisson

  subroutine initializeRandomSeed(rank)
    implicit none
    integer, intent(in) :: rank
    dseed = 123457.D0
    dseed = dseed + rank
    ! integer :: i, n, clock
    ! integer, dimension(:), allocatable :: seed
    !
    ! call random_seed(size = n)
    ! allocate(seed(n))
    !
    ! call system_clock(COUNT = clock)
    !
    ! seed = clock + rank * (/ (i - 1, i = 1, n) /)
    ! call random_seed(PUT = seed)
    ! deallocate(seed)
  end subroutine initializeRandomSeed

  recursive function factorial(n) result(fact)
    implicit none
    integer             :: fact
    integer, intent(in) :: n
    if (n .eq. 0) then
      fact = 1
    else
      fact = n * factorial(n - 1)
    end if
  end function factorial

end module m_aux
