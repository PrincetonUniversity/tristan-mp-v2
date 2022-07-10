module m_aux
  use m_globalnamespace
  implicit none
  real(dprec) :: dseed

  abstract interface
    function spatialDistribution(x_glob, y_glob, z_glob, &
                                 dummy1, dummy2, dummy3)
      real :: spatialDistribution
      real, intent(in), optional :: x_glob, y_glob, z_glob
      real, intent(in), optional :: dummy1, dummy2, dummy3
    end function spatialDistribution
  end interface

  interface STR
    module procedure intToStr
    module procedure realToStr
  end interface STR

  type :: generic_var
    integer :: value_int
    real :: value_real
    logical :: value_bool
  end type generic_var

  type generic_string
    character(len=STR_MAX) :: str
  end type generic_string

  type :: simulation_params
    integer :: count
    integer, allocatable :: param_type(:) ! 1 = int, 2 = float, 3 = bool
    type(generic_string), allocatable :: param_group(:)
    type(generic_string), allocatable :: param_name(:)
    type(generic_var), allocatable :: param_value(:)
  end type simulation_params

  abstract interface
    function getFMT(value, w) result(FMT)
      implicit none
      real, intent(in) :: value
      character(len=STR_MAX) :: FMT
      integer, intent(in), optional :: w
    end function getFMT
  end interface

  type :: warning
    character(len=STR_MAX) :: description
    integer :: counter = 0
  end type warning

  type(warning) :: warnings(100)
  type(simulation_params) :: sim_params

  !--- PRIVATE functions -----------------------------------------!
  private :: intToStr, realToStr
  !...............................................................!
contains
  subroutine initializeSimulationParameters()
    implicit none
    sim_params % count = 0
    allocate (sim_params % param_type(1000))
    allocate (sim_params % param_group(1000))
    allocate (sim_params % param_name(1000))
    allocate (sim_params % param_value(1000))
  end subroutine initializeSimulationParameters

  subroutine initializeWarnings()
    implicit none
    warnings(1) % description = "Synchrotron cooling is too strong"
    warnings(2) % description = "QED probability too large"
    warnings(3) % description = "IC cooling is too strong"
  end subroutine initializeWarnings

  subroutine printDiag(msg, level)
    implicit none
    character(len=*), intent(in) :: msg
    integer, optional, intent(in) :: level
    character(len=STR_MAX) :: dummy
    integer :: sz, i
    if (mpi_rank .eq. 0) then
      open (UNIT_diag, file=diag_file_name, status="old", position="append", form="formatted")
      sz = len(trim(msg))
      dummy = ''
      if (present(level)) then
        sz = sz + level * 3
        do i = 1, level * 3
          dummy(i:i) = '.'
        end do
      end if
      dummy = trim(dummy)//trim(msg)
      write (UNIT_diag, *) trim(dummy)
      close (UNIT_diag)
    end if
  end subroutine printDiag

  subroutine addWarning(id)
    implicit none
    integer, intent(in) :: id
    warnings(id) % counter = warnings(id) % counter + 1
  end subroutine addWarning

  subroutine printWarnings(timestep)
    implicit none
    integer, intent(in) :: timestep
    integer :: ierr, root_rnk = 0, w
    integer :: warnings_global(100)

    do w = 1, 100
      call MPI_REDUCE(warnings(w) % counter, warnings_global(w), 1, MPI_INTEGER, &
                      MPI_SUM, root_rnk, MPI_COMM_WORLD, ierr)
    end do

    if (mpi_rank .eq. root_rnk) then
      open (UNIT_warn, file=warn_file_name, status="old", position="append", form="formatted")
      write (UNIT_warn, *) '=================================================='
      write (UNIT_warn, *) 'Timestep = '//trim(STR(timestep))
      do w = 1, 100
        if (warnings_global(w) .gt. 0) then
          write (UNIT_warn, *) trim(warnings(w) % description)//' -> called '//trim(STR(warnings_global(w)))//' times'
        end if
      end do
      write (UNIT_warn, *) '..................................................'
      close (UNIT_warn)
    end if

    do w = 1, 100
      warnings(w) % counter = 0
    end do
  end subroutine printWarnings

    function getFMTForReal(value, w) result(FMT)
    implicit none
    real, intent(in)              :: value
    character(len=STR_MAX)        :: FMT
    integer, intent(in), optional :: w
    integer                       :: w_
    character(len=10)             :: dummy
    if (.not. present(w)) then
      w_ = 10
    else
      w_ = w
    end if
    write(dummy, '(I10)') w_

    if ((abs(value) .ge. 100000) .or.&
      & ((abs(value) .lt. 1e-2) .and.&
        & (abs(value) .ne. 0.0))) then
      FMT = 'ES' // trim(dummy) // '.3'
    else
      FMT = 'F' // trim(dummy) // '.3'
    end if
  end function getFMTForReal
  
  function getFMTForRealScientific(w) result(FMT)
    implicit none
    character(len=STR_MAX) :: FMT
    integer, intent(in), optional :: w
    integer :: w_
    character(len=10) :: dummy
    if (.not. present(w)) then
      w_ = 10
    else
      w_ = w
    end if
    write (dummy, '(I10)') w_
    FMT = 'ES'//trim(dummy)//'.3'
  end function getFMTForRealScientific

  subroutine printTimeHeader(tstep)
    implicit none
    integer, intent(in) :: tstep
    integer :: i

    do i = 1, 71
      if (i .eq. 1) then
        write (*, '(A2)', advance='no') ' -'
      else if (i .eq. 71) then
        print "(A)", '-'
      else
        write (*, '(A)', advance='no') '-'
      end if
    end do

    write (*, "(1X,A10,I9)", advance='no') "Timestep: ", tstep
    do i = 1, 48
      write (*, "(A)", advance='no') '.'
    end do
    print "(A4)", "[OK]"
    print "(1X,A71)", "[ROUTINE]          [TIME, ms]      [MIN  /  MAX, ms]      [FRACTION, %]"
  end subroutine printTimeHeader

  subroutine printTimeFooter()
    implicit none
    integer :: i

    do i = 1, 71
      if (i .eq. 1) then
        write (*, '(A2)', advance='no') ' .'
      else if (i .eq. 71) then
        print "(A)", '.'
      else
        write (*, '(A)', advance='no') '.'
      end if
    end do
  end subroutine printTimeFooter

  subroutine printTime(dt_arr, msg, fullstep)
    implicit none
    character(len=*), intent(in) :: msg
    character(len=STR_MAX) :: FMT
    character(len=15) :: msg_str
    real(kind=8), intent(in) :: dt_arr(:)
    real, optional, intent(in) :: fullstep
    real :: dt_mean, dt_max, dt_min
    dt_mean = REAL(SUM(dt_arr) * 1000 / mpi_size)
    dt_max = REAL(MAXVAL(dt_arr) * 1000)
    dt_min = REAL(MINVAL(dt_arr) * 1000)
    if (present(fullstep)) then
      if (dt_mean / fullstep .lt. 1e-4) then
        dt_mean = 0; dt_min = 0; dt_max = 0
      end if
    end if

    write (msg_str, '(A15)') msg
    
    if (.not. present(fullstep)) then
        FMT = "(1X,A15,ES14.6,ES12.4,ES11.3)"
        !FMT = "(1X,A15" // &
                !trim(getFMTForReal(dt_mean, 14)) // "," // &
                !trim(getFMTForReal(dt_min, 12)) // "," // &
                !trim(getFMTForReal(dt_max, 11)) // &
               !")"
        print FMT, adjustl(msg_str), dt_mean, dt_min, dt_max
    else
        FMT = "(3X,A13,ES14.6,ES12.4,ES11.3,ES19.2)"
        !FMT = "(3X,A13" // &
                !trim(getFMTForReal(dt_mean, 14)) // "," // &
                !trim(getFMTForReal(dt_min, 12)) // "," // &
                !trim(getFMTForReal(dt_max, 11)) // "," // &
                !trim(getFMTForReal(dt_mean * 100 / fullstep, 19)) // &
               !")"
        print FMT, adjustl(msg_str), dt_mean, dt_min, dt_max, dt_mean * 100 / fullstep
    end if
  end subroutine printTime

  subroutine printNpartHeader()
    implicit none
    print "(1X,A71)", "[NPART per S]       [AVERAGE]      [MIN/MAX per CPU]            [TOTAL]"
  end subroutine printNpartHeader

  subroutine printNpart(npart_arr, msg)
    implicit none
    character(len=*), intent(in) :: msg
    integer(kind=8), intent(in) :: npart_arr(:)
    real :: npart_mean, npart_max, npart_min, npart_sum
    character(len=14) :: msg_str
    character(len=12) :: min_str, max_str
    npart_sum = REAL(SUM(npart_arr))
    npart_mean = npart_sum / mpi_size
    npart_max = REAL(MAXVAL(npart_arr))
    npart_min = REAL(MINVAL(npart_arr))
    write (msg_str, '(A14)') msg
    write (min_str, '(ES12.4)') npart_min
    write (max_str, '(ES12.4)') npart_max
    print "(3X,A14,ES13.6,A12,A1,A12,ES17.10)", adjustl(msg_str), npart_mean, min_str, "/", adjustl(max_str), npart_sum
  end subroutine printNpart

  function intToStr(my_int) result(string)
    implicit none
    integer, intent(in) :: my_int
    character(:), allocatable :: string
    character(len=STR_MAX) :: temp
    write (temp, '(i0)') my_int
    string = trim(temp)
  end function intToStr

  function realToStr(my_real) result(string)
    implicit none
    real, intent(in) :: my_real
    character(:), allocatable :: string
    character(len=STR_MAX) :: temp
    if ((abs(my_real) .ge. 1000.0) .or. ((abs(my_real) .lt. 1e-2) .and. (abs(my_real) .gt. 0.0))) then
      write (temp, '(ES10.2)') my_real
    else
      write (temp, '(F10.2)') my_real
    end if
    string = trim(temp)
  end function realToStr

  function STRtoINT(my_str) result(my_int)
    implicit none
    character(len=*), intent(in) :: my_str
    integer :: my_int
    read (my_str, *) my_int
  end function STRtoINT

  logical function arraysAreEqual(array1, array2)
    integer, dimension(:), intent(in) :: array1, array2
    integer :: i
    arraysAreEqual = (size(array1) .eq. size(array2))
    if (arraysAreEqual) then
      do i = 1, size(array1)
        arraysAreEqual = (array1(i) .eq. array2(i))
        if (.not. arraysAreEqual) exit
      end do
    end if
  end function arraysAreEqual

  real(dprec) function randomNum(DSEED)
    implicit none
    real(dprec) :: DSEED
    real(dprec) :: S2P31, S2P31M, SEED
    DATA S2P31M/2147483647.D0/, S2P31/2147483648.D0/
    SEED = DSEED
    SEED = DMOD(16807.D0 * SEED, S2P31M)
    randomNum = SEED / S2P31
    DSEED = SEED
    return
  end function randomNum

  real function random(DSEED)
    implicit none
    real(dprec) :: DSEED
    real :: rnd
    rnd = 1.0
    do while (rnd .ge. 1.0)
      rnd = REAL(randomNum(DSEED))
    end do
    random = rnd
    return
  end function random

  integer function randomInt(DSEED, amin, amax)
    implicit none
    real(dprec) :: DSEED
    integer, intent(in) :: amin, amax
    randomInt = amin + INT((amax - amin) * random(dseed))
    return
  end function randomInt

  real function poisson(num)
    implicit none
    real, intent(in) :: num
    real(kind=8) :: Lps, pps
    real :: kps, ups
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

    call printDiag("initializeRandomSeed()", 1)
  end subroutine initializeRandomSeed

  subroutine log_normal(n_bins, lognorm)
    integer, intent(in) :: n_bins
    real, allocatable, intent(inout) :: lognorm(:)
    real :: x, y, z, sum
    integer :: i

    allocate (lognorm(n_bins))
    sum = 0.0
    do i = 1, n_bins
      x = random(dseed)
      y = random(dseed)
      z = sqrt(-2.0 * log(x)) * cos(2.0 * M_PI * y) ! now z has standard normal distribution
      z = exp(0.0 + 1.0 * z) ! now z has lognormal distribution with certain sigma=1 and mu=0
      lognorm(i) = z
      sum = sum + z
    end do

    ! this allows having lognorm(max) != 1/
    !   in this case bins are not fixed in upper limit
    lognorm(1) = lognorm(1) / (n_bins + 1.)
    do i = 2, n_bins
      lognorm(i) = lognorm(i) / (n_bins + 1.) + lognorm(i - 1)
    end do
    ! /this allows having lognorm(max) != 1
  end subroutine log_normal

  subroutine lin_normal(n_bins, linnorm)
    integer, intent(in) :: n_bins
    real, allocatable, intent(inout) :: linnorm(:)
    real :: x, sum
    integer :: i

    allocate (linnorm(n_bins))
    sum = 0.0
    do i = 1, n_bins
      x = random(dseed)
      linnorm(i) = x
      sum = sum + x
    end do

    ! this allows having linnorm(max) != 1/
    !   in this case bins are not fixed in upper limit
    linnorm(1) = linnorm(1) / (n_bins + 1.0)
    do i = 2, n_bins
      linnorm(i) = linnorm(i) / (n_bins + 1.0) + linnorm(i - 1)
    end do

    do i = 1, n_bins
      linnorm(i) = 2.0 * linnorm(i)
    end do
    ! /this allows having linnorm(max) != 1
  end subroutine lin_normal

  recursive function factorial(n) result(fact)
    implicit none
    integer :: fact
    integer, intent(in) :: n
    if (n .eq. 0) then
      fact = 1
    else
      fact = n * factorial(n - 1)
    end if
  end function factorial

  subroutine rotateRandomlyIn3D(rx, ry, rz, rnd1, rnd2, rnd3)
    implicit none
    real, intent(inout) :: rx, ry, rz
    real, optional, intent(in) :: rnd1, rnd2, rnd3
    real :: rnd1_, rnd2_, rnd3_, dummy1, dummy2
    real :: rx_, ry_, rz_
    real :: ux, uy, uz, cos_phi, one_m_cos_phi, sin_phi
    ! generate optional arguments
    ! ... each random number is uniform in [0, 1)
    if (.not. present(rnd1)) then
      rnd1_ = random(dseed)
    else
      rnd1_ = rnd1
    end if
    if (.not. present(rnd2)) then
      rnd2_ = random(dseed)
    else
      rnd2_ = rnd2
    end if
    if (.not. present(rnd3)) then
      rnd3_ = random(dseed)
    else
      rnd3_ = rnd3
    end if
    ! generate a random direction in 3d
    dummy1 = 2.0 * rnd1_ - 1.0
    dummy2 = 2.0 * M_PI * rnd2_
    ux = sqrt(1.0 - dummy1**2) * cos(dummy2)
    uy = sqrt(1.0 - dummy1**2) * sin(dummy2)
    uz = dummy1
    ! generate a random angle of rotation
    cos_phi = cos(2.0 * M_PI * rnd3_)
    sin_phi = sin(2.0 * M_PI * rnd3_)
    ! copy old values
    rx_ = rx; ry_ = ry; rz_ = rz

    one_m_cos_phi = (1.0 - cos_phi)

    rx = (one_m_cos_phi * ux**2 + cos_phi) * rx_ + &
         (one_m_cos_phi * ux * uy - sin_phi * uz) * ry_ + &
         (one_m_cos_phi * ux * uz + sin_phi * uy) * rz_

    ry = (one_m_cos_phi * ux * uy + sin_phi * uz) * rx_ + &
         (one_m_cos_phi * uy**2 + cos_phi) * ry_ + &
         (one_m_cos_phi * uy * uz - sin_phi * ux) * rz_

    rz = (one_m_cos_phi * ux * uz - sin_phi * uy) * rx_ + &
         (one_m_cos_phi * uy * uz + sin_phi * ux) * ry_ + &
         (one_m_cos_phi * uz**2 + cos_phi) * rz_
  end subroutine rotateRandomlyIn3D

  subroutine rotateIn3D(rx, ry, rz, ax, ay, az, ang)
    implicit none
    real, intent(inout) :: rx, ry, rz
    real, optional, intent(in) :: ax, ay, az
    real :: ang
    real :: rx_, ry_, rz_
    real :: ux, uy, uz, cos_phi, one_m_cos_phi, sin_phi
    ! generate optional arguments
    ux = ax
    uy = ay
    uz = az
    cos_phi = cos(ang)
    sin_phi = sin(ang)
    ! copy old values
    rx_ = rx; ry_ = ry; rz_ = rz

    one_m_cos_phi = (1.0 - cos_phi)

    rx = (one_m_cos_phi * ux**2 + cos_phi) * rx_ + &
         (one_m_cos_phi * ux * uy - sin_phi * uz) * ry_ + &
         (one_m_cos_phi * ux * uz + sin_phi * uy) * rz_

    ry = (one_m_cos_phi * ux * uy + sin_phi * uz) * rx_ + &
         (one_m_cos_phi * uy**2 + cos_phi) * ry_ + &
         (one_m_cos_phi * uy * uz - sin_phi * ux) * rz_

    rz = (one_m_cos_phi * ux * uz - sin_phi * uy) * rx_ + &
         (one_m_cos_phi * uy * uz + sin_phi * ux) * ry_ + &
         (one_m_cos_phi * uz**2 + cos_phi) * rz_
  end subroutine rotateIn3D

end module m_aux
