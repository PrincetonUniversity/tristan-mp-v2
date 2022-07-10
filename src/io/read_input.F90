#define END_LINE char(10)

module m_readinput
  use m_globalnamespace
  use m_aux
  use m_errors
  implicit none

  interface strToNum
    module procedure strToInt4
    module procedure strToReal4
    module procedure strToInt8
    module procedure strToReal8
  end interface strToNum

  interface getInput
    module procedure getInt4Input
    module procedure getReal4Input
    module procedure getInt8Input
    module procedure getReal8Input
    module procedure getLogicalInput
  end interface getInput

  !--- PRIVATE functions -----------------------------------------!
  private :: strToInt4, strToInt8, strToReal4, strToReal8, &
             getInt4Input, getInt8Input, getReal4Input, &
             getReal8Input, getLogicalInput, parseInput
  !...............................................................!
contains
  ! read input/output filename/directory
  subroutine readCommandlineArgs()
    implicit none
    integer :: i
    character(len=STR_MAX) :: arg, arg1
    do i = 1, command_argument_count()
      call get_command_argument(i, arg)
      select case (arg)
      case ('-i', '--input')
        call get_command_argument(i + 1, arg1)
        input_file_name = trim(arg1)
      case ('-o', '--output')
        call get_command_argument(i + 1, arg1)
        output_dir_name = trim(arg1)
      case ('-r', '--restart')
        call get_command_argument(i + 1, arg1)
        restart_dir_name = trim(arg1)
      case ('-s', '--slice')
        call get_command_argument(i + 1, arg1)
        slice_dir_name = trim(arg1)
      case ('-R', '--RST')
        call get_command_argument(i + 1, arg1)
        restart_from = trim(arg1)
      end select
    end do
  end subroutine readCommandlineArgs

  character(len=STR_MAX) function parseInput(blockname, varname, found)
    implicit none
    character(len=*), intent(in) :: blockname, varname
    integer :: iostatus, k1, k2
    logical, intent(out) :: found
    character(len=STR_MAX) :: istream, value_str

    k1 = 0; k2 = 0
    ! opening the input file
    open (unit=UNIT_input, file=trim(input_file_name), &
          action='read', IOSTAT=iostatus)
    if (iostatus .ne. 0) then
      call throwError('ERROR while reading input file: `'//trim(input_file_name)//'`')
    end if
    find_blockname: do while (.true.)
      read (UNIT_input, *, IOSTAT=iostatus) istream
      if (iostatus .lt. 0) exit
      if (trim(istream) .eq. '<'//trim(blockname)//'>') then
        ! found the right <blockname>
        find_varname: do while (.true.)
          read (UNIT_input, *, IOSTAT=iostatus) istream
          if (iostatus .lt. 0) exit find_blockname
          ! exit if you reach the beginning of a new block:
          if (istream(1:1) .eq. '<') exit find_blockname
          if (trim(istream) .eq. trim(varname)) then
            ! found the right variable
            backspace (UNIT_input)
            read (UNIT_input, '(a)') istream
            k1 = scan(istream, '=')
            k2 = scan(istream, '#')
            if (k2 .eq. 0) k2 = STR_MAX
            value_str = istream(k1 + 1:k2 - 1)
            value_str = adjustl(value_str)
            exit find_blockname
          end if
        end do find_varname
      end if
    end do find_blockname
    close (UNIT_input)
    if (k2 .eq. 0) then
      found = .false.
      parseInput = ''
    else
      found = .true.
      parseInput = value_str
    end if
  end function parseInput

  function isUniqueParam(blockname, varname) result(unique)
    implicit none
    character(len=*), intent(in) :: blockname, varname
    logical :: unique
    integer :: n
    do n = 1, sim_params % count
      if ((trim(blockname) .eq. trim(sim_params % param_group(n) % str)) .and. &
          (trim(varname) .eq. trim(sim_params % param_name(n) % str))) then
        unique = .false.
        return
      end if
    end do
    unique = .true.
    return
  end function

  subroutine getInt4Input(blockname, varname, val, def_val)
    implicit none
    character(len=*), intent(in) :: blockname, varname
    integer(kind=4), optional :: def_val
    integer(kind=4), intent(out) :: val
    character(len=STR_MAX) :: val_str
    logical :: found
    integer :: iostatus

    val_str = parseInput(blockname, varname, found)
    if (.not. found) then
      if (present(def_val)) then
        val = def_val
      else
        call throwError("ERROR variable `"//trim(varname) &
                        //"` not found in <"//trim(blockname)//">"// &
                        END_LINE &
                        //"-- no default value provided")
      end if
    else
      call strToInt4(val_str, val, iostatus)
      if (iostatus .ne. 0) then
        if (present(def_val)) then
          val = def_val
        else
          call throwError("ERROR variable `"//trim(varname)//"`="//trim(val_str) &
                          //" from <"//trim(blockname)//"> not converted"// &
                          END_LINE &
                          //"-- no default value provided")
        end if
      end if
    end if
    if (isUniqueParam(blockname, varname)) then
      sim_params % count = sim_params % count + 1
      sim_params % param_type(sim_params % count) = 1
      sim_params % param_group(sim_params % count) % str = blockname
      sim_params % param_name(sim_params % count) % str = varname
      sim_params % param_value(sim_params % count) % value_int = val
    end if
  end subroutine getInt4Input
  subroutine strToInt4(val_str, val, stat)
    implicit none
    character(len=*), intent(in) :: val_str
    integer(kind=4), intent(out) :: val
    integer, intent(out) :: stat
    real(kind=8) :: val_
    call strToReal8(val_str, val_, stat)
    val = INT(val_)
  end subroutine strToInt4

  subroutine getInt8Input(blockname, varname, val, def_val)
    implicit none
    character(len=*), intent(in) :: blockname, varname
    integer(kind=8), optional :: def_val
    integer(kind=8), intent(out) :: val
    character(len=STR_MAX) :: val_str
    logical :: found
    integer :: iostatus

    val_str = parseInput(blockname, varname, found)
    if (.not. found) then
      if (present(def_val)) then
        val = def_val
      else
        call throwError("ERROR variable `"//trim(varname) &
                        //"` not found in <"//trim(blockname)//">"// &
                        END_LINE &
                        //"-- no default value provided")
      end if
    else
      call strToInt8(val_str, val, iostatus)
      if (iostatus .ne. 0) then
        if (present(def_val)) then
          val = def_val
        else
          call throwError("ERROR variable `"//trim(varname)//"`="//trim(val_str) &
                          //" from <"//trim(blockname)//"> not converted"// &
                          END_LINE &
                          //"-- no default value provided")
        end if
      end if
    end if
    if (isUniqueParam(blockname, varname)) then
      sim_params % count = sim_params % count + 1
      sim_params % param_type(sim_params % count) = 1
      sim_params % param_group(sim_params % count) % str = blockname
      sim_params % param_name(sim_params % count) % str = varname
      sim_params % param_value(sim_params % count) % value_int = INT(val, 4)
    end if
  end subroutine getInt8Input
  subroutine strToInt8(val_str, val, stat)
    implicit none
    character(len=*), intent(in) :: val_str
    integer(kind=8), intent(out) :: val
    integer, intent(out) :: stat
    real(kind=8) :: val_
    call strToReal8(val_str, val_, stat)
    val = INT(val_, KIND(8))
  end subroutine strToInt8

  subroutine getReal4Input(blockname, varname, val, def_val)
    implicit none
    character(len=*), intent(in) :: blockname, varname
    real(kind=4), optional :: def_val
    real(kind=4), intent(out) :: val
    character(len=STR_MAX) :: val_str
    logical :: found
    integer :: iostatus

    val_str = parseInput(blockname, varname, found)
    if (.not. found) then
      if (present(def_val)) then
        val = def_val
      else
        call throwError("ERROR variable `"//trim(varname) &
                        //"` not found in <"//trim(blockname)//">"// &
                        END_LINE &
                        //"-- no default value provided")
      end if
    else
      call strToReal4(val_str, val, iostatus)
      if (iostatus .ne. 0) then
        if (present(def_val)) then
          val = def_val
        else
          call throwError("ERROR variable `"//trim(varname)//"`="//trim(val_str) &
                          //" from <"//trim(blockname)//"> not converted"// &
                          END_LINE &
                          //"-- no default value provided")
        end if
      end if
    end if
    if (isUniqueParam(blockname, varname)) then
      sim_params % count = sim_params % count + 1
      sim_params % param_type(sim_params % count) = 2
      sim_params % param_group(sim_params % count) % str = blockname
      sim_params % param_name(sim_params % count) % str = varname
      sim_params % param_value(sim_params % count) % value_real = val
    end if
  end subroutine getReal4Input
  subroutine strToReal4(val_str, val, stat)
    implicit none
    character(len=*), intent(in) :: val_str
    real(kind=4), intent(out) :: val
    integer, intent(out) :: stat
    real(kind=8) :: val_
    call strToReal8(val_str, val_, stat)
    val = REAL(val_, 4)
  end subroutine strToReal4

  subroutine getReal8Input(blockname, varname, val, def_val)
    implicit none
    character(len=*), intent(in) :: blockname, varname
    real(kind=8), optional :: def_val
    real(kind=8), intent(out) :: val
    character(len=STR_MAX) :: val_str
    logical :: found
    integer :: iostatus

    val_str = parseInput(blockname, varname, found)
    if (.not. found) then
      if (present(def_val)) then
        val = def_val
      else
        call throwError("ERROR variable `"//trim(varname) &
                        //"` not found in <"//trim(blockname)//">"// &
                        END_LINE &
                        //"-- no default value provided")
      end if
    else
      call strToReal8(val_str, val, iostatus)
      if (iostatus .ne. 0) then
        if (present(def_val)) then
          val = def_val
        else
          call throwError("ERROR variable `"//trim(varname)//"`="//trim(val_str) &
                          //" from <"//trim(blockname)//"> not converted"// &
                          END_LINE &
                          //"-- no default value provided")
        end if
      end if
    end if
    if (isUniqueParam(blockname, varname)) then
      sim_params % count = sim_params % count + 1
      sim_params % param_type(sim_params % count) = 2
      sim_params % param_group(sim_params % count) % str = blockname
      sim_params % param_name(sim_params % count) % str = varname
      sim_params % param_value(sim_params % count) % value_real = REAL(val, 4)
    end if
  end subroutine getReal8Input
  subroutine strToReal8(val_str, val, stat)
    implicit none
    character(len=*), intent(in) :: val_str
    real(kind=8), intent(out) :: val
    integer, intent(out) :: stat
    read (val_str, *, iostat=stat) val
  end subroutine strToReal8

  subroutine getLogicalInput(blockname, varname, val, def_val)
    implicit none
    character(len=*), intent(in) :: blockname, varname
    logical, optional :: def_val
    logical, intent(out) :: val
    character(len=STR_MAX) :: val_str
    logical :: found
    integer :: val_int, iostatus

    val_str = parseInput(blockname, varname, found)
    if (.not. found) then
      if (present(def_val)) then
        val = def_val
      else
        call throwError("ERROR variable `"//trim(varname) &
                        //"` not found in <"//trim(blockname)//">"// &
                        END_LINE &
                        //"-- no default value provided")
      end if
    else
      call strToInt4(val_str, val_int, iostatus)
      if (iostatus .ne. 0) then
        if (present(def_val)) then
          val = def_val
        else
          call throwError("ERROR variable `"//trim(varname)//"`="//trim(val_str) &
                          //" from <"//trim(blockname)//"> not converted"// &
                          END_LINE &
                          //"-- no default value provided")
        end if
      else
        if (val_int .eq. 0) then
          val = .false.
        else if (val_int .eq. 1) then
          val = .true.
        else
          call throwError("ERROR variable `"//trim(varname)//"`="//trim(val_str) &
                          //" from <"//trim(blockname)//"> wrong logical value")
        end if
      end if
    end if
    if (isUniqueParam(blockname, varname)) then
      sim_params % count = sim_params % count + 1
      sim_params % param_type(sim_params % count) = 3
      sim_params % param_group(sim_params % count) % str = blockname
      sim_params % param_name(sim_params % count) % str = varname
      sim_params % param_value(sim_params % count) % value_bool = val
    end if
  end subroutine getLogicalInput
end module m_readinput
