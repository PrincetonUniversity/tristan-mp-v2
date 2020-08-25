#include "../defs.F90"

module m_writehelpers
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_helpers
  use m_exchangearray, only: exchangeArray
  implicit none
contains
  ! writes a field specified by `fld_var` from gridcell `i,j,k` ...
  ! ... to `sm_arr(i1, j1, k1)` with proper interpolation etc for the output
  subroutine selectFieldForOutput(fld_var, i1, j1, k1, i, j, k, writing_lgarrQ)
    implicit none
    character(len=STR_MAX), intent(in)  :: fld_var
    integer(kind=2), intent(in)         :: i1, j1, k1, i, j, k
    logical, intent(in)                 :: writing_lgarrQ
    real                                :: ex0, ey0, ez0, bx0, by0, bz0, jx0, jy0, jz0
    select case (trim(fld_var))
    case('ex')
      #ifndef debug
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
      #else
        ex0 = ex(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = ex0 * B_norm
    case('ey')
      #ifndef debug
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
      #else
        ey0 = ey(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = ey0 * B_norm
    case('ez')
      #ifndef debug
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
      #else
        ez0 = ez(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = ez0 * B_norm
    case('bx')
      #ifndef debug
        call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
      #else
        bx0 = bx(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = bx0 * B_norm
    case('by')
      #ifndef debug
        call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
      #else
        by0 = by(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = by0 * B_norm
    case('bz')
      #ifndef debug
        call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
      #else
        bz0 = bz(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = bz0 * B_norm
    case('jx')
      #ifndef debug
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
      #else
        jx0 = jx(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = -jx0 * B_norm
    case('jy')
      #ifndef debug
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
      #else
        jy0 = jy(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = -jy0 * B_norm
    case('jz')
      #ifndef debug
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
      #else
        jz0 = jz(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = -jz0 * B_norm
    case('xx')
      sm_arr(i1, j1, k1) = REAL(this_meshblock%ptr%x0 + i, 4)
    case('yy')
      sm_arr(i1, j1, k1) = REAL(this_meshblock%ptr%y0 + j, 4)
    case('zz')
      sm_arr(i1, j1, k1) = REAL(this_meshblock%ptr%z0 + k, 4)
    case default
      if (((fld_var(1:4) .ne. 'dens') .and.&
         & (fld_var(1:4) .ne. 'enrg') .and.&
         & (fld_var(1:4) .ne. 'dgca')) .or.&
         & (.not. writing_lgarrQ)) then
        call throwError("ERROR: unrecognized `fldname`")
      else
        sm_arr(i1, j1, k1) = lg_arr(i, j, k)
      end if
    end select
  end subroutine selectFieldForOutput

  subroutine prepareFieldForOutput(fldname, writing_lgarrQ)
    implicit none
    character(len=STR_MAX), intent(in)    :: fldname
    logical, intent(out)                  :: writing_lgarrQ
    integer                               :: s

    if (fldname(1:4) .eq. 'dens') then
      writing_lgarrQ = .true.
      s = STRtoINT(fldname(5:5))
      ! fill `lg_arr` with density of species `s`
      #ifndef DEBUG
        call computeDensity(s, reset=.true.)
      #else
        call computeDensity(s, reset=.true., ds=0)
      #endif
      call exchangeArray()
    else if (fldname(1:4) .eq. 'enrg') then
      writing_lgarrQ = .true.
      s = STRtoINT(fldname(5:5))
      ! fill `lg_arr` with energy density of species `s`
      #ifndef DEBUG
        call computeEnergy(s, reset=.true.)
      #else
        call computeEnergy(s, reset=.true., ds=0)
      #endif
      call exchangeArray()
    else if (fldname(1:4) .eq. 'dgca') then
      writing_lgarrQ = .true.
      s = STRtoINT(fldname(5:5))
      call throwError('ERROR: `dgca` not defined without GCA flag.')
    else
      writing_lgarrQ = .false.
    end if
  end subroutine prepareFieldForOutput

end module m_writehelpers
