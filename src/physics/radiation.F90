#include "../defs.F90"

module m_radiation
#ifdef RADIATION

  use m_globalnamespace
  use m_aux
  use m_errors
  use m_particlelogistics
  implicit none

  real              :: emit_gamma_syn, emit_gamma_ic, cool_gamma_syn, cool_gamma_ic, rad_beta_rec
  real              :: rad_dens_lim
  real, allocatable :: rad_spectra(:,:), glob_rad_spectra(:,:)
  integer           :: rad_photon_sp

  !--- PRIVATE variables/functions -------------------------------!
  !...............................................................!
contains
  subroutine particleRadiateSync(s,&
                               & u0, v0, w0, ui, vi, wi,&
                               & dx, dy, dz, xi, yi, zi,&
                               & bx, by, bz, ex, ey, ez)
    implicit none
    real, intent(inout)           :: u0, v0, w0
    real, intent(in)              :: ui, vi, wi
    real, intent(in)              :: bx, by, bz, ex, ey, ez
    real, intent(in)              :: dx, dy, dz
    integer(kind=2), intent(in)   :: xi, yi, zi
    integer, intent(in)           :: s

    ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
    !
    !
    ! THIS ROUTINE IS STILL UNDER CONSTRUCTION
    !   if you need access, please contact developers
    !
    !
    ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  end subroutine particleRadiateSync

  subroutine particleRadiateIC(s,&
                             & u0, v0, w0, ui, vi, wi,&
                             & dx, dy, dz, xi, yi, zi,&
                             & bx, by, bz, ex, ey, ez)
    implicit none
    real, intent(inout)           :: u0, v0, w0
    real, intent(in)              :: ui, vi, wi
    real, intent(in)              :: bx, by, bz, ex, ey, ez
    real, intent(in)              :: dx, dy, dz
    integer(kind=2), intent(in)   :: xi, yi, zi
    integer, intent(in)           :: s

    ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
    !
    !
    ! THIS ROUTINE IS STILL UNDER CONSTRUCTION
    !   if you need access, please contact developers
    !
    !
    ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  end subroutine particleRadiateIC

#endif
end module m_radiation
