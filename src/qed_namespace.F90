module m_qednamespace
  use m_globalnamespace
  implicit none

#ifdef RADIATION
  real :: emit_gamma_syn, emit_gamma_ic, cool_gamma_syn, cool_gamma_ic, rad_beta_rec
  real :: rad_dens_lim, rad_cool_lim
  integer :: rad_photon_sp
  integer :: rad_interval
#endif

end module m_qednamespace
