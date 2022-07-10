module m_outputnamespace
  use m_globalnamespace
  implicit none

  ! history
  logical :: hst_enable
  integer :: hst_interval

  ! user-specified diagnostic output
  logical :: usrout_enable
  integer :: usrout_interval

  ! total output and slice
  integer :: tot_output_index = 0, slice_index = 0

  ! input parameters
  integer :: output_dens_smooth     ! density smoothing window
  ! ... for .tot. outputs
  logical :: tot_output_enable
  logical :: params_enable, prtl_tot_enable
  logical :: flds_tot_enable, spectra_enable, diag_enable
  logical :: flds_at_prtl_enable, xdmf_enable
  logical :: derivatives_enable, npart_enable, T0i_output_enable, Tij_output_enable, Tii_output_enable
  integer :: tot_output_start, tot_output_interval
  integer :: tot_output_stride      ! particle striding
  integer :: output_flds_istep      ! field downsampling for .tot.
  ! ... for `slice` outputs
  logical :: slice_output_enable
  integer :: slice_output_start, slice_output_interval
  ! ... for `spectra`
  integer :: spec_num
  real :: spec_min, spec_max
  logical :: spec_log_bins, spec_dynamic_bins
  integer :: spec_nx, spec_ny, spec_nz
#ifdef RADIATION
  integer :: rad_spec_num
  real :: rad_spec_min, rad_spec_max
#endif

  ! variables
  ! ... for particle/diag output
  integer :: n_prtl_vars, n_dom_vars
  character(len=STR_MAX) :: prtl_vars(100), prtl_var_types(100), dom_vars(100)

  ! ... for spectra
  real :: spec_bin_size
  real, allocatable :: glob_spectra(:, :, :, :, :)

#ifdef RADIATION
  real :: rad_spec_bin_size
  real, allocatable, dimension(:, :) :: rad_spectra, glob_rad_spectra
#endif

  ! ... for `slice` outputs
  integer :: nslices = 0, slice_axes(100), slice_pos(100)

  ! ... for generic field output
  character(len=STR_MAX) :: fld_vars(100)
  integer :: n_fld_vars
end module m_outputnamespace
