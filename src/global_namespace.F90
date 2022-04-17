#include "defs.F90"

!--- GLOBAL_NAMESPACE ------------------------------------------!
! To store predefined(!) parameters and functions/interfaces
!   and share them between modules
!...............................................................!

module m_globalnamespace
  #ifdef MPI08
    use mpi_f08
  #endif

  #ifdef MPI
    include "mpif.h"
  #endif

  integer, parameter      :: dprec = kind(1.0d0)
  integer, parameter      :: sprec = kind(1.0e0)
  integer, parameter      :: UNIT_input = 10, UNIT_output = 20, UNIT_history = 30, UNIT_diag = 40, UNIT_warn = 80
  integer, parameter      :: UNIT_params = 50, UNIT_usrout = 90
  integer, parameter      :: UNIT_restart_fld = 60, UNIT_restart_prtl = 70

  ! algorithm specific parameters
  real                    :: CC, CCINV, CORR

  #ifdef GCA
    real                    :: gca_eoverbmin, gca_rhomin, gca_vperpmax
    logical                 :: gca_enforce_mu0
  #endif

  logical                 :: enable_fieldsolver
  logical                 :: enable_currentdeposit

  ! plasma parameters
  real :: ppc0, c_omp, sigma, B_norm, unit_ch

  ! simulation parameters
  integer                :: start_timestep = 0, final_timestep, warning_count = 0
  logical                :: resize_tiles
  integer                :: min_tile_nprt = 100, max_buffsize = 100
  character(len=STR_MAX) :: input_file_name = 'input',&
                          & output_dir_name = 'output',&
                          & slice_dir_name = 'slices',&
                          & restart_dir_name = 'restart',&
                          & restart_from = 'restart/step_00000',&
                          & diag_file_name = 'diag.log',&
                          & warn_file_name = 'warn.log'

  integer       :: nfilter

  ! mpi variables
  integer       :: mpi_rank, mpi_size, mpi_statsize
  integer       :: sizex, sizey, sizez
  #ifdef HDF5
    integer, parameter  :: UNIT_xdmf = 40
    integer             :: h5comm, h5info
  #endif
  
  ! variables visible globally defined by the user
  real                    :: global_usr_variable_1, global_usr_variable_2, global_usr_variable_3
contains
  subroutine renormalizeUnits()
    implicit none
    CCINV = 1.0 / CC
    B_norm = CC**2 * sqrt(sigma) / c_omp
    unit_ch = CC**2 / (ppc0 * c_omp**2)
  end subroutine renormalizeUnits
end module m_globalnamespace
