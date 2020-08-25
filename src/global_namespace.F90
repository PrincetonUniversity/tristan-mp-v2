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
  integer, parameter      :: UNIT_input = 10, UNIT_output = 20, UNIT_history = 30
  integer, parameter      :: UNIT_params = 50
  integer, parameter      :: UNIT_restart_fld = 60, UNIT_restart_prtl = 70

  ! algorithm specific parameters
  real                    :: CC, CCINV, CORR

  logical                 :: enable_fieldsolver
  logical                 :: enable_currentdeposit

  ! plasma parameters
  real :: ppc0, c_omp, sigma, B_norm, unit_ch

  ! simulation parameters
  integer                :: start_timestep = 0, final_timestep, output_index = 0, slice_index = 0
  logical                :: resize_tiles
  integer                :: min_tile_nprt = 100
  character(len=STR_MAX) :: input_file_name = 'input',&
                          & output_dir_name = 'output',&
                          & restart_dir_name = 'restart',&
                          & slice_dir_name = 'slices'

  integer       :: nfilter, spec_num
  real          :: spec_min, spec_max
  logical       :: spec_log_bins

  ! mpi variables
  integer       :: mpi_rank, mpi_size, mpi_statsize
  integer       :: sizex, sizey, sizez
  #ifdef HDF5
    integer, parameter  :: UNIT_xdmf = 40
    integer             :: h5comm, h5info
  #endif
contains
  subroutine renormalizeUnits()
    implicit none
    CCINV = 1.0 / CC
    B_norm = CC**2 * sqrt(sigma) / c_omp
    unit_ch = CC**2 / (ppc0 * c_omp**2)
  end subroutine renormalizeUnits
end module m_globalnamespace
