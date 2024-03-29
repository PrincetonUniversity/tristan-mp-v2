!--- GLOBAL_NAMESPACE ------------------------------------------!
! To store predefined(!) parameters and functions/interfaces
!   and share them between modules
!...............................................................!

module m_globalnamespace
#ifdef MPI08
  use mpi_f08
#endif

#ifdef HDF5
  use hdf5
#endif

#ifdef MPI
  include "mpif.h"
#endif

  integer, parameter :: dprec = kind(1.0d0)
  integer, parameter :: sprec = kind(1.0e0)
  integer, parameter :: UNIT_input = 10, UNIT_output = 20, UNIT_diag = 40, UNIT_warn = 80
  integer, parameter :: UNIT_params = 50, UNIT_usrout = 90
  integer, parameter :: UNIT_restart = 60

  ! algorithm specific parameters
  real :: CC, CCINV, CORR

  logical :: enable_fieldsolver
  logical :: enable_currentdeposit

  ! plasma parameters
  real :: ppc0, c_omp, sigma, B_norm, unit_ch

  ! simulation parameters
  integer :: start_timestep = 0, final_timestep, warning_count = 0
  logical :: resize_tiles, shrink_tiles
  integer :: min_tile_nprt = 100, max_buffsize = 32
  integer :: t_max_check_interval = 10
  character(len=STR_MAX) :: input_file_name = 'input', &
                            output_dir_name = 'output', &
                            output_dir_spec = 'output/spec', &
                            output_dir_flds = 'output/flds', &
                            output_dir_prtl = 'output/prtl', &
                            slice_dir_name = 'slices', &
                            restart_dir_name = 'restart', &
                            restart_from = 'restart/step_000000', &
                            diag_file_name = 'diag.log', &
                            warn_file_name = 'warn.log'

  integer :: nfilter
  real(kind=8) :: wall_t_max

  ! mpi variables
  integer :: mpi_rank, mpi_size, mpi_statsize
  integer :: sizex, sizey, sizez

#ifdef HDF5
  integer, parameter :: UNIT_xdmf = 40
  integer :: h5comm, h5info
#endif

#ifdef MPI08
  type(MPI_Datatype) :: default_mpi_real
#endif

#ifdef MPI
  integer :: default_mpi_real
#endif

  ! variables visible globally defined by the user
  real :: global_usr_variable_1, global_usr_variable_2, global_usr_variable_3
contains
  subroutine renormalizeUnits()
    implicit none
    CCINV = 1.0 / CC
    B_norm = CC**2 * sqrt(sigma) / c_omp
    unit_ch = CC**2 / (ppc0 * c_omp**2)
  end subroutine renormalizeUnits
end module m_globalnamespace
