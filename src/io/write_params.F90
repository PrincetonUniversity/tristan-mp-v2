module m_writeparams
#ifdef HDF5
  use hdf5
  use m_globalnamespace, only: h5comm, h5info
#endif
  use m_globalnamespace, only: output_dir_name, mpi_rank
  use m_aux
  use m_errors
  use m_domain
  use m_fields
  implicit none

  !--- PRIVATE functions -----------------------------------------!
#ifdef HDF5
  private :: writeParams_hdf5
#else
  private :: writeParams_fmt
#endif
  !...............................................................!
contains

  subroutine writeParams(step, time)
    implicit none
    integer, intent(in) :: step, time
#ifdef HDF5
    call writeParams_hdf5(step, time)
#else
    call writeParams_fmt(step, time)
#endif

    call printDiag("writeParams()", 3)
  end subroutine writeParams

#ifdef HDF5
  subroutine writeParams_hdf5(step, time)
    implicit none
    integer, intent(in) :: step, time
    character(len=STR_MAX) :: stepchar, filename
    integer :: n, error, datarank
    integer(HID_T) :: file_id, dspace_id, dset_id
    integer(HSIZE_T), dimension(1) :: data_dims
    integer, allocatable :: data_int(:)
    real, allocatable :: data_real(:)
    character(len=STR_MAX) :: dsetname

    if (mpi_rank .eq. 0) then
      datarank = 1
      data_dims(1) = 1
      allocate (data_real(1))
      allocate (data_int(1))

      write (stepchar, "(i5.5)") step
      filename = trim(output_dir_name)//'/params.'//trim(stepchar)

      dsetname = 'timestep'
      call h5open_f(error)
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      call h5screate_simple_f(datarank, data_dims, dspace_id, error)
      call h5dcreate_f(file_id, trim(dsetname), H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
      data_int(1) = time
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data_int, data_dims, error)
      call h5dclose_f(dset_id, error)
      call h5sclose_f(dspace_id, error)

      do n = 1, sim_params % count
        dsetname = trim(sim_params % param_group(n) % str)//':'//trim(sim_params % param_name(n) % str)
        call h5screate_simple_f(datarank, data_dims, dspace_id, error)
        if (sim_params % param_type(n) .eq. 1) then
          call h5dcreate_f(file_id, trim(dsetname), H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
          data_int(1) = sim_params % param_value(n) % value_int
          call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data_int, data_dims, error)
        else if (sim_params % param_type(n) .eq. 2) then
          call h5dcreate_f(file_id, trim(dsetname), default_h5_real, dspace_id, dset_id, error)
          data_real(1) = sim_params % param_value(n) % value_real
          call h5dwrite_f(dset_id, default_h5_real, data_real, data_dims, error)
        else if (sim_params % param_type(n) .eq. 3) then
          call h5dcreate_f(file_id, trim(dsetname), H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
          if (sim_params % param_value(n) % value_bool) then
            data_int(1) = 1
          else
            data_int(1) = 0
          end if
          call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data_int, data_dims, error)
        else
          call throwError('ERROR. Unknown `param_type` in `saveAllParameters`.')
        end if
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)
      end do
      call h5fclose_f(file_id, error)
      call h5close_f(error)
    end if
  end subroutine writeParams_hdf5
#else
  subroutine writeParams_fmt(step, time)
    implicit none
    integer, intent(in) :: step, time
    integer :: n
    character(len=STR_MAX) :: FMT
    character(len=STR_MAX) :: filename, stepchar

    if (mpi_rank .eq. 0) then
      write (stepchar, "(i5.5)") step
      filename = trim(output_dir_name)//'/params.'//trim(stepchar)
      open (UNIT_params, file=filename, status="replace", access="stream", form="formatted")
      FMT = '(A52,I10)'
      write (UNIT_params, FMT) 'timestep', time

      do n = 1, sim_params % count
        if (sim_params % param_type(n) .eq. 1) then
          FMT = '(A30,A1,A20,A1,I10)'
          write (UNIT_params, FMT) trim(sim_params % param_group(n) % str), ':', &
            trim(sim_params % param_name(n) % str), ':', &
            sim_params % param_value(n) % value_int
        else if (sim_params % param_type(n) .eq. 2) then
          FMT = getFMTForReal(sim_params % param_value(n) % value_real)
          FMT = '(A30,A1,A20,A1,'//trim(FMT)//')'
          write (UNIT_params, FMT) trim(sim_params % param_group(n) % str), ':', &
            trim(sim_params % param_name(n) % str), ':', &
            sim_params % param_value(n) % value_real
        else if (sim_params % param_type(n) .eq. 3) then
          FMT = '(A30,A1,A20,A1,L10)'
          write (UNIT_params, FMT) trim(sim_params % param_group(n) % str), ':', &
            trim(sim_params % param_name(n) % str), ':', &
            sim_params % param_value(n) % value_bool
        else
          call throwError('ERROR. Unknown `param_type` in `saveAllParameters`.')
        end if
      end do
      close (UNIT_params)
    end if
  end subroutine writeParams_fmt
#endif

end module m_writeparams
