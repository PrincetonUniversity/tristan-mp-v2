module m_writespectra
#ifdef HDF5
  use hdf5
  use m_globalnamespace, only: h5comm, h5info
#endif
  use m_globalnamespace, only: output_dir_name, mpi_rank, output_dir_spec
  use m_outputnamespace, only: glob_spectra, spec_bin_size, spec_nx, spec_ny, spec_nz, &
                               spec_num, spec_min, spec_max, spec_log_bins
  use m_aux
  use m_errors, only: throwError
  use m_domain
  use m_fields
  use m_particles
  use m_exchangearray, only: exchangeArray

#ifdef RADIATION
  use m_outputnamespace, only: glob_rad_spectra, rad_spec_num, rad_spec_min, rad_spec_bin_size
#endif

  implicit none

  !--- PRIVATE functions -----------------------------------------!
#ifdef HDF5
  private :: writeSpectra_hdf5
#endif
  !...............................................................!

contains
  subroutine writeSpectra(step, time)
    implicit none
    integer, intent(in) :: step, time
#ifdef HDF5
    call writeSpectra_hdf5(step, time)
#endif
    call printDiag("writeSpectra()", 3)
  end subroutine writeSpectra

#ifdef HDF5
  subroutine writeSpectra_hdf5(step, time)
    implicit none
    integer, intent(in) :: step, time
    character(len=STR_MAX) :: stepchar, filename
    integer :: error, s, i
    character(len=6) :: dsetname
    integer(HID_T) :: file_id, dset_id, dspace_id
    integer(HSIZE_T) :: bin_dims(1), spec_dims(4)
#ifdef RADIATION
    integer(HSIZE_T) :: rad_dims(1)
#endif
    real, allocatable, dimension(:) :: bin_data
    integer :: root_rnk = 0
    real, allocatable :: xbin_data(:), ybin_data(:), zbin_data(:)

    if (.false.) print *, time

    ! only root rank writes the spectra file
    if (mpi_rank .eq. root_rnk) then
      ! saving the energy bins
      allocate (bin_data(spec_num))
      do i = 1, spec_num
        bin_data(i) = spec_min + REAL(i - 0.5) * spec_bin_size
        if (spec_log_bins) then
          bin_data(i) = exp(bin_data(i))
        end if
      end do
      ! saving the spatial bins
      allocate (xbin_data(spec_nx))
      allocate (ybin_data(spec_ny))
      allocate (zbin_data(spec_nz))
      do i = 1, spec_nx
        xbin_data(i) = (i - 1) * REAL(global_mesh % sx) / REAL(spec_nx)
      end do
      do i = 1, spec_ny
        ybin_data(i) = (i - 1) * REAL(global_mesh % sy) / REAL(spec_ny)
      end do
      do i = 1, spec_nz
        zbin_data(i) = (i - 1) * REAL(global_mesh % sz) / REAL(spec_nz)
      end do

      write (stepchar, "(i5.5)") step
      filename = trim(output_dir_spec)//'/spec.tot.'//trim(stepchar)

      ! Initialize FORTRAN interface
      call h5open_f(error)
      ! Create a new file using default properties
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      ! writing bins:
      bin_dims(1) = spec_num
      call h5screate_simple_f(1, bin_dims, dspace_id, error)
      call h5dcreate_f(file_id, 'ebins', default_h5_real, dspace_id, dset_id, error)
      call h5dwrite_f(dset_id, default_h5_real, bin_data, bin_dims, error)
      call h5dclose_f(dset_id, error)
      call h5sclose_f(dspace_id, error)

      bin_dims(1) = spec_nx
      call h5screate_simple_f(1, bin_dims, dspace_id, error)
      call h5dcreate_f(file_id, 'xbins', default_h5_real, dspace_id, dset_id, error)
      call h5dwrite_f(dset_id, default_h5_real, xbin_data, bin_dims, error)
      call h5dclose_f(dset_id, error)
      call h5sclose_f(dspace_id, error)

      bin_dims(1) = spec_ny
      call h5screate_simple_f(1, bin_dims, dspace_id, error)
      call h5dcreate_f(file_id, 'ybins', default_h5_real, dspace_id, dset_id, error)
      call h5dwrite_f(dset_id, default_h5_real, ybin_data, bin_dims, error)
      call h5dclose_f(dset_id, error)
      call h5sclose_f(dspace_id, error)

      bin_dims(1) = spec_nz
      call h5screate_simple_f(1, bin_dims, dspace_id, error)
      call h5dcreate_f(file_id, 'zbins', default_h5_real, dspace_id, dset_id, error)
      call h5dwrite_f(dset_id, default_h5_real, zbin_data, bin_dims, error)
      call h5dclose_f(dset_id, error)
      call h5sclose_f(dspace_id, error)

#ifdef RADIATION
      rad_dims(1) = rad_spec_num
      if (allocated(bin_data)) deallocate (bin_data)
      allocate (bin_data(rad_spec_num))
      do i = 1, rad_spec_num
        bin_data(i) = rad_spec_min + REAL(i - 0.5) * rad_spec_bin_size
        if (spec_log_bins) then
          bin_data(i) = exp(bin_data(i))
        end if
      end do
      ! writing bins:
      dsetname = 'rbins'
      call h5screate_simple_f(1, rad_dims, dspace_id, error)
      call h5dcreate_f(file_id, dsetname, default_h5_real, dspace_id, &
                       dset_id, error)
      call h5dwrite_f(dset_id, default_h5_real, bin_data, rad_dims, error)
      call h5dclose_f(dset_id, error)
      call h5sclose_f(dspace_id, error)
#endif

      spec_dims(1) = spec_nx
      spec_dims(2) = spec_ny
      spec_dims(3) = spec_nz
      spec_dims(4) = spec_num
      do s = 1, nspec
        ! writing spectra:
        dsetname = 'n'//trim(STR(s))
        call h5screate_simple_f(4, spec_dims, dspace_id, error)
        call h5dcreate_f(file_id, dsetname, default_h5_real, dspace_id, &
                         dset_id, error)
        call h5dwrite_f(dset_id, default_h5_real, glob_spectra(s, :, :, :, :), spec_dims, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

#ifdef RADIATION
        ! writing spectra:
        dsetname = 'nr'//trim(STR(s))
        call h5screate_simple_f(1, rad_dims, dspace_id, error)
        call h5dcreate_f(file_id, dsetname, default_h5_real, dspace_id, &
                         dset_id, error)
        call h5dwrite_f(dset_id, default_h5_real, glob_rad_spectra(s, :), rad_dims, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)
        glob_rad_spectra(s, :) = 0.0
#endif
      end do

      ! Close the file
      call h5fclose_f(file_id, error)
      ! Close FORTRAN interface
      call h5close_f(error)

      if (allocated(bin_data)) deallocate (bin_data)
      if (allocated(xbin_data)) deallocate (xbin_data)
      if (allocated(ybin_data)) deallocate (ybin_data)
      if (allocated(zbin_data)) deallocate (zbin_data)

      if (allocated(glob_spectra)) deallocate (glob_spectra)

#ifdef RADIATION
      if (allocated(glob_rad_spectra)) deallocate (glob_rad_spectra)
#endif
    end if
  end subroutine writeSpectra_hdf5
#endif

end module m_writespectra
