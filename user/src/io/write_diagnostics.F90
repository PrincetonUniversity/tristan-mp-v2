module m_writediagnostics
#ifdef HDF5
  use hdf5
  use m_globalnamespace, only: h5comm, h5info
#endif
  use m_globalnamespace, only: output_dir_name, mpi_rank, mpi_size
  use m_outputnamespace, only: n_dom_vars, dom_vars
  use m_aux
  use m_errors, only: throwError
  use m_domain

  implicit none

  !--- PRIVATE functions -----------------------------------------!
#ifdef HDF5
  private :: writeDomain_hdf5
#endif
  !...............................................................!

contains
  subroutine writeDiagnostics(step, time)
    implicit none
    integer, intent(in) :: step, time
#ifdef HDF5
    call writeDomain_hdf5(step, time)
#endif
    call printDiag("writeDiagnostics()", 3)
  end subroutine writeDiagnostics

#ifdef HDF5
  subroutine writeDomain_hdf5(step, time)
    implicit none
    integer, intent(in) :: step, time
    character(len=STR_MAX) :: stepchar, filename
    integer :: error, datarank, d, rnk
    integer(HID_T) :: file_id, dset_id, dspace_id
    integer(HSIZE_T), dimension(1) :: data_dims
    integer, allocatable, dimension(:) :: domain_data
    if (.false.) print *, time

    datarank = 1
    data_dims(1) = mpi_size

    ! only root rank writes the domain file
    if (mpi_rank .eq. 0) then
      write (stepchar, "(i5.5)") step
      filename = trim(output_dir_name)//'/domain.'//trim(stepchar)

      ! Initialize FORTRAN interface
      call h5open_f(error)
      ! Create a new file using default properties
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      allocate (domain_data(mpi_size))

      do d = 1, n_dom_vars
        select case (trim(dom_vars(d)))
        case ('x0')
          do rnk = 0, mpi_size - 1
            domain_data(rnk + 1) = meshblocks(rnk + 1) % x0
          end do
        case ('y0')
          do rnk = 0, mpi_size - 1
            domain_data(rnk + 1) = meshblocks(rnk + 1) % y0
          end do
        case ('z0')
          do rnk = 0, mpi_size - 1
            domain_data(rnk + 1) = meshblocks(rnk + 1) % z0
          end do
        case ('sx')
          do rnk = 0, mpi_size - 1
            domain_data(rnk + 1) = meshblocks(rnk + 1) % sx
          end do
        case ('sy')
          do rnk = 0, mpi_size - 1
            domain_data(rnk + 1) = meshblocks(rnk + 1) % sy
          end do
        case ('sz')
          do rnk = 0, mpi_size - 1
            domain_data(rnk + 1) = meshblocks(rnk + 1) % sz
          end do
        case default
          call throwError('ERROR: unrecognized `dom_vars`: `'//trim(dom_vars(d))//'`')
        end select

        call h5screate_simple_f(datarank, data_dims, dspace_id, error)
        call h5dcreate_f(file_id, dom_vars(d), H5T_NATIVE_INTEGER, dspace_id, &
                         dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, domain_data, data_dims, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)
      end do

      ! Close the file
      call h5fclose_f(file_id, error)
      ! Close FORTRAN interface
      call h5close_f(error)

      if (allocated(domain_data)) deallocate (domain_data)
    end if
  end subroutine writeDomain_hdf5
#endif

end module m_writediagnostics
