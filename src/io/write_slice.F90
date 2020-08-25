#include "../defs.F90"

module m_writeslice
  #ifdef HDF5
    use hdf5
  #endif
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_helpers
  use m_writehelpers, only: prepareFieldForOutput, selectFieldForOutput

  implicit none

  logical                           :: slice_enable

  integer                           :: nslices = 0
  integer                           :: slice_axes(100)
  integer                           :: slice_pos(100)

  integer                           :: slice_start, slice_interval
  integer, private                  :: n_fld_vars
  character(len=STR_MAX), private   :: fld_vars(100)

  !--- PRIVATE functions -----------------------------------------!
  #ifdef SLICE
    private :: writeSliceX_hdf5, writeSliceY_hdf5, writeSliceZ_hdf5!, writeXDMF_hdf5
  #endif
  private :: initializeSliceOutput
  !...............................................................!

contains
  subroutine writeSlices(time)
    implicit none
    integer, intent(in)         :: time
    integer                     :: step, ierr
    integer                     :: n

    call initializeSliceOutput()

    step = slice_index
    #ifdef SLICE
      do n = 1, nslices
        if (slice_axes(n) .eq. 1) then
          call writeSliceX_hdf5(step, time, slice_pos(n))
            call printDiag((mpi_rank .eq. 0), "...writeSliceX_hdf5()", .true.)
        else if (slice_axes(n) .eq. 2) then
          call writeSliceY_hdf5(step, time, slice_pos(n))
            call printDiag((mpi_rank .eq. 0), "...writeSliceY_hdf5()", .true.)
        else if (slice_axes(n) .eq. 3) then
          call writeSliceZ_hdf5(step, time, slice_pos(n))
            call printDiag((mpi_rank .eq. 0), "...writeSliceZ_hdf5()", .true.)
        else
          call throwError('ERROR. wrong `slice_axes(n)` in `writeSlices`')
        end if
      end do
        call printDiag((mpi_rank .eq. 0), "...writeSlices_hdf5()", .true.)
    #endif
    call printDiag((mpi_rank .eq. 0), "slices()", .true.)
    slice_index = slice_index + 1
  end subroutine writeSlices

  subroutine initializeSliceOutput()
    ! DEP_PRT [particle-dependent]
    implicit none
    integer                   :: s
    integer                   :: ierr, ndown
    ! initialize field variables
    !   total number of fields (excluding particle densities)
    n_fld_vars = 12
    n_fld_vars = n_fld_vars + 2 * nspec
    do s = 1, nspec
      ! hopefully less than 10 species
      fld_vars(s) = 'dens' // STR(s)
    end do
    do s = 1, nspec
      ! hopefully less than 10 species
      fld_vars(nspec + s) = 'enrg' // STR(s)
    end do

    ndown = 2 * nspec + 1

    fld_vars(ndown : n_fld_vars) = (/'ex   ', 'ey   ', 'ez   ',&
                                   & 'bx   ', 'by   ', 'bz   ',&
                                   & 'jx   ', 'jy   ', 'jz   ',&
                                   & 'xx   ', 'yy   ', 'zz   '/)
  end subroutine initializeSliceOutput

  #ifdef SLICE
  subroutine writeSliceX_hdf5(step, time, x_cut)
    implicit none
    integer, intent(in)               :: step, time
    integer, intent(in)               :: x_cut
    character(len=STR_MAX)            :: stepchar, filename, xchar
    integer                           :: ierr, error
    integer(HID_T)                    :: file_id, dspace_id, dset_id
    integer                           :: f, s, dset_rank = 2
    integer(HSIZE_T), dimension(2)    :: global_dims

    real, allocatable                 :: field_data(:,:)
    logical                           :: writing_lgarrQ

    integer                           :: this_x0, this_y0, this_z0, this_sx, this_sy, this_sz
    integer                           :: mb_x0, mb_y0, mb_z0, mb_sx, mb_sy, mb_sz
    integer                           :: root_rnk, rnk
    integer(kind=2)                   :: i, j, k
    real, allocatable                 :: temp_arr(:,:)

    #if defined(MPI08)
      type(MPI_STATUS)                :: istat
    #elif defined(MPI)
      integer                         :: istat(MPI_STATUS_SIZE)
    #endif

    root_rnk = 0

    if (mpi_rank .eq. root_rnk) then
      allocate(field_data(0 : global_mesh%sy - 1, 0 : global_mesh%sz - 1))

      write(xchar, "(i5.5)") x_cut
      write(stepchar, "(i5.5)") step
      filename = trim(slice_dir_name) // '/sliceX=' // trim(xchar) // '.' // trim(stepchar)

      global_dims(1) = global_mesh%sy
      global_dims(2) = global_mesh%sz

      call H5open_f(error)
      call H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      call H5Screate_simple_f(dset_rank, global_dims, dspace_id, error)
    end if

    this_x0 = this_meshblock%ptr%x0
    this_y0 = this_meshblock%ptr%y0
    this_z0 = this_meshblock%ptr%z0
    this_sx = this_meshblock%ptr%sx
    this_sy = this_meshblock%ptr%sy
    this_sz = this_meshblock%ptr%sz

    do f = 1, n_fld_vars
      call prepareFieldForOutput(fld_vars(f), writing_lgarrQ)

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if ((x_cut .ge. this_x0) .and. (x_cut .lt. this_x0 + this_sx)) then
        ! Create dataset by interpolating fields
        i = x_cut - this_x0
        do j = 0, this_sy - 1
          do k = 0, this_sz - 1
            call selectFieldForOutput(fld_vars(f), 0, j, k, i, j, k, writing_lgarrQ)
          end do
        end do
        if (mpi_rank .ne. root_rnk) then
          allocate(temp_arr(this_sy, this_sz))
          temp_arr(:,:) = sm_arr(0,:,:)
          call MPI_SEND(temp_arr(:,:), this_sy * this_sz, MPI_REAL, root_rnk, f, MPI_COMM_WORLD, ierr)
          deallocate(temp_arr)
        else
          field_data(this_y0 : this_y0 + this_sy - 1, this_z0 : this_z0 + this_sz - 1) = sm_arr(0,:,:)
        end if
      end if

      if (mpi_rank .eq. root_rnk) then
        do rnk = 0, mpi_size - 1
          mb_x0 = meshblocks(rnk + 1)%x0
          mb_y0 = meshblocks(rnk + 1)%y0
          mb_z0 = meshblocks(rnk + 1)%z0
          mb_sx = meshblocks(rnk + 1)%sx
          mb_sy = meshblocks(rnk + 1)%sy
          mb_sz = meshblocks(rnk + 1)%sz
          if ((x_cut .ge. mb_x0) .and. (x_cut .lt. mb_x0 + mb_sx) .and. (mpi_rank .ne. rnk)) then
            allocate(temp_arr(mb_sy, mb_sz))
            call MPI_RECV(temp_arr(:, :), mb_sy * mb_sz, MPI_REAL, rnk, f, MPI_COMM_WORLD, istat, ierr)
            field_data(mb_y0 : mb_y0 + mb_sy - 1, mb_z0 : mb_z0 + mb_sz - 1) = temp_arr(:,:)
            deallocate(temp_arr)
          end if
        end do

        call H5Dcreate_f(file_id, fld_vars(f), H5T_NATIVE_REAL, dspace_id, dset_id, error)
        call H5Dwrite_f(dset_id, H5T_NATIVE_REAL, field_data, global_dims, error)
        call H5Dclose_f(dset_id, error)
      end if
    end do ! fld_vars

    if (mpi_rank .eq. root_rnk) then
      call H5Sclose_f(dspace_id, error)

      call H5Fclose_f(file_id, error)
      call H5close_f(error)

      if (allocated(field_data)) deallocate(field_data)
    end if
  end subroutine writeSliceX_hdf5

  subroutine writeSliceY_hdf5(step, time, y_cut)
    implicit none
    integer, intent(in)               :: step, time
    integer, intent(in)               :: y_cut
    character(len=STR_MAX)            :: stepchar, filename, ychar
    integer                           :: ierr, error
    integer(HID_T)                    :: file_id, dspace_id, dset_id
    integer                           :: f, s, dset_rank = 2
    integer(HSIZE_T), dimension(2)    :: global_dims

    real, allocatable                 :: field_data(:,:)
    logical                           :: writing_lgarrQ

    integer                           :: this_x0, this_y0, this_z0, this_sx, this_sy, this_sz
    integer                           :: mb_x0, mb_y0, mb_z0, mb_sx, mb_sy, mb_sz
    integer                           :: root_rnk, rnk
    integer(kind=2)                   :: i, j, k
    real, allocatable                 :: temp_arr(:,:)

    #if defined(MPI08)
      type(MPI_STATUS)                :: istat
    #elif defined(MPI)
      integer                         :: istat(MPI_STATUS_SIZE)
    #endif

    root_rnk = 0

    if (mpi_rank .eq. root_rnk) then
      allocate(field_data(0 : global_mesh%sx - 1, 0 : global_mesh%sz - 1))

      write(ychar, "(i5.5)") y_cut
      write(stepchar, "(i5.5)") step
      filename = trim(slice_dir_name) // '/sliceY=' // trim(ychar) // '.' // trim(stepchar)

      global_dims(1) = global_mesh%sx
      global_dims(2) = global_mesh%sz

      call H5open_f(error)
      call H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      call H5Screate_simple_f(dset_rank, global_dims, dspace_id, error)
    end if

    this_x0 = this_meshblock%ptr%x0
    this_y0 = this_meshblock%ptr%y0
    this_z0 = this_meshblock%ptr%z0
    this_sx = this_meshblock%ptr%sx
    this_sy = this_meshblock%ptr%sy
    this_sz = this_meshblock%ptr%sz

    do f = 1, n_fld_vars
      call prepareFieldForOutput(fld_vars(f), writing_lgarrQ)

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if ((y_cut .ge. this_y0) .and. (y_cut .lt. this_y0 + this_sy)) then
        ! Create dataset by interpolating fields
        j = y_cut - this_y0
        do i = 0, this_sx - 1
          do k = 0, this_sz - 1
            call selectFieldForOutput(fld_vars(f), i, 0, k, i, j, k, writing_lgarrQ)
          end do
        end do
        if (mpi_rank .ne. root_rnk) then
          allocate(temp_arr(this_sx, this_sz))
          temp_arr(:,:) = sm_arr(:,0,:)
          call MPI_SEND(temp_arr(:,:), this_sx * this_sz, MPI_REAL, root_rnk, f, MPI_COMM_WORLD, ierr)
          deallocate(temp_arr)
        else
          field_data(this_x0 : this_x0 + this_sx - 1, this_z0 : this_z0 + this_sz - 1) = sm_arr(:,0,:)
        end if
      end if

      if (mpi_rank .eq. root_rnk) then
        do rnk = 0, mpi_size - 1
          mb_x0 = meshblocks(rnk + 1)%x0
          mb_y0 = meshblocks(rnk + 1)%y0
          mb_z0 = meshblocks(rnk + 1)%z0
          mb_sx = meshblocks(rnk + 1)%sx
          mb_sy = meshblocks(rnk + 1)%sy
          mb_sz = meshblocks(rnk + 1)%sz
          if ((y_cut .ge. mb_y0) .and. (y_cut .lt. mb_y0 + mb_sy) .and. (mpi_rank .ne. rnk)) then
            allocate(temp_arr(mb_sx, mb_sz))
            call MPI_RECV(temp_arr(:,:), mb_sx * mb_sz, MPI_REAL, rnk, f, MPI_COMM_WORLD, istat, ierr)
            field_data(mb_x0 : mb_x0 + mb_sx - 1, mb_z0 : mb_z0 + mb_sz - 1) = temp_arr(:,:)
            deallocate(temp_arr)
          end if
        end do

        call H5Dcreate_f(file_id, fld_vars(f), H5T_NATIVE_REAL, dspace_id, dset_id, error)
        call H5Dwrite_f(dset_id, H5T_NATIVE_REAL, field_data, global_dims, error)
        call H5Dclose_f(dset_id, error)
      end if
    end do ! fld_vars

    if (mpi_rank .eq. root_rnk) then
      call H5Sclose_f(dspace_id, error)

      call H5Fclose_f(file_id, error)
      call H5close_f(error)

      if (allocated(field_data)) deallocate(field_data)
    end if
  end subroutine writeSliceY_hdf5

  subroutine writeSliceZ_hdf5(step, time, z_cut)
    implicit none
    integer, intent(in)               :: step, time
    integer, intent(in)               :: z_cut
    character(len=STR_MAX)            :: stepchar, filename, zchar
    integer                           :: ierr, error
    integer(HID_T)                    :: file_id, dspace_id, dset_id
    integer                           :: f, s, dset_rank = 2
    integer(HSIZE_T), dimension(2)    :: global_dims

    real, allocatable                 :: field_data(:,:)
    logical                           :: writing_lgarrQ

    integer                           :: this_x0, this_y0, this_z0, this_sx, this_sy, this_sz
    integer                           :: mb_x0, mb_y0, mb_z0, mb_sx, mb_sy, mb_sz
    integer                           :: root_rnk, rnk
    integer(kind=2)                   :: i, j, k
    real, allocatable                 :: temp_arr(:,:)

    #if defined(MPI08)
      type(MPI_STATUS)                :: istat
    #elif defined(MPI)
      integer                         :: istat(MPI_STATUS_SIZE)
    #endif

    root_rnk = 0

    if (mpi_rank .eq. root_rnk) then
      allocate(field_data(0 : global_mesh%sx - 1, 0 : global_mesh%sy - 1))

      write(zchar, "(i5.5)") z_cut
      write(stepchar, "(i5.5)") step
      filename = trim(slice_dir_name) // '/sliceZ=' // trim(zchar) // '.' // trim(stepchar)

      global_dims(1) = global_mesh%sx
      global_dims(2) = global_mesh%sy

      call H5open_f(error)
      call H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      call H5Screate_simple_f(dset_rank, global_dims, dspace_id, error)
    end if

    this_x0 = this_meshblock%ptr%x0
    this_y0 = this_meshblock%ptr%y0
    this_z0 = this_meshblock%ptr%z0
    this_sx = this_meshblock%ptr%sx
    this_sy = this_meshblock%ptr%sy
    this_sz = this_meshblock%ptr%sz

    do f = 1, n_fld_vars
      call prepareFieldForOutput(fld_vars(f), writing_lgarrQ)

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if ((z_cut .ge. this_z0) .and. (z_cut .lt. this_z0 + this_sz)) then
        ! Create dataset by interpolating fields
        k = z_cut - this_z0
        do i = 0, this_sx - 1
          do j = 0, this_sy - 1
            call selectFieldForOutput(fld_vars(f), i, j, 0, i, j, k, writing_lgarrQ)
          end do
        end do
        if (mpi_rank .ne. root_rnk) then
          allocate(temp_arr(this_sx, this_sy))
          temp_arr(:,:) = sm_arr(:,:,0)
          call MPI_SEND(temp_arr(:,:), this_sx * this_sy, MPI_REAL, root_rnk, f, MPI_COMM_WORLD, ierr)
          deallocate(temp_arr)
        else
          field_data(this_x0 : this_x0 + this_sx - 1, this_y0 : this_y0 + this_sy - 1) = sm_arr(:,:,0)
        end if
      end if

      if (mpi_rank .eq. root_rnk) then
        do rnk = 0, mpi_size - 1
          mb_x0 = meshblocks(rnk + 1)%x0
          mb_y0 = meshblocks(rnk + 1)%y0
          mb_z0 = meshblocks(rnk + 1)%z0
          mb_sx = meshblocks(rnk + 1)%sx
          mb_sy = meshblocks(rnk + 1)%sy
          mb_sz = meshblocks(rnk + 1)%sz
          if ((z_cut .ge. mb_z0) .and. (z_cut .lt. mb_z0 + mb_sz) .and. (mpi_rank .ne. rnk)) then
            allocate(temp_arr(mb_sx, mb_sy))
            call MPI_RECV(temp_arr(:,:), mb_sx * mb_sy, MPI_REAL, rnk, f, MPI_COMM_WORLD, istat, ierr)
            field_data(mb_x0 : mb_x0 + mb_sx - 1, mb_y0 : mb_y0 + mb_sy - 1) = temp_arr(:,:)
            deallocate(temp_arr)
          end if
        end do

        call H5Dcreate_f(file_id, fld_vars(f), H5T_NATIVE_REAL, dspace_id, dset_id, error)
        call H5Dwrite_f(dset_id, H5T_NATIVE_REAL, field_data, global_dims, error)
        call H5Dclose_f(dset_id, error)
      end if
    end do ! fld_vars

    if (mpi_rank .eq. root_rnk) then
      call H5Sclose_f(dspace_id, error)

      call H5Fclose_f(file_id, error)
      call H5close_f(error)

      if (allocated(field_data)) deallocate(field_data)
    end if
  end subroutine writeSliceZ_hdf5

  #endif

end module m_writeslice
