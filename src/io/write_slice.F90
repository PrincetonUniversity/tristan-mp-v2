module m_writeslice
#ifdef HDF5
  use hdf5
#endif
  use m_globalnamespace, only: mpi_rank
  use m_outputnamespace, only: slice_index, fld_vars, n_fld_vars, &
                               nslices, slice_axes, slice_pos, xdmf_enable
  use m_aux
  use m_errors, only: throwError
  use m_domain
  use m_fields
  use m_helpers
  use m_outputlogistics, only: prepareFieldForOutput, selectFieldForOutput, &
                               defineFieldVarsToOutput
  implicit none

  !--- PRIVATE functions -----------------------------------------!
#ifdef HDF5
  private :: writeSliceX_hdf5, writeSliceY_hdf5, writeSliceZ_hdf5, writeXDMF_hdf5
#endif
  !...............................................................!

contains
  subroutine writeSlices(time)
    implicit none
    integer, intent(in) :: time
    integer :: step
    integer :: n

    call defineFieldVarsToOutput()

    step = slice_index
#ifdef HDF5
    do n = 1, nslices
      if (slice_axes(n) .eq. 1) then
        call writeSliceX_hdf5(step, time, slice_pos(n))
        call printDiag("writeSliceX_hdf5()", 3)
      else if (slice_axes(n) .eq. 2) then
        call writeSliceY_hdf5(step, time, slice_pos(n))
        call printDiag("writeSliceY_hdf5()", 3)
      else if (slice_axes(n) .eq. 3) then
        call writeSliceZ_hdf5(step, time, slice_pos(n))
        call printDiag("writeSliceZ_hdf5()", 3)
      else
        call throwError('ERROR. wrong `slice_axes(n)` in `writeSlices`')
      end if
    end do
    call printDiag("writeSlices_hdf5()", 2)
#endif
    slice_index = slice_index + 1
  end subroutine writeSlices

#ifdef HDF5
  subroutine writeXDMF_hdf5(step, fname, time, ni, nj)
    implicit none
    integer, intent(in) :: step, time, ni, nj
    character(len=*), intent(in) :: fname
    character(len=STR_MAX) :: filename
    integer :: var

    if (.false.) print *, step, time

    filename = trim(slice_dir_name)//'/'//trim(fname)//'.xdmf'

    open (UNIT_xdmf, file=filename, status="replace", access="stream", form="formatted")
    write (UNIT_xdmf, "(A)") &
      '<?xml version="1.0" ?>'
    write (UNIT_xdmf, "(A)") &
      '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
    write (UNIT_xdmf, "(A)") &
      '<Xdmf Version="2.0">'
    write (UNIT_xdmf, "(A)") &
      '  <Domain>'
    write (UNIT_xdmf, "(A)") &
      '    <Grid Name="domain" GridType="Uniform">'
    write (UNIT_xdmf, "(A,A,A,I10,I10,A)") &
      '      <Topology TopologyType=', &
      '"2DCoRectMesh"', ' Dimensions="', &
      nj, ni, '"/>'
    write (UNIT_xdmf, "(A)") &
      '      <Geometry GeometryType="Origin_DxDy">'
    write (UNIT_xdmf, "(A,A)") &
      '        <DataItem Format="XML" Dimensions="2"', &
      ' NumberType="Float" Precision="4">'
    write (UNIT_xdmf, "(A)") &
      '          0.0 0.0'
    write (UNIT_xdmf, "(A)") &
      '        </DataItem>'
    write (UNIT_xdmf, "(A,A)") &
      '        <DataItem Format="XML" Dimensions="2"', &
      ' NumberType="Float" Precision="4">'
    write (UNIT_xdmf, "(A)") &
      '          1.0 1.0'
    write (UNIT_xdmf, "(A)") &
      '        </DataItem>'
    write (UNIT_xdmf, "(A)") &
      '      </Geometry>'

    do var = 1, n_fld_vars
      write (UNIT_xdmf, "(A)") &
        '      <Attribute Name="'//trim(fld_vars(var))//'" Center="Node">'
      write (UNIT_xdmf, "(A,I10,I10,A)") &
        '        <DataItem Format="HDF" Dimensions="', &
        nj, ni, &
        '" NumberType="Float" Precision="4">'
      write (UNIT_xdmf, "(A,A)") &
        trim(fname)//':/', &
        trim(fld_vars(var))
      write (UNIT_xdmf, "(A)") &
        '        </DataItem>'
      write (UNIT_xdmf, "(A)") &
        '      </Attribute>'
    end do

    write (UNIT_xdmf, "(A)") &
      '    </Grid>'
    write (UNIT_xdmf, "(A)") &
      '  </Domain>'
    write (UNIT_xdmf, "(A)") &
      '</Xdmf>'
    close (UNIT_xdmf)
  end subroutine writeXDMF_hdf5

  subroutine writeSliceX_hdf5(step, time, x_cut)
    implicit none
    integer, intent(in) :: step, time
    integer, intent(in) :: x_cut
    character(len=STR_MAX) :: stepchar, filename, xchar
    integer :: ierr, error
    integer(HID_T) :: file_id, dspace_id, dset_id
    integer :: f, dset_rank = 2
    integer(HSIZE_T), dimension(2) :: global_dims

    real, allocatable :: field_data(:, :)
    logical :: writing_lgarrQ

    integer :: this_x0, this_y0, this_z0, this_sx, this_sy, this_sz
    integer :: mb_x0, mb_y0, mb_z0, mb_sx, mb_sy, mb_sz
    integer :: root_rnk, rnk
    integer(kind=2) :: i, j, k
    real, allocatable :: temp_arr(:, :)

#if defined(MPI08)
    type(MPI_STATUS) :: istat
#elif defined(MPI)
    integer :: istat(MPI_STATUS_SIZE)
#endif

    root_rnk = 0

    if (mpi_rank .eq. root_rnk) then
      allocate (field_data(0:global_mesh % sy - 1, 0:global_mesh % sz - 1))

      write (xchar, "(i5.5)") x_cut
      write (stepchar, "(i5.5)") step
      filename = trim(slice_dir_name)//'/sliceX='//trim(xchar)//'.'//trim(stepchar)

      global_dims(1) = global_mesh % sy
      global_dims(2) = global_mesh % sz

      if (xdmf_enable) then
        call writeXDMF_hdf5(step, 'sliceX='//trim(xchar)//'.'//trim(stepchar), &
                            time, INT(global_dims(1)), INT(global_dims(2)))
      end if

      call H5open_f(error)
      call H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      call H5Screate_simple_f(dset_rank, global_dims, dspace_id, error)
    else
      allocate(field_data(0,0))
    end if

    this_x0 = this_meshblock % ptr % x0
    this_y0 = this_meshblock % ptr % y0
    this_z0 = this_meshblock % ptr % z0
    this_sx = this_meshblock % ptr % sx
    this_sy = this_meshblock % ptr % sy
    this_sz = this_meshblock % ptr % sz

    do f = 1, n_fld_vars
      call prepareFieldForOutput(fld_vars(f), writing_lgarrQ)

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if ((x_cut .ge. this_x0) .and. (x_cut .lt. this_x0 + this_sx)) then
        ! Create dataset by interpolating fields
        i = INT(x_cut - this_x0, 2)
        do j = 0, INT(this_sy - 1, 2)
          do k = 0, INT(this_sz - 1, 2)
            call selectFieldForOutput(fld_vars(f), 0_2, j, k, i, j, k, writing_lgarrQ)
          end do
        end do
        if (mpi_rank .ne. root_rnk) then
          allocate (temp_arr(this_sy, this_sz))
          temp_arr(:, :) = sm_arr(0, :, :)
          call MPI_SEND(temp_arr(:, :), this_sy * this_sz, default_mpi_real, root_rnk, f, MPI_COMM_WORLD, ierr)
          deallocate (temp_arr)
        else
          field_data(this_y0:this_y0 + this_sy - 1, this_z0:this_z0 + this_sz - 1) = sm_arr(0, :, :)
        end if
      end if

      if (mpi_rank .eq. root_rnk) then
        do rnk = 0, mpi_size - 1
          mb_x0 = meshblocks(rnk + 1) % x0
          mb_y0 = meshblocks(rnk + 1) % y0
          mb_z0 = meshblocks(rnk + 1) % z0
          mb_sx = meshblocks(rnk + 1) % sx
          mb_sy = meshblocks(rnk + 1) % sy
          mb_sz = meshblocks(rnk + 1) % sz
          if ((x_cut .ge. mb_x0) .and. (x_cut .lt. mb_x0 + mb_sx) .and. (mpi_rank .ne. rnk)) then
            allocate (temp_arr(mb_sy, mb_sz))
            call MPI_RECV(temp_arr(:, :), mb_sy * mb_sz, default_mpi_real, rnk, f, MPI_COMM_WORLD, istat, ierr)
            field_data(mb_y0:mb_y0 + mb_sy - 1, mb_z0:mb_z0 + mb_sz - 1) = temp_arr(:, :)
            deallocate (temp_arr)
          end if
        end do

        call H5Dcreate_f(file_id, fld_vars(f), default_h5_real, dspace_id, dset_id, error)
        call H5Dwrite_f(dset_id, default_h5_real, field_data, global_dims, error)
        call H5Dclose_f(dset_id, error)
      end if
    end do ! fld_vars

    if (mpi_rank .eq. root_rnk) then
      call H5Sclose_f(dspace_id, error)

      call H5Fclose_f(file_id, error)
      call H5close_f(error)

      if (allocated(field_data)) deallocate (field_data)
    end if
  end subroutine writeSliceX_hdf5

  subroutine writeSliceY_hdf5(step, time, y_cut)
    implicit none
    integer, intent(in) :: step, time
    integer, intent(in) :: y_cut
    character(len=STR_MAX) :: stepchar, filename, ychar
    integer :: ierr, error
    integer(HID_T) :: file_id, dspace_id, dset_id
    integer :: f, dset_rank = 2
    integer(HSIZE_T), dimension(2) :: global_dims

    real, allocatable :: field_data(:, :)
    logical :: writing_lgarrQ

    integer :: this_x0, this_y0, this_z0, this_sx, this_sy, this_sz
    integer :: mb_x0, mb_y0, mb_z0, mb_sx, mb_sy, mb_sz
    integer :: root_rnk, rnk
    integer(kind=2) :: i, j, k
    real, allocatable :: temp_arr(:, :)

#if defined(MPI08)
    type(MPI_STATUS) :: istat
#elif defined(MPI)
    integer :: istat(MPI_STATUS_SIZE)
#endif

    root_rnk = 0

    if (mpi_rank .eq. root_rnk) then
      allocate (field_data(0:global_mesh % sx - 1, 0:global_mesh % sz - 1))

      write (ychar, "(i5.5)") y_cut
      write (stepchar, "(i5.5)") step
      filename = trim(slice_dir_name)//'/sliceY='//trim(ychar)//'.'//trim(stepchar)

      global_dims(1) = global_mesh % sx
      global_dims(2) = global_mesh % sz

      if (xdmf_enable) then
        call writeXDMF_hdf5(step, 'sliceY='//trim(ychar)//'.'//trim(stepchar), &
                            time, INT(global_dims(1)), INT(global_dims(2)))
      end if

      call H5open_f(error)
      call H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      call H5Screate_simple_f(dset_rank, global_dims, dspace_id, error)
    else
      allocate(field_data(0, 0))
    end if

    this_x0 = this_meshblock % ptr % x0
    this_y0 = this_meshblock % ptr % y0
    this_z0 = this_meshblock % ptr % z0
    this_sx = this_meshblock % ptr % sx
    this_sy = this_meshblock % ptr % sy
    this_sz = this_meshblock % ptr % sz

    do f = 1, n_fld_vars
      call prepareFieldForOutput(fld_vars(f), writing_lgarrQ)

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if ((y_cut .ge. this_y0) .and. (y_cut .lt. this_y0 + this_sy)) then
        ! Create dataset by interpolating fields
        j = INT(y_cut - this_y0, 2)
        do i = 0, INT(this_sx - 1, 2)
          do k = 0, INT(this_sz - 1, 2)
            call selectFieldForOutput(fld_vars(f), i, 0_2, k, i, j, k, writing_lgarrQ)
          end do
        end do
        if (mpi_rank .ne. root_rnk) then
          allocate (temp_arr(this_sx, this_sz))
          temp_arr(:, :) = sm_arr(:, 0, :)
          call MPI_SEND(temp_arr(:, :), this_sx * this_sz, default_mpi_real, root_rnk, f, MPI_COMM_WORLD, ierr)
          deallocate (temp_arr)
        else
          field_data(this_x0:this_x0 + this_sx - 1, this_z0:this_z0 + this_sz - 1) = sm_arr(:, 0, :)
        end if
      end if

      if (mpi_rank .eq. root_rnk) then
        do rnk = 0, mpi_size - 1
          mb_x0 = meshblocks(rnk + 1) % x0
          mb_y0 = meshblocks(rnk + 1) % y0
          mb_z0 = meshblocks(rnk + 1) % z0
          mb_sx = meshblocks(rnk + 1) % sx
          mb_sy = meshblocks(rnk + 1) % sy
          mb_sz = meshblocks(rnk + 1) % sz
          if ((y_cut .ge. mb_y0) .and. (y_cut .lt. mb_y0 + mb_sy) .and. (mpi_rank .ne. rnk)) then
            allocate (temp_arr(mb_sx, mb_sz))
            call MPI_RECV(temp_arr(:, :), mb_sx * mb_sz, default_mpi_real, rnk, f, MPI_COMM_WORLD, istat, ierr)
            field_data(mb_x0:mb_x0 + mb_sx - 1, mb_z0:mb_z0 + mb_sz - 1) = temp_arr(:, :)
            deallocate (temp_arr)
          end if
        end do

        call H5Dcreate_f(file_id, fld_vars(f), default_h5_real, dspace_id, dset_id, error)
        call H5Dwrite_f(dset_id, default_h5_real, field_data, global_dims, error)
        call H5Dclose_f(dset_id, error)
      end if
    end do ! fld_vars

    if (mpi_rank .eq. root_rnk) then
      call H5Sclose_f(dspace_id, error)

      call H5Fclose_f(file_id, error)
      call H5close_f(error)

      if (allocated(field_data)) deallocate (field_data)
    end if
  end subroutine writeSliceY_hdf5

  subroutine writeSliceZ_hdf5(step, time, z_cut)
    implicit none
    integer, intent(in) :: step, time
    integer, intent(in) :: z_cut
    character(len=STR_MAX) :: stepchar, filename, zchar
    integer :: ierr, error
    integer(HID_T) :: file_id, dspace_id, dset_id
    integer :: f, dset_rank = 2
    integer(HSIZE_T), dimension(2) :: global_dims

    real, allocatable :: field_data(:, :)
    logical :: writing_lgarrQ

    integer :: this_x0, this_y0, this_z0, this_sx, this_sy, this_sz
    integer :: mb_x0, mb_y0, mb_z0, mb_sx, mb_sy, mb_sz
    integer :: root_rnk, rnk
    integer(kind=2) :: i, j, k
    real, allocatable :: temp_arr(:, :)

#if defined(MPI08)
    type(MPI_STATUS) :: istat
#elif defined(MPI)
    integer :: istat(MPI_STATUS_SIZE)
#endif

    root_rnk = 0

    if (mpi_rank .eq. root_rnk) then
      allocate (field_data(0:global_mesh % sx - 1, 0:global_mesh % sy - 1))

      write (zchar, "(i5.5)") z_cut
      write (stepchar, "(i5.5)") step
      filename = trim(slice_dir_name)//'/sliceZ='//trim(zchar)//'.'//trim(stepchar)

      global_dims(1) = global_mesh % sx
      global_dims(2) = global_mesh % sy

      if (xdmf_enable) then
        call writeXDMF_hdf5(step, 'sliceZ='//trim(zchar)//'.'//trim(stepchar), &
                            time, INT(global_dims(1)), INT(global_dims(2)))
      end if

      call H5open_f(error)
      call H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      call H5Screate_simple_f(dset_rank, global_dims, dspace_id, error)
    else
      allocate(field_data(0,0))
    end if

    this_x0 = this_meshblock % ptr % x0
    this_y0 = this_meshblock % ptr % y0
    this_z0 = this_meshblock % ptr % z0
    this_sx = this_meshblock % ptr % sx
    this_sy = this_meshblock % ptr % sy
    this_sz = this_meshblock % ptr % sz

    do f = 1, n_fld_vars
      call prepareFieldForOutput(fld_vars(f), writing_lgarrQ)

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if ((z_cut .ge. this_z0) .and. (z_cut .lt. this_z0 + this_sz)) then
        ! Create dataset by interpolating fields
        k = INT(z_cut - this_z0, 2)
        do i = 0, INT(this_sx - 1, 2)
          do j = 0, INT(this_sy - 1, 2)
            call selectFieldForOutput(fld_vars(f), i, j, 0_2, i, j, k, writing_lgarrQ)
          end do
        end do
        if (mpi_rank .ne. root_rnk) then
          allocate (temp_arr(this_sx, this_sy))
          temp_arr(:, :) = sm_arr(:, :, 0)
          call MPI_SEND(temp_arr(:, :), this_sx * this_sy, default_mpi_real, root_rnk, f, MPI_COMM_WORLD, ierr)
          deallocate (temp_arr)
        else
          field_data(this_x0:this_x0 + this_sx - 1, this_y0:this_y0 + this_sy - 1) = sm_arr(:, :, 0)
        end if
      end if

      if (mpi_rank .eq. root_rnk) then
        do rnk = 0, mpi_size - 1
          mb_x0 = meshblocks(rnk + 1) % x0
          mb_y0 = meshblocks(rnk + 1) % y0
          mb_z0 = meshblocks(rnk + 1) % z0
          mb_sx = meshblocks(rnk + 1) % sx
          mb_sy = meshblocks(rnk + 1) % sy
          mb_sz = meshblocks(rnk + 1) % sz
          if ((z_cut .ge. mb_z0) .and. (z_cut .lt. mb_z0 + mb_sz) .and. (mpi_rank .ne. rnk)) then
            allocate (temp_arr(mb_sx, mb_sy))
            call MPI_RECV(temp_arr(:, :), mb_sx * mb_sy, default_mpi_real, rnk, f, MPI_COMM_WORLD, istat, ierr)
            field_data(mb_x0:mb_x0 + mb_sx - 1, mb_y0:mb_y0 + mb_sy - 1) = temp_arr(:, :)
            deallocate (temp_arr)
          end if
        end do

        call H5Dcreate_f(file_id, fld_vars(f), default_h5_real, dspace_id, dset_id, error)
        call H5Dwrite_f(dset_id, default_h5_real, field_data, global_dims, error)
        call H5Dclose_f(dset_id, error)
      end if
    end do ! fld_vars

    if (mpi_rank .eq. root_rnk) then
      call H5Sclose_f(dspace_id, error)

      call H5Fclose_f(file_id, error)
      call H5close_f(error)

      if (allocated(field_data)) deallocate (field_data)
    end if
  end subroutine writeSliceZ_hdf5
#endif

end module m_writeslice
