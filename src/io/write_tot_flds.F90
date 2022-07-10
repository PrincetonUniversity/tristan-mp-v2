module m_writetotflds
#ifdef HDF5
  use hdf5
  use m_globalnamespace, only: h5comm, h5info
#endif
  use m_globalnamespace, only: output_dir_name, mpi_rank
  use m_outputnamespace, only: xdmf_enable, fld_vars, n_fld_vars, output_flds_istep
  use m_aux
  use m_domain
  use m_fields
  use m_outputlogistics, only: prepareFieldForOutput, selectFieldForOutput
  implicit none

  !--- PRIVATE functions -----------------------------------------!
#ifdef HDF5
  private :: writeFields_hdf5, writeXDMF_hdf5, getBlockDimensions
#endif
  !...............................................................!
contains
  subroutine writeFields(step, time)
    implicit none
    integer, intent(in) :: step, time
#ifdef HDF5
    call writeFields_hdf5(step, time)
#endif
    call printDiag("writeFields()", 3)
  end subroutine writeFields

#ifdef HDF5
  subroutine writeXDMF_hdf5(step, time, ni, nj, nk)
    implicit none
    integer, intent(in) :: step, time, ni, nj, nk
    character(len=STR_MAX) :: stepchar, filename
    integer :: var

    if (.false.) print *, time

    write (stepchar, "(i5.5)") step
    filename = trim(output_dir_name)//'/flds.tot.'//trim(stepchar)//'.xdmf'

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
    write (UNIT_xdmf, "(A,A,A,I10,I10,I10,A)") &
      '      <Topology TopologyType=', &
      '"3DCoRectMesh"', ' Dimensions="', &
      nk, nj, ni, '"/>'
    write (UNIT_xdmf, "(A)") &
      '      <Geometry GeometryType="ORIGIN_DXDYDZ">'
    write (UNIT_xdmf, "(A,A)") &
      '        <DataItem Format="XML" Dimensions="3"', &
      ' NumberType="Float" Precision="4">'
    write (UNIT_xdmf, "(A)") &
      '          0.0 0.0 0.0'
    write (UNIT_xdmf, "(A)") &
      '        </DataItem>'
    write (UNIT_xdmf, "(A,A)") &
      '        <DataItem Format="XML" Dimensions="3"', &
      ' NumberType="Float" Precision="4">'
    write (UNIT_xdmf, "(A)") &
      '          1.0 1.0 1.0'
    write (UNIT_xdmf, "(A)") &
      '        </DataItem>'
    write (UNIT_xdmf, "(A)") &
      '      </Geometry>'

    do var = 1, n_fld_vars
      write (UNIT_xdmf, "(A)") &
        '      <Attribute Name="'//trim(fld_vars(var))//'" Center="Node">'
      write (UNIT_xdmf, "(A,I10,I10,I10,A)") &
        '        <DataItem Format="HDF" Dimensions="', &
        nk, nj, ni, &
        '" NumberType="Float" Precision="4">'
      write (UNIT_xdmf, "(A,A)") &
        '          flds.tot.'//trim(stepchar)//':/', &
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

#ifndef SERIALOUTPUT
  ! parallel field output
  subroutine writeFields_hdf5(step, time)
    implicit none
    integer, intent(in) :: step, time
    character(len=STR_MAX) :: stepchar, filename
    integer(HID_T) :: file_id, dset_id(100), filespace(100), memspace(100), plist_id
    integer :: error, f
    integer(kind=2) :: i, j, k
    logical :: writing_lgarrQ
    integer :: dataset_rank = 3
    integer(HSSIZE_T), dimension(3) :: offsets
    integer(HSIZE_T), dimension(3) :: global_dims, blocks
    integer, dimension(3) :: starts
    integer(kind=2) :: i1, j1, k1

    call getBlockDimensions(this_meshblock % ptr, starts, offsets, blocks, global_dims)

    write (stepchar, "(i5.5)") step
    filename = trim(output_dir_name)//'/flds.tot.'//trim(stepchar)

    if ((mpi_rank .eq. 0) .and. xdmf_enable) then
      call writeXDMF_hdf5(step, time, INT(global_dims(1)), INT(global_dims(2)), INT(global_dims(3)))
    end if

    ! Initialize HDF5 library and Fortran interfaces
    call h5open_f(error)

    ! Setup file access property list with parallel I/O access
    call h5Pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5Pset_fapl_mpio_f(plist_id, h5comm, h5info, error)

    ! Create the file collectively
    call h5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp=plist_id)
    call h5Pclose_f(plist_id, error)

    do f = 1, n_fld_vars
      call h5Screate_simple_f(dataset_rank, global_dims, filespace(f), error)
    end do
    do f = 1, n_fld_vars
      call h5Dcreate_f(file_id, fld_vars(f), default_h5_real, filespace(f), dset_id(f), error)
    end do
    do f = 1, n_fld_vars
      call h5Sselect_hyperslab_f(filespace(f), H5S_SELECT_SET_F, offsets, blocks, error)
    end do
    do f = 1, n_fld_vars
      call h5Screate_simple_f(dataset_rank, blocks, memspace(f), error)
    end do

    call h5Pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

    do f = 1, n_fld_vars
      call prepareFieldForOutput(fld_vars(f), writing_lgarrQ)
      ! Create dataset by interpolating fields
      do i1 = 0, INT(blocks(1) - 1, 2)
        do j1 = 0, INT(blocks(2) - 1, 2)
          do k1 = 0, INT(blocks(3) - 1, 2)
            i = INT(starts(1) + i1 * output_flds_istep, 2)
            j = INT(starts(2) + j1 * output_flds_istep, 2)
            k = INT(starts(3) + k1 * output_flds_istep, 2)
            call selectFieldForOutput(fld_vars(f), i1, j1, k1, i, j, k, writing_lgarrQ)
          end do
        end do
      end do

      ! Write the dataset collectively
      call h5Dwrite_f(dset_id(f), default_h5_real, sm_arr(0:blocks(1) - 1, 0:blocks(2) - 1, 0:blocks(3) - 1), &
                      global_dims, error, file_space_id=filespace(f), mem_space_id=memspace(f), xfer_prp=plist_id)
    end do

    call h5Pclose_f(plist_id, error)
    do f = 1, n_fld_vars
      call h5Sclose_f(filespace(f), error)
    end do
    do f = 1, n_fld_vars
      call h5Sclose_f(memspace(f), error)
    end do
    do f = 1, n_fld_vars
      call h5Dclose_f(dset_id(f), error)
    end do

    ! Close the file
    call h5Fclose_f(file_id, error)
    ! Close FORTRAN interfaces and HDF5 library
    call h5close_f(error)
  end subroutine writeFields_hdf5
#else
  ! serial field output
  subroutine writeFields_hdf5(step, time)
    implicit none
    integer, intent(in) :: step, time
    character(len=STR_MAX) :: stepchar, filename
    integer(HID_T) :: file_id, dset_id(100), filespace(100), memspace(100)
    integer :: error, f, s, root_rnk, rnk, ierr, max_sx, max_sy, max_sz
    integer(kind=2) :: i, j, k
    logical :: writing_intQ, writing_lgarrQ
    integer :: dataset_rank = 3
    integer(HSSIZE_T), dimension(3) :: offsets, rnk_offsets
    integer(HSIZE_T), dimension(3) :: global_dims, blocks, rnk_blocks
    integer, dimension(3) :: starts, rnk_starts
    real :: ex0, ey0, ez0, bx0, by0, bz0
    real :: jx0, jy0, jz0
    integer(kind=2) :: i1, j1, k1

#if defined(MPI08)
    type(MPI_STATUS) :: istat
#elif defined(MPI)
    integer :: istat(MPI_STATUS_SIZE)
#endif

    ! determine root rank
#if defined(SLB) || defined(ALB)
    max_sx = 0; max_sy = 0; max_sz = 0
    do rnk = 0, mpi_size - 1
      if (meshblocks(rnk + 1) % sx .gt. max_sx) max_sx = meshblocks(rnk + 1) % sx
      if (meshblocks(rnk + 1) % sy .gt. max_sy) max_sy = meshblocks(rnk + 1) % sy
      if (meshblocks(rnk + 1) % sz .gt. max_sz) max_sz = meshblocks(rnk + 1) % sz
    end do
    do rnk = 0, mpi_size - 1
      if ((meshblocks(rnk + 1) % sx .eq. max_sx) .and. &
          (meshblocks(rnk + 1) % sz .eq. max_sy) .and. &
          (meshblocks(rnk + 1) % sy .eq. max_sz)) then
        root_rnk = rnk
        exit
      end if
    end do
#else
    root_rnk = 0
#endif

    call getBlockDimensions(this_meshblock % ptr, starts, offsets, blocks, global_dims)

    if (mpi_rank .eq. root_rnk) then
      write (stepchar, "(i5.5)") step
      filename = trim(output_dir_name)//'/flds.tot.'//trim(stepchar)

      if (xdmf_enable) then
        call writeXDMF_hdf5(step, time, INT(global_dims(1)), INT(global_dims(2)), INT(global_dims(3)))
      end if

      call h5open_f(error)
      call h5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, H5P_DEFAULT_F, H5P_DEFAULT_F)
      do f = 1, n_fld_vars
        call h5Screate_simple_f(dataset_rank, global_dims, filespace(f), error)
      end do
    end if

    do f = 1, n_fld_vars
      if (mpi_rank .eq. root_rnk) then
        call h5Dcreate_f(file_id, fld_vars(f), default_h5_real, filespace(f), dset_id(f), error)
        call h5Dclose_f(dset_id(f), error)
      end if

      ! receive and write (or write local)
      do rnk = 0, mpi_size - 1
        if (rnk .eq. mpi_rank) then
          ! Prepare fields to output
          call prepareFieldForOutput(fld_vars(f), writing_lgarrQ)
          do i1 = 0, blocks(1) - 1
            do j1 = 0, blocks(2) - 1
              do k1 = 0, blocks(3) - 1
                i = starts(1) + i1 * output_flds_istep
                j = starts(2) + j1 * output_flds_istep
                k = starts(3) + k1 * output_flds_istep
                call selectFieldForOutput(fld_vars(f), i1, j1, k1, i, j, k, writing_lgarrQ)
              end do
            end do
          end do
          ! now the field is in `sm_arr(0 : blocks(1) - 1, 0 : blocks(2) - 1, 0 : blocks(3) - 1)`
        end if

        if (mpi_rank .eq. root_rnk) then ! root rank does...
          ! get the array to write
          call getBlockDimensions(meshblocks(rnk + 1), rnk_starts, rnk_offsets, rnk_blocks, global_dims)
          if (rnk .ne. root_rnk) then
            call MPI_RECV(sm_arr(0:rnk_blocks(1) - 1, 0:rnk_blocks(2) - 1, 0:rnk_blocks(3) - 1), &
                          INT(rnk_blocks(1) * rnk_blocks(2) * rnk_blocks(3)), default_mpi_real, &
                          rnk, f, MPI_COMM_WORLD, istat, ierr)
          end if

          call h5Sselect_hyperslab_f(filespace(f), H5S_SELECT_SET_F, rnk_offsets, rnk_blocks, error)
          call h5Screate_simple_f(dataset_rank, blocks, memspace(f), error)
          call h5Dopen_f(file_id, fld_vars(f), dset_id(f), error)

          call h5Dwrite_f(dset_id(f), default_h5_real, sm_arr(0:rnk_blocks(1) - 1, 0:rnk_blocks(2) - 1, 0:rnk_blocks(3) - 1), &
                          global_dims, error, file_space_id=filespace(f), mem_space_id=memspace(f))
          call h5Dclose_f(dset_id(f), error)
          call h5Sclose_f(memspace(f), error)
        else if (rnk .eq. mpi_rank) then ! other ranks do...
          ! send to `root_rnk`
          call MPI_SEND(sm_arr(0:blocks(1) - 1, 0:blocks(2) - 1, 0:blocks(3) - 1), &
                        INT(blocks(1) * blocks(2) * blocks(3)), default_mpi_real, &
                        root_rnk, f, MPI_COMM_WORLD, ierr)
        end if
      end do ! ranks

      if (mpi_rank .eq. root_rnk) then
        call h5Sclose_f(filespace(f), error)
      end if
    end do ! field variables

    if (mpi_rank .eq. root_rnk) then
      call h5Fclose_f(file_id, error)
      call h5close_f(error)
    end if
  end subroutine writeFields_hdf5
#endif
#endif

#ifdef HDF5
  subroutine getBlockDimensions(meshblock, starts, offsets, blocks, global_dims)
    implicit none
    type(mesh), intent(in) :: meshblock
    integer(HSSIZE_T), intent(out), dimension(3) :: offsets
    integer(HSIZE_T), intent(out), dimension(3) :: global_dims, blocks
    integer, intent(out), dimension(3) :: starts
    integer :: this_x0, this_y0, this_z0, this_sx, this_sy, this_sz
    integer :: i_start, i_end, j_start, j_end, k_start, k_end
    integer :: offset_i, offset_j, offset_k
    integer :: n_i, n_j, n_k, glob_n_i, glob_n_j, glob_n_k

    ! assuming `global_mesh%{x0,y0,z0} .eq. 0`
    this_x0 = meshblock % x0
    this_y0 = meshblock % y0
    this_z0 = meshblock % z0
    this_sx = meshblock % sx
    this_sy = meshblock % sy
    this_sz = meshblock % sz

    if (output_flds_istep .eq. 1) then
      offset_i = this_x0; offset_j = this_y0; offset_k = this_z0
      n_i = this_sx - 1; n_j = this_sy - 1; n_k = this_sz - 1
      glob_n_i = global_mesh % sx
      glob_n_j = global_mesh % sy
      glob_n_k = global_mesh % sz

      i_start = 0; j_start = 0; k_start = 0
    else
      i_start = 0; i_end = 0
      offset_i = 0; n_i = 0
      glob_n_i = 1

      j_start = 0; j_end = 0
      offset_j = 0; n_j = 0
      glob_n_j = 1

      k_start = 0; k_end = 0
      offset_k = 0; n_k = 0
      glob_n_k = 1
#if defined(oneD) || defined (twoD) || defined (threeD)
      offset_i = CEILING(REAL(this_x0) / REAL(output_flds_istep))
      i_start = CEILING(REAL(this_x0) / REAL(output_flds_istep)) * output_flds_istep - this_x0
      i_end = (CEILING(REAL(this_x0 + this_sx) / REAL(output_flds_istep)) - 1) * output_flds_istep - this_x0
      n_i = (i_end - i_start) / output_flds_istep
      glob_n_i = CEILING(REAL(global_mesh % sx) / REAL(output_flds_istep))
      glob_n_i = MAX(1, glob_n_i)
#endif
#if defined(twoD) || defined (threeD)
      offset_j = CEILING(REAL(this_y0) / REAL(output_flds_istep))
      j_start = CEILING(REAL(this_y0) / REAL(output_flds_istep)) * output_flds_istep - this_y0
      j_end = (CEILING(REAL(this_y0 + this_sy) / REAL(output_flds_istep)) - 1) * output_flds_istep - this_y0
      n_j = (j_end - j_start) / output_flds_istep
      glob_n_j = CEILING(REAL(global_mesh % sy) / REAL(output_flds_istep))
      glob_n_j = MAX(1, glob_n_j)
#endif
#if defined(threeD)
      offset_k = CEILING(REAL(this_z0) / REAL(output_flds_istep))
      k_start = CEILING(REAL(this_z0) / REAL(output_flds_istep)) * output_flds_istep - this_z0
      k_end = (CEILING(REAL(this_z0 + this_sz) / REAL(output_flds_istep)) - 1) * output_flds_istep - this_z0
      n_k = (k_end - k_start) / output_flds_istep
      glob_n_k = CEILING(REAL(global_mesh % sz) / REAL(output_flds_istep))
      glob_n_k = MAX(1, glob_n_k)
#endif
    end if

    starts(1) = i_start
    starts(2) = j_start
    starts(3) = k_start
    offsets(1) = offset_i
    offsets(2) = offset_j
    offsets(3) = offset_k
    blocks(1) = n_i + 1
    blocks(2) = n_j + 1
    blocks(3) = n_k + 1
    global_dims(1) = glob_n_i
    global_dims(2) = glob_n_j
    global_dims(3) = glob_n_k
  end subroutine getBlockDimensions
#endif

end module m_writetotflds
