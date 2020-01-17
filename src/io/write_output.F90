#include "../defs.F90"

module m_writeoutput
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
  use m_exchangearray

  ! extra physics
  #ifdef RADIATION
    use m_radiation
  #endif

  #ifdef BWPAIRPRODUCTION
    use m_bwpairproduction
  #endif

  implicit none

  integer                 :: output_start, output_interval, output_stride, output_istep
  integer                 :: n_fld_vars, n_prtl_vars, n_dom_vars
  character(len=STR_MAX)  :: prtl_vars(100), prtl_var_types(100), fld_vars(100), dom_vars(100)
  integer, allocatable, dimension(:,:) :: glob_spectra
  logical                 :: flds_at_prtl


  !--- PRIVATE functions -----------------------------------------!
  #ifndef HDF5
    private :: writeParticles_Binary, writeFields_Binary,&
             & writeSpectra_Binary
  #else
    private :: writeParticles_hdf5, writeFields_hdf5,&
             & writeSpectra_hdf5, writeDomain_hdf5,&
             & writeXDMF_hdf5
  #endif
  private :: initializeOutput
  !...............................................................!

  !--- PRIVATE variables -----------------------------------------!
  private :: n_fld_vars, n_prtl_vars
  private :: prtl_vars, prtl_var_types
  private :: fld_vars
  !...............................................................!
contains
  subroutine writeOutput(time)
    implicit none
    integer, intent(in)        :: time
    integer                    :: step, ierr

    call initializeOutput()

    step = output_index
    #ifdef HDF5
      call writeParticles_hdf5(step, time)
        call printDiag((mpi_rank .eq. 0), "...writeParticles_hdf5()", .true.)
      call writeFields_hdf5(step, time)
        call printDiag((mpi_rank .eq. 0), "...writeFields_hdf5()", .true.)
      call writeSpectra_hdf5(step, time)
        call printDiag((mpi_rank .eq. 0), "...writeSpectra_hdf5()", .true.)
      call writeDomain_hdf5(step, time)
        call printDiag((mpi_rank .eq. 0), "...writeDomain_hdf5()", .true.)
    #endif
    call printDiag((mpi_rank .eq. 0), "output()", .true.)
    output_index = output_index + 1
  end subroutine writeOutput

  subroutine initializeOutput()
    ! DEP_PRT [particle-dependent]
    implicit none
    real                      :: energy, u_, v_, w_
    integer                   :: s, i, ti, tj, tk, p, spec_index
    integer                   :: ierr
    integer, allocatable, dimension(:,:)  :: spectra
    integer, allocatable, dimension(:)    :: send_spec_int, recv_spec_int
    real, allocatable, dimension(:)       :: send_spec_real, recv_spec_real
    ! initialize particle variables
    if (.not. flds_at_prtl) then
      n_prtl_vars = 8
      prtl_vars(1:n_prtl_vars) = (/'x    ', 'y    ', 'z    ',&
                                 & 'u    ', 'v    ', 'w    ',&
                                 & 'ind  ', 'proc '/)
      prtl_var_types(1:n_prtl_vars) = (/'real ', 'real ', 'real ',&
                                      & 'real ', 'real ', 'real ',&
                                      & 'int  ', 'int  '/)
    else
      n_prtl_vars = 14
      prtl_vars(1:n_prtl_vars) = (/'x    ', 'y    ', 'z    ',&
                                 & 'u    ', 'v    ', 'w    ',&
                                 & 'ind  ', 'proc ',&
                                 & 'ex   ', 'ey   ', 'ez   ',&
                                 & 'bx   ', 'by   ', 'bz   '/)
      prtl_var_types(1:n_prtl_vars) = (/'real ', 'real ', 'real ',&
                                      & 'real ', 'real ', 'real ',&
                                      & 'int  ', 'int  ',&
                                      & 'real ', 'real ', 'real ',&
                                      & 'real ', 'real ', 'real '/)
      do s = 1, nspec
        prtl_vars(n_prtl_vars + s) = 'dens' // STR(s)
        prtl_var_types(n_prtl_vars + s) = 'real '
      end do
      n_prtl_vars = n_prtl_vars + nspec
    end if

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
    fld_vars(2 * nspec + 1 : n_fld_vars) = (/'ex   ', 'ey   ', 'ez   ',&
                                           & 'bx   ', 'by   ', 'bz   ',&
                                           & 'jx   ', 'jy   ', 'jz   ',&
                                           & 'xx   ', 'yy   ', 'zz   '/)

    ! initialize domain output variables
    !   FIX1: maybe add # of particles per domain
    n_dom_vars = 6
    dom_vars(1 : 6) = (/'x0   ', 'y0   ', 'z0   ',&
                      & 'sx   ', 'sy   ', 'sz   '/)

    ! compute spectra
    if (.not. allocated(glob_spectra)) then
      allocate(glob_spectra(nspec, spec_num))
    end if
    allocate(spectra(nspec, spec_num))
    allocate(send_spec_int(spec_num), recv_spec_int(spec_num))
    allocate(send_spec_real(spec_num), recv_spec_real(spec_num))

    spectra(:,:) = 0
    do s = 1, nspec
     do ti = 1, species(s)%tile_nx
       do tj = 1, species(s)%tile_ny
         do tk = 1, species(s)%tile_nz
           ! !$omp simd
           ! !dir$ vector aligned
           do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
             u_ = species(s)%prtl_tile(ti, tj, tk)%u(p)
             v_ = species(s)%prtl_tile(ti, tj, tk)%v(p)
             w_ = species(s)%prtl_tile(ti, tj, tk)%w(p)
             if ((species(s)%m_sp .eq. 0) .and. (species(s)%ch_sp .eq. 0)) then
               energy = sqrt(u_**2 + v_**2 + w_**2)
             else
               energy = sqrt(1.0 + u_**2 + v_**2 + w_**2) - 1.0
             end if
             energy = log(energy)
             if (energy .le. spec_min) then
               spec_index = 1
             else if (energy .ge. spec_max) then
               spec_index = spec_num
             else
               spec_index = INT(CEILING((energy - spec_min) * REAL(spec_num) / (spec_max - spec_min)))
               if (spec_index .lt. 1) spec_index = 1
               if (spec_index .gt. spec_num) spec_index = spec_num
             end if
             spectra(s, spec_index) = spectra(s, spec_index) + 1
           end do
         end do
       end do
     end do
    end do

    ! send to root rank
    do s = 1, nspec
      send_spec_int(:) = spectra(s,:)
      call MPI_REDUCE(send_spec_int, recv_spec_int, spec_num, MPI_INTEGER,&
                    & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      glob_spectra(s,:) = recv_spec_int(:)

      #ifdef RADIATION
        ! compute radiation spectra
        if (allocated(rad_spectra) .and. allocated(glob_rad_spectra)) then
          send_spec_real(:) = rad_spectra(s,:)
          call MPI_REDUCE(send_spec_real, recv_spec_real, spec_num, MPI_REAL,&
                        & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
          glob_rad_spectra(s,:) = recv_spec_real(:)
          rad_spectra(s,:) = 0.0
        end if
      #endif
    end do

    if (allocated(spectra)) deallocate(spectra)
    if (allocated(send_spec_int)) deallocate(send_spec_int)
    if (allocated(recv_spec_int)) deallocate(recv_spec_int)
    if (allocated(send_spec_real)) deallocate(send_spec_real)
    if (allocated(recv_spec_real)) deallocate(recv_spec_real)
  end subroutine initializeOutput

  #ifdef HDF5
  subroutine writeXDMF_hdf5(step, time, ni, nj, nk)
    implicit none
    integer, intent(in)               :: step, time, ni, nj, nk
    character(len=STR_MAX)            :: stepchar, filename
    integer                           :: var

    write(stepchar, "(i5.5)") step
    filename = trim(output_dir_name) // '/flds.tot.' // trim(stepchar) // '.xdmf'

    open (UNIT_xdmf, file=filename, status="replace", access="stream", form="formatted")
    write (UNIT_xdmf, "(A)")&
                           & '<?xml version="1.0" ?>'
    write (UNIT_xdmf, "(A)")&
                           & '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
    write (UNIT_xdmf, "(A)")&
                           & '<Xdmf Version="2.0">'
    write (UNIT_xdmf, "(A)")&
                           & '  <Domain>'
    write (UNIT_xdmf, "(A)")&
                           & '    <Grid Name="domain" GridType="Uniform">'
    write (UNIT_xdmf, "(A,A,A,I10,I10,I10,A)")&
                           & '      <Topology TopologyType=',&
                              & '"3DCoRectMesh"', ' Dimensions="', &
                              & nk, nj, ni, '"/>'
    write (UNIT_xdmf, "(A)")&
                           & '      <Geometry GeometryType="ORIGIN_DXDYDZ">'
    write (UNIT_xdmf, "(A,A)")&
                           & '        <DataItem Format="XML" Dimensions="3"',&
                              & ' NumberType="Float" Precision="4">'
    write (UNIT_xdmf, "(A)")&
                           & '          0.0 0.0 0.0'
    write (UNIT_xdmf, "(A)")&
                           & '        </DataItem>'
    write (UNIT_xdmf, "(A,A)")&
                           & '        <DataItem Format="XML" Dimensions="3"',&
                              & ' NumberType="Float" Precision="4">'
    write (UNIT_xdmf, "(A)")&
                           & '          1.0 1.0 1.0'
    write (UNIT_xdmf, "(A)")&
                           & '        </DataItem>'
    write (UNIT_xdmf, "(A)")&
                           & '      </Geometry>'

    do var = 1, n_fld_vars
      write (UNIT_xdmf, "(A)")&
                           & '      <Attribute Name="' // trim(fld_vars(var)) // '" Center="Node">'
      write (UNIT_xdmf, "(A,I10,I10,I10,A)")&
                           & '        <DataItem Format="HDF" Dimensions="',&
                           & nk, nj, ni,&
                           & '" NumberType="Float" Precision="4">'
      write (UNIT_xdmf, "(A,A)")&
                           & '          flds.tot.' // trim(stepchar) // ':/',&
                           & trim(fld_vars(var))
      write (UNIT_xdmf, "(A)")&
                           & '        </DataItem>'
      write (UNIT_xdmf, "(A)")&
                           & '      </Attribute>'
    end do

    write (UNIT_xdmf, "(A)")&
                           & '    </Grid>'
    write (UNIT_xdmf, "(A)")&
                           & '  </Domain>'
    write (UNIT_xdmf, "(A)")&
                           & '</Xdmf>'
    close (UNIT_xdmf)
  end subroutine writeXDMF_hdf5

  subroutine writeFields_hdf5(step, time)
    implicit none
    integer, intent(in)               :: step, time
    character(len=STR_MAX)            :: stepchar, filename
    integer(HID_T)                    :: file_id, dset_id(40), filespace(40), memspace, plist_id
    integer                           :: error, f, s
    integer(kind=2)                   :: i, j, k
    logical                           :: writing_intQ, writing_lgarrQ
    integer                           :: dataset_rank = 3
    integer(HSSIZE_T), dimension(3)   :: offsets
    integer(HSIZE_T), dimension(3)    :: global_dims, blocks
    real                              :: ex0, ey0, ez0, bx0, by0, bz0
    real                              :: jx0, jy0, jz0

    ! downsampling variables
    integer :: this_x0, this_y0, this_z0, this_sx, this_sy, this_sz
    integer :: i_start, i_end, j_start, j_end, k_start, k_end
    integer :: offset_i, offset_j, offset_k, i1, j1, k1
    integer :: n_i, n_j, n_k, glob_n_i, glob_n_j, glob_n_k

    ! for convenience
    this_x0 = this_meshblock%ptr%x0
    this_y0 = this_meshblock%ptr%y0
    this_z0 = this_meshblock%ptr%z0
    this_sx = this_meshblock%ptr%sx
    this_sy = this_meshblock%ptr%sy
    this_sz = this_meshblock%ptr%sz

    write(stepchar, "(i5.5)") step
    filename = trim(output_dir_name) // '/flds.tot.' // trim(stepchar)

    ! assuming `global_mesh%{x0,y0,z0} .eq. 0`
    if (output_istep .eq. 1) then
      offset_i = this_x0;   offset_j = this_y0;   offset_k = this_z0
      n_i = this_sx - 1;    n_j = this_sy - 1;    n_k = this_sz - 1
      glob_n_i = global_mesh%sx
      glob_n_j = global_mesh%sy
      glob_n_k = global_mesh%sz

      i_start = 0; j_start = 0; k_start = 0
    else
      offset_i = CEILING(REAL(this_x0) / REAL(output_istep))
      offset_j = CEILING(REAL(this_y0) / REAL(output_istep))

      i_start = CEILING(REAL(this_x0) / REAL(output_istep)) * output_istep - this_x0
      i_end = (CEILING(REAL(this_x0 + this_sx) / REAL(output_istep)) - 1) * output_istep - this_x0
      j_start = CEILING(REAL(this_y0) / REAL(output_istep)) * output_istep - this_y0
      j_end = (CEILING(REAL(this_y0 + this_sy) / REAL(output_istep)) - 1) * output_istep - this_y0

      n_i = (i_end - i_start) / output_istep
      n_j = (j_end - j_start) / output_istep

      glob_n_i = CEILING(REAL(global_mesh%sx) / REAL(output_istep))
      glob_n_i = MAX(1, glob_n_i)

      glob_n_j = CEILING(REAL(global_mesh%sy) / REAL(output_istep))
      glob_n_j = MAX(1, glob_n_j)

      #ifndef threeD
        k_start = 0; k_end = 0
        offset_k = 0; n_k = 0
        glob_n_k = 1
      #else
        offset_k = CEILING(REAL(this_z0) / REAL(output_istep))
        k_start = CEILING(REAL(this_z0) / REAL(output_istep)) * output_istep - this_z0
        k_end = (CEILING(REAL(this_z0 + this_sz) / REAL(output_istep)) - 1) * output_istep - this_z0
        n_k = (k_end - k_start) / output_istep
        glob_n_k = CEILING(REAL(global_mesh%sz) / REAL(output_istep))
        glob_n_k = MAX(1, glob_n_k)
      #endif
    end if

    if (mpi_rank .eq. 0) then
      call writeXDMF_hdf5(step, time, glob_n_i, glob_n_j, glob_n_k)
    end if

    offsets(1) = offset_i
    offsets(2) = offset_j
    offsets(3) = offset_k
    blocks(1) = n_i + 1
    blocks(2) = n_j + 1
    blocks(3) = n_k + 1
    global_dims(1) = glob_n_i
    global_dims(2) = glob_n_j
    global_dims(3) = glob_n_k

    ! mpi_f08 thing
    ! call MPI_INFO_CREATE(FILE_INFO_TEMPLATE, error)
    ! call MPI_INFO_SET(FILE_INFO_TEMPLATE, "access_style", "write_once", error)
    ! call MPI_INFO_SET(FILE_INFO_TEMPLATE, "collective_buffering", "true", error)
    ! call MPI_INFO_SET(FILE_INFO_TEMPLATE, "cb_block_size", "4194304", error)
    ! call MPI_INFO_SET(FILE_INFO_TEMPLATE, "cb_buffer_size", "16777216", error)
    ! call MPI_INFO_SET(FILE_INFO_TEMPLATE, "cb_nodes", "1", error)

    call h5open_f(error)
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, h5comm, h5info, error)
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id, error)

    do f = 1, n_fld_vars
      call h5screate_simple_f(dataset_rank, global_dims, filespace(f), error)
    end do

    do f = 1, n_fld_vars
      if (fld_vars(f)(1:4) .eq. 'dens') then
        writing_lgarrQ = .true.
        s = STRtoINT(fld_vars(f)(5:5))
        call computeDensity(s, reset=.true.) ! filled `lg_arr` with density of species `s`
        call exchangeArray()
      else if (fld_vars(f)(1:4) .eq. 'enrg') then
        writing_lgarrQ = .true.
        s = STRtoINT(fld_vars(f)(5:5))
        call computeEnergy(s, reset=.true.) ! filled `lg_arr` with energies of species `s`
        call exchangeArray()
      else
        writing_lgarrQ = .false.
      end if

      call h5dcreate_f(file_id, fld_vars(f), H5T_NATIVE_REAL, filespace(f), &
                     & dset_id(f), error)
      call h5sclose_f(filespace(f), error)
      call h5dget_space_f(dset_id(f), filespace(f), error)

      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

      call h5screate_simple_f(dataset_rank, blocks, memspace, error)
      call h5sselect_hyperslab_f(filespace(f), H5S_SELECT_SET_F, offsets, blocks, error)

      ! Create dataset by interpolating fields
      do i1 = 0, n_i
        do j1 = 0, n_j
          do k1 = 0, n_k
            i = i_start + i1 * output_istep
            j = j_start + j1 * output_istep
            k = k_start + k1 * output_istep
            select case (trim(fld_vars(f)))
            case('ex')
              call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
              sm_arr(i1, j1, k1) = ex0 * B_norm
            case('ey')
              call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
              sm_arr(i1, j1, k1) = ey0 * B_norm
            case('ez')
              call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
              sm_arr(i1, j1, k1) = ez0 * B_norm
            case('bx')
              call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
              sm_arr(i1, j1, k1) = bx0 * B_norm
            case('by')
              call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
              sm_arr(i1, j1, k1) = by0 * B_norm
            case('bz')
              call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
              sm_arr(i1, j1, k1) = bz0 * B_norm
            case('jx')
              call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
              sm_arr(i1, j1, k1) = -jx0 * B_norm
            case('jy')
              call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
              sm_arr(i1, j1, k1) = -jy0 * B_norm
            case('jz')
              call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
              sm_arr(i1, j1, k1) = -jz0 * B_norm
            case('xx')
              sm_arr(i1, j1, k1) = REAL(this_meshblock%ptr%x0 + i, 4)
            case('yy')
              sm_arr(i1, j1, k1) = REAL(this_meshblock%ptr%y0 + j, 4)
            case('zz')
              sm_arr(i1, j1, k1) = REAL(this_meshblock%ptr%z0 + k, 4)
            case default
              if (((fld_vars(f)(1:4) .ne. 'dens') .and. (fld_vars(f)(1:4) .ne. 'enrg')) .or.&
                 & (.not. writing_lgarrQ)) then
                call throwError("ERROR: unrecognized `fld_vars(f)`")
              else
                sm_arr(i1, j1, k1) = lg_arr(i, j, k)
              end if
            end select
          end do
        end do
      end do

      ! Write the dataset collectively
      call h5dwrite_f(dset_id(f), H5T_NATIVE_REAL, sm_arr(0 : n_i, 0 : n_j, 0 : n_k), global_dims, error, &
                    & file_space_id = filespace(f), mem_space_id = memspace, xfer_prp = plist_id)

      call h5dclose_f(dset_id(f), error)
      call h5sclose_f(filespace(f), error)
    end do

    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error)
    call h5fclose_f(file_id, error)
    call h5close_f(error)
  end subroutine writeFields_hdf5

  subroutine writeParticles_hdf5(step, time)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in)               :: step, time
    character(len=STR_MAX)            :: stepchar, filename
    character(len=7)                  :: dsetname
    integer(HID_T)                    :: file_id, dset_id(100), filespace(100), memspace, plist_id
    integer                           :: error, ierr
    integer                           :: rnk, s, dummy_s, p, j, ln_, ti, tj, tk, temp, temp_int
    integer                           :: dataset_rank = 1
    integer(HID_T)                    :: h5type
    integer(HSSIZE_T), dimension(1)   :: offsets
    integer(HSIZE_T), dimension(1)    :: global_dims, blocks
    integer                           :: npart_stride(nspec), npart_stride_global(nspec, mpi_size)
    integer, allocatable, dimension(:):: temp_int_arr, stride_indices_arr, stride_ti_arr, stride_tj_arr, stride_tk_arr
    real, allocatable, dimension(:)   :: temp_real_arr
    real                              :: temp_real1, temp_real2, temp_real3
    logical                           :: writing_intQ

    ! preparation
    ! number of strided particles per each species
    do s = 1, nspec
      npart_stride(s) = 0
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            !$omp simd
            !dir$ vector aligned
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              if (modulo(species(s)%prtl_tile(ti, tj, tk)%ind(p), output_stride) .eq. 0) then
                npart_stride(s) = npart_stride(s) + 1
              end if
            end do ! p
          end do ! tk
        end do ! tj
      end do ! ti
    end do ! s

    call MPI_ALLGATHER(npart_stride, nspec, MPI_INTEGER,&
                    & npart_stride_global, nspec, MPI_INTEGER,&
                    & MPI_COMM_WORLD, ierr)

    write(stepchar, "(i5.5)") step
    filename = trim(output_dir_name) // '/prtl.tot.' // trim(stepchar)

    call h5open_f(error)
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, h5comm, h5info, error)
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id, error)

    do s = 1, nspec
      offsets(1) = 0
      global_dims(1) = 0
      do rnk = 0, mpi_size - 1
        if (rnk .lt. mpi_rank) then
          offsets(1) = offsets(1) + npart_stride_global(s, rnk + 1)
        end if
        global_dims(1) = global_dims(1) + npart_stride_global(s, rnk + 1)
      end do
      blocks(1) = npart_stride(s)

      ! saving the particle (and tile) indices to output for a given species
      allocate(stride_indices_arr(npart_stride(s)))
      allocate(stride_ti_arr(npart_stride(s)))
      allocate(stride_tj_arr(npart_stride(s)))
      allocate(stride_tk_arr(npart_stride(s)))
      j = 1
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            !$omp simd
            !dir$ vector aligned
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              if (modulo(species(s)%prtl_tile(ti, tj, tk)%ind(p), output_stride) .eq. 0) then
                stride_indices_arr(j) = p
                stride_ti_arr(j) = ti
                stride_tj_arr(j) = tj
                stride_tk_arr(j) = tk
                j = j + 1
              end if
            end do ! particles
          end do ! tk
        end do ! tj
      end do ! ti

      do p = 1, n_prtl_vars
        call h5screate_simple_f(dataset_rank, global_dims, filespace(p), error)
      end do

      do p = 1, n_prtl_vars
        ! dataset name `var_name` + '_' + `species #`
        ln_ = len(trim(prtl_vars(p)))
        dsetname(1 : ln_ + 3) = trim(prtl_vars(p)) // '_' // trim(STR(s))
        dsetname(ln_ + 3 : 7) = ' '
        ! creating dataset for a given type
        if (trim(prtl_var_types(p)) .eq. 'int') then
          writing_intQ = .true.
          h5type = H5T_NATIVE_INTEGER
          ! Create dataset
          allocate(temp_int_arr(npart_stride(s)))
          do j = 1, npart_stride(s)
            temp = stride_indices_arr(j)
            ti = stride_ti_arr(j)
            tj = stride_tj_arr(j)
            tk = stride_tk_arr(j)
            select case (trim(prtl_vars(p))) ! select integer variable
              case('ind')
                temp_int = species(s)%prtl_tile(ti, tj, tk)%ind(temp)
                temp_int_arr(j) = temp_int
              case('proc')
                temp_int = species(s)%prtl_tile(ti, tj, tk)%proc(temp)
                temp_int_arr(j) = temp_int
              case default
                call throwError('ERROR: unrecognized `prtl_vars`: `'//trim(prtl_vars(p))//'`')
            end select
          end do
        else if (trim(prtl_var_types(p)) .eq. 'real') then
          writing_intQ = .false.
          h5type = H5T_NATIVE_REAL
          if (prtl_vars(p)(1:4) .eq. 'dens') then
            dummy_s = STRtoINT(prtl_vars(p)(5:5))
            call computeDensity(dummy_s, reset=.true.) ! filled `lg_arr` with density of species `s`
            call exchangeArray()
          end if
          ! Create dataset
          allocate(temp_real_arr(npart_stride(s)))
          do j = 1, npart_stride(s)
            temp = stride_indices_arr(j)
            ti = stride_ti_arr(j)
            tj = stride_tj_arr(j)
            tk = stride_tk_arr(j)
            select case (trim(prtl_vars(p)))
              case('x')
                temp_int = species(s)%prtl_tile(ti, tj, tk)%xi(temp)
                temp_real1 = species(s)%prtl_tile(ti, tj, tk)%dx(temp)
                temp_real_arr(j) = REAL(this_meshblock%ptr%x0 + temp_int) + temp_real1
              case('y')
                temp_int = species(s)%prtl_tile(ti, tj, tk)%yi(temp)
                temp_real1 = species(s)%prtl_tile(ti, tj, tk)%dy(temp)
                temp_real_arr(j) = REAL(this_meshblock%ptr%y0 + temp_int) + temp_real1
              case('z')
                temp_int = species(s)%prtl_tile(ti, tj, tk)%zi(temp)
                temp_real1 = species(s)%prtl_tile(ti, tj, tk)%dz(temp)
                temp_real_arr(j) = REAL(this_meshblock%ptr%z0 + temp_int) + temp_real1
              case('u')
                temp_real1 = species(s)%prtl_tile(ti, tj, tk)%u(temp)
                temp_real_arr(j) = REAL(temp_real1, 4)
              case('v')
                temp_real1 = species(s)%prtl_tile(ti, tj, tk)%v(temp)
                temp_real_arr(j) = REAL(temp_real1, 4)
              case('w')
                temp_real1 = species(s)%prtl_tile(ti, tj, tk)%w(temp)
                temp_real_arr(j) = REAL(temp_real1, 4)
              case('ex')
                call interpFromEdges(species(s)%prtl_tile(ti, tj, tk)%dx(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%dy(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%dz(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%xi(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%yi(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%zi(temp),&
                                   & ex, ey, ez, temp_real1, temp_real2, temp_real3)
                temp_real_arr(j) = REAL(temp_real1 * B_norm, 4)
              case('ey')
                call interpFromEdges(species(s)%prtl_tile(ti, tj, tk)%dx(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%dy(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%dz(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%xi(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%yi(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%zi(temp),&
                                   & ex, ey, ez, temp_real1, temp_real2, temp_real3)
                temp_real_arr(j) = REAL(temp_real2 * B_norm, 4)
              case('ez')
                call interpFromEdges(species(s)%prtl_tile(ti, tj, tk)%dx(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%dy(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%dz(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%xi(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%yi(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%zi(temp),&
                                   & ex, ey, ez, temp_real1, temp_real2, temp_real3)
                temp_real_arr(j) = REAL(temp_real3 * B_norm, 4)
              case('bx')
                call interpFromFaces(species(s)%prtl_tile(ti, tj, tk)%dx(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%dy(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%dz(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%xi(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%yi(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%zi(temp),&
                                   & bx, by, bz, temp_real1, temp_real2, temp_real3)
                temp_real_arr(j) = REAL(temp_real1 * B_norm, 4)
              case('by')
                call interpFromFaces(species(s)%prtl_tile(ti, tj, tk)%dx(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%dy(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%dz(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%xi(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%yi(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%zi(temp),&
                                   & bx, by, bz, temp_real1, temp_real2, temp_real3)
                temp_real_arr(j) = REAL(temp_real2 * B_norm, 4)
              case('bz')
                call interpFromFaces(species(s)%prtl_tile(ti, tj, tk)%dx(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%dy(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%dz(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%xi(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%yi(temp),&
                                   & species(s)%prtl_tile(ti, tj, tk)%zi(temp),&
                                   & bx, by, bz, temp_real1, temp_real2, temp_real3)
                temp_real_arr(j) = REAL(temp_real3 * B_norm, 4)
              case default
                if (prtl_vars(p)(1:4) .ne. 'dens') then
                  call throwError('ERROR: unrecognized `prtl_vars`: `'//trim(prtl_vars(p))//'`')
                else
                  temp_real_arr(j) = lg_arr(species(s)%prtl_tile(ti, tj, tk)%xi(temp),&
                                          & species(s)%prtl_tile(ti, tj, tk)%yi(temp),&
                                          & species(s)%prtl_tile(ti, tj, tk)%zi(temp))
                end if
            end select ! select variable
          end do ! strided prtls
        else ! if unrecognized vartype
          call throwError('ERROR: unrecognized `prtl_var_types`: `'//trim(prtl_var_types(p))//'`')
        end if

        call h5dcreate_f(file_id, dsetname, h5type, filespace(p),&
                       & dset_id(p), error)
        call h5sclose_f(filespace(p), error)
        call h5dget_space_f(dset_id(p), filespace(p), error)

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        call h5screate_simple_f(dataset_rank, blocks, memspace, error)
        call h5sselect_hyperslab_f(filespace(p), H5S_SELECT_SET_F, offsets, blocks, error)

        ! Write the dataset collectively
        if (writing_intQ) then
          call h5dwrite_f(dset_id(p), h5type, temp_int_arr, global_dims, error,&
                        & file_space_id = filespace(p), mem_space_id = memspace,&
                        & xfer_prp = plist_id)
          deallocate(temp_int_arr)
        else
          call h5dwrite_f(dset_id(p), h5type, temp_real_arr, global_dims, error,&
                        & file_space_id = filespace(p), mem_space_id = memspace,&
                        & xfer_prp = plist_id)
          deallocate(temp_real_arr)
        end if

        call h5sclose_f(filespace(p), error)
        call h5dclose_f(dset_id(p), error)
      end do
      deallocate(stride_indices_arr)
      deallocate(stride_ti_arr)
      deallocate(stride_tj_arr)
      deallocate(stride_tk_arr)
    end do

    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error)
    call h5fclose_f(file_id, error)
    call h5close_f(error)
  end subroutine writeParticles_hdf5

  subroutine writeSpectra_hdf5(step, time)
    implicit none
    integer, intent(in)               :: step, time
    character(len=STR_MAX)            :: stepchar, filename
    integer                           :: error, s, i, datarank
    integer(HID_T)                    :: file_id, dset_id, dspace_id
    integer(HSIZE_T), dimension(1)    :: data_dims
    character(len=3)                  :: dsetname
    real, allocatable, dimension(:)   :: bin_data

    datarank = 1
    data_dims(1) = spec_num

    ! only root rank writes the spectra file
    if (mpi_rank .eq. 0) then
      ! saving the energy bins
      allocate(bin_data(spec_num))
      do i = 1, spec_num
        bin_data(i) = spec_min + (REAL(i - 1, 4) / REAL(spec_num, 4)) * (spec_max - spec_min)
      end do

      write(stepchar, "(i5.5)") step
      filename = trim(output_dir_name) // '/spec.tot.' // trim(stepchar)

      ! Initialize FORTRAN interface
      call h5open_f(error)
      ! Create a new file using default properties
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      do s = 1, nspec
        ! writing bins:
        dsetname = 'e' // trim(STR(s))
        call h5screate_simple_f(datarank, data_dims, dspace_id, error)
        call h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, dspace_id, &
                       & dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, bin_data, data_dims, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        ! writing spectra:
        dsetname = 'n' // trim(STR(s))
        call h5screate_simple_f(datarank, data_dims, dspace_id, error)
        call h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, dspace_id, &
                       & dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, glob_spectra(s,:), data_dims, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        #ifdef RADIATION
          if (allocated(glob_rad_spectra)) then
            ! writing bins:
            dsetname = 'er' // trim(STR(s))
            call h5screate_simple_f(datarank, data_dims, dspace_id, error)
            call h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, dspace_id, &
                           & dset_id, error)
            call h5dwrite_f(dset_id, H5T_NATIVE_REAL, bin_data, data_dims, error)
            call h5dclose_f(dset_id, error)
            call h5sclose_f(dspace_id, error)

            ! writing spectra:
            dsetname = 'nr' // trim(STR(s))
            call h5screate_simple_f(datarank, data_dims, dspace_id, error)
            call h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, dspace_id, &
                           & dset_id, error)
            call h5dwrite_f(dset_id, H5T_NATIVE_REAL, glob_rad_spectra(s,:), data_dims, error)
            call h5dclose_f(dset_id, error)
            call h5sclose_f(dspace_id, error)
            glob_rad_spectra(s,:) = 0.0
          end if
        #endif
      end do

      ! Close the file
      call h5fclose_f(file_id, error)
      ! Close FORTRAN interface
      call h5close_f(error)

      if (allocated(bin_data)) deallocate(bin_data)
    end if

    if (allocated(glob_spectra)) deallocate(glob_spectra)
  end subroutine writeSpectra_hdf5

  subroutine writeDomain_hdf5(step, time)
    implicit none
    integer, intent(in)               :: step, time
    character(len=STR_MAX)            :: stepchar, filename
    integer                           :: error, s, i, datarank, d, rnk
    integer(HID_T)                    :: file_id, dset_id, dspace_id
    integer(HSIZE_T), dimension(1)    :: data_dims
    integer, allocatable, dimension(:):: domain_data

    datarank = 1
    data_dims(1) = mpi_size

    ! only root rank writes the domain file
    if (mpi_rank .eq. 0) then
      write(stepchar, "(i5.5)") step
      filename = trim(output_dir_name) // '/domain.' // trim(stepchar)

      ! Initialize FORTRAN interface
      call h5open_f(error)
      ! Create a new file using default properties
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      allocate(domain_data(mpi_size))

      do d = 1, n_dom_vars
        select case (trim(dom_vars(d)))
          case('x0')
            do rnk = 0, mpi_size - 1
              domain_data(rnk + 1) = meshblocks(rnk + 1)%x0
            end do
          case('y0')
            do rnk = 0, mpi_size - 1
              domain_data(rnk + 1) = meshblocks(rnk + 1)%y0
            end do
          case('z0')
            do rnk = 0, mpi_size - 1
              domain_data(rnk + 1) = meshblocks(rnk + 1)%z0
            end do
          case('sx')
            do rnk = 0, mpi_size - 1
              domain_data(rnk + 1) = meshblocks(rnk + 1)%sx
            end do
          case('sy')
            do rnk = 0, mpi_size - 1
              domain_data(rnk + 1) = meshblocks(rnk + 1)%sy
            end do
          case('sz')
            do rnk = 0, mpi_size - 1
              domain_data(rnk + 1) = meshblocks(rnk + 1)%sz
            end do
          case default
            call throwError('ERROR: unrecognized `dom_vars`: `'//trim(dom_vars(d))//'`')
        end select

        call h5screate_simple_f(datarank, data_dims, dspace_id, error)
        call h5dcreate_f(file_id, dom_vars(d), H5T_NATIVE_INTEGER, dspace_id, &
                       & dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, domain_data, data_dims, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)
      end do

      ! Close the file
      call h5fclose_f(file_id, error)
      ! Close FORTRAN interface
      call h5close_f(error)

      if (allocated(domain_data)) deallocate(domain_data)
    end if
  end subroutine writeDomain_hdf5
  #endif

end module m_writeoutput
