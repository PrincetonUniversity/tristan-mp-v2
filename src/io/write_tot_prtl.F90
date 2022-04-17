#include "../defs.F90"

module m_writetotprtl
  #ifdef HDF5
    use hdf5
    use m_globalnamespace, only: h5comm, h5info
  #endif
  use m_globalnamespace, only: output_dir_name, mpi_rank
  use m_outputnamespace, only: tot_output_stride, n_prtl_vars,&
                             & prtl_vars, prtl_var_types
  use m_aux
  use m_errors, only: throwError
  use m_domain
  use m_fields
  use m_particles
  use m_helpers, only: computeDensity, interpFromFaces, interpFromEdges
  use m_exchangearray, only: exchangeArray
  
  #ifdef USROUTPUT
    use m_userfile, only: userExcludeParticles
  #endif

  implicit none

  !--- PRIVATE functions -----------------------------------------!
  #ifdef HDF5
    private :: writeParticles_hdf5
  #endif
  !...............................................................!
contains
  subroutine writeParticles(step, time)
    implicit none
    integer, intent(in)               :: step, time
    #ifdef HDF5
      call writeParticles_hdf5(step, time)
    #endif
    call printDiag("writeParticles()", 3)
  end subroutine writeParticles

  logical function particleIsEligible(s, ti, tj, tk, p)
    implicit none
    integer, intent(in)     :: s, ti, tj, tk, p
    logical                 :: dummy
    dummy = (modulo(species(s)%prtl_tile(ti, tj, tk)%ind(p), tot_output_stride) .eq. 0)
    #ifdef USROUTPUT
      dummy = (dummy .and. userExcludeParticles(s, ti, tj, tk, p))
    #endif
    particleIsEligible = dummy
  end function particleIsEligible

  #ifdef HDF5
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
        if (.not. species(s)%output_sp) cycle
        do ti = 1, species(s)%tile_nx
          do tj = 1, species(s)%tile_ny
            do tk = 1, species(s)%tile_nz
              do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
                if (particleIsEligible(s, ti, tj, tk, p)) then
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

        if (npart_stride(s) .gt. 0) then
          j = 1
          do ti = 1, species(s)%tile_nx
            do tj = 1, species(s)%tile_ny
              do tk = 1, species(s)%tile_nz
                do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
                  if (particleIsEligible(s, ti, tj, tk, p)) then
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
        end if

        do p = 1, n_prtl_vars
          call h5screate_simple_f(dataset_rank, global_dims, filespace(p), error)
        end do

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        do p = 1, n_prtl_vars
          ! dataset name `var_name` + '_' + `species #`
          ln_ = len(trim(prtl_vars(p)))
          dsetname(1 : ln_ + 2) = trim(prtl_vars(p)) // '_' // trim(STR(s))
          dsetname(ln_ + 3 : 7) = ' '
          ! creating dataset for a given type
          if (trim(prtl_var_types(p)) .eq. 'int') then
            writing_intQ = .true.
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
                case('wei')
                  temp_real1 = species(s)%prtl_tile(ti, tj, tk)%weight(temp)
                  temp_real_arr(j) = temp_real1
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
                  if (prtl_vars(p)(1:4) .eq. 'dens') then
                    temp_real_arr(j) = lg_arr(species(s)%prtl_tile(ti, tj, tk)%xi(temp),&
                                            & species(s)%prtl_tile(ti, tj, tk)%yi(temp),&
                                            & species(s)%prtl_tile(ti, tj, tk)%zi(temp))
                  #ifdef PRTLPAYLOADS
                    else if (prtl_vars(p)(1:3) .eq. 'pld') then
                      dummy_s = STRtoINT(prtl_vars(p)(4:4))
                      if (dummy_s .eq. 1) then
                        temp_real_arr(j) = species(s)%prtl_tile(ti, tj, tk)%payload1(temp)
                      else if (dummy_s .eq. 2) then
                        temp_real_arr(j) = species(s)%prtl_tile(ti, tj, tk)%payload2(temp)
                      else if (dummy_s .eq. 3) then
                        temp_real_arr(j) = species(s)%prtl_tile(ti, tj, tk)%payload3(temp)
                      else
                        call throwError('ERROR: only 3 payloads are allowed.')
                      end if
                  #endif
                  else
                    call throwError('ERROR: unrecognized `prtl_vars`: `'//trim(prtl_vars(p))//'`')
                  end if
              end select ! select variable
            end do ! strided prtls
          else ! if unrecognized vartype
            call throwError('ERROR: unrecognized `prtl_var_types`: `'//trim(prtl_var_types(p))//'`')
          end if

          if (writing_intQ) then
            call h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, filespace(p), dset_id(p), error)
          else
            call h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, filespace(p), dset_id(p), error)
          end if
          call h5sclose_f(filespace(p), error)
          call h5dget_space_f(dset_id(p), filespace(p), error)

          call h5screate_simple_f(dataset_rank, blocks, memspace, error)
          call h5dget_space_f(dset_id(p), filespace(p), error)
          call h5sselect_hyperslab_f(filespace(p), H5S_SELECT_SET_F, offsets, blocks, error)

          ! Write the dataset collectively
          if (writing_intQ) then
            call h5dwrite_f(dset_id(p), H5T_NATIVE_INTEGER, temp_int_arr, global_dims, error,&
                          & file_space_id = filespace(p), mem_space_id = memspace,&
                          & xfer_prp = h5p_default_f)
            deallocate(temp_int_arr)
          else
            call h5dwrite_f(dset_id(p), H5T_NATIVE_REAL, temp_real_arr, global_dims, error,&
                          & file_space_id = filespace(p), mem_space_id = memspace,&
                          & xfer_prp = h5p_default_f)
            deallocate(temp_real_arr)
          end if
          call h5sclose_f(memspace, error)
          call h5dclose_f(dset_id(p), error)
          call h5sclose_f(filespace(p), error)
        end do
        deallocate(stride_indices_arr)
        deallocate(stride_ti_arr)
        deallocate(stride_tj_arr)
        deallocate(stride_tk_arr)
      end do

      call h5pclose_f(plist_id, error)
      call h5fclose_f(file_id, error)
      call h5close_f(error)
    end subroutine writeParticles_hdf5

  #endif

end module m_writetotprtl
