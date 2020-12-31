#include "../defs.F90"

module m_outputlogistics
  use m_globalnamespace
  use m_outputnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_readinput, only: getInput
  use m_helpers, only: computeDensity, computeMomentum, computeNpart
  use m_helpers, only: interpFromFaces, interpFromEdges
  use m_exchangearray, only: exchangeArray

  implicit none
contains

  subroutine initializeOutput()
    implicit none
    call getInput('output', 'enable', tot_output_enable, .true.)
    ! individual `.tot.`, `diag` and `spec` outputs
    call getInput('output', 'params_enable', params_enable, .true.)
    call getInput('output', 'prtl_enable', prtl_tot_enable, .true.)
    call getInput('output', 'flds_enable', flds_tot_enable, .true.)
    call getInput('output', 'spec_enable', spectra_enable, .true.)

    call getInput('output', 'diag_enable', diag_enable, .false.)

    call getInput('output', 'start', tot_output_start, 0)
    call getInput('output', 'interval', tot_output_interval, 10)
    call getInput('output', 'stride', tot_output_stride, 10)
    call getInput('output', 'istep', output_flds_istep, 1)
    call getInput('output', 'smooth_window', output_dens_smooth, 2)

    call getInput('output', 'spec_log_bins', spec_log_bins, .true.)
    call getInput('output', 'spec_min', spec_min, 1e-2)
    call getInput('output', 'spec_max', spec_max, 1e2)
    call getInput('output', 'spec_num', spec_num, 100)

    call getInput('output', 'spec_dynamic_bins', spec_dynamic_bins, .false.)

    #ifdef oneD
      call getInput('output', 'spec_nx', spec_nx, 1)
      spec_ny = 1; spec_nz = 1
    #elif twoD
      call getInput('output', 'spec_nx', spec_nx, 1)
      call getInput('output', 'spec_ny', spec_ny, 1)
      spec_nz = 1
    #elif threeD
      call getInput('output', 'spec_nx', spec_nx, 1)
      call getInput('output', 'spec_ny', spec_ny, 1)
      call getInput('output', 'spec_nz', spec_nz, 1)
    #endif

    if (spec_log_bins) then
      spec_min = log(spec_min)
      spec_max = log(spec_max)
    endif
    spec_bin_size = (spec_max - spec_min) / spec_num

    call getInput('output', 'flds_at_prtl', flds_at_prtl_enable, .false.)
    call getInput('output', 'write_xdmf', xdmf_enable, .true.)
    call getInput('output', 'write_nablas', derivatives_enable, .false.)
    call getInput('output', 'write_momenta', momenta_enable, .false.)
    call getInput('output', 'write_npart', npart_enable, .false.)

    #if defined(HDF5) && defined(MPI08)
      h5comm = MPI_COMM_WORLD%MPI_VAL
      h5info = MPI_INFO_NULL%MPI_VAL
    #elif defined(HDF5) && defined(MPI)
      h5comm = MPI_COMM_WORLD
      h5info = MPI_INFO_NULL
    #endif
  end subroutine initializeOutput

  subroutine initializeSlice()
    implicit none
    integer                 :: i
    character(len=STR_MAX)  :: var_name
    call getInput('slice_output', 'enable', slice_output_enable, .false.)
    call getInput('slice_output', 'start', slice_output_start, 0)
    call getInput('slice_output', 'interval', slice_output_interval, 10)

    #ifndef threeD
      slice_output_enable = .false.
    #endif

    slice_axes(:) = -1
    slice_pos(:) = -1

    do i = 1, 100
      write (var_name, "(A7,I1)") "sliceX_", i
      call getInput('slice_output', var_name, slice_pos(nslices + 1), -1)
      if (slice_pos(nslices + 1) .ne. -1) then
        nslices = nslices + 1
        slice_axes(nslices) = 1
      else
        exit
      end if
    end do

    do i = 1, 100
      write (var_name, "(A7,I1)") "sliceY_", i
      call getInput('slice_output', var_name, slice_pos(nslices + 1), -1)
      if (slice_pos(nslices + 1) .ne. -1) then
        nslices = nslices + 1
        slice_axes(nslices) = 2
      else
        exit
      end if
    end do

    do i = 1, 100
      write (var_name, "(A7,I1)") "sliceZ_", i
      call getInput('slice_output', var_name, slice_pos(nslices + 1), -1)
      if (slice_pos(nslices + 1) .ne. -1) then
        nslices = nslices + 1
        slice_axes(nslices) = 3
      else
        exit
      end if
    end do

  end subroutine initializeSlice

  subroutine prepareOutput()
    ! DEP_PRT [particle-dependent]
    implicit none
    integer                   :: s
    integer                   :: ierr, ndown, pid
    ! initialize particle variables
    n_prtl_vars = 9
    prtl_vars(1:n_prtl_vars) = (/'x    ', 'y    ', 'z    ',&
                               & 'u    ', 'v    ', 'w    ',&
                               & 'wei  ', 'ind  ', 'proc '/)
    prtl_var_types(1:n_prtl_vars) = (/'real ', 'real ', 'real ',&
                                    & 'real ', 'real ', 'real ',&
                                    & 'real ', 'int  ', 'int  '/)
    if (flds_at_prtl_enable) then
      n_prtl_vars = n_prtl_vars + 6
      prtl_vars(10:n_prtl_vars) = (/'ex   ', 'ey   ', 'ez   ',&
                                  & 'bx   ', 'by   ', 'bz   '/)
      prtl_var_types(10:n_prtl_vars) = (/'real ', 'real ', 'real ',&
                                       & 'real ', 'real ', 'real '/)
      do s = 1, nspec
        prtl_vars(n_prtl_vars + s) = 'dens' // STR(s)
        prtl_var_types(n_prtl_vars + s) = 'real '
      end do
      n_prtl_vars = n_prtl_vars + nspec
    end if

    #ifdef PRTLPAYLOADS
      do pid = 1, 3
        prtl_vars(n_prtl_vars + pid) = 'pld' // STR(pid)
        prtl_var_types(n_prtl_vars + pid) = 'real '
      end do
      n_prtl_vars = n_prtl_vars + 3
    #endif

    call prepareSpectraForOutput()
    call defineFieldVarsToOutput()

    ! initialize domain output variables
    !   FIX1: maybe add # of particles per domain
    n_dom_vars = 6
    dom_vars(1 : 6) = (/'x0   ', 'y0   ', 'z0   ',&
                      & 'sx   ', 'sy   ', 'sz   '/)

  end subroutine prepareOutput

  subroutine prepareSpectraForOutput()
    implicit none
    real                      :: energy, u_, v_, w_, x_g, y_g, z_g, emax, glob_emax
    integer                   :: s, i, ti, tj, tk, p, spec_index
    integer                   :: spec_x_index, spec_y_index, spec_z_index
    integer                   :: ierr, root_rnk
    real, allocatable         :: spectra(:,:,:,:,:)
    real, allocatable         :: send_spec(:,:,:,:), recv_spec(:,:,:,:)

    root_rnk = 0

    if (spec_dynamic_bins) then
      ! find global max energy
      emax = 0.0
      if (spec_log_bins) emax = -10.0
      do s = 1, nspec
        do ti = 1, species(s)%tile_nx
          do tj = 1, species(s)%tile_ny
            do tk = 1, species(s)%tile_nz
              do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
                u_ = species(s)%prtl_tile(ti, tj, tk)%u(p)
                v_ = species(s)%prtl_tile(ti, tj, tk)%v(p)
                w_ = species(s)%prtl_tile(ti, tj, tk)%w(p)
                if ((species(s)%m_sp .eq. 0) .and. (species(s)%ch_sp .eq. 0)) then
                  energy = sqrt(u_**2 + v_**2 + w_**2)
                else
                  energy = sqrt(1.0 + u_**2 + v_**2 + w_**2) - 1.0
                end if
                if (spec_log_bins) energy = log(energy + 1e-8)
                if (energy .gt. emax) emax = energy
              end do
            end do
          end do
        end do
      end do
      call MPI_ALLREDUCE(emax, glob_emax, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr)
      if (glob_emax .gt. spec_max) then
        spec_max = glob_emax
        spec_num = INT(CEILING((spec_max - spec_min) / spec_bin_size))
      end if
    end if

    ! compute spectra
    if (mpi_rank .eq. root_rnk) then
      if (.not. allocated(glob_spectra)) then
        allocate(glob_spectra(nspec, spec_nx, spec_ny, spec_nz, spec_num))
      end if
    end if

    allocate(spectra(nspec, spec_nx, spec_ny, spec_nz, spec_num))
    allocate(send_spec(spec_nx, spec_ny, spec_nz, spec_num), recv_spec(spec_nx, spec_ny, spec_nz, spec_num))
    spectra(:,:,:,:,:) = 0

    do s = 1, nspec
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              x_g = REAL(this_meshblock%ptr%x0 + species(s)%prtl_tile(ti, tj, tk)%xi(p)) +&
                      & species(s)%prtl_tile(ti, tj, tk)%dx(p)
              y_g = REAL(this_meshblock%ptr%y0 + species(s)%prtl_tile(ti, tj, tk)%yi(p)) +&
                      & species(s)%prtl_tile(ti, tj, tk)%dy(p)
              z_g = REAL(this_meshblock%ptr%z0 + species(s)%prtl_tile(ti, tj, tk)%zi(p)) +&
                      & species(s)%prtl_tile(ti, tj, tk)%dz(p)

              ! find energy bin
              u_ = species(s)%prtl_tile(ti, tj, tk)%u(p)
              v_ = species(s)%prtl_tile(ti, tj, tk)%v(p)
              w_ = species(s)%prtl_tile(ti, tj, tk)%w(p)
              if ((species(s)%m_sp .eq. 0) .and. (species(s)%ch_sp .eq. 0)) then
                energy = sqrt(u_**2 + v_**2 + w_**2)
              else
                energy = sqrt(1.0 + u_**2 + v_**2 + w_**2) - 1.0
              end if
              if (spec_log_bins) energy = log(energy + 1e-8)
              if (energy .le. spec_min) then
                spec_index = 1
              else if (energy .ge. spec_max) then
                spec_index = spec_num
              else
                spec_index = INT(CEILING((energy - spec_min) * REAL(spec_num) / (spec_max - spec_min)))
                if (spec_index .lt. 1) spec_index = 1
                if (spec_index .gt. spec_num) spec_index = spec_num
              end if

              ! find spatial bin
              spec_x_index = CEILING(x_g / (REAL(global_mesh%sx) / REAL(spec_nx)))
              spec_y_index = CEILING(y_g / (REAL(global_mesh%sy) / REAL(spec_ny)))
              spec_z_index = CEILING(z_g / (REAL(global_mesh%sz) / REAL(spec_nz)))

              spectra(s, spec_x_index, spec_y_index, spec_z_index, spec_index) =&
                    & spectra(s, spec_x_index, spec_y_index, spec_z_index, spec_index) + species(s)%prtl_tile(ti, tj, tk)%weight(p)
            end do
          end do
        end do
      end do
    end do

    ! send to root rank
    do s = 1, nspec
      send_spec(:,:,:,:) = spectra(s,:,:,:,:)
      call MPI_REDUCE(send_spec, recv_spec, spec_nx * spec_ny * spec_nz * spec_num, MPI_REAL,&
                    & MPI_SUM, root_rnk, MPI_COMM_WORLD, ierr)
      if (mpi_rank .eq. root_rnk) then
        glob_spectra(s,:,:,:,:) = recv_spec(:,:,:,:)
      end if

    end do

    if (allocated(spectra)) deallocate(spectra)
    if (allocated(send_spec)) deallocate(send_spec)
    if (allocated(recv_spec)) deallocate(recv_spec)
  end subroutine prepareSpectraForOutput

  subroutine defineFieldVarsToOutput()
    implicit none
    integer     :: s
    ! initialize field variables
    !   total number of fields (excluding particle densities)
    n_fld_vars = 0
    do s = 1, nspec
      fld_vars(n_fld_vars + 1) = 'dens' // STR(s)
      fld_vars(n_fld_vars + 2) = 'enrg' // STR(s)
      n_fld_vars = n_fld_vars + 2
      if (momenta_enable) then
        fld_vars(n_fld_vars + 1) = 'momX' // STR(s)
        fld_vars(n_fld_vars + 2) = 'momY' // STR(s)
        fld_vars(n_fld_vars + 3) = 'momZ' // STR(s)
        n_fld_vars = n_fld_vars + 3
      end if
      if (npart_enable) then
        fld_vars(n_fld_vars + 1) = 'nprt' // STR(s)
        n_fld_vars = n_fld_vars + 1
      end if
    end do

    fld_vars(n_fld_vars + 1 : n_fld_vars + 1 + 12) =&
                                 & (/'ex   ', 'ey   ', 'ez   ',&
                                   & 'bx   ', 'by   ', 'bz   ',&
                                   & 'jx   ', 'jy   ', 'jz   ',&
                                   & 'xx   ', 'yy   ', 'zz   '/)
    n_fld_vars = n_fld_vars + 12
    if (derivatives_enable) then
      fld_vars(n_fld_vars + 1 : n_fld_vars + 1 + 4) = (/'curlBx', 'curlBy', 'curlBz', 'divE'/)
      n_fld_vars = n_fld_vars + 4
    end if
  end subroutine defineFieldVarsToOutput

  ! writes a field specified by `fld_var` from gridcell `i,j,k` ...
  ! ... to `sm_arr(i1, j1, k1)` with proper interpolation etc for the output
  subroutine selectFieldForOutput(fld_var, i1, j1, k1, i, j, k, writing_lgarrQ)
    implicit none
    character(len=STR_MAX), intent(in)  :: fld_var
    integer(kind=2), intent(in)         :: i1, j1, k1, i, j, k
    logical, intent(in)                 :: writing_lgarrQ
    real                                :: ex0, ey0, ez0, bx0, by0, bz0, jx0, jy0, jz0
    real                                :: dx1, dx2, dy1, dy2, dz1, dz2, divE
    select case (trim(fld_var))
    case('ex')
      #ifndef DEBUG
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
      #else
        ex0 = ex(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = ex0 * B_norm
    case('ey')
      #ifndef DEBUG
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
      #else
        ey0 = ey(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = ey0 * B_norm
    case('ez')
      #ifndef DEBUG
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
      #else
        ez0 = ez(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = ez0 * B_norm
    case('bx')
      #ifndef DEBUG
        call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
      #else
        bx0 = bx(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = bx0 * B_norm
    case('by')
      #ifndef DEBUG
        call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
      #else
        by0 = by(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = by0 * B_norm
    case('bz')
      #ifndef DEBUG
        call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
      #else
        bz0 = bz(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = bz0 * B_norm
    case('jx')
      #ifndef DEBUG
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
      #else
        jx0 = jx(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = -jx0 * B_norm
    case('jy')
      #ifndef DEBUG
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
      #else
        jy0 = jy(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = -jy0 * B_norm
    case('jz')
      #ifndef DEBUG
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
      #else
        jz0 = jz(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = -jz0 * B_norm
    case('curlBx')
      #ifdef oneD
        dx1 = 0.0; dx2 = 0.0
      #elif twoD
        dx1 = (bz(    i,    j,k) - bz(    i,j - 1,    k))
        dx2 = (bz(i - 1,    j,k) - bz(i - 1,j - 1,    k))
      #elif threeD
        dx1 = (bz(    i,    j,    k) - bz(    i,j - 1,    k)) - (by(    i,    j,    k) - by(    i,    j,k - 1))
        dx2 = (bz(i - 1,    j,    k) - bz(i - 1,j - 1,    k)) - (by(i - 1,    j,    k) - by(i - 1,    j,k - 1))
      #endif
      sm_arr(i1, j1, k1) = B_norm * 0.5 * (dx1 + dx2)
    case('curlBy')
      #ifdef oneD
        dy1 = -(bz(i, j, k) - bz(i - 1, j, k))
        dy2 = dy1
      #elif twoD
        dy1 = -(bz(   i,    j,    k) - bz(i - 1,    j,    k))
        dy2 = -(bz(   i,j - 1,    k) - bz(i - 1,j - 1,    k))
      #elif threeD
        dy1 = (bx(    i,    j,    k) - bx(    i,    j,k - 1)) - (bz(    i,    j,    k) - bz(i - 1,    j,    k))
        dy2 = (bx(    i,j - 1,    k) - bx(    i,j - 1,k - 1)) - (bz(    i,j - 1,    k) - bz(i - 1,j - 1,    k))
      #endif
      sm_arr(i1, j1, k1) = B_norm * 0.5 * (dy1 + dy2)
    case('curlBz')
      #ifdef oneD
        dz1 = (by(i, j, k) - by(i - 1, j, k))
        dz2 = dz1
      #elif twoD
        dz1 = (by(i, j, k) - by(i - 1, j, k)) - (bx(i, j, k) - bx(i, j - 1, k))
        dz2 = dz1
      #elif threeD
        dz1 = (by(    i,    j,    k) - by(i - 1,    j,    k)) - (bx(    i,    j,    k) - bx(    i,j - 1,    k))
        dz2 = (by(    i,    j,k - 1) - by(i - 1,    j,k - 1)) - (bx(    i,    j,k - 1) - bx(    i,j - 1,k - 1))
      #endif
      sm_arr(i1, j1, k1) = B_norm * 0.5 * (dz1 + dz2)
    case('divE')
      divE = 0.0
      #if defined(oneD) || defined (twoD) || defined (threeD)
        divE = divE + (ex(i, j, k) - ex(i - 1, j, k))
      #endif
      #if defined(twoD) || defined (threeD)
        divE = divE + (ey(i, j, k) - ey(i, j - 1, k))
      #endif
      #if defined(threeD)
        divE = divE + (ez(i, j, k) - ez(i, j, k - 1))
      #endif
      sm_arr(i1, j1, k1) = B_norm * divE
    case('xx')
      sm_arr(i1, j1, k1) = REAL(this_meshblock%ptr%x0 + i, 4)
    case('yy')
      sm_arr(i1, j1, k1) = REAL(this_meshblock%ptr%y0 + j, 4)
    case('zz')
      sm_arr(i1, j1, k1) = REAL(this_meshblock%ptr%z0 + k, 4)
    case default
      if (((fld_var(1:4) .ne. 'dens') .and.&
         & (fld_var(1:4) .ne. 'enrg') .and.&
         & (fld_var(1:3) .ne. 'mom') .and.&
         & (fld_var(1:4) .ne. 'nprt')) .and.&
         & (.not. writing_lgarrQ)) then
        call throwError("ERROR: unrecognized `fldname`: " // trim(fld_var))
      else
        ! interpolating cell-centered values to nodes (when smoothing is large enough)
        if (output_dens_smooth .gt. 1) then
          #ifdef oneD
            sm_arr(i1, j1, k1) = 0.5 * (lg_arr(i, j, k) + lg_arr(i - 1, j, k))
          #elif twoD
            sm_arr(i1, j1, k1) = 0.25 * (lg_arr(i, j, k) + lg_arr(i - 1, j, k) +&
                                       & lg_arr(i - 1, j - 1, k) + lg_arr(i, j - 1, k))
          #elif threeD
            sm_arr(i1, j1, k1) = 0.125 * (lg_arr(i, j, k) + lg_arr(i - 1, j - 1, k - 1) +&
                                        & lg_arr(i - 1, j, k) + lg_arr(i, j - 1, k) + lg_arr(i, j, k - 1) +&
                                        & lg_arr(i - 1, j - 1, k) + lg_arr(i, j - 1, k - 1) + lg_arr(i - 1, j, k - 1))
          #endif
        else
          sm_arr(i1, j1, k1) = lg_arr(i, j, k)
        end if
      end if
    end select
  end subroutine selectFieldForOutput

  subroutine prepareFieldForOutput(fldname, writing_lgarrQ)
    implicit none
    character(len=STR_MAX), intent(in)    :: fldname
    logical, intent(out)                  :: writing_lgarrQ
    integer                               :: s

    if (fldname(1:4) .eq. 'dens') then
      writing_lgarrQ = .true.
      s = STRtoINT(fldname(5:5))
      ! fill `lg_arr` with density of species `s`
      call computeDensity(s, reset=.true., ds=output_dens_smooth)
      call exchangeArray()
    else if (fldname(1:4) .eq. 'enrg') then
      writing_lgarrQ = .true.
      s = STRtoINT(fldname(5:5))
      ! fill `lg_arr` with energy density of species `s`
      call computeMomentum(s, 0, reset=.true., ds=output_dens_smooth)
      call exchangeArray()
    else if (fldname(1:4) .eq. 'nprt') then
      writing_lgarrQ = .true.
      s = STRtoINT(fldname(5:5))
      call computeNpart(s, reset=.true., ds=output_dens_smooth)
      call exchangeArray()
    else if (fldname(1:4) .eq. 'momX') then
      writing_lgarrQ = .true.
      s = STRtoINT(fldname(5:5))
      call computeMomentum(s, 1, reset=.true., ds=output_dens_smooth)
      call exchangeArray()
    else if (fldname(1:4) .eq. 'momY') then
      writing_lgarrQ = .true.
      s = STRtoINT(fldname(5:5))
      call computeMomentum(s, 2, reset=.true., ds=output_dens_smooth)
      call exchangeArray()
    else if (fldname(1:4) .eq. 'momZ') then
      writing_lgarrQ = .true.
      s = STRtoINT(fldname(5:5))
      call computeMomentum(s, 3, reset=.true., ds=output_dens_smooth)
      call exchangeArray()
    else
      writing_lgarrQ = .false.
    end if
  end subroutine prepareFieldForOutput

end module m_outputlogistics
