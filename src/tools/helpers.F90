module m_helpers
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_qednamespace
  implicit none
contains
  logical function numbersAreClose(number1, number2)
    implicit none
    real, intent(in) :: number1, number2
    real :: abs1, abs2
    real :: diff

    abs1 = abs(number1); abs2 = abs(number2)
    diff = abs(number1 - number2)

    if ((number1 .le. number2) .and. (number1 .ge. number2)) then
      numbersAreClose = .true.
    else if ((number1 .le. 0.0) .or. (number2 .le. 0.0) .or. &
             (abs1 + abs2 .lt. TINYREAL)) then
      numbersAreClose = (diff .lt. TINYREAL)
    else
      numbersAreClose = (diff / (abs1 + abs2) .lt. 0.5 * TINYREAL)
    end if
  end function numbersAreClose

  subroutine checkNpart(msg)
    implicit none
    integer :: ti, tj, tk, s, nprt, ierr
    character(len=*), intent(in) :: msg
    integer :: glob_nprt
    do s = 1, nspec
      nprt = 0
      do tk = 1, species(s) % tile_nz
        do tj = 1, species(s) % tile_ny
          do ti = 1, species(s) % tile_nx
            nprt = nprt + species(s) % prtl_tile(ti, tj, tk) % npart_sp
          end do
        end do
      end do
      call MPI_REDUCE(nprt, glob_nprt, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (mpi_rank .eq. 0) then
        print *, "npart of", s, msg, ":", glob_nprt
      end if
    end do
  end subroutine checkNpart

  subroutine globalToLocalCoords(x_glob, y_glob, z_glob, &
                                 x_loc, y_loc, z_loc, adjustQ, containedQ)
    implicit none
    real, intent(in) :: x_glob, y_glob, z_glob
    real, intent(out) :: x_loc, y_loc, z_loc
    logical, optional, intent(in) :: adjustQ
    logical, optional, intent(out) :: containedQ
    logical :: adjustQ_
    if (present(adjustQ)) then
      adjustQ_ = adjustQ
    else
      adjustQ_ = .false.
    end if

    if (adjustQ_) then
      x_loc = x_glob; y_loc = y_glob; z_loc = z_glob
#if defined(oneD) || defined (twoD) || defined (threeD)
      x_loc = MAX(0.0, MIN(x_glob - REAL(this_meshblock % ptr % x0), REAL(this_meshblock % ptr % sx)))
#endif
#if defined (twoD) || defined (threeD)
      y_loc = MAX(0.0, MIN(y_glob - REAL(this_meshblock % ptr % y0), REAL(this_meshblock % ptr % sy)))
#endif
#if defined(threeD)
      z_loc = MAX(0.0, MIN(z_glob - REAL(this_meshblock % ptr % z0), REAL(this_meshblock % ptr % sz)))
#endif
    else
      x_loc = x_glob; y_loc = y_glob; z_loc = z_glob
#if defined(oneD) || defined (twoD) || defined (threeD)
      x_loc = x_glob - REAL(this_meshblock % ptr % x0)
#endif
#if defined (twoD) || defined (threeD)
      y_loc = y_glob - REAL(this_meshblock % ptr % y0)
#endif
#if defined(threeD)
      z_loc = z_glob - REAL(this_meshblock % ptr % z0)
#endif
      if (present(containedQ)) then
#ifdef oneD
        containedQ = ((x_loc .ge. 0) .and. &
                      (x_loc .lt. REAL(this_meshblock % ptr % sx)))
#elif defined(twoD)
        containedQ = ((x_loc .ge. 0) .and. &
                      (x_loc .lt. REAL(this_meshblock % ptr % sx)) .and. &
                      (y_loc .ge. 0) .and. &
                      (y_loc .lt. REAL(this_meshblock % ptr % sy)))
#elif defined(threeD)
        containedQ = ((x_loc .ge. 0) .and. &
                      (x_loc .lt. REAL(this_meshblock % ptr % sx)) .and. &
                      (y_loc .ge. 0) .and. &
                      (y_loc .lt. REAL(this_meshblock % ptr % sy)) .and. &
                      (z_loc .ge. 0) .and. &
                      (z_loc .lt. REAL(this_meshblock % ptr % sz)))
#endif
      end if
    end if
  end subroutine globalToLocalCoords

  subroutine localToCellBasedCoords(x_loc, y_loc, z_loc, &
                                    xi, yi, zi, dx, dy, dz)
    implicit none
    real, intent(in) :: x_loc, y_loc, z_loc
    real, intent(out) :: dx, dy, dz
    integer(kind=2), intent(out) :: xi, yi, zi

    xi = INT(FLOOR(x_loc), 2); dx = x_loc - FLOOR(x_loc)
    yi = INT(FLOOR(y_loc), 2); dy = y_loc - FLOOR(y_loc)
    zi = INT(FLOOR(z_loc), 2); dz = z_loc - FLOOR(z_loc)
  end subroutine localToCellBasedCoords

  subroutine generateCoordInRegion(xmin, xmax, ymin, ymax, zmin, zmax, &
                                   x_, y_, z_, xi_, yi_, zi_, dx_, dy_, dz_)
    implicit none
    real :: rnd
    real, intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax
    real, intent(out) :: x_, y_, z_, dx_, dy_, dz_
    integer(kind=2), intent(out) :: xi_, yi_, zi_

    x_ = 0.5; xi_ = 0; dx_ = 0.5
    y_ = 0.5; yi_ = 0; dy_ = 0.5
    z_ = 0.5; zi_ = 0; dz_ = 0.5
#if defined(oneD) || defined (twoD) || defined (threeD)
    rnd = random(dseed)
    x_ = xmin + rnd * (xmax - xmin)
    xi_ = INT(FLOOR(x_), 2); dx_ = x_ - FLOOR(x_)
    if (xi_ .eq. this_meshblock % ptr % sx) then
      xi_ = xi_ - 1_2; dx_ = dx_ + 1.0
    end if
#endif
#if defined (twoD) || defined (threeD)
    rnd = random(dseed)
    y_ = ymin + rnd * (ymax - ymin)
    yi_ = INT(FLOOR(y_), 2); dy_ = y_ - FLOOR(y_)
    if (yi_ .eq. this_meshblock % ptr % sy) then
      yi_ = yi_ - 1_2; dy_ = dy_ + 1.0
    end if
#endif
#if defined(threeD)
    rnd = random(dseed)
    z_ = zmin + rnd * (zmax - zmin)
    zi_ = INT(FLOOR(z_), 2); dz_ = z_ - FLOOR(z_)
    if (zi_ .eq. this_meshblock % ptr % sz) then
      zi_ = zi_ - 1_2; dz_ = dz_ + 1.0
    end if
#endif
  end subroutine generateCoordInRegion

  function rnkToInd(rnk)
    implicit none
    integer, intent(in) :: rnk
    integer, dimension(3) :: rnkToInd
    if ((rnk .lt. 0) .or. (rnk .ge. mpi_size)) then
      rnkToInd = (/-1, -1, -1/)
    else
      rnkToInd(3) = rnk / (sizex * sizey)
      rnkToInd(2) = (rnk - sizex * sizey * rnkToInd(3)) / sizex
      rnkToInd(1) = rnk - sizex * sizey * rnkToInd(3) - sizex * rnkToInd(2)
    end if
  end function rnkToInd

  function indToRnk(ind)
    implicit none
    integer, intent(in) :: ind(3)
    integer :: ind_t(3), indToRnk
    ind_t = ind
    if (boundary_x .eq. 1) then
      if (ind_t(1) .lt. 0) ind_t(1) = ind_t(1) + sizex
      ind_t(1) = modulo(ind_t(1), sizex)
    end if
    if (boundary_y .eq. 1) then
      if (ind_t(2) .lt. 0) ind_t(2) = ind_t(2) + sizey
      ind_t(2) = modulo(ind_t(2), sizey)
    end if
    if (boundary_z .eq. 1) then
      if (ind_t(3) .lt. 0) ind_t(3) = ind_t(3) + sizez
      ind_t(3) = modulo(ind_t(3), sizez)
    end if
    if ((ind_t(1) .lt. 0) .or. (ind_t(1) .ge. sizex) .or. &
        (ind_t(2) .lt. 0) .or. (ind_t(2) .ge. sizey) .or. &
        (ind_t(3) .lt. 0) .or. (ind_t(3) .ge. sizez)) then
      indToRnk = -1
    else
      indToRnk = ind_t(3) * sizex * sizey + ind_t(2) * sizex + ind_t(1)
    end if
  end function indToRnk

  subroutine reassignNeighborsForAll(mblocks)
    implicit none
    type(mesh), allocatable, target, intent(inout) :: mblocks(:)
    integer :: rnk
    integer :: ind1, ind2, ind3, inds(3)
    do rnk = 0, mpi_size - 1
      do ind3 = -1, 1
        do ind2 = -1, 1
          do ind1 = -1, 1
            inds(1) = ind1; inds(2) = ind2; inds(3) = ind3
            call assignNeighbor(rnk, inds, mblocks)
          end do
        end do
      end do
    end do
    call computeNumberOfNeighbors()
  end subroutine reassignNeighborsForAll

  subroutine assignNeighbor(rnk, inds1, mblocks)
    implicit none
    type(mesh), allocatable, target, intent(inout) :: mblocks(:)
    integer, intent(in) :: rnk, inds1(3)
    integer :: rnk2, inds0(3)
    integer :: inds_temp(3), i
    inds0 = rnkToInd(rnk)
    do i = 1, 3
      inds_temp(i) = inds0(i) + inds1(i)
    end do
    rnk2 = indToRnk(inds_temp)
    if (rnk2 .eq. -1) then
      mblocks(rnk + 1) % neighbor(inds1(1), inds1(2), inds1(3)) % ptr => null()
    else
      mblocks(rnk + 1) % neighbor(inds1(1), inds1(2), inds1(3)) % ptr => mblocks(rnk2 + 1)
    end if
  end subroutine assignNeighbor

  subroutine computeNumberOfNeighbors()
    integer :: ind1, ind2, ind3
    integer :: cntr
    cntr = 0
    do ind3 = -1, 1
      do ind2 = -1, 1
        do ind1 = -1, 1
          if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
#ifdef oneD
          if ((ind2 .ne. 0) .or. (ind3 .ne. 0)) cycle
#elif defined(twoD)
          if (ind3 .ne. 0) cycle
#endif
          if (.not. associated(this_meshblock % ptr % neighbor(ind1, ind2, ind3) % ptr)) cycle
          cntr = cntr + 1
        end do
      end do
    end do
    sendrecv_neighbors = cntr
  end subroutine computeNumberOfNeighbors

  subroutine computeDensity(s, reset, ds, charge)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in) :: s
    logical, intent(in) :: reset
    logical, optional, intent(in) :: charge
    integer, optional, intent(in) :: ds
    integer :: p, ti, tj, tk
    integer(kind=2), pointer, contiguous :: pt_xi(:), pt_yi(:), pt_zi(:)
    real, pointer, contiguous :: pt_wei(:)
    logical :: charge_
    integer(kind=2) :: i, j, k, i1, i2, j1, j2, k1, k2, ds_, pow
    real :: contrib

    if (.not. present(ds)) then
      ds_ = 2_2
    else
      ds_ = INT(ds, 2)
    end if

    if (.not. present(charge)) then
      charge_ = .false.
    else
      charge_ = charge
    end if

#ifdef oneD
    pow = 1_2
#elif defined(twoD)
    pow = 2_2
#elif defined(threeD)
    pow = 3_2
#endif

    if (species(s) % m_sp .le. 0) then
      contrib = 1.0 / (2.0 * REAL(ds_) + 1.0)**pow
    else
      if (charge_) then
        contrib = species(s) % ch_sp / (2.0 * REAL(ds_) + 1.0)**pow
      else
        contrib = species(s) % m_sp / (2.0 * REAL(ds_) + 1.0)**pow
      end if
    end if

    if (reset) then
      lg_arr(:, :, :) = 0
    end if

    do tk = 1, species(s) % tile_nz
      do tj = 1, species(s) % tile_ny
        do ti = 1, species(s) % tile_nx
          pt_xi => species(s) % prtl_tile(ti, tj, tk) % xi
          pt_yi => species(s) % prtl_tile(ti, tj, tk) % yi
          pt_zi => species(s) % prtl_tile(ti, tj, tk) % zi
          pt_wei => species(s) % prtl_tile(ti, tj, tk) % weight
          do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
            i = pt_xi(p); j = pt_yi(p); k = pt_zi(p)

            i1 = 0; i2 = 0
            j1 = 0; j2 = 0
            k1 = 0; k2 = 0
#if defined(oneD) || defined (twoD) || defined (threeD)
            i1 = max(i - ds_, -NGHOST)
            i2 = min(i + ds_, INT(this_meshblock % ptr % sx + NGHOST - 1, 2))
#endif
#if defined (twoD) || defined (threeD)
            j1 = max(j - ds_, -NGHOST)
            j2 = min(j + ds_, INT(this_meshblock % ptr % sy + NGHOST - 1, 2))
#endif
#if defined (threeD)
            k1 = max(k - ds_, -NGHOST)
            k2 = min(k + ds_, INT(this_meshblock % ptr % sz + NGHOST - 1, 2))
#endif

            do k = k1, k2
              do j = j1, j2
                do i = i1, i2
                  lg_arr(i, j, k) = lg_arr(i, j, k) + pt_wei(p) * contrib
                end do
              end do
            end do

          end do
          pt_xi => null(); pt_yi => null(); pt_zi => null(); pt_wei => null()
        end do
      end do
    end do
  end subroutine computeDensity

  subroutine computePrtCurr(s, component, reset, ds)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in) :: s, component
    logical, intent(in) :: reset
    integer, optional, intent(in) :: ds
    integer :: p, ti, tj, tk
    integer(kind=2), pointer, contiguous :: pt_xi(:), pt_yi(:), pt_zi(:)
    real, pointer, contiguous :: pt_u(:), pt_v(:), pt_w(:), pt_wei(:)
    integer(kind=2) :: i, j, k
    integer :: i1, i2, j1, j2, k1, k2, ds_
    integer :: pow
    logical :: charged
    real :: comp
    real :: contrib

    if (.not. present(ds)) then
      ds_ = 2
    else
      ds_ = ds
    end if

#ifdef oneD
    pow = 1
#elif defined(twoD)
    pow = 2
#elif defined(threeD)
    pow = 3
#endif

    contrib = species(s) % ch_sp / (2.0 * REAL(ds_) + 1.0)**pow

    if (reset) then
      lg_arr(:, :, :) = 0
    end if
    do tk = 1, species(s) % tile_nz
      do tj = 1, species(s) % tile_ny
        do ti = 1, species(s) % tile_nx
          pt_xi => species(s) % prtl_tile(ti, tj, tk) % xi
          pt_yi => species(s) % prtl_tile(ti, tj, tk) % yi
          pt_zi => species(s) % prtl_tile(ti, tj, tk) % zi
          pt_wei => species(s) % prtl_tile(ti, tj, tk) % weight
          pt_u => species(s) % prtl_tile(ti, tj, tk) % u
          pt_v => species(s) % prtl_tile(ti, tj, tk) % v
          pt_w => species(s) % prtl_tile(ti, tj, tk) % w
          do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
            i = pt_xi(p); j = pt_yi(p); k = pt_zi(p)
            if (component .eq. 1) then
              comp = pt_u(p) / sqrt(1.0 + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
            else if (component .eq. 2) then
              comp = pt_v(p) / sqrt(1.0 + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
            else if (component .eq. 3) then
              comp = pt_w(p) / sqrt(1.0 + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
            end if

            i1 = 0; i2 = 0
            j1 = 0; j2 = 0
            k1 = 0; k2 = 0
#if defined(oneD) || defined (twoD) || defined (threeD)
            i1 = max(i - ds_, -NGHOST)
            i2 = min(i + ds_, this_meshblock % ptr % sx + NGHOST - 1)
#endif
#if defined (twoD) || defined (threeD)
            j1 = max(j - ds_, -NGHOST)
            j2 = min(j + ds_, this_meshblock % ptr % sy + NGHOST - 1)
#endif
#if defined (threeD)
            k1 = max(k - ds_, -NGHOST)
            k2 = min(k + ds_, this_meshblock % ptr % sz + NGHOST - 1)
#endif

            do k = k1, k2
              do j = j1, j2
                do i = i1, i2
                  lg_arr(i, j, k) = lg_arr(i, j, k) + comp * pt_wei(p) * contrib
                end do
              end do
            end do

          end do
          pt_xi => null(); pt_yi => null(); pt_zi => null()
          pt_u => null(); pt_v => null(); pt_w => null()
          pt_wei => null()
        end do
      end do
    end do
  end subroutine computePrtCurr

  subroutine computeMomentum(s, component, reset, ds)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in) :: s, component
    logical, intent(in) :: reset
    integer, optional, intent(in) :: ds
    integer :: p, ti, tj, tk
    integer(kind=2), pointer, contiguous :: pt_xi(:), pt_yi(:), pt_zi(:)
    real, pointer, contiguous :: pt_u(:), pt_v(:), pt_w(:), pt_wei(:)
    integer(kind=2) :: i, j, k, i1, i2, j1, j2, k1, k2, ds_, pow
    logical :: massive
    real :: comp = 0.0
    real :: contrib

    if (.not. present(ds)) then
      ds_ = 2_2
    else
      ds_ = INT(ds, 2)
    end if

#ifdef oneD
    pow = 1_2
#elif defined(twoD)
    pow = 2_2
#elif defined(threeD)
    pow = 3_2
#endif

    massive = (species(s) % m_sp .gt. 0)
    if (.not. massive) then
      contrib = 1.0 / (2.0 * REAL(ds_) + 1.0)**pow
    else
      contrib = species(s) % m_sp / (2.0 * REAL(ds_) + 1.0)**pow
    end if

    if (reset) then
      lg_arr(:, :, :) = 0
    end if
    do tk = 1, species(s) % tile_nz
      do tj = 1, species(s) % tile_ny
        do ti = 1, species(s) % tile_nx
          pt_xi => species(s) % prtl_tile(ti, tj, tk) % xi
          pt_yi => species(s) % prtl_tile(ti, tj, tk) % yi
          pt_zi => species(s) % prtl_tile(ti, tj, tk) % zi
          pt_wei => species(s) % prtl_tile(ti, tj, tk) % weight
          pt_u => species(s) % prtl_tile(ti, tj, tk) % u
          pt_v => species(s) % prtl_tile(ti, tj, tk) % v
          pt_w => species(s) % prtl_tile(ti, tj, tk) % w
          do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
            i = pt_xi(p); j = pt_yi(p); k = pt_zi(p)
            if (component .eq. 0) then
              if (massive) then
                comp = sqrt(1.0 + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
              else
                comp = sqrt(pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
              end if
            else if (component .eq. 1) then
              comp = pt_u(p)
            else if (component .eq. 2) then
              comp = pt_v(p)
            else if (component .eq. 3) then
              comp = pt_w(p)
            else if (component .eq. 4) then
              comp = pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2  ! momentum squared
              if (massive) comp = comp * species(s) % m_sp
            end if

            i1 = 0; i2 = 0
            j1 = 0; j2 = 0
            k1 = 0; k2 = 0
#if defined(oneD) || defined (twoD) || defined (threeD)
            i1 = max(i - ds_, -NGHOST)
            i2 = min(i + ds_, INT(this_meshblock % ptr % sx + NGHOST - 1, 2))
#endif
#if defined (twoD) || defined (threeD)
            j1 = max(j - ds_, -NGHOST)
            j2 = min(j + ds_, INT(this_meshblock % ptr % sy + NGHOST - 1, 2))
#endif
#if defined (threeD)
            k1 = max(k - ds_, -NGHOST)
            k2 = min(k + ds_, INT(this_meshblock % ptr % sz + NGHOST - 1, 2))
#endif

            do k = k1, k2
              do j = j1, j2
                do i = i1, i2
                  lg_arr(i, j, k) = lg_arr(i, j, k) + comp * pt_wei(p) * contrib
                end do
              end do
            end do

          end do
          pt_xi => null(); pt_yi => null(); pt_zi => null()
          pt_u => null(); pt_v => null(); pt_w => null()
          pt_wei => null()
        end do
      end do
    end do
  end subroutine computeMomentum

  subroutine computeEnergyMomentum(s, component1, component2, reset, ds)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in) :: s, component1, component2
    logical, intent(in) :: reset
    integer, optional, intent(in) :: ds
    integer :: p, ti, tj, tk
    integer(kind=2), pointer, contiguous :: pt_xi(:), pt_yi(:), pt_zi(:)
    real, pointer, contiguous :: pt_u(:), pt_v(:), pt_w(:), pt_wei(:)
    integer(kind=2) :: i, j, k
    integer(kind=2) :: i1, i2, j1, j2, k1, k2, ds_
    integer :: pow
    logical :: massive
    real :: comp1, comp2
    real :: contrib

    real :: p0inv

    if (.not. present(ds)) then
      ds_ = 2_2
    else
      ds_ = INT(ds, 2)
    end if

#ifdef oneD
    pow = 1
#elif defined(twoD)
    pow = 2
#elif defined(threeD)
    pow = 3
#endif

    massive = (species(s) % m_sp .ne. 0)
    if (.not. massive) then
      contrib = 1.0 / (2.0 * REAL(ds_) + 1.0)**pow
    else
      contrib = species(s) % m_sp / (2.0 * REAL(ds_) + 1.0)**pow
    end if

    if (reset) then
      lg_arr(:, :, :) = 0
    end if
    do tk = 1, species(s) % tile_nz
      do tj = 1, species(s) % tile_ny
        do ti = 1, species(s) % tile_nx
          pt_xi => species(s) % prtl_tile(ti, tj, tk) % xi
          pt_yi => species(s) % prtl_tile(ti, tj, tk) % yi
          pt_zi => species(s) % prtl_tile(ti, tj, tk) % zi
          pt_wei => species(s) % prtl_tile(ti, tj, tk) % weight
          pt_u => species(s) % prtl_tile(ti, tj, tk) % u
          pt_v => species(s) % prtl_tile(ti, tj, tk) % v
          pt_w => species(s) % prtl_tile(ti, tj, tk) % w
          do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
            i = pt_xi(p); j = pt_yi(p); k = pt_zi(p)

            if (massive) then
              p0inv = 1.0 / sqrt(1.0 + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
            else
              p0inv = 1.0 / sqrt(pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
            end if

            if (component1 .eq. 0) then
              comp1 = 1.0 / p0inv
            else if (component1 .eq. 1) then
              comp1 = pt_u(p)
            else if (component1 .eq. 2) then
              comp1 = pt_v(p)
            else if (component1 .eq. 3) then
              comp1 = pt_w(p)
            end if

            if (component2 .eq. 0) then
              comp2 = 1.0 / p0inv
            else if (component2 .eq. 1) then
              comp2 = pt_u(p)
            else if (component2 .eq. 2) then
              comp2 = pt_v(p)
            else if (component2 .eq. 3) then
              comp2 = pt_w(p)
            end if

            i1 = 0; i2 = 0
            j1 = 0; j2 = 0
            k1 = 0; k2 = 0
#if defined(oneD) || defined (twoD) || defined (threeD)
            i1 = max(i - ds_, -NGHOST)
            i2 = min(i + ds_, this_meshblock % ptr % sx + NGHOST - 1)
#endif
#if defined (twoD) || defined (threeD)
            j1 = max(j - ds_, -NGHOST)
            j2 = min(j + ds_, INT(this_meshblock % ptr % sy + NGHOST - 1, 2))
#endif
#if defined (threeD)
            k1 = max(k - ds_, -NGHOST)
            k2 = min(k + ds_, INT(this_meshblock % ptr % sz + NGHOST - 1, 2))
#endif

            do k = k1, k2
              do j = j1, j2
                do i = i1, i2
                  lg_arr(i, j, k) = lg_arr(i, j, k) + pt_wei(p) * contrib * comp1 * comp2 * p0inv
                end do
              end do
            end do

          end do
          pt_xi => null(); pt_yi => null(); pt_zi => null()
          pt_u => null(); pt_v => null(); pt_w => null()
          pt_wei => null()
        end do
      end do
    end do
  end subroutine computeEnergyMomentum

  subroutine computeFluidVelocity(component, ds)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in) :: component
    integer, optional, intent(in) :: ds
    integer :: s
    integer :: p, ti, tj, tk
    integer(kind=2), pointer, contiguous :: pt_xi(:), pt_yi(:), pt_zi(:)
    real, pointer, contiguous :: pt_u(:), pt_v(:), pt_w(:), pt_wei(:)
    integer(kind=2) :: pow, ds_, i, j, k, i1, i2, j1, j2, k1, k2
    real :: contrib, inv_gamma, comp = 0.0

    if (.not. present(ds)) then
      ds_ = 2_2
    else
      ds_ = INT(ds, 2)
    end if

#ifdef oneD
    pow = 1_2
#elif defined(twoD)
    pow = 2_2
#elif defined(threeD)
    pow = 3_2
#endif

    lg_arr(:, :, :) = 0.0
    do s = 1, nspec
      if (.not. species(s) % move_sp) cycle
      if (.not. species(s) % deposit_sp) cycle
      if (species(s) % m_sp .le. 0.0) cycle

      ! this factor takes into account smoothing and mass of the species
      contrib = species(s) % m_sp / (2.0 * REAL(ds_) + 1.0)**pow

      do tk = 1, species(s) % tile_nz
        do tj = 1, species(s) % tile_ny
          do ti = 1, species(s) % tile_nx
            pt_xi => species(s) % prtl_tile(ti, tj, tk) % xi
            pt_yi => species(s) % prtl_tile(ti, tj, tk) % yi
            pt_zi => species(s) % prtl_tile(ti, tj, tk) % zi
            pt_wei => species(s) % prtl_tile(ti, tj, tk) % weight
            pt_u => species(s) % prtl_tile(ti, tj, tk) % u
            pt_v => species(s) % prtl_tile(ti, tj, tk) % v
            pt_w => species(s) % prtl_tile(ti, tj, tk) % w
            do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
              ! compute 3-velocity
              inv_gamma = 1.0 / sqrt(1.0 + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
              if (component .eq. 0) then
                comp = pt_u(p) * inv_gamma
              else if (component .eq. 1) then
                comp = pt_v(p) * inv_gamma
              else if (component .eq. 2) then
                comp = pt_w(p) * inv_gamma
              end if

              i = pt_xi(p); j = pt_yi(p); k = pt_zi(p)

              i1 = 0; i2 = 0
              j1 = 0; j2 = 0
              k1 = 0; k2 = 0
#if defined(oneD) || defined (twoD) || defined (threeD)
              i1 = max(i - ds_, -NGHOST)
              i2 = min(i + ds_, INT(this_meshblock % ptr % sx + NGHOST - 1, 2))
#endif
#if defined (twoD) || defined (threeD)
              j1 = max(j - ds_, -NGHOST)
              j2 = min(j + ds_, INT(this_meshblock % ptr % sy + NGHOST - 1, 2))
#endif
#if defined (threeD)
              k1 = max(k - ds_, -NGHOST)
              k2 = min(k + ds_, INT(this_meshblock % ptr % sz + NGHOST - 1, 2))
#endif

              do k = k1, k2
                do j = j1, j2
                  do i = i1, i2
                    lg_arr(i, j, k) = lg_arr(i, j, k) + comp * pt_wei(p) * contrib
                  end do
                end do
              end do
            end do ! particles on tile
            pt_xi => null(); pt_yi => null(); pt_zi => null()
            pt_u => null(); pt_v => null(); pt_w => null()
            pt_wei => null()
          end do
        end do
      end do ! tiles
    end do ! species
  end subroutine computeFluidVelocity

  subroutine computeFluidDensity(ds)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, optional, intent(in) :: ds
    integer :: s
    integer :: p, ti, tj, tk
    integer(kind=2), pointer, contiguous :: pt_xi(:), pt_yi(:), pt_zi(:)
    real, pointer, contiguous :: pt_wei(:)
    integer(kind=2) :: pow, ds_, i, j, k, i1, i2, j1, j2, k1, k2
    real :: contrib

    if (.not. present(ds)) then
      ds_ = 2_2
    else
      ds_ = INT(ds, 2)
    end if

#ifdef oneD
    pow = 1_2
#elif defined(twoD)
    pow = 2_2
#elif defined(threeD)
    pow = 3_2
#endif

    lg_arr(:, :, :) = 0.0
    do s = 1, nspec
      if (.not. species(s) % move_sp) cycle
      if (.not. species(s) % deposit_sp) cycle
      if (species(s) % m_sp .le. 0) cycle

      ! this factor takes into account smoothing and mass of the species
      contrib = species(s) % m_sp / (2.0 * REAL(ds_) + 1.0)**pow

      do tk = 1, species(s) % tile_nz
        do tj = 1, species(s) % tile_ny
          do ti = 1, species(s) % tile_nx
            pt_xi => species(s) % prtl_tile(ti, tj, tk) % xi
            pt_yi => species(s) % prtl_tile(ti, tj, tk) % yi
            pt_zi => species(s) % prtl_tile(ti, tj, tk) % zi
            pt_wei => species(s) % prtl_tile(ti, tj, tk) % weight
            do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
              ! compute density
              i = pt_xi(p); j = pt_yi(p); k = pt_zi(p)

              i1 = 0; i2 = 0
              j1 = 0; j2 = 0
              k1 = 0; k2 = 0
#if defined(oneD) || defined (twoD) || defined (threeD)
              i1 = max(i - ds_, -NGHOST)
              i2 = min(i + ds_, INT(this_meshblock % ptr % sx + NGHOST - 1, 2))
#endif
#if defined (twoD) || defined (threeD)
              j1 = max(j - ds_, -NGHOST)
              j2 = min(j + ds_, INT(this_meshblock % ptr % sy + NGHOST - 1, 2))
#endif
#if defined (threeD)
              k1 = max(k - ds_, -NGHOST)
              k2 = min(k + ds_, INT(this_meshblock % ptr % sz + NGHOST - 1, 2))
#endif

              do k = k1, k2
                do j = j1, j2
                  do i = i1, i2
                    lg_arr(i, j, k) = lg_arr(i, j, k) + pt_wei(p) * contrib
                  end do
                end do
              end do
            end do ! particles on tile
            pt_xi => null(); pt_yi => null(); pt_zi => null()
            pt_wei => null()
          end do
        end do
      end do ! tiles
    end do ! species
  end subroutine computeFluidDensity

  subroutine computeNpart(s, reset, ds)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in) :: s
    logical, intent(in) :: reset
    integer, optional, intent(in) :: ds
    integer :: p, ti, tj, tk
    integer(kind=2), pointer, contiguous :: pt_xi(:), pt_yi(:), pt_zi(:)
    integer(kind=2) :: i, j, k, i1, i2, j1, j2, k1, k2, ds_, pow
    real :: contrib

    if (.not. present(ds)) then
      ds_ = 2_2
    else
      ds_ = INT(ds, 2)
    end if

#ifdef oneD
    pow = 1_2
#elif defined(twoD)
    pow = 2_2
#elif defined(threeD)
    pow = 3_2
#endif

    contrib = 1.0 / (2.0 * REAL(ds_) + 1.0)**pow

    if (reset) then
      lg_arr(:, :, :) = 0
    end if
    do tk = 1, species(s) % tile_nz
      do tj = 1, species(s) % tile_ny
        do ti = 1, species(s) % tile_nx
          pt_xi => species(s) % prtl_tile(ti, tj, tk) % xi
          pt_yi => species(s) % prtl_tile(ti, tj, tk) % yi
          pt_zi => species(s) % prtl_tile(ti, tj, tk) % zi
          do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
            i = pt_xi(p); j = pt_yi(p); k = pt_zi(p)
            i1 = 0; i2 = 0
            j1 = 0; j2 = 0
            k1 = 0; k2 = 0
#if defined(oneD) || defined (twoD) || defined (threeD)
            i1 = max(i - ds_, -NGHOST)
            i2 = min(i + ds_, INT(this_meshblock % ptr % sx + NGHOST - 1, 2))
#endif
#if defined (twoD) || defined (threeD)
            j1 = max(j - ds_, -NGHOST)
            j2 = min(j + ds_, INT(this_meshblock % ptr % sy + NGHOST - 1, 2))
#endif
#if defined (threeD)
            k1 = max(k - ds_, -NGHOST)
            k2 = min(k + ds_, INT(this_meshblock % ptr % sz + NGHOST - 1, 2))
#endif
            do k = k1, k2
              do j = j1, j2
                do i = i1, i2
                  lg_arr(i, j, k) = lg_arr(i, j, k) + contrib
                end do
              end do
            end do
          end do
          pt_xi => null(); pt_yi => null(); pt_zi => null()
        end do
      end do
    end do
  end subroutine computeNpart

  subroutine interpFromEdges(dx, dy, dz, i, j, k, &
                             fx, fy, fz, &
                             intfx, intfy, intfz, comp)
    implicit none
    integer(kind=2), intent(in) :: i, j, k
    real, intent(in) :: dx, dy, dz
    real, intent(in) :: fx(this_meshblock % ptr % i1:this_meshblock % ptr % i2, &
                           this_meshblock % ptr % j1:this_meshblock % ptr % j2, &
                           this_meshblock % ptr % k1:this_meshblock % ptr % k2)
    real, intent(in) :: fy(this_meshblock % ptr % i1:this_meshblock % ptr % i2, &
                           this_meshblock % ptr % j1:this_meshblock % ptr % j2, &
                           this_meshblock % ptr % k1:this_meshblock % ptr % k2)
    real, intent(in) :: fz(this_meshblock % ptr % i1:this_meshblock % ptr % i2, &
                           this_meshblock % ptr % j1:this_meshblock % ptr % j2, &
                           this_meshblock % ptr % k1:this_meshblock % ptr % k2)
    real, intent(out) :: intfx, intfy, intfz
    real :: c000, c100, c001, c101, c010, c110, c011, c111, &
            c00, c01, c10, c11, c0, c1
    integer, optional, intent(in) :: comp
    integer :: comp_
    if (.not. present(comp)) then
      comp_ = 0
    else
      comp_ = comp
    end if
    ! f_x
    if ((comp_ .eq. 1) .or. (comp_ .le. 0)) then
#ifdef oneD
      c0 = 0.5 * (fx(i, j, k) + fx(i - 1, j, k))
      c1 = 0.5 * (fx(i, j, k) + fx(i + 1, j, k))
      intfx = c0 * (1 - dx) + c1 * dx
#elif defined(twoD) || defined(threeD)
      c000 = 0.5 * (fx(i, j, k) + fx(i - 1, j, k))
      c100 = 0.5 * (fx(i, j, k) + fx(i + 1, j, k))
      c010 = 0.5 * (fx(i, j + 1, k) + fx(i - 1, j + 1, k))
      c110 = 0.5 * (fx(i, j + 1, k) + fx(i + 1, j + 1, k))
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c0 = c00 * (1 - dy) + c10 * dy
#ifdef twoD
      intfx = c0
#else
      c001 = 0.5 * (fx(i, j, k + 1) + fx(i - 1, j, k + 1))
      c101 = 0.5 * (fx(i, j, k + 1) + fx(i + 1, j, k + 1))
      c011 = 0.5 * (fx(i, j + 1, k + 1) + fx(i - 1, j + 1, k + 1))
      c111 = 0.5 * (fx(i, j + 1, k + 1) + fx(i + 1, j + 1, k + 1))
      c01 = c001 * (1 - dx) + c101 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c1 = c01 * (1 - dy) + c11 * dy
      intfx = c0 * (1 - dz) + c1 * dz
#endif
#endif
    end if

    ! f_y
    if ((comp_ .eq. 2) .or. (comp_ .le. 0)) then
#ifdef oneD
      c0 = 0.5 * (fy(i, j, k) + fy(i, j, k))
      c1 = 0.5 * (fy(i + 1, j, k) + fy(i + 1, j, k))
      intfy = c0 * (1 - dx) + c1 * dx
#elif defined(twoD) || defined(threeD)
      c000 = 0.5 * (fy(i, j, k) + fy(i, j - 1, k))
      c100 = 0.5 * (fy(i + 1, j, k) + fy(i + 1, j - 1, k))
      c010 = 0.5 * (fy(i, j, k) + fy(i, j + 1, k))
      c110 = 0.5 * (fy(i + 1, j, k) + fy(i + 1, j + 1, k))
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c0 = c00 * (1 - dy) + c10 * dy
#ifdef twoD
      intfy = c0
#else
      c001 = 0.5 * (fy(i, j, k + 1) + fy(i, j - 1, k + 1))
      c101 = 0.5 * (fy(i + 1, j, k + 1) + fy(i + 1, j - 1, k + 1))
      c011 = 0.5 * (fy(i, j, k + 1) + fy(i, j + 1, k + 1))
      c111 = 0.5 * (fy(i + 1, j, k + 1) + fy(i + 1, j + 1, k + 1))
      c01 = c001 * (1 - dx) + c101 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c1 = c01 * (1 - dy) + c11 * dy
      intfy = c0 * (1 - dz) + c1 * dz
#endif
#endif
    end if

    ! f_z
    if ((comp_ .eq. 3) .or. (comp_ .le. 0)) then
#ifdef oneD
      c0 = fz(i, j, k)
      c1 = fz(i + 1, j, k)
      intfz = c0 * (1 - dx) + c1 * dx
#elif defined(twoD)
      c000 = fz(i, j, k)
      c100 = fz(i + 1, j, k)
      c010 = fz(i, j + 1, k)
      c110 = fz(i + 1, j + 1, k)
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      intfz = c00 * (1 - dy) + c10 * dy
#elif defined(threeD)
      c000 = 0.5 * (fz(i, j, k) + fz(i, j, k - 1))
      c100 = 0.5 * (fz(i + 1, j, k) + fz(i + 1, j, k - 1))
      c010 = 0.5 * (fz(i, j + 1, k) + fz(i, j + 1, k - 1))
      c110 = 0.5 * (fz(i + 1, j + 1, k) + fz(i + 1, j + 1, k - 1))
      c001 = 0.5 * (fz(i, j, k) + fz(i, j, k + 1))
      c101 = 0.5 * (fz(i + 1, j, k) + fz(i + 1, j, k + 1))
      c011 = 0.5 * (fz(i, j + 1, k) + fz(i, j + 1, k + 1))
      c111 = 0.5 * (fz(i + 1, j + 1, k) + fz(i + 1, j + 1, k + 1))
      c00 = c000 * (1 - dx) + c100 * dx
      c01 = c001 * (1 - dx) + c101 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c0 = c00 * (1 - dy) + c10 * dy
      c1 = c01 * (1 - dy) + c11 * dy
      intfz = c0 * (1 - dz) + c1 * dz
#endif
    end if
  end subroutine interpFromEdges

  subroutine interpFromFaces(dx, dy, dz, i, j, k, &
                             fx, fy, fz, &
                             intfx, intfy, intfz, comp)
    implicit none
    integer(kind=2), intent(in) :: i, j, k
    real, intent(in) :: dx, dy, dz
    real, intent(in) :: fx(this_meshblock % ptr % i1:this_meshblock % ptr % i2, &
                           this_meshblock % ptr % j1:this_meshblock % ptr % j2, &
                           this_meshblock % ptr % k1:this_meshblock % ptr % k2)
    real, intent(in) :: fy(this_meshblock % ptr % i1:this_meshblock % ptr % i2, &
                           this_meshblock % ptr % j1:this_meshblock % ptr % j2, &
                           this_meshblock % ptr % k1:this_meshblock % ptr % k2)
    real, intent(in) :: fz(this_meshblock % ptr % i1:this_meshblock % ptr % i2, &
                           this_meshblock % ptr % j1:this_meshblock % ptr % j2, &
                           this_meshblock % ptr % k1:this_meshblock % ptr % k2)
    real, intent(out) :: intfx, intfy, intfz
    real :: c000, c100, c001, c101, c010, c110, c011, c111, &
            c00, c01, c10, c11, c0, c1
    integer, optional, intent(in) :: comp
    integer :: comp_
    if (.not. present(comp)) then
      comp_ = 0
    else
      comp_ = comp
    end if
    ! f_x
    if ((comp_ .eq. 1) .or. (comp_ .le. 0)) then
#ifdef oneD
      c0 = 0.5 * (fx(i, j, k) + fx(i, j, k))
      c1 = 0.5 * (fx(i + 1, j, k) + fx(i + 1, j, k))
      intfx = c0 * (1 - dx) + c1 * dx
#elif defined(twoD)
      c000 = 0.5 * (fx(i, j, k) + fx(i, j - 1, k))
      c100 = 0.5 * (fx(i + 1, j, k) + fx(i + 1, j - 1, k))
      c010 = 0.5 * (fx(i, j, k) + fx(i, j + 1, k))
      c110 = 0.5 * (fx(i + 1, j, k) + fx(i + 1, j + 1, k))
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      intfx = c00 * (1 - dy) + c10 * dy
#elif defined(threeD)
      c000 = 0.25 * (fx(i, j, k) + fx(i, j - 1, k) + &
                     fx(i, j, k - 1) + fx(i, j - 1, k - 1))
      c100 = 0.25 * (fx(i + 1, j, k) + fx(i + 1, j - 1, k) + &
                     fx(i + 1, j, k - 1) + fx(i + 1, j - 1, k - 1))
      c001 = 0.25 * (fx(i, j, k) + fx(i, j, k + 1) + &
                     fx(i, j - 1, k) + fx(i, j - 1, k + 1))
      c101 = 0.25 * (fx(i + 1, j, k) + fx(i + 1, j, k + 1) + &
                     fx(i + 1, j - 1, k) + fx(i + 1, j - 1, k + 1))
      c010 = 0.25 * (fx(i, j, k) + fx(i, j + 1, k) + &
                     fx(i, j, k - 1) + fx(i, j + 1, k - 1))
      c110 = 0.25 * (fx(i + 1, j, k) + fx(i + 1, j, k - 1) + &
                     fx(i + 1, j + 1, k - 1) + fx(i + 1, j + 1, k))
      c011 = 0.25 * (fx(i, j, k) + fx(i, j + 1, k) + &
                     fx(i, j + 1, k + 1) + fx(i, j, k + 1))
      c111 = 0.25 * (fx(i + 1, j, k) + fx(i + 1, j + 1, k) + &
                     fx(i + 1, j + 1, k + 1) + fx(i + 1, j, k + 1))
      c00 = c000 * (1 - dx) + c100 * dx
      c01 = c001 * (1 - dx) + c101 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c0 = c00 * (1 - dy) + c10 * dy
      c1 = c01 * (1 - dy) + c11 * dy
      intfx = c0 * (1 - dz) + c1 * dz
#endif
    end if

    ! b_y
    if ((comp_ .eq. 2) .or. (comp_ .le. 0)) then
#ifdef oneD
      c0 = 0.5 * (fy(i - 1, j, k) + fy(i, j, k))
      c1 = 0.5 * (fy(i, j, k) + fy(i + 1, j, k))
      intfy = c0 * (1 - dx) + c1 * dx
#elif defined(twoD)
      c000 = 0.5 * (fy(i - 1, j, k) + fy(i, j, k))
      c100 = 0.5 * (fy(i, j, k) + fy(i + 1, j, k))
      c010 = 0.5 * (fy(i - 1, j + 1, k) + fy(i, j + 1, k))
      c110 = 0.5 * (fy(i, j + 1, k) + fy(i + 1, j + 1, k))
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      intfy = c00 * (1 - dy) + c10 * dy
#elif defined(threeD)
      c000 = 0.25 * (fy(i - 1, j, k - 1) + fy(i - 1, j, k) + &
                     fy(i, j, k - 1) + fy(i, j, k))
      c100 = 0.25 * (fy(i, j, k - 1) + fy(i, j, k) + &
                     fy(i + 1, j, k - 1) + fy(i + 1, j, k))
      c001 = 0.25 * (fy(i - 1, j, k) + fy(i - 1, j, k + 1) + &
                     fy(i, j, k) + fy(i, j, k + 1))
      c101 = 0.25 * (fy(i, j, k) + fy(i, j, k + 1) + &
                     fy(i + 1, j, k) + fy(i + 1, j, k + 1))
      c010 = 0.25 * (fy(i - 1, j + 1, k - 1) + fy(i - 1, j + 1, k) + &
                     fy(i, j + 1, k - 1) + fy(i, j + 1, k))
      c110 = 0.25 * (fy(i, j + 1, k - 1) + fy(i, j + 1, k) + &
                     fy(i + 1, j + 1, k - 1) + fy(i + 1, j + 1, k))
      c011 = 0.25 * (fy(i - 1, j + 1, k) + fy(i - 1, j + 1, k + 1) + &
                     fy(i, j + 1, k) + fy(i, j + 1, k + 1))
      c111 = 0.25 * (fy(i, j + 1, k) + fy(i, j + 1, k + 1) + &
                     fy(i + 1, j + 1, k) + fy(i + 1, j + 1, k + 1))
      c00 = c000 * (1 - dx) + c100 * dx
      c01 = c001 * (1 - dx) + c101 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c0 = c00 * (1 - dy) + c10 * dy
      c1 = c01 * (1 - dy) + c11 * dy
      intfy = c0 * (1 - dz) + c1 * dz
#endif
    end if

    ! b_z
    if ((comp_ .eq. 3) .or. (comp_ .le. 0)) then
#ifdef oneD
      c0 = 0.5 * (fz(i - 1, j, k) + fz(i, j, k))
      c1 = 0.5 * (fz(i, j, k) + fz(i + 1, j, k))
      intfz = c0 * (1 - dx) + c1 * dx
#elif defined(twoD)
      c000 = 0.25 * (fz(i - 1, j - 1, k) + fz(i - 1, j, k) + &
                     fz(i, j - 1, k) + fz(i, j, k))
      c100 = 0.25 * (fz(i, j - 1, k) + fz(i, j, k) + &
                     fz(i + 1, j - 1, k) + fz(i + 1, j, k))
      c010 = 0.25 * (fz(i - 1, j, k) + fz(i - 1, j + 1, k) + &
                     fz(i, j, k) + fz(i, j + 1, k))
      c110 = 0.25 * (fz(i, j, k) + fz(i, j + 1, k) + &
                     fz(i + 1, j, k) + fz(i + 1, j + 1, k))
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      intfz = c00 * (1 - dy) + c10 * dy
#elif defined(threeD)
      c000 = 0.25 * (fz(i - 1, j - 1, k) + fz(i - 1, j, k) + &
                     fz(i, j - 1, k) + fz(i, j, k))
      c100 = 0.25 * (fz(i, j - 1, k) + fz(i, j, k) + &
                     fz(i + 1, j - 1, k) + fz(i + 1, j, k))
      c001 = 0.25 * (fz(i - 1, j - 1, k + 1) + fz(i - 1, j, k + 1) + &
                     fz(i, j - 1, k + 1) + fz(i, j, k + 1))
      c101 = 0.25 * (fz(i, j - 1, k + 1) + fz(i, j, k + 1) + &
                     fz(i + 1, j - 1, k + 1) + fz(i + 1, j, k + 1))
      c010 = 0.25 * (fz(i - 1, j, k) + fz(i - 1, j + 1, k) + &
                     fz(i, j, k) + fz(i, j + 1, k))
      c110 = 0.25 * (fz(i, j, k) + fz(i, j + 1, k) + &
                     fz(i + 1, j, k) + fz(i + 1, j + 1, k))
      c011 = 0.25 * (fz(i - 1, j, k + 1) + fz(i - 1, j + 1, k + 1) + &
                     fz(i, j, k + 1) + fz(i, j + 1, k + 1))
      c111 = 0.25 * (fz(i, j, k + 1) + fz(i, j + 1, k + 1) + &
                     fz(i + 1, j, k + 1) + fz(i + 1, j + 1, k + 1))
      c00 = c000 * (1 - dx) + c100 * dx
      c01 = c001 * (1 - dx) + c101 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c0 = c00 * (1 - dy) + c10 * dy
      c1 = c01 * (1 - dy) + c11 * dy
      intfz = c0 * (1 - dz) + c1 * dz
#endif
    end if
  end subroutine interpFromFaces

  subroutine depositCurrentsFromSingleParticle(s, tile, p, x1, y1, z1, &
                                               x2, y2, z2, multiplier)
    implicit none
    type(particle_tile), intent(in) :: tile
    integer, intent(in) :: s, p
    real, intent(in) :: x1, y1, z1, x2, y2, z2
    real, optional, intent(in) :: multiplier
    real :: xr, yr, zr
    integer(kind=2) :: i1, i2, j1, j2, k1, k2
    integer(kind=2) :: i1p1, i2p1, j1p1, j2p1, k1p1, k2p1
    real :: Wx1, Wy1, Wz1, Wx2, Wy2, Wz2
    real :: onemWx1, onemWy1, onemWz1, onemWx2, onemWy2, onemWz2
    real :: Fx1, Fy1, Fz1, Fx2, Fy2, Fz2
    real :: weighted_charge

    if (present(multiplier)) then
      weighted_charge = multiplier * tile % weight(p) * species(s) % ch_sp * unit_ch / B_norm
    else
      weighted_charge = tile % weight(p) * species(s) % ch_sp * unit_ch / B_norm
    end if

#ifdef oneD
    i1 = INT(FLOOR(x1), 2); i2 = INT(FLOOR(x2), 2)
    j1 = 0; j2 = 0
    k1 = 0; k2 = 0
    i1p1 = i1 + 1_2; i2p1 = i2 + 1_2
#elif defined(twoD)
    i1 = INT(FLOOR(x1), 2); i2 = INT(FLOOR(x2), 2)
    j1 = INT(FLOOR(y1), 2); j2 = INT(FLOOR(y2), 2)
    k1 = 0; k2 = 0
    i1p1 = i1 + 1_2; i2p1 = i2 + 1_2
    j1p1 = j1 + 1_2; j2p1 = j2 + 1_2
#elif defined(threeD)
    i1 = INT(FLOOR(x1), 2); i2 = INT(FLOOR(x2), 2)
    j1 = INT(FLOOR(y1), 2); j2 = INT(FLOOR(y2), 2)
    k1 = INT(FLOOR(z1), 2); k2 = INT(FLOOR(z2), 2)
    i1p1 = i1 + 1_2; i2p1 = i2 + 1_2
    j1p1 = j1 + 1_2; j2p1 = j2 + 1_2
    k1p1 = k1 + 1_2; k2p1 = k2 + 1_2
#endif

    ! this "function" takes
    ! ... the start and end coordinates: `x1`, `x2`, `y1`, `y2`, `z1`, `z2` ...
    ! ... the start and end cells: `i1`, `i2`, `j1`, `j2`, `k1`, `k2` ...
    ! ... the start and end cells + 1: `i1p1`, `i2p1` etc ...
    ! ... the weighted_chargeed charge: `weighted_charge = weight * charge_sp * unit_charge / Bnorm`
    ! ... and deposits proper currents to corresponding components
#include "zigzag_deposit.F08"
  end subroutine depositCurrentsFromSingleParticle

  subroutine checkEverything()
    implicit none
    ! this is an ultimate function that checks that all the implicit assertions are satisfied
    integer, dimension(3) :: field_shape
    integer :: s, ti, tj, tk, p, i, j, k, dummy1, dummy2, ierr, shape_x
    integer(kind=2), pointer, contiguous :: pt_xi(:), pt_yi(:), pt_zi(:)
    real, pointer, contiguous :: pt_dx(:), pt_dy(:), pt_dz(:)
    real, pointer, contiguous :: pt_ux(:), pt_uy(:), pt_uz(:)
    integer, pointer, contiguous :: pt_proc(:)

    ! 1. check that the domain size is larger than the number of ghost zones
#ifdef oneD
    if (this_meshblock % ptr % sx .lt. NGHOST) then
      call throwError('ERROR: ghost zones overflow the domain size in '//trim(STR(mpi_rank)))
    end if
#elif defined(twoD)
    if ((this_meshblock % ptr % sx .lt. NGHOST) .or. &
        (this_meshblock % ptr % sy .lt. NGHOST)) then
      call throwError('ERROR: ghost zones overflow the domain size in '//trim(STR(mpi_rank)))
    end if
#elif defined(threeD)
    if ((this_meshblock % ptr % sx .lt. NGHOST) .or. &
        (this_meshblock % ptr % sy .lt. NGHOST) .or. &
        (this_meshblock % ptr % sz .lt. NGHOST)) then
      call throwError('ERROR: ghost zones overflow the domain size in '//trim(STR(mpi_rank)))
    end if
#endif

    ! 2. check field dimensions
    shape_x = this_meshblock % ptr % sx + 2 * NGHOST
    shape_x = ((shape_x + VEC_LEN - 1) / VEC_LEN) * VEC_LEN
#ifdef oneD
    field_shape = (/shape_x, 1, 1/)
#elif defined(twoD)
    field_shape = (/shape_x, this_meshblock % ptr % sy + 2 * NGHOST, 1/)
#elif defined(threeD)
    field_shape = (/shape_x, this_meshblock % ptr % sy + 2 * NGHOST, this_meshblock % ptr % sz + 2 * NGHOST/)
#endif

    if (.not. arraysAreEqual(shape(ex), field_shape)) then
      print *, shape(ex), field_shape
      call throwError('ERROR: incorrect shape of `ex`.')
    end if
    if (.not. arraysAreEqual(shape(ey), field_shape)) then
      call throwError('ERROR: incorrect shape of `ey`.')
    end if
    if (.not. arraysAreEqual(shape(ez), field_shape)) then
      call throwError('ERROR: incorrect shape of `ez`.')
    end if
    if (.not. arraysAreEqual(shape(bx), field_shape)) then
      call throwError('ERROR: incorrect shape of `bx`.')
    end if
    if (.not. arraysAreEqual(shape(by), field_shape)) then
      call throwError('ERROR: incorrect shape of `by`.')
    end if
    if (.not. arraysAreEqual(shape(bz), field_shape)) then
      call throwError('ERROR: incorrect shape of `bz`.')
    end if
    if (.not. arraysAreEqual(shape(jx), field_shape)) then
      call throwError('ERROR: incorrect shape of `jx`.')
    end if
    if (.not. arraysAreEqual(shape(jy), field_shape)) then
      call throwError('ERROR: incorrect shape of `jy`.')
    end if
    if (.not. arraysAreEqual(shape(jz), field_shape)) then
      call throwError('ERROR: incorrect shape of `jz`.')
    end if
    if (.not. arraysAreEqual(shape(jx_buff), field_shape)) then
      call throwError('ERROR: incorrect shape of `jx_buff`.')
    end if
    if (.not. arraysAreEqual(shape(jy_buff), field_shape)) then
      call throwError('ERROR: incorrect shape of `jy_buff`.')
    end if
    if (.not. arraysAreEqual(shape(jz_buff), field_shape)) then
      call throwError('ERROR: incorrect shape of `jz_buff`.')
    end if
    if (.not. arraysAreEqual(shape(lg_arr), field_shape)) then
      call throwError('ERROR: incorrect shape of `lg_arr`.')
    end if
    field_shape = (/this_meshblock % ptr % sx, this_meshblock % ptr % sy, this_meshblock % ptr % sz/)
    if (.not. arraysAreEqual(shape(sm_arr), field_shape)) then
      call throwError('ERROR: incorrect shape of `sm_arr`.')
    end if

    ! 3. check tiles & particles
    do s = 1, nspec
      do ti = 1, species(s) % tile_nx
        dummy1 = species(s) % prtl_tile(ti, 1, 1) % x1
        dummy2 = species(s) % prtl_tile(ti, 1, 1) % x2
        do tj = 2, species(s) % tile_ny
          do tk = 2, species(s) % tile_nz
            if ((species(s) % prtl_tile(ti, tj, tk) % x1 .ne. dummy1) .or. &
                (species(s) % prtl_tile(ti, tj, tk) % x2 .ne. dummy2)) then
              call throwError('ERROR: tiles are misaligned in x.')
            end if
          end do
        end do
      end do

      do tj = 1, species(s) % tile_ny
        dummy1 = species(s) % prtl_tile(1, tj, 1) % y1
        dummy2 = species(s) % prtl_tile(1, tj, 1) % y2
        do ti = 2, species(s) % tile_nx
          do tk = 2, species(s) % tile_nz
            if ((species(s) % prtl_tile(ti, tj, tk) % y1 .ne. dummy1) .or. &
                (species(s) % prtl_tile(ti, tj, tk) % y2 .ne. dummy2)) then
              call throwError('ERROR: tiles are misaligned in y.')
            end if
          end do
        end do
      end do

      do tk = 1, species(s) % tile_nz
        dummy1 = species(s) % prtl_tile(1, 1, tk) % z1
        dummy2 = species(s) % prtl_tile(1, 1, tk) % z2
        do ti = 2, species(s) % tile_nx
          do tj = 2, species(s) % tile_ny
            if ((species(s) % prtl_tile(ti, tj, tk) % z1 .ne. dummy1) .or. &
                (species(s) % prtl_tile(ti, tj, tk) % z2 .ne. dummy2)) then
              call throwError('ERROR: tiles are misaligned in z.')
            end if
          end do
        end do
      end do

      do tk = 1, species(s) % tile_nz
        do tj = 1, species(s) % tile_ny
          do ti = 1, species(s) % tile_nx

            if ((species(s) % prtl_tile(ti, tj, tk) % x1 .lt. 0) .or. &
                (species(s) % prtl_tile(ti, tj, tk) % x2 .gt. this_meshblock % ptr % sx) .or. &
                (species(s) % prtl_tile(ti, tj, tk) % y1 .lt. 0) .or. &
                (species(s) % prtl_tile(ti, tj, tk) % y2 .gt. this_meshblock % ptr % sy) .or. &
                (species(s) % prtl_tile(ti, tj, tk) % z1 .lt. 0) .or. &
                (species(s) % prtl_tile(ti, tj, tk) % z2 .gt. this_meshblock % ptr % sz)) then
              call throwError('ERROR: incorrect tile sizes.')
            end if

            pt_xi => species(s) % prtl_tile(ti, tj, tk) % xi
            pt_yi => species(s) % prtl_tile(ti, tj, tk) % yi
            pt_zi => species(s) % prtl_tile(ti, tj, tk) % zi
            pt_dx => species(s) % prtl_tile(ti, tj, tk) % dx
            pt_dy => species(s) % prtl_tile(ti, tj, tk) % dy
            pt_dz => species(s) % prtl_tile(ti, tj, tk) % dz
            pt_ux => species(s) % prtl_tile(ti, tj, tk) % u
            pt_uy => species(s) % prtl_tile(ti, tj, tk) % v
            pt_uz => species(s) % prtl_tile(ti, tj, tk) % w
            pt_proc => species(s) % prtl_tile(ti, tj, tk) % proc
            do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
              i = pt_xi(p); j = pt_yi(p); k = pt_zi(p)
#if defined(oneD) || defined (twoD) || defined (threeD)
              if ((i .lt. species(s) % prtl_tile(ti, tj, tk) % x1) .or. &
                  (i .ge. species(s) % prtl_tile(ti, tj, tk) % x2)) then
                call throwError('ERROR: particle in wrong tile in i.')
              end if
#endif
#if defined (twoD) || defined (threeD)
              if ((j .lt. species(s) % prtl_tile(ti, tj, tk) % y1) .or. &
                  (j .ge. species(s) % prtl_tile(ti, tj, tk) % y2)) then
                call throwError('ERROR: particle in wrong tile in j.')
              end if
#endif
#if defined(threeD)
              if ((k .lt. species(s) % prtl_tile(ti, tj, tk) % z1) .or. &
                  (k .ge. species(s) % prtl_tile(ti, tj, tk) % z2)) then
                call throwError('ERROR: particle in wrong tile in k.')
              end if
#endif

#if defined (oneD)
              pt_yi(p) = 0; pt_zi(p) = 0; 
#elif defined (twoD)
              pt_zi(p) = 0; 
#endif

#if defined(oneD) || defined (twoD) || defined (threeD)
              if ((pt_dx(p) .lt. 0) .or. (pt_dx(p) .gt. 1)) then
                call throwError('ERROR: invalid particle coordinate in x.')
              end if
#endif
#if defined (twoD) || defined (threeD)
              if ((pt_dy(p) .lt. 0) .or. (pt_dy(p) .gt. 1)) then
                call throwError('ERROR: invalid particle coordinate in y.')
              end if
#endif
#if defined(threeD)
              if ((pt_dz(p) .lt. 0) .or. (pt_dz(p) .gt. 1)) then
                call throwError('ERROR: invalid particle coordinate in z.')
              end if
#endif

              ! safeguard for crazy velocities
              if (abs(pt_ux(p)) .gt. 1e10) then
                call throwError('ERROR: particles have crazy velocities in x.')
              end if
              if (abs(pt_uy(p)) .gt. 1e10) then
                call throwError('ERROR: particles have crazy velocities in y.')
              end if
              if (abs(pt_uz(p)) .gt. 1e10) then
                call throwError('ERROR: particles have crazy velocities in z.')
              end if

#if defined (oneD)
              pt_dy(p) = 0.5; pt_dz(p) = 0.5; 
#elif defined (twoD)
              pt_dz(p) = 0.5; 
#endif

              if (abs(pt_proc(p)) .gt. 10 * mpi_size) then
                call throwError('ERROR: particle proc wrong.')
              end if
            end do
            pt_ux => null(); pt_uy => null(); pt_uz => null()
            pt_xi => null(); pt_yi => null(); pt_zi => null()
            pt_dx => null(); pt_dy => null(); pt_dz => null()
            pt_proc => null()
          end do
        end do
      end do
    end do

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    call printDiag("checkEverything()", 1)
  end subroutine checkEverything
end module m_helpers
