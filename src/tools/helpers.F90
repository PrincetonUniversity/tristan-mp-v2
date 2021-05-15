#include "../defs.F90"

module m_helpers
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  implicit none
contains
  logical function numbersAreClose(number1, number2)
    implicit none
    real, intent(in)  :: number1, number2
    real              :: abs1, abs2
    real              :: diff

    abs1 = abs(number1); abs2 = abs(number2)
    diff = abs(number1 - number2)

    if (number1 .eq. number2) then
      numbersAreClose = .true.
    else if ((number1 .eq. 0.0) .or. (number2 .eq. 0.0) .or.&
           & (abs1 + abs2 .lt. TINYREAL)) then
      numbersAreClose = (diff .lt. TINYREAL)
    else
      numbersAreClose = (diff / (abs1 + abs2) .lt. 0.5 * TINYREAL)
    end if
  end function numbersAreClose

  subroutine checkNpart(msg)
    implicit none
    integer             :: ti, tj, tk, s, nprt
    character(len=*), intent(in) :: msg
    integer :: glob_nprt
    do s = 1, nspec
      nprt = 0
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            nprt = nprt + species(s)%prtl_tile(ti, tj, tk)%npart_sp
          end do
        end do
      end do
      call MPI_REDUCE(nprt, glob_nprt, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD)
      if (mpi_rank .eq. 0) then
        print *, "npart of", s, msg, ":", glob_nprt
      end if
    end do
  end subroutine checkNpart

  subroutine globalToLocalCoords(x_glob, y_glob, z_glob,&
                               & x_loc, y_loc, z_loc, adjustQ, containedQ)
    implicit none
    real, intent(in)              :: x_glob, y_glob, z_glob
    real, intent(out)             :: x_loc, y_loc, z_loc
    logical, optional, intent(in) :: adjustQ
    logical, optional, intent(out):: containedQ
    logical                       :: adjustQ_
    if (present(adjustQ)) then
      adjustQ_ = adjustQ
    else
      adjustQ_ = .false.
    end if

    if (adjustQ_) then
      x_loc = x_glob; y_loc = y_glob; z_loc = z_glob
      #if defined(oneD) || defined (twoD) || defined (threeD)
        x_loc = MAX(0.0, MIN(x_glob - REAL(this_meshblock%ptr%x0), REAL(this_meshblock%ptr%sx)))
      #endif

      #if defined (twoD) || defined (threeD)
        y_loc = MAX(0.0, MIN(y_glob - REAL(this_meshblock%ptr%y0), REAL(this_meshblock%ptr%sy)))
      #endif

      #if defined(threeD)
        z_loc = MAX(0.0, MIN(z_glob - REAL(this_meshblock%ptr%z0), REAL(this_meshblock%ptr%sz)))
      #endif
    else
      x_loc = x_glob; y_loc = y_glob; z_loc = z_glob
      #if defined(oneD) || defined (twoD) || defined (threeD)
        x_loc = x_glob - REAL(this_meshblock%ptr%x0)
      #endif

      #if defined (twoD) || defined (threeD)
        y_loc = y_glob - REAL(this_meshblock%ptr%y0)
      #endif

      #if defined(threeD)
        z_loc = z_glob - REAL(this_meshblock%ptr%z0)
      #endif
      if (present(containedQ)) then
        #ifdef oneD
          containedQ = ((x_glob .ge. REAL(this_meshblock%ptr%x0)) .and.&
                      & (x_glob .lt. REAL(this_meshblock%ptr%x0 + this_meshblock%ptr%sx)))
        #elif twoD
          containedQ = ((x_glob .ge. REAL(this_meshblock%ptr%x0)) .and.&
                      & (x_glob .lt. REAL(this_meshblock%ptr%x0 + this_meshblock%ptr%sx)) .and.&
                      & (y_glob .ge. REAL(this_meshblock%ptr%y0)) .and.&
                      & (y_glob .lt. REAL(this_meshblock%ptr%y0 + this_meshblock%ptr%sy)))
        #elif threeD
          containedQ = ((x_glob .ge. REAL(this_meshblock%ptr%x0)) .and.&
                      & (x_glob .lt. REAL(this_meshblock%ptr%x0 + this_meshblock%ptr%sx)) .and.&
                      & (y_glob .ge. REAL(this_meshblock%ptr%y0)) .and.&
                      & (y_glob .lt. REAL(this_meshblock%ptr%y0 + this_meshblock%ptr%sy)) .and.&
                      & (z_glob .ge. REAL(this_meshblock%ptr%z0)) .and.&
                      & (z_glob .lt. REAL(this_meshblock%ptr%z0 + this_meshblock%ptr%sz)))
        #endif
      end if
    end if
  end subroutine globalToLocalCoords

  subroutine localToCellBasedCoords(x_loc, y_loc, z_loc,&
                                  & xi, yi, zi, dx, dy, dz)
    implicit none
    real, intent(in)              :: x_loc, y_loc, z_loc
    real, intent(out)             :: dx, dy, dz
    integer(kind=2), intent(out)  :: xi, yi, zi

    xi = FLOOR(x_loc); dx = x_loc - FLOOR(x_loc)
    yi = FLOOR(y_loc); dy = y_loc - FLOOR(y_loc)
    zi = FLOOR(z_loc); dz = z_loc - FLOOR(z_loc)
  end subroutine localToCellBasedCoords

  subroutine generateCoordInRegion(xmin, xmax, ymin, ymax, zmin, zmax,&
                                 & x_, y_, z_, xi_, yi_, zi_, dx_, dy_, dz_)
    implicit none
    real                          :: rnd
    real, intent(in)              :: xmin, xmax, ymin, ymax, zmin, zmax
    real, intent(out)             :: x_, y_, z_, dx_, dy_, dz_
    integer(kind=2), intent(out)  :: xi_, yi_, zi_

    x_ = 0.5; xi_ = 0; dx_ = 0.5
    y_ = 0.5; yi_ = 0; dy_ = 0.5
    z_ = 0.5; zi_ = 0; dz_ = 0.5
    #if defined(oneD) || defined (twoD) || defined (threeD)
      rnd = random(dseed)
      x_ = xmin + rnd * (xmax - xmin)
      xi_ = FLOOR(x_); dx_ = x_ - FLOOR(x_)
      if (xi_ .eq. this_meshblock%ptr%sx) then
        xi_ = xi_ - 1; dx_ = dx_ + 1.0
      end if
    #endif

    #if defined (twoD) || defined (threeD)
      rnd = random(dseed)
      y_ = ymin + rnd * (ymax - ymin)
      yi_ = FLOOR(y_); dy_ = y_ - FLOOR(y_)
      if (yi_ .eq. this_meshblock%ptr%sy) then
        yi_ = yi_ - 1; dy_ = dy_ + 1.0
      end if
    #endif

    #if defined(threeD)
      rnd = random(dseed)
      z_ = zmin + rnd * (zmax - zmin)
      zi_ = FLOOR(z_); dz_ = z_ - FLOOR(z_)
      if (zi_ .eq. this_meshblock%ptr%sz) then
        zi_ = zi_ - 1; dz_ = dz_ + 1.0
      end if
    #endif
  end subroutine generateCoordInRegion

  function rnkToInd(rnk)
    implicit none
    integer, intent(in)   :: rnk
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
    integer, intent(in)                 :: ind(3)
    integer                             :: ind_t(3), indToRnk
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
    if ((ind_t(1) .lt. 0) .or. (ind_t(1) .ge. sizex) .or.&
      & (ind_t(2) .lt. 0) .or. (ind_t(2) .ge. sizey) .or.&
      & (ind_t(3) .lt. 0) .or. (ind_t(3) .ge. sizez)) then
      indToRnk = -1
    else
      indToRnk = ind_t(3) * sizex * sizey + ind_t(2) * sizex + ind_t(1)
    end if
  end function indToRnk

  subroutine reassignNeighborsForAll()
    implicit none
    integer   :: rnk
    integer   :: ind1, ind2, ind3
    do rnk = 0, mpi_size - 1
      do ind1 = -1, 1
        do ind2 = -1, 1
          do ind3 = -1, 1
            call assignNeighbor(rnk, (/ ind1, ind2, ind3/))
          end do
        end do
      end do
    end do
    call computeNumberOfNeighbors()
  end subroutine reassignNeighborsForAll

  subroutine assignNeighbor(rnk, inds1)
    implicit none
    integer, intent(in)       :: rnk, inds1(3)
    integer                   :: rnk2, inds0(3)
    inds0 = rnkToInd(rnk)
    rnk2 = indToRnk([inds0(1) + inds1(1), inds0(2) + inds1(2), inds0(3) + inds1(3)])
    if (rnk2 .eq. -1) then
      meshblocks(rnk + 1)%neighbor(inds1(1), inds1(2), inds1(3))%ptr => null()
    else
      meshblocks(rnk + 1)%neighbor(inds1(1), inds1(2), inds1(3))%ptr => meshblocks(rnk2 + 1)
    end if
  end subroutine assignNeighbor

  subroutine computeNumberOfNeighbors()
    integer   :: ind1, ind2, ind3
    integer   :: cntr
    cntr = 0
    do ind1 = -1, 1
      do ind2 = -1, 1
        do ind3 = -1, 1
          if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
          #ifdef oneD
            if ((ind2 .ne. 0) .or. (ind3 .ne. 0)) cycle
          #elif twoD
            if (ind3 .ne. 0) cycle
          #endif
          if (.not. associated(this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr)) cycle
          cntr = cntr + 1
        end do
      end do
    end do
    sendrecv_neighbors = cntr
  end subroutine computeNumberOfNeighbors

  subroutine computeDensity(s, reset, ds, charge)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in)                   :: s
    logical, intent(in)                   :: reset
    logical, optional, intent(in)         :: charge
    integer, optional, intent(in)         :: ds
    integer                               :: p, ti, tj, tk
    integer(kind=2), pointer, contiguous  :: pt_xi(:), pt_yi(:), pt_zi(:)
    real, pointer, contiguous             :: pt_wei(:)
    logical                               :: charge_
    integer(kind=2) :: i, j, k
    integer :: i1, i2, j1, j2, k1, k2, ds_
    integer :: pow
    real    :: contrib

    if (.not. present(ds)) then
      ds_ = 2
    else
      ds_ = ds
    end if

    if (.not. present(charge)) then
      charge_ = .false.
    else
      charge_ = charge
    end if

    #ifdef oneD
      pow = 1
    #elif twoD
      pow = 2
    #elif threeD
      pow = 3
    #endif

    if (species(s)%m_sp .eq. 0) then
      contrib = 1.0 / (2.0 * REAL(ds_) + 1.0)**pow
    else
      if (charge_) then
        contrib = species(s)%ch_sp / (2.0 * REAL(ds_) + 1.0)**pow
      else
        contrib = species(s)%m_sp / (2.0 * REAL(ds_) + 1.0)**pow
      end if
    end if

    if (reset) then
      lg_arr(:,:,:) = 0
    end if

    do ti = 1, species(s)%tile_nx
      do tj = 1, species(s)%tile_ny
        do tk = 1, species(s)%tile_nz
          pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
          pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
          pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi
          pt_wei => species(s)%prtl_tile(ti, tj, tk)%weight
          do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
            i = pt_xi(p); j = pt_yi(p); k = pt_zi(p)

            i1 = 0; i2 = 0
            j1 = 0; j2 = 0
            k1 = 0; k2 = 0
            #if defined(oneD) || defined (twoD) || defined (threeD)
              i1 = max(i - ds_, -NGHOST)
              i2 = min(i + ds_, this_meshblock%ptr%sx + NGHOST - 1)
            #endif

            #if defined (twoD) || defined (threeD)
              j1 = max(j - ds_, -NGHOST)
              j2 = min(j + ds_, this_meshblock%ptr%sy + NGHOST - 1)
            #endif

            #if defined (threeD)
              k1 = max(k - ds_, -NGHOST)
              k2 = min(k + ds_, this_meshblock%ptr%sz + NGHOST - 1)
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

  subroutine computeMomentum(s, component, reset, ds)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in)                   :: s, component
    logical, intent(in)                   :: reset
    integer, optional, intent(in)         :: ds
    integer                               :: p, ti, tj, tk
    integer(kind=2), pointer, contiguous  :: pt_xi(:), pt_yi(:), pt_zi(:)
    real, pointer, contiguous             :: pt_u(:), pt_v(:), pt_w(:), pt_wei(:)
    integer(kind=2) :: i, j, k
    integer :: i1, i2, j1, j2, k1, k2, ds_
    integer :: pow
    logical :: massive
    real    :: comp
    real    :: contrib

    if (.not. present(ds)) then
      ds_ = 2
    else
      ds_ = ds
    end if

    #ifdef oneD
      pow = 1
    #elif twoD
      pow = 2
    #elif threeD
      pow = 3
    #endif

    massive = (species(s)%m_sp .ne. 0)
    if (.not. massive) then
      contrib = 1.0 / (2.0 * REAL(ds_) + 1.0)**pow
    else
      contrib = species(s)%m_sp / (2.0 * REAL(ds_) + 1.0)**pow
    end if

    if (reset) then
      lg_arr(:,:,:) = 0
    end if
    do ti = 1, species(s)%tile_nx
      do tj = 1, species(s)%tile_ny
        do tk = 1, species(s)%tile_nz
          pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
          pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
          pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi
          pt_wei => species(s)%prtl_tile(ti, tj, tk)%weight
          pt_u => species(s)%prtl_tile(ti, tj, tk)%u
          pt_v => species(s)%prtl_tile(ti, tj, tk)%v
          pt_w => species(s)%prtl_tile(ti, tj, tk)%w
          do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
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
            end if

            i1 = 0; i2 = 0
            j1 = 0; j2 = 0
            k1 = 0; k2 = 0
            #if defined(oneD) || defined (twoD) || defined (threeD)
              i1 = max(i - ds_, -NGHOST)
              i2 = min(i + ds_, this_meshblock%ptr%sx + NGHOST - 1)
            #endif

            #if defined (twoD) || defined (threeD)
              j1 = max(j - ds_, -NGHOST)
              j2 = min(j + ds_, this_meshblock%ptr%sy + NGHOST - 1)
            #endif

            #if defined (threeD)
              k1 = max(k - ds_, -NGHOST)
              k2 = min(k + ds_, this_meshblock%ptr%sz + NGHOST - 1)
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

  subroutine computeNpart(s, reset, ds)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in)                   :: s
    logical, intent(in)                   :: reset
    integer, optional, intent(in)         :: ds
    integer                               :: p, ti, tj, tk
    integer(kind=2), pointer, contiguous  :: pt_xi(:), pt_yi(:), pt_zi(:)
    real, pointer, contiguous             :: pt_u(:), pt_v(:), pt_w(:), pt_wei(:)
    integer(kind=2) :: i, j, k
    integer :: i1, i2, j1, j2, k1, k2, ds_
    integer :: pow
    real    :: contrib

    if (.not. present(ds)) then
      ds_ = 2
    else
      ds_ = ds
    end if

    #ifdef oneD
      pow = 1
    #elif twoD
      pow = 2
    #elif threeD
      pow = 3
    #endif

    contrib = 1.0 / (2.0 * REAL(ds_) + 1.0)**pow

    if (reset) then
      lg_arr(:,:,:) = 0
    end if
    do ti = 1, species(s)%tile_nx
      do tj = 1, species(s)%tile_ny
        do tk = 1, species(s)%tile_nz
          pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
          pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
          pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi
          do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
            i = pt_xi(p); j = pt_yi(p); k = pt_zi(p)
            i1 = 0; i2 = 0
            j1 = 0; j2 = 0
            k1 = 0; k2 = 0
            #if defined(oneD) || defined (twoD) || defined (threeD)
              i1 = max(i - ds_, -NGHOST)
              i2 = min(i + ds_, this_meshblock%ptr%sx + NGHOST - 1)
            #endif

            #if defined (twoD) || defined (threeD)
              j1 = max(j - ds_, -NGHOST)
              j2 = min(j + ds_, this_meshblock%ptr%sy + NGHOST - 1)
            #endif

            #if defined (threeD)
              k1 = max(k - ds_, -NGHOST)
              k2 = min(k + ds_, this_meshblock%ptr%sz + NGHOST - 1)
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
                           & fx, fy, fz, &
                           & intfx, intfy, intfz)
    implicit none
    integer(kind=2), intent(in)   :: i, j, k
    real, intent(in)              :: dx, dy, dz
    #ifdef oneD
      real, intent(inout) :: fx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0)
      real, intent(inout) :: fy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0)
      real, intent(inout) :: fz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0)
    #elif twoD
      real, intent(inout) :: fx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                               & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0)
      real, intent(inout) :: fy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0)
      real, intent(inout) :: fz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0)
    #elif threeD
      real, intent(inout) :: fx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                               & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                               & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST)
      real, intent(inout) :: fy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                              & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST)
      real, intent(inout) :: fz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                               & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                               & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST)
    #endif
    real, intent(out)             :: intfx, intfy, intfz
    real                          :: c000, c100, c001, c101, c010, c110, c011, c111,&
                                   & c00, c01, c10, c11, c0, c1
    ! f_x
    #ifdef oneD
      c0 = 0.5 * (fx(     i,      j,      k) + fx(  i - 1,      j,      k))
      c1 = 0.5 * (fx(     i,      j,      k) + fx(  i + 1,      j,      k))
      intfx = c0 * (1 - dx) + c1 * dx
    #elif defined(twoD) || defined(threeD)
      c000 = 0.5 * (fx(     i,      j,      k) + fx(  i - 1,      j,      k))
      c100 = 0.5 * (fx(     i,      j,      k) + fx(  i + 1,      j,      k))
      c010 = 0.5 * (fx(     i,  j + 1,      k) + fx(  i - 1,  j + 1,      k))
      c110 = 0.5 * (fx(     i,  j + 1,      k) + fx(  i + 1,  j + 1,      k))
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c0 = c00 * (1 - dy) + c10 * dy
      #ifdef twoD
        intfx = c0
      #else
        c001 = 0.5 * (fx(     i,      j,  k + 1) + fx(  i - 1,      j,  k + 1))
        c101 = 0.5 * (fx(     i,      j,  k + 1) + fx(  i + 1,      j,  k + 1))
        c011 = 0.5 * (fx(     i,  j + 1,  k + 1) + fx(  i - 1,  j + 1,  k + 1))
        c111 = 0.5 * (fx(     i,  j + 1,  k + 1) + fx(  i + 1,  j + 1,  k + 1))
        c01 = c001 * (1 - dx) + c101 * dx
        c11 = c011 * (1 - dx) + c111 * dx
        c1 = c01 * (1 - dy) + c11 * dy
        intfx = c0 * (1 - dz) + c1 * dz
      #endif

    #endif

    ! f_y
    #ifdef oneD
      c0 = 0.5 * (fy(     i,      j,      k) + fy(      i,      j,      k))
      c1 = 0.5 * (fy( i + 1,      j,      k) + fy(  i + 1,      j,      k))
      intfy = c0 * (1 - dx) + c1 * dx
    #elif defined(twoD) || defined(threeD)
      c000 = 0.5 * (fy(     i,      j,      k) + fy(      i,  j - 1,      k))
      c100 = 0.5 * (fy( i + 1,      j,      k) + fy(  i + 1,  j - 1,      k))
      c010 = 0.5 * (fy(     i,      j,      k) + fy(      i,  j + 1,      k))
      c110 = 0.5 * (fy( i + 1,      j,      k) + fy(  i + 1,  j + 1,      k))
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c0 = c00 * (1 - dy) + c10 * dy
      #ifdef twoD
        intfy = c0
      #else
        c001 = 0.5 * (fy(     i,      j,  k + 1) + fy(      i,  j - 1,  k + 1))
        c101 = 0.5 * (fy( i + 1,      j,  k + 1) + fy(  i + 1,  j - 1,  k + 1))
        c011 = 0.5 * (fy(     i,      j,  k + 1) + fy(      i,  j + 1,  k + 1))
        c111 = 0.5 * (fy( i + 1,      j,  k + 1) + fy(  i + 1,  j + 1,  k + 1))
        c01 = c001 * (1 - dx) + c101 * dx
        c11 = c011 * (1 - dx) + c111 * dx
        c1 = c01 * (1 - dy) + c11 * dy
        intfy = c0 * (1 - dz) + c1 * dz
      #endif

    #endif

    ! f_z
    #ifdef oneD
      c0 = fz(     i,      j,      k)
      c1 = fz( i + 1,      j,      k)
      intfz = c0 * (1 - dx) + c1 * dx
    #elif twoD
      c000 = fz(     i,      j,      k)
      c100 = fz( i + 1,      j,      k)
      c010 = fz(     i,  j + 1,      k)
      c110 = fz( i + 1,  j + 1,      k)
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      intfz = c00 * (1 - dy) + c10 * dy
    #elif threeD
      c000 = 0.5 * (fz(     i,      j,      k) + fz(      i,      j,  k - 1))
      c100 = 0.5 * (fz( i + 1,      j,      k) + fz(  i + 1,      j,  k - 1))
      c010 = 0.5 * (fz(     i,  j + 1,      k) + fz(      i,  j + 1,  k - 1))
      c110 = 0.5 * (fz( i + 1,  j + 1,      k) + fz(  i + 1,  j + 1,  k - 1))
      c001 = 0.5 * (fz(     i,      j,      k) + fz(      i,      j,  k + 1))
      c101 = 0.5 * (fz( i + 1,      j,      k) + fz(  i + 1,      j,  k + 1))
      c011 = 0.5 * (fz(     i,  j + 1,      k) + fz(      i,  j + 1,  k + 1))
      c111 = 0.5 * (fz( i + 1,  j + 1,      k) + fz(  i + 1,  j + 1,  k + 1))
      c00 = c000 * (1 - dx) + c100 * dx
      c01 = c001 * (1 - dx) + c101 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c0 = c00 * (1 - dy) + c10 * dy
      c1 = c01 * (1 - dy) + c11 * dy
      intfz = c0 * (1 - dz) + c1 * dz
    #endif
  end subroutine interpFromEdges

  subroutine interpFromFaces(dx, dy, dz, i, j, k, &
                           & fx, fy, fz, &
                           & intfx, intfy, intfz)
    implicit none
    integer(kind=2), intent(in)   :: i, j, k
    real, intent(in)              :: dx, dy, dz
    #ifdef oneD
      real, intent(inout) :: fx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0)
      real, intent(inout) :: fy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0)
      real, intent(inout) :: fz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0)
    #elif twoD
      real, intent(inout) :: fx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                               & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0)
      real, intent(inout) :: fy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0)
      real, intent(inout) :: fz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0)
    #elif threeD
      real, intent(inout) :: fx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                               & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                               & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST)
      real, intent(inout) :: fy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                              & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST)
      real, intent(inout) :: fz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                               & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                               & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST)
    #endif
    real, intent(out)             :: intfx, intfy, intfz
    real                          :: c000, c100, c001, c101, c010, c110, c011, c111,&
                                   & c00, c01, c10, c11, c0, c1
    ! f_x
    #ifdef oneD
      c0 = 0.5 * (fx(      i,      j,      k) + fx(      i,  j,      k))
      c1 = 0.5 * (fx(  i + 1,      j,      k) + fx(  i + 1,  j,      k))
      intfx = c0 * (1 - dx) + c1 * dx
    #elif twoD
      c000 = 0.5 * (fx(      i,      j,      k) + fx(      i,  j - 1,      k))
      c100 = 0.5 * (fx(  i + 1,      j,      k) + fx(  i + 1,  j - 1,      k))
      c010 = 0.5 * (fx(      i,      j,      k) + fx(      i,  j + 1,      k))
      c110 = 0.5 * (fx(  i + 1,      j,      k) + fx(  i + 1,  j + 1,      k))
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      intfx = c00 * (1 - dy) + c10 * dy
    #elif threeD
      c000 = 0.25 * (fx(      i,      j,      k) + fx(      i,  j - 1,      k) +&
                   & fx(      i,      j,  k - 1) + fx(      i,  j - 1,  k - 1))
      c100 = 0.25 * (fx(  i + 1,      j,      k) + fx(  i + 1,  j - 1,      k) +&
                   & fx(  i + 1,      j,  k - 1) + fx(  i + 1,  j - 1,  k - 1))
      c001 = 0.25 * (fx(      i,      j,      k) + fx(      i,      j,  k + 1) +&
                   & fx(      i,  j - 1,      k) + fx(      i,  j - 1,  k + 1))
      c101 = 0.25 * (fx(  i + 1,      j,      k) + fx(  i + 1,      j,  k + 1) +&
                   & fx(  i + 1,  j - 1,      k) + fx(  i + 1,  j - 1,  k + 1))
      c010 = 0.25 * (fx(      i,      j,      k) + fx(      i,  j + 1,      k) +&
                   & fx(      i,      j,  k - 1) + fx(      i,  j + 1,  k - 1))
      c110 = 0.25 * (fx(  i + 1,      j,      k) + fx(  i + 1,      j,  k - 1) +&
                   & fx(  i + 1,  j + 1,  k - 1) + fx(  i + 1,  j + 1,      k))
      c011 = 0.25 * (fx(      i,      j,      k) + fx(      i,  j + 1,      k) +&
                   & fx(      i,  j + 1,  k + 1) + fx(      i,      j,  k + 1))
      c111 = 0.25 * (fx(  i + 1,      j,      k) + fx(  i + 1,  j + 1,      k) +&
                   & fx(  i + 1,  j + 1,  k + 1) + fx(  i + 1,      j,  k + 1))
      c00 = c000 * (1 - dx) + c100 * dx
      c01 = c001 * (1 - dx) + c101 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c0 = c00 * (1 - dy) + c10 * dy
      c1 = c01 * (1 - dy) + c11 * dy
      intfx = c0 * (1 - dz) + c1 * dz
    #endif

    ! b_y
    #ifdef oneD
      c0 = 0.5 * (fy(  i - 1,      j,      k) + fy(      i,      j,      k))
      c1 = 0.5 * (fy(      i,      j,      k) + fy(  i + 1,      j,      k))
      intfy = c0 * (1 - dx) + c1 * dx
    #elif twoD
      c000 = 0.5 * (fy(  i - 1,      j,      k) + fy(      i,      j,      k))
      c100 = 0.5 * (fy(      i,      j,      k) + fy(  i + 1,      j,      k))
      c010 = 0.5 * (fy(  i - 1,  j + 1,      k) + fy(      i,  j + 1,      k))
      c110 = 0.5 * (fy(      i,  j + 1,      k) + fy(  i + 1,  j + 1,      k))
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      intfy = c00 * (1 - dy) + c10 * dy
    #elif threeD
      c000 = 0.25 * (fy(  i - 1,      j,  k - 1) + fy(  i - 1,      j,      k) +&
                   & fy(      i,      j,  k - 1) + fy(      i,      j,      k))
      c100 = 0.25 * (fy(      i,      j,  k - 1) + fy(      i,      j,      k) +&
                   & fy(  i + 1,      j,  k - 1) + fy(  i + 1,      j,      k))
      c001 = 0.25 * (fy(  i - 1,      j,      k) + fy(  i - 1,      j,  k + 1) +&
                   & fy(      i,      j,      k) + fy(      i,      j,  k + 1))
      c101 = 0.25 * (fy(      i,      j,      k) + fy(      i,      j,  k + 1) +&
                   & fy(  i + 1,      j,      k) + fy(  i + 1,      j,  k + 1))
      c010 = 0.25 * (fy(  i - 1,  j + 1,  k - 1) + fy(  i - 1,  j + 1,      k) +&
                   & fy(      i,  j + 1,  k - 1) + fy(      i,  j + 1,      k))
      c110 = 0.25 * (fy(      i,  j + 1,  k - 1) + fy(      i,  j + 1,      k) +&
                   & fy(  i + 1,  j + 1,  k - 1) + fy(  i + 1,  j + 1,      k))
      c011 = 0.25 * (fy(  i - 1,  j + 1,      k) + fy(  i - 1,  j + 1,  k + 1) +&
                   & fy(      i,  j + 1,      k) + fy(      i,  j + 1,  k + 1))
      c111 = 0.25 * (fy(      i,  j + 1,      k) + fy(      i,  j + 1,  k + 1) +&
                   & fy(  i + 1,  j + 1,      k) + fy(  i + 1,  j + 1,  k + 1))
      c00 = c000 * (1 - dx) + c100 * dx
      c01 = c001 * (1 - dx) + c101 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c0 = c00 * (1 - dy) + c10 * dy
      c1 = c01 * (1 - dy) + c11 * dy
      intfy = c0 * (1 - dz) + c1 * dz
    #endif

    ! b_z
    #ifdef oneD
      c0 = 0.5 * (fz(  i - 1,      j,      k) + fz(      i,      j,      k))
      c1 = 0.5 * (fz(      i,      j,      k) + fz(  i + 1,      j,      k))
      intfz = c0 * (1 - dx) + c1 * dx
    #elif twoD
      c000 = 0.25 * (fz(  i - 1,  j - 1,      k) + fz(  i - 1,      j,      k) +&
                   & fz(      i,  j - 1,      k) + fz(      i,      j,      k))
      c100 = 0.25 * (fz(      i,  j - 1,      k) + fz(      i,      j,      k) +&
                   & fz(  i + 1,  j - 1,      k) + fz(  i + 1,      j,      k))
      c010 = 0.25 * (fz(  i - 1,      j,      k) + fz(  i - 1,  j + 1,      k) +&
                   & fz(      i,      j,      k) + fz(      i,  j + 1,      k))
      c110 = 0.25 * (fz(      i,      j,      k) + fz(      i,  j + 1,      k) +&
                   & fz(  i + 1,      j,      k) + fz(  i + 1,  j + 1,      k))
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      intfz = c00 * (1 - dy) + c10 * dy
    #elif threeD
      c000 = 0.25 * (fz(  i - 1,  j - 1,      k) + fz(  i - 1,      j,      k) +&
                   & fz(      i,  j - 1,      k) + fz(      i,      j,      k))
      c100 = 0.25 * (fz(      i,  j - 1,      k) + fz(      i,      j,      k) +&
                   & fz(  i + 1,  j - 1,      k) + fz(  i + 1,      j,      k))
      c001 = 0.25 * (fz(  i - 1,  j - 1,  k + 1) + fz(  i - 1,      j,  k + 1) +&
                   & fz(      i,  j - 1,  k + 1) + fz(      i,      j,  k + 1))
      c101 = 0.25 * (fz(      i,  j - 1,  k + 1) + fz(      i,      j,  k + 1) +&
                   & fz(  i + 1,  j - 1,  k + 1) + fz(  i + 1,      j,  k + 1))
      c010 = 0.25 * (fz(  i - 1,      j,      k) + fz(  i - 1,  j + 1,      k) +&
                   & fz(      i,      j,      k) + fz(      i,  j + 1,      k))
      c110 = 0.25 * (fz(      i,      j,      k) + fz(      i,  j + 1,      k) +&
                   & fz(  i + 1,      j,      k) + fz(  i + 1,  j + 1,      k))
      c011 = 0.25 * (fz(  i - 1,      j,  k + 1) + fz(  i - 1,  j + 1,  k + 1) +&
                   & fz(      i,      j,  k + 1) + fz(      i,  j + 1,  k + 1))
      c111 = 0.25 * (fz(      i,      j,  k + 1) + fz(      i,  j + 1,  k + 1) +&
                   & fz(  i + 1,      j,  k + 1) + fz(  i + 1,  j + 1,  k + 1))
      c00 = c000 * (1 - dx) + c100 * dx
      c01 = c001 * (1 - dx) + c101 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c0 = c00 * (1 - dy) + c10 * dy
      c1 = c01 * (1 - dy) + c11 * dy
      intfz = c0 * (1 - dz) + c1 * dz
    #endif
  end subroutine interpFromFaces

  subroutine depositCurrentsFromSingleParticle(s, tile, p, x1, y1, z1,&
                                                         & x2, y2, z2, multiplier)
    implicit none
    type(particle_tile), intent(in) :: tile
    integer, intent(in)       :: s, p
    real, intent(in)          :: x1, y1, z1, x2, y2, z2
    real, optional, intent(in):: multiplier
    real                      :: xr, yr, zr
    integer(kind=2)           :: i1, i2, j1, j2, k1, k2
    integer(kind=2)           :: i1p1, i2p1, j1p1, j2p1, k1p1, k2p1
    real                      :: Wx1, Wy1, Wz1, Wx2, Wy2, Wz2
    real                      :: onemWx1, onemWy1, onemWz1, onemWx2, onemWy2, onemWz2
    real                      :: Fx1, Fy1, Fz1, Fx2, Fy2, Fz2
    real                      :: weighted_charge

    if (present(multiplier)) then
      weighted_charge = multiplier * tile%weight(p) * species(s)%ch_sp * unit_ch / B_norm
    else
      weighted_charge = tile%weight(p) * species(s)%ch_sp * unit_ch / B_norm
    end if

    #ifdef oneD
      i1 = FLOOR(x1);   i2 = FLOOR(x2)
      j1 = 0;           j2 = 0
      k1 = 0;           k2 = 0
      i1p1 = i1 + 1_2;  i2p1 = i2 + 1_2
    #elif twoD
      i1 = FLOOR(x1);  i2 = FLOOR(x2)
      j1 = FLOOR(y1);  j2 = FLOOR(y2)
      k1 = 0;          k2 = 0
      i1p1 = i1 + 1_2;  i2p1 = i2 + 1_2
      j1p1 = j1 + 1_2;  j2p1 = j2 + 1_2
    #elif threeD
      i1 = FLOOR(x1);  i2 = FLOOR(x2)
      j1 = FLOOR(y1);  j2 = FLOOR(y2)
      k1 = FLOOR(z1);  k2 = FLOOR(z2)
      i1p1 = i1 + 1_2;  i2p1 = i2 + 1_2
      j1p1 = j1 + 1_2;  j2p1 = j2 + 1_2
      k1p1 = k1 + 1_2;  k2p1 = k2 + 1_2
    #endif

    ! this "function" takes
    ! ... the start and end coordinates: `x1`, `x2`, `y1`, `y2`, `z1`, `z2` ...
    ! ... the start and end cells: `i1`, `i2`, `j1`, `j2`, `k1`, `k2` ...
    ! ... the start and end cells + 1: `i1p1`, `i2p1` etc ...
    ! ... the weighted_chargeed charge: `weighted_charge = weight * charge_sp * unit_charge / Bnorm`
    ! ... and deposits proper currents to corresponding components
    include "zigzag_deposit.F"
  end subroutine depositCurrentsFromSingleParticle
end module m_helpers
