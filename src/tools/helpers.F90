#include "../defs.F90"

module m_helpers
  use m_globalnamespace
  use m_aux
  use m_domain
  use m_particles
  use m_fields
  implicit none
contains
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
                               & x_loc, y_loc, z_loc, adjustQ_)
    implicit none
    real, intent(in)              :: x_glob, y_glob, z_glob
    real, intent(out)             :: x_loc, y_loc, z_loc
    logical, optional, intent(in) :: adjustQ_
    logical                       :: adjustQ
    if (present(adjustQ_)) then
      adjustQ = adjustQ_
    else
      adjustQ = .false.
    end if

    if (adjustQ) then
      x_loc = MAX(0.0, MIN(x_glob - REAL(this_meshblock%ptr%x0), REAL(this_meshblock%ptr%sx)))
      y_loc = MAX(0.0, MIN(y_glob - REAL(this_meshblock%ptr%y0), REAL(this_meshblock%ptr%sy)))
      #ifdef threeD
        z_loc = MAX(0.0, MIN(z_glob - REAL(this_meshblock%ptr%z0), REAL(this_meshblock%ptr%sz)))
      #else
        z_loc = z_glob
      #endif
    else
      x_loc = x_glob - REAL(this_meshblock%ptr%x0)
      y_loc = y_glob - REAL(this_meshblock%ptr%y0)
      #ifdef threeD
        z_loc = z_glob - REAL(this_meshblock%ptr%z0)
      #else
        z_loc = z_glob
      #endif
    end if
  end subroutine globalToLocalCoords

  subroutine generateCoordInRegion(xmin, xmax, ymin, ymax, zmin, zmax,&
                                 & x_, y_, z_, xi_, yi_, zi_, dx_, dy_, dz_)
    implicit none
    real                          :: rnd
    real, intent(in)              :: xmin, xmax, ymin, ymax, zmin, zmax
    real, intent(out)             :: x_, y_, z_, dx_, dy_, dz_
    integer(kind=2), intent(out)  :: xi_, yi_, zi_

    rnd = random(dseed)
    x_ = xmin + rnd * (xmax - xmin)
    xi_ = INT(FLOOR(x_), 2); dx_ = x_ - FLOOR(x_)
    if (xi_ .eq. this_meshblock%ptr%sx) then
      xi_ = xi_ - 1; dx_ = dx_ + 1.0
    end if
    rnd = random(dseed)
    y_ = ymin + rnd * (ymax - ymin)
    yi_ = INT(FLOOR(y_), 2); dy_ = y_ - FLOOR(y_)
    if (yi_ .eq. this_meshblock%ptr%sy) then
      yi_ = yi_ - 1; dy_ = dy_ + 1.0
    end if
    #ifdef threeD
      rnd = random(dseed)
      z_ = zmin + rnd * (zmax - zmin)
      zi_ = INT(FLOOR(z_), 2); dz_ = z_ - FLOOR(z_)
      if (zi_ .eq. this_meshblock%ptr%sz) then
        zi_ = zi_ - 1; dz_ = dz_ + 1.0
      end if
    #else
      z_ = 0.5
      zi_ = 0; dz_ = 0.5
    #endif
  end subroutine

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
          #ifndef threeD
            if (ind3 .ne. 0) cycle
          #endif
          if (.not. associated(this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr)) cycle
          cntr = cntr + 1
        end do
      end do
    end do
    sendrecv_neighbors = cntr
  end subroutine computeNumberOfNeighbors

  subroutine computeDensity(s, reset)
    implicit none
    integer, intent(in)                   :: s
    logical, intent(in)                   :: reset
    integer                               :: p, ti, tj, tk
    integer(kind=2), pointer, contiguous  :: pt_xi(:), pt_yi(:), pt_zi(:)
    integer(kind=2) :: i, j, k
    integer :: i1, i2, j1, j2, k1, k2, ds
    integer :: pow
    ds = 2
    #ifndef threeD
      pow = 2
    #else
      pow = 3
    #endif
    if (reset) then
      lg_arr(:,:,:) = 0
    end if
    do ti = 1, species(s)%tile_nx
      do tj = 1, species(s)%tile_ny
        do tk = 1, species(s)%tile_nz
          pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
          pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
          pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi
          ! FIX1 vectorize/align
          do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
            i = pt_xi(p); j = pt_yi(p); k = pt_zi(p)

            i1 = max(i - ds, -NGHOST)
            i2 = min(i + ds, this_meshblock%ptr%sx + NGHOST - 1)

            j1 = max(j - ds, -NGHOST)
            j2 = min(j + ds, this_meshblock%ptr%sy + NGHOST - 1)

            #ifndef threeD
              k1 = 0; k2 = 0
            #else
              k1 = max(k - ds, -NGHOST)
              k2 = min(k + ds, this_meshblock%ptr%sz + NGHOST - 1)
            #endif

            do k = k1, k2
              do j = j1, j2
                do i = i1, i2
                  lg_arr(i, j, k) = lg_arr(i, j, k) + 1.0 / (2 * ds + 1.0)**pow
                end do
              end do
            end do

          end do
          pt_xi => null(); pt_yi => null(); pt_zi => null()
        end do
      end do
    end do
  end subroutine computeDensity

  subroutine computeEnergy(s, reset)
    implicit none
    integer, intent(in)                   :: s
    logical, intent(in)                   :: reset
    integer                               :: p, ti, tj, tk
    integer(kind=2), pointer, contiguous  :: pt_xi(:), pt_yi(:), pt_zi(:)
    real, pointer, contiguous             :: pt_u(:), pt_v(:), pt_w(:)
    integer(kind=2) :: i, j, k
    integer :: i1, i2, j1, j2, k1, k2, ds
    integer :: pow
    logical :: massive
    real    :: energy

    if (species(s)%m_sp .gt. 0) then
      massive = .true.
    else
      massive = .false.
    end if

    ds = 2
    #ifndef threeD
      pow = 2
    #else
      pow = 3
    #endif
    if (reset) then
      lg_arr(:,:,:) = 0
    end if
    do ti = 1, species(s)%tile_nx
      do tj = 1, species(s)%tile_ny
        do tk = 1, species(s)%tile_nz
          pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
          pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
          pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi
          pt_u => species(s)%prtl_tile(ti, tj, tk)%u
          pt_v => species(s)%prtl_tile(ti, tj, tk)%v
          pt_w => species(s)%prtl_tile(ti, tj, tk)%w
          ! FIX1 vectorize/align
          do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
            i = pt_xi(p); j = pt_yi(p); k = pt_zi(p)
            if (massive) then
              energy = sqrt(1.0 + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
            else
              energy = sqrt(pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
            end if

            i1 = max(i - ds, -NGHOST)
            i2 = min(i + ds, this_meshblock%ptr%sx + NGHOST - 1)

            j1 = max(j - ds, -NGHOST)
            j2 = min(j + ds, this_meshblock%ptr%sy + NGHOST - 1)

            #ifndef threeD
              k1 = 0; k2 = 0
            #else
              k1 = max(k - ds, -NGHOST)
              k2 = min(k + ds, this_meshblock%ptr%sz + NGHOST - 1)
            #endif

            do k = k1, k2
              do j = j1, j2
                do i = i1, i2
                  lg_arr(i, j, k) = lg_arr(i, j, k) + energy / (2 * ds + 1.0)**pow
                end do
              end do
            end do

          end do
          pt_xi => null(); pt_yi => null(); pt_zi => null()
          pt_u => null(); pt_v => null(); pt_w => null()
        end do
      end do
    end do
  end subroutine computeEnergy

  subroutine interpFromEdges(dx, dy, dz, i, j, k, &
                           & fx, fy, fz, &
                           & intfx, intfy, intfz)
    !$omp declare simd(interpFromEdges)
    implicit none
    integer(kind=2), intent(in)   :: i, j, k
    real, intent(in)              :: dx, dy, dz
    real, intent(in)              :: fx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                                      & fldBoundZ)
    real, intent(in)              :: fy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                                      & fldBoundZ)
    real, intent(in)              :: fz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                                      & fldBoundZ)
    real, intent(out)             :: intfx, intfy, intfz
    real                          :: c000, c100, c001, c101, c010, c110, c011, c111,&
                                   & c00, c01, c10, c11, c0, c1
    ! f_x
    c000 = 0.5 * (fx(     i,      j,      k) + fx(  i - 1,      j,      k))
    c100 = 0.5 * (fx(     i,      j,      k) + fx(  i + 1,      j,      k))
    c010 = 0.5 * (fx(     i,  j + 1,      k) + fx(  i - 1,  j + 1,      k))
    c110 = 0.5 * (fx(     i,  j + 1,      k) + fx(  i + 1,  j + 1,      k))
    c00 = c000 * (1 - dx) + c100 * dx
    c10 = c010 * (1 - dx) + c110 * dx
    c0 = c00 * (1 - dy) + c10 * dy
    #ifndef threeD
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

    ! f_y
    c000 = 0.5 * (fy(     i,      j,      k) + fy(      i,  j - 1,      k))
    c100 = 0.5 * (fy( i + 1,      j,      k) + fy(  i + 1,  j - 1,      k))
    c010 = 0.5 * (fy(     i,      j,      k) + fy(      i,  j + 1,      k))
    c110 = 0.5 * (fy( i + 1,      j,      k) + fy(  i + 1,  j + 1,      k))
    c00 = c000 * (1 - dx) + c100 * dx
    c10 = c010 * (1 - dx) + c110 * dx
    c0 = c00 * (1 - dy) + c10 * dy
    #ifndef threeD
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

    ! f_z
    #ifndef threeD
      c000 = fz(     i,      j,      k)
      c100 = fz( i + 1,      j,      k)
      c010 = fz(     i,  j + 1,      k)
      c110 = fz( i + 1,  j + 1,      k)
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      intfz = c00 * (1 - dy) + c10 * dy
    #else
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
    !$omp declare simd(interpFromFaces)
    implicit none
    integer(kind=2), intent(in)   :: i, j, k
    real, intent(in)              :: dx, dy, dz
    real, intent(in)              :: fx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                                      & fldBoundZ)
    real, intent(in)              :: fy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                                      & fldBoundZ)
    real, intent(in)              :: fz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                                      & fldBoundZ)
    real, intent(out)             :: intfx, intfy, intfz
    real                          :: c000, c100, c001, c101, c010, c110, c011, c111,&
                                   & c00, c01, c10, c11, c0, c1
    ! f_x
    #ifndef threeD
      c000 = 0.5 * (fx(      i,      j,      k) + fx(      i,  j - 1,      k))
      c100 = 0.5 * (fx(  i + 1,      j,      k) + fx(  i + 1,  j - 1,      k))
      c010 = 0.5 * (fx(      i,      j,      k) + fx(      i,  j + 1,      k))
      c110 = 0.5 * (fx(  i + 1,      j,      k) + fx(  i + 1,  j + 1,      k))
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      intfx = c00 * (1 - dy) + c10 * dy
    #else
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
    #ifndef threeD
      c000 = 0.5 * (fy(  i - 1,      j,      k) + fy(      i,      j,      k))
      c100 = 0.5 * (fy(      i,      j,      k) + fy(  i + 1,      j,      k))
      c010 = 0.5 * (fy(  i - 1,  j + 1,      k) + fy(      i,  j + 1,      k))
      c110 = 0.5 * (fy(      i,  j + 1,      k) + fy(  i + 1,  j + 1,      k))
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      intfy = c00 * (1 - dy) + c10 * dy
    #else
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
    #ifndef threeD
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
    #else
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
end module m_helpers
