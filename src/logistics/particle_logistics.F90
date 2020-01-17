#include "../defs.F90"

module m_particlelogistics
  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_errors
  use m_domain
  use m_particles
contains
  subroutine copyParticleFromTo(s, p_from, p_to, ti, tj, tk)
    ! DEP_PRT [particle-dependent]
    implicit none
    ! within a single tile
    integer, intent(in)   :: s, p_from, p_to, ti, tj, tk
    species(s)%prtl_tile(ti, tj, tk)%xi(p_to) = species(s)%prtl_tile(ti, tj, tk)%xi(p_from)
    species(s)%prtl_tile(ti, tj, tk)%yi(p_to) = species(s)%prtl_tile(ti, tj, tk)%yi(p_from)
    species(s)%prtl_tile(ti, tj, tk)%zi(p_to) = species(s)%prtl_tile(ti, tj, tk)%zi(p_from)

    species(s)%prtl_tile(ti, tj, tk)%dx(p_to) = species(s)%prtl_tile(ti, tj, tk)%dx(p_from)
    species(s)%prtl_tile(ti, tj, tk)%dy(p_to) = species(s)%prtl_tile(ti, tj, tk)%dy(p_from)
    species(s)%prtl_tile(ti, tj, tk)%dz(p_to) = species(s)%prtl_tile(ti, tj, tk)%dz(p_from)

    species(s)%prtl_tile(ti, tj, tk)%u(p_to) = species(s)%prtl_tile(ti, tj, tk)%u(p_from)
    species(s)%prtl_tile(ti, tj, tk)%v(p_to) = species(s)%prtl_tile(ti, tj, tk)%v(p_from)
    species(s)%prtl_tile(ti, tj, tk)%w(p_to) = species(s)%prtl_tile(ti, tj, tk)%w(p_from)

    species(s)%prtl_tile(ti, tj, tk)%ind(p_to) = species(s)%prtl_tile(ti, tj, tk)%ind(p_from)
    species(s)%prtl_tile(ti, tj, tk)%proc(p_to) = species(s)%prtl_tile(ti, tj, tk)%proc(p_from)
  end subroutine copyParticleFromTo

  subroutine removeParticleFromTile(s, ti, tj, tk, p)
    implicit none
    integer, intent(in)   :: s, ti, tj, tk, p
    if (p .ne. species(s)%prtl_tile(ti, tj, tk)%npart_sp) then
      call copyParticleFromTo(s, species(s)%prtl_tile(ti, tj, tk)%npart_sp, p, ti, tj, tk)
    end if
    species(s)%prtl_tile(ti, tj, tk)%npart_sp = species(s)%prtl_tile(ti, tj, tk)%npart_sp - 1
  end subroutine removeParticleFromTile

  subroutine createParticle(s, xi, yi, zi, dx, dy, dz, u, v, w, &
                          & ind, proc)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in)           :: s
    integer(kind=2), intent(in)   :: xi, yi, zi
    real, intent(in)              :: dx, dy, dz, u, v, w
    integer                       :: p
    integer                       :: ti, tj, tk
    integer, optional, intent(in) :: ind, proc
    ti = INT(FLOOR(REAL(xi) / REAL(species(s)%tile_sx))) + 1
    tj = INT(FLOOR(REAL(yi) / REAL(species(s)%tile_sy))) + 1
    tk = INT(FLOOR(REAL(zi) / REAL(species(s)%tile_sz))) + 1
    #ifdef DEBUG
      if ((ti .gt. species(s)%tile_nx) .or. &
        & (tj .gt. species(s)%tile_ny) .or. &
        & (tk .gt. species(s)%tile_nz)) then
        print *, mpi_rank, xi, yi, zi, ti, tj, tk
        print *, species(s)%tile_nx, species(s)%tile_ny, species(s)%tile_nz
        call throwError('ERROR: wrong ti, tj, tk in `createParticle`')
      end if
      if ((xi .lt. species(s)%prtl_tile(ti, tj, tk)%x1) .or. &
        & (xi .ge. species(s)%prtl_tile(ti, tj, tk)%x2) .or. &
        & (yi .lt. species(s)%prtl_tile(ti, tj, tk)%y1) .or. &
        & (yi .ge. species(s)%prtl_tile(ti, tj, tk)%y2) .or. &
        & (zi .lt. species(s)%prtl_tile(ti, tj, tk)%z1) .or. &
        & (zi .ge. species(s)%prtl_tile(ti, tj, tk)%z2)) then
        print *, xi, yi, zi, dx, dy, dz
        print *, species(s)%prtl_tile(ti, tj, tk)%x1,&
               & species(s)%prtl_tile(ti, tj, tk)%x2,&
               & species(s)%prtl_tile(ti, tj, tk)%y1,&
               & species(s)%prtl_tile(ti, tj, tk)%y2,&
               & species(s)%prtl_tile(ti, tj, tk)%z1,&
               & species(s)%prtl_tile(ti, tj, tk)%z2
        print *, ti, tj, tk
        print *, species(s)%tile_nx, species(s)%tile_ny, species(s)%tile_nz
        print *, species(s)%tile_sx, species(s)%tile_sy, species(s)%tile_sz
        call throwError('ERROR: wrong ti, tj, tk in `createParticle` according to x1,x2,etc')
      end if
    #endif
    species(s)%prtl_tile(ti, tj, tk)%npart_sp = species(s)%prtl_tile(ti, tj, tk)%npart_sp + 1
    p = species(s)%prtl_tile(ti, tj, tk)%npart_sp

    species(s)%prtl_tile(ti, tj, tk)%xi(p) = xi
    species(s)%prtl_tile(ti, tj, tk)%dx(p) = dx

    species(s)%prtl_tile(ti, tj, tk)%yi(p) = yi
    species(s)%prtl_tile(ti, tj, tk)%dy(p) = dy

    species(s)%prtl_tile(ti, tj, tk)%zi(p) = zi
    species(s)%prtl_tile(ti, tj, tk)%dz(p) = dz

    species(s)%prtl_tile(ti, tj, tk)%u(p) = u
    species(s)%prtl_tile(ti, tj, tk)%v(p) = v
    species(s)%prtl_tile(ti, tj, tk)%w(p) = w

    if (present(ind) .and. present(proc)) then
      ! moving particle from one tile/meshblock to another
      species(s)%prtl_tile(ti, tj, tk)%ind(p) = ind
      #ifdef DEBUG
        species(s)%prtl_tile(ti, tj, tk)%proc(p) = mpi_rank
      #else
        species(s)%prtl_tile(ti, tj, tk)%proc(p) = proc
      #endif
    else
      ! create a very new particle
      species(s)%prtl_tile(ti, tj, tk)%ind(p) = species(s)%cntr_sp
      species(s)%prtl_tile(ti, tj, tk)%proc(p) = mpi_rank
      species(s)%cntr_sp = species(s)%cntr_sp + 1
    end if
  end subroutine createParticle

  subroutine allocateParticles(prt, sz)
    ! DEP_PRT [particle-dependent]
    implicit none
    type(particle_tile), intent(inout)    :: prt
    integer, intent(in)                   :: sz
    if (allocated(prt%xi)) deallocate(prt%xi)
    if (allocated(prt%yi)) deallocate(prt%yi)
    if (allocated(prt%zi)) deallocate(prt%zi)
    if (allocated(prt%dx)) deallocate(prt%dx)
    if (allocated(prt%dy)) deallocate(prt%dy)
    if (allocated(prt%dz)) deallocate(prt%dz)
    if (allocated(prt%u)) deallocate(prt%u)
    if (allocated(prt%v)) deallocate(prt%v)
    if (allocated(prt%w)) deallocate(prt%w)
    allocate(prt%xi(sz)); allocate(prt%yi(sz)); allocate(prt%zi(sz))
    allocate(prt%dx(sz)); allocate(prt%dy(sz)); allocate(prt%dz(sz))
    allocate(prt%u(sz)); allocate(prt%v(sz)); allocate(prt%w(sz))
    allocate(prt%ind(sz)); allocate(prt%proc(sz))
  end subroutine allocateParticles

  subroutine checkTileSizes()
    implicit none
    integer     :: s, ti, tj, tk
    do s = 1, nspec ! loop over species
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            if (species(s)%prtl_tile(ti, tj, tk)%npart_sp .ge. species(s)%prtl_tile(ti, tj, tk)%maxptl_sp * 0.7) then
              ! increase the tile size
              if (resize_tiles) then
                call reallocTileSize(species(s)%prtl_tile(ti, tj, tk), .true.)
              else
                print *, "DANGER: `maxptl_sp` in tiles too low, consider increasing `maxptl` or turning on `resize_tiles`."
              end if
            else if ((species(s)%prtl_tile(ti, tj, tk)%npart_sp .lt. (species(s)%prtl_tile(ti, tj, tk)%maxptl_sp * 0.3)) .and.&
                   & ((species(s)%prtl_tile(ti, tj, tk)%maxptl_sp * 0.5) .gt. min_tile_nprt)) then
              ! decrease the tile size
              if (resize_tiles) then
                call reallocTileSize(species(s)%prtl_tile(ti, tj, tk), .false.)
              end if
            end if
          end do
        end do
      end do
    end do
  end subroutine checkTileSizes

  subroutine reallocTileSize(tile, increase_flag)
    ! DEP_PRT [particle-dependent]
    implicit none
    type(particle_tile), intent(inout)          :: tile
    ! `.true.` if need to increase, otherwise `.false.`
    logical, intent(in)                         :: increase_flag
    integer(kind=2), allocatable, dimension(:)  :: dummy_int2
    integer, allocatable, dimension(:)          :: dummy_int
    real, allocatable, dimension(:)             :: dummy_real
    integer                                     :: current_npart

    if (increase_flag) then
      ! increase twice
      tile%maxptl_sp = INT(tile%maxptl_sp * 1.5)
    else
      ! decrease twice
      tile%maxptl_sp = INT(tile%maxptl_sp * 0.5)
    end if

    if (tile%npart_sp .gt. tile%maxptl_sp) then
      call throwError('ERROR: `npart > maxptl` in `reallocTileSize`')
    end if

    allocate(dummy_int2(tile%maxptl_sp))
    allocate(dummy_int(tile%maxptl_sp))
    allocate(dummy_real(tile%maxptl_sp))

    current_npart = tile%npart_sp

    dummy_int2(1 : current_npart) = tile%xi(1 : current_npart)
    deallocate(tile%xi); allocate(tile%xi(tile%maxptl_sp))
    tile%xi(1 : current_npart) = dummy_int2(1 : current_npart)

    dummy_int2(1 : current_npart) = tile%yi(1 : current_npart)
    deallocate(tile%yi); allocate(tile%yi(tile%maxptl_sp))
    tile%yi(1 : current_npart) = dummy_int2(1 : current_npart)

    dummy_int2(1 : current_npart) = tile%zi(1 : current_npart)
    deallocate(tile%zi); allocate(tile%zi(tile%maxptl_sp))
    tile%zi(1 : current_npart) = dummy_int2(1 : current_npart)

    dummy_real(1 : current_npart) = tile%dx(1 : current_npart)
    deallocate(tile%dx); allocate(tile%dx(tile%maxptl_sp))
    tile%dx(1 : current_npart) = dummy_real(1 : current_npart)

    dummy_real(1 : current_npart) = tile%dy(1 : current_npart)
    deallocate(tile%dy); allocate(tile%dy(tile%maxptl_sp))
    tile%dy(1 : current_npart) = dummy_real(1 : current_npart)

    dummy_real(1 : current_npart) = tile%dz(1 : current_npart)
    deallocate(tile%dz); allocate(tile%dz(tile%maxptl_sp))
    tile%dz(1 : current_npart) = dummy_real(1 : current_npart)

    dummy_real(1 : current_npart) = tile%u(1 : current_npart)
    deallocate(tile%u); allocate(tile%u(tile%maxptl_sp))
    tile%u(1 : current_npart) = dummy_real(1 : current_npart)

    dummy_real(1 : current_npart) = tile%v(1 : current_npart)
    deallocate(tile%v); allocate(tile%v(tile%maxptl_sp))
    tile%v(1 : current_npart) = dummy_real(1 : current_npart)

    dummy_real(1 : current_npart) = tile%w(1 : current_npart)
    deallocate(tile%w); allocate(tile%w(tile%maxptl_sp))
    tile%w(1 : current_npart) = dummy_real(1 : current_npart)

    dummy_int(1 : current_npart) = tile%ind(1 : current_npart)
    deallocate(tile%ind); allocate(tile%ind(tile%maxptl_sp))
    tile%ind(1 : current_npart) = dummy_int(1 : current_npart)

    dummy_int(1 : current_npart) = tile%proc(1 : current_npart)
    deallocate(tile%proc); allocate(tile%proc(tile%maxptl_sp))
    tile%proc(1 : current_npart) = dummy_int(1 : current_npart)

    deallocate(dummy_int2)
    deallocate(dummy_int)
    deallocate(dummy_real)
    if (increase_flag) then
      call printDiag(.true., "...reallocTileSize(+).."//trim(STR(tile%npart_sp))&
                    & //".."//trim(STR(tile%maxptl_sp)),&
                    & .true.)
    else
      call printDiag(.true., "...reallocTileSize(-).."//trim(STR(tile%npart_sp))&
                    & //".."//trim(STR(tile%maxptl_sp)),&
                    & .true.)
    end if
  end subroutine reallocTileSize

  subroutine clearGhostParticles()
    implicit none
    integer                       :: s, p, ti, tj, tk
    do s = 1, nspec ! loop over species
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            ! FIX1 try to vectorize this
            p = 1
            do while (p .le. species(s)%prtl_tile(ti, tj, tk)%npart_sp)
              if (species(s)%prtl_tile(ti, tj, tk)%proc(p) .lt. 0) then
                call removeParticleFromTile(s, ti, tj, tk, p)
              else
                p = p + 1
              end if
            end do ! p
          end do ! tk
        end do ! tj
      end do ! ti
    end do ! s
    call printDiag((mpi_rank .eq. 0), "clearGhostParticles()", .true.)
  end subroutine clearGhostParticles

  subroutine injectParticleGlobally(s, x_glob, y_glob, z_glob, u, v, w)
    implicit none
    integer, intent(in) :: s
    real, intent(in)    :: x_glob, y_glob, z_glob
    real, intent(in)    :: u, v, w
    real                :: x_loc, y_loc, z_loc
    real                :: dx_, dy_, dz_
    integer(kind=2)     :: xi_, yi_, zi_

    ! convert local to global
    call globalToLocalCoords(x_glob, y_glob, z_glob,&
                           & x_loc, y_loc, z_loc)

    x_loc = x_loc + TINYXYZ
    y_loc = y_loc + TINYXYZ
    z_loc = z_loc + TINYXYZ
    ! check if the coordinate is within the current MPI domain
    if ((x_loc .ge. 0.0) .and. (x_loc .lt. REAL(this_meshblock%ptr%sx)) .and.&
      & (y_loc .ge. 0.0) .and. (y_loc .lt. REAL(this_meshblock%ptr%sy)) .and.&
      & (z_loc .ge. 0.0) .and. (z_loc .lt. REAL(this_meshblock%ptr%sz))) then
      ! transform coordinates
      xi_ = INT(FLOOR(x_loc), 2); dx_ = x_loc - FLOOR(x_loc)
      yi_ = INT(FLOOR(y_loc), 2); dy_ = y_loc - FLOOR(y_loc)
      zi_ = INT(FLOOR(z_loc), 2); dz_ = z_loc - FLOOR(z_loc)

      call createParticle(s, xi_, yi_, zi_, dx_, dy_, dz_, u, v, w)
    end if
  end subroutine
end module m_particlelogistics
