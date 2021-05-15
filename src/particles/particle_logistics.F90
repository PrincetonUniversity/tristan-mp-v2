#include "../defs.F90"

module m_particlelogistics
  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_errors
  use m_domain
  use m_particles
contains
  ! Subroutine to move particles around
  subroutine createParticleFromAttributes(s, xi, yi, zi, dx, dy, dz,&
                                           & u, v, w,&
                                           #ifdef PRTLPAYLOADS
                                            & payload1, payload2, payload3,&
                                           #endif
                                           & ind, proc, weight)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in)                     :: s
    integer(kind=2), intent(in)             :: xi, yi, zi
    real, intent(in)                        :: dx, dy, dz, u, v, w

    #ifdef PRTLPAYLOADS
      real                                    :: payload1, payload2, payload3
    #endif

    integer                                 :: p
    integer                                 :: ti, tj, tk
    integer, intent(in)                     :: ind, proc
    real                                    :: weight
    ti = 1; tj = 1; tk = 1
    #if defined(oneD) || defined (twoD) || defined (threeD)
      ti = FLOOR(REAL(xi) / REAL(species(s)%tile_sx)) + 1
    #endif

    #if defined (twoD) || defined (threeD)
      tj = FLOOR(REAL(yi) / REAL(species(s)%tile_sy)) + 1
    #endif

    #if defined (threeD)
      tk = FLOOR(REAL(zi) / REAL(species(s)%tile_sz)) + 1
    #endif
    ! if the debug flag is enabled ...
    ! ... check that the particle is within the boundaries ...
    ! ... of that tile and that the tile exists
    #ifdef DEBUG
      if ((s .le. 0) .or. (s .gt. nspec)) then
        call throwError('Wrong species in `createParticleFromAttributes`.')
      end if
      if ((ti .gt. species(s)%tile_nx) .or. &
        & (tj .gt. species(s)%tile_ny) .or. &
        & (tk .gt. species(s)%tile_nz)) then
        print *, mpi_rank, xi, yi, zi, ti, tj, tk
        print *, species(s)%tile_nx, species(s)%tile_ny, species(s)%tile_nz
        call throwError('ERROR: wrong ti, tj, tk in `createParticleFromAttributes`')
      end if
      if ((xi .lt. species(s)%prtl_tile(ti, tj, tk)%x1) .or. &
        & (xi .ge. species(s)%prtl_tile(ti, tj, tk)%x2) .or. &
        & (yi .lt. species(s)%prtl_tile(ti, tj, tk)%y1) .or. &
        & (yi .ge. species(s)%prtl_tile(ti, tj, tk)%y2) .or. &
        & (zi .lt. species(s)%prtl_tile(ti, tj, tk)%z1) .or. &
        & (zi .ge. species(s)%prtl_tile(ti, tj, tk)%z2)) then
        print *, s, xi, yi, zi, dx, dy, dz
        print *, species(s)%prtl_tile(ti, tj, tk)%x1,&
               & species(s)%prtl_tile(ti, tj, tk)%x2,&
               & species(s)%prtl_tile(ti, tj, tk)%y1,&
               & species(s)%prtl_tile(ti, tj, tk)%y2,&
               & species(s)%prtl_tile(ti, tj, tk)%z1,&
               & species(s)%prtl_tile(ti, tj, tk)%z2
        print *, ti, tj, tk
        print *, species(s)%tile_nx, species(s)%tile_ny, species(s)%tile_nz
        print *, species(s)%tile_sx, species(s)%tile_sy, species(s)%tile_sz
        call throwError('ERROR: wrong ti, tj, tk in `createParticleFromAttributes` according to x1,x2,etc')
      end if
    #endif
    if (species(s)%prtl_tile(ti, tj, tk)%npart_sp .eq. species(s)%prtl_tile(ti, tj, tk)%maxptl_sp) then
      call throwError('ERROR: npart_sp > maxptl_sp in createParticleFromAttributes')
    end if
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

    species(s)%prtl_tile(ti, tj, tk)%ind(p) = ind
    species(s)%prtl_tile(ti, tj, tk)%proc(p) = proc

    species(s)%prtl_tile(ti, tj, tk)%weight(p) = weight

    #ifdef PRTLPAYLOADS
      species(s)%prtl_tile(ti, tj, tk)%payload1(p) = payload1
      species(s)%prtl_tile(ti, tj, tk)%payload2(p) = payload2
      species(s)%prtl_tile(ti, tj, tk)%payload3(p) = payload3
    #endif
  end subroutine createParticleFromAttributes

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

    species(s)%prtl_tile(ti, tj, tk)%weight(p_to) = species(s)%prtl_tile(ti, tj, tk)%weight(p_from)

    #ifdef PRTLPAYLOADS
      species(s)%prtl_tile(ti, tj, tk)%payload1(p_to) = species(s)%prtl_tile(ti, tj, tk)%payload1(p_from)
      species(s)%prtl_tile(ti, tj, tk)%payload2(p_to) = species(s)%prtl_tile(ti, tj, tk)%payload2(p_from)
      species(s)%prtl_tile(ti, tj, tk)%payload3(p_to) = species(s)%prtl_tile(ti, tj, tk)%payload3(p_from)
    #endif
  end subroutine copyParticleFromTo

  subroutine createEmptyTile(s, ti, tj, tk, maxptl)
    implicit none
    integer, intent(in) :: s, ti, tj, tk, maxptl
    integer             :: maxptl_on_tile
    maxptl_on_tile = maxptl / (species(s)%tile_nx * species(s)%tile_ny * species(s)%tile_nz)

    species(s)%prtl_tile(ti, tj, tk)%spec = s

    species(s)%prtl_tile(ti, tj, tk)%x1 = (ti - 1) * species(s)%tile_sx
    species(s)%prtl_tile(ti, tj, tk)%x2 = min(ti * species(s)%tile_sx, this_meshblock%ptr%sx)
    species(s)%prtl_tile(ti, tj, tk)%y1 = (tj - 1) * species(s)%tile_sy
    species(s)%prtl_tile(ti, tj, tk)%y2 = min(tj * species(s)%tile_sy, this_meshblock%ptr%sy)
    species(s)%prtl_tile(ti, tj, tk)%z1 = (tk - 1) * species(s)%tile_sz
    species(s)%prtl_tile(ti, tj, tk)%z2 = min(tk * species(s)%tile_sz, this_meshblock%ptr%sz)
    #ifdef DEBUG
      if ((species(s)%prtl_tile(ti, tj, tk)%x1 .eq. 0) .and.&
        & (species(s)%prtl_tile(ti, tj, tk)%x2 .eq. 0) .and.&
        & (species(s)%prtl_tile(ti, tj, tk)%y1 .eq. 0) .and.&
        & (species(s)%prtl_tile(ti, tj, tk)%y2 .eq. 0) .and.&
        & (species(s)%prtl_tile(ti, tj, tk)%z1 .eq. 0) .and.&
        & (species(s)%prtl_tile(ti, tj, tk)%z2 .eq. 0)) then
        print *, ti, tj, tk
        print *, species(s)%prtl_tile(ti, tj, tk)%x1,&
         & species(s)%prtl_tile(ti, tj, tk)%x2,&
         & species(s)%prtl_tile(ti, tj, tk)%y1,&
         & species(s)%prtl_tile(ti, tj, tk)%y2,&
         & species(s)%prtl_tile(ti, tj, tk)%z1,&
         & species(s)%prtl_tile(ti, tj, tk)%z2
       call throwError('ERROR: in `initializeParticles`')
      end if
    #endif
    call allocateParticlesOnEmptyTile(s, species(s)%prtl_tile(ti, tj, tk), maxptl_on_tile)
  end subroutine createEmptyTile

  subroutine allocateParticlesOnEmptyTile(s, tile, sz)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in)                   :: s, sz
    type(particle_tile), intent(inout)    :: tile

    tile%npart_sp = 0
    tile%maxptl_sp = sz

    if (allocated(tile%xi)) deallocate(tile%xi)
    if (allocated(tile%yi)) deallocate(tile%yi)
    if (allocated(tile%zi)) deallocate(tile%zi)
    if (allocated(tile%dx)) deallocate(tile%dx)
    if (allocated(tile%dy)) deallocate(tile%dy)
    if (allocated(tile%dz)) deallocate(tile%dz)
    if (allocated(tile%u)) deallocate(tile%u)
    if (allocated(tile%v)) deallocate(tile%v)
    if (allocated(tile%w)) deallocate(tile%w)
    if (allocated(tile%ind)) deallocate(tile%ind)
    if (allocated(tile%proc)) deallocate(tile%proc)
    if (allocated(tile%weight)) deallocate(tile%weight)
    allocate(tile%xi(sz)); allocate(tile%yi(sz)); allocate(tile%zi(sz))
    allocate(tile%dx(sz)); allocate(tile%dy(sz)); allocate(tile%dz(sz))
    allocate(tile%u(sz)); allocate(tile%v(sz)); allocate(tile%w(sz))
    allocate(tile%ind(sz)); allocate(tile%proc(sz)); allocate(tile%weight(sz))

    #ifdef PRTLPAYLOADS
      if (allocated(tile%payload1)) deallocate(tile%payload1)
      if (allocated(tile%payload2)) deallocate(tile%payload2)
      if (allocated(tile%payload3)) deallocate(tile%payload3)
      allocate(tile%payload1(sz))
      allocate(tile%payload2(sz))
      allocate(tile%payload3(sz))
    #endif
  end subroutine allocateParticlesOnEmptyTile

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

    dummy_real(1 : current_npart) = tile%weight(1 : current_npart)
    deallocate(tile%weight); allocate(tile%weight(tile%maxptl_sp))
    tile%weight(1 : current_npart) = dummy_real(1 : current_npart)

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

    #ifdef PRTLPAYLOADS
      dummy_real(1 : current_npart) = tile%payload1(1 : current_npart)
      deallocate(tile%payload1); allocate(tile%payload1(tile%maxptl_sp))
      tile%payload1(1 : current_npart) = dummy_real(1 : current_npart)

      dummy_real(1 : current_npart) = tile%payload2(1 : current_npart)
      deallocate(tile%payload2); allocate(tile%payload2(tile%maxptl_sp))
      tile%payload2(1 : current_npart) = dummy_real(1 : current_npart)

      dummy_real(1 : current_npart) = tile%payload3(1 : current_npart)
      deallocate(tile%payload3); allocate(tile%payload3(tile%maxptl_sp))
      tile%payload3(1 : current_npart) = dummy_real(1 : current_npart)
    #endif

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

  ! Subroutine to create brand new particles
  subroutine createParticle(s, xi, yi, zi, dx, dy, dz, u, v, w, &
                          & ind, proc, weight)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in)           :: s
    integer(kind=2), intent(in)   :: xi, yi, zi
    real, intent(in)              :: dx, dy, dz, u, v, w
    integer                       :: p
    integer                       :: ti, tj, tk
    integer, optional, intent(in) :: ind, proc
    real, optional                :: weight
    integer                       :: ind_, proc_
    real                          :: weight_
    #ifdef PRTLPAYLOADS
      real                          :: payload1, payload2, payload3
      payload1 = 0.0; payload2 = 0.0; payload3 = 0.0
    #endif
    if (present(ind) .and. present(proc)) then
      ! moving particle from one tile/meshblock to another
      ind_ = ind
      proc_ = proc
    else
      ! create a very new particle
      ind_ = species(s)%cntr_sp
      proc_ = mpi_rank
      species(s)%cntr_sp = species(s)%cntr_sp + 1
    end if

    if (present(weight)) then
      weight_ = weight
    else
      weight_ = 1
    end if
    call createParticleFromAttributes(s, xi=xi, yi=yi, zi=zi, dx=dx, dy=dy, dz=dz,&
                                       & u=u, v=v, w=w,&
                                       #ifdef PRTLPAYLOADS
                                        & payload1=payload1, payload2=payload2, payload3=payload3,&
                                       #endif
                                       & ind=ind_, proc=proc_, weight=weight_)
  end subroutine createParticle

  subroutine injectParticleGlobally(s, x_glob, y_glob, z_glob, u, v, w, weight)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in) :: s
    real, intent(in)    :: x_glob, y_glob, z_glob
    real, intent(in)    :: u, v, w
    real                :: x_g, y_g, z_g, x_loc, y_loc, z_loc
    real                :: dx_, dy_, dz_
    integer(kind=2)     :: xi_, yi_, zi_
    real                :: weight_
    real, optional, intent(in) :: weight
    logical             :: contained_flag

    if (present(weight)) then
      weight_ = weight
    else
      weight_ = 1.0
    end if

    #ifdef oneD
      x_g = x_glob
      y_g = 0.5
      z_g = 0.5
    #elif twoD
      x_g = x_glob
      y_g = y_glob
      z_g = 0.5
    #elif threeD
      x_g = x_glob
      y_g = y_glob
      z_g = z_glob
    #endif

    call globalToLocalCoords(x_g, y_g, z_g, x_loc, y_loc, z_loc, containedQ=contained_flag)
    x_loc = x_loc + TINYXYZ
    y_loc = y_loc + TINYXYZ
    z_loc = z_loc + TINYXYZ
    ! check if the coordinate is within the current MPI domain
    if (contained_flag) then
      ! transform coordinates
      call localToCellBasedCoords(x_loc, y_loc, z_loc, xi_, yi_, zi_, dx_, dy_, dz_)
      call createParticle(s, xi_, yi_, zi_, dx_, dy_, dz_, u, v, w, weight=weight_)
    end if
  end subroutine injectParticleGlobally

  subroutine injectParticleLocally(s, x_loc, y_loc, z_loc,&
                                 & u, v, w, ind, proc, weight)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in)                   :: s
    real, intent(in)                      :: x_loc, y_loc, z_loc, u, v, w
    integer, optional, intent(in)         :: ind, proc
    real, optional, intent(in)            :: weight
    integer(kind=2)                       :: xi, yi, zi
    real                                  :: dx, dy, dz

    call localToCellBasedCoords(x_loc, y_loc, z_loc, xi, yi, zi, dx, dy, dz)
    if (present(ind) .and. present(proc) .and. present(weight)) then
      call createParticle(s, xi, yi, zi, dx, dy, dz, u, v, w, ind=ind, proc=proc, weight=weight)
    else if (present(ind) .and. present(proc)) then
      call createParticle(s, xi, yi, zi, dx, dy, dz, u, v, w, ind=ind, proc=proc)
    else if (present(weight)) then
      call createParticle(s, xi, yi, zi, dx, dy, dz, u, v, w, weight=weight)
    else
      call createParticle(s, xi, yi, zi, dx, dy, dz, u, v, w)
    end if
  end subroutine injectParticleLocally

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

  subroutine removeParticleFromTile(s, ti, tj, tk, p)
    implicit none
    integer, intent(in)   :: s, ti, tj, tk, p
    if (p .ne. species(s)%prtl_tile(ti, tj, tk)%npart_sp) then
      call copyParticleFromTo(s, species(s)%prtl_tile(ti, tj, tk)%npart_sp, p, ti, tj, tk)
    end if
    species(s)%prtl_tile(ti, tj, tk)%npart_sp = species(s)%prtl_tile(ti, tj, tk)%npart_sp - 1
  end subroutine removeParticleFromTile

end module m_particlelogistics
