module m_particlelogistics
  use m_globalnamespace
  use m_readinput, only: getInput
  use m_aux
  use m_helpers
  use m_errors
  use m_domain
  use m_particles
contains
  subroutine initializeParticles()
    implicit none
    integer :: s, ti, tj, tk
    character(len=STR_MAX) :: var_name

    allocate (species(nspec))
    do s = 1, nspec
#ifdef oneD
      call getInput('grid', 'tileX', species(s) % tile_sx)
      species(s) % tile_sy = 1
      species(s) % tile_sz = 1
#elif twoD
      call getInput('grid', 'tileX', species(s) % tile_sx)
      call getInput('grid', 'tileY', species(s) % tile_sy)
      species(s) % tile_sz = 1
#elif threeD
      call getInput('grid', 'tileX', species(s) % tile_sx)
      call getInput('grid', 'tileY', species(s) % tile_sy)
      call getInput('grid', 'tileZ', species(s) % tile_sz)
#endif
      species(s) % tile_nx = ceiling(real(this_meshblock % ptr % sx) / real(species(s) % tile_sx))
      species(s) % tile_ny = ceiling(real(this_meshblock % ptr % sy) / real(species(s) % tile_sy))
      species(s) % tile_nz = ceiling(real(this_meshblock % ptr % sz) / real(species(s) % tile_sz))
      allocate (species(s) % prtl_tile(species(s) % tile_nx, &
                                       species(s) % tile_ny, &
                                       species(s) % tile_nz))
    end do

    allocate (maxptl_array(nspec))

    do s = 1, nspec
      write (var_name, "(A6,I1)") "maxptl", s
      call getInput('particles', var_name, maxptl_array(s))
      write (var_name, "(A1,I1)") "m", s
      call getInput('particles', var_name, species(s) % m_sp)
      write (var_name, "(A2,I1)") "ch", s
      call getInput('particles', var_name, species(s) % ch_sp)

      write (var_name, "(A7,I1)") "deposit", s
      call getInput('particles', var_name, species(s) % deposit_sp, (species(s) % ch_sp .ne. 0))
      write (var_name, "(A4,I1)") "move", s
      call getInput('particles', var_name, species(s) % move_sp, .true.)
      write (var_name, "(A6,I1)") "output", s
      call getInput('particles', var_name, species(s) % output_sp, .true.)

      if ((species(s) % m_sp .eq. 0) .and. (species(s) % ch_sp .ne. 0)) then
        call throwError('ERROR: massless charged particles are not allowed')
      end if
      if ((species(s) % m_sp .ne. 0) .and. (species(s) % ch_sp .eq. 0)) then
        call throwError('ERROR: massive zero-charge particles are not allowed')
      end if
      if ((species(s) % ch_sp .eq. 0) .and. (species(s) % deposit_sp)) then
        call throwError('ERROR: zero-charged particles cannot deposit current')
      end if

#ifdef DOWNSAMPLING
      write (var_name, "(A3,I1)") "dwn", s
      call getInput('particles', var_name, species(s) % dwn_sp, .false.)
#endif

      ! extra physics properties
#ifdef RADIATION
      write (var_name, "(A4,I1)") "cool", s
      call getInput('particles', var_name, species(s) % cool_sp, .false.)
      if ((species(s) % cool_sp) .and. (species(s) % m_sp .eq. 0)) then
        call throwError('Unable to cool `m=0` particles.')
      end if
#endif

      do ti = 1, species(s) % tile_nx
        do tj = 1, species(s) % tile_ny
          do tk = 1, species(s) % tile_nz
            call createEmptyTile(s, ti, tj, tk, maxptl_array(s))
          end do
        end do
      end do
      species(s) % cntr_sp = 0
    end do
    call printDiag("initializeParticles()", 1)
  end subroutine initializeParticles

  ! Subroutine to move particles around
  subroutine createParticleFromAttributes(s, xi, yi, zi, dx, dy, dz, &
                                          u, v, w, &
#ifdef PRTLPAYLOADS
                                          payload1, payload2, payload3, &
#endif
                                          ind, proc, weight)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in) :: s
    integer(kind=2), intent(in) :: xi, yi, zi
    real, intent(in) :: dx, dy, dz, u, v, w

#ifdef PRTLPAYLOADS
    real, intent(in) :: payload1, payload2, payload3
#endif

    integer, intent(in) :: ind, proc
    real, intent(in) :: weight

    integer :: ti, tj, tk

    character(len=STR_MAX) :: dummy_string

    ti = 1; tj = 1; tk = 1
#if defined(oneD) || defined (twoD) || defined (threeD)
    ti = FLOOR(REAL(xi) / REAL(species(s) % tile_sx)) + 1
#endif
#if defined (twoD) || defined (threeD)
    tj = FLOOR(REAL(yi) / REAL(species(s) % tile_sy)) + 1
#endif
#if defined (threeD)
    tk = FLOOR(REAL(zi) / REAL(species(s) % tile_sz)) + 1
#endif
    ! if the debug flag is enabled ...
    ! ... check that the particle is within the boundaries ...
    ! ... of that tile and that the tile exists
#ifdef DEBUG
    if ((s .le. 0) .or. (s .gt. nspec)) then
      call throwError('Wrong species in `createParticleFromAttributes`.')
    end if
    if ((ti .gt. species(s) % tile_nx) .or. &
        (tj .gt. species(s) % tile_ny) .or. &
        (tk .gt. species(s) % tile_nz)) then
      print *, mpi_rank, xi, yi, zi, ti, tj, tk
      print *, species(s) % tile_nx, species(s) % tile_ny, species(s) % tile_nz
      call throwError('ERROR: wrong ti, tj, tk in `createParticleFromAttributes`')
    end if
    if ((xi .lt. species(s) % prtl_tile(ti, tj, tk) % x1) .or. &
        (xi .ge. species(s) % prtl_tile(ti, tj, tk) % x2) .or. &
        (yi .lt. species(s) % prtl_tile(ti, tj, tk) % y1) .or. &
        (yi .ge. species(s) % prtl_tile(ti, tj, tk) % y2) .or. &
        (zi .lt. species(s) % prtl_tile(ti, tj, tk) % z1) .or. &
        (zi .ge. species(s) % prtl_tile(ti, tj, tk) % z2)) then
      print *, s, xi, yi, zi, dx, dy, dz
      print *, species(s) % prtl_tile(ti, tj, tk) % x1, &
        species(s) % prtl_tile(ti, tj, tk) % x2, &
        species(s) % prtl_tile(ti, tj, tk) % y1, &
        species(s) % prtl_tile(ti, tj, tk) % y2, &
        species(s) % prtl_tile(ti, tj, tk) % z1, &
        species(s) % prtl_tile(ti, tj, tk) % z2
      print *, ti, tj, tk
      print *, species(s) % tile_nx, species(s) % tile_ny, species(s) % tile_nz
      print *, species(s) % tile_sx, species(s) % tile_sy, species(s) % tile_sz
      call throwError('ERROR: wrong ti, tj, tk in `createParticleFromAttributes` according to x1,x2,etc')
    end if
#endif
    if (species(s) % prtl_tile(ti, tj, tk) % npart_sp .eq. species(s) % prtl_tile(ti, tj, tk) % maxptl_sp) then
      write (dummy_string, '(I5)') s
      if (resize_tiles) then
        call reallocTileSize(species(s) % prtl_tile(ti, tj, tk), .true.)
      else
        call throwError('ERROR: npart_sp > maxptl_sp in createParticleFromAttributes for species #'//trim(dummy_string))
      end if
    end if

    call putParticleOnTile(s, ti, tj, tk, &
                           xi, yi, zi, dx, dy, dz, &
                           u, v, w, &
#ifdef PRTLPAYLOADS
                           payload1, payload2, payload3, &
#endif
                           ind, proc, weight)
  end subroutine createParticleFromAttributes

  subroutine putParticleOnTile(s, ti, tj, tk, &
                               xi, yi, zi, dx, dy, dz, &
                               u, v, w, &
#ifdef PRTLPAYLOADS
                               payload1, payload2, payload3, &
#endif
                               ind, proc, weight)
    implicit none
    integer, intent(in) :: s
    integer, intent(in) :: ti, tj, tk
    integer(kind=2), intent(in) :: xi, yi, zi
    real, intent(in) :: dx, dy, dz, u, v, w

#ifdef PRTLPAYLOADS
    real, intent(in) :: payload1, payload2, payload3
#endif

    integer, intent(in) :: ind, proc
    real, intent(in) :: weight

    integer :: p

    species(s) % prtl_tile(ti, tj, tk) % npart_sp = species(s) % prtl_tile(ti, tj, tk) % npart_sp + 1
    p = species(s) % prtl_tile(ti, tj, tk) % npart_sp

    species(s) % prtl_tile(ti, tj, tk) % xi(p) = xi
    species(s) % prtl_tile(ti, tj, tk) % dx(p) = dx

    species(s) % prtl_tile(ti, tj, tk) % yi(p) = yi
    species(s) % prtl_tile(ti, tj, tk) % dy(p) = dy

    species(s) % prtl_tile(ti, tj, tk) % zi(p) = zi
    species(s) % prtl_tile(ti, tj, tk) % dz(p) = dz

    species(s) % prtl_tile(ti, tj, tk) % u(p) = u
    species(s) % prtl_tile(ti, tj, tk) % v(p) = v
    species(s) % prtl_tile(ti, tj, tk) % w(p) = w

    species(s) % prtl_tile(ti, tj, tk) % ind(p) = ind
    species(s) % prtl_tile(ti, tj, tk) % proc(p) = proc

    species(s) % prtl_tile(ti, tj, tk) % weight(p) = weight

#ifdef PRTLPAYLOADS
    species(s) % prtl_tile(ti, tj, tk) % payload1(p) = payload1
    species(s) % prtl_tile(ti, tj, tk) % payload2(p) = payload2
    species(s) % prtl_tile(ti, tj, tk) % payload3(p) = payload3
#endif
  end subroutine putParticleOnTile

  subroutine putEnrouteParticleOnTile(s, ti, tj, tk, enroute)
    implicit none
    integer, intent(in) :: s, ti, tj, tk
    type(prtl_enroute), intent(in) :: enroute
    call putParticleOnTile(s, ti, tj, tk, &
                           enroute % xi, enroute % yi, enroute % zi, enroute % dx, enroute % dy, enroute % dz, &
                           enroute % u, enroute % v, enroute % w, &
#ifdef PRTLPAYLOADS
                           enroute % payload1, enroute % payload2, enroute % payload3, &
#endif
                           enroute % ind, enroute % proc, enroute % weight)
  end subroutine putEnrouteParticleOnTile

  subroutine copyParticleFromTo(s, p_from, p_to, ti, tj, tk)
    ! DEP_PRT [particle-dependent]
    implicit none
    ! within a single tile
    integer, intent(in) :: s, p_from, p_to, ti, tj, tk
    species(s) % prtl_tile(ti, tj, tk) % xi(p_to) = species(s) % prtl_tile(ti, tj, tk) % xi(p_from)
    species(s) % prtl_tile(ti, tj, tk) % yi(p_to) = species(s) % prtl_tile(ti, tj, tk) % yi(p_from)
    species(s) % prtl_tile(ti, tj, tk) % zi(p_to) = species(s) % prtl_tile(ti, tj, tk) % zi(p_from)

    species(s) % prtl_tile(ti, tj, tk) % dx(p_to) = species(s) % prtl_tile(ti, tj, tk) % dx(p_from)
    species(s) % prtl_tile(ti, tj, tk) % dy(p_to) = species(s) % prtl_tile(ti, tj, tk) % dy(p_from)
    species(s) % prtl_tile(ti, tj, tk) % dz(p_to) = species(s) % prtl_tile(ti, tj, tk) % dz(p_from)

    species(s) % prtl_tile(ti, tj, tk) % u(p_to) = species(s) % prtl_tile(ti, tj, tk) % u(p_from)
    species(s) % prtl_tile(ti, tj, tk) % v(p_to) = species(s) % prtl_tile(ti, tj, tk) % v(p_from)
    species(s) % prtl_tile(ti, tj, tk) % w(p_to) = species(s) % prtl_tile(ti, tj, tk) % w(p_from)

    species(s) % prtl_tile(ti, tj, tk) % ind(p_to) = species(s) % prtl_tile(ti, tj, tk) % ind(p_from)
    species(s) % prtl_tile(ti, tj, tk) % proc(p_to) = species(s) % prtl_tile(ti, tj, tk) % proc(p_from)

    species(s) % prtl_tile(ti, tj, tk) % weight(p_to) = species(s) % prtl_tile(ti, tj, tk) % weight(p_from)

#ifdef PRTLPAYLOADS
    species(s) % prtl_tile(ti, tj, tk) % payload1(p_to) = species(s) % prtl_tile(ti, tj, tk) % payload1(p_from)
    species(s) % prtl_tile(ti, tj, tk) % payload2(p_to) = species(s) % prtl_tile(ti, tj, tk) % payload2(p_from)
    species(s) % prtl_tile(ti, tj, tk) % payload3(p_to) = species(s) % prtl_tile(ti, tj, tk) % payload3(p_from)
#endif
  end subroutine copyParticleFromTo

  subroutine createEmptyTile(s, ti, tj, tk, maxptl, meshblock)
    implicit none
    integer(kind=8), intent(in) :: maxptl
    integer, intent(in) :: s, ti, tj, tk
    type(mesh), optional, intent(in) :: meshblock
    integer :: maxptl_on_tile
    type(mesh) :: meshblock_
    if (present(meshblock)) then
      meshblock_ = meshblock
    else
      meshblock_ = this_meshblock % ptr
    end if

    maxptl_on_tile = INT(maxptl / INT(species(s) % tile_nx * species(s) % tile_ny * species(s) % tile_nz, 8), 4)

    species(s) % prtl_tile(ti, tj, tk) % spec = s

    species(s) % prtl_tile(ti, tj, tk) % x1 = (ti - 1) * species(s) % tile_sx
    species(s) % prtl_tile(ti, tj, tk) % x2 = min(ti * species(s) % tile_sx, meshblock_ % sx)
    species(s) % prtl_tile(ti, tj, tk) % y1 = (tj - 1) * species(s) % tile_sy
    species(s) % prtl_tile(ti, tj, tk) % y2 = min(tj * species(s) % tile_sy, meshblock_ % sy)
    species(s) % prtl_tile(ti, tj, tk) % z1 = (tk - 1) * species(s) % tile_sz
    species(s) % prtl_tile(ti, tj, tk) % z2 = min(tk * species(s) % tile_sz, meshblock_ % sz)
#ifdef DEBUG
    if ((species(s) % prtl_tile(ti, tj, tk) % x1 .eq. 0) .and. &
        (species(s) % prtl_tile(ti, tj, tk) % x2 .eq. 0) .and. &
        (species(s) % prtl_tile(ti, tj, tk) % y1 .eq. 0) .and. &
        (species(s) % prtl_tile(ti, tj, tk) % y2 .eq. 0) .and. &
        (species(s) % prtl_tile(ti, tj, tk) % z1 .eq. 0) .and. &
        (species(s) % prtl_tile(ti, tj, tk) % z2 .eq. 0)) then
      print *, ti, tj, tk
      print *, species(s) % prtl_tile(ti, tj, tk) % x1, &
        species(s) % prtl_tile(ti, tj, tk) % x2, &
        species(s) % prtl_tile(ti, tj, tk) % y1, &
        species(s) % prtl_tile(ti, tj, tk) % y2, &
        species(s) % prtl_tile(ti, tj, tk) % z1, &
        species(s) % prtl_tile(ti, tj, tk) % z2
      call throwError('ERROR: in `initializeParticles`')
    end if
#endif
    call allocateParticlesOnEmptyTile(species(s) % prtl_tile(ti, tj, tk), maxptl_on_tile)
  end subroutine createEmptyTile

  subroutine allocateParticlesOnEmptyTile(tile, sz)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in) :: sz
    type(particle_tile), intent(inout) :: tile

    tile % npart_sp = 0
    tile % maxptl_sp = sz

    if (allocated(tile % xi)) deallocate (tile % xi)
    if (allocated(tile % yi)) deallocate (tile % yi)
    if (allocated(tile % zi)) deallocate (tile % zi)
    if (allocated(tile % dx)) deallocate (tile % dx)
    if (allocated(tile % dy)) deallocate (tile % dy)
    if (allocated(tile % dz)) deallocate (tile % dz)
    if (allocated(tile % u)) deallocate (tile % u)
    if (allocated(tile % v)) deallocate (tile % v)
    if (allocated(tile % w)) deallocate (tile % w)
    if (allocated(tile % ind)) deallocate (tile % ind)
    if (allocated(tile % proc)) deallocate (tile % proc)
    if (allocated(tile % weight)) deallocate (tile % weight)
    allocate (tile % xi(sz)); allocate (tile % yi(sz)); allocate (tile % zi(sz))
    allocate (tile % dx(sz)); allocate (tile % dy(sz)); allocate (tile % dz(sz))
    allocate (tile % u(sz)); allocate (tile % v(sz)); allocate (tile % w(sz))
    allocate (tile % ind(sz)); allocate (tile % proc(sz)); allocate (tile % weight(sz))

#ifdef PRTLPAYLOADS
    if (allocated(tile % payload1)) deallocate (tile % payload1)
    if (allocated(tile % payload2)) deallocate (tile % payload2)
    if (allocated(tile % payload3)) deallocate (tile % payload3)
    allocate (tile % payload1(sz))
    allocate (tile % payload2(sz))
    allocate (tile % payload3(sz))
#endif
  end subroutine allocateParticlesOnEmptyTile

  subroutine checkTileSizes()
    implicit none
    integer :: s, ti, tj, tk
    do s = 1, nspec ! loop over species
      do ti = 1, species(s) % tile_nx
        do tj = 1, species(s) % tile_ny
          do tk = 1, species(s) % tile_nz
            if (species(s) % prtl_tile(ti, tj, tk) % npart_sp .ge. species(s) % prtl_tile(ti, tj, tk) % maxptl_sp * 0.7) then
              ! increase the tile size
              if (resize_tiles) then
                call reallocTileSize(species(s) % prtl_tile(ti, tj, tk), .true.)
              else
                print *, "DANGER: `maxptl_sp` in tiles too low, consider increasing `maxptl` or turning on `resize_tiles`."
              end if
            else if ((species(s) % prtl_tile(ti, tj, tk) % npart_sp .lt. (species(s) % prtl_tile(ti, tj, tk) % maxptl_sp * 0.3)) .and. &
                     ((species(s) % prtl_tile(ti, tj, tk) % maxptl_sp * 0.5) .gt. min_tile_nprt)) then
              ! decrease the tile size
              if (resize_tiles) then
                call reallocTileSize(species(s) % prtl_tile(ti, tj, tk), .false.)
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
    type(particle_tile), intent(inout) :: tile
    ! `.true.` if need to increase, otherwise `.false.`
    logical, intent(in) :: increase_flag
    integer(kind=2), allocatable, dimension(:) :: dummy_int2
    integer, allocatable, dimension(:) :: dummy_int
    real, allocatable, dimension(:) :: dummy_real
    integer :: current_npart
    if (increase_flag) then
      ! increase twice
      tile % maxptl_sp = INT(tile % maxptl_sp * 1.5)
    else
      ! decrease twice
      tile % maxptl_sp = INT(tile % maxptl_sp * 0.5)
    end if

#ifndef LOWMEM
    if (tile % npart_sp .gt. tile % maxptl_sp) then
      call throwError('ERROR: `npart > maxptl` in `reallocTileSize`')
    end if
#else
    do while (tile % npart_sp .gt. tile % maxptl_sp)
      tile % maxptl_sp = INT(tile % maxptl_sp * 1.5) + 1
#ifdef DEBUG
      print *, 'increasing', tile % npart_sp, tile % maxptl_sp
#endif
    end do
#endif

    allocate (dummy_int2(tile % maxptl_sp))
    allocate (dummy_int(tile % maxptl_sp))
    allocate (dummy_real(tile % maxptl_sp))

    current_npart = tile % npart_sp

    dummy_int2(1:current_npart) = tile % xi(1:current_npart)
    deallocate (tile % xi); allocate (tile % xi(tile % maxptl_sp))
    tile % xi(1:current_npart) = dummy_int2(1:current_npart)

    dummy_int2(1:current_npart) = tile % yi(1:current_npart)
    deallocate (tile % yi); allocate (tile % yi(tile % maxptl_sp))
    tile % yi(1:current_npart) = dummy_int2(1:current_npart)

    dummy_int2(1:current_npart) = tile % zi(1:current_npart)
    deallocate (tile % zi); allocate (tile % zi(tile % maxptl_sp))
    tile % zi(1:current_npart) = dummy_int2(1:current_npart)

    dummy_real(1:current_npart) = tile % weight(1:current_npart)
    deallocate (tile % weight); allocate (tile % weight(tile % maxptl_sp))
    tile % weight(1:current_npart) = dummy_real(1:current_npart)

    dummy_real(1:current_npart) = tile % dx(1:current_npart)
    deallocate (tile % dx); allocate (tile % dx(tile % maxptl_sp))
    tile % dx(1:current_npart) = dummy_real(1:current_npart)

    dummy_real(1:current_npart) = tile % dy(1:current_npart)
    deallocate (tile % dy); allocate (tile % dy(tile % maxptl_sp))
    tile % dy(1:current_npart) = dummy_real(1:current_npart)

    dummy_real(1:current_npart) = tile % dz(1:current_npart)
    deallocate (tile % dz); allocate (tile % dz(tile % maxptl_sp))
    tile % dz(1:current_npart) = dummy_real(1:current_npart)

    dummy_real(1:current_npart) = tile % u(1:current_npart)
    deallocate (tile % u); allocate (tile % u(tile % maxptl_sp))
    tile % u(1:current_npart) = dummy_real(1:current_npart)

    dummy_real(1:current_npart) = tile % v(1:current_npart)
    deallocate (tile % v); allocate (tile % v(tile % maxptl_sp))
    tile % v(1:current_npart) = dummy_real(1:current_npart)

    dummy_real(1:current_npart) = tile % w(1:current_npart)
    deallocate (tile % w); allocate (tile % w(tile % maxptl_sp))
    tile % w(1:current_npart) = dummy_real(1:current_npart)

    dummy_int(1:current_npart) = tile % ind(1:current_npart)
    deallocate (tile % ind); allocate (tile % ind(tile % maxptl_sp))
    tile % ind(1:current_npart) = dummy_int(1:current_npart)

    dummy_int(1:current_npart) = tile % proc(1:current_npart)
    deallocate (tile % proc); allocate (tile % proc(tile % maxptl_sp))
    tile % proc(1:current_npart) = dummy_int(1:current_npart)

#ifdef PRTLPAYLOADS
    dummy_real(1:current_npart) = tile % payload1(1:current_npart)
    deallocate (tile % payload1); allocate (tile % payload1(tile % maxptl_sp))
    tile % payload1(1:current_npart) = dummy_real(1:current_npart)

    dummy_real(1:current_npart) = tile % payload2(1:current_npart)
    deallocate (tile % payload2); allocate (tile % payload2(tile % maxptl_sp))
    tile % payload2(1:current_npart) = dummy_real(1:current_npart)

    dummy_real(1:current_npart) = tile % payload3(1:current_npart)
    deallocate (tile % payload3); allocate (tile % payload3(tile % maxptl_sp))
    tile % payload3(1:current_npart) = dummy_real(1:current_npart)
#endif

    deallocate (dummy_int2)
    deallocate (dummy_int)
    deallocate (dummy_real)
  end subroutine reallocTileSize

  ! Subroutine to create brand new particles
  subroutine createParticle(s, xi, yi, zi, dx, dy, dz, u, v, w, &
                            ind, proc, weight)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in) :: s
    integer(kind=2), intent(in) :: xi, yi, zi
    real, intent(in) :: dx, dy, dz, u, v, w
    integer, optional, intent(in) :: ind, proc
    real, optional :: weight
    integer :: ind_, proc_
    real :: weight_
#ifdef PRTLPAYLOADS
    real :: payload1, payload2, payload3
    payload1 = 0.0; payload2 = 0.0; payload3 = 0.0
#endif
    if (present(ind) .and. present(proc)) then
      ! moving particle from one tile/meshblock to another
      ind_ = ind
      proc_ = proc
    else
      ! create a very new particle
      ind_ = species(s) % cntr_sp
      proc_ = mpi_rank
      species(s) % cntr_sp = species(s) % cntr_sp + 1
    end if

    if (present(weight)) then
      weight_ = weight
    else
      weight_ = 1
    end if
    call createParticleFromAttributes(s, xi=xi, yi=yi, zi=zi, dx=dx, dy=dy, dz=dz, &
                                      u=u, v=v, w=w, &
#ifdef PRTLPAYLOADS
                                      payload1=payload1, payload2=payload2, payload3=payload3, &
#endif
                                      ind=ind_, proc=proc_, weight=weight_)
  end subroutine createParticle

  subroutine injectParticleGlobally(s, x_glob, y_glob, z_glob, u, v, w, weight)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in) :: s
    real, intent(in) :: x_glob, y_glob, z_glob
    real, intent(in) :: u, v, w
    real :: x_g, y_g, z_g, x_loc, y_loc, z_loc
    real :: dx_, dy_, dz_
    integer(kind=2) :: xi_, yi_, zi_
    real :: weight_
    real, optional, intent(in) :: weight
    logical :: contained_flag

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

  subroutine injectParticleLocally(s, x_loc, y_loc, z_loc, &
                                   u, v, w, ind, proc, weight)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in) :: s
    real, intent(in) :: x_loc, y_loc, z_loc, u, v, w
    integer, optional, intent(in) :: ind, proc
    real, optional, intent(in) :: weight
    integer(kind=2) :: xi, yi, zi
    real :: dx, dy, dz

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
    integer :: s, p, ti, tj, tk
    do s = 1, nspec ! loop over species
      do ti = 1, species(s) % tile_nx
        do tj = 1, species(s) % tile_ny
          do tk = 1, species(s) % tile_nz
            ! FIX1 try to vectorize this
            p = 1
            do while (p .le. species(s) % prtl_tile(ti, tj, tk) % npart_sp)
              if (species(s) % prtl_tile(ti, tj, tk) % proc(p) .lt. 0) then
                call removeParticleFromTile(s, ti, tj, tk, p)
              else
                p = p + 1
              end if
            end do ! p
          end do ! tk
        end do ! tj
      end do ! ti
    end do ! s
    call printDiag("clearGhostParticles()", 2)
  end subroutine clearGhostParticles

  subroutine removeParticleFromTile(s, ti, tj, tk, p)
    implicit none
    integer, intent(in) :: s, ti, tj, tk, p
    if (p .ne. species(s) % prtl_tile(ti, tj, tk) % npart_sp) then
      call copyParticleFromTo(s, species(s) % prtl_tile(ti, tj, tk) % npart_sp, p, ti, tj, tk)
    end if
    species(s) % prtl_tile(ti, tj, tk) % npart_sp = species(s) % prtl_tile(ti, tj, tk) % npart_sp - 1
  end subroutine removeParticleFromTile

  subroutine copyToEnroute(spec_id, ti, tj, tk, prtl_id, enroute)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in) :: spec_id, prtl_id, ti, tj, tk
    type(prtl_enroute), intent(inout) :: enroute
    enroute % weight = species(spec_id) % prtl_tile(ti, tj, tk) % weight(prtl_id)
    enroute % xi = species(spec_id) % prtl_tile(ti, tj, tk) % xi(prtl_id)
    enroute % yi = species(spec_id) % prtl_tile(ti, tj, tk) % yi(prtl_id)
    enroute % zi = species(spec_id) % prtl_tile(ti, tj, tk) % zi(prtl_id)
    enroute % dx = species(spec_id) % prtl_tile(ti, tj, tk) % dx(prtl_id)
    enroute % dy = species(spec_id) % prtl_tile(ti, tj, tk) % dy(prtl_id)
    enroute % dz = species(spec_id) % prtl_tile(ti, tj, tk) % dz(prtl_id)
    enroute % u = species(spec_id) % prtl_tile(ti, tj, tk) % u(prtl_id)
    enroute % v = species(spec_id) % prtl_tile(ti, tj, tk) % v(prtl_id)
    enroute % w = species(spec_id) % prtl_tile(ti, tj, tk) % w(prtl_id)
    enroute % ind = species(spec_id) % prtl_tile(ti, tj, tk) % ind(prtl_id)
    enroute % proc = species(spec_id) % prtl_tile(ti, tj, tk) % proc(prtl_id)

#ifdef PRTLPAYLOADS
    enroute % payload1 = species(spec_id) % prtl_tile(ti, tj, tk) % payload1(prtl_id)
    enroute % payload2 = species(spec_id) % prtl_tile(ti, tj, tk) % payload2(prtl_id)
    enroute % payload3 = species(spec_id) % prtl_tile(ti, tj, tk) % payload3(prtl_id)
#endif
  end subroutine copyToEnroute

  subroutine copyFromEnroute(enroute, spec_id)
    implicit none
    type(prtl_enroute), intent(in) :: enroute
    integer, intent(in) :: spec_id
    ! DEP_PRT [particle-dependent]
    call createParticleFromAttributes(spec_id, enroute % xi, enroute % yi, enroute % zi, &
                                      enroute % dx, enroute % dy, enroute % dz, &
                                      enroute % u, enroute % v, enroute % w, &
#ifdef PRTLPAYLOADS
                                      enroute % payload1, enroute % payload2, enroute % payload3, &
#endif
                                      enroute % ind, enroute % proc, enroute % weight)
  end subroutine copyFromEnroute

  subroutine shiftParticlesX(shift)
    implicit none
    integer, intent(in) :: shift
    integer :: s, ti, tj, tk, p
    do s = 1, nspec
      do ti = 1, species(s) % tile_nx
        do tj = 1, species(s) % tile_ny
          do tk = 1, species(s) % tile_nz
            do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
              species(s) % prtl_tile(ti, tj, tk) % xi(p) = species(s) % prtl_tile(ti, tj, tk) % xi(p) + INT(shift, 2)
            end do
          end do
        end do
      end do
    end do
  end subroutine shiftParticlesX

  subroutine shiftParticlesY(shift)
    implicit none
    integer, intent(in) :: shift
    integer :: s, ti, tj, tk, p
    do s = 1, nspec
      do ti = 1, species(s) % tile_nx
        do tj = 1, species(s) % tile_ny
          do tk = 1, species(s) % tile_nz
            do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
              species(s) % prtl_tile(ti, tj, tk) % yi(p) = species(s) % prtl_tile(ti, tj, tk) % yi(p) + INT(shift, 2)
            end do
          end do
        end do
      end do
    end do
  end subroutine shiftParticlesY

  subroutine shiftParticlesZ(shift)
    implicit none
    integer, intent(in) :: shift
    integer :: s, ti, tj, tk, p
    do s = 1, nspec
      do ti = 1, species(s) % tile_nx
        do tj = 1, species(s) % tile_ny
          do tk = 1, species(s) % tile_nz
            do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
              species(s) % prtl_tile(ti, tj, tk) % zi(p) = species(s) % prtl_tile(ti, tj, tk) % zi(p) + INT(shift, 2)
            end do
          end do
        end do
      end do
    end do
  end subroutine shiftParticlesZ

  subroutine extractParticlesFromEnroute(cnt, spec_id)
    implicit none
    integer, intent(in) :: cnt, spec_id
    integer :: p
    do p = 1, cnt
      call copyFromEnroute(recv_enroute % enroute(p), spec_id)
    end do
  end subroutine extractParticlesFromEnroute

  subroutine reallocateParticles(meshblock)
    implicit none
    type(mesh), intent(in) :: meshblock
    integer :: s, p, ti, tj, tk, ti_, tj_, tk_
    integer(kind=8) :: maxptl_
    integer, dimension(:, :, :), allocatable :: maxptl_tmpar

    call deallocateParticles()

    do s = 1, nspec
      species(s) % tile_nx = ceiling(real(meshblock % sx) / real(species(s) % tile_sx))
      species(s) % tile_ny = ceiling(real(meshblock % sy) / real(species(s) % tile_sy))
      species(s) % tile_nz = ceiling(real(meshblock % sz) / real(species(s) % tile_sz))
      allocate (species(s) % prtl_tile(species(s) % tile_nx, &
                                       species(s) % tile_ny, &
                                       species(s) % tile_nz))

      if (resize_tiles) then

        allocate (maxptl_tmpar(1:species(s) % tile_nx, 1:species(s) % tile_ny, 1:species(s) % tile_nz))

        maxptl_tmpar = 0

        do p = 1, prtl_backup(s) % cnt
          ti_ = 1; tj_ = 1; tk_ = 1
#if defined(oneD) || defined (twoD) || defined (threeD)
          ti_ = FLOOR(REAL(prtl_backup(s) % enroute(p) % xi) / REAL(species(s) % tile_sx)) + 1
#endif
#if defined (twoD) || defined (threeD)
          tj_ = FLOOR(REAL(prtl_backup(s) % enroute(p) % yi) / REAL(species(s) % tile_sy)) + 1
#endif
#if defined (threeD)
          tk_ = FLOOR(REAL(prtl_backup(s) % enroute(p) % zi) / REAL(species(s) % tile_sz)) + 1
#endif
          ti_ = MIN(MAX(ti_, 1), species(s) % tile_nx)
          tj_ = MIN(MAX(tj_, 1), species(s) % tile_ny)
          tk_ = MIN(MAX(tk_, 1), species(s) % tile_nz)
          maxptl_tmpar(ti_, tj_, tk_) = maxptl_tmpar(ti_, tj_, tk_) + 1
        end do

        do ti = 1, species(s) % tile_nx
          do tj = 1, species(s) % tile_ny
            do tk = 1, species(s) % tile_nz

              maxptl_ = max(maxptl_tmpar(ti, tj, tk), min_tile_nprt)
              maxptl_ = maxptl_ * (species(s) % tile_nx * species(s) % tile_ny * species(s) % tile_nz)

              call createEmptyTile(s, ti, tj, tk, maxptl_, meshblock)

            end do
          end do
        end do

        deallocate (maxptl_tmpar)

      else

        do ti = 1, species(s) % tile_nx
          do tj = 1, species(s) % tile_ny
            do tk = 1, species(s) % tile_nz

              call createEmptyTile(s, ti, tj, tk, maxptl_array(s), meshblock)

            end do
          end do
        end do

      end if

      species(s) % cntr_sp = 0
    end do

    call reallocateEnrouteArray(meshblock)
  end subroutine reallocateParticles

  subroutine reallocateEnrouteArray(meshblock)
    implicit none
    type(mesh), intent(in) :: meshblock
    integer :: buffsize, buffsize_x, buffsize_y
    integer :: buffsize_xy
    integer :: buffsize_z, buffsize_xz, buffsize_yz
    integer :: buffsize_xyz, old_buffsize, min_buffsize
    integer :: multiplier, ind1, ind2, ind3
    type(prtl_enroute), allocatable :: enroute_temp(:)

    multiplier = max(INT(ppc0), 1) * max_buffsize

    buffsize_x = 0
    buffsize_y = 0; buffsize_xy = 0
    buffsize_z = 0; buffsize_xz = 0; buffsize_yz = 0; buffsize_xyz = 0
#if defined (oneD) || defined (twoD) || defined (threeD)
    buffsize_x = meshblock % sy * meshblock % sz * multiplier
#endif
#if defined (twoD) || defined (threeD)
    buffsize_y = meshblock % sx * meshblock % sz * multiplier
    buffsize_xy = meshblock % sz * multiplier
#endif
#if defined(threeD)
    buffsize_z = meshblock % sx * meshblock % sy * multiplier
    buffsize_xz = meshblock % sz * multiplier
    buffsize_yz = meshblock % sx * multiplier
    buffsize_xyz = multiplier
#endif

#ifdef oneD
    buffsize = multiplier
#elif twoD
    buffsize = MAX0(meshblock % sx, meshblock % sy, meshblock % sz) * multiplier
#elif threeD
    buffsize = MAX0(meshblock % sx, meshblock % sy, meshblock % sz)**2 * multiplier
#endif

    if (allocated(recv_enroute % enroute)) then
      old_buffsize = recv_enroute % max
      min_buffsize = MIN(old_buffsize, buffsize)

      ! reallocate
      allocate (enroute_temp(1:buffsize))
      enroute_temp(1:min_buffsize) = recv_enroute % enroute(1:min_buffsize)
      deallocate (recv_enroute % enroute)
      allocate (recv_enroute % enroute(1:buffsize))
      recv_enroute % enroute(1:min_buffsize) = enroute_temp(1:min_buffsize)
      deallocate (enroute_temp)
    else
      ! allocate from scratch
      allocate (recv_enroute % enroute(1:buffsize))
    end if

    recv_enroute % max = buffsize
    recv_enroute % cnt = 0

#ifndef LOWMEM

    do ind1 = -1, 1
      do ind2 = -1, 1
        do ind3 = -1, 1
          if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
#ifdef oneD
          if ((ind2 .ne. 0) .or. (ind3 .ne. 0)) cycle
#elif twoD
          if (ind3 .ne. 0) cycle
#endif
          if ((ind2 .eq. 0) .and. (ind3 .eq. 0)) then
            buffsize = buffsize_x
          else if ((ind1 .eq. 0) .and. (ind3 .eq. 0)) then
            buffsize = buffsize_y
          else if ((ind1 .eq. 0) .and. (ind2 .eq. 0)) then
            buffsize = buffsize_z
          else if (ind3 .eq. 0) then
            buffsize = buffsize_xy
          else if (ind2 .eq. 0) then
            buffsize = buffsize_xz
          else if (ind1 .eq. 0) then
            buffsize = buffsize_yz
          else
            buffsize = buffsize_xyz
          end if
          if (allocated(enroute_bot % get(ind1, ind2, ind3) % enroute)) then
            ! reallocate
            old_buffsize = enroute_bot % get(ind1, ind2, ind3) % max
            min_buffsize = MIN(old_buffsize, buffsize)
            allocate (enroute_temp(1:buffsize))
            enroute_temp(1:min_buffsize) = enroute_bot % get(ind1, ind2, ind3) % enroute(1:min_buffsize)

            call reallocateEnroute(ind1, ind2, ind3, buffsize)

            enroute_bot % get(ind1, ind2, ind3) % enroute(1:min_buffsize) = enroute_temp(1:min_buffsize)
            deallocate (enroute_temp)
          else
            ! allocate from scratch
            call reallocateEnroute(ind1, ind2, ind3, buffsize)

          end if
          enroute_bot % get(ind1, ind2, ind3) % max = buffsize
          enroute_bot % get(ind1, ind2, ind3) % cnt = 0
        end do
      end do
    end do

#endif

  end subroutine reallocateEnrouteArray

  subroutine reallocateEnroute(ind1, ind2, ind3, buffsize)
    implicit none
    integer, intent(in) :: ind1, ind2, ind3, buffsize
    enroute_bot % get(ind1, ind2, ind3) % max = buffsize
    if (allocated(enroute_bot % get(ind1, ind2, ind3) % enroute)) then
      deallocate (enroute_bot % get(ind1, ind2, ind3) % enroute)
    end if
    allocate (enroute_bot % get(ind1, ind2, ind3) % enroute(buffsize))
  end subroutine reallocateEnroute

  subroutine backupParticles()
    implicit none
    integer :: s, p, ti, tj, tk
    integer :: npart
    call deallocateParticleBackup()
    allocate (prtl_backup(nspec))

    do s = 1, nspec
      ! count number of particles
      npart = 0
      do ti = 1, species(s) % tile_nx
        do tj = 1, species(s) % tile_ny
          do tk = 1, species(s) % tile_nz
            npart = npart + species(s) % prtl_tile(ti, tj, tk) % npart_sp
          end do
        end do
      end do

      ! allocate backup array
      prtl_backup(s) % max = npart
      allocate (prtl_backup(s) % enroute(npart))

      ! copy particles to backup array
      prtl_backup(s) % cnt = 0
      do ti = 1, species(s) % tile_nx
        do tj = 1, species(s) % tile_ny
          do tk = 1, species(s) % tile_nz
            do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
              prtl_backup(s) % cnt = prtl_backup(s) % cnt + 1
              call copyToEnroute(s, ti, tj, tk, p, prtl_backup(s) % enroute(prtl_backup(s) % cnt))
#ifdef DEBUG
              if (prtl_backup(s) % cnt .gt. prtl_backup(s) % max) then
                call throwError('ERROR: something went wrong in `backupParticles`: cnt > max.')
              end if
#endif
            end do
          end do
        end do
      end do

    end do
  end subroutine backupParticles

  subroutine restoreParticlesFromBackup()
    implicit none
    integer :: s, p, ti, tj, tk
    do s = 1, nspec
      do p = 1, prtl_backup(s) % cnt
        ti = 1; tj = 1; tk = 1
#if defined(oneD) || defined (twoD) || defined (threeD)
        ti = FLOOR(REAL(prtl_backup(s) % enroute(p) % xi) / REAL(species(s) % tile_sx)) + 1
#endif
#if defined (twoD) || defined (threeD)
        tj = FLOOR(REAL(prtl_backup(s) % enroute(p) % yi) / REAL(species(s) % tile_sy)) + 1
#endif
#if defined (threeD)
        tk = FLOOR(REAL(prtl_backup(s) % enroute(p) % zi) / REAL(species(s) % tile_sz)) + 1
#endif
        ti = MIN(MAX(ti, 1), species(s) % tile_nx)
        tj = MIN(MAX(tj, 1), species(s) % tile_ny)
        tk = MIN(MAX(tk, 1), species(s) % tile_nz)
        call putEnrouteParticleOnTile(s, ti, tj, tk, prtl_backup(s) % enroute(p))
      end do
    end do
  end subroutine restoreParticlesFromBackup

  subroutine deallocateParticles()
    implicit none
    integer :: s

    do s = 1, nspec
      if (allocated(species(s) % prtl_tile)) deallocate (species(s) % prtl_tile)
    end do
  end subroutine deallocateParticles

  subroutine deallocateParticleBackup()
    implicit none
    if (allocated(prtl_backup)) deallocate (prtl_backup)
  end subroutine deallocateParticleBackup

end module m_particlelogistics
