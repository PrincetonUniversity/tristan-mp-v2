#include "../defs.F90"

module m_particledownsampling
#ifdef DOWNSAMPLING

  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_errors
  use m_domain
  use m_particles
  use m_particlelogistics
  use m_particlebinning
  implicit none

  ! Auxiliary type for a group of merging particles
  type :: particleDwnGroup
    ! total momentum components and energy of the group
    real                  :: tot_px, tot_py, tot_pz, tot_en
    ! total weight of the group
    real                  :: tot_wei
    ! indices of particles on a tile contained in the group
    integer, allocatable  :: indices(:)
    ! size of the group
    integer               :: size
    ! average momentum unit vector of the bin
    real                  :: bin_px, bin_py, bin_pz
  end type particleDwnGroup

  integer           :: dwn_start, dwn_interval
  real              :: dwn_maxweight, dwn_mom_spread
  logical           :: dwn_cartesian_bins, dwn_int_weights, dwn_dynamic_bins

  !--- PRIVATE variables/functions -------------------------------!
  private :: downsampleParticles,&
           & downsampleOnTile_Spherical, downsampleAllBins_Spherical, downsampleBin_Spherical,&
           & downsampleOnTile_Cartesian, downsampleAllBins_Cartesian, downsampleBin_Cartesian,&
           & mergeParticlesInGroup
  !...............................................................!
contains
  subroutine downsamplingStep(timestep)
    implicit none
    integer, intent(in) :: timestep
    if ((timestep .ge. dwn_start) .and.&
      & (modulo(timestep, dwn_interval) .eq. 0)) then
      call downsampleParticles()
    end if
    call printDiag("downsamplingStep()", 2)
  end subroutine downsamplingStep

  subroutine downsampleParticles()
    implicit none
    integer :: s, ti, tj, tk
    integer :: bin_limit

    type(positionBin_XYZ), allocatable  :: position_grid(:,:,:)
    type(particle_tile), allocatable    :: downsampling_tile
    integer                             :: nx_bin, ny_bin, nz_bin
    integer                             :: pi, pj, pk, p_ind, p
    #ifdef DEBUG
      real                                        :: nbinned, tmpnpart
    #endif

    do s = 1, nspec
      if (species(s)%dwn_sp .and. species(s)%ch_sp .ne. 0) then
        ! merging of charged particles based on cells

        do ti = 1, species(s)%tile_nx
          do tj = 1, species(s)%tile_ny
            do tk = 1, species(s)%tile_nz

              #ifdef DEBUG
                nbinned = 0
                tmpnpart = real(species(s)%prtl_tile(ti, tj, tk)%npart_sp)
              #endif

              nx_bin = INT(species(s)%tile_sx)
              ny_bin = INT(species(s)%tile_sy)
              nz_bin = INT(species(s)%tile_sz)

              call initializePositionBins(species(s)%prtl_tile(ti, tj, tk), position_grid, nx_bin, ny_bin, nz_bin)
              call binParticlePositions(species(s)%prtl_tile(ti, tj, tk), position_grid, nx_bin, ny_bin, nz_bin)

              if (allocated(downsampling_tile)) deallocate(downsampling_tile)
              allocate(downsampling_tile)

                do pi = 1, nx_bin
                  do pj = 1, ny_bin
                    do pk = 1, nz_bin

                    call allocateParticlesOnEmptyTile(s, downsampling_tile, position_grid(pi, pj, pk)%npart)

                      downsampling_tile%spec = s
                      downsampling_tile%npart_sp = position_grid(pi, pj, pk)%npart

                      downsampling_tile%x1 = (ti - 1) * species(s)%tile_sx
                      downsampling_tile%x2 = min(ti * species(s)%tile_sx, this_meshblock%ptr%sx)
                      downsampling_tile%y1 = (tj - 1) * species(s)%tile_sy
                      downsampling_tile%y2 = min(tj * species(s)%tile_sy, this_meshblock%ptr%sy)
                      downsampling_tile%z1 = (tk - 1) * species(s)%tile_sz
                      downsampling_tile%z2 = min(tk * species(s)%tile_sz, this_meshblock%ptr%sz)

                    #ifdef DEBUG
                      nbinned = nbinned + position_grid(pi, pj, pk)%npart
                    #endif

                    do p = 1, position_grid(pi, pj, pk)%npart
                      p_ind = position_grid(pi, pj, pk)%indices(p)
                      call fillDownsamplingTile(species(s)%prtl_tile(ti, tj, tk), downsampling_tile, p, p_ind)
                    enddo

                      if (downsampling_tile%npart_sp .gt. 5) then
                        ! decide whether to use cartesian OR spherical binning
                        if (dwn_cartesian_bins) then
                          call downsampleOnTile_Cartesian(downsampling_tile)
                        else
                          call downsampleOnTile_Spherical(downsampling_tile)
                        end if
                      end if

                    do p = 1, position_grid(pi, pj, pk)%npart
                      p_ind = position_grid(pi, pj, pk)%indices(p)
                      species(s)%prtl_tile(ti, tj, tk)%proc(p_ind) = downsampling_tile%proc(p)
                    enddo

                  enddo
                enddo
              enddo

              #ifdef DEBUG
                if (nbinned.ne.tmpnpart) then
                  print *, nbinned, tmpnpart
                  call throwError('[downsampleParticles] Unequal number of particles in tile and position bins.')
                endif
              #endif

              if (allocated(downsampling_tile)) deallocate(downsampling_tile)

            end do ! loop tk
          end do ! loop tj
        end do ! loop ti
      else if (species(s)%dwn_sp .and. species(s)%ch_sp .eq. 0) then
        ! merging of photons based on tiles
        do ti = 1, species(s)%tile_nx
          do tj = 1, species(s)%tile_ny
            do tk = 1, species(s)%tile_nz
              if (species(s)%prtl_tile(ti, tj, tk)%npart_sp .gt. 5) then
                ! decide whether to use cartesian OR spherical binning
                if (dwn_cartesian_bins) then
                  call downsampleOnTile_Cartesian(species(s)%prtl_tile(ti, tj, tk))
                else
                  call downsampleOnTile_Spherical(species(s)%prtl_tile(ti, tj, tk))
                end if
              end if
            end do
          end do
        end do
      end if
    end do

  end subroutine downsampleParticles

  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  ! Spherical momenta binning.
  subroutine downsampleOnTile_Spherical(tile)
    implicit none
    type(particle_tile), intent(inout)  :: tile
    type(momentumBin_Sph), allocatable  :: momentum_bins(:)
    real                                :: rot_ax_1, rot_ax_2, rot_ang
    ! generate a random rotation axis and a random rotation angle for a tile
    rot_ax_1 = random(dseed)
    rot_ax_2 = random(dseed)
    rot_ang = random(dseed)

    call initializeMomentumBins_Spherical(momentum_bins, tile%npart_sp)
    call binParticlesOnTile_Spherical(momentum_bins, tile, rot_ax_1, rot_ax_2, rot_ang)
    call downsampleAllBins_Spherical(momentum_bins, tile, rot_ax_1, rot_ax_2, rot_ang)
  end subroutine downsampleOnTile_Spherical

  subroutine downsampleAllBins_Spherical(momentum_bins, tile, ax1, ax2, ang)
    implicit none
    type(particle_tile), intent(inout)            :: tile
    type(momentumBin_Sph), allocatable, intent(inout) :: momentum_bins(:)
    integer                                       :: e_b, th_b, ph_b, p_ind, p, npart
    real, intent(in)                              :: ax1, ax2, ang
    #ifdef DEBUG
      real                                        :: en, u, v, w, theta, phi
    #endif

    ! loop through all the bins...
    ! ... in all 3 values (energy, theta, phi)
    do e_b = 0, dwn_n_energy_bins - 1
      do th_b = 0, momentum_bins(e_b)%n_theta_bins + 1
        do ph_b = 0, momentum_bins(e_b)%theta_bins(th_b)%n_phi_bins - 1
          npart = momentum_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%npart
          #ifdef DEBUG
            do p_ind = 1, npart
              p = momentum_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%indices(p_ind)
              if ((p .le. 0) .or. (p .gt. tile%npart_sp)) then
                call throwError('Wrong index in `downsampleAllBinsSpherical()`.')
              end if
            end do
          #endif
          if (npart .gt. 2) then
            call downsampleBin_Spherical(tile,&
                        & momentum_bins(e_b)%theta_bins(th_b)%theta_mid,&
                        & momentum_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%phi_mid,&
                        & momentum_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%indices, npart,&
                        & ax1, ax2, ang)
          end if
        end do
      end do
    end do
  end subroutine downsampleAllBins_Spherical

  ! on each bin we are forming groups of particles...
  ! ... with cumulative weights less than `dwn_maxweight`...
  ! ... and sending them to merge into separate routine
  subroutine downsampleBin_Spherical(tile, theta_mid, phi_mid, indices, npart,&
                         & ax1, ax2, ang)
    implicit none
    type(particle_tile), intent(inout)  :: tile
    integer, allocatable, intent(inout) :: indices(:)
    integer, intent(inout)              :: npart
    real, intent(in)                    :: theta_mid, phi_mid
    real, intent(in)                    :: ax1, ax2, ang
    type(particleDwnGroup)              :: group
    integer       :: p_ind, p, s
    real          :: en
    logical       :: masslessQ

    s = tile%spec
    if ((species(s)%m_sp .eq. 0) .and. (species(s)%ch_sp .eq. 0)) then
      masslessQ = .true.
    else
      masslessQ = .false.
    end if

    allocate(group%indices(npart))
    group%indices(:) = -1
    group%size = 0
    group%tot_px = 0.0; group%tot_py = 0.0; group%tot_pz = 0.0
    group%tot_en = 0.0; group%tot_wei = 0.0

    group%bin_px = cos(theta_mid) * cos(phi_mid)
    group%bin_py = cos(theta_mid) * sin(phi_mid)
    group%bin_pz = sin(theta_mid)
    ! rotate the bin center back to match the binned particles ...
    ! ... notice that angle is now `-ang` since we are rotating back
    call rotateRandomlyIn3D(group%bin_px, group%bin_py, group%bin_pz,&
                          & ax1, ax2, -ang)

    p_ind = 1
    do while (p_ind .le. npart)
      p = indices(p_ind)
      if (tile%weight(p) .gt. dwn_maxweight) then
        ! particle too heavy to merge
        indices(p_ind) = indices(npart)
        npart = npart - 1
        cycle
      else
        group%tot_wei = group%tot_wei + tile%weight(p)
        group%tot_px = group%tot_px + tile%weight(p) * tile%u(p)
        group%tot_py = group%tot_py + tile%weight(p) * tile%v(p)
        group%tot_pz = group%tot_pz + tile%weight(p) * tile%w(p)
        if (masslessQ) then
          en = sqrt(tile%u(p)**2 + tile%v(p)**2 + tile%w(p)**2)
        else
          en = sqrt(1.0 + tile%u(p)**2 + tile%v(p)**2 + tile%w(p)**2)
        end if
        group%tot_en = group%tot_en + tile%weight(p) * en

        group%indices(group%size + 1) = p
        group%size = group%size + 1

        if ((group%tot_wei .ge. dwn_maxweight) .or. (p_ind .eq. npart)) then
          ! once there are enough particles in the group...
          ! ... send a group of these particles to merge...
          ! ... then reset the quantities
          if (group%size .gt. 5) then
            call mergeParticlesInGroup(group, tile)
          end if
          group%indices(:) = -1
          group%size = 0
          group%tot_px = 0.0; group%tot_py = 0.0; group%tot_pz = 0.0
          group%tot_en = 0.0; group%tot_wei = 0.0
        end if

        p_ind = p_ind + 1
      end if
    end do
  end subroutine downsampleBin_Spherical
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  ! Cartesian momenta binning.
  subroutine downsampleOnTile_Cartesian(tile)
    implicit none
    type(particle_tile), intent(inout)  :: tile
    type(momentumBin_XYZ), allocatable  :: momentum_bins(:,:,:)
    integer                             :: p
    real                                :: rot_ax_1, rot_ax_2, rot_ang
    real                                :: px_min, px_max, py_min, py_max, pz_min, pz_max
    real                                :: px_mid, py_mid, pz_mid
    ! generate a random rotation axis and a random rotation angle for a tile
    rot_ax_1 = random(dseed)
    rot_ax_2 = random(dseed)

    if (.not. dwn_dynamic_bins) then
      rot_ang = random(dseed)
      px_min = -dwn_energy_max
      px_max = dwn_energy_max
      py_min = -dwn_energy_max
      py_max = dwn_energy_max
      pz_min = -dwn_energy_max
      pz_max = dwn_energy_max
    else

      rot_ang = 0.0

      px_min = MINVAL(tile%u(1 : tile%npart_sp))
      py_min = MINVAL(tile%v(1 : tile%npart_sp))
      pz_min = MINVAL(tile%w(1 : tile%npart_sp))
      px_max = MAXVAL(tile%u(1 : tile%npart_sp))
      py_max = MAXVAL(tile%v(1 : tile%npart_sp))
      pz_max = MAXVAL(tile%w(1 : tile%npart_sp))

      if (px_min.eq.px_max) then
        px_min = px_min - 1e-5
        px_max = px_max + 1e-5
      endif

      if (py_min.eq.py_max) then
        py_min = py_min - 1e-5
        py_max = py_max + 1e-5
      endif

        if (pz_min.eq.pz_max) then
        pz_min = pz_min - 1e-5
        pz_max = pz_max + 1e-5
      endif

      px_min = px_min * (1.0 - sign(1.0, px_min) * 0.01)
      py_min = py_min * (1.0 - sign(1.0, py_min) * 0.01)
      pz_min = pz_min * (1.0 - sign(1.0, pz_min) * 0.01)
      px_max = px_max * (1.0 + sign(1.0, px_max) * 0.01)
      py_max = py_max * (1.0 + sign(1.0, py_max) * 0.01)
      pz_max = pz_max * (1.0 + sign(1.0, pz_max) * 0.01)

      px_mid = SUM(tile%u(1 : tile%npart_sp)) / tile%npart_sp
      py_mid = SUM(tile%v(1 : tile%npart_sp)) / tile%npart_sp
      pz_mid = SUM(tile%w(1 : tile%npart_sp)) / tile%npart_sp

      px_min = MAX(px_min, px_mid - 0.5 * dwn_mom_spread)
      px_max = MIN(px_max, px_mid + 0.5 * dwn_mom_spread)
      py_min = MAX(py_min, py_mid - 0.5 * dwn_mom_spread)
      py_max = MIN(py_max, py_mid + 0.5 * dwn_mom_spread)
      pz_min = MAX(pz_min, pz_mid - 0.5 * dwn_mom_spread)
      pz_max = MIN(pz_max, pz_mid + 0.5 * dwn_mom_spread)
    end if

    call initializeMomentumBins_Cartesian(momentum_bins, tile%npart_sp,&
                                        & px_min, px_max, py_min, py_max, pz_min, pz_max)
    call binParticlesOnTile_Cartesian(momentum_bins, tile, rot_ax_1, rot_ax_2, rot_ang,&
                                    & px_min, px_max, py_min, py_max, pz_min, pz_max)

    call downsampleAllBins_Cartesian(momentum_bins, tile, rot_ax_1, rot_ax_2, rot_ang)

  end subroutine downsampleOnTile_Cartesian

  subroutine downsampleAllBins_Cartesian(momentum_bins, tile, ax1, ax2, ang)
    implicit none
    type(particle_tile), intent(inout)            :: tile
    type(momentumBin_XYZ), allocatable, intent(inout) :: momentum_bins(:,:,:)
    real, intent(in)                              :: ax1, ax2, ang
    integer :: p, pi, pj, pk, npart
    real    :: px_mid, py_mid, pz_mid

    ! loop through all the bins
    do pi = 1, dwn_n_mom_bins
      do pj = 1, dwn_n_mom_bins
        do pk = 1, dwn_n_mom_bins
          npart = momentum_bins(pi, pj, pk)%npart

          if (npart .gt. 5) then
            px_mid = 0.5 * (momentum_bins(pi, pj, pk)%px_max + momentum_bins(pi, pj, pk)%px_min)
            py_mid = 0.5 * (momentum_bins(pi, pj, pk)%py_max + momentum_bins(pi, pj, pk)%py_min)
            pz_mid = 0.5 * (momentum_bins(pi, pj, pk)%pz_max + momentum_bins(pi, pj, pk)%pz_min)
            call downsampleBin_Cartesian(tile, px_mid, py_mid, pz_mid,&
                                       & momentum_bins(pi, pj, pk)%indices, npart,&
                                       & ax1, ax2, ang)
          end if
        end do
      end do
    end do
  end subroutine downsampleAllBins_Cartesian

  subroutine downsampleBin_Cartesian(tile, px_mid, py_mid, pz_mid,&
                                   & indices, npart, ax1, ax2, ang)
    implicit none
    type(particle_tile), intent(inout)  :: tile
    integer, allocatable, intent(inout) :: indices(:)
    integer, intent(inout)              :: npart
    real, intent(in)                    :: px_mid, py_mid, pz_mid
    real, intent(in)                    :: ax1, ax2, ang
    type(particleDwnGroup)              :: group
    integer       :: p_ind, p, s
    real          :: en
    logical       :: masslessQ

    s = tile%spec
    if ((species(s)%m_sp .eq. 0) .and. (species(s)%ch_sp .eq. 0)) then
      masslessQ = .true.
    else
      masslessQ = .false.
    end if

    allocate(group%indices(npart))
    group%indices(:) = -1
    group%size = 0
    group%tot_px = 0.0; group%tot_py = 0.0; group%tot_pz = 0.0
    group%tot_en = 0.0; group%tot_wei = 0.0

    group%bin_px = px_mid
    group%bin_py = py_mid
    group%bin_pz = pz_mid
    ! rotate the bin center back to match the binned particles ...
    ! ... notice that angle is now `-ang` since we are rotating back
    call rotateRandomlyIn3D(group%bin_px, group%bin_py, group%bin_pz,&
                          & ax1, ax2, -ang)

    p_ind = 1
    do while (p_ind .le. npart)
      p = indices(p_ind)
      if (tile%weight(p) .gt. dwn_maxweight) then
        ! particle too heavy to merge
        indices(p_ind) = indices(npart)
        npart = npart - 1
        cycle
      else
        group%tot_wei = group%tot_wei + tile%weight(p)
        group%tot_px = group%tot_px + tile%weight(p) * tile%u(p)
        group%tot_py = group%tot_py + tile%weight(p) * tile%v(p)
        group%tot_pz = group%tot_pz + tile%weight(p) * tile%w(p)
        if (masslessQ) then
          en = sqrt(tile%u(p)**2 + tile%v(p)**2 + tile%w(p)**2)
        else
          en = sqrt(1.0 + tile%u(p)**2 + tile%v(p)**2 + tile%w(p)**2)
        end if
        group%tot_en = group%tot_en + tile%weight(p) * en

        group%indices(group%size + 1) = p
        group%size = group%size + 1

        if ((group%tot_wei .ge. dwn_maxweight) .or. (p_ind .eq. npart)) then
          ! once there are enough particles in the group...
          ! ... send a group of these particles to merge...
          ! ... then reset the quantities

          if (group%size .gt. 5) then
            call mergeParticlesInGroup(group, tile)
          end if
          group%indices(:) = -1
          group%size = 0
          group%tot_px = 0.0; group%tot_py = 0.0; group%tot_pz = 0.0
          group%tot_en = 0.0; group%tot_wei = 0.0
        end if

        p_ind = p_ind + 1
      end if
    end do

  end subroutine downsampleBin_Cartesian
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !   The algorithm is adapted from Vranic et al. 2014 [1411.2248v1]
  subroutine mergeParticlesInGroup(group, tile)
    implicit none
    type(particleDwnGroup), intent(in)  :: group
    type(particle_tile), intent(inout)  :: tile
    real            :: wA, wB, dw
    integer         :: p, p_ind, s
    real            :: pxA, pxB, pyA, pyB, pzA, pzB, rnd
    real            :: dxA, dyA, dzA, dxB, dyB, dzB
    integer(kind=2) :: xAi, yAi, zAi, xBi, yBi, zBi
    real            :: pA, pB, enA, enB, tot_p, cos_th, sin_th
    real            :: temp1_x, temp1_y, temp1_z, temp1
    real            :: temp2_x, temp2_y, temp2_z, temp2
    logical         :: masslessQ
    real            :: x2, y2, z2, totweight
    real            :: x1, y1, z1

    s = tile%spec
    if ((species(s)%m_sp .eq. 0) .and. (species(s)%ch_sp .eq. 0)) then
      masslessQ = .true.
    else
      masslessQ = .false.
    end if

    ! There are two directions in this problem...
    ! ... parallel to the total momentum...
    ! ... and perpendicular (which is arbitrary)

    ! We choose the perpendicular direction,..
    ! ... so that particles are created in the plane...
    ! ... containing both the total momentum...
    ! ... and the average momentum of the bin

    temp1_x = group%bin_pz * group%tot_py - group%bin_py * group%tot_pz
    temp1_y = -group%bin_pz * group%tot_px + group%bin_px * group%tot_pz
    temp1_z = group%bin_py * group%tot_px - group%bin_px * group%tot_py

    temp2_x = group%tot_pz * temp1_y - group%tot_py * temp1_z
    temp2_y = -group%tot_pz * temp1_x + group%tot_px * temp1_z
    temp2_z = group%tot_py * temp1_x - group%tot_px * temp1_y
    temp2 = sqrt(temp2_x**2 + temp2_y**2 + temp2_z**2)

    ! Now vector `temp2` is perpendicular to `tot_p`...
    ! ... but also lies in the plane containing the...
    ! ... average momentum of the bin

    tot_p = sqrt(group%tot_px**2 + group%tot_py**2 + group%tot_pz**2)

    ! new weights
    ! ... weights of A and B are a matter of choice ...
    wA = group%tot_wei / 2.0
    wB = group%tot_wei / 2.0
    ! ... but the code below works even if wA != wB ...
    ! ... (as long as their sum is tot_wei)

    ! merge into 2 particles w/ integer weights. only possible if the total ...
    ! ... group weight is an integer.
    if (dwn_int_weights .and. (abs(INT(group%tot_wei) - group%tot_wei) .eq. 0.0)) then
      dw = wA - INT(wA)
      wA = wA - dw
      wB = wB + dw
    endif

    if (masslessQ) then
      ! new energies & momenta (magnitudes)
      enA = group%tot_en / (2.0 * wA)
      pA = enA
      enB = group%tot_en / (2.0 * wB)
      pB = enB

      ! direction between new particle momenta...
      ! ... and the total group momentum
      cos_th = tot_p / group%tot_en
      sin_th = sqrt(1.0 - cos_th**2)
    else
      ! new energies & momenta (magnitudes)
      enA = (group%tot_en**2 + wA**2 - wB**2) /&
            & (2.0 * group%tot_en * wA)
      pA = sqrt(enA**2 - 1.0)
      enB = (group%tot_en**2 + wB**2 - wA**2) /&
            & (2.0 * group%tot_en * wB)
      pB = sqrt(enB**2 - 1.0)

      ! direction between new particle momenta...
      ! ... and the total group momentum
      cos_th = MIN(tot_p / (2.0 * wA * pA), 1.0)
      sin_th = sqrt(1.0 - cos_th**2)
    end if

    ! new momenta
    pxA = (group%tot_px / tot_p) * cos_th * pA +& ! parallel component
        & (temp2_x / temp2) * sin_th * pA   ! perp component
    pyA = (group%tot_py / tot_p) * cos_th * pA +&
        & (temp2_y / temp2) * sin_th * pA
    pzA = (group%tot_pz / tot_p) * cos_th * pA +&
        & (temp2_z / temp2) * sin_th * pA

    pxB = (group%tot_px / tot_p) * cos_th * pB -& ! parallel component
        & (temp2_x / temp2) * sin_th * pB   ! perp component
    pyB = (group%tot_py / tot_p) * cos_th * pB -&
        & (temp2_y / temp2) * sin_th * pB
    pzB = (group%tot_pz / tot_p) * cos_th * pB -&
        & (temp2_z / temp2) * sin_th * pB

    #ifdef DEBUG
      ! check weight conservation
      if (wA + wB .ne. group%tot_wei) then
        print *, wA, wB, group%tot_wei
        call throwError('Weight is not conserved in `mergeParticlesInGroup()`')
      end if
      if (.not. numbersAreClose(pxA * wA + pxB * wB, group%tot_px)) then
        call throwError('Px is not conserved in `mergeParticlesInGroup()`')
      end if
      if (.not. numbersAreClose(pyA * wA + pyB * wB, group%tot_py)) then
        call throwError('Py is not conserved in `mergeParticlesInGroup()`')
      end if
      if (.not. numbersAreClose(pzA * wA + pzB * wB, group%tot_pz)) then
        call throwError('Pz is not conserved in `mergeParticlesInGroup()`')
      end if
      if (.not. numbersAreClose(enA * wA + enB * wB, group%tot_en)) then
        call throwError('Energy is not conserved in `mergeParticlesInGroup()`')
      end if
      if ((.not. numbersAreClose(pxA**2 + pyA**2 + pzA**2, pA**2)) .or.&
        & (.not. numbersAreClose(pxB**2 + pyB**2 + pzB**2, pB**2))) then
        call throwError('Wrong momenta projections in `mergeParticlesInGroup()`')
      end if
    #endif

    ! Determine center of mass of particles to-be-merged
    x2 = 0.0
    y2 = 0.0
    z2 = 0.0
    totweight = 0.0

    do p_ind = 1, group%size
      p = group%indices(p_ind)

      x2 = x2 + tile%weight(p) * (real(tile%xi(p)) + tile%dx(p))
      y2 = y2 + tile%weight(p) * (real(tile%yi(p)) + tile%dy(p))
      z2 = z2 + tile%weight(p) * (real(tile%zi(p)) + tile%dz(p))

      totweight = real(totweight + tile%weight(p))
    enddo

    x2 = x2 / totweight
    y2 = y2 / totweight
    z2 = z2 / totweight

    xBi = floor(x2); dxB = x2 - real(xBi)
    yBi = floor(y2); dyB = y2 - real(yBi)
    zBi = floor(z2); dzB = z2 - real(zBi)

    xAi = xBi; dxA = dxB
    yAi = yBi; dyA = dyB
    zAi = zBi; dzA = dzB

    ! Extra current deposit for center-of-mass shift
    do p_ind = 1, group%size
      p = group%indices(p_ind)

      x1 = real(tile%xi(p)) + tile%dx(p)
      y1 = real(tile%yi(p)) + tile%dy(p)
      z1 = real(tile%zi(p)) + tile%dz(p)

      call depositCurrentsFromSingleParticle(s, tile, p, x1, y1, z1, x2, y2, z2)
    end do

    ! "nullify" merged particles
    do p_ind = 1, group%size
      p = group%indices(p_ind)
      tile%proc(p) = -1
    end do

    ! inject new particles
    call createParticle(s, xAi, yAi, zAi, dxA, dyA, dzA, pxA, pyA, pzA, weight=wA)
    call createParticle(s, xBi, yBi, zBi, dxB, dyB, dzB, pxB, pyB, pzB, weight=wB)

  end subroutine mergeParticlesInGroup
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
end module m_particledownsampling
