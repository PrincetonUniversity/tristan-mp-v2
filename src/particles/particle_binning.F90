#include "../defs.F90"

module m_particlebinning
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  implicit none

  ! min energy for merging particles
  real              :: dwn_energy_min, dwn_energy_max

  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  ! Spherical momenta binning...
  ! ... auxiliary types
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! binning is logarithmic or linear in energy,..
  ! ... uniform in `theta`,..
  ! ... and approximately uniform in `phi`
  ! each `energy` bin contains `theta`-bins...
  ! ... which further contains `phi`-bins
  type :: phiBin_Sph
    ! left and right bounds of the sub-bin (for debugging)
    real                 :: phi_min, phi_max, phi_mid
    ! number of particles in the sub-bin
    integer              :: npart
    ! particle indices in this sub-bin
    integer, allocatable :: indices(:)
  end type phiBin_Sph

  type :: thetaBin_Sph
    ! upper and lower boundes of the bin (for debugging)
    real                       :: theta_min, theta_max, theta_mid
    ! number of sub-bins in phi
    integer                    :: n_phi_bins
    ! sub-bins in phi
    type(phiBin_Sph), allocatable  :: phi_bins(:)
  end type thetaBin_Sph

  type :: momentumBin_Sph
    ! upper and lower bounds for the energy bin
    real                          :: e_min, e_max
    ! polar bin opening angle (half of it)
    real                          :: th0_bin
    ! number of sub-bins in theta
    integer                       :: n_theta_bins
    ! sub-bins in theta
    type(thetaBin_Sph), allocatable   :: theta_bins(:)
  end type momentumBin_Sph

  ! number of bins
  integer                         :: dwn_n_angular_bins, dwn_n_energy_bins
  ! log or linearly distributed energy bins
  logical                         :: dwn_log_e_bins

  !--- PRIVATE variables/functions -------------------------------!
  private :: findEnergyBin_Sph, findThetaBin_Sph, findPhiBin_Sph
  private :: initializeThetaBins_Sph, initializePhiBins_Sph
  !...............................................................!
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  ! Cartesian momenta binning...
  ! ... auxiliary types
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  type :: momentumBin_XYZ
    ! upper and lower bounds for the momentum components
    real                    :: px_min, px_max
    real                    :: py_min, py_max
    real                    :: pz_min, pz_max
    ! # of particles in the bin
    integer                 :: npart
    ! indices of all the particles in a given bin
    integer, allocatable    :: indices(:)
  end type momentumBin_XYZ

  ! number of bins
  integer                   :: dwn_n_mom_bins
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  ! Position binning ...
  ! ... auxiliary types
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  type :: positionBin_XYZ
    ! # of particles in the bin
    integer                 :: npart
    ! indices of all the particles in a given bin
    integer, allocatable    :: indices(:)
  end type positionBin_XYZ
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
contains

  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  ! Spherical momenta binning.
  subroutine binParticlesOnTile_Spherical(momentum_bins, tile, ax1, ax2, ang)
    implicit none
    ! FIX: this is for photons only
    type(momentumBin_Sph), allocatable, intent(inout) :: momentum_bins(:)
    type(particle_tile), intent(in)               :: tile
    real, intent(in)                              :: ax1, ax2, ang
    integer                                       :: s
    logical                                       :: masslessQ

    integer :: p
    integer :: energy_ind, theta_ind, phi_ind
    real    :: prtl_ux, prtl_uy, prtl_uz, prtl_energy, prtl_gamma
    real    :: u_theta, u_phi
    integer :: dummy_int

    s = tile%spec
    if ((species(s)%m_sp .eq. 0) .and. (species(s)%ch_sp .eq. 0)) then
      masslessQ = .true.
    else
      masslessQ = .false.
    end if

    do p = 1, tile%npart_sp
      prtl_ux = tile%u(p); prtl_uy = tile%v(p); prtl_uz = tile%w(p)
      ! rotation (not really random, because axis and angles are passed)
      call rotateRandomlyIn3D(prtl_ux, prtl_uy, prtl_uz, ax1, ax2, ang)
      prtl_energy = sqrt(prtl_ux**2 + prtl_uy**2 + prtl_uz**2)
      if (.not. masslessQ) then
        prtl_gamma = sqrt(1.0 + prtl_ux**2 + prtl_uy**2 + prtl_uz**2)
      end if
      prtl_ux = prtl_ux / prtl_energy
      prtl_uy = prtl_uy / prtl_energy
      prtl_uz = prtl_uz / prtl_energy
      if (.not. masslessQ) then
        ! binning by `lorentz-factor - 1` (kinetic energy)
        prtl_energy = prtl_gamma - 1.0
      end if
      if ((prtl_energy .ge. momentum_bins(0)%e_min) .and.&
        & (prtl_energy .lt. momentum_bins(dwn_n_energy_bins - 1)%e_max)) then
        u_theta = asin(prtl_uz)
        u_phi = atan2(prtl_uy, prtl_ux)
        if (u_phi .lt. 0) u_phi = u_phi + 2 * M_PI
        call findEnergyBin_Sph(momentum_bins,&
                        & prtl_energy, energy_ind)
        call findThetaBin_Sph(momentum_bins(energy_ind),&
                        & u_theta, theta_ind)
        dummy_int = momentum_bins(energy_ind)%theta_bins(theta_ind)%n_phi_bins
        call findPhiBin_Sph(momentum_bins(energy_ind)%theta_bins(theta_ind),&
                      & u_phi, phi_ind)
        dummy_int =&
            & momentum_bins(energy_ind)%theta_bins(theta_ind)%phi_bins(phi_ind)%npart
        momentum_bins(energy_ind)%theta_bins(theta_ind)%phi_bins(phi_ind)%npart = dummy_int + 1
        momentum_bins(energy_ind)%theta_bins(theta_ind)%phi_bins(phi_ind)%&
                                  &indices(dummy_int + 1) = p
      end if
    end do
  end subroutine binParticlesOnTile_Spherical

  ! - - - finding bins - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine findEnergyBin_Sph(momentum_bins, energy, en_ind)
    implicit none
    type(momentumBin_Sph), intent(in), allocatable  :: momentum_bins(:)
    real, intent(in)                            :: energy
    integer, intent(out)                        :: en_ind
    integer                                     :: e_b
    en_ind = -1
    do e_b = 0, dwn_n_energy_bins - 1
      if ((energy .lt. momentum_bins(e_b)%e_max) .and.&
         &(energy .ge. momentum_bins(e_b)%e_min))then
        en_ind = e_b
        exit
      end if
    end do
    #ifdef DEBUG
      if ((en_ind .lt. 0) .or. (en_ind .ge. dwn_n_energy_bins)) then
        call throwError('Something is wrong in `findEnergyBin_Sph()`')
      end if
    #endif
  end subroutine findEnergyBin_Sph

  subroutine findThetaBin_Sph(momentum_bin, u_theta, th_ind)
    implicit none
    type(momentumBin_Sph), intent(in) :: momentum_bin
    real, intent(in)              :: u_theta
    integer, intent(out)          :: th_ind
    real                          :: d_theta

    d_theta = (M_PI - 2 * momentum_bin%th0_bin) / momentum_bin%n_theta_bins
    if (u_theta .le. -0.5 * M_PI + momentum_bin%th0_bin) then
      ! if on the southern pole bin
      th_ind = 0
    else if (u_theta .gt. 0.5 * M_PI - momentum_bin%th0_bin) then
      ! if on the northern pole bin
      th_ind = momentum_bin%n_theta_bins + 1
    else
      th_ind = INT((u_theta + 0.5 * M_PI - momentum_bin%th0_bin) / d_theta) + 1
      th_ind = MIN(th_ind, momentum_bin%n_theta_bins)
      #ifdef DEBUG
        if ((th_ind .gt. momentum_bin%n_theta_bins) .or. (th_ind .le. 0)) then
          call throwError('Something is wrong in `findThetaBin_Sph()`')
        end if
      #endif
    end if
  end subroutine findThetaBin_Sph

  subroutine findPhiBin_Sph(theta_bin, u_phi, ph_bin)
    implicit none
    type(thetaBin_Sph), intent(in)  :: theta_bin
    real, intent(in)            :: u_phi
    integer, intent(out)        :: ph_bin
    ph_bin = INT(u_phi * theta_bin%n_phi_bins / (2 * M_PI))
    ph_bin = MAX(0, MIN(ph_bin, theta_bin%n_phi_bins - 1))
    #ifdef DEBUG
      if ((ph_bin .lt. 0) .or. (ph_bin .ge. theta_bin%n_phi_bins)) then
        print *, u_phi, ph_bin, theta_bin%theta_min, theta_bin%theta_max
        call throwError('Something is wrong in `findPhiBin_Sph()`')
      end if
    #endif
  end subroutine findPhiBin_Sph

  ! - - - initializing bins - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine initializeMomentumBins_Spherical(momentum_bins, nparts_in_tile)
    implicit none
    integer, intent(in)                         :: nparts_in_tile
    type(momentumBin_Sph), intent(out), allocatable :: momentum_bins(:)
    integer                                     :: e_b
    real, allocatable                           :: normbins(:)
    real                                        :: e_min, e_max

    if (dwn_log_e_bins) then
      ! initializing randomized energy bins ...
      ! ... intervals of which have lognormal distribution
      call log_normal(dwn_n_energy_bins + 1, normbins)
      ! randomize maximum energy bin from `E` to `10 * E`...
      ! ... for more randomness
      e_min = dwn_energy_min
      e_max = 10**log10(dwn_energy_max) + random(dseed)
    else
      ! initializing randomized energy bins ...
      ! ... intervals of which have linear distribution
      call lin_normal(dwn_n_energy_bins + 1, normbins)
      ! randomize maximum energy bin from `E` to `2 * E`...
      ! ... for more randomness
      e_min = dwn_energy_min
      e_max = dwn_energy_max * (1.0 + random(dseed))
    end if

    allocate(momentum_bins(0 : dwn_n_energy_bins - 1))
    do e_b = 0, dwn_n_energy_bins - 1
      momentum_bins(e_b)%e_min = (e_min + (e_max - e_min) * normbins(e_b + 1))
      momentum_bins(e_b)%e_max = (e_min + (e_max - e_min) * normbins(e_b + 2))
      momentum_bins(e_b)%th0_bin = 0.5 * M_PI / dwn_n_angular_bins
      momentum_bins(e_b)%n_theta_bins = dwn_n_angular_bins
      call initializeThetaBins_Sph(momentum_bins(e_b), nparts_in_tile)
    end do
  end subroutine initializeMomentumBins_Spherical

  subroutine initializeThetaBins_Sph(momentum_bin, nparts_in_tile)
    implicit none
    type(momentumBin_Sph), intent(inout)        :: momentum_bin
    integer, intent(in)                         :: nparts_in_tile
    real                                        :: d_theta
    integer                                     :: th_b, ph_b, dummy4
    real                                        :: dummy1, dummy2, dummy3

    ! bins go from `0 -> n_theta_bins + 1`...
    ! ... with `n_theta_bins + 2` bins overall
    allocate(momentum_bin%theta_bins(0 : momentum_bin%n_theta_bins + 1))

    d_theta = (M_PI - 2 * momentum_bin%th0_bin) / momentum_bin%n_theta_bins

    do th_b = 0, momentum_bin%n_theta_bins + 1
      if (th_b .eq. 0) then
        ! south polar bin
        dummy1 = -0.5 * M_PI
        dummy2 = -0.5 * M_PI + momentum_bin%th0_bin
        dummy3 = -0.5 * M_PI
        dummy4 = 1
      else if (th_b .eq. momentum_bin%n_theta_bins + 1) then
        ! north polar bin
        dummy1 = 0.5 * M_PI - momentum_bin%th0_bin
        dummy2 = 0.5 * M_PI
        dummy3 = 0.5 * M_PI
        dummy4 = 1
      else
        dummy1 = -0.5 * M_PI + momentum_bin%th0_bin + d_theta * (th_b - 1)
        dummy2 = -0.5 * M_PI + momentum_bin%th0_bin + d_theta * th_b
        dummy3 = 0.5 * (dummy1 + dummy2)
        dummy4 = INT(2 * momentum_bin%n_theta_bins * cos(dummy3))
      end if
      momentum_bin%theta_bins(th_b)%theta_min = dummy1
      momentum_bin%theta_bins(th_b)%theta_max = dummy2
      momentum_bin%theta_bins(th_b)%theta_mid = dummy3
      momentum_bin%theta_bins(th_b)%n_phi_bins = dummy4
      call initializePhiBins_Sph(momentum_bin, momentum_bin%theta_bins(th_b), nparts_in_tile)
    end do
  end subroutine initializeThetaBins_Sph

  subroutine initializePhiBins_Sph(momentum_bin, theta_bin, nparts_in_tile)
    implicit none
    type(momentumBin_Sph), intent(inout)        :: momentum_bin
    type(thetaBin_Sph), intent(inout)           :: theta_bin
    integer, intent(in)                     :: nparts_in_tile
    integer                                 :: ph_b
    real                                    :: d_phi

    ! `phi` bins go from `0 -> n_ph_bins - 1`...
    ! ... with `n_ph_bins` bins overall
    allocate(theta_bin%phi_bins(0 : theta_bin%n_phi_bins - 1))
    if (abs(theta_bin%theta_mid) .ne. M_PI * 0.5) then
      d_phi = 2 * M_PI / REAL(theta_bin%n_phi_bins)
    else
      d_phi = 2 * M_PI
    end if
    do ph_b = 0, theta_bin%n_phi_bins - 1
      theta_bin%phi_bins(ph_b)%phi_min = d_phi * ph_b
      theta_bin%phi_bins(ph_b)%phi_max = d_phi * (ph_b + 1)
      theta_bin%phi_bins(ph_b)%phi_mid = d_phi * (ph_b + 0.5)
      theta_bin%phi_bins(ph_b)%npart = 0
      allocate(theta_bin%phi_bins(ph_b)%indices(nparts_in_tile))
    end do
  end subroutine initializePhiBins_Sph
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  ! Cartesian momenta binning.
  ! - - - initializing bins - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine initializeMomentumBins_Cartesian(momentum_bins, nparts_in_tile,&
                                            & px_min, px_max, py_min, py_max, pz_min, pz_max)
    implicit none
    integer, intent(in)                               :: nparts_in_tile
    type(momentumBin_XYZ), intent(out), allocatable   :: momentum_bins(:,:,:)
    real, intent(in)    :: px_min, px_max, py_min, py_max, pz_min, pz_max
    integer             :: pi, pj, pk
    real                :: del_ex, del_ey, del_ez

    allocate(momentum_bins(dwn_n_mom_bins, dwn_n_mom_bins, dwn_n_mom_bins))
    del_ex = (px_max - px_min) / dwn_n_mom_bins
    del_ey = (py_max - py_min) / dwn_n_mom_bins
    del_ez = (pz_max - pz_min) / dwn_n_mom_bins
    do pi = 1, dwn_n_mom_bins
      do pj = 1, dwn_n_mom_bins
        do pk = 1, dwn_n_mom_bins
          momentum_bins(pi, pj, pk)%px_min = px_min + (pi - 1) * del_ex
          momentum_bins(pi, pj, pk)%px_max = px_min + pi * del_ex
          momentum_bins(pi, pj, pk)%py_min = py_min + (pj - 1) * del_ey
          momentum_bins(pi, pj, pk)%py_max = py_min + pj * del_ey
          momentum_bins(pi, pj, pk)%pz_min = pz_min + (pk - 1) * del_ez
          momentum_bins(pi, pj, pk)%pz_max = pz_min + pk * del_ez
          momentum_bins(pi, pj, pk)%npart = 0
          allocate(momentum_bins(pi, pj, pk)%indices(nparts_in_tile))
        end do
      end do
    end do
  end subroutine initializeMomentumBins_Cartesian

  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  ! - - - initializing bins - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine initializePositionBins(tile, position_bins, n_tile_sx, n_tile_sy, n_tile_sz)
    implicit none
    type(particle_tile), intent(in)                   :: tile
    type(positionBin_XYZ), intent(out), allocatable   :: position_bins(:,:,:)
    integer, intent(in)                               :: n_tile_sx, n_tile_sy, n_tile_sz
    integer                                           :: pi, pj, pk
    integer                                           :: s

    allocate(position_bins(n_tile_sx, n_tile_sy, n_tile_sz))

    do pi = 1, n_tile_sx
      do pj = 1, n_tile_sy
        do pk = 1, n_tile_sz
          position_bins(pi, pj, pk)%npart = 0
          allocate(position_bins(pi, pj, pk)%indices(tile%npart_sp))
        end do
      end do
    end do
  end subroutine initializePositionBins

  subroutine fillDownsamplingTile(tile, dwn_tile, p, p_ind)
    implicit none
    ! DEP_PRT [particle-dependent]
    type(particle_tile), intent(in)                   :: tile
    type(particle_tile), intent(inout)                :: dwn_tile
    integer, intent(in)                               :: p
    integer, intent(in)                               :: p_ind

    dwn_tile%xi(p) = tile%xi(p_ind)
    dwn_tile%dx(p) = tile%dx(p_ind)

    dwn_tile%yi(p) = tile%yi(p_ind)
    dwn_tile%dy(p) = tile%dy(p_ind)

    dwn_tile%zi(p) = tile%zi(p_ind)
    dwn_tile%dz(p) = tile%dz(p_ind)

    dwn_tile%u(p) = tile%u(p_ind)
    dwn_tile%v(p) = tile%v(p_ind)
    dwn_tile%w(p) = tile%w(p_ind)

    dwn_tile%ind(p) = tile%ind(p_ind)
    dwn_tile%proc(p) = tile%proc(p_ind)

    dwn_tile%weight(p) = tile%weight(p_ind)

    #ifdef PRTLPAYLOADS
      dwn_tile%payload1(p) = tile%payload1(p_ind)
      dwn_tile%payload2(p) = tile%payload2(p_ind)
      dwn_tile%payload3(p) = tile%payload3(p_ind)
    #endif

  end subroutine fillDownsamplingTile

  subroutine binParticlePositions(tile, position_bins, n_tile_sx, n_tile_sy, n_tile_sz)
    implicit none

    type(positionBin_XYZ), allocatable, intent(inout) :: position_bins(:,:,:)
    integer, intent(in)                               :: n_tile_sx, n_tile_sy, n_tile_sz
    type(particle_tile), intent(in)   :: tile
    integer :: p, pi, pj, pk

    do p = 1, tile%npart_sp
      pi = tile%xi(p) - tile%x1 + 1
      pj = tile%yi(p) - tile%y1 + 1
      pk = tile%zi(p) - tile%z1 + 1

      pi = floor(REAL(mod(tile%xi(p),n_tile_sx))) + 1
      pj = floor(REAL(mod(tile%yi(p),n_tile_sy))) + 1
      pk = floor(REAL(mod(tile%zi(p),n_tile_sz))) + 1

      position_bins(pi, pj, pk)%npart = position_bins(pi, pj, pk)%npart + 1
      position_bins(pi, pj, pk)%indices(position_bins(pi, pj, pk)%npart) = p
    end do
  end subroutine binParticlePositions
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  subroutine binParticlesOnTile_Cartesian(momentum_bins, tile, ax1, ax2, ang,&
                                        & px_min, px_max, py_min, py_max, pz_min, pz_max)
    implicit none
    ! FIX: this is for photons only
    type(momentumBin_XYZ), allocatable, intent(inout) :: momentum_bins(:,:,:)
    type(particle_tile), intent(in)   :: tile
    real, intent(in)                  :: px_min, px_max, py_min, py_max, pz_min, pz_max
    real, intent(in)                  :: ax1, ax2, ang
    integer                           :: s
    logical                           :: masslessQ

    integer :: p, pi, pj, pk
    real    :: prtl_ux, prtl_uy, prtl_uz, del_ex, del_ey, del_ez, prtl_energy
    integer :: dummy_int

    s = tile%spec
    if ((species(s)%m_sp .eq. 0) .and. (species(s)%ch_sp .eq. 0)) then
      masslessQ = .true.
    else
      masslessQ = .false.
    end if

    del_ex = (px_max - px_min) / dwn_n_mom_bins
    del_ey = (py_max - py_min) / dwn_n_mom_bins
    del_ez = (pz_max - pz_min) / dwn_n_mom_bins

    do p = 1, tile%npart_sp
      prtl_ux = tile%u(p); prtl_uy = tile%v(p); prtl_uz = tile%w(p)

      if (.not. masslessQ) then
        prtl_energy = sqrt(1.0 + prtl_ux**2 + prtl_uy**2 + prtl_uz**2)
      else
        prtl_energy = sqrt(prtl_ux**2 + prtl_uy**2 + prtl_uz**2)
      end if

      if ((prtl_energy .ge. dwn_energy_min) .and. (prtl_energy .lt. dwn_energy_max)) then
        ! rotation (not really random, because axis and angles are passed)
        call rotateRandomlyIn3D(prtl_ux, prtl_uy, prtl_uz, ax1, ax2, ang)

        if ((prtl_ux .ge. px_min) .and. (prtl_ux .lt. px_max) .and.&
          & (prtl_uy .ge. py_min) .and. (prtl_uy .lt. py_max) .and.&
          & (prtl_uz .ge. pz_min) .and. (prtl_uz .lt. pz_max)) then
          pi = INT(CEILING((prtl_ux - px_min) / del_ex))
          pj = INT(CEILING((prtl_uy - py_min) / del_ey))
          pk = INT(CEILING((prtl_uz - pz_min) / del_ez))
          pi = MIN(MAX(1, pi), dwn_n_mom_bins)
          pj = MIN(MAX(1, pj), dwn_n_mom_bins)
          pk = MIN(MAX(1, pk), dwn_n_mom_bins)
          momentum_bins(pi, pj, pk)%npart = momentum_bins(pi, pj, pk)%npart + 1
          momentum_bins(pi, pj, pk)%indices(momentum_bins(pi, pj, pk)%npart) = p
        end if
      end if
    end do
  end subroutine binParticlesOnTile_Cartesian
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
end module m_particlebinning


! #ifdef DEBUG
!   ! if ((pi .le. 0) .or. (pi .gt. dwn_n_mom_bins) .or.&
!   !   & (pj .le. 0) .or. (pj .gt. dwn_n_mom_bins) .or.&
!   !   & (pk .le. 0) .or. (pk .gt. dwn_n_mom_bins)) then
!   !   print *, pi, pj, pk
!   !   print *, prtl_ux, px_min, del_ex
!   !   print *, prtl_uy, py_min, del_ey
!   !   print *, prtl_uz, pz_min, del_ey
!   !   call throwError('ERROR: Wrong bin in `binParticlesOnTile_Cartesian()`.')
!   ! end if
!   ! if ((prtl_ux .le. momentum_bins(pi, pj, pk)%px_min) .or.&
!   !   & (prtl_ux .gt. momentum_bins(pi, pj, pk)%px_max) .or.&
!   !   & (prtl_uy .le. momentum_bins(pi, pj, pk)%py_min) .or.&
!   !   & (prtl_uy .gt. momentum_bins(pi, pj, pk)%py_max) .or.&
!   !   & (prtl_uz .le. momentum_bins(pi, pj, pk)%pz_min) .or.&
!   !   & (prtl_uz .gt. momentum_bins(pi, pj, pk)%pz_max)) then
!   !   print *, prtl_ux, momentum_bins(pi, pj, pk)%px_min, momentum_bins(pi, pj, pk)%px_max, px_min, px_max
!   !   print *, prtl_uy, momentum_bins(pi, pj, pk)%py_min, momentum_bins(pi, pj, pk)%py_max, py_min, py_max
!   !   print *, prtl_uz, momentum_bins(pi, pj, pk)%pz_min, momentum_bins(pi, pj, pk)%pz_max, pz_min, pz_max
!   !   call throwError('ERROR: Wrong bin by PXYZ in `binParticlesOnTile_Cartesian()`.')
!   ! end if
! #endif
