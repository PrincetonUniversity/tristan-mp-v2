module m_userfile
  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
  use m_exchangefields
  use m_exchangeparts
#ifdef USROUTPUT
  use m_writeusroutput
#endif
  use m_qednamespace

  implicit none

  !--- PRIVATE variables -----------------------------------------!

  integer, parameter, private :: n_modes = 8
  complex(kind=4), parameter, private :: ii = (0.0, 1.0)
  real, dimension(n_modes), private :: kx_ant, ky_ant, kz_ant
  complex(kind=4), dimension(n_modes), private :: b_k
  real, private :: omega0, gamma0, deltaB, T0, eph0, L0, ppc_ph
  integer, private :: esc_interval
  real, allocatable, private :: bx_ant(:, :, :), by_ant(:, :, :)
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution, advanceBext
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    real :: Pmax
    call getInput('problem', 'omega0', omega0, 0.01)
    call getInput('problem', 'gamma0', gamma0, 0.001)
    call getInput('problem', 'deltaB', deltaB, 0.8)
    call getInput('problem', 'T0', T0, 2e-3)
    call getInput('problem', 'eph0', eph0, T0)
    call getInput('problem', 'L0', L0, 1.0)
    call getInput('problem', 'ppc_ph', ppc_ph, ppc0)
    call getInput('problem', 'esc_interval', esc_interval, 4)
  end subroutine userReadInput

#ifdef PRTLPAYLOADS
  elemental subroutine usrSetPhPld(u0, v0, w0, over_e_temp, incr_pld1, incr_pld2, incr_pld3)
    !$omp declare simd(usrSetPhPld)
    real, intent(in) :: u0, v0, w0, over_e_temp
    real, intent(out) :: incr_pld1, incr_pld2, incr_pld3
    incr_pld1 = CC * u0 * over_e_temp
    incr_pld2 = CC * v0 * over_e_temp
    incr_pld3 = CC * w0 * over_e_temp
  end subroutine
  elemental subroutine usrSetElPld(q_over_m, u0, v0, w0, over_e_temp, ex0, &
                                   ey0, ez0, bx0, by0, bz0, incr_pld1, incr_pld2, incr_pld3)
    !$omp declare simd(usrSetElPld)
    real, intent(in) :: q_over_m, u0, v0, w0, over_e_temp, ex0, ey0, ez0, bx0, by0, bz0
    real, intent(out) :: incr_pld1, incr_pld2, incr_pld3
    incr_pld1 = CC * u0 * over_e_temp
    incr_pld2 = CC * v0 * over_e_temp
    incr_pld3 = CCINV * q_over_m * B_norm * (ex0 * bx0 + ey0 * by0 + ez0 * bz0) * &
               (u0 * bx0 + v0 * by0 + w0 * bz0) * over_e_temp / (bx0**2 + by0**2 + bz0**2)
  end subroutine
#endif

  function userSpatialDistribution(x_glob, y_glob, z_glob, &
                                   dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3

    return
  end function

  subroutine userDeallocate()
    implicit none
    ! external antenna fields
    if (allocated(bx_ant)) deallocate (bx_ant)
    if (allocated(by_ant)) deallocate (by_ant)
  end subroutine userDeallocate

  function userSLBload(x_glob, y_glob, z_glob, &
                       dummy1, dummy2, dummy3)
    real :: userSLBload
    ! global coordinates
    real, intent(in), optional :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional :: dummy1, dummy2, dummy3
    return
  end function

  ! random sampling of a Planck spectrum
  real function planck_sample(rnd)
    implicit none
    real, intent(in) :: rnd
    real :: prob, n
    n = 0.0
    prob = 0.0
    do while ((prob .lt. rnd) .and. (n .lt. 40.0))
      n = n + 1.0
      prob = prob + 1.0 / (1.20206 * n**3)
    end do
    planck_sample = -log(random(dseed) * random(dseed) * random(dseed) + 1e-16) / n
    return
  end function planck_sample

  subroutine advanceBext()
    implicit none
    integer :: mode, ierr
    complex(kind=4) :: u
    if (mpi_rank .eq. 0) then
      do mode = 1, n_modes / 2
        u = CMPLX(random(dseed) - 0.5, random(dseed) - 0.5, kind=4)
        b_k(mode) = b_k(mode) * cexp(-(ii * CMPLX(omega0, kind=4) + CMPLX(gamma0, kind=4)))
        b_k(mode) = b_k(mode) + CMPLX(deltaB * sqrt(12.0 * gamma0), kind=4) * u
        b_k(mode + n_modes / 2) = CONJG(b_k(mode))
      end do
    end if
    call MPI_BCAST(b_k, n_modes, MPI_COMPLEX, 0, MPI_COMM_WORLD, ierr)
  end subroutine advanceBext

  subroutine userInitParticles()
    implicit none
    real :: background_n
    type(region) :: back_region
    integer :: nph, n
    real :: xg, yg, zg, costh, phi, kx, ky, kz, eph
    integer :: s, ti, tj, tk, p
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution
    background_n = REAL(ppc0) * 0.5
    back_region % x_min = 0.0
    back_region % x_max = REAL(global_mesh % sx)
    back_region % y_min = 0.0
    back_region % y_max = REAL(global_mesh % sy)
    back_region % z_min = 0.0
    back_region % z_max = REAL(global_mesh % sz)
    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, background_n, T0)
    nph = INT(this_meshblock % ptr % sx * this_meshblock % ptr % sy * this_meshblock % ptr % sz * ppc_ph)
    do n = 1, nph
      xg = random(dseed) * REAL(this_meshblock % ptr % sx) + REAL(this_meshblock % ptr % x0)
      yg = random(dseed) * REAL(this_meshblock % ptr % sy) + REAL(this_meshblock % ptr % y0)
      zg = random(dseed) * REAL(this_meshblock % ptr % sz) + REAL(this_meshblock % ptr % z0)
      costh = 2.0 * random(dseed) - 1.0
      phi = 2.0 * M_PI * random(dseed)
      eph = eph0 * planck_sample(random(dseed))
      kx = eph * sqrt(1.0 - costh**2) * cos(phi)
      ky = eph * sqrt(1.0 - costh**2) * sin(phi)
      kz = eph * costh
      ! isotropic photons:
      call injectParticleGlobally(3, xg, yg, zg, kx, ky, kz)
    end do
#ifdef PRTLPAYLOADS
    do s = 1, nspec
      do tk = 1, species(s) % tile_nz
        do tj = 1, species(s) % tile_ny
          do ti = 1, species(s) % tile_nx
            do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
              ! initialize uniform displacements for photons at t=0
              if (s .eq. 3) then
                xg = (random(dseed) - 0.5) * REAL(global_mesh % sx)
                yg = (random(dseed) - 0.5) * REAL(global_mesh % sy)
                zg = (random(dseed) - 0.5) * REAL(global_mesh % sz)
                species(s) % prtl_tile(ti, tj, tk) % payload1(p) = xg
                species(s) % prtl_tile(ti, tj, tk) % payload2(p) = yg
                species(s) % prtl_tile(ti, tj, tk) % payload3(p) = zg
              else
                species(s) % prtl_tile(ti, tj, tk) % payload1(p) = 0.0
                species(s) % prtl_tile(ti, tj, tk) % payload2(p) = 0.0
                species(s) % prtl_tile(ti, tj, tk) % payload3(p) = 0.0
              end if
            end do
          end do
        end do
      end do
    end do
#endif
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: mode
    integer :: i1, i2, j1, j2, k1, k2
    real :: lx, ly, lz
    real :: ierr
#ifdef DEBUG
    integer :: i, j, k, im1, jm1, km1
    real :: check_x, check_y, check_z, db_sq, db_rms
#endif
    i1 = this_meshblock % ptr % i1
    i2 = this_meshblock % ptr % i2
    j1 = this_meshblock % ptr % j1
    j2 = this_meshblock % ptr % j2
    k1 = this_meshblock % ptr % k1
    k2 = this_meshblock % ptr % k2
    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = 0; by(:, :, :) = 0; bz(:, :, :) = 1.0 ! guide field in z direction
    jx(:, :, :) = 0; jy(:, :, :) = 0; jz(:, :, :) = 0
    call userDeallocate()
    allocate (bx_ant(i1:i2, j1:j2, k1:k2))
    allocate (by_ant(i1:i2, j1:j2, k1:k2))
    lx = REAL(global_mesh % sx)
    ly = REAL(global_mesh % sy)
    lz = REAL(global_mesh % sz)
    kx_ant(1) = 2.0 * M_PI / lx; ky_ant(1) = 0.0; kz_ant(1) = 2.0 * M_PI / lz
    kx_ant(2) = 2.0 * M_PI / lx; ky_ant(2) = 0.0; kz_ant(2) = -2.0 * M_PI / lz
    kx_ant(3) = 0.0; ky_ant(3) = 2.0 * M_PI / ly; kz_ant(3) = 2.0 * M_PI / lz
    kx_ant(4) = 0.0; ky_ant(4) = 2.0 * M_PI / ly; kz_ant(4) = -2.0 * M_PI / lz
    kx_ant(5) = -2.0 * M_PI / lx; ky_ant(5) = 0.0; kz_ant(5) = -2.0 * M_PI / lz
    kx_ant(6) = -2.0 * M_PI / lx; ky_ant(6) = 0.0; kz_ant(6) = 2.0 * M_PI / lz
    kx_ant(7) = 0.0; ky_ant(7) = -2.0 * M_PI / ly; kz_ant(7) = -2.0 * M_PI / lz
    kx_ant(8) = 0.0; ky_ant(8) = -2.0 * M_PI / ly; kz_ant(8) = 2.0 * M_PI / lz
    if (mpi_rank .eq. 0) then
      do mode = 1, n_modes / 2
        b_k(mode) = CMPLX(deltaB, kind=4) * cexp(ii * CMPLX(2.0 * M_PI * random(dseed), kind=4))
        b_k(mode + n_modes / 2) = CONJG(b_k(mode))
      end do
    end if
    call MPI_BCAST(b_k, n_modes, MPI_COMPLEX, 0, MPI_COMM_WORLD, ierr)
    call userCurrentDeposit(-1)
    ! initialize the fields to match the external current at t=0:
    bx(i1:i2, j1:j2, k1:k2) = bx_ant(i1:i2, j1:j2, k1:k2)
    by(i1:i2, j1:j2, k1:k2) = by_ant(i1:i2, j1:j2, k1:k2)
#ifdef DEBUG
    db_sq = 0.0
    do k = 0, this_meshblock % ptr % sz - 1
      km1 = k - 1
      do j = 0, this_meshblock % ptr % sy - 1
        jm1 = j - 1
        do i = 0, this_meshblock % ptr % sx - 1
          im1 = i - 1
          check_x = jx(i, j, k) + CORR * CC *&
                    & (by(i, j, km1) - by(i, j, k) - bz(i, jm1, k) + bz(i, j, k))
          check_y = jy(i, j, k) + CORR * CC *&
                    & (bz(im1, j, k) - bz(i, j, k) - bx(i, j, km1) + bx(i, j, k))
          check_z = jz(i, j, k) + CORR * CC *&
                    & (bx(i, jm1, k) - bx(i, j, k) - by(im1, j, k) + by(i, j, k))
          if ((abs(check_x) .gt. 1e-7) .or. &
              (abs(check_y) .gt. 1e-7) .or. &
              (abs(check_z) .gt. 1e-7)) then
            print *, "check_x, check_y, check_z = ", check_x, check_y, check_z
            call throwError('Curl(B) != J at initialization!')
          end if
          db_sq = db_sq + (bx(i, j, k)**2 + by(i, j, k)**2) / &
                  (REAL(global_mesh % sx) * REAL(global_mesh % sy) * REAL(global_mesh % sz))
        end do
      end do
    end do
    db_rms = 0.0
    call MPI_REDUCE(db_sq, db_rms, 1, default_mpi_real, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    db_rms = sqrt(db_rms)
    if (mpi_rank .eq. 0) then
      print *, "antenna amplitude at t=0, dB_rms / B0 = ", db_rms
    end if
#endif
  end subroutine userInitFields
  !............................................................!

  !--- driving ------------------------------------------------!
  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
    integer :: i, j, k, im1, jm1, km1, mode
    real :: x_, y_, z_shift, x_shift, y_shift
    real :: k_dot_r
    real, dimension(n_modes) :: coef_x, coef_y
    do mode = 1, n_modes
      coef_x(mode) = ky_ant(mode) / sqrt(REAL(n_modes) * (kx_ant(mode)**2 + ky_ant(mode)**2))
      coef_y(mode) = -kx_ant(mode) / sqrt(REAL(n_modes) * (kx_ant(mode)**2 + ky_ant(mode)**2))
    end do
    ! update the coefficients:
    if (step .gt. 0) call advanceBext()
    bx_ant(:, :, :) = 0.0; by_ant(:, :, :) = 0.0
    ! initialize the external field for this time step:
    do k = this_meshblock % ptr % k1, this_meshblock % ptr % k2
      z_shift = REAL(k + this_meshblock % ptr % z0) + 0.5
      do j = this_meshblock % ptr % j1, this_meshblock % ptr % j2
        y_ = REAL(j + this_meshblock % ptr % y0)
        y_shift = y_ + 0.5
        do i = this_meshblock % ptr % i1, this_meshblock % ptr % i2
          x_ = REAL(i + this_meshblock % ptr % x0)
          x_shift = x_ + 0.5
          do mode = 1, n_modes
            k_dot_r = kx_ant(mode) * x_ + ky_ant(mode) * y_shift + kz_ant(mode) * z_shift
            bx_ant(i, j, k) = bx_ant(i, j, k) + coef_x(mode) * REAL(ii * b_k(mode) * cexp(ii * CMPLX(k_dot_r, kind=4)))
            k_dot_r = kx_ant(mode) * x_shift + ky_ant(mode) * y_ + kz_ant(mode) * z_shift
            by_ant(i, j, k) = by_ant(i, j, k) + coef_y(mode) * REAL(ii * b_k(mode) * cexp(ii * CMPLX(k_dot_r, kind=4)))
          end do
        end do
      end do
    end do
    ! now setup the current as J = CC curl(B)
    do k = 0, this_meshblock % ptr % sz - 1
      km1 = k - 1
      do j = 0, this_meshblock % ptr % sy - 1
        jm1 = j - 1
        do i = 0, this_meshblock % ptr % sx - 1
          im1 = i - 1
          ! minus sign because this is how the current is deposited:
          jx(i, j, k) = jx(i, j, k) - CORR * CC * (by_ant(i, j, km1) - by_ant(i, j, k))
          jy(i, j, k) = jy(i, j, k) - CORR * CC * (bx_ant(i, j, k) - bx_ant(i, j, km1))
          jz(i, j, k) = jz(i, j, k) - CORR * CC * (bx_ant(i, jm1, k) - bx_ant(i, j, k) - by_ant(im1, j, k) + by_ant(i, j, k))
        end do
      end do
    end do
    call printDiag("userCurrentDeposit()", 2)
  end subroutine userCurrentDeposit

  subroutine userDriveParticles(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userDriveParticles

  subroutine userExternalFields(xp, yp, zp, &
                                ex_ext, ey_ext, ez_ext, &
                                bx_ext, by_ext, bz_ext)
    implicit none
    real, intent(in) :: xp, yp, zp
    real, intent(out) :: ex_ext, ey_ext, ez_ext
    real, intent(out) :: bx_ext, by_ext, bz_ext
    ! some functions of xp, yp, zp
    ex_ext = 0.0; ey_ext = 0.0; ez_ext = 0.0
    bx_ext = 0.0; by_ext = 0.0; bz_ext = 0.0
  end subroutine userExternalFields

  !............................................................!

  !--- boundaries ---------------------------------------------!
  subroutine userParticleBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
    integer :: s, ti, tj, tk, p, p_esc, s_esc
    real :: costh, phi, eph, xg, yg, zg
#ifdef PRTLPAYLOADS
    ! check escape condition every esc_interval steps:
    if (modulo(step, esc_interval) .eq. 0) then
      s = 3 ! photon species
      do tk = 1, species(s) % tile_nz
        do tj = 1, species(s) % tile_ny
          do ti = 1, species(s) % tile_nx
            s_esc = 4  ! the escaping photon species
            ! reset escaping species particle number:
            species(s_esc) % prtl_tile(ti, tj, tk) % npart_sp = 0
            do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
              ! NOTE: payloads update is done from inside the particle mover
              ! check if distance traveled matches escape condition:
              if (((abs(species(s) % prtl_tile(ti, tj, tk) % payload1(p)) .gt. (L0 * 0.5 * REAL(global_mesh % sx))) .or. &
                   (abs(species(s) % prtl_tile(ti, tj, tk) % payload2(p)) .gt. (L0 * 0.5 * REAL(global_mesh % sy))) .or. &
                   (abs(species(s) % prtl_tile(ti, tj, tk) % payload3(p)) .gt. (L0 * 0.5 * REAL(global_mesh % sz))))) then

                ! "move" photon to species=s_esc that represents the escaping population:
                species(s_esc) % prtl_tile(ti, tj, tk) % npart_sp = species(s_esc) % prtl_tile(ti, tj, tk) % npart_sp + 1
                species(s_esc) % cntr_sp = species(s_esc) % cntr_sp + 1
                p_esc = species(s_esc) % prtl_tile(ti, tj, tk) % npart_sp
                species(s_esc) % prtl_tile(ti, tj, tk) % xi(p_esc) = species(s) % prtl_tile(ti, tj, tk) % xi(p)
                species(s_esc) % prtl_tile(ti, tj, tk) % dx(p_esc) = species(s) % prtl_tile(ti, tj, tk) % dx(p)
                species(s_esc) % prtl_tile(ti, tj, tk) % yi(p_esc) = species(s) % prtl_tile(ti, tj, tk) % yi(p)
                species(s_esc) % prtl_tile(ti, tj, tk) % dy(p_esc) = species(s) % prtl_tile(ti, tj, tk) % dy(p)
                species(s_esc) % prtl_tile(ti, tj, tk) % zi(p_esc) = species(s) % prtl_tile(ti, tj, tk) % zi(p)
                species(s_esc) % prtl_tile(ti, tj, tk) % dz(p_esc) = species(s) % prtl_tile(ti, tj, tk) % dz(p)
                species(s_esc) % prtl_tile(ti, tj, tk) % u(p_esc) = species(s) % prtl_tile(ti, tj, tk) % u(p)
                species(s_esc) % prtl_tile(ti, tj, tk) % v(p_esc) = species(s) % prtl_tile(ti, tj, tk) % v(p)
                species(s_esc) % prtl_tile(ti, tj, tk) % w(p_esc) = species(s) % prtl_tile(ti, tj, tk) % w(p)
                species(s_esc) % prtl_tile(ti, tj, tk) % ind(p_esc) = species(s) % prtl_tile(ti, tj, tk) % ind(p)
                species(s_esc) % prtl_tile(ti, tj, tk) % proc(p_esc) = species(s) % prtl_tile(ti, tj, tk) % proc(p)
                species(s_esc) % prtl_tile(ti, tj, tk) % weight(p_esc) = species(s) % prtl_tile(ti, tj, tk) % weight(p)
                species(s_esc) % prtl_tile(ti, tj, tk) % payload1(p_esc) = species(s) % prtl_tile(ti, tj, tk) % payload1(p)
                species(s_esc) % prtl_tile(ti, tj, tk) % payload2(p_esc) = species(s) % prtl_tile(ti, tj, tk) % payload2(p)
                species(s_esc) % prtl_tile(ti, tj, tk) % payload3(p_esc) = species(s) % prtl_tile(ti, tj, tk) % payload3(p)

                ! particle escapes and is replaced with a new one from thermal bath:
                costh = 2.0 * random(dseed) - 1.0
                phi = 2.0 * M_PI * random(dseed)
                eph = eph0 * planck_sample(random(dseed))
                species(s) % prtl_tile(ti, tj, tk) % w(p) = eph * costh
                costh = sqrt(1.0 - costh**2)
                species(s) % prtl_tile(ti, tj, tk) % u(p) = eph * costh * cos(phi)
                species(s) % prtl_tile(ti, tj, tk) % v(p) = eph * costh * sin(phi)
                ! reset position:
                !xg = (random(dseed) - 0.5) * REAL(global_mesh % sx)
                !yg = (random(dseed) - 0.5) * REAL(global_mesh % sy)
                !zg = (random(dseed) - 0.5) * REAL(global_mesh % sz)
                !species(s)%prtl_tile(ti, tj, tk)%payload1(p) = xg
                !species(s)%prtl_tile(ti, tj, tk)%payload2(p) = yg
                !species(s)%prtl_tile(ti, tj, tk)%payload3(p) = zg
                species(s) % prtl_tile(ti, tj, tk) % payload1(p) = 0
                species(s) % prtl_tile(ti, tj, tk) % payload2(p) = 0
                species(s) % prtl_tile(ti, tj, tk) % payload3(p) = 0
                ! change ind and proc to imply this is a "new" particle:
                species(s) % prtl_tile(ti, tj, tk) % ind(p) = species(s) % cntr_sp
                species(s) % prtl_tile(ti, tj, tk) % proc(p) = mpi_rank
                species(s) % cntr_sp = species(s) % cntr_sp + 1
              end if
            end do
          end do
        end do
      end do
    end if
#endif
    call printDiag("userParticleBoundaryConditions()", 2)
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
    logical :: updateE_, updateB_

    if (present(updateE)) then
      updateE_ = updateE
    else
      updateE_ = .true.
    end if

    if (present(updateB)) then
      updateB_ = updateB
    else
      updateB_ = .true.
    end if
  end subroutine userFieldBoundaryConditions
  !............................................................!

  subroutine writeUsrRestart(rst_file)
    implicit none
    integer, intent(in) :: rst_file
    write (rst_file) kx_ant, ky_ant, kz_ant
    write (rst_file) b_k
  end subroutine writeUsrRestart

  subroutine readUsrRestart(rst_file)
    implicit none
    integer, intent(in) :: rst_file
    integer :: i1, i2, j1, j2, k1, k2
    character(len=STR_MAX) :: filename, mpichar
    i1 = this_meshblock % ptr % i1
    i2 = this_meshblock % ptr % i2
    j1 = this_meshblock % ptr % j1
    j2 = this_meshblock % ptr % j2
    k1 = this_meshblock % ptr % k1
    k2 = this_meshblock % ptr % k2
    read (rst_file) kx_ant, ky_ant, kz_ant
    read (rst_file) b_k
    ! make arrays for the external fields:
    call userDeallocate()
    allocate (bx_ant(i1:i2, j1:j2, k1:k2))
    allocate (by_ant(i1:i2, j1:j2, k1:k2))
  end subroutine readUsrRestart

  !--- user-specific output -----------------------------------!
#ifdef USROUTPUT
  subroutine userOutput(step)
    implicit none
    integer, optional, intent(in) :: step
    ! ...
  end subroutine userOutput

  logical function userExcludeParticles(s, ti, tj, tk, p)
    implicit none
    integer, intent(in) :: s, ti, tj, tk, p
    userExcludeParticles = .true.
  end function userExcludeParticles
#endif
  !............................................................!

end module m_userfile
