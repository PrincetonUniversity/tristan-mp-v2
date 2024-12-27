! Configuration for this userfile:
! ```
!   $ python configure.py -3d --user=user_psr_simple -slb --gca=2 --radiation=sync -absorb
! ```

module m_userfile
  use m_globalnamespace
  use m_qednamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
  use m_helpers
#ifdef USROUTPUT
  use m_writeusroutput
#endif
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real, private :: xc_g, yc_g, zc_g, psr_bstar, cooling_on
  real, private :: psr_angle, psr_period, psr_omega, psr_radius
  real, private :: inj_mult
  real, private :: shell_width, prtl_kick, rmin_dr
  real, private :: sigma_nGJ, nGJ, inj_dr
  real, private :: sigGJ_limiter, edotb_thr_closed
  real, private :: omegaB0, psr_rlc, alpha_inj

  real, private :: b_at_LC
  real, private :: gammarad_LC, gammarad_dummy
  real, private :: eph_at_LC, gammac_dummy
  integer, private :: vol_inj_start
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution, randomPointInSphericalShell
  private :: distance3D_sq, dumpParticlesHere
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'vol_inj_start', vol_inj_start)
    call getInput('problem', 'alpha_inj', alpha_inj)

    call getInput('problem', 'psr_radius', psr_radius)
    call getInput('problem', 'psr_angle', psr_angle)
    psr_angle = psr_angle * M_PI / 180.0
    call getInput('problem', 'psr_period', psr_period)
    psr_omega = 2.0 * M_PI / psr_period
    call getInput('problem', 'psr_bstar', psr_bstar)
    call getInput('problem', 'prtl_kick', prtl_kick, 0.0)
    if (prtl_kick .lt. 1.0) then
      ! interpreting as beta * gamma 
      prtl_kick = sqrt(1.0 + prtl_kick**2)
    end if

    call getInput('problem', 'inj_mult', inj_mult, 1.0)
    call getInput('problem', 'sigGJ_limiter', sigGJ_limiter, 1000.0)
    call getInput('problem', 'edotb_thr_closed', edotb_thr_closed, 0.002)

    call getInput('problem', 'cooling_on', cooling_on, 0.0)

    psr_rlc = CC / psr_omega
    omegaB0 = sqrt(sigma) * CC / c_omp

#ifdef RADIATION
    b_at_LC = 0.1 * psr_bstar * (psr_radius / psr_rlc)**3

    ! gamma_rad for the field near LC
    call getInput('problem', 'grad_LC', gammarad_LC, 0.5)
    gammarad_dummy = gammarad_LC * (b_at_LC)**0.5

    ! gamma of particle radiating me c^2 photon at LC
    call getInput('problem', 'eph_at_LC', eph_at_LC, 10.0)
    gammac_dummy = eph_at_LC * (b_at_LC)**0.5

    ! redefine `gamma_syn`, `emit_gamma_syn`
    cool_gamma_syn = gammarad_dummy
    emit_gamma_syn = gammac_dummy
#endif

    call getInput('problem', 'rmin_dr', rmin_dr, 1.0)
    call getInput('problem', 'shell_width', shell_width, 2.0)
    call getInput('problem', 'inj_dr', inj_dr, 1.0)

    xc_g = 0.5 * global_mesh % sx
    yc_g = 0.5 * global_mesh % sy
    zc_g = 0.5 * global_mesh % sz

    global_usr_variable_1 = psr_radius
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob, &
                                   dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3

    return
  end function

  function userSLBload(x_glob, y_glob, z_glob, &
                       dummy1, dummy2, dummy3)
    real :: userSLBload
    ! global coordinates
    real, intent(in), optional :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional :: dummy1, dummy2, dummy3
    real :: radius2, psrrad
    psrrad = global_usr_variable_1
    radius2 = (dummy1 * 0.5 - x_glob)**2 + (dummy2 * 0.5 - y_glob)**2 + (dummy3 * 0.5 - z_glob)**2 + 1.0
    userSLBload = psrrad**2 / radius2 + exp(-(dummy3 * 0.5 - z_glob)**2 / (psrrad * 0.5)**2)
    if (radius2 .lt. psrrad**2) then
      userSLBload = 1.0 / exp((psrrad**2 - radius2) / psrrad**2)
    end if
    return
  end function

  subroutine userInitParticles()
    implicit none
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    real :: bx0, by0, bz0, x_, y_, z_
    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = 0; by(:, :, :) = 0; bz(:, :, :) = 0
    jx(:, :, :) = 0; jy(:, :, :) = 0; jz(:, :, :) = 0

    do i = 0, this_meshblock % ptr % sx - 1
      i_glob = i + this_meshblock % ptr % x0
      do j = 0, this_meshblock % ptr % sy - 1
        j_glob = j + this_meshblock % ptr % y0
        do k = 0, this_meshblock % ptr % sz - 1
          k_glob = k + this_meshblock % ptr % z0

          x_ = REAL(i_glob); y_ = REAL(j_glob) + 0.5; z_ = REAL(k_glob) + 0.5
          call getBfield(0, 0.0, x_, y_, z_, bx0, by0, bz0)
          bx(i, j, k) = bx0

          x_ = REAL(i_glob) + 0.5; y_ = REAL(j_glob); z_ = REAL(k_glob) + 0.5
          call getBfield(0, 0.0, x_, y_, z_, bx0, by0, bz0)
          by(i, j, k) = by0

          x_ = REAL(i_glob) + 0.5; y_ = REAL(j_glob) + 0.5; z_ = REAL(k_glob)
          call getBfield(0, 0.0, x_, y_, z_, bx0, by0, bz0)
          bz(i, j, k) = bz0
        end do
      end do
    end do
  end subroutine userInitFields
  !............................................................!

  !--- driving ------------------------------------------------!
  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
    ! called after particles move and deposit ...
    ! ... and before the currents are added to the electric field
  end subroutine userCurrentDeposit

  subroutine userDriveParticles(step)
    implicit none
    integer, optional, intent(in) :: step
    ! ... dummy loop ...
    integer :: s, ti, tj, tk, p
    !real      :: gamma, nx, ny, nz, norm
    !real      :: gamma_max = 50.0
    !real      :: beta_max, u_max
    !beta_max = sqrt(1.0 - 1.0 / gamma_max**2)
    !u_max = gamma_max * beta_max
    !! limit particle gamma-factor
    !do s = 1, nspec
    !do ti = 1, species(s)%tile_nx
    !do tj = 1, species(s)%tile_ny
    !do tk = 1, species(s)%tile_nz
    !do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              !! for boris particles only
    !if (species(s)%prtl_tile(ti, tj, tk)%proc(p) .lt. mpi_size) then
    !gamma = sqrt(1.0 + species(s)%prtl_tile(ti, tj, tk)%u(p)**2 +&
    !& species(s)%prtl_tile(ti, tj, tk)%v(p)**2 +&
    !& species(s)%prtl_tile(ti, tj, tk)%w(p)**2)
    !if (gamma .gt. gamma_max) then
    !norm = species(s)%prtl_tile(ti, tj, tk)%u(p)**2 +&
    !& species(s)%prtl_tile(ti, tj, tk)%v(p)**2 +&
    !& species(s)%prtl_tile(ti, tj, tk)%w(p)**2
    !species(s)%prtl_tile(ti, tj, tk)%u(p) = u_max *&
    !& (species(s)%prtl_tile(ti, tj, tk)%u(p) / norm)
    !species(s)%prtl_tile(ti, tj, tk)%v(p) = u_max *&
    !& (species(s)%prtl_tile(ti, tj, tk)%v(p) / norm)
    !species(s)%prtl_tile(ti, tj, tk)%w(p) = u_max *&
    !& (species(s)%prtl_tile(ti, tj, tk)%w(p) / norm)
    !end if
    !end if
    !end do
    !end do
    !end do
    !end do
    !end do
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
  subroutine dumpParticlesHere(xi, yi, zi, dx, dy, dz, weight)
    implicit none
    integer(kind=2), intent(in) :: xi, yi, zi
    real, intent(in) :: dx, dy, dz, weight
    call createParticle(1, xi, yi, zi, dx, dy, dz, 0.0, 0.0, 0.0, weight=weight)
    call createParticle(2, xi, yi, zi, dx, dy, dz, 0.0, 0.0, 0.0, weight=weight)
  end subroutine dumpParticlesHere

  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  ! VOLUME INJECTOR
  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  !subroutine userParticleBoundaryConditions(step)
  !implicit none
  !integer, optional, intent(in) :: step
  !integer                       :: s, ti, tj, tk, p
  !real                          :: x_g, y_g, z_g, r_g
  !integer                       :: n_part, n, sign
  !real                          :: x_loc, y_loc, z_loc, dx, dy, dz
  !integer(kind=2)               :: xi, yi, zi, xi_eb, yi_eb, zi_eb
  !real                          :: x_glob, y_glob, z_glob, weight, ppc, dens, sig, rmax, rmin
  !real                          :: e_dot_b, b_sqr, delta_er, bx0, by0, bz0, ex0, ey0, ez0
  !real                          :: u_, v_, w_, nx, ny, nz, rr, vx, vy, vz, gamma, rmax_sph, smax_car
  !real                          :: dens_GJ, e_b_scale, j_dot_b, density, jx0, jy0, jz0
  !real                          :: bz_at_surface, mu_x, mu_y, mu_z, mu_dot_n, phase, dummy_, cos_theta
  !logical                       :: dummy_flag
  !#ifdef GCA
  !real                          :: vE_x, vE_y, vE_z
  !real                          :: e0_SQR, b0_SQR
  !#endif
  !real        :: dr_i, rmin_inj, rmax_inj
  !integer     :: ri
  !dr_i = 2.0

  !#ifdef RADIATION
      !! disable cooling until certain point
  !if (step .lt. cooling_on * psr_period) then
  !species(1)%cool_sp = .false.
  !species(2)%cool_sp = .false.
  !else
  !species(1)%cool_sp = .true.
  !species(2)%cool_sp = .true.
  !end if
  !#endif

  !if (sigGJ_limiter .ne. 0) then
      !! compute number density and write to `lg_arr`
  !call computeDensity(1, reset=.true., ds=0, charge=.false.)
  !call computeDensity(2, reset=.false., ds=0, charge=.false.)
  !end if

  !nGJ = 2 * psr_omega * B_norm * psr_bstar / (CC * abs(unit_ch))
  !ppc = 0.5 * ppc0
  !weight = inj_mult * nGJ / ppc
  !do ri = 0, INT(shell_width / dr_i)
      !! iterate over spherical shells
  !rmin_inj = psr_radius + inj_dr + REAL(ri) * dr_i
  !rmax_inj = psr_radius + inj_dr + REAL(ri + 1) * dr_i
  !n_part = INT(2.0 * (4.0 * M_PI / 3.0) * (rmax_inj**3 - rmin_inj**3) * ppc)
  !do n = 1, n_part
        !! for each spherical shell iterate over the given number of particles to inject
  !call randomPointInSphericalShell(rmin_inj, rmax_inj, x_glob, y_glob, z_glob)
  !rr = sqrt(x_glob**2 + y_glob**2 + z_glob**2)
  !nx = x_glob / rr
  !ny = y_glob / rr
  !nz = z_glob / rr
  !cos_theta = (nz)
  !x_glob = x_glob + xc_g
  !y_glob = y_glob + yc_g
  !z_glob = z_glob + zc_g
  !call globalToLocalCoords(x_glob, y_glob, z_glob, x_loc, y_loc, z_loc, containedQ=dummy_flag)
        !! if particle is within the current MPI meshblock
  !if (dummy_flag) then
          !! open zone injection
  !call localToCellBasedCoords(x_loc, y_loc, z_loc, xi, yi, zi, dx, dy, dz)
  !dummy_flag = .true.
          !! limiter on sigma
  !if (sigGJ_limiter .ne. 0) then
  !density = lg_arr(xi, yi, zi)
  !b_sqr = bx(xi, yi, zi)**2 + by(xi, yi, zi)**2 + bz(xi, yi, zi)**2
  !if (density .gt. 0) then
  !sig = b_sqr * sigma * ppc0 / density
  !else
  !sig = sigma
  !end if
  !dummy_flag = (sig .gt. sigGJ_limiter * (psr_radius / rr)**3)
  !end if
  !if (dummy_flag) then
  !call dumpParticlesHere(xi, yi, zi, dx, dy, dz, weight)
  !end if
  !end if
  !end do ! npart
  !end do ! ri

    !! max radius for spherical absorption and absorbing layer size in cartesian absorption
  !rmax_sph = 0.5 * MIN(global_mesh%sx, global_mesh%sy, global_mesh%sz) - ds_abs + 4.0
  !smax_car = ds_abs / 4.0
    !! remove particles falling into the star
  !do s = 1, nspec
  !do ti = 1, species(s)%tile_nx
  !do tj = 1, species(s)%tile_ny
  !do tk = 1, species(s)%tile_nz
  !do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
  !x_g = REAL(species(s)%prtl_tile(ti, tj, tk)%xi(p) + this_meshblock%ptr%x0)&
  !& + species(s)%prtl_tile(ti, tj, tk)%dx(p)
  !y_g = REAL(species(s)%prtl_tile(ti, tj, tk)%yi(p) + this_meshblock%ptr%y0)&
  !& + species(s)%prtl_tile(ti, tj, tk)%dy(p)
  !z_g = REAL(species(s)%prtl_tile(ti, tj, tk)%zi(p) + this_meshblock%ptr%z0)&
  !& + species(s)%prtl_tile(ti, tj, tk)%dz(p)
  !r_g = sqrt((x_g - xc_g)**2 + (y_g - yc_g)**2 + (z_g - zc_g)**2)
  !if (&
  !& (r_g .lt. (psr_radius - rmin_dr)) .or.&
  !& ((boundary_x .eq. 2) .and. (r_g .gt. rmax_sph)) .or.&
  !& ((boundary_x .eq. 0) .and. ((x_g .lt. smax_car) .or.&
  !& (x_g .gt. global_mesh%sx - smax_car) .or.&
  !& (y_g .lt. smax_car) .or.&
  !& (y_g .gt. global_mesh%sy - smax_car) .or.&
  !& (z_g .lt. smax_car) .or.&
  !& (z_g .gt. global_mesh%sz - smax_car))&
  !& )&
  !& ) then
  !species(s)%prtl_tile(ti, tj, tk)%proc(p) = -1
  !end if
  !end do
  !end do
  !end do
  !end do
  !end do
  !end subroutine userParticleBoundaryConditions

  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  ! SURFACE INJECTOR EVERYWHERE
  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  !subroutine userParticleBoundaryConditions(step)
  !implicit none
  !integer, optional, intent(in) :: step
  !integer                       :: s, ti, tj, tk, p
  !real                          :: x_g, y_g, z_g, r_g
  !integer                       :: n_part, n, sign
  !real                          :: x_loc, y_loc, z_loc, dx, dy, dz
  !integer(kind=2)               :: xi, yi, zi, xi_eb, yi_eb, zi_eb
  !real                          :: x_glob, y_glob, z_glob, weight, ppc, dens, sig, rmax, rmin
  !real                          :: e_dot_b, b_sqr, delta_er, bx0, by0, bz0, ex0, ey0, ez0
  !real                          :: u_, v_, w_, nx, ny, nz, rr, vx, vy, vz, gamma, rmax_sph, smax_car
  !real                          :: dens_GJ, e_b_scale, j_dot_b, density, jx0, jy0, jz0
  !real                          :: bz_at_surface, mu_x, mu_y, mu_z, mu_dot_n, phase, dummy_
  !logical                       :: dummy_flag
  !#ifdef GCA
  !real                          :: vE_x, vE_y, vE_z
  !real                          :: e0_SQR, b0_SQR
  !#endif

  !phase = step * psr_omega
  !mu_x = sin(psr_angle) * cos(phase)
  !mu_y = sin(psr_angle) * sin(phase)
  !mu_z = cos(psr_angle)

  !#ifdef RADIATION
      !! disable cooling until certain point
  !if (step .lt. cooling_on * psr_period) then
  !species(1)%cool_sp = .false.
  !species(2)%cool_sp = .false.
  !else
  !species(1)%cool_sp = .true.
  !species(2)%cool_sp = .true.
  !end if
  !#endif

  !nGJ = 2 * psr_omega * B_norm * psr_bstar / (CC * abs(unit_ch))

  !ppc = 0.5 * ppc0
  !if (sigGJ_limiter .ne. 0) then
      !! compute number density and write to `lg_arr`
  !call computeDensity(1, reset=.true., ds=0, charge=.false.)
  !call computeDensity(2, reset=.false., ds=0, charge=.false.)
  !end if
  !n_part = INT(2.0 * (4.0 * M_PI / 3.0) * ((psr_radius + shell_width + inj_dr)**3 - (psr_radius + inj_dr)**3) * ppc)
  !do n = 1, n_part
  !call randomPointInSphericalShell(psr_radius + inj_dr, psr_radius + inj_dr + shell_width, x_glob, y_glob, z_glob)
  !rr = sqrt(x_glob**2 + y_glob**2 + z_glob**2)
  !nx = x_glob / rr
  !ny = y_glob / rr
  !nz = z_glob / rr

  !mu_dot_n = mu_x * nx + mu_y * ny + mu_z * nz

  !x_glob = x_glob + xc_g
  !y_glob = y_glob + yc_g
  !z_glob = z_glob + zc_g
  !call globalToLocalCoords(x_glob, y_glob, z_glob, x_loc, y_loc, z_loc, containedQ=dummy_flag)
      !! if particle is within the current MPI meshblock
  !if (dummy_flag) then
  !call localToCellBasedCoords(x_loc, y_loc, z_loc, xi, yi, zi, dx, dy, dz)
        !! limiter on sigma
  !dummy_flag = .true.
  !if (sigGJ_limiter .ne. 0) then
  !density = lg_arr(xi, yi, zi)
  !b_sqr = bx(xi, yi, zi)**2 + by(xi, yi, zi)**2 + bz(xi, yi, zi)**2
  !if (density .gt. 0) then
  !sig = b_sqr * sigma * ppc0 / density
  !else
  !sig = sigma
  !end if
  !dummy_flag = (sig .gt. sigGJ_limiter * (psr_radius / rr)**3)
  !end if ! sigma limiter enabled
  !if (dummy_flag) then
          !! initialize at rest
  !call interpFromEdges(dx, dy, dz, xi, yi, zi, ex, ey, ez, ex0, ey0, ez0)
  !call interpFromFaces(dx, dy, dz, xi, yi, zi, bx, by, bz, bx0, by0, bz0)
  !b_sqr = sqrt(bx0**2 + by0**2 + bz0**2)
  !sign = 1
  !if (bx0 * nx + by0 * ny + bz0 * nz .lt. 0) then
  !sign = -1
  !bx0 = -bx0; by0 = -by0; bz0 = -bz0
  !end if
  !nx = bx0 / b_sqr
  !ny = by0 / b_sqr
  !nz = bz0 / b_sqr
  !dummy_ = abs(ex0 * bx0 + ey0 * by0 + ez0 * bz0) / (bx0**2 + by0**2 + bz0**2)
  !if (dummy_ .gt. edotb_thr_closed) then
  !weight = inj_mult * dummy_ * nGJ / ppc
  !if (prtl_kick .eq. 0) then
  !call createParticle(1, xi, yi, zi, dx, dy, dz, 0.0, 0.0, 0.0, weight=weight)
  !call createParticle(2, xi, yi, zi, dx, dy, dz, 0.0, 0.0, 0.0, weight=weight)
  !else
  !#ifndef GCA
  !u_ = nx * prtl_kick
  !v_ = ny * prtl_kick
  !w_ = nz * prtl_kick
  !call createParticle(1, xi, yi, zi, dx, dy, dz, u_, v_, w_, weight=weight)
  !call createParticle(2, xi, yi, zi, dx, dy, dz, u_, v_, w_, weight=weight)
  !#else
                !! kick along ExB:
  !ex0 = ex(xi, yi, zi); ey0 = ey(xi, yi, zi); ez0 = ez(xi, yi, zi)
  !bx0 = bx(xi, yi, zi); by0 = by(xi, yi, zi); bz0 = bz(xi, yi, zi)
  !b0_SQR = bx0**2 + by0**2 + bz0**2
  !e0_SQR = ex0**2 + ey0**2 + ez0**2

  !dummy_ = 1.0 / (b0_SQR + TINYFLD)

  !vE_x = (bz0 * ey0 - by0 * ez0) * dummy_
  !vE_y = (-bz0 * ex0 + bx0 * ez0) * dummy_
  !vE_z = (by0 * ex0 - bx0 * ey0) * dummy_

  !dummy_ = 1.0 / sqrt(abs(1.0 - vE_x**2 - vE_y**2 - vE_z**2) + TINYFLD)

  !u_ = vE_x * dummy_
  !v_ = vE_y * dummy_
  !w_ = vE_z * dummy_

  !dummy_ = sqrt(prtl_kick**2 - 1.0)

  !u_ = u_ + nx * dummy_
  !v_ = v_ + ny * dummy_
  !w_ = w_ + nz * dummy_

  !call createParticleFromAttributes(1, xi=xi, yi=yi, zi=zi, dx=dx, dy=dy, dz=dz,&
  !& xi_past=xi, yi_past=yi, zi_past=zi,&
  !& dx_past=dx, dy_past=dy, dz_past=dz,&
  !& u=u_, v=v_, w=w_,&
  !& u_eff=u_, v_eff=v_, w_eff=w_, u_par=sign*dummy_, u_perp=0.0,&
  !& ind=species(1)%cntr_sp, proc=mpi_rank + 2 * mpi_size, weight=weight)
  !species(1)%cntr_sp = species(1)%cntr_sp + 1
  !call createParticleFromAttributes(2, xi=xi, yi=yi, zi=zi, dx=dx, dy=dy, dz=dz,&
  !& xi_past=xi, yi_past=yi, zi_past=zi,&
  !& dx_past=dx, dy_past=dy, dz_past=dz,&
  !& u=u_, v=v_, w=w_,&
  !& u_eff=u_, v_eff=v_, w_eff=w_, u_par=sign*dummy_, u_perp=0.0,&
  !& ind=species(2)%cntr_sp, proc=mpi_rank + 2 * mpi_size, weight=weight)
  !species(2)%cntr_sp = species(2)%cntr_sp + 1
  !#endif
  !end if
  !end if ! E.B limiter
  !end if ! SIGMA limiter
  !end if ! current MPI block
  !end do

    !! max radius for spherical absorption and absorbing layer size in cartesian absorption
  !rmax_sph = 0.5 * MIN(global_mesh%sx, global_mesh%sy, global_mesh%sz) - ds_abs + 4.0
  !smax_car = ds_abs / 4.0
    !! remove particles falling into the star
  !do s = 1, nspec
  !do ti = 1, species(s)%tile_nx
  !do tj = 1, species(s)%tile_ny
  !do tk = 1, species(s)%tile_nz
  !do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
  !x_g = REAL(species(s)%prtl_tile(ti, tj, tk)%xi(p) + this_meshblock%ptr%x0)&
  !& + species(s)%prtl_tile(ti, tj, tk)%dx(p)
  !y_g = REAL(species(s)%prtl_tile(ti, tj, tk)%yi(p) + this_meshblock%ptr%y0)&
  !& + species(s)%prtl_tile(ti, tj, tk)%dy(p)
  !z_g = REAL(species(s)%prtl_tile(ti, tj, tk)%zi(p) + this_meshblock%ptr%z0)&
  !& + species(s)%prtl_tile(ti, tj, tk)%dz(p)
  !r_g = sqrt((x_g - xc_g)**2 + (y_g - yc_g)**2 + (z_g - zc_g)**2)
  !if (&
  !& (r_g .lt. (psr_radius - rmin_dr)) .or.&
  !& ((boundary_x .eq. 2) .and. (r_g .gt. rmax_sph)) .or.&
  !& ((boundary_x .eq. 0) .and. ((x_g .lt. smax_car) .or.&
  !& (x_g .gt. global_mesh%sx - smax_car) .or.&
  !& (y_g .lt. smax_car) .or.&
  !& (y_g .gt. global_mesh%sy - smax_car) .or.&
  !& (z_g .lt. smax_car) .or.&
  !& (z_g .gt. global_mesh%sz - smax_car))&
  !& )&
  !& ) then
  !species(s)%prtl_tile(ti, tj, tk)%proc(p) = -1
  !end if
  !end do
  !end do
  !end do
  !end do
  !end do
  !end subroutine userParticleBoundaryConditions

  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  ! SURFACE INJECTOR ON POLAR CAP + CLOSED ZONE INJECTOR
  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  subroutine userParticleBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
    integer :: s, ti, tj, tk, p, i, j, k, i_glob, j_glob, k_glob
    real :: x_g, y_g, z_g, r_g
    integer :: n_part, n, sign, fi
    real :: x_loc, y_loc, z_loc, dx, dy, dz
    integer(kind=2) :: xi, yi, zi, xi_eb, yi_eb, zi_eb
    real :: x_glob, y_glob, z_glob, weight, ppc, dens, sig, rmax, rmin
    real :: e_dot_b, b_sqr, delta_er, bx0, by0, bz0, ex0, ey0, ez0
    real :: u_, v_, w_, nx, ny, nz, rr, vx, vy, vz, gamma, rmax_sph, smax_car
    real :: dens_GJ, e_b_scale, j_dot_b, density, jx0, jy0, jz0
    real :: bz_at_surface, mu_x, mu_y, mu_z, mu_dot_n, phase, dummy_
    logical :: dummy_flag
    real :: cos_theta

    real :: dr_i, rmin_inj, rmax_inj

    phase = step * psr_omega
    mu_x = sin(psr_angle) * cos(phase)
    mu_y = sin(psr_angle) * sin(phase)
    mu_z = cos(psr_angle)

#ifdef RADIATION
    ! disable cooling until certain point
    if (step .lt. cooling_on * psr_period) then
      species(1) % cool_sp = .false.
      species(2) % cool_sp = .false.
    else
      species(1) % cool_sp = .true.
      species(2) % cool_sp = .true.
    end if
#endif

    nGJ = 2 * psr_omega * B_norm * psr_bstar / (CC * abs(unit_ch))

    ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    ppc = 0.5 * ppc0
    ! if (sigGJ_limiter .ne. 0) then
    ! compute number density and write to `lg_arr`
    call computeDensity(1, reset=.true., ds=0, charge=.false.)
    call computeDensity(2, reset=.false., ds=0, charge=.false.)
    ! end if
    n_part = INT(2.0 * (4.0 * M_PI / 3.0) * ((psr_radius + shell_width + inj_dr)**3 - (psr_radius + inj_dr)**3) * ppc)
    do n = 1, n_part
      call randomPointInSphericalShell(psr_radius + inj_dr, psr_radius + inj_dr + shell_width, x_glob, y_glob, z_glob)
      rr = sqrt(x_glob**2 + y_glob**2 + z_glob**2)
      nx = x_glob / rr
      ny = y_glob / rr
      nz = z_glob / rr

      cos_theta = mu_x * nx + mu_y * ny + mu_z * nz

      bz_at_surface = (3.0 * nz * cos_theta - mu_z)

      dummy_flag = (((random(dseed) * 2.0) .lt. abs(bz_at_surface)) .and. &
                    (abs(cos_theta) .gt. 0.5))
      ! injection in the open zone
      if (dummy_flag) then
        x_glob = x_glob + xc_g
        y_glob = y_glob + yc_g
        z_glob = z_glob + zc_g
        call globalToLocalCoords(x_glob, y_glob, z_glob, x_loc, y_loc, z_loc, containedQ=dummy_flag)
        ! if particle is within the current MPI meshblock
        if (dummy_flag) then
          call localToCellBasedCoords(x_loc, y_loc, z_loc, xi, yi, zi, dx, dy, dz)
          dummy_flag = .true.

          ! limiter on sigma
          if (sigGJ_limiter .ne. 0) then
            density = lg_arr(xi, yi, zi)
            b_sqr = bx(xi, yi, zi)**2 + by(xi, yi, zi)**2 + bz(xi, yi, zi)**2
            if (density .gt. 0) then
              sig = b_sqr * sigma * ppc0 / density
            else
              sig = sigma
            end if
            dummy_flag = (sig .gt. sigGJ_limiter)
          end if

          if (dummy_flag) then
            ! kick along local b-field:
            bx0 = bx(xi, yi, zi); by0 = by(xi, yi, zi); bz0 = bz(xi, yi, zi)
            b_sqr = sqrt(bx0**2 + by0**2 + bz0**2)
            sign = 1
            if (bx0 * nx + by0 * ny + bz0 * nz .lt. 0) then
              sign = -1
              bx0 = -bx0; by0 = -by0; bz0 = -bz0
            end if
            nx = bx0 / b_sqr
            ny = by0 / b_sqr
            nz = bz0 / b_sqr

            weight = inj_mult * nGJ / ppc

            u_ = nx * prtl_kick
            v_ = ny * prtl_kick
            w_ = nz * prtl_kick
            call createParticle(1, xi, yi, zi, dx, dy, dz, u_, v_, w_, weight=weight)
            call createParticle(2, xi, yi, zi, dx, dy, dz, u_, v_, w_, weight=weight)
          end if
        end if
      end if ! open zone injection

      ! ! injection in the closed zone
      ! if (abs(cos_theta) .lt. 0.5) then
      !   x_glob = x_glob + xc_g
      !   y_glob = y_glob + yc_g
      !   z_glob = z_glob + zc_g
      !   call globalToLocalCoords(x_glob, y_glob, z_glob, x_loc, y_loc, z_loc, containedQ=dummy_flag)
      !   ! if particle is within the current MPI meshblock
      !   if (dummy_flag) then
      !     call localToCellBasedCoords(x_loc, y_loc, z_loc, xi, yi, zi, dx, dy, dz)
      !     ! initialize at rest
      !     call interpFromEdges(dx, dy, dz, xi, yi, zi, ex, ey, ez, ex0, ey0, ez0)
      !     call interpFromFaces(dx, dy, dz, xi, yi, zi, bx, by, bz, bx0, by0, bz0)
      !     dummy_ = abs(ex0 * bx0 + ey0 * by0 + ez0 * bz0) / (bx0**2 + by0**2 + bz0**2)
      !     if (dummy_ .gt. edotb_thr_closed) then
      !       weight = dummy_ * nGJ / ppc
      !       call createParticle(1, xi, yi, zi, dx, dy, dz, 0.0, 0.0, 0.0, weight=weight)
      !       call createParticle(2, xi, yi, zi, dx, dy, dz, 0.0, 0.0, 0.0, weight=weight)
      !     end if ! E.B limiter
      !   end if ! current MPI block
      ! end if ! closed zone

    end do

    ! volume injection
    nGJ = 2 * psr_omega * B_norm * psr_bstar / (CC * abs(unit_ch))
    ppc = 0.5 * ppc0
    weight = inj_mult * nGJ / ppc
    rmin_inj = 2.0 * psr_radius
    rmax_inj = 4.0 * psr_radius
    n_part = INT(2.0 * (4.0 * M_PI / 3.0) * (rmax_inj**3 - rmin_inj**3) * ppc) / 1000
    do n = 1, n_part
      ! for each spherical shell iterate over the given number of particles to inject
      call randomPointInSphericalShell(rmin_inj, rmax_inj, x_glob, y_glob, z_glob)
      rr = sqrt(x_glob**2 + y_glob**2 + z_glob**2)
      nx = x_glob / rr
      ny = y_glob / rr
      nz = z_glob / rr
      cos_theta = (nz)
      if (abs(cos_theta) .lt. 0.6) then
        x_glob = x_glob + xc_g
        y_glob = y_glob + yc_g
        z_glob = z_glob + zc_g
        call globalToLocalCoords(x_glob, y_glob, z_glob, x_loc, y_loc, z_loc, containedQ=dummy_flag)
        ! if particle is within the current MPI meshblock
        if (dummy_flag) then
          ! open zone injection
          call localToCellBasedCoords(x_loc, y_loc, z_loc, xi, yi, zi, dx, dy, dz)
          density = lg_arr(xi, yi, zi)
          if (density .lt. 20.0 * nGJ * (psr_radius / rr)**2) then
            call dumpParticlesHere(xi, yi, zi, dx, dy, dz, weight)
          end if
        end if
      end if
    end do ! npart

    ! max radius for spherical absorption and absorbing layer size in cartesian absorption
    rmax_sph = 0.5 * MIN(global_mesh % sx, global_mesh % sy, global_mesh % sz) - ds_abs + 4.0
    smax_car = ds_abs / 4.0

    !if (step .ge. vol_inj_start) then
      !! inject in low density regions
    !call computeDensity(1, reset=.true., ds=0, charge=.false.)
    !call computeDensity(2, reset=.false., ds=0, charge=.false.)
    !do i = 0, this_meshblock%ptr%sx - 1
    !i_glob = i + this_meshblock%ptr%x0
    !do j = 0, this_meshblock%ptr%sy - 1
    !j_glob = j + this_meshblock%ptr%y0
    !do k = 0, this_meshblock%ptr%sz - 1
    !k_glob = k + this_meshblock%ptr%z0
    !r_g = sqrt(REAL(i_glob - xc_g)**2 + REAL(j_glob - yc_g)**2 + REAL(k_glob - zc_g)**2) + 0.0001
    !density = lg_arr(i, j, k)
    !if ((r_g .gt. 0.5 * psr_radius) .and.&
    !& (density .lt. 5.0 * (nGJ * (psr_radius / r_g)**3)) .and.&
    !& (r_g .lt. 0.5 * rmax_sph)) then
    !ppc = 0.5 * ppc0
    !weight = alpha_inj * (5.0 * nGJ * (psr_radius / r_g)**3) / ppc
    !ppc = ppc * tanh((r_g - psr_radius) / (0.5 * psr_radius))
    !do while (ppc .gt. 0)
    !if (random(dseed) .lt. ppc) then
    !xi = INT(i, 2); yi = INT(j, 2); zi = INT(k, 2)
    !dx = random(dseed); dy = random(dseed); dz = random(dseed)
    !call createParticle(1, xi, yi, zi, dx, dy, dz, 0.0, 0.0, 0.0, weight=weight)
    !call createParticle(2, xi, yi, zi, dx, dy, dz, 0.0, 0.0, 0.0, weight=weight)
                  !!end if ! E.B limiter
    !end if
    !ppc = ppc - 1.0
    !end do
    !end if
    !end do
    !end do
    !end do
    !end if

    ! remove particles falling into the star
    do s = 1, nspec
      do ti = 1, species(s) % tile_nx
        do tj = 1, species(s) % tile_ny
          do tk = 1, species(s) % tile_nz
            do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
              x_g = REAL(species(s) % prtl_tile(ti, tj, tk) % xi(p) + this_meshblock % ptr % x0) &
                    + species(s) % prtl_tile(ti, tj, tk) % dx(p)
              y_g = REAL(species(s) % prtl_tile(ti, tj, tk) % yi(p) + this_meshblock % ptr % y0) &
                    + species(s) % prtl_tile(ti, tj, tk) % dy(p)
              z_g = REAL(species(s) % prtl_tile(ti, tj, tk) % zi(p) + this_meshblock % ptr % z0) &
                    + species(s) % prtl_tile(ti, tj, tk) % dz(p)
              r_g = sqrt((x_g - xc_g)**2 + (y_g - yc_g)**2 + (z_g - zc_g)**2)
              if ( &
                (r_g .lt. (psr_radius - rmin_dr)) .or. &
                ((boundary_x .eq. 2) .and. (r_g .gt. rmax_sph)) .or. &
                ((boundary_x .eq. 0) .and. ((x_g .lt. smax_car) .or. &
                                            (x_g .gt. global_mesh % sx - smax_car) .or. &
                                            (y_g .lt. smax_car) .or. &
                                            (y_g .gt. global_mesh % sy - smax_car) .or. &
                                            (z_g .lt. smax_car) .or. &
                                            (z_g .gt. global_mesh % sz - smax_car)) &
                 ) &
                ) then
                species(s) % prtl_tile(ti, tj, tk) % proc(p) = -1
              end if
            end do
          end do
        end do
      end do
    end do
  end subroutine userParticleBoundaryConditions

  !! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  !! DISK-DOME INJECTOR
  !! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  !subroutine userParticleBoundaryConditions(step)
  !implicit none
  !integer, optional, intent(in) :: step
  !integer                       :: s, ti, tj, tk, p
  !real                          :: x_g, y_g, z_g, r_g
  !integer                       :: n_part, n, sign
  !real                          :: x_loc, y_loc, z_loc, dx, dy, dz
  !integer(kind=2)               :: xi, yi, zi, xi_eb, yi_eb, zi_eb
  !real                          :: x_glob, y_glob, z_glob, weight, ppc, dens, sig, rmax, rmin
  !real                          :: e_dot_b, b_sqr, delta_er, bx0, by0, bz0, ex0, ey0, ez0
  !real                          :: u_, v_, w_, nx, ny, nz, rr, vx, vy, vz, gamma, rmax_sph, smax_car
  !real                          :: dens_GJ, e_b_scale, j_dot_b, density, jx0, jy0, jz0
  !real                          :: bz_at_surface, mu_x, mu_y, mu_z, mu_dot_n, phase
  !logical                       :: dummy_flag
  !real                          :: dummy_

  !phase = step * psr_omega
  !mu_x = sin(psr_angle) * cos(phase)
  !mu_y = sin(psr_angle) * sin(phase)
  !mu_z = cos(psr_angle)
  !nGJ = 2 * psr_omega * B_norm * psr_bstar / (CC * abs(unit_ch))

    !! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  !ppc = 0.5 * ppc0
  !n_part = INT(2.0 * (4.0 * M_PI / 3.0) * ((psr_radius + shell_width + inj_dr)**3 - (psr_radius + inj_dr)**3) * ppc)
  !do n = 1, n_part
  !call randomPointInSphericalShell(psr_radius + inj_dr, psr_radius + inj_dr + shell_width, x_glob, y_glob, z_glob)
  !x_glob = x_glob + xc_g
  !y_glob = y_glob + yc_g
  !z_glob = z_glob + zc_g
  !call globalToLocalCoords(x_glob, y_glob, z_glob, x_loc, y_loc, z_loc, containedQ=dummy_flag)
      !! if particle is within the current MPI meshblock
  !if (dummy_flag) then
  !call localToCellBasedCoords(x_loc, y_loc, z_loc, xi, yi, zi, dx, dy, dz)
        !! initialize at rest
  !call interpFromEdges(dx, dy, dz, xi, yi, zi, ex, ey, ez, ex0, ey0, ez0)
  !call interpFromFaces(dx, dy, dz, xi, yi, zi, bx, by, bz, bx0, by0, bz0)
  !dummy_ = abs(ex0 * bx0 + ey0 * by0 + ez0 * bz0) / (bx0**2 + by0**2 + bz0**2)
  !if (dummy_ .gt. edotb_thr_closed) then
  !weight = dummy_ * nGJ / ppc
  !call createParticle(1, xi, yi, zi, dx, dy, dz, 0.0, 0.0, 0.0, weight=weight)
  !call createParticle(2, xi, yi, zi, dx, dy, dz, 0.0, 0.0, 0.0, weight=weight)
  !end if ! E.B limiter
  !end if ! current MPI block
  !end do

    !! max radius for spherical absorption and absorbing layer size in cartesian absorption
  !rmax_sph = 0.5 * MIN(global_mesh%sx, global_mesh%sy, global_mesh%sz) - ds_abs + 4.0
  !smax_car = ds_abs / 4.0
    !! remove particles falling into the star
  !do s = 1, nspec
  !do ti = 1, species(s)%tile_nx
  !do tj = 1, species(s)%tile_ny
  !do tk = 1, species(s)%tile_nz
  !do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
  !x_g = REAL(species(s)%prtl_tile(ti, tj, tk)%xi(p) + this_meshblock%ptr%x0)&
  !& + species(s)%prtl_tile(ti, tj, tk)%dx(p)
  !y_g = REAL(species(s)%prtl_tile(ti, tj, tk)%yi(p) + this_meshblock%ptr%y0)&
  !& + species(s)%prtl_tile(ti, tj, tk)%dy(p)
  !z_g = REAL(species(s)%prtl_tile(ti, tj, tk)%zi(p) + this_meshblock%ptr%z0)&
  !& + species(s)%prtl_tile(ti, tj, tk)%dz(p)
  !r_g = sqrt((x_g - xc_g)**2 + (y_g - yc_g)**2 + (z_g - zc_g)**2)
  !if (&
  !& (r_g .lt. (psr_radius - rmin_dr)) .or.&
  !& ((boundary_x .eq. 2) .and. (r_g .gt. rmax_sph)) .or.&
  !& ((boundary_x .eq. 0) .and. ((x_g .lt. smax_car) .or.&
  !& (x_g .gt. global_mesh%sx - smax_car) .or.&
  !& (y_g .lt. smax_car) .or.&
  !& (y_g .gt. global_mesh%sy - smax_car) .or.&
  !& (z_g .lt. smax_car) .or.&
  !& (z_g .gt. global_mesh%sz - smax_car))&
  !& )&
  !& ) then
  !species(s)%prtl_tile(ti, tj, tk)%proc(p) = -1
  !end if
  !end do
  !end do
  !end do
  !end do
  !end do
  !end subroutine userParticleBoundaryConditions

  subroutine randomPointInSphericalShell(rmin, rmax, x, y, z)
    implicit none
    real, intent(in) :: rmin, rmax
    real, intent(out) :: x, y, z
    real :: TH, X0, Y0, Z0, T, R
    TH = random(dseed) * 2.0 * M_PI
    Z0 = 2.0 * (random(dseed) - 0.5)
    X0 = sqrt(1.0 - Z0**2) * cos(TH)
    Y0 = sqrt(1.0 - Z0**2) * sin(TH)
    T = random(dseed) * (rmax**3 - rmin**3) + rmin**3
    R = T**(1.0 / 3.0)
    x = X0 * R
    y = Y0 * R
    z = Z0 * R
  end subroutine randomPointInSphericalShell

  real function distance3D_sq(xA, yA, zA, xB, yB, zB)
    real, intent(in) :: xA, yA, zA, xB, yB, zB
    distance3D_sq = (xA - xB)**2 + (yA - yB)**2 + (zA - zB)**2
  end function distance3D_sq

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
    logical :: updateE_, updateB_
    integer :: i, j, k, i_glob, j_glob, k_glob, mx, my, mz
    real :: supersph_radius_sq
    real, allocatable :: ex_new(:, :, :), ey_new(:, :, :), ez_new(:, :, :), &
                         bx_new(:, :, :), by_new(:, :, :), bz_new(:, :, :)
    real :: rx, ry, rz, rr_sqr, shift_B, s
    real :: bx_dip, by_dip, bz_dip, b_dip_dot_r, b_int_dot_r
    real :: scaleEpar, scaleEperp, scaleBperp, scaleBpar, scale
    real :: vx, vy, vz, ex_dip, ey_dip, ez_dip
    real :: shift_E, e_int_dot_r, e_dip_dot_r
    real :: rr, x_, y_, z_, rlimit
    logical :: dummy_log

    ! update only within this "supersphere"
    supersph_radius_sq = (psr_radius * 2)**2
    ! test if meshblock "touches" the supersphere
    x_ = max(REAL(this_meshblock % ptr % x0 - NGHOST), min(REAL(this_meshblock % ptr % x0 + this_meshblock % ptr % sx - 1 + NGHOST), xc_g))
    y_ = max(REAL(this_meshblock % ptr % y0 - NGHOST), min(REAL(this_meshblock % ptr % y0 + this_meshblock % ptr % sy - 1 + NGHOST), yc_g))
    z_ = max(REAL(this_meshblock % ptr % z0 - NGHOST), min(REAL(this_meshblock % ptr % z0 + this_meshblock % ptr % sz - 1 + NGHOST), zc_g))
    dummy_log = (distance3D_sq(x_, y_, z_, xc_g, yc_g, zc_g) .lt. supersph_radius_sq)

    if (dummy_log) then

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

      if (updateE_) then
        allocate (ex_new(-NGHOST:this_meshblock % ptr % sx - 1 + NGHOST, &
                         -NGHOST:this_meshblock % ptr % sy - 1 + NGHOST, &
                         -NGHOST:this_meshblock % ptr % sz - 1 + NGHOST))
        allocate (ey_new(-NGHOST:this_meshblock % ptr % sx - 1 + NGHOST, &
                         -NGHOST:this_meshblock % ptr % sy - 1 + NGHOST, &
                         -NGHOST:this_meshblock % ptr % sz - 1 + NGHOST))
        allocate (ez_new(-NGHOST:this_meshblock % ptr % sx - 1 + NGHOST, &
                         -NGHOST:this_meshblock % ptr % sy - 1 + NGHOST, &
                         -NGHOST:this_meshblock % ptr % sz - 1 + NGHOST))
      end if
      if (updateB_) then
        allocate (bx_new(-NGHOST:this_meshblock % ptr % sx - 1 + NGHOST, &
                         -NGHOST:this_meshblock % ptr % sy - 1 + NGHOST, &
                         -NGHOST:this_meshblock % ptr % sz - 1 + NGHOST))
        allocate (by_new(-NGHOST:this_meshblock % ptr % sx - 1 + NGHOST, &
                         -NGHOST:this_meshblock % ptr % sy - 1 + NGHOST, &
                         -NGHOST:this_meshblock % ptr % sz - 1 + NGHOST))
        allocate (bz_new(-NGHOST:this_meshblock % ptr % sx - 1 + NGHOST, &
                         -NGHOST:this_meshblock % ptr % sy - 1 + NGHOST, &
                         -NGHOST:this_meshblock % ptr % sz - 1 + NGHOST))
      end if

      ! B fields are set few cells below the E fields
      shift_B = 4.0
      shift_E = 0.0

      scaleEpar = 0.5
      scaleEperp = 0.25
      scaleBperp = 0.5
      scaleBpar = 0.25

      if (updateE_ .or. updateB_) then
        ! saving boundaries into `e/b_new` arrays
        do i = 0, this_meshblock % ptr % sx - 1
          i_glob = i + this_meshblock % ptr % x0
          do j = 0, this_meshblock % ptr % sy - 1
            j_glob = j + this_meshblock % ptr % y0
            do k = 0, this_meshblock % ptr % sz - 1
              k_glob = k + this_meshblock % ptr % z0
              if ((i_glob - xc_g)**2 + (j_glob - yc_g)**2 + (k_glob - zc_g)**2 .lt. supersph_radius_sq) then
                ! ... setting B-field
                if (updateB_) then
                  ! setting `B_x`
                  rx = REAL(i_glob) - xc_g
                  ry = REAL(j_glob) + 0.5 - yc_g
                  rz = REAL(k_glob) + 0.5 - zc_g
                  rr_sqr = rx**2 + ry**2 + rz**2
                  b_int_dot_r = bx(i, j, k) * rx + &
                                0.25 * (by(i, j, k) + by(i, j + 1, k) + by(i - 1, j, k) + by(i - 1, j + 1, k)) * ry + &
                                0.25 * (bz(i, j, k) + bz(i, j, k + 1) + bz(i - 1, j, k) + bz(i - 1, j, k + 1)) * rz
                  call getBfield(step, 0.0, rx + xc_g, ry + yc_g, rz + zc_g, bx_dip, by_dip, bz_dip)
                  b_dip_dot_r = bx_dip * rx + by_dip * ry + bz_dip * rz

                  scale = scaleBperp
                  s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_B) / scale)
                  bx_new(i, j, k) = (b_int_dot_r - b_dip_dot_r) * (rx / rr_sqr) * (1.0 - s)
                  scale = scaleBpar
                  s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_B) / scale)
                  bx_new(i, j, k) = bx_new(i, j, k) + bx_dip + &
                                    ((bx(i, j, k) - b_int_dot_r * rx / rr_sqr) - &
                                     (bx_dip - b_dip_dot_r * rx / rr_sqr)) * (1.0 - s)

                  ! setting `B_y`
                  rx = REAL(i_glob) + 0.5 - xc_g
                  ry = REAL(j_glob) - yc_g
                  rz = REAL(k_glob) + 0.5 - zc_g
                  rr_sqr = rx**2 + ry**2 + rz**2
                  b_int_dot_r = 0.25 * (bx(i, j, k) + bx(i + 1, j, k) + bx(i, j - 1, k) + bx(i + 1, j - 1, k)) * rx + &
                                by(i, j, k) * ry + &
                                0.25 * (bz(i, j, k) + bz(i, j, k + 1) + bz(i, j - 1, k) + bz(i, j - 1, k + 1)) * rz
                  call getBfield(step, 0.0, rx + xc_g, ry + yc_g, rz + zc_g, bx_dip, by_dip, bz_dip)
                  b_dip_dot_r = bx_dip * rx + by_dip * ry + bz_dip * rz

                  scale = scaleBperp
                  s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_B) / scale)
                  by_new(i, j, k) = (b_int_dot_r - b_dip_dot_r) * (ry / rr_sqr) * (1.0 - s)
                  scale = scaleBpar
                  s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_B) / scale)
                  by_new(i, j, k) = by_new(i, j, k) + by_dip + &
                                    ((by(i, j, k) - b_int_dot_r * ry / rr_sqr) - &
                                     (by_dip - b_dip_dot_r * ry / rr_sqr)) * (1.0 - s)

                  ! setting `B_z`
                  rx = REAL(i_glob) + 0.5 - xc_g
                  ry = REAL(j_glob) + 0.5 - yc_g
                  rz = REAL(k_glob) - zc_g
                  rr_sqr = rx**2 + ry**2 + rz**2
                  b_int_dot_r = 0.25 * (bx(i, j, k) + bx(i + 1, j, k) + bx(i, j, k - 1) + bx(i + 1, j, k - 1)) * rx + &
                                0.25 * (by(i, j, k) + by(i, j + 1, k) + by(i, j, k - 1) + by(i, j + 1, k - 1)) * ry + &
                                bz(i, j, k) * rz
                  call getBfield(step, 0.0, rx + xc_g, ry + yc_g, rz + zc_g, bx_dip, by_dip, bz_dip)
                  b_dip_dot_r = bx_dip * rx + by_dip * ry + bz_dip * rz

                  scale = scaleBperp
                  s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_B) / scale)
                  bz_new(i, j, k) = (b_int_dot_r - b_dip_dot_r) * (rz / rr_sqr) * (1.0 - s)
                  scale = scaleBpar
                  s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_B) / scale)
                  bz_new(i, j, k) = bz_new(i, j, k) + bz_dip + &
                                    ((bz(i, j, k) - b_int_dot_r * rz / rr_sqr) - &
                                     (bz_dip - b_dip_dot_r * rz / rr_sqr)) * (1.0 - s)
                end if
                ! ... setting E-field
                if (updateE_) then
                  ! setting `E_x`
                  rx = REAL(i_glob) + 0.5 - xc_g
                  ry = REAL(j_glob) - yc_g
                  rz = REAL(k_glob) - zc_g
                  rr_sqr = rx**2 + ry**2 + rz**2
                  e_int_dot_r = ex(i, j, k) * rx + &
                                0.25 * (ey(i, j, k) + ey(i + 1, j, k) + ey(i, j - 1, k) + ey(i + 1, j - 1, k)) * ry + &
                                0.25 * (ez(i, j, k) + ez(i + 1, j, k) + ez(i, j, k - 1) + ez(i + 1, j, k - 1)) * rz
                  call getBfield(step, 0.0, rx + xc_g, ry + yc_g, rz + zc_g, bx_dip, by_dip, bz_dip)
                  vx = -psr_omega * ry
                  vy = psr_omega * rx
                  vz = 0.0
                  ex_dip = -(vy * bz_dip - vz * by_dip) * CCINV
                  ey_dip = (vx * bz_dip - vz * bx_dip) * CCINV
                  ez_dip = -(vx * by_dip - vy * bx_dip) * CCINV
                  e_dip_dot_r = ex_dip * rx + ey_dip * ry + ez_dip * rz

                  scale = scaleEpar
                  s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_E) / scale)
                  ex_new(i, j, k) = ex_dip + &
                                    ((ex(i, j, k) - e_int_dot_r * rx / rr_sqr) - &
                                     (ex_dip - e_dip_dot_r * rx / rr_sqr)) * (1.0 - s)
                  scale = scaleEperp
                  s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_E) / scale)
                  ex_new(i, j, k) = ex_new(i, j, k) + &
                                    (e_int_dot_r - e_dip_dot_r) * (rx / rr_sqr) * (1.0 - s)

                  ! setting `E_y`
                  rx = REAL(i_glob) - xc_g
                  ry = REAL(j_glob) + 0.5 - yc_g
                  rz = REAL(k_glob) - zc_g
                  rr_sqr = rx**2 + ry**2 + rz**2
                  e_int_dot_r = 0.25 * (ex(i, j, k) + ex(i, j + 1, k) + ex(i - 1, j, k) + ex(i - 1, j + 1, k)) * rx + &
                                ey(i, j, k) * ry + &
                                0.25 * (ez(i, j, k) + ez(i, j + 1, k) + ez(i, j, k - 1) + ez(i, j + 1, k - 1)) * rz
                  call getBfield(step, 0.0, rx + xc_g, ry + yc_g, rz + zc_g, bx_dip, by_dip, bz_dip)
                  vx = -psr_omega * ry
                  vy = psr_omega * rx
                  vz = 0.0
                  ex_dip = -(vy * bz_dip - vz * by_dip) * CCINV
                  ey_dip = (vx * bz_dip - vz * bx_dip) * CCINV
                  ez_dip = -(vx * by_dip - vy * bx_dip) * CCINV
                  e_dip_dot_r = ex_dip * rx + ey_dip * ry + ez_dip * rz

                  scale = scaleEpar
                  s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_E) / scale)
                  ey_new(i, j, k) = ey_dip + &
                                    ((ey(i, j, k) - e_int_dot_r * ry / rr_sqr) - &
                                     (ey_dip - e_dip_dot_r * ry / rr_sqr)) * (1.0 - s)
                  scale = scaleEperp
                  s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_E) / scale)
                  ey_new(i, j, k) = ey_new(i, j, k) + &
                                    (e_int_dot_r - e_dip_dot_r) * (ry / rr_sqr) * (1.0 - s)

                  ! setting `E_z`
                  rx = REAL(i_glob) - xc_g
                  ry = REAL(j_glob) - yc_g
                  rz = REAL(k_glob) + 0.5 - zc_g
                  rr_sqr = rx**2 + ry**2 + rz**2
                  e_int_dot_r = 0.25 * (ex(i, j, k) + ex(i, j, k + 1) + ex(i - 1, j, k) + ex(i - 1, j, k + 1)) * rx + &
                                0.25 * (ey(i, j, k) + ey(i, j - 1, k) + ey(i, j, k + 1) + ey(i, j - 1, k + 1)) * ry + &
                                ez(i, j, k) * rz
                  call getBfield(step, 0.0, rx + xc_g, ry + yc_g, rz + zc_g, bx_dip, by_dip, bz_dip)
                  vx = -psr_omega * ry
                  vy = psr_omega * rx
                  vz = 0.0
                  ex_dip = -(vy * bz_dip - vz * by_dip) * CCINV
                  ey_dip = (vx * bz_dip - vz * bx_dip) * CCINV
                  ez_dip = -(vx * by_dip - vy * bx_dip) * CCINV
                  e_dip_dot_r = ex_dip * rx + ey_dip * ry + ez_dip * rz

                  scale = scaleEpar
                  s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_E) / scale)
                  ez_new(i, j, k) = ez_dip + &
                                    ((ez(i, j, k) - e_int_dot_r * rz / rr_sqr) - &
                                     (ez_dip - e_dip_dot_r * rz / rr_sqr)) * (1.0 - s)
                  scale = scaleEperp
                  s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_E) / scale)
                  ez_new(i, j, k) = ez_new(i, j, k) + &
                                    (e_int_dot_r - e_dip_dot_r) * (rz / rr_sqr) * (1.0 - s)
                end if
              end if
            end do
          end do
        end do
        ! copying to real array
        do i = 0, this_meshblock % ptr % sx - 1
          i_glob = i + this_meshblock % ptr % x0
          do j = 0, this_meshblock % ptr % sy - 1
            j_glob = j + this_meshblock % ptr % y0
            do k = 0, this_meshblock % ptr % sz - 1
              k_glob = k + this_meshblock % ptr % z0
              if ((i_glob - xc_g)**2 + (j_glob - yc_g)**2 + (k_glob - zc_g)**2 .lt. supersph_radius_sq) then
                if (updateE_) then
                  ex(i, j, k) = ex_new(i, j, k)
                  ey(i, j, k) = ey_new(i, j, k)
                  ez(i, j, k) = ez_new(i, j, k)
                end if
                if (updateB_) then
                  bx(i, j, k) = bx_new(i, j, k)
                  by(i, j, k) = by_new(i, j, k)
                  bz(i, j, k) = bz_new(i, j, k)
                end if
              end if
            end do
          end do
        end do
      end if

      if (updateE_) deallocate (ex_new, ey_new, ez_new)
      if (updateB_) deallocate (bx_new, by_new, bz_new)

    end if
  end subroutine userFieldBoundaryConditions
  !............................................................!

  !--- auxiliary functions ------------------------------------!
  subroutine getBfield(step, offset, x_g, y_g, z_g, &
                       obx, oby, obz)
    integer, intent(in) :: step
    real, intent(in) :: x_g, y_g, z_g, offset
    real, intent(out) :: obx, oby, obz
    call getDipole(step, offset, x_g, y_g, z_g, obx, oby, obz)
  end subroutine getBfield

  subroutine getDipole(step, offset, x_g, y_g, z_g, &
                       obx, oby, obz)
    implicit none
    integer, intent(in) :: step
    real, intent(in) :: x_g, y_g, z_g, offset
    real, intent(out) :: obx, oby, obz
    real :: phase, nx, ny, nz, rr, mux, muy, muz, mu_dot_n

    phase = psr_omega * step
    nx = x_g - xc_g
    ny = y_g - yc_g
    nz = z_g - zc_g

    rr = 1.0 / sqrt(nx**2 + ny**2 + nz**2)
    nx = nx * rr
    ny = ny * rr
    nz = nz * rr
    rr = rr**3

    mux = psr_radius**3 * sin(psr_angle) * cos(phase + offset)
    muy = psr_radius**3 * sin(psr_angle) * sin(phase + offset)
    muz = psr_radius**3 * cos(psr_angle)

    mu_dot_n = mux * nx + muy * ny + muz * nz

    obx = psr_bstar * (3.0 * nx * mu_dot_n - mux) * rr
    oby = psr_bstar * (3.0 * ny * mu_dot_n - muy) * rr
    obz = psr_bstar * (3.0 * nz * mu_dot_n - muz) * rr
  end subroutine getDipole

  real function shape(rad, rad0)
    implicit none
    real, intent(in) :: rad, rad0
    real :: del
    del = 1.0
    shape = 0.5 * (1.0 - tanh((rad - rad0) / del))
  end function shape

  !............................................................!

  !--- user-specific output -----------------------------------!
#ifdef USROUTPUT
  subroutine userOutput(step)
    implicit none
    integer, optional, intent(in) :: step
    integer :: root_rank = 0
    real, allocatable :: r_bins(:)
    real :: dr, x_glob, y_glob, z_glob, r_glob, fr_factor
    real :: dummy_x, dummy_y, dummy_z, dummy1, dummy2
    real, allocatable :: sum_ExBr_f(:), sum_f(:), sum_ExBr_f_global(:), sum_f_global(:)
    integer :: ri, rnum = 50, i, j, k, ierr

    allocate (r_bins(rnum))
    allocate (sum_ExBr_f(rnum), sum_f(rnum))
    allocate (sum_ExBr_f_global(rnum), sum_f_global(rnum))
    sum_ExBr_f(:) = 0.0
    sum_f(:) = 0.0

    dr = (REAL(global_mesh % sx) * 0.5 - psr_radius) / REAL(rnum)
    do ri = 1, rnum
      r_bins(ri) = psr_radius + ri * dr
    end do

    do i = 0, this_meshblock % ptr % sx - 1
      x_glob = REAL(i + this_meshblock % ptr % x0)
      do j = 0, this_meshblock % ptr % sy - 1
        y_glob = REAL(j + this_meshblock % ptr % y0)
        do k = 0, this_meshblock % ptr % sz - 1
          z_glob = REAL(k + this_meshblock % ptr % z0)
          r_glob = sqrt((x_glob - xc_g)**2 + (y_glob - yc_g)**2 + (z_glob - zc_g)**2)
          ! compute ExB_r
          dummy_x = -(ez(i, j, k) * by(i, j, k)) + ey(i, j, k) * bz(i, j, k)
          dummy_y = ez(i, j, k) * bx(i, j, k) - ex(i, j, k) * bz(i, j, k)
          dummy_z = -(ey(i, j, k) * bx(i, j, k)) + ex(i, j, k) * by(i, j, k)
          if (r_glob .gt. psr_radius / 2.0) then
            dummy1 = (dummy_x * (x_glob - xc_g) + &
                      dummy_y * (y_glob - yc_g) + &
                      dummy_z * (z_glob - zc_g)) / r_glob
          else
            dummy1 = 0.0
          end if
          do ri = 1, rnum
            fr_factor = exp(-(r_glob - r_bins(ri))**2 / (dr * 0.5)**2)
            sum_ExBr_f(ri) = sum_ExBr_f(ri) + dummy1 * fr_factor
            sum_f(ri) = sum_f(ri) + fr_factor
          end do
        end do
      end do
    end do

    call MPI_REDUCE(sum_ExBr_f, sum_ExBr_f_global, rnum, default_mpi_real, MPI_SUM, root_rank, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(sum_f, sum_f_global, rnum, default_mpi_real, MPI_SUM, root_rank, MPI_COMM_WORLD, ierr)

    if (mpi_rank .eq. root_rank) then
      ! normalizations
      sum_ExBr_f_global(:) = sum_ExBr_f_global(:) * r_bins(:)**2 * CC * B_norm**2 / sum_f_global(:)
      call writeUsrOutputTimestep(step)
      call writeUsrOutputArray('r_bins', r_bins)
      call writeUsrOutputArray('ExB_flux', sum_ExBr_f_global)
      call writeUsrOutputEnd()
    end if
  end subroutine userOutput

  logical function userExcludeParticles(s, ti, tj, tk, p)
    implicit none
    integer, intent(in) :: s, ti, tj, tk, p
    real :: xx, yy, zz, rr
    real :: uu, vv, ww, gg

    xx = REAL(this_meshblock % ptr % x0 + species(s) % prtl_tile(ti, tj, tk) % xi(p))
    xx = xx + species(s) % prtl_tile(ti, tj, tk) % dx(p)
    xx = xx - REAL(global_mesh % sx) * 0.5

    yy = REAL(this_meshblock % ptr % y0 + species(s) % prtl_tile(ti, tj, tk) % yi(p))
    yy = yy + species(s) % prtl_tile(ti, tj, tk) % dy(p)
    yy = yy - REAL(global_mesh % sy) * 0.5

    zz = REAL(this_meshblock % ptr % z0 + species(s) % prtl_tile(ti, tj, tk) % zi(p))
    zz = zz + species(s) % prtl_tile(ti, tj, tk) % dz(p)
    zz = zz - REAL(global_mesh % sz) * 0.5

    rr = sqrt(xx**2 + yy**2 + zz**2)

    uu = REAL(species(s) % prtl_tile(ti, tj, tk) % u(p))
    vv = REAL(species(s) % prtl_tile(ti, tj, tk) % v(p))
    ww = REAL(species(s) % prtl_tile(ti, tj, tk) % w(p))
    gg = sqrt(1.0 + uu**2 + vv**2 + ww**2)

    userExcludeParticles = ((rr .gt. global_usr_variable_1 + 40) .and. &
                            (rr .lt. REAL(global_mesh % sx) * 0.5 - 100) .and. &
                            (gg .gt. 20.0))
  end function userExcludeParticles
#endif
  !............................................................!

#include "optional.F"
end module m_userfile
