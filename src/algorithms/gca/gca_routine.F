#include "../../defs.F90"

! save interpolation of fields at `t = n`
ex0_n = ex0;  ey0_n = ey0;  ez0_n = ez0
bx0_n = bx0;  by0_n = by0;  bz0_n = bz0
bb0_n = sqrt(bx0_n**2 + by0_n**2 + bz0_n**2)

! compute position at `t = n-1/2`
! and save position at `t = n`
#if defined(oneD) || defined(twoD) || defined(threeD)
  x_ = (REAL(pt_xi(p) + pt_xi_past(p)) + (pt_dx(p) + pt_dx_past(p))) * 0.5
  xi_ = FLOOR(x_);   dx_ = x_ - REAL(xi_);
  x_n = REAL(pt_xi(p)) + pt_dx(p)
  pt_xi_past(p) = pt_xi(p);     pt_dx_past(p) = pt_dx(p)
#endif

#if defined(twoD) || defined(threeD)
  y_ = (REAL(pt_yi(p) + pt_yi_past(p)) + (pt_dy(p) + pt_dy_past(p))) * 0.5
  yi_ = FLOOR(y_);   dy_ = y_ - REAL(yi_)
  y_n = REAL(pt_yi(p)) + pt_dy(p)
  pt_yi_past(p) = pt_yi(p);     pt_dy_past(p) = pt_dy(p)
#endif

#if defined(threeD)
  z_ = (REAL(pt_zi(p) + pt_zi_past(p)) + (pt_dz(p) + pt_dz_past(p))) * 0.5
  zi_ = FLOOR(z_);   dz_ = z_ - REAL(zi_)
  z_n = REAL(pt_zi(p)) + pt_dz(p)
  pt_zi_past(p) = pt_zi(p);     pt_dz_past(p) = pt_dz(p)
#endif

#ifdef oneD
  yi_ = 0_2; dy_ = 0.5
  zi_ = 0_2; dz_ = 0.5
#elif twoD
  zi_ = 0_2; dz_ = 0.5
#endif

! save velocities at `t = n-1/2`
u0 = pt_u(p); v0 = pt_v(p); w0 = pt_w(p)
Gamma_ = sqrt(1.0 + u0**2 + v0**2 + w0**2)

! compute fields at `t = n-1/2`
#ifdef oneD
  lind = xi_
#elif twoD
  lind = xi_ + (NGHOST + yi_) * iy
#elif threeD
  lind = xi_ + (NGHOST + yi_) * iy + (NGHOST + zi_) * iz
#endif
dx = dx_; dy = dy_; dz = dz_

#ifndef DEBUG
  include "interp_efield.F"
  include "interp_bfield.F"
#else
  call interpFromEdges(dx, dy, dz, pt_xi(p), pt_yi(p), pt_zi(p),&
                     & ex, ey, ez, ex0, ey0, ez0)
  call interpFromFaces(dx, dy, dz, pt_xi(p), pt_yi(p), pt_zi(p),&
                     & bx, by, bz, bx0, by0, bz0)
#endif

e0_SQR = ex0**2 + ey0**2 + ez0**2
b0_SQR = bx0**2 + by0**2 + bz0**2
#include "gca_ecrossb_drift.F"

bb0 = sqrt(b0_SQR)
dummy_ = 1.0 / sqrt(b0_SQR + TINYFLD)

! at time `t = n-1/2`
bpar_x = bx0 * dummy_
bpar_y = by0 * dummy_
bpar_z = bz0 * dummy_

if (pt_proc(p) .ge. mpi_size) then
  ! particle was doing GCA last timestep
  u_par = pt_u_par(p)
  if (.not. gca_enforce_mu0) then
    uperp_g = pt_u_perp(p)

    ! to recover the direction later
    uperp_g_x = u0 - u_par * bpar_x - vE_x * Gamma_
    uperp_g_y = v0 - u_par * bpar_y - vE_y * Gamma_
    uperp_g_z = w0 - u_par * bpar_z - vE_z * Gamma_
    dummy_ = sqrt(uperp_g_x**2 + uperp_g_y**2 + uperp_g_z**2)
    uperp_g_x = uperp_g_x * uperp_g / dummy_
    uperp_g_y = uperp_g_y * uperp_g / dummy_
    uperp_g_z = uperp_g_z * uperp_g / dummy_
  else
    uperp_g = 0.0
  end if
else
  ! particle was doing normal push last timestep
  u_par = bpar_x * u0 + bpar_y * v0 + bpar_z * w0

  uperp_g_x = u0 - u_par * bpar_x - vE_x * Gamma_
  uperp_g_y = v0 - u_par * bpar_y - vE_y * Gamma_
  uperp_g_z = w0 - u_par * bpar_z - vE_z * Gamma_
  uperp_g = sqrt(uperp_g_x**2 + uperp_g_y**2 + uperp_g_z**2)

  if (gca_enforce_mu0) then
    ! update u_par to keep v_par constant (cooling in one timestep)
    ! ... but set u_perp_g = 0
    u_par = u_par * (gammaE / Gamma_) / sqrt(abs(1.0 - (u_par * gammaE / Gamma_)**2) + TINYFLD)
  end if
end if

doNormalPushQ = ((.not. species(s)%gca_sp) .or.&                                              ! if species cannot do GCA
          & (bb0 * bb0_n .eq. 0.0) .or.&                                                      ! if B = 0
          & ((CC**2 * Gamma_ / (abs(q_over_m) * bb0 * B_norm)) .ge. gca_rhomin) .or.&         ! if Larmor too large
          & (uperp_g / Gamma_ .gt. gca_vperpmax) .or.&                                        ! if uperp_g too large
          & (sqrt(e0_SQR) .ge. gca_eoverbmin * bb0))                                          ! of E > B * fraction

doNormalPushQ = (doNormalPushQ .and. (.not. userEnforceGCA(xi_, yi_, zi_, dx_, dy_, dz_,&     ! GCA enforced by the user ...
                                               & u0, v0, w0, pt_wei(p))))                     ! ... (defaults to .false.)

#ifdef DEBUG
  if (sqrt(1.0 + u0**2 + v0**2 + w0**2) .gt. 1e8) then
    print *, 'bad particle!!!!', s, pt_ind(p), pt_proc(p), sqrt(1.0 + u0**2 + v0**2 + w0**2)
    call throwError('ERROR: runaway particle found!')
  end if
  ! if ((pt_ind(p) .eq. 31731) .and. ((pt_proc(p) .eq. 283) .or. (pt_proc(p) .eq. 283 + 2 * mpi_size))) then
  !   dummy_ = sqrt(((x_ + this_meshblock%ptr%x0) - global_mesh%sx * 0.5)**2 +&
  !               & ((y_ + this_meshblock%ptr%y0) - global_mesh%sy * 0.5)**2 +&
  !               & ((z_ + this_meshblock%ptr%z0) - global_mesh%sz * 0.5)**2)
  !   print *, 'prtl:', pt_u(p), pt_v(p), pt_w(p), dummy_, uperp_g, pt_proc(p) .gt. mpi_size
  ! end if
#endif

! GCA/normal push condition . . . . . . . . . . . . . . . . . . . . . . . . . . . .
if (doNormalPushQ) then
  if (pt_proc(p) .ge. mpi_size) pt_proc(p) = pt_proc(p) - 2 * mpi_size
  ! use normal push here
  ex0 = ex0_n; ey0 = ey0_n; ez0 = ez0_n
  bx0 = bx0_n; by0 = by0_n; bz0 = bz0_n
  #ifndef VAY
    include "boris_push.F"
  #else
    include "vay_push.F"
  #endif
  pt_u(p) = u0;     pt_v(p) = v0;     pt_w(p) = w0
  over_e_temp = 1.0 / sqrt(1.0 + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
  ! compute "effective" velocities for the deposit
  include "position_update.F"
  pt_u_eff(p) = u0 * CC * over_e_temp
  pt_v_eff(p) = v0 * CC * over_e_temp
  pt_w_eff(p) = w0 * CC * over_e_temp
else
  if (pt_proc(p) .lt. mpi_size) pt_proc(p) = pt_proc(p) + 2 * mpi_size
  ! use GCA here
  if (.not. gca_enforce_mu0) then
    ! at time `t = n-1/2`
    dummy_ = 1.0 / (uperp_g + TINYFLD)
    bperp_x = uperp_g_x * dummy_
    bperp_y = uperp_g_y * dummy_
    bperp_z = uperp_g_z * dummy_

    ! `mu` at time `t = n-1/2` (in units of mc^2/2)
    mu_ = uperp_g**2 / (bb0 * gammaE)
  else
    mu_ = 0.0
  end if

  ! compute E_|| at `t = n`
  epar0_n = (ex0_n * bx0_n + ey0_n * by0_n + ez0_n * bz0_n) / bb0_n
  ! update u_||
  u_par = u_par + q_over_m * B_norm * epar0_n * CCINV

  ! . . . . Solve nonlinear equation on coordinates at `t = n+1` . . . . . . . . . . . . .
  ! compute `vE` and `gammaE` at time `t = n` - - >
  ex0 = ex0_n; ey0 = ey0_n; ez0 = ez0_n
  bx0 = bx0_n; by0 = by0_n; bz0 = bz0_n
  #include "gca_ecrossb_drift.F"

  vE_x_n = vE_x; vE_y_n = vE_y; vE_z_n = vE_z; gammaE_n = gammaE
  ! < - - compute `vE` and `gammaE` at time `t = n`

  ! `Gamma` at `t = n`
  Gamma_n = gammaE_n * sqrt(1.0 + u_par**2 + mu_ * bb0_n * gammaE_n)

  ! initial guess
  #if defined(oneD) || defined(twoD) || defined(threeD)
    x_n1 = x_n
  #endif

  #if defined(twoD) || defined(threeD)
    y_n1 = y_n
  #endif

  #if defined(threeD)
    z_n1 = z_n
  #endif

  do iter = 1, GCAITER
    ! convert coordinates at `t = n+1` to `xi` and `dx`
    #if defined(oneD) || defined(twoD) || defined(threeD)
      xi_ = FLOOR(x_n1); dx_ = x_n1 - REAL(xi_)
    #endif

    #if defined(twoD) || defined(threeD)
      yi_ = FLOOR(y_n1); dy_ = y_n1 - REAL(yi_)
    #endif

    #if defined(threeD)
      zi_ = FLOOR(z_n1); dz_ = z_n1 - REAL(zi_)
    #endif

    ! compute fields at `t = n+1`
    #ifdef oneD
      lind = xi_
    #elif twoD
      lind = xi_ + (NGHOST + yi_) * iy
    #elif threeD
      lind = xi_ + (NGHOST + yi_) * iy + (NGHOST + zi_) * iz
    #endif
    dx = dx_; dy = dy_; dz = dz_
    #ifndef DEBUG
      include "interp_efield.F"
      include "interp_bfield.F"
    #else
      call interpFromEdges(dx, dy, dz, pt_xi(p), pt_yi(p), pt_zi(p),&
                         & ex, ey, ez, ex0, ey0, ez0)
      call interpFromFaces(dx, dy, dz, pt_xi(p), pt_yi(p), pt_zi(p),&
                         & bx, by, bz, bx0, by0, bz0)
    #endif

    e0_SQR = ex0**2 + ey0**2 + ez0**2
    b0_SQR = bx0**2 + by0**2 + bz0**2
    #include "gca_ecrossb_drift.F"

    vE_x_n1 = vE_x;   vE_y_n1 = vE_y;   vE_z_n1 = vE_z;   gammaE_n1 = gammaE
    bx0_n1 = bx0;     by0_n1 = by0;     bz0_n1 = bz0
    bb0_n1 = sqrt(b0_SQR)

    ! `Gamma` at `t = n`
    Gamma_n1 = gammaE_n1 * sqrt(1.0 + u_par**2 + mu_ * bb0_n1 * gammaE_n1)

    ! update coordinates
    #if defined(oneD) || defined(twoD) || defined(threeD)
      x_n1 = x_n + 0.5 * ((bx0_n / bb0_n) / Gamma_n + (bx0_n1 / bb0_n1) / Gamma_n1) * u_par * CC +&
                 & 0.5 * (vE_x_n + vE_x_n1) * CC
    #endif

    #if defined(twoD) || defined(threeD)
      y_n1 = y_n + 0.5 * ((by0_n / bb0_n) / Gamma_n + (by0_n1 / bb0_n1) / Gamma_n1) * u_par * CC +&
                 & 0.5 * (vE_y_n + vE_y_n1) * CC
    #endif

    #if defined(threeD)
      z_n1 = z_n + 0.5 * ((bz0_n / bb0_n) / Gamma_n + (bz0_n1 / bb0_n1) / Gamma_n1) * u_par * CC +&
                 & 0.5 * (vE_z_n + vE_z_n1) * CC
    #endif
  end do
  ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  ! compute "effective" velocities for the deposit
  #ifdef oneD
    pt_u_eff(p) = x_n1 - x_n
    pt_v_eff(p) = 0.5 * ((by0_n / bb0_n) / Gamma_n + (by0_n1 / bb0_n1) / Gamma_n1) * u_par * CC +&
                & 0.5 * (vE_y_n + vE_y_n1) * CC
    pt_w_eff(p) = 0.5 * ((bz0_n / bb0_n) / Gamma_n + (bz0_n1 / bb0_n1) / Gamma_n1) * u_par * CC +&
                & 0.5 * (vE_z_n + vE_z_n1) * CC
  #elif twoD
    pt_u_eff(p) = x_n1 - x_n
    pt_v_eff(p) = y_n1 - y_n
    pt_w_eff(p) = 0.5 * ((bz0_n / bb0_n) / Gamma_n + (bz0_n1 / bb0_n1) / Gamma_n1) * u_par * CC +&
                & 0.5 * (vE_z_n + vE_z_n1) * CC
  #elif threeD
    pt_u_eff(p) = x_n1 - x_n
    pt_v_eff(p) = y_n1 - y_n
    pt_w_eff(p) = z_n1 - z_n
  #endif

  ! convert coordinates at `t = n+1` to `xi` and `dx` ...
  ! ... compute position at `t = n+1/2`
  ! ... and update particle coordinates
  #if defined(oneD) || defined(twoD) || defined(threeD)
    xi_ = FLOOR(x_n1);   dx_ = x_n1 - REAL(xi_)
    pt_xi(p) = xi_;       pt_dx(p) = dx_
    x_ = (REAL(pt_xi(p) + pt_xi_past(p)) + (pt_dx(p) + pt_dx_past(p))) * 0.5
    xi_ = FLOOR(x_);     dx_ = x_ - REAL(xi_)
  #endif

  #if defined(twoD) || defined(threeD)
    yi_ = FLOOR(y_n1);   dy_ = y_n1 - REAL(yi_)
    pt_yi(p) = yi_;       pt_dy(p) = dy_
    y_ = (REAL(pt_yi(p) + pt_yi_past(p)) + (pt_dy(p) + pt_dy_past(p))) * 0.5
    yi_ = FLOOR(y_);     dy_ = y_ - REAL(yi_)
  #endif

  #if defined(threeD)
    zi_ = FLOOR(z_n1);   dz_ = z_n1 - REAL(zi_)
    pt_zi(p) = zi_;       pt_dz(p) = dz_
    z_ = (REAL(pt_zi(p) + pt_zi_past(p)) + (pt_dz(p) + pt_dz_past(p))) * 0.5
    zi_ = FLOOR(z_);     dz_ = z_ - REAL(zi_)
  #endif

  ! interpolate fields at `t = n+1/2`
  #ifdef oneD
    lind = xi_
  #elif twoD
    lind = xi_ + (NGHOST + yi_) * iy
  #elif threeD
    lind = xi_ + (NGHOST + yi_) * iy + (NGHOST + zi_) * iz
  #endif
  dx = dx_; dy = dy_; dz = dz_
  #ifndef DEBUG
    include "interp_efield.F"
    include "interp_bfield.F"
  #else
    call interpFromEdges(dx, dy, dz, pt_xi(p), pt_yi(p), pt_zi(p),&
                       & ex, ey, ez, ex0, ey0, ez0)
    call interpFromFaces(dx, dy, dz, pt_xi(p), pt_yi(p), pt_zi(p),&
                       & bx, by, bz, bx0, by0, bz0)
  #endif

  e0_SQR = ex0**2 + ey0**2 + ez0**2
  b0_SQR = bx0**2 + by0**2 + bz0**2
  #include "gca_ecrossb_drift.F"

  ! at time `t = n+1/2`
  bb0 = sqrt(b0_SQR)
  uperp_g = sqrt(mu_ * bb0 * gammaE)

  dummy_ = 1.0 / sqrt(b0_SQR + TINYFLD)

  bpar_x = bx0 * dummy_
  bpar_y = by0 * dummy_
  bpar_z = bz0 * dummy_

  Gamma_ = gammaE * sqrt(1.0 + u_par**2 + uperp_g**2)

  if (.not. gca_enforce_mu0) then
    ! recover cartesian velocities
    pt_u(p) = bpar_x * u_par + vE_x * Gamma_ + uperp_g * bperp_x
    pt_v(p) = bpar_y * u_par + vE_y * Gamma_ + uperp_g * bperp_y
    pt_w(p) = bpar_z * u_par + vE_z * Gamma_ + uperp_g * bperp_z
  else
    ! recover cartesian velocities
    pt_u(p) = bpar_x * u_par + vE_x * Gamma_
    pt_v(p) = bpar_y * u_par + vE_y * Gamma_
    pt_w(p) = bpar_z * u_par + vE_z * Gamma_
  end if

  ! save par/perp velocities
  pt_u_par(p) = u_par
  pt_u_perp(p) = uperp_g
end if ! normal pusher / gca criterion
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
