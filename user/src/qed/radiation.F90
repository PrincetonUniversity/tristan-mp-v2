module m_radiation
  use m_globalnamespace
  use m_outputnamespace
  use m_qednamespace
  use m_readinput, only: getInput

  use m_aux
  use m_errors
  use m_particlelogistics
#ifdef RADIATION
  implicit none

  !--- PRIVATE variables/functions -------------------------------!
  !...............................................................!
contains
  subroutine initializeRadiation()
    implicit none
    call getInput('radiation', 'interval', rad_interval, 1)
    call getInput('radiation', 'emit_gamma_syn', emit_gamma_syn, 10.0)
    call getInput('radiation', 'emit_gamma_ic', emit_gamma_ic, 10.0)
    call getInput('radiation', 'gamma_syn', cool_gamma_syn, 10.0)
    call getInput('radiation', 'gamma_ic', cool_gamma_ic, 10.0)
    call getInput('radiation', 'beta_rec', rad_beta_rec, 0.1)
    call getInput('radiation', 'dens_limit', rad_dens_lim, 0.0)
    call getInput('radiation', 'cool_limit', rad_cool_lim, 0.0)
#ifdef EMIT
    call getInput('radiation', 'photon_sp', rad_photon_sp, 3)
    if ((rad_photon_sp .le. 0) .or. &
        (nspec .lt. rad_photon_sp) .or. &
        (species(rad_photon_sp) % ch_sp .ne. 0) .or. &
        (species(rad_photon_sp) % m_sp .ne. 0)) then
      call throwError('Wrong choice of `photon_sp`.')
    end if
#endif
    call printDiag("initializeRadiation()", 1)
  end subroutine initializeRadiation

  subroutine particleRadiateSync(timestep, s, &
                                 u0, v0, w0, ui, vi, wi, &
                                 dx, dy, dz, xi, yi, zi, &
                                 weight, &
                                 bx, by, bz, ex, ey, ez, &
                                 index, proc)
    implicit none
    integer, intent(in) :: timestep
    real, intent(inout) :: u0, v0, w0
    real, intent(in) :: ui, vi, wi
    real, intent(in) :: bx, by, bz, ex, ey, ez
    real, intent(in) :: dx, dy, dz
    integer(kind=2), intent(in) :: xi, yi, zi
    real, intent(in) :: weight
    integer, intent(in) :: s, index, proc

    real :: uci, vci, wci, g0, gci, betaci, over_gci
#ifdef EMIT
    real :: kx, ky, kz
#endif

    real :: e_bar_x, e_bar_y, e_bar_z, e_bar_sq, beta_dot_e
    real :: chiR, chiR_sq, kappaR_x, kappaR_y, kappaR_z
    real :: tau_emit, eph_emit, dummy_

    integer :: spec_index

    g0 = sqrt(1.0 + u0**2 + v0**2 + w0**2)
    if ((g0 .gt. 1.5) .and. &
        ((rad_dens_lim .eq. 0) .or. (lg_arr(xi, yi, zi) / ppc0 .lt. rad_dens_lim)) .and. &
        (cool_gamma_syn .gt. 0.0) &
        ) then

      uci = 0.5 * (u0 + ui)
      vci = 0.5 * (v0 + vi)
      wci = 0.5 * (w0 + wi)

      gci = sqrt(1.0 + uci**2 + vci**2 + wci**2)
      over_gci = 1.0 / gci
      betaci = sqrt(1.0 - over_gci**2)

      e_bar_x = ex + (vci * bz - wci * by) * over_gci
      e_bar_y = ey + (wci * bx - uci * bz) * over_gci
      e_bar_z = ez + (uci * by - vci * bx) * over_gci
      e_bar_sq = e_bar_x**2 + e_bar_y**2 + e_bar_z**2
      beta_dot_e = (ex * uci + ey * vci + ez * wci) * over_gci

      chiR_sq = abs(e_bar_sq - beta_dot_E**2)
      chiR = sqrt(chiR_sq)

      kappaR_x = (bz * e_bar_y - by * e_bar_z) + (ex * beta_dot_e)
      kappaR_y = (-bz * e_bar_x + bx * e_bar_z) + (ey * beta_dot_e)
      kappaR_z = (by * e_bar_x - bx * e_bar_y) + (ez * beta_dot_e)

      dummy_ = B_norm * rad_beta_rec * CCINV / cool_gamma_syn**2

      if ((rad_cool_lim .gt. 0) .and. (dummy_ * chiR_sq * gci .gt. rad_cool_lim)) then
        dummy_ = rad_cool_lim / (chiR_sq * gci)
        call addWarning(1)
      end if

      tau_emit = dummy_ * betaci * emit_gamma_syn**2 * chiR
      eph_emit = (gci / emit_gamma_syn)**2 * chiR

#ifndef EMIT
      u0 = u0 + dummy_ * (kappaR_x - chiR_sq * gci * uci)
      v0 = v0 + dummy_ * (kappaR_y - chiR_sq * gci * vci)
      w0 = w0 + dummy_ * (kappaR_z - chiR_sq * gci * wci)
#else
      u0 = u0 + dummy_ * kappaR_x
      v0 = v0 + dummy_ * kappaR_y
      w0 = w0 + dummy_ * kappaR_z

      over_gci = 1.0 / sqrt(uci**2 + vci**2 + wci**2)
      kx = uci * over_gci; ky = vci * over_gci; kz = wci * over_gci
      u0 = u0 - tau_emit * kx * eph_emit
      v0 = v0 - tau_emit * ky * eph_emit
      w0 = w0 - tau_emit * kz * eph_emit

      if ((random(dseed) .lt. tau_emit) .and. (modulo(timestep + index, rad_interval) .eq. 0)) then
        call createParticle(rad_photon_sp, xi, yi, zi, dx, dy, dz, &
                            kx * eph_emit, ky * eph_emit, kz * eph_emit, weight=(weight * rad_interval))
      end if
#endif

      if (eph_emit .gt. TINYFLD) then
        if (spec_log_bins) eph_emit = log(eph_emit)
        if (eph_emit .le. rad_spec_min) then
          spec_index = 1
        else if (eph_emit .ge. rad_spec_max) then
          spec_index = rad_spec_num
        else
          spec_index = INT(CEILING((eph_emit - rad_spec_min) * REAL(rad_spec_num) / (rad_spec_max - rad_spec_min)))
          if (spec_index .lt. 1) spec_index = 1
          if (spec_index .gt. rad_spec_num) spec_index = rad_spec_num
        end if
        rad_spectra(s, spec_index) = rad_spectra(s, spec_index) + tau_emit * weight
      end if
    end if
  end subroutine particleRadiateSync

  subroutine particleRadiateIC(timestep, s, &
                               u0, v0, w0, ui, vi, wi, &
                               dx, dy, dz, xi, yi, zi, &
                               weight, &
                               bx, by, bz, ex, ey, ez, &
                               index, proc)
    implicit none
    integer, intent(in) :: timestep
    real, intent(inout) :: u0, v0, w0
    real, intent(in) :: ui, vi, wi
    real, intent(in) :: bx, by, bz, ex, ey, ez
    real, intent(in) :: dx, dy, dz
    integer(kind=2), intent(in) :: xi, yi, zi
    real, intent(in) :: weight
    integer, intent(in) :: s, index, proc

    real :: uci, vci, wci, g0, gci, betaci, over_gci
#ifdef EMIT
    real :: kx, ky, kz
#endif

    real :: tau_emit, eph_emit, dummy_
    integer :: spec_index

    g0 = sqrt(1.0 + u0**2 + v0**2 + w0**2)
    if ((g0 .gt. 1.5) .and. &
        ((rad_dens_lim .eq. 0) .or. (lg_arr(xi, yi, zi) / ppc0 .lt. rad_dens_lim)) .and. &
        (cool_gamma_ic .gt. 0.0) &
        ) then

      uci = 0.5 * (u0 + ui)
      vci = 0.5 * (v0 + vi)
      wci = 0.5 * (w0 + wi)

      gci = sqrt(1.0 + uci**2 + vci**2 + wci**2)
      over_gci = 1.0 / gci
      betaci = sqrt(1.0 - over_gci**2)

      dummy_ = B_norm * rad_beta_rec * CCINV / cool_gamma_ic**2

      if ((rad_cool_lim .gt. 0) .and. (dummy_ * gci .gt. rad_cool_lim)) then
        dummy_ = rad_cool_lim / gci
        call addWarning(2)
      end if

      tau_emit = dummy_ * betaci * emit_gamma_ic**2
      eph_emit = (gci / emit_gamma_ic)**2

#ifndef EMIT
      u0 = u0 - dummy_ * gci * uci
      v0 = v0 - dummy_ * gci * vci
      w0 = w0 - dummy_ * gci * wci
#else
      over_gci = 1.0 / sqrt(uci**2 + vci**2 + wci**2)
      kx = uci * over_gci; ky = vci * over_gci; kz = wci * over_gci
      u0 = u0 - tau_emit * kx * eph_emit
      v0 = v0 - tau_emit * ky * eph_emit
      w0 = w0 - tau_emit * kz * eph_emit

      if ((random(dseed) .lt. tau_emit) .and. (modulo(timestep + index, rad_interval) .eq. 0)) then
        call createParticle(rad_photon_sp, xi, yi, zi, dx, dy, dz, &
                            kx * eph_emit, ky * eph_emit, kz * eph_emit, weight=(weight * rad_interval))
      end if
#endif
      if (eph_emit .gt. TINYFLD) then
        if (spec_log_bins) eph_emit = log(eph_emit)
        if (eph_emit .le. rad_spec_min) then
          spec_index = 1
        else if (eph_emit .ge. rad_spec_max) then
          spec_index = rad_spec_num
        else
          spec_index = INT(CEILING((eph_emit - rad_spec_min) * REAL(rad_spec_num) / (rad_spec_max - rad_spec_min)))
          if (spec_index .lt. 1) spec_index = 1
          if (spec_index .gt. rad_spec_num) spec_index = rad_spec_num
        end if
        rad_spectra(s, spec_index) = rad_spectra(s, spec_index) + tau_emit * weight
      end if
    end if
  end subroutine particleRadiateIC

#endif
end module m_radiation
