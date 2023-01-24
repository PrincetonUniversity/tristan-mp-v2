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
  implicit none

  !--- PRIVATE variables -----------------------------------------!

  real, private :: background_T, beta_kick
  logical, private :: init_J, init_rho, init_ufl, smooth_flds, single_tube
  integer, private :: nsmooth
  real, private :: ppc_buff
  real, private :: r_j, alpha, c_param, x1, y1, x2, y2 ! flux tube params
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution, filterInX, filterInY, filterFields
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'background_T', background_T, 1e-3)
    call getInput('problem', 'beta_kick', beta_kick, 0.1)
    call getInput('problem', 'c_param', c_param, 0.01)
    call getInput('problem', 'init_J', init_J, .true.)
    call getInput('problem', 'init_rho', init_rho, .true.)
    call getInput('problem', 'init_ufl', init_ufl, .true.)
    call getInput('problem', 'single_tube', single_tube, .false.)
    call getInput('problem', 'smooth_flds', smooth_flds, .true.)
    call getInput('problem', 'nsmooth', nsmooth, 8)
    call getInput('problem', 'ppc_buff', ppc_buff, ppc0)
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob, &
                                   dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3

    return
  end function

  subroutine userDeallocate()
    implicit none
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

  !
  ! This assumes userInitFields() is called *before*
  ! userInitParticles(), so that the fields are already properly set up.
  !
  subroutine userInitParticles()
    implicit none
    real :: background_n
    type(region) :: back_region
    integer :: s, ti, tj, tk, p
    integer(kind=2) :: i, j, k
    real :: ux, uy, uz, gam, wei, wei_new
    real :: beta_x, beta_y, beta_z, beta_sq
    real :: ux_boost, uy_boost, uz_boost, gam_boost, boost, beta_dot_u
    real :: dx, dy, dz, jx0, jy0, jz0, rho0
    real :: ex0, ey0, ez0, bx0, by0, bz0
    real :: x_, y_, r1, r2
    real :: dwn_ppc_buff, wei_buff
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    dwn_ppc_buff = REAL(ppc_buff) / REAL(ppc0)
    ! high ppc where the flux tubes are:
    background_n = REAL(ppc0) * 0.5
    !back_region%x_min = 0.0
    !back_region%x_max = REAL(global_mesh%sx)
    back_region % x_min = REAL(global_mesh % sx) * 0.2
    back_region % x_max = REAL(global_mesh % sx) * 0.8
    back_region % y_min = 0.0
    back_region % y_max = REAL(global_mesh % sy)
    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, background_n, background_T)

    ! less ppc in the buffer zones:
    background_n = background_n * dwn_ppc_buff
    wei_buff = 1.0 / dwn_ppc_buff
    back_region % x_min = 0.0
    back_region % x_max = REAL(global_mesh % sx) * 0.2
    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, background_n, background_T, weights=wei_buff)
    back_region % x_min = REAL(global_mesh % sx) * 0.8
    back_region % x_max = REAL(global_mesh % sx)
    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, background_n, background_T, weights=wei_buff)

    ! loop over all particles and give them a boost to setup current profile:
    do s = 1, nspec
      do ti = 1, species(s) % tile_nx
        do tj = 1, species(s) % tile_ny
          do tk = 1, species(s) % tile_nz
            do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp

              i = species(s) % prtl_tile(ti, tj, tk) % xi(p)
              j = species(s) % prtl_tile(ti, tj, tk) % yi(p)
              k = species(s) % prtl_tile(ti, tj, tk) % zi(p)
              dx = species(s) % prtl_tile(ti, tj, tk) % dx(p)
              dy = species(s) % prtl_tile(ti, tj, tk) % dy(p)
              dz = species(s) % prtl_tile(ti, tj, tk) % dz(p)
              wei = species(s) % prtl_tile(ti, tj, tk) % weight(p)
              ux = species(s) % prtl_tile(ti, tj, tk) % u(p)
              uy = species(s) % prtl_tile(ti, tj, tk) % v(p)
              uz = species(s) % prtl_tile(ti, tj, tk) % w(p)
              gam = sqrt(1.0 + ux**2 + uy**2 + uz**2)
              if ((dx .lt. 0) .or. (dy .lt. 0) .or. (dz .lt. 0)) then
                call throwError('ERROR: particle `dx`, `dy` or `dz` < 0 !')
              end if

              beta_x = 0.0; beta_y = 0.0; beta_z = 0.0
              if (init_J) then
                call interpFromEdges(dx, dy, dz, i, j, k, &
                                     jx, jy, jz, jx0, jy0, jz0)
                ! the conversion factors come from the normalization
                ! in the current deposit
                beta_x = jx0 * sqrt(sigma) * c_omp * sign(1.0, species(s) % ch_sp) / CC
                beta_y = jy0 * sqrt(sigma) * c_omp * sign(1.0, species(s) % ch_sp) / CC
                beta_z = jz0 * sqrt(sigma) * c_omp * sign(1.0, species(s) % ch_sp) / CC
              end if
              if (init_ufl) then
                !x_ = REAL(i + this_meshblock%ptr%x0) + dx
                !y_ = REAL(j + this_meshblock%ptr%y0) + dy
                !r1 = sqrt((x_ - x1)**2 + (y_ - y1)**2) / r_j
                !r2 = sqrt((x_ - x2)**2 + (y_ - y2)**2) / r_j
                !if (r1 .lt. 1.0) beta_y = beta_y + beta_kick
                !if ((r2 .lt. 1.0) .and. (.not. single_tube)) beta_y = beta_y - beta_kick
                call interpFromEdges(dx, dy, dz, i, j, k, &
                                     ex, ey, ez, ex0, ey0, ez0)
                call interpFromFaces(dx, dy, dz, i, j, k, &
                                     bx, by, bz, bx0, by0, bz0)
                beta_y = beta_y + (ez0 * bx0 - ex0 * bz0) / (bx0**2 + bz0**2)
              end if

              beta_sq = beta_x**2 + beta_y**2 + beta_z**2
              if (beta_sq .ge. 1) then
                call throwError('ERROR: `beta_sq` >= 1 in `userInitParticles()`')
              end if
              beta_dot_u = ux * beta_x + uy * beta_y + uz * beta_z
              if (-beta_dot_u / gam .gt. random(dseed)) then
                ux = ux - 2.0 * beta_dot_u * beta_x / beta_sq
                uy = uy - 2.0 * beta_dot_u * beta_y / beta_sq
                uz = uz - 2.0 * beta_dot_u * beta_z / beta_sq
              end if
              gam_boost = 1.0 / sqrt(1.0 - beta_sq)
              ux_boost = gam_boost * beta_x
              uy_boost = gam_boost * beta_y
              uz_boost = gam_boost * beta_z
              boost = (ux * ux_boost + uy * uy_boost + uz * uz_boost) / (gam_boost + 1) + gam
              ux = ux + boost * ux_boost
              uy = uy + boost * uy_boost
              uz = uz + boost * uz_boost
              species(s) % prtl_tile(ti, tj, tk) % u(p) = ux
              species(s) % prtl_tile(ti, tj, tk) % v(p) = uy
              species(s) % prtl_tile(ti, tj, tk) % w(p) = uz

              if (init_rho) then
                ! setup the (tiny) charge density perturbation:
                rho0 = lg_arr(i, j, 0) * (1.0 - dx) * (1.0 - dy)
                rho0 = rho0 + lg_arr(i + 1, j, 0) * dx * (1.0 - dy)
                rho0 = rho0 + lg_arr(i, j + 1, 0) * (1.0 - dx) * dy
                rho0 = rho0 + lg_arr(i + 1, j + 1, 0) * dx * dy
                wei_new = wei + rho0 * sqrt(sigma) * c_omp * sign(1.0, species(s) % ch_sp)
                if (abs(wei_new - wei) .gt. 0.25) then
                  call throwError('ERROR: |`wei_new` - `wei`| above the 0.25 safety threshold in `userInitParticles()`')
                end if
                if (wei_new .le. 0) then
                  call throwError('ERROR: `wei_new` < 0 in `userInitParticles()`')
                end if
                species(s) % prtl_tile(ti, tj, tk) % weight(p) = wei_new
              end if
            end do
          end do
        end do
      end do
    end do
    lg_arr(:, :, :) = 0.0

  end subroutine userInitParticles

  ! copy from (private) subroutine in m_filtering:
  subroutine filterInX(arr, do_n_times)
    implicit none
#ifdef twoD
    real, intent(inout) :: arr(-NGHOST:this_meshblock % ptr % sx - 1 + NGHOST, &
                               -NGHOST:this_meshblock % ptr % sy - 1 + NGHOST, 0:0)
#elif defined(threeD)
    real, intent(inout) :: arr(-NGHOST:this_meshblock % ptr % sx - 1 + NGHOST, &
                               -NGHOST:this_meshblock % ptr % sy - 1 + NGHOST, &
                               -NGHOST:this_meshblock % ptr % sz - 1 + NGHOST)
#endif
    integer, intent(in) :: do_n_times
    real :: tmp2, tmp1
    integer :: i, j, k, n_pass
    integer :: imin, imax, jmin, jmax, kmin, kmax
    if (do_n_times .gt. NGHOST) then
      call throwError('ERROR: `filterInX()` called with `do_n_times` > NGHOST.')
    end if
    do n_pass = 1, do_n_times
#ifdef twoD
      imin = -NGHOST + n_pass
      imax = this_meshblock % ptr % sx - 1 + NGHOST - n_pass
      jmin = -NGHOST + n_pass
      jmax = this_meshblock % ptr % sy - 1 + NGHOST - n_pass
      kmin = 0; kmax = 0
#elif defined(threeD)
      imin = -NGHOST + n_pass
      imax = this_meshblock % ptr % sx - 1 + NGHOST - n_pass
      jmin = -NGHOST + n_pass
      jmax = this_meshblock % ptr % sy - 1 + NGHOST - n_pass
      kmin = -NGHOST + n_pass
      kmax = this_meshblock % ptr % sz - 1 + NGHOST - n_pass
#endif
      if (modulo(imax - imin, 2) .eq. 0) then
        do k = kmin, kmax
          do j = jmin, jmax
            tmp2 = arr(imin - 1, j, k)
            do i = imin, imax - 1, 2
              tmp1 = 0.25 * arr(i - 1, j, k) + 0.5 * arr(i, j, k) + 0.25 * arr(i + 1, j, k)
              arr(i - 1, j, k) = tmp2
              tmp2 = 0.25 * arr(i, j, k) + 0.5 * arr(i + 1, j, k) + 0.25 * arr(i + 2, j, k)
              arr(i, j, k) = tmp1
            end do
            tmp1 = 0.25 * arr(imax - 1, j, k) + 0.5 * arr(imax, j, k) + 0.25 * arr(imax + 1, j, k)
            arr(imax - 1, j, k) = tmp2
            arr(imax, j, k) = tmp1
          end do
        end do
      else
        do k = kmin, kmax
          do j = jmin, jmax
            tmp2 = arr(imin - 1, j, k)
            do i = imin, imax, 2
              tmp1 = 0.25 * arr(i - 1, j, k) + 0.5 * arr(i, j, k) + 0.25 * arr(i + 1, j, k)
              arr(i - 1, j, k) = tmp2
              tmp2 = 0.25 * arr(i, j, k) + 0.5 * arr(i + 1, j, k) + 0.25 * arr(i + 2, j, k)
              arr(i, j, k) = tmp1
            end do
            arr(imax, j, k) = tmp2
          end do
        end do
      end if
    end do
  end subroutine filterInX

  ! copy from (private) subroutine in m_filtering:
  subroutine filterInY(arr, do_n_times)
    implicit none
#ifdef twoD
    real, intent(inout) :: arr(-NGHOST:this_meshblock % ptr % sx - 1 + NGHOST,&
                             & -NGHOST:this_meshblock % ptr % sy - 1 + NGHOST, 0:0)
#elif defined(threeD)
    real, intent(inout) :: arr(-NGHOST:this_meshblock % ptr % sx - 1 + NGHOST,&
                             & -NGHOST:this_meshblock % ptr % sy - 1 + NGHOST,&
                             & -NGHOST:this_meshblock % ptr % sz - 1 + NGHOST)
#endif
    integer, intent(in) :: do_n_times
    real :: tmp2, tmp1
    integer :: i, j, k, n_pass
    integer :: imin, imax, jmin, jmax, kmin, kmax
    if (do_n_times .gt. NGHOST) then
      call throwError('ERROR: `filterInY()` called with `do_n_times` > NGHOST.')
    end if
    do n_pass = 1, do_n_times
#ifdef twoD
      imin = -NGHOST + n_pass
      imax = this_meshblock % ptr % sx - 1 + NGHOST - n_pass
      jmin = -NGHOST + n_pass
      jmax = this_meshblock % ptr % sy - 1 + NGHOST - n_pass
      kmin = 0; kmax = 0
#elif defined(threeD)
      imin = -NGHOST + n_pass
      imax = this_meshblock % ptr % sx - 1 + NGHOST - n_pass
      jmin = -NGHOST + n_pass
      jmax = this_meshblock % ptr % sy - 1 + NGHOST - n_pass
      kmin = -NGHOST + n_pass
      kmax = this_meshblock % ptr % sz - 1 + NGHOST - n_pass
#endif
      if (modulo(jmax - jmin, 2) .eq. 0) then
        do k = kmin, kmax
          do i = imin, imax
            tmp2 = arr(i, jmin - 1, k)
            do j = jmin, jmax - 1, 2
              tmp1 = 0.25 * arr(i, j - 1, k) + 0.5 * arr(i, j, k) + 0.25 * arr(i, j + 1, k)
              arr(i, j - 1, k) = tmp2
              tmp2 = 0.25 * arr(i, j, k) + 0.5 * arr(i, j + 1, k) + 0.25 * arr(i, j + 2, k)
              arr(i, j, k) = tmp1
            end do
            tmp1 = 0.25 * arr(i, jmax - 1, k) + 0.5 * arr(i, jmax, k) + 0.25 * arr(i, jmax + 1, k)
            arr(i, jmax - 1, k) = tmp2
            arr(i, jmax, k) = tmp1
          end do
        end do
      else
        do k = kmin, kmax
          do i = imin, imax
            tmp2 = arr(i, jmin - 1, k)
            do j = jmin, jmax, 2
              tmp1 = 0.25 * arr(i, j - 1, k) + 0.5 * arr(i, j, k) + 0.25 * arr(i, j + 1, k)
              arr(i, j - 1, k) = tmp2
              tmp2 = 0.25 * arr(i, j, k) + 0.5 * arr(i, j + 1, k) + 0.25 * arr(i, j + 2, k)
              arr(i, j, k) = tmp1
            end do
            arr(i, jmax, k) = tmp2
          end do
        end do
      end if
    end do
  end subroutine filterInY

  subroutine filterFields(nsmooth)
    implicit none
    integer, intent(in) :: nsmooth
    integer :: n_pass, iter, i
    n_pass = NGHOST
    iter = 0
    do while (.true.)
      if (n_pass .ge. nsmooth) then
        if (nsmooth .gt. 0) then
          i = nsmooth - iter * NGHOST
          call filterInX(ex, i)
          call filterInX(ey, i)
          call filterInX(ez, i)
          call filterInX(bz, i)
          call filterInY(ex, i)
          call filterInY(ey, i)
          call filterInY(ez, i)
          call filterInY(bz, i)
          call exchangeFields(.true., .true.)
        end if
        exit
      else
        i = NGHOST
        call filterInX(ex, i)
        call filterInX(ey, i)
        call filterInX(ez, i)
        call filterInX(bz, i)
        call filterInY(ex, i)
        call filterInY(ey, i)
        call filterInY(ez, i)
        call filterInY(bz, i)
        call exchangeFields(.true., .true.)
        n_pass = n_pass + NGHOST
        iter = iter + 1
      end if
    end do
  end subroutine filterFields

  subroutine userInitFields()
    implicit none
    integer :: i, j
    real :: x_, y_, r1, r2, x_shift, y_shift, r1_shift, r2_shift
    real :: bx0, bz0
    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = 0; by(:, :, :) = 0; bz(:, :, :) = 0
    jx(:, :, :) = 0; jy(:, :, :) = 0; jz(:, :, :) = 0

    if (NGHOST .lt. 3) then
      call throwError('ERROR: `NGHOST` needs to be >= 3.')
    end if

    ! set tube radius to 1/4 of box size along x:
    r_j = REAL(global_mesh % sx) / 4
    alpha = 3.8317059702075  ! 1st zero of J_1(x) bessel function
    x1 = 0.5 * REAL(global_mesh % sx)
    x2 = x1
    y1 = 0.5 * REAL(global_mesh % sy) - r_j
    y2 = 0.5 * REAL(global_mesh % sy) + r_j

    ! first setup the az vector potential (temporarily stored into ey array)
    ! and the bz field
    do j = -NGHOST, this_meshblock % ptr % sy - 1 + NGHOST
      do i = -NGHOST, this_meshblock % ptr % sx - 1 + NGHOST

        x_ = REAL(i + this_meshblock % ptr % x0)
        y_ = REAL(j + this_meshblock % ptr % y0)
        x_shift = x_ + 0.5
        y_shift = y_ + 0.5
        ! 1st flux tube radius
        r1 = sqrt((x_ - x1)**2 + (y_ - y1)**2) / r_j
        r1_shift = sqrt((x_shift - x1)**2 + (y_shift - y1)**2) / r_j
        ! 2nd flux tube radius:
        r2 = sqrt((x_ - x2)**2 + (y_ - y2)**2) / r_j
        r2_shift = sqrt((x_shift - x2)**2 + (y_shift - y2)**2) / r_j

        if (r1 .lt. 1.0) then
          ey(i, j, :) = bessel_jn(0, r1 * alpha) * r_j / alpha
        else if ((r2 .lt. 1.0) .and. (.not. single_tube)) then
          ey(i, j, :) = bessel_jn(0, r2 * alpha) * r_j / alpha
        else
          ey(i, j, :) = bessel_jn(0, alpha) * r_j / alpha
        end if
        if (r1_shift .lt. 1.0) then
          bz(i, j, :) = sqrt(bessel_jn(0, r1_shift * alpha)**2 + c_param)
        else if ((r2_shift .lt. 1.0) .and. (.not. single_tube)) then
          bz(i, j, :) = sqrt(bessel_jn(0, r2_shift * alpha)**2 + c_param)
        else
          bz(i, j, :) = sqrt(bessel_jn(0, alpha)**2 + c_param)
        end if

      end do
    end do

    ! setup (the non-smooth) bx and by from the curl of vector potential:
    do j = -NGHOST, this_meshblock % ptr % sy - 1 + NGHOST - 1
      do i = -NGHOST, this_meshblock % ptr % sx - 1 + NGHOST - 1
        bx(i, j, :) = -ey(i, j, 0) + ey(i, j + 1, 0)
        by(i, j, :) = ey(i, j, 0) - ey(i + 1, j, 0)
      end do
    end do

    ! setup the motional E field for each flux tube:
    do j = -NGHOST + 1, this_meshblock % ptr % sy - 1 + NGHOST - 1
      do i = -NGHOST + 1, this_meshblock % ptr % sx - 1 + NGHOST - 1

        x_ = REAL(i + this_meshblock % ptr % x0)   ! ez coord
        y_ = REAL(j + this_meshblock % ptr % y0)   ! ez coord
        x_shift = x_ + 0.5                     ! ex coord
        y_shift = y_                           ! ex coord
        ! 1st flux tube radius
        r1 = sqrt((x_ - x1)**2 + (y_ - y1)**2) / r_j
        r1_shift = sqrt((x_shift - x1)**2 + (y_shift - y1)**2) / r_j
        ! 2nd flux tube radius:
        r2 = sqrt((x_ - x2)**2 + (y_ - y2)**2) / r_j
        r2_shift = sqrt((x_shift - x2)**2 + (y_shift - y2)**2) / r_j

        bx0 = 0.5 * (bx(i, j, 0) + bx(i, j - 1, 0)) ! bx at location of ez
        bz0 = 0.5 * (bz(i, j, 0) + bz(i, j - 1, 0)) ! bz at location of ex

        if (r1 .lt. 1.0) then
          ez(i, j, :) = beta_kick * bx0
        else if ((r2 .lt. 1.0) .and. (.not. single_tube)) then
          ez(i, j, :) = -beta_kick * bx0
        else
          ez(i, j, :) = 0.0
        end if
        if (r1_shift .lt. 1.0) then
          ex(i, j, :) = -beta_kick * bz0
        else if ((r2_shift .lt. 1.0) .and. (.not. single_tube)) then
          ex(i, j, :) = beta_kick * bz0
        else
          ex(i, j, :) = 0.0
        end if

      end do
    end do

    ! filter on the fields:
    if (smooth_flds) call filterFields(nsmooth)

    bx(:, :, :) = 0.0; by(:, :, :) = 0.0
    ! setup the smooth bx and by from the smooth vector potential:
    ! analytically B_phi = J_1( alpha * r / r_j)  (in cylindrical coords)
    ! the maximum of B_phi inside the tube is = 0.581865
    ! so, sigma_{in-plane} = 0.581865**2 sigma_{code} when
    ! based on max in-plane field.
    ! for sigma_in defined using the *mean* in-plane field
    ! inside the tube the conversion is sigma_in = 0.162215 * sigma
    do j = -NGHOST + 1, this_meshblock % ptr % sy - 1 + NGHOST - 1
      do i = -NGHOST + 1, this_meshblock % ptr % sx - 1 + NGHOST - 1
        bx(i, j, :) = -ey(i, j, 0) + ey(i, j + 1, 0)
        by(i, j, :) = ey(i, j, 0) - ey(i + 1, j, 0)
      end do
    end do
    ey(:, :, :) = 0.0
    ! now setup the current as J = CC curl(B)
    do j = -NGHOST + 1, this_meshblock % ptr % sy - 1 + NGHOST - 1
      do i = -NGHOST + 1, this_meshblock % ptr % sx - 1 + NGHOST - 1
        jx(i, j, :) = CORR * CC * (-bz(i, j - 1, 0) + bz(i, j, 0))
        jy(i, j, :) = CORR * CC * (bz(i - 1, j, 0) - bz(i, j, 0))
        jz(i, j, :) = CORR * CC * (bx(i, j - 1, 0) - bx(i, j, 0) - by(i - 1, j, 0) + by(i, j, 0))
      end do
    end do
    ! calculate the (tiny) charge density perturbation and store it into lg_arr:
    do j = -NGHOST + 1, this_meshblock % ptr % sy - 1 + NGHOST - 1
      do i = -NGHOST + 1, this_meshblock % ptr % sx - 1 + NGHOST - 1
        lg_arr(i, j, :) = ex(i, j, 0) - ex(i - 1, j, 0)
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
    ! integer :: s, ti, tj, tk, p
    ! do s = 1, nspec
    !   do ti = 1, species(s)%tile_nx
    !     do tj = 1, species(s)%tile_ny
    !       do tk = 1, species(s)%tile_nz
    !         do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
    !           ...
    !         end do
    !       end do
    !     end do
    !   end do
    ! end do
  end subroutine userDriveParticles

  subroutine userExternalFields(xp, yp, zp,&
                              & ex_ext, ey_ext, ez_ext,&
                              & bx_ext, by_ext, bz_ext)
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
  end subroutine writeUsrRestart

  subroutine readUsrRestart(rst_file)
    implicit none
    integer, intent(in) :: rst_file
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

#include "optional.F"
end module m_userfile
