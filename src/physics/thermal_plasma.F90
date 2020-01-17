#include "../defs.F90"

module m_thermalplasma
  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_particlelogistics
  implicit none

  logical :: maxwell_generated = .false.

  type :: maxwellian
    ! tabulated maxwellian is only used with `T < 0.1 [me c^2]`
    !   so non-relativistic maxwellian is used with the parameter
    !     `beta_wave = beta / sqrt(T)` for simplicity
    real                            :: temperature, shift_gamma
    real, allocatable, dimension(:) :: DF_table
    real, allocatable, dimension(:) :: beta_table
    logical                         :: generated, shift_flag
    integer                         :: npoints, shift_dir
  end type maxwellian

  !--- PRIVATE functions -----------------------------------------!
  private :: tabulateMaxwellian, generateFromMaxwellian, deallocateMaxwellian
  !...............................................................!
contains
  ! See more details in Zenitani 2015
  !   arXiv:1504.03910v1

  ! this should only be used for `T << 1`
  subroutine tabulateMaxwellian(maxw, n)
    implicit none
    type(maxwellian), intent(inout) :: maxw
    integer, intent(in)     :: n
    integer                 :: iter
    real                    :: beta_wave1, beta_wave2, beta_wave_max, df
    ! `beta_wave` is `beta / sqrt(T)`

    if (maxw%generated) then
      call throwError('ERROR: maxwell table already generated.')
    else
      maxw%generated = .true.
    end if

    maxw%npoints = n
    allocate(maxw%DF_table(maxw%npoints))
    allocate(maxw%beta_table(maxw%npoints))

    beta_wave_max = 0.9 / sqrt(maxw%temperature)
    do iter = 1, maxw%npoints
      beta_wave1 = beta_wave_max * REAL(iter - 1) / REAL(maxw%npoints)
      beta_wave2 = beta_wave_max * REAL(iter) / REAL(maxw%npoints)
      df = (beta_wave2 - beta_wave1) * 0.5 *&
            & (exp(-beta_wave1**2 * 0.5) * beta_wave1**2 + exp(-beta_wave2**2 * 0.5) * beta_wave2**2)
      maxw%beta_table(iter) = beta_wave2
      if (iter .eq. 1) then
        maxw%DF_table(iter) = df
      else
        maxw%DF_table(iter) = maxw%DF_table(iter - 1) + df
      end if
    end do
    do iter = 1, maxw%npoints
      maxw%DF_table(iter) = maxw%DF_table(iter) / maxw%DF_table(maxw%npoints)
    end do
  end subroutine tabulateMaxwellian

  subroutine generateFromMaxwellian(maxw, u_, v_, w_)
    implicit none
    type(maxwellian), intent(in) :: maxw
    real, intent(out)            :: u_, v_, w_
    real                         :: U, ETA, X1, X2, X3, X4, X5, X6, X7, X8, dx1, dx2, BETA, gamma, gamma1
    logical                      :: flag
    integer                      :: iter
    if (maxw%temperature .lt. 0.1) then
      ! using tabulated Maxwellian
      if (.not. maxw%generated) then
        call throwError('ERROR: maxwell table not generated yet.')
      end if
      X3 = random(dseed)
      do iter = 1, maxw%npoints
        if (maxw%DF_table(iter) .ge. X3) then
          if (iter .gt. 1) then
            dx1 = (maxw%DF_table(iter) - X3) / (maxw%DF_table(iter) - maxw%DF_table(iter - 1))
            dx2 = (X3 - maxw%DF_table(iter - 1)) / (maxw%DF_table(iter) - maxw%DF_table(iter - 1))
            U = maxw%beta_table(iter) * dx2 +&
              & maxw%beta_table(iter - 1) * dx1
          else
            dx2 = X3 / maxw%DF_table(iter)
            U = maxw%beta_table(iter) * dx2
          end if
          ! `U` now is in terms of `beta_wave = beta / sqrt(T)`
          !   step 1 (convert to `beta`):
          U = U * sqrt(maxw%temperature)
           !   step 2 (convert to 4-velocity):
          U = U / sqrt(1.0 - U**2)
          exit
        end if
      end do
    else
      ! using Sobol method
      flag = .false.
      do while (.not. flag)
        X4 = random(dseed); X5 = random(dseed)
        X6 = random(dseed); X7 = random(dseed)
        if (X4 * X5 * X6 * X7 .eq. 0) cycle
        U = -maxw%temperature * log(X4 * X5 * X6)
        ETA = -maxw%temperature * log(X4 * X5 * X6 * X7)
        if (ETA**2 - U**2 .gt. 1) then
          flag = .true.
        end if
      end do
    end if

    ! generate projections
    X1 = random(dseed); X2 = random(dseed)
    u_ = U * (2.0 * X1 - 1.0)
    v_ = 2.0 * U * sqrt(X1 * (1.0 - X1)) * cos(2.0 * M_PI * X2)
    w_ = 2.0 * U * sqrt(X1 * (1.0 - X1)) * sin(2.0 * M_PI * X2)

    ! shift maxwellian
    if (maxw%shift_flag) then
      X8 = random(dseed)
      gamma = sqrt(1.0 + U**2)
      BETA = sqrt(1.0 - 1.0 / maxw%shift_gamma**2)
      select case (maxw%shift_dir)
        case (+1) ! +x
          if (-BETA * u_ / gamma .gt. X8) u_ = -u_
          u_ = maxw%shift_gamma * (u_ + BETA * sqrt(1 + U**2))
        case (-1) ! -x
          BETA = -BETA
          if (-BETA * u_ / gamma .gt. X8) u_ = -u_
          u_ = maxw%shift_gamma * (u_ + BETA * sqrt(1 + U**2))
        case (+2) ! +y
          if (-BETA * v_ / gamma .gt. X8) v_ = -v_
          v_ = maxw%shift_gamma * (v_ + BETA * sqrt(1 + U**2))
        case (-2) ! -y
          BETA = -BETA
          if (-BETA * v_ / gamma .gt. X8) v_ = -v_
          v_ = maxw%shift_gamma * (v_ + BETA * sqrt(1 + U**2))
        case (+3) ! +z
          if (-BETA * w_ / gamma .gt. X8) w_ = -w_
          w_ = maxw%shift_gamma * (w_ + BETA * sqrt(1 + U**2))
        case (-3) ! -z
          BETA = -BETA
          if (-BETA * w_ / gamma .gt. X8) w_ = -w_
          w_ = maxw%shift_gamma * (w_ + BETA * sqrt(1 + U**2))
        case default
      end select
    end if
  end subroutine generateFromMaxwellian

  subroutine deallocateMaxwellian(maxw)
    implicit none
    type(maxwellian), intent(inout) :: maxw
    if (allocated(maxw%DF_table)) deallocate(maxw%DF_table)
    if (allocated(maxw%beta_table)) deallocate(maxw%beta_table)
  end subroutine deallocateMaxwellian

  subroutine fillRegionWithThermalPlasma(fill_region, fill_species, num_species, ndens_sp,&
                                       & temperature, shift_gamma, shift_dir,&
                                       & spat_distr_ptr,&
                                       & dummy1, dummy2, dummy3)
    implicit none
    ! assuming that the charges of all species given in `fill_species` add up to `0`
    type(region), intent(in)         :: fill_region
    integer, intent(in)              :: num_species
    integer, intent(in)              :: fill_species(num_species)
    real, intent(in)                 :: ndens_sp, temperature
    real, optional, intent(in)       :: shift_gamma
    integer, optional, intent(in)    :: shift_dir
    type(maxwellian)                 :: fill_maxwellian
    integer                          :: num_part, n, s, spec_
    integer(kind=2)                  :: xi_, yi_, zi_
    real                             :: fill_xmin, fill_xmax,&
                                      & fill_ymin, fill_ymax,&
                                      & fill_zmin, fill_zmax
    real                             :: u_, v_, w_, dx_, dy_, dz_
    real                             :: x_, y_, z_, rnd, num_part_r
    real                             :: x_glob, y_glob, z_glob

    procedure (spatialDistribution), pointer, intent(in), optional :: spat_distr_ptr
    real, intent(in), optional                                     :: dummy1, dummy2, dummy3
    real                                                           :: dummy1_, dummy2_, dummy3_

    if (present(dummy1)) then
      dummy1_ = dummy1
    else
      dummy1_ = 0.0
    end if
    if (present(dummy2)) then
      dummy2_ = dummy2
    else
      dummy2_ = 0.0
    end if
    if (present(dummy3)) then
      dummy3_ = dummy3
    else
      dummy3_ = 0.0
    end if

    fill_maxwellian%temperature = temperature
    fill_maxwellian%generated = .false.
    if (present(shift_gamma)) then
      fill_maxwellian%shift_gamma = shift_gamma
      fill_maxwellian%shift_flag = .true.
    else
      fill_maxwellian%shift_flag = .false.
    end if
    if (temperature .lt. 0.1) then
      call tabulateMaxwellian(fill_maxwellian, 2000)
    end if

    ! global to local coordinates
    #ifndef threeD
      call globalToLocalCoords(fill_region%x_min, fill_region%y_min, 0.0,&
                             & fill_xmin, fill_ymin, fill_zmin, adjustQ_ = .true.)
      call globalToLocalCoords(fill_region%x_max, fill_region%y_max, 0.0,&
                             & fill_xmax, fill_ymax, fill_zmax, adjustQ_ = .true.)
    #else
      call globalToLocalCoords(fill_region%x_min, fill_region%y_min, fill_region%z_min,&
                             & fill_xmin, fill_ymin, fill_zmin, adjustQ_ = .true.)
      call globalToLocalCoords(fill_region%x_max, fill_region%y_max, fill_region%z_max,&
                             & fill_xmax, fill_ymax, fill_zmax, adjustQ_ = .true.)
    #endif

    #ifndef threeD
      num_part_r = REAL(ndens_sp) * (fill_xmax - fill_xmin)&
                                & * (fill_ymax - fill_ymin)
    #else
      num_part_r = REAL(ndens_sp) * (fill_xmax - fill_xmin)&
                                & * (fill_ymax - fill_ymin)&
                                & * (fill_zmax - fill_zmin)
    #endif
    if (num_part_r .lt. 10.0) then
      if (num_part_r .ne. 0.0) then
        num_part_r = poisson(num_part_r)
      else
        num_part_r = 0.0
      end if
    else
      num_part_r = CEILING(num_part_r)
    end if
    num_part = INT(num_part_r)

    n = 0
    do while (n .lt. num_part)
      ! generate coords for all species
      call generateCoordInRegion(fill_xmin, fill_xmax, fill_ymin, fill_ymax, fill_zmin, fill_zmax,&
                               & x_, y_, z_, xi_, yi_, zi_, dx_, dy_, dz_)

      ! if spatial distribution function is present, compute it
      !   otherwise use uniform distribution
      if (present(spat_distr_ptr)) then
        x_glob = REAL(this_meshblock%ptr%x0) + x_
        y_glob = REAL(this_meshblock%ptr%y0) + y_
        z_glob = REAL(this_meshblock%ptr%z0) + z_
        rnd = spat_distr_ptr(x_glob = x_glob, y_glob = y_glob, z_glob = z_glob,&
                           & dummy1 = dummy1_, dummy2 = dummy2_, dummy3 = dummy3_)
      else
        rnd = 1.0
      end if
      if (random(dseed) .lt. rnd) then
        do s = 1, num_species
          ! generate momenta for every species individually
          spec_ = fill_species(s)
          !   shift direction is opposite for opposite signed species
          if (present(shift_gamma)) then
            fill_maxwellian%shift_dir = INT(SIGN(1.0, species(spec_)%ch_sp)) * shift_dir
          end if
          call generateFromMaxwellian(fill_maxwellian, u_, v_, w_)
          call createParticle(spec_, xi_, yi_, zi_, dx_, dy_, dz_, u_, v_, w_)
        end do
      end if
      n = n + 1
    end do
    call deallocateMaxwellian(fill_maxwellian)
  end subroutine fillRegionWithThermalPlasma
end module m_thermalplasma
