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

  type :: maxwellian
    ! tabulated maxwellian is used 1d/2d for any T, and 3d for T < t_crit
    ! ... in all cases if T < t_crit -- we use non-relativistic maxwellian
    ! ... and `u_table` contains `beta` instead of `4-velocity`
    real                            :: temperature, shift_gamma
    real, allocatable, dimension(:) :: DF_table
    real, allocatable, dimension(:) :: u_table
    logical                         :: generated = .false., shift_flag = .false.
    integer                         :: npoints, shift_dir
    integer                         :: dimension = 3
  end type maxwellian

  !--- PRIVATE variables -----------------------------------------!
  real    :: t_crit = 0.01
  private :: t_crit
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: tabulateMaxwellian
  private :: deallocateMaxwellian
  !...............................................................!
contains
  ! See more details in Zenitani 2015
  !   arXiv:1504.03910v1

  subroutine tabulateMaxwellian(maxw)
    implicit none
    type(maxwellian), intent(inout) :: maxw
    integer                         :: iter
    real                            :: u_1, u_2, u_max, df, temp

    if (maxw%generated) then
      call throwError('ERROR: maxwell table already generated.')
    else
      maxw%generated = .true.
    end if

    temp = maxw%temperature

    if (temp .lt. t_crit) then
      maxw%npoints = 2000
      u_max = 5 * sqrt(2 * temp)
    else
      maxw%npoints = 10000
      u_max = 50 * temp
    end if

    allocate(maxw%DF_table(maxw%npoints))
    allocate(maxw%u_table(maxw%npoints))

    do iter = 1, maxw%npoints
      u_1 = u_max * REAL(iter - 1) / REAL(maxw%npoints)
      u_2 = u_max * REAL(iter) / REAL(maxw%npoints)
      if (temp .ge. t_crit) then
        if (maxw%dimension .eq. 1) then
          df = (u_2 - u_1) * 0.5 * (exp(-sqrt(1.0 + u_1**2) / temp) + exp(-sqrt(1.0 + u_2**2) / temp))
        else if (maxw%dimension .eq. 2) then
          df = (u_2 - u_1) * 0.5 * (exp(-sqrt(1.0 + u_1**2) / temp) * u_1 + exp(-sqrt(1.0 + u_2**2) / temp) * u_2)
        else if (maxw%dimension .eq. 3) then
          df = (u_2 - u_1) * 0.5 * (exp(-sqrt(1.0 + u_1**2) / temp) * u_1**2 + exp(-sqrt(1.0 + u_2**2) / temp) * u_2**2)
        else
          call throwError('ERROR: Unknown dimension in `tabulateMaxwellian`.')
        end if
      else
        if (maxw%dimension .eq. 1) then
          df = (u_2 - u_1) * 0.5 * (exp(-u_1**2 * 0.5 / temp) + exp(-u_2**2 * 0.5 / temp))
        else if (maxw%dimension .eq. 2) then
          df = (u_2 - u_1) * 0.5 * (exp(-u_1**2 * 0.5 / temp) * u_1 + exp(-u_2**2 * 0.5 / temp) * u_2)
        else if (maxw%dimension .eq. 3) then
          df = (u_2 - u_1) * 0.5 * (exp(-u_1**2 * 0.5 / temp) * u_1**2 + exp(-u_2**2 * 0.5 / temp) * u_2**2)
        else
          call throwError('ERROR: Unknown dimension in `tabulateMaxwellian`.')
        end if
      end if
      maxw%u_table(iter) = u_2
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
    type(maxwellian), intent(inout) :: maxw
    real, intent(out)               :: u_, v_, w_
    real                            :: U, ETA, X1, X2, X3, X4, X5, X6, X7, X8, dx1, dx2, BETA, gamma, gamma1
    logical                         :: flag
    integer                         :: iter

    if ((maxw%temperature .lt. t_crit) .or. (maxw%dimension .ne. 3)) then
      ! using tabulated Maxwellian
      if (.not. maxw%generated) then
        call tabulateMaxwellian(maxw)
      end if
      X3 = random(dseed)
      do iter = 1, maxw%npoints
        if (maxw%DF_table(iter) .ge. X3) then
          if (iter .gt. 1) then
            dx1 = (maxw%DF_table(iter) - X3) / (maxw%DF_table(iter) - maxw%DF_table(iter - 1))
            dx2 = (X3 - maxw%DF_table(iter - 1)) / (maxw%DF_table(iter) - maxw%DF_table(iter - 1))
            U = maxw%u_table(iter) * dx2 +&
              & maxw%u_table(iter - 1) * dx1
          else
            dx2 = X3 / maxw%DF_table(iter)
            U = maxw%u_table(iter) * dx2
          end if
          ! if T < t_crit -> convert from `beta` to `4-velocity`
          if (maxw%temperature .lt. t_crit) then
            U = U / sqrt(1.0 - U**2)
          end if
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
    if (maxw%dimension .eq. 1) then
      X1 = random(dseed)
      u_ = SIGN(U, X1 - 0.5)
      v_ = 0.0
      w_ = 0.0
    else if (maxw%dimension .eq. 2) then
      X1 = random(dseed)
      u_ = U * sin(2.0 * M_PI * X1)
      v_ = U * cos(2.0 * M_PI * X1)
      w_ = 0.0
    else if (maxw%dimension .eq. 3) then
      X1 = random(dseed); X2 = random(dseed)
      u_ = U * (2.0 * X1 - 1.0)
      v_ = 2.0 * U * sqrt(X1 * (1.0 - X1)) * cos(2.0 * M_PI * X2)
      w_ = 2.0 * U * sqrt(X1 * (1.0 - X1)) * sin(2.0 * M_PI * X2)
    else
      call throwError('ERROR: Unknown dimension in `generateFromMaxwellian`.')
    end if

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
    if (allocated(maxw%u_table)) deallocate(maxw%u_table)
  end subroutine deallocateMaxwellian

  subroutine fillRegionWithThermalPlasma(fill_region, fill_species, num_species, ndens_sp,&
                                       & temperature, shift_gamma, shift_dir, zero_current,&
                                       & dimension, weights, spat_distr_ptr,&
                                       & dummy1, dummy2, dummy3)
    implicit none
    ! assuming that the charges of all species given in `fill_species` add up to `0`
    type(region), intent(in)          :: fill_region
    integer, intent(in)               :: num_species
    integer, intent(in)               :: fill_species(num_species)
    real, intent(in)                  :: ndens_sp, temperature
    real, optional, intent(in)        :: shift_gamma
    integer, optional, intent(in)     :: shift_dir, dimension
    type(maxwellian)                  :: fill_maxwellian
    integer                           :: num_part, n, s, spec_, dimension_
    integer(kind=2)                   :: xi_, yi_, zi_
    real                              :: fill_xmin, fill_xmax,&
                                       & fill_ymin, fill_ymax,&
                                       & fill_zmin, fill_zmax
    real                              :: u_, v_, w_, dx_, dy_, dz_
    real                              :: x_, y_, z_, rnd, num_part_r
    real                              :: x_glob, y_glob, z_glob

    real, intent(in), optional        :: weights
    real                              :: weights_
    logical, intent(in), optional     :: zero_current
    logical                           :: zero_current_

    procedure (spatialDistribution), pointer, intent(in), optional :: spat_distr_ptr
    real, intent(in), optional                                     :: dummy1, dummy2, dummy3
    real                                                           :: dummy1_, dummy2_, dummy3_

    if (.not. present(dimension)) then
      dimension_ = 3
    else
      dimension_ = dimension
    end if

    if (.not. present(zero_current)) then
      zero_current_ = .false.
    else
      zero_current_ = zero_current
    end if

    if (.not. present(weights)) then
      weights_ = 1.0
    else
      weights_ = weights
    end if

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

    fill_maxwellian%dimension = dimension_
    fill_maxwellian%generated = .false.
    if (present(shift_gamma)) then
      fill_maxwellian%shift_gamma = abs(shift_gamma)
      fill_maxwellian%shift_flag = .true.
    else
      fill_maxwellian%shift_flag = .false.
    end if

    ! global to local coordinates
    #ifdef oneD
      call globalToLocalCoords(fill_region%x_min, 0.0, 0.0,&
                             & fill_xmin, fill_ymin, fill_zmin, adjustQ = .true.)
      call globalToLocalCoords(fill_region%x_max, 0.0, 0.0,&
                             & fill_xmax, fill_ymax, fill_zmax, adjustQ = .true.)
      num_part_r = REAL(ndens_sp) * (fill_xmax - fill_xmin)
    #elif twoD
      call globalToLocalCoords(fill_region%x_min, fill_region%y_min, 0.0,&
                             & fill_xmin, fill_ymin, fill_zmin, adjustQ = .true.)
      call globalToLocalCoords(fill_region%x_max, fill_region%y_max, 0.0,&
                             & fill_xmax, fill_ymax, fill_zmax, adjustQ = .true.)
      num_part_r = REAL(ndens_sp) * (fill_xmax - fill_xmin)&
                                & * (fill_ymax - fill_ymin)
    #elif threeD
      call globalToLocalCoords(fill_region%x_min, fill_region%y_min, fill_region%z_min,&
                             & fill_xmin, fill_ymin, fill_zmin, adjustQ = .true.)
      call globalToLocalCoords(fill_region%x_max, fill_region%y_max, fill_region%z_max,&
                           & fill_xmax, fill_ymax, fill_zmax, adjustQ = .true.)
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
      if ((.not. present(spat_distr_ptr)) .or. (random(dseed) .lt. rnd)) then
        do s = 1, num_species
          ! generate momenta for every species individually
          spec_ = fill_species(s)
          if ((spec_ .le. 0) .or. (spec_ .gt. nspec)) then
            call throwError('Wrong species specified in fillRegionWithThermalPlasma.')
          end if
          !   shift direction is opposite for opposite signed species
          if (present(shift_gamma)) then
            if (zero_current_) then
              fill_maxwellian%shift_dir = shift_dir
            else
              fill_maxwellian%shift_dir = INT(SIGN(1.0, species(spec_)%ch_sp)) * shift_dir
            end if
          end if
          fill_maxwellian%temperature = temperature / species(spec_)%m_sp
          call generateFromMaxwellian(fill_maxwellian, u_, v_, w_)
          call createParticle(spec_, xi_, yi_, zi_, dx_, dy_, dz_, u_, v_, w_,&
                            & weight = weights_)
        end do
      end if
      n = n + 1
    end do
    call deallocateMaxwellian(fill_maxwellian)
  end subroutine fillRegionWithThermalPlasma
end module m_thermalplasma
