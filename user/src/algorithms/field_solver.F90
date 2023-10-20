module m_fldsolver
  use m_globalnamespace
  use m_aux
  use m_domain
  use m_fields
  implicit none

  !--- PRIVATE functions -----------------------------------------!
  private :: lambdaAbsorb
  !...............................................................!
contains
  subroutine advanceBHalfstep()
    implicit none
    integer :: i, j, k, ip1, jp1, kp1
    real :: const
#ifdef ABSORB
    real :: lam, lam1, lam2, xg, yg, zg
#endif

#ifdef BLINNE
    integer :: ip2, jp2, im1, jm1
#endif
    const = CORR * 0.5 * CC

#ifndef ABSORB
#  ifndef BLINNE
#    include "default_faraday.F08"
#  else
#    include "blinne_faraday.F08"
#  endif
#else
#  include "absorb_faraday.F08"
#endif
    call printDiag("advanceBHalfstep()", 2)
  end subroutine advanceBHalfstep

  subroutine advanceEFullstep()
    implicit none
    integer :: i, j, k, im1, jm1, km1
    real :: const
#ifdef ABSORB
    real :: lam, lam1, lam2, xg, yg, zg
#endif
    const = CORR * CC

#ifndef ABSORB
#  include "default_ampere.F08"
#else
#  include "absorb_ampere.F08"
#endif

    call printDiag("advanceEFullstep()", 2)
  end subroutine advanceEFullstep

  subroutine addCurrents()
    implicit none
    integer :: xmin, xmax, ymin, ymax, zmin, zmax
    xmin = 0
    ymin = 0
    zmin = 0

#ifdef oneD
    xmax = this_meshblock % ptr % sx - 1
    ymax = 0
    zmax = 0
#elif defined(twoD)
    xmax = this_meshblock % ptr % sx - 1
    ymax = this_meshblock % ptr % sy - 1
    zmax = 0
#elif defined(threeD)
    xmax = this_meshblock % ptr % sx - 1
    ymax = this_meshblock % ptr % sy - 1
    zmax = this_meshblock % ptr % sz - 1
#endif
    ! "-" sign is taken care of in the deposit
    ex(xmin:xmax, ymin:ymax, zmin:zmax) = &
      ex(xmin:xmax, ymin:ymax, zmin:zmax) + &
      jx(xmin:xmax, ymin:ymax, zmin:zmax)
    ey(xmin:xmax, ymin:ymax, zmin:zmax) = &
      ey(xmin:xmax, ymin:ymax, zmin:zmax) + &
      jy(xmin:xmax, ymin:ymax, zmin:zmax)
    ez(xmin:xmax, ymin:ymax, zmin:zmax) = &
      ez(xmin:xmax, ymin:ymax, zmin:zmax) + &
      jz(xmin:xmax, ymin:ymax, zmin:zmax)
    call printDiag("addCurrents()", 2)
  end subroutine addCurrents

  real function lambdaAbsorb(x0, y0, z0)
    implicit none
    real, intent(in) :: x0, y0, z0 ! global coordinates
    real :: K_abs
    real :: gr_max, gr_bound, gc_x, gc_y, gc_z, radius
    lambdaAbsorb = 0.0
    K_abs = CC / 3.0
    if (absorb_x .eq. 2) then
      ! radial open boundaries (in all directions)
#ifdef oneD
      gc_x = global_mesh % sx * 0.5
      gr_max = gc_x
      gr_bound = (gr_max - ds_abs)**2
      radius = (x0 - gc_x)**2
#elif defined(twoD)
      gc_x = global_mesh % sx * 0.5
      gc_y = global_mesh % sy * 0.5
      gr_max = MIN(gc_x, gc_y)
      gr_bound = (gr_max - ds_abs)**2
      radius = (x0 - gc_x)**2 + (y0 - gc_y)**2
#elif defined(threeD)
      gc_x = global_mesh % sx * 0.5
      gc_y = global_mesh % sy * 0.5
      gc_z = global_mesh % sz * 0.5
      gr_max = MIN(gc_x, gc_y, gc_z)
      gr_bound = (gr_max - ds_abs)**2
      radius = (x0 - gc_x)**2 + (y0 - gc_y)**2 + (z0 - gc_z)**2
#endif
      if (radius .gt. gr_bound) then
        radius = sqrt(radius)
        gr_bound = sqrt(gr_bound)
        lambdaAbsorb = -MIN(K_abs * ((radius - gr_bound) / ds_abs)**3, K_abs)
      end if
    else
      ! simple cartesian open boundaries
#if defined(oneD) || defined (twoD) || defined (threeD)
      if (absorb_x .eq. 1) then
        ! open boundaries in x direction
        if (x0 .lt. ds_abs) then
          lambdaAbsorb = -K_abs * ((ds_abs - x0) / ds_abs)**3
        else if (x0 .gt. global_mesh % sx - ds_abs) then
          lambdaAbsorb = -K_abs * ((x0 - (global_mesh % sx - 1.0 - ds_abs)) / ds_abs)**3
        end if
      end if
#endif
#if defined(twoD) || defined (threeD)
      if (absorb_y .eq. 1) then
        ! open boundaries in y direction
        if (y0 .lt. ds_abs) then
          lambdaAbsorb = -K_abs * ((ds_abs - y0) / ds_abs)**3
        else if (y0 .gt. global_mesh % sy - ds_abs) then
          lambdaAbsorb = -K_abs * ((y0 - (global_mesh % sy - 1.0 - ds_abs)) / ds_abs)**3
        end if
      end if
#endif
#if defined(threeD)
      if (absorb_z .eq. 0) then
        ! open boundaries in z direction
        if (z0 .lt. ds_abs) then
          lambdaAbsorb = -K_abs * ((ds_abs - z0) / ds_abs)**3
        else if (z0 .gt. global_mesh % sz - ds_abs) then
          lambdaAbsorb = -K_abs * ((z0 - (global_mesh % sz - 1.0 - ds_abs)) / ds_abs)**3
        end if
      end if
#endif
    end if
  end function lambdaAbsorb
end module m_fldsolver

