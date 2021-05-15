#include "../defs.F90"

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
    real :: lam, lam1, lam2, xg, yg, zg
    const = CORR * 0.5 * CC

    #ifndef ABSORB

      #ifdef oneD
        k = 0
        j = 0
        do i = 0, this_meshblock%ptr%sx - 1
          ip1 = i + 1
          by(i, j, k) = by(i, j, k) + const *&
                    & (ez(ip1, j, k) - ez(i, j, k))
          bz(i, j, k) = bz(i, j, k) + const *&
                    & (-ey(ip1, j, k) + ey(i, j, k))
        enddo
      #elif twoD
        k = 0
        do j = 0, this_meshblock%ptr%sy - 1
          jp1 = j + 1
          do i = 0, this_meshblock%ptr%sx - 1
            ip1 = i + 1
            bx(i, j, k) = bx(i, j, k) + const *&
                      & (-ez(i, jp1, k) + ez(i, j, k))
            by(i, j, k) = by(i, j, k) + const *&
                      & (ez(ip1, j, k) - ez(i, j, k))
            bz(i, j, k) = bz(i, j, k) + const *&
                      & (ex(i, jp1, k) - ex(i, j, k) - ey(ip1, j, k) + ey(i, j, k))
          enddo
        enddo
      #elif threeD
        do k = 0, this_meshblock%ptr%sz - 1
          kp1 = k + 1
          do j = 0, this_meshblock%ptr%sy - 1
            jp1 = j + 1
            do i = 0, this_meshblock%ptr%sx - 1
              ip1 = i + 1
              bx(i, j, k) = bx(i, j, k) + const *&
                        & (ey(i, j, kp1) - ey(i, j, k) - ez(i, jp1, k) + ez(i, j, k))
              by(i, j, k) = by(i, j, k) + const *&
                        & (ez(ip1, j, k) - ez(i, j, k) - ex(i, j, kp1) + ex(i, j, k))
              bz(i, j, k) = bz(i, j, k) + const *&
                        & (ex(i, jp1, k) - ex(i, j, k) - ey(ip1, j, k) + ey(i, j, k))
            enddo
          enddo
        enddo
      #endif

    #else

      #ifdef oneD
        k = 0
        zg = 0.0
        j = 0
        yg = 0.0
        do i = 0, this_meshblock%ptr%sx - 1
          ip1 = i + 1
          xg = REAL(i + this_meshblock%ptr%x0)

          lam = 0.25 * lambdaAbsorb(xg, yg, zg)
          lam1 = (1.0 + lam) / (1.0 - lam)
          bx(i, j, k) = lam1 * bx(i, j, k)

          lam = 0.25 * lambdaAbsorb(xg + 0.5, yg, zg)
          lam1 = (1.0 + lam) / (1.0 - lam)
          lam2 = 1.0 / (1.0 - lam)
          by(i, j, k) = lam1 * by(i, j, k) + lam2 * const *&
                    & (ez(ip1, j, k) - ez(i, j, k))

          lam = 0.25 * lambdaAbsorb(xg + 0.5, yg, zg)
          lam1 = (1.0 + lam) / (1.0 - lam)
          lam2 = 1.0 / (1.0 - lam)
          bz(i, j, k) = lam1 * bz(i, j, k) + lam2 * const *&
                    & (-ey(ip1, j, k) + ey(i, j, k))
        enddo
      #elif twoD
        k = 0
        zg = 0.0
        do j = 0, this_meshblock%ptr%sy - 1
          jp1 = j + 1
          yg = REAL(j + this_meshblock%ptr%y0)
          do i = 0, this_meshblock%ptr%sx - 1
            ip1 = i + 1
            xg = REAL(i + this_meshblock%ptr%x0)

            lam = 0.25 * lambdaAbsorb(xg, yg + 0.5, zg)
            lam1 = (1.0 + lam) / (1.0 - lam)
            lam2 = 1.0 / (1.0 - lam)
            bx(i, j, k) = lam1 * bx(i, j, k) + lam2 * const *&
                      & (-ez(i, jp1, k) + ez(i, j, k))

            lam = 0.25 * lambdaAbsorb(xg + 0.5, yg, zg)
            lam1 = (1.0 + lam) / (1.0 - lam)
            lam2 = 1.0 / (1.0 - lam)
            by(i, j, k) = lam1 * by(i, j, k) + lam2 * const *&
                      & (ez(ip1, j, k) - ez(i, j, k))

            lam = 0.25 * lambdaAbsorb(xg + 0.5, yg + 0.5, zg)
            lam1 = (1.0 + lam) / (1.0 - lam)
            lam2 = 1.0 / (1.0 - lam)
            bz(i, j, k) = lam1 * bz(i, j, k) + lam2 * const *&
                      & (ex(i, jp1, k) - ex(i, j, k) - ey(ip1, j, k) + ey(i, j, k))
          enddo
        enddo
      #elif threeD
        do k = 0, this_meshblock%ptr%sz - 1
          kp1 = k + 1
          zg = REAL(k + this_meshblock%ptr%z0)
          do j = 0, this_meshblock%ptr%sy - 1
            jp1 = j + 1
            yg = REAL(j + this_meshblock%ptr%y0)
            do i = 0, this_meshblock%ptr%sx - 1
              ip1 = i + 1
              xg = REAL(i + this_meshblock%ptr%x0)

              lam = 0.25 * lambdaAbsorb(xg, yg + 0.5, zg + 0.5)
              lam1 = (1.0 + lam) / (1.0 - lam)
              lam2 = 1.0 / (1.0 - lam)
              bx(i, j, k) = lam1 * bx(i, j, k) + lam2 * const *&
                        & (ey(i, j, kp1) - ey(i, j, k) - ez(i, jp1, k) + ez(i, j, k))

              lam = 0.25 * lambdaAbsorb(xg + 0.5, yg, zg + 0.5)
              lam1 = (1.0 + lam) / (1.0 - lam)
              lam2 = 1.0 / (1.0 - lam)
              by(i, j, k) = lam1 * by(i, j, k) + lam2 * const *&
                        & (ez(ip1, j, k) - ez(i, j, k) - ex(i, j, kp1) + ex(i, j, k))

              lam = 0.25 * lambdaAbsorb(xg + 0.5, yg + 0.5, zg)
              lam1 = (1.0 + lam) / (1.0 - lam)
              lam2 = 1.0 / (1.0 - lam)
              bz(i, j, k) = lam1 * bz(i, j, k) + lam2 * const *&
                        & (ex(i, jp1, k) - ex(i, j, k) - ey(ip1, j, k) + ey(i, j, k))
            enddo
          enddo
        enddo
      #endif

    #endif
    call printDiag((mpi_rank .eq. 0), "advanceBHalfstep()", .true.)
  end subroutine advanceBHalfstep

  subroutine advanceEFullstep()
    implicit none
    integer :: i, j, k, im1, jm1, km1
    real :: const
    real :: lam, lam1, lam2, xg, yg, zg
    const = CORR * CC

    #ifndef ABSORB

      #ifdef oneD
        k = 0
        j = 0
        do i = 0, this_meshblock%ptr%sx - 1
          im1 = i - 1
          ey(i, j, k) = ey(i, j, k) + const *&
                    & (bz(im1, j, k) - bz(i, j, k))
          ez(i, j, k) = ez(i, j, k) + const *&
                    & (-by(im1, j, k) + by(i, j, k))
        enddo
      #elif twoD
        k = 0
        do j = 0, this_meshblock%ptr%sy - 1
          jm1 = j - 1
          do i = 0, this_meshblock%ptr%sx - 1
            im1 = i - 1
            ex(i, j, k) = ex(i, j, k) + const *&
                      & (-bz(i, jm1, k) + bz(i, j, k))
            ey(i, j, k) = ey(i, j, k) + const *&
                      & (bz(im1, j, k) - bz(i, j, k))
            ez(i, j, k) = ez(i, j, k) + const *&
                      & (bx(i, jm1, k) - bx(i, j, k) - by(im1, j, k) + by(i, j, k))
          enddo
        enddo
      #elif threeD
        do k = 0, this_meshblock%ptr%sz - 1
          km1 = k - 1
          do j = 0, this_meshblock%ptr%sy - 1
            jm1 = j - 1
            do i = 0, this_meshblock%ptr%sx - 1
              im1 = i - 1
              ex(i, j, k) = ex(i, j, k) + const *&
                        & (by(i, j, km1) - by(i, j, k) - bz(i, jm1, k) + bz(i, j, k))
              ey(i, j, k) = ey(i, j, k) + const *&
                        & (bz(im1, j, k) - bz(i, j, k) - bx(i, j, km1) + bx(i, j, k))
              ez(i, j, k) = ez(i, j, k) + const *&
                        & (bx(i, jm1, k) - bx(i, j, k) - by(im1, j, k) + by(i, j, k))
            enddo
          enddo
        enddo
      #endif

    #else

      #ifdef oneD
        k = 0
        zg = 0.0
        j = 0
        yg = 0.0
        do i = 0, this_meshblock%ptr%sx - 1
          im1 = i - 1
          xg = REAL(i + this_meshblock%ptr%x0)

          lam = 0.5 * lambdaAbsorb(xg + 0.5, yg, zg)
          lam1 = (1.0 + lam) / (1.0 - lam)
          ex(i, j, k) = lam1 * ex(i, j, k)

          lam = 0.5 * lambdaAbsorb(xg, yg, zg)
          lam1 = (1.0 + lam) / (1.0 - lam)
          lam2 = 1.0 / (1.0 - lam)
          ey(i, j, k) = lam1 * ey(i, j, k) + lam2 * const *&
                    & (bz(im1, j, k) - bz(i, j, k))

          lam = 0.5 * lambdaAbsorb(xg, yg, zg)
          lam1 = (1.0 + lam) / (1.0 - lam)
          lam2 = 1.0 / (1.0 - lam)
          ez(i, j, k) = lam1 * ez(i, j, k) + lam2 * const *&
                    & (-by(im1, j, k) + by(i, j, k))
        enddo
      #elif twoD
        k = 0
        zg = 0.0
        do j = 0, this_meshblock%ptr%sy - 1
          jm1 = j - 1
          yg = REAL(j + this_meshblock%ptr%y0)
          do i = 0, this_meshblock%ptr%sx - 1
            im1 = i - 1
            xg = REAL(i + this_meshblock%ptr%x0)

            lam = 0.5 * lambdaAbsorb(xg + 0.5, yg, zg)
            lam1 = (1.0 + lam) / (1.0 - lam)
            lam2 = 1.0 / (1.0 - lam)
            ex(i, j, k) = lam1 * ex(i, j, k) + lam2 * const *&
                      & (-bz(i, jm1, k) + bz(i, j, k))

            lam = 0.5 * lambdaAbsorb(xg, yg + 0.5, zg)
            lam1 = (1.0 + lam) / (1.0 - lam)
            lam2 = 1.0 / (1.0 - lam)
            ey(i, j, k) = lam1 * ey(i, j, k) + lam2 * const *&
                      & (bz(im1, j, k) - bz(i, j, k))

            lam = 0.5 * lambdaAbsorb(xg, yg, zg)
            lam1 = (1.0 + lam) / (1.0 - lam)
            lam2 = 1.0 / (1.0 - lam)
            ez(i, j, k) = lam1 * ez(i, j, k) + lam2 * const *&
                      & (bx(i, jm1, k) - bx(i, j, k) - by(im1, j, k) + by(i, j, k))
          enddo
        enddo
      #elif threeD
        do k = 0, this_meshblock%ptr%sz - 1
          km1 = k - 1
          zg = REAL(k + this_meshblock%ptr%z0)
          do j = 0, this_meshblock%ptr%sy - 1
            jm1 = j - 1
            yg = REAL(j + this_meshblock%ptr%y0)
            do i = 0, this_meshblock%ptr%sx - 1
              im1 = i - 1
              xg = REAL(i + this_meshblock%ptr%x0)

              lam = 0.5 * lambdaAbsorb(xg + 0.5, yg, zg)
              lam1 = (1.0 + lam) / (1.0 - lam)
              lam2 = 1.0 / (1.0 - lam)
              ex(i, j, k) = lam1 * ex(i, j, k) + lam2 * const *&
                        & (by(i, j, km1) - by(i, j, k) - bz(i, jm1, k) + bz(i, j, k))

              lam = 0.5 * lambdaAbsorb(xg, yg + 0.5, zg)
              lam1 = (1.0 + lam) / (1.0 - lam)
              lam2 = 1.0 / (1.0 - lam)
              ey(i, j, k) = lam1 * ey(i, j, k) + lam2 * const *&
                        & (bz(im1, j, k) - bz(i, j, k) - bx(i, j, km1) + bx(i, j, k))

              lam = 0.5 * lambdaAbsorb(xg, yg, zg + 0.5)
              lam1 = (1.0 + lam) / (1.0 - lam)
              lam2 = 1.0 / (1.0 - lam)
              ez(i, j, k) = lam1 * ez(i, j, k) + lam2 * const *&
                        & (bx(i, jm1, k) - bx(i, j, k) - by(im1, j, k) + by(i, j, k))
            enddo
          enddo
        enddo
      #endif

    #endif
    call printDiag((mpi_rank .eq. 0), "advanceEFullstep()", .true.)
  end subroutine advanceEFullstep

  subroutine addCurrents()
    implicit none
    integer :: xmin, xmax, ymin, ymax, zmin, zmax
    xmin = 0
    ymin = 0
    zmin = 0

    #ifdef oneD
      xmax = this_meshblock%ptr%sx - 1
      ymax = 0
      zmax = 0
    #elif twoD
      xmax = this_meshblock%ptr%sx - 1
      ymax = this_meshblock%ptr%sy - 1
      zmax = 0
    #elif threeD
      xmax = this_meshblock%ptr%sx - 1
      ymax = this_meshblock%ptr%sy - 1
      zmax = this_meshblock%ptr%sz - 1
    #endif
    ! "-" sign is taken care of in the deposit
    ex(xmin:xmax, ymin:ymax, zmin:zmax) = &
        & ex(xmin:xmax, ymin:ymax, zmin:zmax) + &
        & jx(xmin:xmax, ymin:ymax, zmin:zmax)
    ey(xmin:xmax, ymin:ymax, zmin:zmax) = &
        & ey(xmin:xmax, ymin:ymax, zmin:zmax) + &
        & jy(xmin:xmax, ymin:ymax, zmin:zmax)
    ez(xmin:xmax, ymin:ymax, zmin:zmax) = &
        & ez(xmin:xmax, ymin:ymax, zmin:zmax) + &
        & jz(xmin:xmax, ymin:ymax, zmin:zmax)
    call printDiag((mpi_rank .eq. 0), "addCurrents()", .true.)
  end subroutine addCurrents

  real function lambdaAbsorb(x0, y0, z0)
    implicit none
    real, intent(in)  :: x0, y0, z0 ! global coordinates
    real              :: K_abs
    real              :: gr_max, gr_bound, gc_x, gc_y, gc_z, radius
    lambdaAbsorb = 0.0
    K_abs = CC / 3.0
    if (boundary_x .eq. 2) then
      ! radial open boundaries (in all directions)
      #ifdef oneD
        gc_x = global_mesh%sx * 0.5
        gr_max = gc_x
        gr_bound = (gr_max - ds_abs)**2
        radius = (x0 - gc_x)**2
      #elif twoD
        gc_x = global_mesh%sx * 0.5
        gc_y = global_mesh%sy * 0.5
        gr_max = MIN(gc_x, gc_y)
        gr_bound = (gr_max - ds_abs)**2
        radius = (x0 - gc_x)**2 + (y0 - gc_y)**2
      #elif threeD
        gc_x = global_mesh%sx * 0.5
        gc_y = global_mesh%sy * 0.5
        gc_z = global_mesh%sz * 0.5
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
        if (boundary_x .eq. 0) then
          ! open boundaries in x direction
          if (x0 .lt. ds_abs) then
            lambdaAbsorb = -K_abs * ((ds_abs - x0) / ds_abs)**3
          else if (x0 .gt. global_mesh%sx - ds_abs) then
            lambdaAbsorb = -K_abs * ((x0 - (global_mesh%sx - 1.0 - ds_abs)) / ds_abs)**3
          end if
        end if
      #endif

      #if defined(twoD) || defined (threeD)
        if (boundary_y .eq. 0) then
          ! open boundaries in y direction
          if (y0 .lt. ds_abs) then
            lambdaAbsorb = -K_abs * ((ds_abs - y0) / ds_abs)**3
          else if (y0 .gt. global_mesh%sy - ds_abs) then
            lambdaAbsorb = -K_abs * ((y0 - (global_mesh%sy - 1.0 - ds_abs)) / ds_abs)**3
          end if
        end if
      #endif

      #if defined(threeD)
        if (boundary_z .eq. 0) then
          ! open boundaries in z direction
          if (z0 .lt. ds_abs) then
            lambdaAbsorb = -K_abs * ((ds_abs - z0) / ds_abs)**3
          else if (z0 .gt. global_mesh%sz - ds_abs) then
            lambdaAbsorb = -K_abs * ((z0 - (global_mesh%sz - 1.0 - ds_abs)) / ds_abs)**3
          end if
        end if
      #endif
    end if
  end function lambdaAbsorb
end module m_fldsolver
