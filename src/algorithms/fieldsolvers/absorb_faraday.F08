#ifdef oneD
    k = 0
    zg = 0.0
    j = 0
    yg = 0.0
    do i = 0, this_meshblock % ptr % sx - 1
      ip1 = i + 1
      xg = REAL(i + this_meshblock % ptr % x0)

      lam = 0.25 * lambdaAbsorb(xg, yg, zg)
      lam1 = (1.0 + lam) / (1.0 - lam)
      bx(i, j, k) = lam1 * bx(i, j, k)

      lam = 0.25 * lambdaAbsorb(xg + 0.5, yg, zg)
      lam1 = (1.0 + lam) / (1.0 - lam)
      lam2 = 1.0 / (1.0 - lam)
      by(i, j, k) = lam1 * by(i, j, k) + lam2 * const * &
                    (ez(ip1, j, k) - ez(i, j, k))

      lam = 0.25 * lambdaAbsorb(xg + 0.5, yg, zg)
      lam1 = (1.0 + lam) / (1.0 - lam)
      lam2 = 1.0 / (1.0 - lam)
      bz(i, j, k) = lam1 * bz(i, j, k) + lam2 * const * &
                    (-ey(ip1, j, k) + ey(i, j, k))
    end do
#elif defined(twoD)
    k = 0
    zg = 0.0
    do j = 0, this_meshblock % ptr % sy - 1
      jp1 = j + 1
      yg = REAL(j + this_meshblock % ptr % y0)
      do i = 0, this_meshblock % ptr % sx - 1
        ip1 = i + 1
        xg = REAL(i + this_meshblock % ptr % x0)

        lam = 0.25 * lambdaAbsorb(xg, yg + 0.5, zg)
        lam1 = (1.0 + lam) / (1.0 - lam)
        lam2 = 1.0 / (1.0 - lam)
        bx(i, j, k) = lam1 * bx(i, j, k) + lam2 * const * &
                      (-ez(i, jp1, k) + ez(i, j, k)) - &
                      2.0 * lam * lam2 * bx_target

        lam = 0.25 * lambdaAbsorb(xg + 0.5, yg, zg)
        lam1 = (1.0 + lam) / (1.0 - lam)
        lam2 = 1.0 / (1.0 - lam)
        by(i, j, k) = lam1 * by(i, j, k) + lam2 * const * &
                      (ez(ip1, j, k) - ez(i, j, k)) - &
                      2.0 * lam * lam2 * by_target

        lam = 0.25 * lambdaAbsorb(xg + 0.5, yg + 0.5, zg)
        lam1 = (1.0 + lam) / (1.0 - lam)
        lam2 = 1.0 / (1.0 - lam)
        bz(i, j, k) = lam1 * bz(i, j, k) + lam2 * const * &
                      (ex(i, jp1, k) - ex(i, j, k) - ey(ip1, j, k) + ey(i, j, k)) - &
                      2.0 * lam * lam2 * bz_target
      end do
    end do
#elif defined(threeD)
    do k = 0, this_meshblock % ptr % sz - 1
      kp1 = k + 1
      zg = REAL(k + this_meshblock % ptr % z0)
      do j = 0, this_meshblock % ptr % sy - 1
        jp1 = j + 1
        yg = REAL(j + this_meshblock % ptr % y0)
        do i = 0, this_meshblock % ptr % sx - 1
          ip1 = i + 1
          xg = REAL(i + this_meshblock % ptr % x0)

          lam = 0.25 * lambdaAbsorb(xg, yg + 0.5, zg + 0.5)
          lam1 = (1.0 + lam) / (1.0 - lam)
          lam2 = 1.0 / (1.0 - lam)
          bx(i, j, k) = lam1 * bx(i, j, k) + lam2 * const * &
                        (ey(i, j, kp1) - ey(i, j, k) - ez(i, jp1, k) + ez(i, j, k))

          lam = 0.25 * lambdaAbsorb(xg + 0.5, yg, zg + 0.5)
          lam1 = (1.0 + lam) / (1.0 - lam)
          lam2 = 1.0 / (1.0 - lam)
          by(i, j, k) = lam1 * by(i, j, k) + lam2 * const * &
                        (ez(ip1, j, k) - ez(i, j, k) - ex(i, j, kp1) + ex(i, j, k))

          lam = 0.25 * lambdaAbsorb(xg + 0.5, yg + 0.5, zg)
          lam1 = (1.0 + lam) / (1.0 - lam)
          lam2 = 1.0 / (1.0 - lam)
          bz(i, j, k) = lam1 * bz(i, j, k) + lam2 * const * &
                        (ex(i, jp1, k) - ex(i, j, k) - ey(ip1, j, k) + ey(i, j, k))
        end do
      end do
    end do
#endif
