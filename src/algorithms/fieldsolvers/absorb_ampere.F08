#ifdef oneD
    k = 0
    zg = 0.0
    j = 0
    yg = 0.0
    do i = 0, this_meshblock % ptr % sx - 1
      im1 = i - 1
      xg = REAL(i + this_meshblock % ptr % x0)

      lam = 0.5 * lambdaAbsorb(xg + 0.5, yg, zg)
      lam1 = (1.0 + lam) / (1.0 - lam)
      ex(i, j, k) = lam1 * ex(i, j, k)

      lam = 0.5 * lambdaAbsorb(xg, yg, zg)
      lam1 = (1.0 + lam) / (1.0 - lam)
      lam2 = 1.0 / (1.0 - lam)
      ey(i, j, k) = lam1 * ey(i, j, k) + lam2 * const * &
                    (bz(im1, j, k) - bz(i, j, k))

      lam = 0.5 * lambdaAbsorb(xg, yg, zg)
      lam1 = (1.0 + lam) / (1.0 - lam)
      lam2 = 1.0 / (1.0 - lam)
      ez(i, j, k) = lam1 * ez(i, j, k) + lam2 * const * &
                    (-by(im1, j, k) + by(i, j, k))
    end do
#elif defined(twoD)
    k = 0
    zg = 0.0
    do j = 0, this_meshblock % ptr % sy - 1
      jm1 = j - 1
      yg = REAL(j + this_meshblock % ptr % y0)
      do i = 0, this_meshblock % ptr % sx - 1
        im1 = i - 1
        xg = REAL(i + this_meshblock % ptr % x0)

        lam = 0.5 * lambdaAbsorb(xg + 0.5, yg, zg)
        lam1 = (1.0 + lam) / (1.0 - lam)
        lam2 = 1.0 / (1.0 - lam)
        ex(i, j, k) = lam1 * ex(i, j, k) + lam2 * const * &
                      (-bz(i, jm1, k) + bz(i, j, k))

        lam = 0.5 * lambdaAbsorb(xg, yg + 0.5, zg)
        lam1 = (1.0 + lam) / (1.0 - lam)
        lam2 = 1.0 / (1.0 - lam)
        ey(i, j, k) = lam1 * ey(i, j, k) + lam2 * const * &
                      (bz(im1, j, k) - bz(i, j, k))

        lam = 0.5 * lambdaAbsorb(xg, yg, zg)
        lam1 = (1.0 + lam) / (1.0 - lam)
        lam2 = 1.0 / (1.0 - lam)
        ez(i, j, k) = lam1 * ez(i, j, k) + lam2 * const * &
                      (bx(i, jm1, k) - bx(i, j, k) - by(im1, j, k) + by(i, j, k))
      end do
    end do
#elif defined(threeD)
    do k = 0, this_meshblock % ptr % sz - 1
      km1 = k - 1
      zg = REAL(k + this_meshblock % ptr % z0)
      do j = 0, this_meshblock % ptr % sy - 1
        jm1 = j - 1
        yg = REAL(j + this_meshblock % ptr % y0)
        do i = 0, this_meshblock % ptr % sx - 1
          im1 = i - 1
          xg = REAL(i + this_meshblock % ptr % x0)

          lam = 0.5 * lambdaAbsorb(xg + 0.5, yg, zg)
          lam1 = (1.0 + lam) / (1.0 - lam)
          lam2 = 1.0 / (1.0 - lam)
          ex(i, j, k) = lam1 * ex(i, j, k) + lam2 * const * &
                        (by(i, j, km1) - by(i, j, k) - bz(i, jm1, k) + bz(i, j, k))

          lam = 0.5 * lambdaAbsorb(xg, yg + 0.5, zg)
          lam1 = (1.0 + lam) / (1.0 - lam)
          lam2 = 1.0 / (1.0 - lam)
          ey(i, j, k) = lam1 * ey(i, j, k) + lam2 * const * &
                        (bz(im1, j, k) - bz(i, j, k) - bx(i, j, km1) + bx(i, j, k))

          lam = 0.5 * lambdaAbsorb(xg, yg, zg + 0.5)
          lam1 = (1.0 + lam) / (1.0 - lam)
          lam2 = 1.0 / (1.0 - lam)
          ez(i, j, k) = lam1 * ez(i, j, k) + lam2 * const * &
                        (bx(i, jm1, k) - bx(i, j, k) - by(im1, j, k) + by(i, j, k))
        end do
      end do
    end do
#endif
