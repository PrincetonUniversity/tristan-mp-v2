#ifdef oneD
    k = 0
    j = 0
    do i = 0, this_meshblock % ptr % sx - 1
      im1 = i - 1
      ey(i, j, k) = ey(i, j, k) + const * &
                    (bz(im1, j, k) - bz(i, j, k))
      ez(i, j, k) = ez(i, j, k) + const * &
                    (-by(im1, j, k) + by(i, j, k))
    end do
#elif defined(twoD)
    k = 0
    do j = 0, this_meshblock % ptr % sy - 1
      jm1 = j - 1
      do i = 0, this_meshblock % ptr % sx - 1
        im1 = i - 1
        ex(i, j, k) = ex(i, j, k) + const * &
                      (-bz(i, jm1, k) + bz(i, j, k))
        ey(i, j, k) = ey(i, j, k) + const * &
                      (bz(im1, j, k) - bz(i, j, k))
        ez(i, j, k) = ez(i, j, k) + const * &
                      (bx(i, jm1, k) - bx(i, j, k) - by(im1, j, k) + by(i, j, k))
      end do
    end do
#elif defined(threeD)
    do k = 0, this_meshblock % ptr % sz - 1
      km1 = k - 1
      do j = 0, this_meshblock % ptr % sy - 1
        jm1 = j - 1
        do i = 0, this_meshblock % ptr % sx - 1
          im1 = i - 1
          ex(i, j, k) = ex(i, j, k) + const * &
                        (by(i, j, km1) - by(i, j, k) - bz(i, jm1, k) + bz(i, j, k))
          ey(i, j, k) = ey(i, j, k) + const * &
                        (bz(im1, j, k) - bz(i, j, k) - bx(i, j, km1) + bx(i, j, k))
          ez(i, j, k) = ez(i, j, k) + const * &
                        (bx(i, jm1, k) - bx(i, j, k) - by(im1, j, k) + by(i, j, k))
        end do
      end do
    end do
#endif
