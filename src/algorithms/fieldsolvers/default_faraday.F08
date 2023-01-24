#ifdef oneD
    k = 0
    j = 0
    do i = 0, this_meshblock % ptr % sx - 1
      ip1 = i + 1
      by(i, j, k) = by(i, j, k) + const * &
                    (ez(ip1, j, k) - ez(i, j, k))
      bz(i, j, k) = bz(i, j, k) + const * &
                    (-ey(ip1, j, k) + ey(i, j, k))
    end do
#elif defined(twoD)
    k = 0
    do j = 0, this_meshblock % ptr % sy - 1
      jp1 = j + 1
      do i = 0, this_meshblock % ptr % sx - 1
        ip1 = i + 1
        bx(i, j, k) = bx(i, j, k) + const * &
                      (-ez(i, jp1, k) + ez(i, j, k))
        by(i, j, k) = by(i, j, k) + const * &
                      (ez(ip1, j, k) - ez(i, j, k))
        bz(i, j, k) = bz(i, j, k) + const * &
                      (ex(i, jp1, k) - ex(i, j, k) - ey(ip1, j, k) + ey(i, j, k))
      end do
    end do
#elif defined(threeD)
    do k = 0, this_meshblock % ptr % sz - 1
      kp1 = k + 1
      do j = 0, this_meshblock % ptr % sy - 1
        jp1 = j + 1
        do i = 0, this_meshblock % ptr % sx - 1
          ip1 = i + 1
          bx(i, j, k) = bx(i, j, k) + const * &
                        (ey(i, j, kp1) - ey(i, j, k) - ez(i, jp1, k) + ez(i, j, k))
          by(i, j, k) = by(i, j, k) + const * &
                        (ez(ip1, j, k) - ez(i, j, k) - ex(i, j, kp1) + ex(i, j, k))
          bz(i, j, k) = bz(i, j, k) + const * &
                        (ex(i, jp1, k) - ex(i, j, k) - ey(ip1, j, k) + ey(i, j, k))
        end do
      end do
    end do
#endif
