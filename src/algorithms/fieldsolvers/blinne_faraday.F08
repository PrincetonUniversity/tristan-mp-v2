#ifdef oneD
  call throwError("1D Blinne not implemented")
#elif defined(twoD)
  k = 0
  do j = 0, this_meshblock % ptr % sy - 1
    jp1 = j + 1
    jp2 = j + 2
    jm1 = j - 1
    do i = 0, this_meshblock % ptr % sx - 1
      ip1 = i + 1
      ip2 = i + 2
      im1 = i - 1
      bx(i, j, k) = bx(i, j, k) - const * (1.325 * (ez(i, jp1, k) - ez(i, j, k)) &
                                           - 0.065 * (ez(ip1, jp1, k) - ez(ip1, j, k) &
                                                      + ez(im1, jp1, k) - ez(im1, j, k) &
                                                      + ez(i, jp2, k) - ez(i, jm1, k)))
      by(i, j, k) = by(i, j, k) + const * (1.325 * (ez(ip1, j, k) - ez(i, j, k)) &
                                           - 0.065 * (ez(ip1, jp1, k) - ez(i, jp1, k) &
                                                      + ez(ip1, jm1, k) - ez(i, jm1, k) &
                                                      + ez(ip2, j, k) - ez(im1, j, k)))
      bz(i, j, k) = bz(i, j, k) - const * (1.325 * (ey(ip1, j, k) - ey(i, j, k)) &
                                           - 0.065 * (ey(ip1, jp1, k) - ey(i, jp1, k) &
                                                      + ey(ip1, jm1, k) - ey(i, jm1, k) &
                                                      + ey(ip2, j, k) - ey(im1, j, k))) &
                    + const * (1.325 * (ex(i, jp1, k) - ex(i, j, k)) &
                               - 0.065 * (ex(ip1, jp1, k) - ex(ip1, j, k) &
                                          + ex(im1, jp1, k) - ex(im1, j, k) &
                                          + ex(i, jp2, k) - ex(i, jm1, k)))
    end do
  end do
#elif defined(threeD)
  call throwError("1D Blinne not implemented")
#endif
