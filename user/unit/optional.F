subroutine writeUsrRestart(rst_file)
  implicit none
  integer, intent(in) :: rst_file
end subroutine writeUsrRestart

subroutine readUsrRestart(rst_file)
  implicit none
  integer, intent(in) :: rst_file
end subroutine readUsrRestart

subroutine userDeallocate()
  implicit none
end subroutine userDeallocate

elemental subroutine usrSetPhPld(u0, v0, w0, over_e_temp, &
                                 incr_pld1, incr_pld2, incr_pld3)
  !$omp declare simd(usrSetPhPld)
  real, intent(in) :: u0, v0, w0, over_e_temp
  real, intent(out) :: incr_pld1, incr_pld2, incr_pld3
  incr_pld1 = 0.0; incr_pld2 = 0.0; incr_pld3 = 0.0
end subroutine

elemental subroutine usrSetElPld(q_over_m, u0, v0, w0, over_e_temp, &
                                 ex0, ey0, ez0, bx0, by0, bz0, &
                                 incr_pld1, incr_pld2, incr_pld3)
  !$omp declare simd(usrSetElPld)
  real, intent(in) :: q_over_m, u0, v0, w0, over_e_temp, &
                      ex0, ey0, ez0, &
                      bx0, by0, bz0
  real, intent(out) :: incr_pld1, incr_pld2, incr_pld3
  incr_pld1 = 0.0; incr_pld2 = 0.0; incr_pld3 = 0.0
end subroutine
