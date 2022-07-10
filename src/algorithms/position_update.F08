! this "function" takes
! ... inverse energy: `over_e_temp` ...
! ... reads the velocities from: `pt_*(p)` ...
! ... and updates the particle position `pt_*(p)`

#if defined(oneD) || defined(twoD) || defined(threeD)
  pt_dx(p) = pt_dx(p) + CC * pt_u(p) * over_e_temp
  temp_i = INT(pt_dx(p), 2)
  temp_r = MAX(SIGN(1.0, pt_dx(p)) + temp_i, REAL(temp_i)) - 1
  temp_i = INT(temp_r, 2)
  pt_xi(p) = pt_xi(p) + temp_i
  pt_dx(p) = pt_dx(p) - temp_r
#endif

#if defined(twoD) || defined(threeD)
  pt_dy(p) = pt_dy(p) + CC * pt_v(p) * over_e_temp
  temp_i = INT(pt_dy(p), 2)
  temp_r = MAX(SIGN(1.0, pt_dy(p)) + temp_i, REAL(temp_i)) - 1
  temp_i = INT(temp_r, 2)
  pt_yi(p) = pt_yi(p) + temp_i
  pt_dy(p) = pt_dy(p) - temp_r
#endif

#if defined(threeD)
  pt_dz(p) = pt_dz(p) + CC * pt_w(p) * over_e_temp
  temp_i = INT(pt_dz(p), 2)
  temp_r = MAX(SIGN(1.0, pt_dz(p)) + temp_i, REAL(temp_i)) - 1
  temp_i = INT(temp_r, 2)
  pt_zi(p) = pt_zi(p) + temp_i
  pt_dz(p) = pt_dz(p) - temp_r
#endif
