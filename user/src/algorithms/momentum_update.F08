! this "function" takes
! ... the field quantities: `bx0`, `by0`, `bz0`, `ex0`, `ey0`, `ez0` ...
! ... and the velocities: `u0`, `v0`, `w0` ...
! ... and returns the updated velocities `u0`, `v0`, `w0`
u0 = pt_u(p); v0 = pt_v(p); w0 = pt_w(p)
#ifndef VAY
#include "boris_push.F08"
#else
#include "vay_push.F08"
#endif
pt_u(p) = u0; pt_v(p) = v0; pt_w(p) = w0
! RADIATION >
#ifdef RADIATION
#ifdef SYNCHROTRON
if (species(s) % cool_sp) then
  call particleRadiateSync(timestep, s, &
                           pt_u(p), pt_v(p), pt_w(p), u_init, v_init, w_init, &
                           pt_dx(p), pt_dy(p), pt_dz(p), pt_xi(p), pt_yi(p), pt_zi(p), pt_wei(p), &
                           bx_rad, by_rad, bz_rad, ex_rad, ey_rad, ez_rad, &
                           index=pt_ind(p), proc=pt_proc(p))
end if
#endif
#ifdef INVERSECOMPTON
if (species(s) % cool_sp) then
  call particleRadiateIC(timestep, s, &
                         pt_u(p), pt_v(p), pt_w(p), u_init, v_init, w_init, &
                         pt_dx(p), pt_dy(p), pt_dz(p), pt_xi(p), pt_yi(p), pt_zi(p), pt_wei(p), &
                         bx_rad, by_rad, bz_rad, ex_rad, ey_rad, ez_rad, &
                         index=pt_ind(p), proc=pt_proc(p))
end if
#endif
#endif
! </ RADIATION
