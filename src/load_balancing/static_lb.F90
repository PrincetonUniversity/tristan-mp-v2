#include "../defs.F90"

module m_staticlb
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_helpers
  implicit none

  !--- PRIVATE variables/functions -------------------------------!

  !...............................................................!
contains
  ! the following subroutine computes the load
  !   (e.g. # of particles) for any given domain
  ! this can be arbitrary function for which the algorithm
  !   balances the domain distribution
  subroutine computeLoadSLB(load_SLB, spat_load_ptr)
    implicit none
    integer, intent(out) :: load_SLB
    integer              :: i, j, k
    real                 :: load_real
    real                 :: x_glob, y_glob, z_glob
    real                 :: sx_glob, sy_glob, sz_glob

    procedure (spatialDistribution), pointer, intent(in) :: spat_load_ptr
    
    sx_glob = REAL(global_mesh%sx)
    sy_glob = REAL(global_mesh%sy)
    sz_glob = REAL(global_mesh%sz)

    load_real = 0.0
    do i = 0, this_meshblock%ptr%sx - 1
      x_glob = REAL(i + this_meshblock%ptr%x0)
      do j = 0, this_meshblock%ptr%sy - 1
        y_glob = REAL(j + this_meshblock%ptr%y0)
        do k = 0, this_meshblock%ptr%sz - 1
          z_glob = REAL(k + this_meshblock%ptr%z0)
          load_real = load_real +&
            & spat_load_ptr(x_glob = x_glob, y_glob = y_glob, z_glob = z_glob,&
                          & dummy1 = sx_glob, dummy2 = sy_glob, dummy3 = sz_glob)
        end do
      end do
    end do
    load_SLB = INT(load_real)
  end subroutine computeLoadSLB

end module m_staticlb
