#include "../defs.F90"

module m_adaptivelb
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
  subroutine computeLoadALB(load_ALB)
    implicit none
    integer, intent(out) :: load_ALB
    integer              :: s, ti, tj, tk
    load_ALB = 0
    do s = 1, nspec
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            load_ALB = load_ALB + species(s)%prtl_tile(ti, tj, tk)%npart_sp
          end do
        end do
      end do
    end do
  end subroutine computeLoadALB

end module m_adaptivelb
