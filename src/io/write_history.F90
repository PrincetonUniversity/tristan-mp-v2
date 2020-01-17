#include "../defs.F90"

module m_writehistory
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_helpers
  implicit none
contains
  ! FIX2 total E^2, total B^2, total E_kin
  ! time, time step, total mass, momenta, kinetic energies in three directions, total energy, and magnetic energies in three directions
  subroutine writeHistory()
    implicit none
  end subroutine writeHistory

end module m_writehistory
