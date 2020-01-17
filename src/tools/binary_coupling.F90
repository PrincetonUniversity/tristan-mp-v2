#include "../defs.F90"

module m_bincoupling
  use m_globalnamespace
  use m_aux
  use m_domain
  use m_particles
  use m_errors
  implicit none

  type :: spec_ind_pair
    ! object which identifies a particle ...
    !     ... by its species and index on a given tile
    integer :: spec, index
  end type spec_ind_pair

  type :: couple
    ! object that contains two particles ...
    !     ... given by their species and index ...
    !     ... on a single tile
    type(spec_ind_pair) :: part_1
    type(spec_ind_pair) :: part_2
  end type

  !--- PRIVATE functions -----------------------------------------!
  !...............................................................!
contains

  ! this routine pairs particles in two sets #1 and #2 randomly
  !   - two sets can either be equal, or have no intersection
  !   - sets are passed as an array of species (indices) at a given tile
  subroutine coupleParticlesOnTile(ti, tj, tk,&
                                 & sp_arr_1, n_sp_1,&
                                 & sp_arr_2, n_sp_2,&
                                 & coupled_pairs, num_couples,&
                                 & num_group_1, num_group_2)
    implicit none
    integer, intent(in)                     :: ti, tj, tk
    integer, intent(in)                     :: n_sp_1, n_sp_2 ! # of species in set #1 and #2
    integer, intent(in)                     :: sp_arr_1(n_sp_1), sp_arr_2(n_sp_2) ! array of species in set #1 and #2
    integer, intent(out)                    :: num_couples
    type(couple), allocatable, intent(out)  :: coupled_pairs(:)
    integer, optional, intent(out)          :: num_group_1, num_group_2

    ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
    !
    !
    ! THIS ROUTINE IS STILL UNDER CONSTRUCTION
    !   if you need access, please contact developers
    !
    !
    ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  end subroutine

end module m_bincoupling
