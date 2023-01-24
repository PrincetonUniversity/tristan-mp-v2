module m_bincoupling
  ! DEP_PRT [particle-dependent]
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
    real :: wei ! records weight to allow arbitrary-weight particle splitting
  end type spec_ind_pair

  type :: couple
    ! object that contains two particles ...
    !     ... given by their species and index ...
    !     ... on a single tile
    type(spec_ind_pair) :: part_1
    type(spec_ind_pair) :: part_2
  end type couple

  !--- PRIVATE functions -----------------------------------------!
  !...............................................................!
contains
  ! to treat certain species equivalently ...
  !   ... we convert them into a "set" of `spec_ind_pair` objects
  subroutine prtlToSet(ti, tj, tk, &
                       sp_arr, n_sp, &
                       set, set_size)
    implicit none
    integer, intent(in) :: ti, tj, tk
    integer, intent(in) :: n_sp ! # of species in set
    integer, intent(in) :: sp_arr(n_sp)
    integer, intent(out) :: set_size
    type(spec_ind_pair), allocatable, intent(out) :: set(:)
    integer :: s, si, i, p

    ! computing number of particles in the set
    set_size = 0
    do si = 1, n_sp
      s = sp_arr(si)
      set_size = set_size + species(s) % prtl_tile(ti, tj, tk) % npart_sp
    end do
    allocate (set(set_size))
    ! assigning particles in the set
    i = 1
    do si = 1, n_sp
      s = sp_arr(si)
      do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
        set(i) % spec = s
        set(i) % index = p
        set(i) % wei = 1.0
        i = i + 1
      end do
    end do
  end subroutine prtlToSet

  subroutine prtlToSetWeighted(ti, tj, tk, &
                               sp_arr, n_sp, &
                               set, set_size, set_weight)
    implicit none
    integer, intent(in) :: ti, tj, tk
    integer, intent(in) :: n_sp ! # of species in set
    integer, intent(in) :: sp_arr(n_sp)
    integer, intent(out) :: set_size
    real, intent(out) :: set_weight
    type(spec_ind_pair), allocatable, intent(out) :: set(:)
    integer :: set_size_
    type(spec_ind_pair), allocatable :: set_(:)
    integer :: s, si, i, p, q
    real :: wei, wei_split

    ! computing number of particles in the set
    set_size_ = 0
    do si = 1, n_sp
      s = sp_arr(si)
      do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
        set_size_ = set_size_ + CEILING(species(s) % prtl_tile(ti, tj, tk) % weight(p))
      end do
    end do
    allocate (set_(set_size_))
    ! assigning particles in the set
    i = 1
    set_weight = 0.0
    do si = 1, n_sp
      s = sp_arr(si)
      do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
        wei = species(s) % prtl_tile(ti, tj, tk) % weight(p)
        ! evenly distribute the weight (this avoids creating tiny weight particles via splitting):
        wei_split = wei / CEILING(wei)
        set_weight = set_weight + wei
        do q = 1, (CEILING(wei))
          set_(i) % spec = s
          set_(i) % index = p
          set_(i) % wei = wei_split
          i = i + 1
        end do
      end do
    end do
    set_size = i - 1
    if (allocated(set)) deallocate (set)
    allocate (set(set_size))
    set(1:set_size) = set_(1:set_size)
  end subroutine prtlToSetWeighted

  ! Knuth's algorithm to randomly shuffle a set
  subroutine shuffleSet(set, set_size)
    implicit none
    integer, intent(in) :: set_size
    type(spec_ind_pair), intent(inout) :: set(set_size)
    type(spec_ind_pair) :: temp
    integer :: i, j
    do i = 1, set_size - 1
      j = randomInt(dseed, i, set_size + 1)
      temp = set(i)
      set(i) = set(j)
      set(j) = temp
    end do
  end subroutine shuffleSet

  ! this routine pairs particles in two sets #1 and #2 randomly
  !   - two sets can either be equal, or have no intersection
  !   - sets are passed as an array of species (indices) at a given tile
  subroutine coupleParticlesOnTile(ti, tj, tk, &
                                   sp_arr_1, n_sp_1, &
                                   sp_arr_2, n_sp_2, &
                                   coupled_pairs, num_couples, &
                                   num_group_1, num_group_2, &
                                   wei_group_1, wei_group_2)
    implicit none
    integer, intent(in) :: ti, tj, tk
    integer, intent(in) :: n_sp_1, n_sp_2 ! # of species in set #1 and #2
    integer, intent(in) :: sp_arr_1(n_sp_1), sp_arr_2(n_sp_2) ! array of species in set #1 and #2
    integer :: num_1, num_2 ! total # of particles in sets #1 and #2
    real :: wei_1, wei_2 ! total weight of sets #1 and #2
    integer :: i, j
    type(spec_ind_pair), allocatable :: set_1(:), set_2(:) ! set #1 and #2 saved as "tuples" of species and index
    integer, intent(out) :: num_couples
    type(couple), allocatable, intent(out) :: coupled_pairs(:)
    integer, optional, intent(out) :: num_group_1, num_group_2
    real, optional, intent(out) :: wei_group_1, wei_group_2

    ! auxiliary variables
    integer :: common_species
    logical :: pairing_correctQ, same_setsQ

    ! check the number of common elements
    common_species = 0
    do i = 1, n_sp_1
      do j = 1, n_sp_2
        if (sp_arr_1(i) .eq. sp_arr_2(j)) then
          common_species = common_species + 1
        end if
      end do
    end do

    if (common_species .eq. 0) then
      same_setsQ = .false.
    else
      if ((n_sp_1 .eq. n_sp_2) .and. (n_sp_1 .eq. common_species)) then
        same_setsQ = .true.
      else
        same_setsQ = .false.
        call throwError("Sets should be equal or have no intersection in `coupleParticlesOnTile()`")
      end if
    end if

    if (same_setsQ) then ! if two sets are exactly the same (e.g. gamma+gamma)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call prtlToSetWeighted(ti, tj, tk, sp_arr_1, n_sp_1, set_1, num_1, wei_1)
      ! shuffle the set
      call shuffleSet(set_1, num_1)
      num_couples = INT(num_1 / 2)
      allocate (coupled_pairs(num_couples))
      ! assign pairs
      do i = 1, num_couples
        coupled_pairs(i) % part_1 = set_1(i)
        coupled_pairs(i) % part_2 = set_1(num_1 - i + 1)
      end do
      ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    else ! if two sets have no common elements (e.g. compton scattering)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call prtlToSetWeighted(ti, tj, tk, sp_arr_1, n_sp_1, set_1, num_1, wei_1)
      call prtlToSetWeighted(ti, tj, tk, sp_arr_2, n_sp_2, set_2, num_2, wei_2)

      num_couples = min(num_1, num_2)
      if (num_couples .gt. 0) then
        allocate (coupled_pairs(max(num_1, num_2)))
        call shuffleSet(set_1, num_1)
        call shuffleSet(set_2, num_2)
        do i = 1, max(num_1, num_2)
          coupled_pairs(i) % part_1 = set_1(modulo(i-1, num_1)+1)
          coupled_pairs(i) % part_2 = set_2(modulo(i-1, num_2)+1)
        end do
      end if
      ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    end if

    if (present(num_group_1)) num_group_1 = num_1
    if (present(num_group_2)) num_group_2 = num_2
    if (present(wei_group_1)) wei_group_1 = wei_1
    if (present(wei_group_2)) wei_group_2 = wei_2

    if (allocated(set_1)) deallocate (set_1)
    if (allocated(set_2)) deallocate (set_2)
  end subroutine coupleParticlesOnTile

end module m_bincoupling
