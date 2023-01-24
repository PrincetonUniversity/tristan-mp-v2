module m_particles
  use m_globalnamespace
  implicit none

  type :: particle_tile
    ! DEP_PRT [particle-dependent]
    integer :: npart_sp, maxptl_sp
    ! species index for a given tile
    integer :: spec
    ! tile boundaries in local coordinates
    integer :: x1, x2, y1, y2, z1, z2
    integer(kind=2), allocatable, dimension(:) :: xi, yi, zi
    real, allocatable, dimension(:) :: dx, dy, dz
    real, allocatable, dimension(:) :: u, v, w
    real, allocatable, dimension(:) :: weight
    integer, allocatable, dimension(:) :: ind, proc
    !dir$ attributes align: 64 :: xi, yi, zi, dx, dy, dz, u, v, w, weight, ind, proc

#ifdef PRTLPAYLOADS
    real, allocatable, dimension(:) :: payload1, payload2, payload3
    !dir$ attributes align: 64 :: payload1, payload2, payload3
#endif
    ! > `proc < 0` means the particle will be deleted once the `clearGhostParticles()` is called
  end type particle_tile

  type :: particle_species
    integer :: cntr_sp
    real :: m_sp, ch_sp
    ! sizes and boundaries of the tiles in each direction
    integer :: tile_sx, tile_sy, tile_sz
    ! numbers of the tiles in each direction
    integer :: tile_nx, tile_ny, tile_nz
    type(particle_tile), allocatable, dimension(:, :, :) :: prtl_tile
    ! `true/false` - whether this species deposit currents or not
    logical :: deposit_sp
    ! `true/false` - whether this species moves or not
    logical :: move_sp
    ! `true/false` - whether this species is saved into hist, fld, prtl output
    logical :: output_sp_hist
    logical :: output_sp_fld
    logical :: output_sp_prtl
    logical :: flds_at_prtl_sp
    logical :: dens_at_prtl_sp
    integer :: n_prtl_vars_sp, n_dom_vars_sp
    character(len=STR_MAX) :: prtl_vars_sp(64), prtl_var_types_sp(64)

    ! extra physics properties
#ifdef RADIATION
    ! `true/false` - either apply cooling to species or not
    logical :: cool_sp
#endif

#ifdef DOWNSAMPLING
    ! `true/false` - either downsample species or not
    logical :: dwn_sp
#endif

  end type particle_species

  ! particle types for exchange between processors />
  type :: prtl_enroute
    sequence

    ! DEP_PRT [particle-dependent]
    integer(kind=2) :: xi, yi, zi

    integer(kind=2) :: dummy2
    real :: dx, dy, dz

    real :: u, v, w

    real :: weight

#ifdef PRTLPAYLOADS
    real :: payload1, payload2, payload3
    real :: dummy3
#endif

    real :: dummy4

    integer :: ind, proc
  end type prtl_enroute

  type :: enroute_array
    type(prtl_enroute), allocatable :: enroute(:)
    integer :: cnt
    integer :: max
  end type enroute_array

  type :: enroute_handler
    type(enroute_array), dimension(-1:1, -1:1, -1:1) :: get
  end type enroute_handler
  ! </ particle types for exchange between processors

  ! main container for particles
  type(particle_species), target, allocatable :: species(:)
  ! number of species
  integer :: nspec
  integer(kind=8), allocatable :: maxptl_array(:)

  ! type(prtl_enroute), allocatable, dimension(:)    :: recv_enroute
  type(enroute_array) :: recv_enroute
  type(enroute_handler) :: enroute_bot

  type(enroute_array), allocatable, dimension(:) :: prtl_backup

#ifdef MPI08
  type(MPI_DATATYPE) :: myMPI_ENROUTE
#endif

#ifdef MPI
  integer :: myMPI_ENROUTE
#endif

end module m_particles
