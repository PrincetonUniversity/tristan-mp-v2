#include "../defs.F90"

module m_particles
  use m_globalnamespace
  implicit none

  type :: particle_tile
    ! DEP_PRT [particle-dependent]
    integer                                     :: npart_sp, maxptl_sp
    ! species index for a given tile
    integer                                     :: spec
    ! tile boundaries in local coordinates
    integer                                     :: x1, x2, y1, y2, z1, z2
    integer(kind=2), allocatable, dimension(:)  :: xi, yi, zi
    real, allocatable, dimension(:)             :: dx, dy, dz
    real, allocatable, dimension(:)             :: u, v, w
    real, allocatable, dimension(:)             :: weight
    integer, allocatable, dimension(:)          :: ind, proc
    !dir$ attributes align: 64 :: xi, yi, zi, dx, dy, dz, u, v, w, weight, ind, proc
    #ifdef GCA
      ! GCA specific variables
      integer(kind=2), allocatable, dimension(:)  :: xi_past, yi_past, zi_past
      real, allocatable, dimension(:)             :: dx_past, dy_past, dz_past
      real, allocatable, dimension(:)             :: u_eff, v_eff, w_eff
      real, allocatable, dimension(:)             :: u_par, u_perp
      !dir$ attributes align: 64 :: xi_past, yi_past, zi_past, dx_past, dy_past, dz_past, u_eff, v_eff, w_eff, u_par, u_perp
    #endif

    #ifdef PRTLPAYLOADS
      real, allocatable, dimension(:)             :: payload1, payload2, payload3
    #endif
    ! > `proc < 0` means the particle will be deleted once the `clearGhostParticles()` is called
  end type particle_tile

  type :: particle_species
    integer     :: cntr_sp
    real        :: m_sp, ch_sp
    ! sizes and boundaries of the tiles in each direction
    integer     :: tile_sx, tile_sy, tile_sz
    ! numbers of the tiles in each direction
    integer     :: tile_nx, tile_ny, tile_nz
    type (particle_tile), allocatable, dimension(:,:,:) :: prtl_tile
    ! `true/false` - whether this species deposit currents or not
    logical     :: deposit_sp
    ! `true/false` - whether this species moves or not
    logical     :: move_sp
    ! `true/false` - whether this species is saved into particle output
    logical     :: output_sp

    #ifdef GCA
      ! `true/false` - either this species can be treated in a GCA mover, or not
      logical     :: gca_sp
    #endif

    #ifdef DOWNSAMPLING
      ! `true/false` - either downsample species or not
      logical     :: dwn_sp
    #endif

  end type particle_species

  ! particle types for exchange between processors />
  type :: prtl_enroute
    ! DEP_PRT [particle-dependent]
    integer(kind=2)   :: xi, yi, zi

    #ifdef GCA
      integer(kind=2)   :: xi_past, yi_past, zi_past
    #endif

    real              :: dx, dy, dz

    #ifdef GCA
      real              :: dx_past, dy_past, dz_past
    #endif

    real              :: u, v, w

    #ifdef GCA
      real              :: u_eff, v_eff, w_eff
      real              :: u_perp, u_par
    #endif

    real              :: weight

    #ifdef PRTLPAYLOADS
      real              :: payload1, payload2, payload3
    #endif

    integer           :: ind, proc
  end type prtl_enroute

  type :: enroute_array
    type(prtl_enroute), allocatable     :: enroute(:)
    integer                             :: cnt
    integer                             :: max
  end type enroute_array

  type :: enroute_handler
    type(enroute_array), dimension(-1:1,-1:1,-1:1)   :: get
  end type enroute_handler
  ! </ particle types for exchange between processors

  ! main container for particles
  type(particle_species), target, allocatable   :: species(:)
  ! number of species
  integer                                       :: nspec
  integer(kind=8), allocatable                  :: maxptl_array(:)

  ! type(prtl_enroute), allocatable, dimension(:)    :: recv_enroute
  type(enroute_array)                              :: recv_enroute
  type(enroute_handler)                            :: enroute_bot

  type(enroute_array), allocatable, dimension(:)   :: prtl_backup

  #ifdef MPI08
    type(MPI_DATATYPE)                             :: myMPI_ENROUTE
  #endif

  #ifdef MPI
    integer                                        :: myMPI_ENROUTE
  #endif

end module m_particles
