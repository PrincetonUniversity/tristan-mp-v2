#include "../defs.F90"

!--- DOMAIN ----------------------------------------------------!
! To store all domain related variables and constructions
!   - all grid meshblocks are stored as `mesh` type
!...............................................................!

module m_domain
  use m_globalnamespace
  implicit none

  type :: meshptr
    type(mesh), pointer :: ptr
  end type meshptr

  type :: mesh
    integer             :: rnk          ! rank of the cpu that takes care of the current meshblock
    integer             :: x0, y0, z0   ! global coordinates (in cells) of the corner
    integer             :: sx, sy, sz   ! # of cells in each dimension
    ! pointers to the neighboring meshblocks
    type(meshptr), dimension(-1:1,-1:1,-1:1) :: neighbor
  end type mesh

  type :: region
    #if defined(oneD) || defined (twoD) || defined (threeD)
      real     :: x_min, x_max
    #endif

    #if defined(twoD) || defined (threeD)
      real     :: y_min, y_max
    #endif

    #if defined(threeD)
      real     :: z_min, z_max
    #endif
  end type region

  type(meshptr)                      :: this_meshblock  ! pointer to current (rank) meshblock
  type(mesh)                         :: global_mesh     ! global mesh parameters
  type(mesh), allocatable, target    :: meshblocks(:)   ! meshblocks for all cpus

  ! boundary conditions for all dimensions
  !   - boundary = 1: periodic
  !   - boundary = 0: open
  integer                            :: boundary_x, boundary_y, boundary_z
  integer                            :: sendrecv_neighbors

  ! loadbalancing variables
  type(mesh), allocatable :: new_meshblocks(:)

  ! global constants for SLB
  logical :: slb_x, slb_y, slb_z
  integer :: slb_sxmin, slb_symin, slb_szmin

  ! global constants for ALB
  logical :: alb_x, alb_y, alb_z
  integer :: alb_sxmin, alb_symin, alb_szmin
  integer :: alb_int_x, alb_int_y, alb_int_z
  integer :: alb_start_x, alb_start_y, alb_start_z

  ! load per each MPI process
  integer :: lb_load

  ! array of loads per each MPI process
  integer, allocatable :: lb_load_glob(:)
  ! array of loads per each slab in each direction
  integer, allocatable :: lb_group_x0(:), lb_group_x1(:),&
                        & lb_group_y0(:), lb_group_y1(:),&
                        & lb_group_z0(:), lb_group_z1(:)

end module m_domain
