module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
  use m_helpers
  use m_exchangearray, only: exchangeArray
#ifdef USROUTPUT
    use m_writeusroutput
#endif
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real, private    :: nCS_over_nUP, width, upstream_T, cs_x, cs_x1, cs_x2
  real, private    :: Rc, radloop
  real, private    :: radius, current_width, a
  integer, private :: injector_reset_interval, open_boundaries, no_cooling
  integer, private :: cs_lecs, cs_ions, cs_heavy, up_lecs, up_ions, up_heavy

  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    integer    :: istat(MPI_STATUS_SIZE)
    character(len=STR_MAX)            :: filename
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer    :: amode, file_handle, ierrmpi

    call getInput('problem', 'upstream_T', upstream_T)
    call getInput('problem', 'up_lecs', up_lecs, 1)
    call getInput('problem', 'up_ions', up_ions, 2)
    call getInput('problem', 'open_boundaries', open_boundaries, -1)
    call getInput('problem', 'RmaxOverRc', radloop, 10.0)
    Rc=REAL(global_mesh%sx)/(2*radloop)
    if (mpi_rank .eq. 1) then
       print *, "Rc=", Rc
    end if
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob,&
                                 & dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional  :: x_glob, y_glob, z_glob
    real, intent(in), optional  :: dummy1, dummy2, dummy3
    real                :: sx_glob, sy_glob
    real                        :: rad2, rad
        
    sx_glob = REAL(global_mesh%sx)
    sy_glob = REAL(global_mesh%sy)
    rad2 = (x_glob - 0.5*sx_glob)**2 + (y_glob - sy_glob)**2
    rad = sqrt(rad2)
    userSpatialDistribution = 0.3+0.7/cosh(2*rad/Rc)
    return
  end function userSpatialDistribution

  function findbphi(rad)
    real, intent(in) :: rad
    real             :: findbphi, radnorm    
    
    radnorm = rad/Rc
    findbphi = rad/Rc * exp(1-rad/Rc)
    return
  end function findbphi
  
  function findT(rad)
    real, intent(in) :: rad
    real             :: findT, radnorm, localDens

    radnorm = rad/Rc
    localDens = REAL(ppc0)*(0.3+0.7/cosh(2*radnorm))
    findT = CC * (exp(2-2*radnorm)*(1+2*radnorm-2*radnorm**2)+5*exp(-2.0)) / 4
    findT = findT / localDens

    return
  end function findT

  function findV(rad)
    real, intent(in) :: rad
    real             :: findV, radnorm, jz, localDens

    radnorm = rad/Rc
    jz = CC * 1/Rc * (2-rad/Rc) * exp(1-rad/Rc)
    localDens = REAL(ppc0)*(0.3+0.7/cosh(2*radnorm))
    findV = jz/localDens
  end function findV

  subroutine userInitParticles()
    implicit none
    type(region)        :: back_region
    integer             :: i, j
    real                :: vshift, rad2, rad, temp, shift_gamma 
    real                :: nUP_elec, nUP_pos, sx_glob, sy_glob
    procedure (spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution


    nUP_elec = 0.5 * ppc0
    nUP_pos = 0.5 * ppc0

    sx_glob = REAL(global_mesh%sx)
    sy_glob = REAL(global_mesh%sy)

    do i = 0, this_meshblock%ptr%sx - 1
       do j = 0, this_meshblock%ptr%sy - 1
          back_region%x_min = this_meshblock%ptr%x0+i
          back_region%x_max = this_meshblock%ptr%y0+j
          back_region%y_min = this_meshblock%ptr%x0+i+1
          back_region%y_max = this_meshblock%ptr%y0+j+1
          rad2=(this_meshblock%ptr%x0+i+0.5-0.5*sx_glob)**2+(this_meshblock%ptr%y0+j+0.5-0.5*sy_glob)**2
          rad=sqrt(rad2) 
          vshift=findV(rad)*sqrt(sigma)*c_omp
          if (vshift .ge. 1) then
             call throwError('ERROR: `shift_beta` >= 1 in `userInitParticles()`')
          end if
          shift_gamma = 1.0 / sqrt(1.0 - vshift**2)
          temp = findT(rad)
          call fillRegionWithThermalPlasma(back_region, (/up_lecs,up_ions/), 2, nUP_pos, temp,&
               & shift_gamma = shift_gamma, shift_dir = 3,&
               & spat_distr_ptr = spat_distr_ptr, dummy1 = (0.5 * sy_glob), dummy2 = (0.5 * sy_glob))

       end do
    end do

  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    real    :: rad2, rad, Bedge, const
    integer :: i_glob, j_glob
    real    :: bphi, cosphi, sinphi, flux, a, a2
    real    :: x_glob, sx_glob, sy_glob, y_glob
    real    :: va, vrec


    ex(:,:,:) = 0; ey(:,:,:) = 0; ez(:,:,:) = 0
    bx(:,:,:) = 0; by(:,:,:) = 0; bz(:,:,:) = 0
    jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0

    sx_glob = REAL(global_mesh%sx)
    sy_glob = REAL(global_mesh%sy)  

    do i = -NGHOST, this_meshblock%ptr%sx - 1 + NGHOST
       do j = -NGHOST, this_meshblock%ptr%sy - 1 + NGHOST  
          i_glob = i + this_meshblock%ptr%x0
          x_glob = REAL(i_glob)
          j_glob = j + this_meshblock%ptr%y0 
          y_glob = REAL(j_glob)

          ! B_x : i, j+1/2
          ! B_y : i+1/2, j
          ! B_z : i+1/2, j+1/2

          ! B_x : i, j+1/2 
          rad2 = (x_glob - sx_glob * 0.5)**2 + (y_glob + 0.5 - 0.5 * sy_glob)**2
          rad = sqrt(rad2)
          bphi = findbphi(rad)
          sinphi = (y_glob - 0.5 * (sy_glob + 0.5)) / (rad+1.0d-10)
          bx(i, j, :) = -bphi * sinphi

          ! i+1/2 j
          rad2 = (x_glob + 0.5 - sx_glob * 0.5)**2 + (y_glob - 0.5 * sy_glob)**2
          rad = sqrt(rad2)
          bphi = findbphi(rad)
          cosphi = (x_glob + 0.5 - 0.5 * sx_glob) / (rad + 1e-10)
          by(i, j, :) = bphi * cosphi
       end do
    end do
  end subroutine userInitFields
  !............................................................!

  !--- driving ------------------------------------------------!
  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
    ! called after particles move and deposit ...
    ! ... and before the currents are added to the electric field
  end subroutine userCurrentDeposit

  subroutine userDriveParticles(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userDriveParticles

  subroutine userExternalFields(xp, yp, zp,&
                              & ex_ext, ey_ext, ez_ext,&
                              & bx_ext, by_ext, bz_ext)
    implicit none
    real, intent(in)  :: xp, yp, zp
    real, intent(out) :: ex_ext, ey_ext, ez_ext
    real, intent(out) :: bx_ext, by_ext, bz_ext
    ! some functions of xp, yp, zp
    ex_ext = 0.0; ey_ext = 0.0; ez_ext = 0.0
    bx_ext = 0.0; by_ext = 0.0; bz_ext = 0.0
  end subroutine userExternalFields
  !............................................................!

  !--- boundaries ---------------------------------------------!
  subroutine userParticleBoundaryConditions(step)
    implicit none
    integer, optional, intent(in)             :: step

  end subroutine userParticleBoundaryConditions
    
  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB

  end subroutine userFieldBoundaryConditions
  !............................................................!

#include "optional.F"

  !--- user-specific output -----------------------------------!
#ifdef USROUTPUT
    subroutine userOutput(step)
      implicit none
      integer, optional, intent(in) :: step
      integer                       :: root_rank = 0
      real, allocatable             :: y_bins(:), ExB_arr(:), ExB_arr_global(:)
      ! real                          :: dr, x_glob, y_glob, z_glob, r_glob
      real                          :: dummy_x, dummy_y, dummy_z, dummy
      ! real, allocatable             :: sum_ExBr_f(:), sum_f(:), sum_ExBr_f_global(:), sum_f_global(:)
      ! integer                       :: ri, rnum = 50, i, j, k, ierr
      integer                       :: x_bin, yi, ynum, i, j, k, ierr

      ynum = INT(global_mesh%sy)

      ! allocate(y_bins(ynum))
      allocate(ExB_arr(ynum))
      allocate(ExB_arr_global(ynum))
      ExB_arr(:) = 0.0

      ! do yi = 0, ynum - 1
      !   y_bins(yi + 1) = 0.5 + REAL(yi)
      ! end do

      if (this_meshblock%ptr%x0 .eq. 0) then
        x_bin = INT(measure_x * global_mesh%sx)
        i = x_bin; k = 0
        do j = 0, this_meshblock%ptr%sy - 1
          dummy_x = -(ez(i,j,k) * by(i,j,k)) + ey(i,j,k) * bz(i,j,k)
          dummy = bx(i,j,k)**2 + by(i,j,k)**2 + bz(i,j,k)**2

          yi = j + this_meshblock%ptr%y0
          ExB_arr(yi + 1) = dummy_x / dummy
        end do
      end if

      call MPI_REDUCE(ExB_arr, ExB_arr_global, ynum, MPI_REAL, MPI_SUM, root_rank, MPI_COMM_WORLD, ierr)

      if (mpi_rank .eq. root_rank) then
        call writeUsrOutputTimestep(step)
        ! call writeUsrOutputArray('y', y_bins)
        call writeUsrOutputArray('ExB', ExB_arr_global)
        call writeUsrOutputEnd()
      end if
    end subroutine userOutput

    logical function userExcludeParticles(s, ti, tj, tk, p)
      implicit none
      integer, intent(in)       :: s, ti, tj, tk, p
      userExcludeParticles = .true.
    end function userExcludeParticles
#endif
  !............................................................!
end module m_userfile
