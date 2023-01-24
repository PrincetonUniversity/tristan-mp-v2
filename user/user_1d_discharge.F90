module m_userfile
  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
  use m_exchangefields
  use m_exchangeparts
#ifdef USROUTPUT
  use m_writeusroutput
#endif
  use m_qednamespace

  implicit none

  !--- PRIVATE variables -----------------------------------------!

  real, private :: contrast, temp, chi
  real, private :: force, jGJ
  real, private :: xGJ, x05GJ, ThEndGrav, xEndGrav
  real(kind=8), private    :: injector_x1
  integer, private :: starWidth, height, open_boundaries
  integer :: ncellsInjector
  logical, private :: movingInj
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput() 
    implicit none
    real :: Pmax
    call getInput('problem', 'contrast', contrast, 10.0)
    call getInput('problem', 'height', height, 10)
    call getInput('problem', 'temperature', temp, 1e-2)
    call getInput('problem', 'width_star', starWidth, 10)
    call getInput('problem', 'open_boundaries', open_boundaries, -1)
    call getInput('problem', 'nInjector', ncellsInjector, height) 
    call getInput('problem', 'chi', chi, 0.5) !! j_magnetosphere / j_GJ
    force = temp * CC**2 / height !! gravitational force
    jGJ = -CC**3 / c_omp**2 / B_norm!! goldreich-julian current
    xGJ = height*log(contrast)+starWidth
    x05GJ = height*log(2*contrast)+starWidth
    xEndGrav = (x05GJ+xGJ)/2
    ThEndGrav = (x05GJ-xGJ)/4
    if (mpi_rank .eq. 1) then
       print *, "injector", ncellsInjector
    end if
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob, &
                                   dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional :: x_glob,  y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3
    real ::glob_sx, center, width

    if (x_glob>=starWidth) then
       userSpatialDistribution = exp(-(x_glob-REAL(starWidth))/REAL(height))
    else
       userSpatialDistribution = 0
    end if
    
    return
  end function userSpatialDistribution

  function userCloudDistribution(x_glob, y_glob, z_glob, &
         dummy1, dummy2, dummy3)
    real :: userCloudDistribution
    real, intent(in), optional :: x_glob,  y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3
    real ::glob_sx, center, width

    glob_sx=REAL(global_mesh % sx)
    center = glob_sx/2
    width=100.0
    userCloudDistribution = exp(-(x_glob-center)**2/width**2)
    
    return
  end function userCloudDistribution
  
  subroutine userDeallocate()
    implicit none
  end subroutine userDeallocate

#ifdef PRTLPAYLOADS    
  subroutine userPayload(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userPayload
#endif
  
  function userSLBload(x_glob, y_glob, z_glob, &
                       dummy1, dummy2, dummy3)
    real :: userSLBload
    ! global coordinates
    real, intent(in), optional :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional :: dummy1, dummy2, dummy3
    return
  end function

#ifdef PRTLPAYLOADS
  elemental subroutine usrSetPhPld(u0, v0, w0, over_e_temp, incr_pld1, incr_pld2, incr_pld3, incr_pld4, incr_pld5)
    !$omp declare simd(usrSetPhPld)
    real, intent(in) :: u0, v0, w0, over_e_temp
    real, intent(out) :: incr_pld1, incr_pld2, incr_pld3, incr_pld4, incr_pld5
    incr_pld1 = 0.0
    incr_pld2 = 0.0
    incr_pld3 = 0.0
    incr_pld4 = 0.0
    incr_pld5 = 0.0
  end subroutine
  elemental subroutine usrSetElPld(q_over_m, u0, v0, w0, over_e_temp, ex0, &
                                   ey0, ez0, bx0, by0, bz0, incr_pld1, incr_pld2, incr_pld3, incr_pld4, incr_pld5)
    !$omp declare simd(usrSetElPld)
    real, intent(in) :: q_over_m, u0, v0, w0, over_e_temp, ex0, ey0, ez0, bx0, by0, bz0
    real, intent(out) :: incr_pld1, incr_pld2, incr_pld3, incr_pld4, incr_pld5
    incr_pld1 = 0.0
    incr_pld2 = 0.0
    incr_pld3 = 0.0
    incr_pld4 = 0.0
    incr_pld5 = 0.0
  end subroutine
#endif
  
  subroutine userInitParticles()
    implicit none
    type(region) :: back_region
    real :: baseDensity, chargeDensity, charge, x_glob
    integer :: i, j
    integer :: s, ti, tj, tk, p
    real :: LocalWeight
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()

    
    spat_distr_ptr => userSpatialDistribution
    baseDensity = 0.5 * contrast * ppc0

    back_region % x_min = REAL(starWidth) 
    back_region % x_max = REAL(global_mesh % sx)
    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, baseDensity, &
         & temp, spat_distr_ptr = spat_distr_ptr,dimension=1)

    do i = 1, this_meshblock%ptr%sx
       charge = ex(i,0,0) - ex(i-1,0,0)
       x_glob = REAL(i + this_meshblock%ptr%x0)
       back_region % x_min = x_glob-0.5
       back_region % x_max = x_glob+0.5
       LocalWeight = charge*B_norm/(ppc0*unit_ch)-1.0
       if (abs(LocalWeight) .gt. 1e-6) then
          if (LocalWeight > 0) then
             do j=1, int(ppc0)
                call createParticle(2, int(i-1,2), int(0,2), int(0,2), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, weight=LocalWeight)
             end do
          else
             do j=1, int(ppc0)
                call createParticle(1, int(i-1,2), int(0,2), int(0,2), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, weight=-LocalWeight)
             end do
          end if
       end if
    end do

  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, i_glob
    real :: x_glob
    ey(:,:,:) = 0; ez(:,:,:) = 0
    
    bx(:,:,:) = 1; by(:,:,:) = 0; bz(:,:,:) = 0
    
    jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0

#ifdef CURVATURE
    curvature(:,:,:) = 0
    BabsCurv(:, :, :) = 0
#endif
    
    do i = -NGHOST, this_meshblock%ptr%sx - 1 + NGHOST    
       i_glob = i + this_meshblock%ptr%x0
       x_glob = REAL(i + this_meshblock%ptr%x0)+0.5
       if (i_glob>=starWidth) then
          if (((x_glob-xGJ)/(0.5*height)) < 5) then
             ex(i,:,:) = 0.5*x_glob-0.25*height*log(cosh(-xGJ/(0.5*height)))+0.25*height*log(cosh((x_glob-xGJ)/(0.5*height)))
          else
             ex(i,:,:) = 0.5*x_glob-0.25*height*log(cosh(-xGJ/(0.5*height)))+0.25*height*((x_glob-xGJ)/(0.5*height)-log(2.0))
          end if
       else
          ex(i,:,:) = 0.0
       end if
    end do

    ex(:,:,:) = unit_ch * ppc0 * ex(:,:,:) / B_norm
  end subroutine userInitFields
  !............................................................!

  !--- driving ------------------------------------------------!
  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
    integer :: i, j, k
    real :: x_glob
    integer :: s, ti, tj, tk, p
    real :: glob_sx
    
    glob_sx = REAL(global_mesh % sx)
    do k = 0, this_meshblock % ptr % sz - 1
       do j = 0, this_meshblock % ptr % sy - 1
          do i = 0, this_meshblock % ptr % sx - 1
             x_glob=REAL(i+this_meshblock%ptr%x0)
             jx(i, j, k) = jx(i, j, k) + chi * jGJ
             if (x_glob > glob_sx-ds_abs) then
                if (chi<1.0 .and. chi>0.0) then
                   jx(i, j, k) = jx(i, j, k) + jGJ*(1.01-chi)*(x_glob - glob_sx + ds_abs)/ds_abs
                end if
             end if
          end do
       end do
    end do

#ifdef  PRTLPAYLOADS
    do s = 1, nspec
       do ti = 1, species(s)%tile_nx
          do tj = 1, species(s)%tile_ny
             do tk = 1, species(s)%tile_nz
                do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
                   x_glob=REAL(species(s)%prtl_tile(ti, tj, tk)%xi(p) + this_meshblock%ptr%x0)&
                        & + species(s)%prtl_tile(ti, tj, tk)%dx(p)
                   species(s)%prtl_tile(ti, tj, tk)%payload2(p) = chi
                   if (x_glob > glob_sx-ds_abs) then
                      if (chi<1.0 .and. chi>0.0) then
                         species(s)%prtl_tile(ti, tj, tk)%payload2(p) = species(s)%prtl_tile(ti, tj, tk)%payload2(p) + (1.01-chi)*(x_glob - glob_sx + ds_abs)/ds_abs
                      end if
                   end if
                end do
             end do
          end do
       end do
    end do
#endif
    
  end subroutine userCurrentDeposit

  subroutine userDriveParticles(step)
    implicit none
    integer, optional, intent(in) :: step
    integer :: s, ti, tj, tk, p, x_glob
    do s = 1, nspec
       do ti = 1, species(s)%tile_nx
          do tj = 1, species(s)%tile_ny
             do tk = 1, species(s)%tile_nz
                do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
                   x_glob=REAL(species(s)%prtl_tile(ti, tj, tk)%xi(p) + this_meshblock%ptr%x0)&
                        & + species(s)%prtl_tile(ti, tj, tk)%dx(p)
#ifdef  PRTLPAYLOADS
                   species(s)%prtl_tile(ti, tj, tk)%payload1(p) = (0.25*tanh((x_glob-0.7*starWidth)/(0.15*starWidth))+0.5-0.25*tanh((x_glob-xEndGrav)/ThEndGrav))
#endif
                   species(s)%prtl_tile(ti, tj, tk)%u(p) = species(s)%prtl_tile(ti, tj, tk)%u(p)&
                        - force* (0.25*tanh((x_glob-0.7*starWidth)/(0.15*starWidth))+0.5-0.25*tanh((x_glob-xEndGrav)/ThEndGrav))                   
                end do
             end do
          end do
       end do
    end do
    call printDiag("userDriveParticles()", 2)
  end subroutine userDriveParticles

  subroutine userExternalFields(xp, yp, zp, &
                                ex_ext, ey_ext, ez_ext, &
                                bx_ext, by_ext, bz_ext)
    implicit none
    real, intent(in) :: xp, yp, zp
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
    type(region) :: back_region
    integer, optional, intent(in) :: step
    integer :: s, ti, tj, tk, p, i
    real :: xcur, partsneeded
    integer :: index, i_glob
    real :: x_glob
    integer :: ierr
    integer :: root_rank = 0
    real :: glob_sx
    real :: baseDensity
    real, dimension(ncellsInjector) :: density, lackOfDensity, densityGlobal
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    procedure(spatialCloudDistribution), pointer :: spat_cloud_distr_ptr => null()

    spat_cloud_distr_ptr => userCloudDistribution
    
    density(:) = 0
    lackOfDensity(:) = 0
    densityGlobal(:) = 0 
    
    
    do s = 1, nspec
       do ti = 1, species(s)%tile_nx
          do tj = 1, species(s)%tile_ny
             do tk = 1, species(s)%tile_nz
                do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
                   i_glob = species(s)%prtl_tile(ti, tj, tk)%xi(p) + this_meshblock%ptr%x0
                   if (i_glob .ge. starWidth .and. i_glob .lt. starWidth+ncellsInjector) then
                      index=i_glob-starWidth+1
                      density(index)=density(index)+1
                   end if
                end do
             end do
          end do
        end do
    end do

    call MPI_REDUCE(density, densityGlobal, ncellsInjector, MPI_INTEGER, MPI_SUM, root_rank, MPI_COMM_WORLD, ierr)

    call MPI_BCAST(densityGlobal, ncellsInjector, MPI_INTEGER, root_rank, MPI_COMM_WORLD, ierr)

    do i=1,ncellsInjector
       x_glob=REAL(i+starWidth-1, 8)
       lackOfDensity(i) = ppc0*contrast*exp(-(x_glob-REAL(starWidth,8))/REAL(height,8)) - densityGlobal(i)
    end do

    do i=1,ncellsInjector
       if (lackOfDensity(i)>0) then
          x_glob=REAL(i+starWidth-1,8)
          partsneeded = lackOfDensity(i)/2.0
          back_region % x_min = REAL(x_glob, 8)
          back_region % x_max = REAL(x_glob+1.0, 8)
          call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, partsneeded, temp,dimension=1)
       end if
    end do
       
    if(step .ge. 2000 .and. (int(step/1000) .ne. int((step-200)/1000))) then
       glob_sx=REAL(global_mesh % sx)       
       back_region % x_min = glob_sx/2.0-200
       back_region % x_max = glob_sx/2.0+200
       baseDensity=0.2*ppc0*contrast/200
       if (mpi_rank .eq. 1) then
          print *, "call cloud"
       end if
       call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, baseDensity, &
            & temp, spat_distr_ptr = spat_cloud_distr_ptr,dimension=1, weights=-1.0)
       do s = 1, nspec
          do ti = 1, species(s)%tile_nx
             do tj = 1, species(s)%tile_ny
                do tk = 1, species(s)%tile_nz
                   do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp                      
                      if (species(s)%prtl_tile(ti, tj, tk)%weight(p) < 0) then
                         species(s)%prtl_tile(ti, tj, tk)%weight(p) = 1 
                         species(s)%prtl_tile(ti, tj, tk)%u(p) = species(s)%prtl_tile(ti, tj, tk)%u(p) + 0.5
                      end if
                   end do
                end do
             end do
          end do
       end do
    end if
    call printDiag("userParticlesBCend()", 2)
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
    integer :: i_glob, i
    logical :: updateE_, updateB_
    real :: glob_sx
    call printDiag("userFieldBCbegin()", 2) 
    if ((step .ge. open_boundaries) .and. (open_boundaries .ge. 0)) then
      boundary_x = 0
    end if
    if (step .eq. open_boundaries) then
      call reassignNeighborsForAll(meshblocks)
    end if

    
    if (present(updateE)) then
       updateE_ = updateE
    else
       updateE_ = .true.
    end if
    
    if (present(updateB)) then
       updateB_ = updateB
    else
       updateB_ = .true.
    end if

    if(updateB_) then
       do i = -NGHOST, this_meshblock%ptr%sx - 1 + NGHOST
          i_glob = i + this_meshblock%ptr%x0
          if (i_glob < starWidth) then
             by(i,:,:)=0.0
             bz(i,:,:)=0.0
          end if
          if (i_glob <= starWidth) then
             bx(i,:,:)=1.0
          end if
       end do
    end if

    if(updateE_) then
       do i = -NGHOST, this_meshblock%ptr%sx - 1 + NGHOST
          i_glob = i + this_meshblock%ptr%x0
          if (i_glob < starWidth) then
             ex(i,:,:)=0.0
          end if
          
          if (i_glob <= starWidth) then
             ey(i,:,:)=0.0
             ez(i,:,:)=0.0
          end if
       end do
    end if

    glob_sx = global_mesh % sx
    if(updateB_) then
       do i = -NGHOST, this_meshblock%ptr%sx - 1 + NGHOST
          i_glob = i + this_meshblock%ptr%x0
          if (i_glob > (glob_sx-10)) then
             bx(i,:,:)=1.0
          end if
       end do
    end if

    call printDiag("userFieldBCend()", 2) 
  end subroutine userFieldBoundaryConditions
  !............................................................!

  subroutine writeUsrRestart(rst_file)
    implicit none
    integer, intent(in) :: rst_file
  end subroutine writeUsrRestart

  subroutine readUsrRestart(rst_file)
    implicit none
    integer, intent(in) :: rst_file
    integer :: i1, i2, j1, j2, k1, k2
    character(len=STR_MAX) :: filename, mpichar
  end subroutine readUsrRestart

  !--- user-specific output -----------------------------------!
#ifdef USROUTPUT
  subroutine userOutput(step)
    implicit none
    integer, optional, intent(in) :: step
    ! ...
  end subroutine userOutput

  logical function userExcludeParticles(s, ti, tj, tk, p)
    implicit none
    integer, intent(in) :: s, ti, tj, tk, p
    userExcludeParticles = .true.
  end function userExcludeParticles
#endif
  !............................................................!

end module m_userfile

