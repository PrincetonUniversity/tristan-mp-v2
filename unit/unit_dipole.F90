#include "../src/defs.F90"

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
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  integer, private  :: fld_geometry
  real, private     :: xc_g, yc_g, zc_g, radius
  real, private     :: posx, posy, posz, velocity, angle
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution, randomPointInSphericalShell
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    ! B-field geometry: 1 = monopole, 2 = dipole
    call getInput('problem', 'fld_geometry', fld_geometry)
    call getInput('problem', 'radius', radius)
    call getInput('problem', 'posx', posx)
    call getInput('problem', 'posy', posy)
    call getInput('problem', 'posz', posz, 0.5)
    call getInput('problem', 'angle', angle)
    call getInput('problem', 'velocity', velocity)

    xc_g = 0.5 * global_mesh%sx
    yc_g = 0.5 * global_mesh%sy
    zc_g = 0.5 * global_mesh%sz
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob,&
                                 & dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional  :: x_glob, y_glob, z_glob
    real, intent(in), optional  :: dummy1, dummy2, dummy3

    return
  end function

  function userSLBload(x_glob, y_glob, z_glob,&
                     & dummy1, dummy2, dummy3)
    real :: userSLBload
    ! global coordinates
    real, intent(in), optional  :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional  :: dummy1, dummy2, dummy3
    real                        :: radius2
    radius2 = (dummy1 * 0.5 - x_glob)**2 + (dummy2 * 0.5 - y_glob)**2 + (dummy3 * 0.5 - z_glob)**2 + 1.0
    userSLBload = 40**2 / radius2
    if (radius2 .lt. 40**2) then
      userSLBload = 1.0 / exp((40**2 - radius2) / 40**2)
    end if
    return
  end function

  subroutine userInitParticles()
    implicit none
    procedure (spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    real    :: bx0, by0, bz0
    ex(:,:,:) = 0; ey(:,:,:) = 0; ez(:,:,:) = 0
    bx(:,:,:) = 0; by(:,:,:) = 0; bz(:,:,:) = 0
    jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0

    #if defined(twoD)
      k = 0
      do i = 0, this_meshblock%ptr%sx - 1
        i_glob = i + this_meshblock%ptr%x0
        do j = 0, this_meshblock%ptr%sy - 1
          j_glob = j + this_meshblock%ptr%y0
          call getBfield(REAL(i_glob), REAL(j_glob) + 0.5, REAL(k), bx0, by0, bz0)
          bx(i, j, k) = bx0
          call getBfield(REAL(i_glob) + 0.5, REAL(j_glob), REAL(k), bx0, by0, bz0)
          by(i, j, k) = by0
          call getBfield(REAL(i_glob) + 0.5, REAL(j_glob) + 0.5, REAL(k), bx0, by0, bz0)
          bz(i, j, k) = bz0
        end do
      end do
    #elif defined(threeD)
      do i = 0, this_meshblock%ptr%sx - 1
        i_glob = i + this_meshblock%ptr%x0
        do j = 0, this_meshblock%ptr%sy - 1
          j_glob = j + this_meshblock%ptr%y0
          do k = 0, this_meshblock%ptr%sz - 1
            k_glob = k + this_meshblock%ptr%z0
            call getBfield(REAL(i_glob), REAL(j_glob) + 0.5, REAL(k_glob) + 0.5, bx0, by0, bz0)
            bx(i, j, k) = bx0
            call getBfield(REAL(i_glob) + 0.5, REAL(j_glob), REAL(k_glob) + 0.5, bx0, by0, bz0)
            by(i, j, k) = by0
            call getBfield(REAL(i_glob) + 0.5, REAL(j_glob) + 0.5, REAL(k_glob), bx0, by0, bz0)
            bz(i, j, k) = bz0
          end do
        end do
      end do
    #endif
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
    ! ... dummy loop ...
    ! integer :: s, ti, tj, tk, p
    ! do s = 1, nspec
    !   do ti = 1, species(s)%tile_nx
    !     do tj = 1, species(s)%tile_ny
    !       do tk = 1, species(s)%tile_nz
    !         do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
    !           ...
    !         end do
    !       end do
    !     end do
    !   end do
    ! end do
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
    integer, optional, intent(in) :: step
    integer                       :: s, ti, tj, tk, p
    real                          :: x_g, y_g, z_g, r_g
    real                          :: Ux, Uy, Uz

    #if defined(twoD)
      Ux = velocity * cos(angle * M_PI / 180.0)
      Uy = velocity * sin(angle * M_PI / 180.0)
      Uz = 0.0
    #elif defined(threeD)
      Ux = velocity * cos(angle * M_PI / 180.0)
      Uy = 0.0
      Uz = velocity * sin(angle * M_PI / 180.0)
    #endif
    if (step .eq. 0) call injectParticleGlobally(1, posx, posy, posz, Ux, Uy, Uz)
    if (step .eq. 0) call injectParticleGlobally(2, posx, posy, posz, Ux, Uy, Uz)

    ! remove particles falling into the star
    do s = 1, nspec
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              #if defined (twoD)
                x_g = REAL(species(s)%prtl_tile(ti, tj, tk)%xi(p) + this_meshblock%ptr%x0)&
                    & + species(s)%prtl_tile(ti, tj, tk)%dx(p)
                y_g = REAL(species(s)%prtl_tile(ti, tj, tk)%yi(p) + this_meshblock%ptr%y0)&
                    & + species(s)%prtl_tile(ti, tj, tk)%dy(p)
                r_g = sqrt((x_g - xc_g)**2 + (y_g - yc_g)**2)
                if (r_g .lt. (radius - 2.0)) then
                  species(s)%prtl_tile(ti, tj, tk)%proc(p) = -1
                end if
              #elif defined(threeD)
                x_g = REAL(species(s)%prtl_tile(ti, tj, tk)%xi(p) + this_meshblock%ptr%x0)&
                    & + species(s)%prtl_tile(ti, tj, tk)%dx(p)
                y_g = REAL(species(s)%prtl_tile(ti, tj, tk)%yi(p) + this_meshblock%ptr%y0)&
                    & + species(s)%prtl_tile(ti, tj, tk)%dy(p)
                z_g = REAL(species(s)%prtl_tile(ti, tj, tk)%zi(p) + this_meshblock%ptr%z0)&
                    & + species(s)%prtl_tile(ti, tj, tk)%dz(p)
                r_g = sqrt((x_g - xc_g)**2 + (y_g - yc_g)**2 + (z_g - zc_g)**2)
                if (r_g .lt. (radius - 2.0)) then
                  species(s)%prtl_tile(ti, tj, tk)%proc(p) = -1
                end if
              #endif
            end do
          end do
        end do
      end do
    end do
  end subroutine userParticleBoundaryConditions

  subroutine randomPointInSphericalShell(rmin, rmax, x, y, z)
    implicit none
    real, intent(in)    :: rmin, rmax
    real, intent(out)   :: x, y, z
    real                :: TH, X0, Y0, Z0, T, R
    TH = random(dseed) * 2.0 * M_PI
    Z0 = 2.0 * (random(dseed) - 0.5)
    X0 = sqrt(1.0 - Z0**2) * cos(TH)
    Y0 = sqrt(1.0 - Z0**2) * sin(TH)
    T = random(dseed) * (rmax**3 - rmin**3) + rmin**3
    R = T**(1.0 / 3.0)
    x = X0 * R
    y = Y0 * R
    z = Z0 * R
  end subroutine randomPointInSphericalShell

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
    logical                       :: updateE_, updateB_
    integer                       :: i, j, k, i_glob, j_glob, k_glob, mx, my, mz
    real                          :: supersph_radius_sq
    real, allocatable             :: ex_new(:,:,:), ey_new(:,:,:), ez_new(:,:,:),&
                                   & bx_new(:,:,:), by_new(:,:,:), bz_new(:,:,:)
    real                          :: rx, ry, rz, rr_sqr, shift_B, s
    real                          :: bx_dip, by_dip, bz_dip, b_dip_dot_r, b_int_dot_r
    real                          :: scaleEpar, scaleEperp, scaleBperp, scaleBpar, scale
    real                          :: vx, vy, vz, ex_dip, ey_dip, ez_dip
    real                          :: shift_E, e_int_dot_r, e_dip_dot_r

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
  end subroutine userFieldBoundaryConditions
  !............................................................!

  !--- auxiliary functions ------------------------------------!
  subroutine getBfield(x_g, y_g, z_g,&
                     & obx, oby, obz)
    real, intent(in)    :: x_g, y_g, z_g
    real, intent(out)   :: obx, oby, obz
    if (fld_geometry .eq. 1) then
      call getMonopole(x_g, y_g, z_g, obx, oby, obz)
    else if (fld_geometry .eq. 2) then
      call getDipole(x_g, y_g, z_g, obx, oby, obz)
    else
      print *, "Something went wrong in `usr_psr`."
      stop
    end if
  end subroutine getBfield

  subroutine getDipole(x_g, y_g, z_g,&
                     & obx, oby, obz)
    implicit none
    real, intent(in)    :: x_g, y_g, z_g
    real, intent(out)   :: obx, oby, obz
    real                :: phase, nx, ny, nz, rr, mux, muy, muz, mu_dot_n

    nx = x_g - xc_g
    ny = y_g - yc_g
    nz = z_g - zc_g

    rr = sqrt(nx**2 + ny**2 + nz**2)
    nx = nx / rr
    ny = ny / rr
    nz = nz / rr

    rr = 1.0 / rr**3

    #if defined(twoD)
      mux = 0.0
      muy = radius**3
      muz = 0.0
    #elif defined(threeD)
      mux = 0.0
      muy = 0.0
      muz = radius**3
    #endif

    mu_dot_n = mux * nx + muy * ny + muz * nz

    obx = (3.0 * nx * mu_dot_n - mux) * rr
    oby = (3.0 * ny * mu_dot_n - muy) * rr
    obz = (3.0 * nz * mu_dot_n - muz) * rr
  end subroutine getDipole

  subroutine getMonopole(x_g, y_g, z_g,&
                       & obx, oby, obz)
    implicit none
    real, intent(in)    :: x_g, y_g, z_g
    real, intent(out)   :: obx, oby, obz
    real                :: nx, ny, nz, rr
    nx = x_g - xc_g
    ny = y_g - yc_g
    nz = z_g - zc_g

    rr = sqrt(nx**2 + ny**2 + nz**2)
    rr = 1.0 / rr**3

    obx = radius**2 * nx * rr
    oby = radius**2 * ny * rr
    obz = radius**2 * nz * rr
  end subroutine getMonopole

  real function shape(rad, rad0)
    implicit none
    real, intent(in)  :: rad, rad0
    real              :: del
    del = 1.0
    shape = 0.5 * (1.0 - tanh((rad - rad0) / del))
  end function shape

  !............................................................!
end module m_userfile
