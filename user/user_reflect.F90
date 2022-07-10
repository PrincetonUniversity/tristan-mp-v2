module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
  use m_helpers, only: depositCurrentsFromSingleParticle
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  integer, private :: wall_x_location
  real, private :: shift_gamma, background_T
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'wall_x', wall_x_location)
    call getInput('problem', 'gamma_launch', shift_gamma)
    call getInput('problem', 'temperature', background_T)
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob, &
                                   dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3

    return
  end function

  function userSLBload(x_glob, y_glob, z_glob, &
                       dummy1, dummy2, dummy3)
    real :: userSLBload
    ! global coordinates
    real, intent(in), optional :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional :: dummy1, dummy2, dummy3
    return
  end function

  subroutine userInitParticles()
    implicit none
    real :: background_n
    type(region) :: back_region
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    background_n = REAL(ppc0) * 0.5
    back_region % x_min = wall_x_location + 1.0
    back_region % y_min = 0.0
    back_region % x_max = REAL(global_mesh % sx)
    back_region % y_max = REAL(global_mesh % sy)
    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, background_n, background_T, &
                                     shift_gamma=shift_gamma, shift_dir=-1, zero_current=.true.)
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = 0; by(:, :, :) = 0; bz(:, :, :) = 0
    jx(:, :, :) = 0; jy(:, :, :) = 0; jz(:, :, :) = 0
    ! ... dummy loop ...
    ! do i = 0, this_meshblock%ptr%sx - 1
    !   i_glob = i + this_meshblock%ptr%x0
    !   do j = 0, this_meshblock%ptr%sy - 1
    !     j_glob = j + this_meshblock%ptr%y0
    !     do k = 0, this_meshblock%ptr%sz - 1
    !       k_glob = k + this_meshblock%ptr%z0
    !       ...
    !     end do
    !   end do
    ! end do
  end subroutine userInitFields
  !............................................................!

  !--- driving ------------------------------------------------!
  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
    ! called after particles move and deposit ...
    ! ... and before the currents are added to the electric field
    integer :: s, ti, tj, tk, p
    real :: u_, v_, w_, x_l, y_l, z_l, x_g, inv_gamma, gamma
    real :: tfrac, xcolis, xnew, x0

    do s = 1, nspec
      do ti = 1, species(s) % tile_nx
        do tj = 1, species(s) % tile_ny
          do tk = 1, species(s) % tile_nz
            do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
              x_l = REAL(species(s) % prtl_tile(ti, tj, tk) % xi(p)) + species(s) % prtl_tile(ti, tj, tk) % dx(p)
              x_g = x_l + REAL(this_meshblock % ptr % x0)
              if (x_g .lt. wall_x_location) then
                u_ = species(s) % prtl_tile(ti, tj, tk) % u(p)
                v_ = species(s) % prtl_tile(ti, tj, tk) % v(p)
                w_ = species(s) % prtl_tile(ti, tj, tk) % w(p)
                gamma = sqrt(1.0 + u_**2 + v_**2 + w_**2)
                inv_gamma = 1.0 / gamma
                x0 = x_l - u_ * inv_gamma * CC
                tfrac = abs((x0 - wall_x_location) / (u_ * inv_gamma * CC))
                if (tfrac .lt. 1) then
                  xcolis = x0 + u_ * inv_gamma * CC * tfrac
                end if
                y_l = REAL(species(s) % prtl_tile(ti, tj, tk) % yi(p)) + species(s) % prtl_tile(ti, tj, tk) % dy(p)
                z_l = REAL(species(s) % prtl_tile(ti, tj, tk) % zi(p)) + species(s) % prtl_tile(ti, tj, tk) % dz(p)
                call depositCurrentsFromSingleParticle(s, species(s) % prtl_tile(ti, tj, tk), p, &
                                                       xcolis, y_l, z_l, x_l, y_l, z_l, -1.0)
                ! reflecting particle
                species(s) % prtl_tile(ti, tj, tk) % u(p) = -u_
                u_ = species(s) % prtl_tile(ti, tj, tk) % u(p)
                tfrac = min(abs((x_l - xcolis) / max(abs(x_l - x0), 1e-9)), 1.0)
                xnew = xcolis + u_ * inv_gamma * CC * tfrac
                species(s) % prtl_tile(ti, tj, tk) % xi(p) = INT(FLOOR(xnew), 2)
                species(s) % prtl_tile(ti, tj, tk) % dx(p) = xnew - REAL(species(s) % prtl_tile(ti, tj, tk) % xi(p))
                call depositCurrentsFromSingleParticle(s, species(s) % prtl_tile(ti, tj, tk), p, &
                                                       xcolis, y_l, z_l, xnew, y_l, z_l, 1.0)
              end if
            end do
          end do
        end do
      end do
    end do
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
    integer, optional, intent(in) :: step
    real :: shift_beta, background_n
    type(region) :: back_region

    shift_beta = sqrt(1.0 - shift_gamma**-2)
    background_n = REAL(ppc0) * 0.5
    back_region % x_min = REAL(global_mesh % sx) - CC * shift_beta
    back_region % y_min = 0.0
    back_region % x_max = REAL(global_mesh % sx)
    back_region % y_max = REAL(global_mesh % sy)
    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, background_n, background_T, &
                                     shift_gamma=shift_gamma, shift_dir=-1, zero_current=.true.)
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
    logical :: updateE_, updateB_
    integer :: i, j, k, i_glob
    real :: delta_x

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

    delta_x = 0

    if (updateE_) then
      if (this_meshblock % ptr % x0 .le. wall_x_location + delta_x) then
        do i = 0, this_meshblock % ptr % sx - 1
          i_glob = i + this_meshblock % ptr % x0
          if (i_glob .le. FLOOR(wall_x_location + delta_x)) then
            ey(i, :, :) = 0.0
            ez(i, :, :) = 0.0
          end if
        end do
      end if
    end if

    ! if (updateB_) then
    !   if (this_meshblock%ptr%x0 .le. wall_x_location + delta_x) then
    !     do i = 0, this_meshblock%ptr%sx - 1
    !       i_glob = i + this_meshblock%ptr%x0
    !       if (i_glob .le. FLOOR(wall_x_location + delta_x)) then
    !         bx(i, :, :) = 0.0
    !       end if
    !     end do
    !   end if
    ! end if
  end subroutine userFieldBoundaryConditions
  !............................................................!
end module m_userfile
