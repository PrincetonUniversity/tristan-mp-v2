! Configuration for this userfile:
! ```
!   $ python configure.py --nghosts=5 --user=unit_bw -qed -bwpp
! ```

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
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real, private :: background_T, photon_to_lepton, photon_energy, photon_mfp, photon_weight
  logical, private :: photon_outflow
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'temperature', background_T)
    call getInput('problem', 'photon_to_lepton', photon_to_lepton)
    call getInput('problem', 'photon_energy', photon_energy)
    call getInput('problem', 'photon_outflow', photon_outflow)
    call getInput('problem', 'photon_mfp', photon_mfp)
    call getInput('problem', 'photon_weight', photon_weight)
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
    type(region) :: back_region
    integer :: n, npart, s, ti, tj, tk
    integer(kind=2) :: xi_, yi_, zi_
    real :: x_, y_, z_, u_, v_, w_, vol, dx_, dy_, dz_, U_rand, TH_rand
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    back_region % x_min = 0.0
    back_region % x_max = REAL(global_mesh % sx)
#if defined(twoD) || defined (threeD)
    back_region % y_min = 0.0
    back_region % y_max = REAL(global_mesh % sy)
#endif
#if defined(threeD)
    back_region % z_min = 0.0
    back_region % z_max = REAL(global_mesh % sz)
#endif
    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, ppc0, background_T)

    if (.not. photon_outflow) then
      ! initializing photons if not outflow boundaries
      vol = (back_region % x_max - back_region % x_min) * (back_region % y_max - back_region % y_min)
      npart = vol * ppc0 * photon_to_lepton / photon_weight
      do n = 1, npart
        x_ = random(dseed) * (back_region % x_max - back_region % x_min) + back_region % x_min
        y_ = random(dseed) * (back_region % y_max - back_region % y_min) + back_region % y_min
        z_ = 0.5
        U_rand = random(dseed) * 2.0 - 1.0
        TH_rand = 2.0 * M_PI * random(dseed)
        u_ = photon_energy * sqrt(1.0 - U_rand**2) * cos(TH_rand)
        v_ = photon_energy * sqrt(1.0 - U_rand**2) * sin(TH_rand)
        w_ = photon_energy * U_rand
        call injectParticleGlobally(3, x_, y_, z_, u_, v_, w_, weight=photon_weight)
      end do
    end if
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
    integer :: s, ti, tj, tk, p, n, npart
    integer :: xi_glob, yi_glob
    real :: u_, v_, w_, x_, y_, z_, U_rand, TH_rand, vol
    type(region) :: back_region

    back_region % x_min = 0.0
    back_region % x_max = REAL(global_mesh % sx)
#if defined(twoD) || defined (threeD)
    back_region % y_min = 0.0
    back_region % y_max = REAL(global_mesh % sy)
#endif
#if defined(threeD)
    back_region % z_min = 0.0
    back_region % z_max = REAL(global_mesh % sz)
#endif

    if (photon_outflow) then
      ! inject photons
      vol = (back_region % x_max - back_region % x_min) * (back_region % y_max - back_region % y_min)
      ! normalize by mfp
      npart = vol * ppc0 * photon_to_lepton / photon_mfp
      do n = 1, npart
        x_ = random(dseed) * (back_region % x_max - back_region % x_min) + back_region % x_min
        y_ = random(dseed) * (back_region % y_max - back_region % y_min) + back_region % y_min
        z_ = 0.5
        U_rand = random(dseed) * 2.0 - 1.0
        TH_rand = 2.0 * M_PI * random(dseed)
        u_ = photon_energy * sqrt(1.0 - U_rand**2) * cos(TH_rand)
        v_ = photon_energy * sqrt(1.0 - U_rand**2) * sin(TH_rand)
        w_ = photon_energy * U_rand
        ! payloads are set to 0
        call injectParticleGlobally(3, x_, y_, z_, u_, v_, w_, weight=1.0)
      end do
      ! clear photons that existed for too long
      s = 3
      do ti = 1, species(s) % tile_nx
        do tj = 1, species(s) % tile_ny
          do tk = 1, species(s) % tile_nz
            do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
              species(s) % prtl_tile(ti, tj, tk) % payload1(p) = species(s) % prtl_tile(ti, tj, tk) % payload1(p) + 1.0
              if (species(s) % prtl_tile(ti, tj, tk) % payload1(p) .gt. photon_mfp) then
                species(s) % prtl_tile(ti, tj, tk) % proc(p) = -1
              end if
            end do
          end do
        end do
      end do

    end if

    ! if (photon_outflow) then
    !   s = 3
    !   do ti = 1, species(s)%tile_nx
    !     do tj = 1, species(s)%tile_ny
    !       do tk = 1, species(s)%tile_nz
    !         do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
    !           xi_glob = species(s)%prtl_tile(ti, tj, tk)%xi(p) + this_meshblock%ptr%x0
    !           yi_glob = species(s)%prtl_tile(ti, tj, tk)%yi(p) + this_meshblock%ptr%y0
    !           if ((xi_glob .eq. 0) .or. (xi_glob .eq. global_mesh%sx - 1) .or.&
    !             & (yi_glob .eq. 0) .or. (yi_glob .eq. global_mesh%sy - 1)) then
    !             species(s)%prtl_tile(ti, tj, tk)%proc(p) = -1
    !           end if
    !         end do
    !       end do
    !     end do
    !   end do
    ! end if
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
    logical :: updateE_, updateB_

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

#include "optional.F"
end module m_userfile
