module m_loadbalancing
  use m_globalnamespace
  use m_readinput, only: getInput
  use m_aux
  use m_errors
  use m_domain
  implicit none

contains
  subroutine initializeLB()
    implicit none
    ! initializing static LB variables
    slb_x = .false.; slb_sxmin = -1
    slb_y = .false.; slb_symin = -1
    slb_z = .false.; slb_szmin = -1

    alb_x = .false.; alb_sxmin = -1; alb_int_x = -1; alb_start_x = -1
    alb_y = .false.; alb_symin = -1; alb_int_y = -1; alb_start_y = -1
    alb_z = .false.; alb_szmin = -1; alb_int_z = -1; alb_start_z = -1

#if defined(oneD) || defined (twoD) || defined (threeD)
#ifdef SLB
    call getInput('static_load_balancing', 'in_x', slb_x, .false.)
    call getInput('static_load_balancing', 'sx_min', slb_sxmin, 10)
#endif

#ifdef ALB
    call getInput('adaptive_load_balancing', 'in_x', alb_x, .false.)
    call getInput('adaptive_load_balancing', 'sx_min', alb_sxmin, 10)
    call getInput('adaptive_load_balancing', 'interval_x', alb_int_x, 1000)
    call getInput('adaptive_load_balancing', 'start_x', alb_start_x, 0)
    call getInput('adaptive_load_balancing', 'slab_x', alb_slab_x, 1)
#endif
#endif
#if defined(twoD) || defined (threeD)
#ifdef SLB
    call getInput('static_load_balancing', 'in_y', slb_y, .false.)
    call getInput('static_load_balancing', 'sy_min', slb_symin, 10)
#endif

#ifdef ALB
    call getInput('adaptive_load_balancing', 'in_y', alb_y, .false.)
    call getInput('adaptive_load_balancing', 'sy_min', alb_symin, 10)
    call getInput('adaptive_load_balancing', 'interval_y', alb_int_y, 1000)
    call getInput('adaptive_load_balancing', 'start_y', alb_start_y, 1)
    call getInput('adaptive_load_balancing', 'slab_y', alb_slab_y, 1)
#endif
#endif
#if defined(threeD)
#ifdef SLB
    call getInput('static_load_balancing', 'in_z', slb_z, .false.)
    call getInput('static_load_balancing', 'sz_min', slb_szmin, 10)
#endif

#ifdef ALB
    call getInput('adaptive_load_balancing', 'in_z', alb_z, .false.)
    call getInput('adaptive_load_balancing', 'sz_min', alb_szmin, 10)
    call getInput('adaptive_load_balancing', 'interval_z', alb_int_z, 1000)
    call getInput('adaptive_load_balancing', 'start_z', alb_start_z, 1)
    call getInput('adaptive_load_balancing', 'slab_z', alb_slab_z, 1)
#endif
#endif

    call printDiag("initializeLB()", 1)
  end subroutine initializeLB

end module m_loadbalancing
