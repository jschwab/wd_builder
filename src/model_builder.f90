module model_builder

  use const_def
  use crlibm_lib
  use star_def
  use star_lib

  implicit none

  private
  public :: build_wd

  logical, parameter :: dbg = .true.

  integer, parameter :: i_Tc = 1
  integer, parameter :: i_file = 1
  integer, parameter :: i_qxq = 1

  ! data for composition from file
  integer :: num_pts
  real(dp), allocatable :: xq_data(:), xa_data(:,:)

  abstract interface
     subroutine get_xa_interface(s, q, xa)
       use const_def, only: dp
       use star_def, only: star_info
       type (star_info), pointer :: s
       real(dp), intent(in) :: q
       real(dp) :: xa(:)
     end subroutine get_xa_interface
  end interface

  procedure (get_xa_interface), pointer :: get_xa
  
contains

  
  subroutine build_wd(id, ierr)

    use num_lib, only: look_for_brackets, safe_root

    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s

    integer, parameter :: wd_lipar = 1
    integer, pointer :: ipar(:)
    integer :: wd_lrpar
    real(dp), pointer :: rpar(:)

    real(dp), pointer :: xh(:,:), q(:), dq(:)

    integer :: i, j, k, nz

    real(dp) :: mstar, mstar1, rstar, T_c, rho_c, L_core

    real(dp) :: lnd, dlnd, lnd1, lnd3, y1, y3, epsx, epsy

    integer, parameter :: imax = 100

    include 'formats'

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    if (s% x_logical_ctrl(i_file)) then
       call read_composition_from_file(s, ierr)
       get_xa => get_xa_from_file
    else
       get_xa => get_xa_from_code
    end if
    
    mstar = s% initial_mass * Msun
    s% mstar = mstar
    s% star_mass = mstar/Msun
    s% xmstar = mstar

    s% M_center = 0
    s% L_center = 0
    s% R_center = 0
    s% v_center = 0

    ! fix central temperature
    T_c = s% x_ctrl(i_Tc)

    ! rough guess for initial central density
    rho_c = 1e7 * pow_cr(s% initial_mass / 0.8d0, 5d0)

    ! pick a luminosity for the core; will impose L(m) = Lcore * m
    ! this won't be the final L, MESA is just happier if gradT = 0
    L_core = 0.1 * Lsun

    if (dbg) then
       write(*,1) 'T_c', T_c
       write(*,1) 'rho_c', rho_c
       write(*,1) 'L_core', L_core
       write(*,1) 'mstar/Msun', mstar/Msun
    end if

    wd_lrpar = 3
    allocate(rpar(wd_lrpar))
    i = 1 ! rpar(1) for mstar result
    rpar(i+1) = T_c; i = i+1
    rpar(i+1) = L_core; i = i+1

    if (i /= wd_lrpar) then
       write(*,*) 'i /= wd_lrpar', i, wd_lrpar
       write(*,*) 'wd'
       ierr = -1
       return
    end if

    allocate(ipar(wd_lipar))
    ipar(1) = id

    lnd = log_cr(rho_c)
    dlnd = 0.1d0

    call look_for_brackets(lnd, dlnd, lnd1, lnd3, wd_f, y1, y3, &
         imax, wd_lrpar, rpar, wd_lipar, ipar, ierr)
    if (ierr /= 0) then
       if (dbg) then
          if (dbg) write(*,*) 'look_for_brackets ierr', ierr
          write(*,1) 'lnd1', lnd1
          write(*,1) 'lnd3', lnd3
          write(*,1) 'y1', y1
          write(*,1) 'y3', y3
       end if
       return
    end if


    epsx = 1d-3 ! limit for variation in lnd
    epsy = 1d-3 ! limit for matching desired mass as fraction of total mass

    lnd = safe_root(wd_f, lnd1, lnd3, y1, y3, imax, epsx, epsy, &
         wd_lrpar, rpar, wd_lipar, ipar, ierr)
    if (ierr /= 0) then
       if (dbg) write(*,*) 'safe_root ierr', ierr
       return
    end if

    mstar1 = rpar(1)

    xh => s% xh
    q => s% q
    dq => s% dq
    nz = s% nz

    if (dbg) then
       write(*,*)
       write(*,*) 'finished build_wd model'
       write(*,1) 'mstar1/Msun', mstar1/Msun
       write(*,1) '(mstar-mstar1)/mstar', (mstar-mstar1)/mstar
       write(*,1) 'log10(r/Rsun)', log10_cr(exp_cr(xh(s% i_lnR,1))/Rsun)
       if (s% i_lum /= 0) write(*,1) 'log10(L/Lsun)', log10_cr(xh(s% i_lum,1)/Lsun)
       write(*,1) 'log10(Tsurf)', xh(s% i_lnT,1)/ln10
       write(*,1) 'Tsurf', exp_cr(xh(s% i_lnT,1))
       write(*,*) 'nz', nz
       write(*,*)
    end if

    ! The following deallocations deal with arrays which were
    ! needed in the root bracket/solve above, but will get
    ! overwritten during the call to allocate_star_info_arrays

    if (ASSOCIATED(s% xh_old)) deallocate(s% xh_old)
    if (ASSOCIATED(s% xh_older)) deallocate(s% xh_older)
    if (ASSOCIATED(s% equ1)) deallocate(s% equ1)
    if (ASSOCIATED(s% xh_pre)) deallocate(s% xh_pre)

    call star_allocate_arrays(id, ierr)
    if (ierr /= 0) then
       return
    end if

    do k=1,nz
       do j=1,s% nvar_hydro
          s% xh(j,k) = xh(j,k)
       end do
       s% q(k) = q(k)
       s% dq(k) = dq(k)
       call get_xa(s, s% q(k), s% xa(:, k))
    end do

    deallocate(xh, q, dq, ipar, rpar)

  end subroutine build_wd

  real(dp) function wd_f(lnd, dfdx, lrpar, rpar, lipar, ipar, ierr)
    integer, intent(in) :: lrpar, lipar
    real(dp), intent(in) :: lnd
    real(dp), intent(out) :: dfdx
    integer, intent(inout), pointer :: ipar(:) ! (lipar)
    real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
    integer, intent(out) :: ierr

    type (star_info), pointer :: s
    real(dp) :: rho_c, T_c, L_core, d_log10_P
    real(dp), pointer :: xa(:)
    integer :: i, nz
    real(dp) :: mstar, mstar1

    logical, parameter :: dbg = .true.

    include 'formats'

    ierr = 0
    wd_f = 0
    if (lipar <= 0) then
       write(*,*) 'lipar', lipar
       write(*,*) 'wd f'
       ierr = -1
       return
    end if

    call star_ptr(ipar(1), s, ierr)
    if (ierr /= 0) return

    if (associated(s% xh)) deallocate(s% xh)
    if (associated(s% q)) deallocate(s% q)
    if (associated(s% dq)) deallocate(s% dq)

    rho_c = exp_cr(lnd)

    i = 1 ! rpar(1) for mstar result
    T_c = rpar(i+1); i = i+1
    L_core = rpar(i+1); i = i+1
    if (i > lrpar) then
       write(*,*) 'i > lrpar', i, lrpar
       write(*,*) 'wd f'
       ierr = -1
       return
    end if

    mstar = s% mstar ! desired value
    mstar1 = mstar ! to keep gfortran quiet

    if (dbg) write(*,*) 'call build1_wd_model'
    call build1_wd_model(s, T_c, rho_c, L_core, nz, mstar1, ierr)
    if (ierr /= 0) then
       write(*,*) 'failed in build1_wd_model'
       return
    end if

    s% nz = nz

    rpar(1) = mstar1 ! return the actual mass

    wd_f = (mstar - mstar1) / mstar
    dfdx = 0

    if (dbg) then
       write(*,1) 'rho_c', rho_c
       write(*,1) 'mstar1', mstar1 / msun
       write(*,1) 'wd_f', wd_f
       write(*,*)
    end if

  end function wd_f


  subroutine build1_wd_model(s, T_c, rho_c, L_core, nz, mstar, ierr)
    use chem_def
    use eos_def
    use kap_lib
    use chem_lib
    use eos_lib, only: Radiation_Pressure
    use eos_support, only: solve_eos_given_PgasT
    type (star_info), pointer :: s
    real(dp), intent(in) :: T_c, rho_c, L_core
    real(dp) :: x, z, abar, zbar
    real(dp), allocatable :: xa(:)
    integer, intent(out) :: nz
    real(dp), intent(out) :: mstar ! the mass of the constructed model
    integer, intent(out) :: ierr

    real(dp) :: logRho_guess
    real(dp), parameter :: LOGRHO_TOL = 1E-6_dp
    real(dp), parameter :: LOGPGAS_TOL = 1E-6_dp

    integer :: i, ii, k, j, i_lnd, i_lnT, i_lnR, prune, max_retries
    real(dp), parameter :: &
         dlogPgas = 0.01d0, q_at_nz = 1d-5
    real(dp) :: &
         P_surf_limit, y, logPgas, Prad, Pgas, try_dlogPgas, logPgas0, &
         res(num_eos_basic_results), P_c, logP, m, &
         d_eos_dlnd(num_eos_basic_results), d_eos_dlnT(num_eos_basic_results), &
         d_eos_dabar(num_eos_basic_results), d_eos_dzbar(num_eos_basic_results), &
         lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
         cgrav, r, rmid, rho, logRho, T, lnT, L, P, P0, dm, m0, L0, r0, lnT0, T0, &
         rho0, rho_mid, Pmid, chiRho0, chiRho_mid, chiT0, chiT_mid, Cp0, Cp_mid, &
         grada0, grada_mid, mmid, Tmid, Lmid, &
         chiRho, chiT, Cp, grada, gradT, logT_surf_limit, logP_surf_limit
    real(dp), pointer :: xh(:,:), q(:), dq(:) ! model structure info

    real(dp) :: opacity, kap_cond, kap_rad, dlnkap_dlnd, dlnkap_dlnT, frac_Type2
    real(dp) :: Zbase, XC, XN, XO, XNe

    real(dp) :: z2bar, ye, mass_correction, sumx
    integer :: species

    logical, parameter :: dbg = .false.

    include 'formats'

    ierr = 0

    logP_surf_limit = s% job% pre_ms_logP_surf_limit
    P_surf_limit = exp10_cr(logP_surf_limit)
    if (dbg) write(*,1) 'logP_surf_limit', logP_surf_limit

    
    logT_surf_limit = s% job% pre_ms_logT_surf_limit
    if (logT_surf_limit <= 0) logT_surf_limit = 3.7d0
    if (dbg) write(*,1) 'logT_surf_limit', logT_surf_limit

    i_lnd = s% i_lnd
    i_lnT = s% i_lnT
    i_lnR = s% i_lnR

    if (i_lnd == 0) then
       write(*,*) 'Sorry: require lnPgas_flag be .false. for build1_wd_model'
       ierr = -1
       return
    end if

    if (i_lnT == 0) then
       write(*,*) 'Sorry: require E_flag be .false. for build1_wd_model'
       ierr = -1
       return
    end if

    cgrav = standard_cgrav

    allocate(xa(s% species))

    call get_xa(s, 0d0, xa)

    call set_composition_info

    call star_get_eos( &
         s, 0, z, x, abar, zbar, xa, &
         rho_c, log10_cr(rho_c), T_c, log10_cr(T_c), &
         res, d_eos_dlnd, d_eos_dlnT, &
         d_eos_dabar, d_eos_dzbar, ierr)
    if (ierr /= 0) then
       write(*,*) 'failed in get_eos'
       return
    end if
    call unpack_eos_results

    logPgas = res(i_lnPgas)/ln10
    Pgas = exp10_cr(logPgas)
    P_c = Pgas + Radiation_Pressure(T_c) ! center pressure

    mstar = s% mstar ! desired total mass
    m = q_at_nz*mstar ! mass at nz
    ! pressure at innermost point using K&W 10.6
    P = P_c - 3*cgrav/(8*pi)*pow_cr(pi4*rho_c/3,4d0/3d0)*pow_cr(m,two_thirds)
    logP = log10_cr(P)

    ! estimate nz from lgP
    nz = 1 + (logP - logP_surf_limit)/dlogPgas

    ! temperature at nz assuming isothermal
    lnT = log_cr(T_c)
    T = exp_cr(lnT)

    ! density at nz
    logRho_guess = log10_cr(rho_c)
    call solve_eos_given_PgasT( &
         s, 0, z, xa(s% net_iso(ih1)), abar, zbar, xa, &
         lnT/ln10, log10_cr(Pgas), logRho_guess, LOGRHO_TOL, LOGPGAS_TOL, &
         logRho, res, d_eos_dlnd, d_eos_dlnT, d_eos_dabar, d_eos_dzbar, &
         ierr)
    if (ierr /= 0) return
    rho = exp10_cr(logRho)
    call unpack_eos_results

    r = pow_cr(m/(pi4*rho/3),one_third) ! radius at nz

    y = 1 - (x+z)

    L = L_core * (m/mstar)


    call kap_get( &
         s% kap_handle, zbar, X, Z, Zbase, XC, XN, XO, XNe, log10_cr(rho), log10_Cr(T), &
         lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
         frac_Type2, opacity, dlnkap_dlnd, dlnkap_dlnT, ierr)

    if (ierr /= 0) then
       write(*,*) 'failed in kap_get'
       return
    end if

    call eval_gradT( &
         s, zbar, x, y, xa, rho, m, mstar, r, T, lnT, L, P, &
         chiRho, chiT, Cp, opacity, grada, &
         lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
         gradT, ierr )
    if (ierr /= 0) return



    allocate(xh(s% nvar_hydro,nz), q(nz), dq(nz), stat=ierr)
    if (ierr /= 0) return
    s% xh => xh
    s% dq => dq
    s% q => q

    xh(i_lnd, nz) = logRho*ln10
    xh(i_lnT, nz) = lnT
    xh(i_lnR, nz) = log_cr(r)
    if (s% i_lum /= 0) xh(s% i_lum,nz) = L

    q(nz) = q_at_nz
    dq(nz) = q_at_nz

    if (dbg) write(*,*) 'nz', nz

    max_retries = 10
    prune = 0
    step_loop: do k = nz-1, 1, -1

       try_dlogPgas = dlogPgas
       logPgas0 = logPgas
       P0 = P
       m0 = m
       L0 = L
       r0 = r
       lnT0 = lnT
       T0 = T
       rho0 = rho
       chiRho0 = chiRho
       chiT0 = chiT
       Cp0 = Cp
       grada0 = grada
       dm = 0 ! for gfortran

       if (dbg) write(*,3) 'step', k, nz, logPgas0

       retry_loop: do j = 1, max_retries

          logPgas = logPgas0 - try_dlogPgas
          Pgas = exp10_cr(logPgas)

          if (j > 1) write(*,2) 'retry', j, logPgas

          do i = 1, 2

             Prad = Radiation_Pressure(T)
             P = Pgas + Prad

             rho_mid = (rho+rho0)/2

             do ii = 1, 10 ! repeat to get hydrostatic balance
                rmid = pow_cr((r*r*r + r0*r0*r0)/2,one_third)
                mmid = (m + m0)/2
                if (ii == 10) exit
                dm = -pi4*pow4(rmid)*(P-P0)/(cgrav*mmid)
                m = m0 + dm ! mass at point k
                r = pow_cr(r0*r0*r0 + dm/((4*pi/3)*rho_mid),one_third)
                if (dbg) write(*,2) 'r', ii, r, m, dm
             end do

             L = L0 + dm*(L_core/mstar) ! luminosity at point k
             Lmid = (L0+L)/2

             Pmid = (P+P0)/2

             chiRho_mid = (chiRho0 + chiRho)/2
             chiT_mid = (chiT0 + chiT)/2
             Cp_mid = (Cp0 + Cp)/2
             grada_mid = (grada0 + grada)/2

             do ii = 1, 2
                Tmid = (T+T0)/2

                call get_xa(s, m/mstar, xa)

                call set_composition_info

                call kap_get( &
                     s% kap_handle, zbar, X, Z, Zbase, XC, XN, XO, XNe, log10_cr(rho_mid), log10_Cr(Tmid), &
                     lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                     frac_Type2, opacity, dlnkap_dlnd, dlnkap_dlnT, ierr)
                if (ierr /= 0) then
                   write(*,*) 'failed in kap_get'
                   return
                end if

                call kap_get_elect_cond_opacity( &
                     zbar, log10_cr(rho_mid), log10_Cr(Tmid), &
                     kap_cond, dlnkap_dlnd, dlnkap_dlnT, ierr)

                kap_rad = 1d0/(1d0/opacity - 1d0/kap_cond)

                if (kap_cond .lt. kap_rad) then ! suggested by Evan
                   call eval_gradT( &
                        s, zbar, x, y, xa, rho_mid, mmid, mstar, rmid, Tmid, log_cr(Tmid), Lmid, Pmid, &
                        chiRho_mid, chiT_mid, Cp_mid, opacity, grada_mid, &
                        lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                        gradT, ierr )
                   if (ierr /= 0) return
                else
                   gradT = 0.235 ! suggested by Lars
                end if

                T = T0 + Tmid*gradT*(P-P0)/Pmid
                lnT = log_cr(T)
                if (dbg) write(*,2) 'T', ii, T
             end do

             if (i == 2) exit

             logRho_guess = logRho
             call solve_eos_given_PgasT( &
                  s, 0, z, x, abar, zbar, xa, &
                  lnT/ln10, logPgas, logRho_guess, LOGRHO_TOL, LOGPGAS_TOL, &
                  logRho, res, d_eos_dlnd, d_eos_dlnT, d_eos_dabar, d_eos_dzbar, &
                  ierr)
             rho = exp10_cr(logRho)
             if (ierr /= 0) return
             call unpack_eos_results

          end do

          if (lnT <= logT_surf_limit*ln10) then
             if (dbg) write(*,*) 'have reached lgT_surf_limit', lnT/ln10, logT_surf_limit
             prune = k
             exit step_loop
          end if

          if (P <= P_surf_limit) then
             if (dbg) write(*,1) 'have reached P_surf limit', P, P_surf_limit
             prune = k
             exit step_loop
          end if

          if (logRho <= -12) then
             if (dbg) write(*,1) 'have reached logRho_surf limit', logRho, -12d0
             prune = k
             exit step_loop
          end if


          xh(i_lnd, k) = logRho*ln10
          xh(i_lnT, k) = lnT
          xh(i_lnR, k) = log_cr(r)
          if (s% i_lum /= 0) xh(s% i_lum,k) = L
          q(k) = m/mstar
          dq(k) = dm/mstar

          if (dbg) then
             write(*,2) 'xh(i_lnd, k)', k, xh(i_lnd, k)
             write(*,2) 'xh(i_lnT, k)', k, xh(i_lnT, k)
             write(*,2) 'xh(i_lnR, k)', k, xh(i_lnR, k)
             write(*,2) 'L', k, L
             write(*,2) 'q(k)', k, q(k)
             write(*,2) 'dq(k)', k, dq(k)
          end if

          exit retry_loop

       end do retry_loop

    end do step_loop

    if (prune > 0) then ! move stuff and reduce nz
       if (dbg) write(*,*) 'prune', prune
       do k=1,nz-prune
          xh(:,k) = xh(:,k+prune)
          q(k) = q(k+prune)
          dq(k) = dq(k+prune)
       end do
       m = mstar*q(1)
       nz = nz-prune
       if (dbg) write(*,*) 'final nz', nz
    end if

    mstar = m ! actual total mass

    call star_normalize_dqs(nz, dq, ierr)
    if (ierr /= 0) then
       if (s% report_ierr) write(*,*) 'set_qs failed in pre ms model'
       return
    end if
    call star_set_qs(nz, q, dq, ierr)
    if (ierr /= 0) then
       if (s% report_ierr) write(*,*) 'set_qs failed in pre ms model'
       return
    end if

  contains

    subroutine set_composition_info

      call basic_composition_info( &
           s% species, s% chem_id, xa, x, y, z, abar, zbar, z2bar, ye, &
           mass_correction, sumx)

      xc = xa(s% net_iso(ic12))
      xn = xa(s% net_iso(in14))
      xo = xa(s% net_iso(io16))
      xne = xa(s% net_iso(ine20))
      Zbase = s% Zbase

    end subroutine set_composition_info


    subroutine unpack_eos_results
      chiRho = res(i_chiRho)
      chiT = res(i_chiT)
      Cp = res(i_cp)
      grada = res(i_grad_ad)
      lnfree_e = res(i_lnfree_e)
      d_lnfree_e_dlnRho = d_eos_dlnd(i_lnfree_e)
      d_lnfree_e_dlnT = d_eos_dlnT(i_lnfree_e)
    end subroutine unpack_eos_results


  end subroutine build1_wd_model


  subroutine eval_gradT( &
       s, zbar, x, y, xa, rho, m, mstar, r, T, lnT, L, P, &
       chiRho, chiT, Cp, opacity, grada, &
       lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
       gradT, ierr )
    use chem_def, only: ih1
    use kap_lib, only : kap_get

    type (star_info), pointer :: s
    real(dp), intent(in) :: &
         zbar, x, y, xa(:), rho, m, mstar, r, T, lnT, L, P, &
         chiRho, chiT, Cp, opacity, grada, &
         lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
    real(dp), intent(out) :: gradT
    integer, intent(out) :: ierr

    real(dp) :: Z, Zbase, XC, XN, XO, XNe, frac_Type2

    real(dp) :: &
         cgrav, Cv, csound, &
         prev_conv_vel, max_conv_vel, dt, gradL_composition_term, tau
    integer :: mixing_type
    real(dp) :: mlt_basics(num_mlt_results)
    real(dp), target :: mlt_partials1_ary(num_mlt_partials*num_mlt_results)
    real(dp), pointer :: mlt_partials1(:), mlt_partials(:,:)
    real(dp), parameter :: alpha_semiconvection = 0, thermohaline_coeff = 0, &
         gradr_factor = 1
    real(dp) :: gradT_smooth_low, gradT_smooth_high, alfa, beta, &
         d_alfa_dq00, d_alfa_dqm1, d_alfa_dqp1, &
         normal_mlt_gradT_factor, &
         T_00, T_m1, rho_00, rho_m1, P_00, P_m1, &
         chiRho_00, d_chiRho_00_dlnd, d_chiRho_00_dlnT, &
         chiRho_m1, d_chiRho_m1_dlnd, d_chiRho_m1_dlnT, &
         chiT_00, d_chiT_00_dlnd, d_chiT_00_dlnT, &
         chiT_m1, d_chiT_m1_dlnd, d_chiT_m1_dlnT, &
         Cp_00, d_Cp_00_dlnd, d_Cp_00_dlnT, &
         Cp_m1, d_Cp_m1_dlnd, d_Cp_m1_dlnT, &
         opacity_00, d_opacity_00_dlnd, d_opacity_00_dlnT, &
         opacity_m1, d_opacity_m1_dlnd, d_opacity_m1_dlnT, &
         grada_00, d_grada_00_dlnd, d_grada_00_dlnT, &
         grada_m1, d_grada_m1_dlnd, d_grada_m1_dlnT
    logical :: smooth_gradT
    smooth_gradT = .false.
    normal_mlt_gradT_factor = 1d0

    ierr = 0
    mlt_partials1 => mlt_partials1_ary
    mlt_partials(1:num_mlt_partials,1:num_mlt_results) => &
         mlt_partials1(1:num_mlt_partials*num_mlt_results)

    cgrav = standard_cgrav
    gradL_composition_term = 0
    Cv = Cp
    tau = 1
    prev_conv_vel = -1
    max_conv_vel = 1d99
    dt = -1
    csound = 0 ! not used when dt <= 0
    ! not used
    alfa=0d0; beta=0d0; d_alfa_dq00=0d0; d_alfa_dqm1=0d0; d_alfa_dqp1=0d0
    T_00=0d0; T_m1=0d0; rho_00=0d0; rho_m1=0d0; P_00=0d0; P_m1=0d0
    chiRho_00=0d0; d_chiRho_00_dlnd=0d0; d_chiRho_00_dlnT=0d0
    chiRho_m1=0d0; d_chiRho_m1_dlnd=0d0; d_chiRho_m1_dlnT=0d0
    chiT_00=0d0; d_chiT_00_dlnd=0d0; d_chiT_00_dlnT=0d0
    chiT_m1=0d0; d_chiT_m1_dlnd=0d0; d_chiT_m1_dlnT=0d0
    Cp_00=0d0; d_Cp_00_dlnd=0d0; d_Cp_00_dlnT=0d0
    Cp_m1=0d0; d_Cp_m1_dlnd=0d0; d_Cp_m1_dlnT=0d0
    opacity_00=0d0; d_opacity_00_dlnd=0d0; d_opacity_00_dlnT=0d0
    opacity_m1=0d0; d_opacity_m1_dlnd=0d0; d_opacity_m1_dlnT=0d0
    grada_00=0d0; d_grada_00_dlnd=0d0; d_grada_00_dlnT=0d0
    grada_m1=0d0; d_grada_m1_dlnd=0d0; d_grada_m1_dlnT=0d0

    call star_mlt_eval( &
         s% id, 0, cgrav, m, mstar, r, L, x, T, rho, P, &
         chiRho, chiT, Cp, opacity, grada, &

                                ! not used
         alfa, beta, d_alfa_dq00, d_alfa_dqm1, d_alfa_dqp1, &
         T_00, T_m1, rho_00, rho_m1, P_00, P_m1, &
         chiRho_00, d_chiRho_00_dlnd, d_chiRho_00_dlnT, &
         chiRho_m1, d_chiRho_m1_dlnd, d_chiRho_m1_dlnT, &
         chiT_00, d_chiT_00_dlnd, d_chiT_00_dlnT, &
         chiT_m1, d_chiT_m1_dlnd, d_chiT_m1_dlnT, &
         Cp_00, d_Cp_00_dlnd, d_Cp_00_dlnT, &
         Cp_m1, d_Cp_m1_dlnd, d_Cp_m1_dlnT, &
         opacity_00, d_opacity_00_dlnd, d_opacity_00_dlnT, &
         opacity_m1, d_opacity_m1_dlnd, d_opacity_m1_dlnT, &
         grada_00, d_grada_00_dlnd, d_grada_00_dlnT, &
         grada_m1, d_grada_m1_dlnd, d_grada_m1_dlnT, &

         gradr_factor, gradL_composition_term, &
         alpha_semiconvection, s% semiconvection_option, &
         thermohaline_coeff, s% thermohaline_option, ih1, &
         s% mixing_length_alpha, s% alt_scale_height_flag, s% remove_small_D_limit, &
         s% MLT_option, s% Henyey_MLT_y_param, s% Henyey_MLT_nu_param, &
         gradT_smooth_low, gradT_smooth_high, smooth_gradT, &
         normal_mlt_gradT_factor, &
         prev_conv_vel, max_conv_vel, s% mlt_accel_g_theta, dt, tau, .false., &
         mixing_type, mlt_basics, mlt_partials1, ierr)
    if (ierr /= 0) return

    gradT = mlt_basics(mlt_gradT) ! actual temperature gradient dlnT/dlnP

  end subroutine eval_gradT

  subroutine get_xa_from_code(s, q, xa)

    use chem_def

    type (star_info), pointer :: s
    real(dp), intent(in) :: q
    real(dp) :: xa(:)

    xa = 0

    if (q < 0.50d0) then
       xa(s% net_iso(ihe4))  = 0.0d0
       xa(s% net_iso(ic12))  = 0.23d0
       xa(s% net_iso(in14))  = 0.0d0
       xa(s% net_iso(io16))  = 0.75d0
       xa(s% net_iso(ine22)) = 0.02d0

    else if (0.50 <= q .and. q < 0.90) then
       xa(s% net_iso(ihe4))  = 0.0d0
       xa(s% net_iso(ic12))  = 0.48d0
       xa(s% net_iso(in14))  = 0.0d0
       xa(s% net_iso(io16))  = 0.50d0
       xa(s% net_iso(ine22)) = 0.02d0

    else if (0.90 <= q .and. q < 0.999) then
       xa(s% net_iso(ihe4))  = 0.0d0
       xa(s% net_iso(ic12))  = 0.65d0
       xa(s% net_iso(in14))  = 0.02d0
       xa(s% net_iso(io16))  = 0.33d0
       xa(s% net_iso(ine22)) = 0.00d0

    else if (0.999 < q) then
       xa(s% net_iso(ihe4))  = 0.98d0
       xa(s% net_iso(ic12))  = 0.0d0
       xa(s% net_iso(in14))  = 0.02d0
       xa(s% net_iso(io16))  = 0.0d0
       xa(s% net_iso(ine22)) = 0.0d0
    end if


  end subroutine get_xa_from_code


  subroutine read_composition_from_file(s, ierr)

    type (star_info), pointer :: s
    
    integer :: num_species
    integer :: i, iounit, ierr
    
    open(newunit=iounit, file=trim(s% job% relax_composition_filename), &
         status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
       write(*,*) 'open failed', ierr, iounit
       write(*, '(a)') 'failed to open ' // trim(s% job% relax_composition_filename)
       close(iounit)
       return
    end if
    
    read(iounit, *, iostat=ierr) num_pts, num_species
    if (ierr /= 0) then
       close(iounit)
       write(*, '(a)') 'failed while trying to read 1st line of ' // &
            trim(s% job% relax_composition_filename)
       return
    end if
    
    if(num_species .ne. s% species) then
       write(*,*) 'Error in ',trim(s% job% relax_composition_filename)
       write(*,'(a,I4,a)') 'got ',num_species,' species'
       write(*,'(a,I4,a)') 'expected ', s% species,' species'
       write(*,*)
       ierr=-1
       return
    end if
    
    allocate(xq_data(num_pts), xa_data(num_species,num_pts))
    do i = 1, num_pts
       read(iounit,*,iostat=ierr) xq_data(i), xa_data(1:num_species,i)
       if (ierr /= 0) then
          close(iounit)
          write(*, '(a)') &
               'failed while trying to read ' // trim(s% job% relax_composition_filename)
          write(*,*) 'line', i+1
          write(*,*) 'perhaps wrong info in 1st line?'
          write(*,*) '1st line must have num_pts and num_species in that order'
          deallocate(xq_data, xa_data)
          return
       end if
    end do
    close(iounit)

  end subroutine read_composition_from_file
  


  subroutine get_xa_from_file(s, q, xa)
  
    use interp_1d_def, only: pm_work_size
    use interp_1d_lib, only: interpolate_vector, interp_pm
  
    type (star_info), pointer :: s
    real(dp), intent(in) :: q
    real(dp) :: xa(:)
  
    real(dp), pointer :: work(:)
  
    real(dp) :: x_new(1), v_new(1)
  
    integer :: j, ierr
  
    allocate(work(pm_work_size*num_pts))

    if (s% x_logical_ctrl(i_qxq)) then
       x_new(1) = q
    else
       x_new(1) = 1d0 - q
    end if
    
    do j = 1, s% species
       call interpolate_vector( &
            num_pts, xq_data, 1, x_new, xa_data(j,:), v_new, &
            interp_pm, pm_work_size, work, 'get_xa_from_target', ierr)
       xa(j) = v_new(1)
    end do
  
    deallocate(work)
  
  end subroutine get_xa_from_file
  


end module model_builder
