module model_builder

  use const_def
  use crlibm_lib
  use star_def
  use star_lib

  implicit none

  private
  public :: build_wd

  integer, parameter :: i_xo16 = 1
  
  logical, parameter :: dbg = .true.

contains

  subroutine build_wd(id, ierr)

    use chem_def
    use chem_lib, only : basic_composition_info
    use num_lib, only: look_for_brackets, safe_root
    
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s

    real(dp) :: x, y, z, abar, zbar, z2bar, ye, mass_correction, sumx
    integer :: species
    real(dp), dimension(:), allocatable :: xa

    integer, parameter :: pre_ms_lipar = 1
    integer, pointer :: ipar(:)
    integer :: pre_ms_lrpar
    real(dp), pointer :: rpar(:)

    real(dp), pointer :: xh(:,:), q(:), dq(:)

    integer :: i, j, k, nz
     
  !    real(dp) :: mu_eff ! effective mean molecular weight for gas particles
  !    real(dp) :: initial_y, initial_h1, initial_h2, initial_he3, initial_he4, &
  !       xsol_he3, xsol_he4
  !    integer :: initial_zfracs
  !    real(dp), parameter :: max_mass_to_create = 30, min_mass_to_create = 0.03d0

  !    include 'formats'

  !    ipar => ipar_ary
  !    rpar => rpar_ary

  !    pre_ms_lrpar = rpar_init+nmet+species
    
    
    !       abar, zbar, mass_correction, z2bar, ye, sumx, &
  !       lgL, eps_grav, metals(nmet), lnd, dlnd, lnd1, lnd3, y1, y3, epsx, epsy, &
  !       guess_rho_c, d_log10_P

   ! real(mstar, mstar1, lgM,  &

        ! lgL, eps_grav, metals(nmet), &
        ! T_c, guess_rho_c, d_log10_P

    real(dp) :: mstar, mstar1, rstar, mu_eff, T_c, rho_c, lgL, eps_grav

    real(dp) :: xc, xn, xo, xne, xheavy, Zbase

    real(dp) :: lnd, dlnd, lnd1, lnd3, y1, y3, epsx, epsy

    integer :: funit
    
    integer, parameter :: imax = 100

    include 'formats'
    
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    species = s% species

    mstar = s% initial_mass * Msun
    s% mstar = mstar
    s% star_mass = mstar/Msun
    s% xmstar = mstar

    s% M_center = 0
    s% L_center = 0
    s% R_center = 0
    s% v_center = 0

    allocate(xa(species))
    xa = 0
    xa(s% net_iso(ic12)) = 1d0 - s% x_ctrl(i_xo16)
    xa(s% net_iso(io16)) = s% x_ctrl(i_xo16)
    
    call basic_composition_info( &
        species, s% chem_id, xa, x, y, z, abar, zbar, z2bar, ye, &
        mass_correction, sumx)

    if (dbg) then
       write(*,1) 'xc', 1d0 - xa(s% net_iso(io16))
       write(*,1) 'xo', xa(s% net_iso(io16))
    end if

    ! guess for Tc and Rhoc
    ! fix Tc, iterate Rhoc
    T_c = 2e8

    ! all low mass WDs are around here
    rho_c = 1e7

    ! pick a luminosity
    lgL = -1d0

    ! use uniform eps_grav to give that luminosity
    eps_grav = exp10_cr(lgL)*Lsun/mstar

    
    if (dbg) then
       write(*,1) 'T_c', T_c
       write(*,1) 'rho_c', rho_c
       write(*,1) 'lgL', lgL
       write(*,1) 'eps_grav', eps_grav
       write(*,1) 'mstar/Msun', mstar/Msun
    end if

    pre_ms_lrpar = 7 + species 
    allocate(rpar(pre_ms_lrpar))
    i = 1 ! rpar(1) for mstar result
    rpar(i+1) = T_c; i = i+1
    rpar(i+1) = eps_grav; i = i+1
    rpar(i+1) = x; i = i+1
    rpar(i+1) = z; i = i+1
    rpar(i+1) = abar; i = i+1
    rpar(i+1) = zbar; i = i+1

    rpar(i+1:i+species) = xa(1:species); i = i+species

    if (i /= pre_ms_lrpar) then
       write(*,*) 'i /= pre_ms_lrpar', i, pre_ms_lrpar
       write(*,*) 'pre ms'
       ierr = -1
       return
    end if

    allocate(ipar(pre_ms_lipar))
    ipar(1) = id

     lnd = log_cr(rho_c)
     dlnd = 0.1d0

     call look_for_brackets(lnd, dlnd, lnd1, lnd3, pre_ms_f, y1, y3, &
           imax, pre_ms_lrpar, rpar, pre_ms_lipar, ipar, ierr)
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

     lnd = safe_root(pre_ms_f, lnd1, lnd3, y1, y3, imax, epsx, epsy, &
              pre_ms_lrpar, rpar, pre_ms_lipar, ipar, ierr)
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
        write(*,*) 'finished pre-MS model'
        write(*,1) 'mstar1/Msun', mstar1/Msun
        write(*,1) '(mstar-mstar1)/mstar', (mstar-mstar1)/mstar
        write(*,1) 'log10(r/Rsun)', log10_cr(exp_cr(xh(s% i_lnR,1))/Rsun)
        if (s% i_lum /= 0) write(*,1) 'log10(L/Lsun)', log10_cr(xh(s% i_lum,1)/Lsun)
        write(*,1) 'log10(Tsurf)', xh(s% i_lnT,1)/ln10
        write(*,1) 'Tsurf', exp_cr(xh(s% i_lnT,1))
        write(*,*) 'nz', nz
        write(*,*)
        !stop 'debug: pre ms'
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

     open(newunit=funit, file='initial_model.dat')

     do j=1,s% nvar_hydro
        write(funit, *) '! ', j, trim(s% nameofvar(j))
     end do
     
     do k=1,nz
        do j=1,s% nvar_hydro
           s% xh(j,k) = xh(j,k)
        end do
        write(funit,*) xh(:,k)
        do j=1,species
           s% xa(j,k) = xa(j)
        end do
        s% q(k) = q(k)
        s% dq(k) = dq(k)
     end do

     close(funit)

     deallocate(xa, xh, q, dq)

  end subroutine build_wd


  ! subroutine build_pre_ms_model(id, s, nvar_hydro, species, ierr)
  !    use chem_def
  !    use chem_lib, only: basic_composition_info, chem_Xsol
  !    use adjust_xyz, only: get_xa_for_standard_metals
  !    use alloc, only: allocate_star_info_arrays
  !    use num_lib, only: look_for_brackets, safe_root

  !    type (star_info), pointer :: s
  !    integer, intent(in) :: id, nvar_hydro, species
  !    integer, intent(out) :: ierr

  !    real(dp) :: &
  !       initial_z, x, y, z, xa(species), mstar, mstar1, lgM, rstar, rho_c, &
  !       abar, zbar, mass_correction, z2bar, ye, sumx, &
  !       lgL, eps_grav, metals(nmet), lnd, dlnd, lnd1, lnd3, y1, y3, epsx, epsy, &
  !       T_c, guess_rho_c, d_log10_P
  !    integer :: i, j, k, nz, pre_ms_lrpar
  !    real(dp), pointer :: xh(:,:), q(:), dq(:)
  !    integer, parameter :: pre_ms_lipar = 1, imax = 100, rpar_init = 8
  !    integer, target :: ipar_ary(pre_ms_lipar)
  !    integer, pointer :: ipar(:)
  !    real(dp), target :: rpar_ary(rpar_init+nmet+species) ! (pre_ms_lrpar)
  !    real(dp), pointer :: rpar(:)
  !    real(dp) :: mu_eff ! effective mean molecular weight for gas particles
  !    real(dp) :: initial_y, initial_h1, initial_h2, initial_he3, initial_he4, &
  !       xsol_he3, xsol_he4
  !    integer :: initial_zfracs
  !    real(dp), parameter :: max_mass_to_create = 30, min_mass_to_create = 0.03d0

  !    include 'formats'

  !    ipar => ipar_ary
  !    rpar => rpar_ary

  !    pre_ms_lrpar = rpar_init+nmet+species

  !    if (nvar_hydro > 4) then
  !       write(*,*) 'sorry, build_pre_ms_model only supports the basic 4 vars.'
  !       ierr = -1
  !       return
  !    end if
  !    mstar = min(max_mass_to_create*Msun, max(min_mass_to_create*Msun, s% mstar))
  !    s% mstar = mstar
  !    s% star_mass = mstar/Msun
  !    s% xmstar = mstar

  !    T_c = s% pre_ms_T_c
  !    if (T_c <= 0) T_c = 9d5
  !    if (T_c >= 1d6) then
  !       write(*,1) 'log T_c', log10_cr(T_c)
  !       write(*,*) 'center temperature is too high for a pre main sequence model.'
  !       write(*,*) 'please pick T_c below 10^6 so code can ignore nuclear reactions.'
  !       ierr = -1
  !       return
  !    end if

  !    ierr = 0
  !    initial_z = s% initial_z
  !    initial_y = s% initial_y
  !    if (initial_y < 0) initial_y = max(0d0, min(1d0, 0.24d0 + 2*initial_z))
  !    s% M_center = 0
  !    s% L_center = 0
  !    s% R_center = 0
  !    s% v_center = 0

  !    initial_h1 = max(0d0, min(1d0, 1d0 - (initial_z + initial_y)))
  !    initial_h2 = 0
  !    if (s% initial_he3 < 0d0) then
  !       xsol_he3 = chem_Xsol('he3')
  !       xsol_he4 = chem_Xsol('he4')
  !       initial_he3 = initial_y*xsol_he3/(xsol_he3 + xsol_he4)
  !       initial_he4 = initial_y*xsol_he4/(xsol_he3 + xsol_he4)
  !    else if (s% initial_he3 < initial_y) then
  !       initial_he3 = s% initial_he3
  !       initial_he4 = initial_y - s% initial_he3
  !    else
  !       write(*,*) "ERROR: initial_he3 is larger than initial_y"
  !       ierr = -1
  !    end if
  !    initial_zfracs = s% pre_ms_initial_zfracs

  !    call get_xa_for_standard_metals(s, &
  !       species, s% chem_id, s% net_iso, &
  !       initial_h1, initial_h2, initial_he3, initial_he4, &
  !       initial_zfracs, s% pre_ms_dump_missing_heaviest, &
  !       xa, ierr)
  !    if (ierr /= 0) return

  !    if (dbg) then
  !       write(*,*) 'abundances'
  !       do i=1,species
  !          write(*,1) trim(chem_isos% name(s% chem_id(i))), xa(i)
  !       end do
  !       write(*,*)
  !    end if

  !    if (abs(1-sum(xa(:))) > 1d-8) then
  !       write(*,1) 'initial_h1', initial_h1
  !       write(*,1) 'initial_h2', initial_h2
  !       write(*,1) 'initial_he3', initial_he3
  !       write(*,1) 'initial_he4', initial_he4
  !       write(*,1) 'initial_y', initial_y
  !       write(*,1) 'initial_z', initial_z
  !       write(*,1) 'initial_h1+h2+he3+he4+z', &
  !          initial_h1 + initial_h2 + initial_he3 + initial_he4 + initial_z
  !       write(*,1) 'sum(xa(:))', sum(xa(:))
  !       write(*,*) 'build_pre_ms_model'
  !       ierr = -1
  !       return
  !    end if

  !    call basic_composition_info( &
  !       species, s% chem_id, xa(:), x, y, z, abar, zbar, z2bar, ye, &
  !       mass_correction, sumx)

  !    mu_eff = 4 / (3 + 5*x) 
  !       ! estimate mu_eff assuming complete ionization and Z << 1

  !    call set_metals ! for opacities

  !    guess_rho_c = s% pre_ms_guess_rho_c
  !    if (guess_rho_c <= 0) then ! use n=3/2 polytrope
  !       rstar = Rsun*7.41d6*(mu_eff/0.6d0)*(mstar/Msun)/T_c ! Ushomirsky et al, 6
  !       rho_c = 8.44d0*(mstar/Msun)*pow3(Rsun/rstar) ! Ushomirsky et al, 5
  !       rho_c = 1.2d0*rho_c ! polytrope value is usually too small
  !    else
  !       rho_c = guess_rho_c
  !    end if

  !    ! pick a luminosity that is above the zams level
  !    lgM = log10_cr(mstar/Msun)
  !    if (lgM > 1) then
  !       lgL = 4.5d0 + 2.5d0*(lgM - 1)
  !    else if (lgM > 0.5d0) then
  !       lgL = 3d0 + 3*(lgM - 0.5d0)
  !    else if (lgM > 0) then
  !       lgL = 4d0*lgM + 0.5d0
  !    else
  !       lgL = 0.5d0
  !    end if

  !    ! use uniform eps_grav to give that luminosity
  !    eps_grav = exp10_cr(lgL)*Lsun/mstar

  !    if (dbg) then
  !       write(*,1) 'initial_z', initial_z
  !       write(*,1) 'T_c', T_c
  !       write(*,1) 'rho_c', rho_c
  !       write(*,1) 'lgL', lgL
  !       write(*,1) 'eps_grav', eps_grav
  !       write(*,1) 'mstar/Msun', mstar/Msun
  !    end if

  !    if (ASSOCIATED(s% xh)) deallocate(s% xh)
  !    nullify(s% q)
  !    nullify(s% dq)

  !    d_log10_P = s% pre_ms_d_log10_P

  !    i = 1 ! rpar(1) for mstar result
  !    rpar(i+1) = T_c; i = i+1
  !    rpar(i+1) = eps_grav; i = i+1
  !    rpar(i+1) = x; i = i+1
  !    rpar(i+1) = initial_z; i = i+1
  !    rpar(i+1) = abar; i = i+1
  !    rpar(i+1) = zbar; i = i+1
  !    rpar(i+1) = d_log10_P; i = i+1

  !    rpar(i+1:i+species) = xa(1:species); i = i+species
  !    rpar(i+1:i+nmet) = metals(1:nmet); i = i+nmet

  !    if (i /= pre_ms_lrpar) then
  !       write(*,*) 'i /= pre_ms_lrpar', i, pre_ms_lrpar
  !       write(*,*) 'pre ms'
  !       ierr = -1
  !       return
  !    end if

  !    ipar(1) = id

  !    lnd = log_cr(rho_c)
  !    dlnd = 0.01d0

  !    call look_for_brackets(lnd, dlnd, lnd1, lnd3, pre_ms_f, y1, y3, &
  !          imax, pre_ms_lrpar, rpar, pre_ms_lipar, ipar, ierr)
  !    if (ierr /= 0) then
  !       if (dbg) then
  !          if (dbg) write(*,*) 'look_for_brackets ierr', ierr
  !          write(*,1) 'lnd1', lnd1
  !          write(*,1) 'lnd3', lnd3
  !          write(*,1) 'y1', y1
  !          write(*,1) 'y3', y3
  !       end if
  !       return
  !    end if

  !    epsx = 1d-3 ! limit for variation in lnd
  !    epsy = 1d-3 ! limit for matching desired mass as fraction of total mass

  !    lnd = safe_root(pre_ms_f, lnd1, lnd3, y1, y3, imax, epsx, epsy, &
  !             pre_ms_lrpar, rpar, pre_ms_lipar, ipar, ierr)
  !    if (ierr /= 0) then
  !       if (dbg) write(*,*) 'safe_root ierr', ierr
  !       return
  !    end if

  !    mstar1 = rpar(1)

  !    xh => s% xh
  !    q => s% q
  !    dq => s% dq
  !    nz = s% nz

  !    if (dbg) then
  !       write(*,*)
  !       write(*,*) 'finished pre-MS model'
  !       write(*,1) 'mstar1/Msun', mstar1/Msun
  !       write(*,1) '(mstar-mstar1)/mstar', (mstar-mstar1)/mstar
  !       write(*,1) 'log10(r/Rsun)', log10_cr(exp_cr(xh(s% i_lnR,1))/Rsun)
  !       if (s% i_lum /= 0) write(*,1) 'log10(L/Lsun)', log10_cr(xh(s% i_lum,1)/Lsun)
  !       write(*,1) 'log10(Tsurf)', xh(s% i_lnT,1)/ln10
  !       write(*,1) 'Tsurf', exp_cr(xh(s% i_lnT,1))
  !       write(*,*) 'nz', nz
  !       write(*,*)
  !       stop 'debug: pre ms'
  !    end if

  !    ! The following deallocations deal with arrays which were
  !    ! needed in the root bracket/solve above, but will get
  !    ! overwritten during the call to allocate_star_info_arrays

  !    if (ASSOCIATED(s% xh_old)) deallocate(s% xh_old)
  !    if (ASSOCIATED(s% xh_older)) deallocate(s% xh_older)
  !    if (ASSOCIATED(s% equ1)) deallocate(s% equ1)
  !    if (ASSOCIATED(s% xh_pre)) deallocate(s% xh_pre)

  !    call allocate_star_info_arrays(s, ierr)
  !    if (ierr /= 0) then
  !       call dealloc
  !       return
  !    end if

  !    do k=1,nz
  !       do j=1,nvar_hydro
  !          s% xh(j,k) = xh(j,k)
  !          !write(*,3) trim(s% nameofvar(j)), j, k, xh(j,k)
  !       end do
  !       do j=1,species
  !          s% xa(j,k) = xa(j)
  !       end do
  !       s% q(k) = q(k)
  !       s% dq(k) = dq(k)
  !    end do

  !    call dealloc


  !    subroutine dealloc
  !       deallocate(xh, q, dq)
  !    end subroutine dealloc

  ! end subroutine build_pre_ms_model


  real(dp) function pre_ms_f(lnd, dfdx, lrpar, rpar, lipar, ipar, ierr)
     integer, intent(in) :: lrpar, lipar
     real(dp), intent(in) :: lnd
     real(dp), intent(out) :: dfdx
     integer, intent(inout), pointer :: ipar(:) ! (lipar)
     real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
     integer, intent(out) :: ierr

     type (star_info), pointer :: s
     real(dp) :: rho_c, T_c, eps_grav, x, z, abar, zbar, d_log10_P
     real(dp), pointer :: xa(:)
     integer :: i, nz, species
     real(dp) :: mstar, mstar1

     logical, parameter :: dbg = .true.

     include 'formats'

     ierr = 0
     pre_ms_f = 0
     if (lipar <= 0) then
        write(*,*) 'lipar', lipar
        write(*,*) 'pre_ms f'
        ierr = -1
        return
     end if

     call star_ptr(ipar(1), s, ierr)
     if (ierr /= 0) return

     species = s% species

     if (associated(s% xh)) deallocate(s% xh)
     if (associated(s% q)) deallocate(s% q)
     if (associated(s% dq)) deallocate(s% dq)

     rho_c = exp_cr(lnd)

     i = 1 ! rpar(1) for mstar result
     T_c = rpar(i+1); i = i+1
     eps_grav = rpar(i+1); i = i+1
     x = rpar(i+1); i = i+1
     z = rpar(i+1); i = i+1
     abar = rpar(i+1); i = i+1
     zbar = rpar(i+1); i = i+1
     xa => rpar(i+1:i+species); i = i+species
     if (i > lrpar) then
        write(*,*) 'i > lrpar', i, lrpar
        write(*,*) 'pre_ms f'
        ierr = -1
        return
     end if

     mstar = s% mstar ! desired value
     mstar1 = mstar ! to keep gfortran quiet

     if (dbg) write(*,*) 'call build1_pre_ms_model'
     call build1_pre_ms_model( &
           s, T_c, rho_c, 0d0, eps_grav, &
           x, s% initial_z, abar, zbar, &
           xa, nz, mstar1, ierr)
     if (ierr /= 0) then
        write(*,*) 'failed in build1_pre_ms_model'
        !ierr = 0
        return
     end if

     s% nz = nz

     rpar(1) = mstar1 ! return the actual mass

     pre_ms_f = (mstar - mstar1) / mstar
     dfdx = 0

     if (dbg) then
        write(*,1) 'rho_c', rho_c
        write(*,1) 'mstar1', mstar1 / msun
        write(*,1) 'pre_ms_f', pre_ms_f
        write(*,*)
     end if

  end function pre_ms_f


  subroutine build1_pre_ms_model( &
           s, T_c, rho_c, d_log10_P_in, eps_grav_in, &
           x, z, abar, zbar, xa, nz, mstar, ierr)
     use chem_def
     use eos_def
     use kap_lib
     use chem_lib
     use eos_lib, only: Radiation_Pressure
     use eos_support, only: get_eos, solve_eos_given_PgasT
!     use star_utils, only: normalize_dqs, set_qs
     type (star_info), pointer :: s
     real(dp), intent(in) :: &
        T_c, rho_c, d_log10_P_in, eps_grav_in, &
        x, z, abar, zbar, xa(:)
     integer, intent(out) :: nz
     real(dp), intent(out) :: mstar ! the mass of the constructed model
     integer, intent(out) :: ierr

     real(dp) :: logRho_guess
     real(dp), parameter :: LOGRHO_TOL = 1E-6_dp
     real(dp), parameter :: LOGPGAS_TOL = 1E-6_dp

     integer :: i, ii, k, j, i_lnd, i_lnT, i_lnR, prune, max_retries
     real(dp), parameter :: &
        delta_logPgas = 0.01d0, q_at_nz = 1d-5
     real(dp) :: &
        P_surf_limit, y, dlogPgas, logPgas, Prad, Pgas, try_dlogPgas, logPgas0, &
        res(num_eos_basic_results), eps_grav, P_c, logP, m, &
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
     
     logical, parameter :: dbg = .false.

     include 'formats'

     ierr = 0

     logP_surf_limit = s% pre_ms_logP_surf_limit
     if (logP_surf_limit <= 0) logP_surf_limit = 3.5d0
     P_surf_limit = exp10_cr(logP_surf_limit)

     logT_surf_limit = s% pre_ms_logT_surf_limit
     if (logT_surf_limit <= 0) logT_surf_limit = 3.7d0

     if (dbg) write(*,1) 'logT_surf_limit', logT_surf_limit

     i_lnd = s% i_lnd
     i_lnT = s% i_lnT
     i_lnR = s% i_lnR

     if (i_lnd == 0) then
        write(*,*) 'Sorry: require lnPgas_flag be .false. for create_pre_main_sequence'
        ierr = -1
        return
     end if

     if (i_lnT == 0) then
        write(*,*) 'Sorry: require E_flag be .false. for create_pre_main_sequence'
        ierr = -1
        return
     end if

     cgrav = standard_cgrav

     eps_grav = eps_grav_in
     if (dbg) write(*,1) 'eps_grav', eps_grav

     if (d_log10_P_in == 0) then
        dlogPgas = delta_logPgas
     else
        dlogPgas = abs(d_log10_P_in)
     end if
     if (dbg) write(*,1) 'dlogPgas', dlogPgas

     call get_eos( &
           s, 0, z, x, abar, zbar, xa, &
           rho_c, log10_cr(rho_c), T_c, log10_cr(T_c), &
           res, d_eos_dlnd, d_eos_dlnT, &
           d_eos_dabar, d_eos_dzbar, ierr)
     if (ierr /= 0) return
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
     nz = 1 + (logP - s% pre_ms_logP_surf_limit)/dlogPgas

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

     L = eps_grav*m ! L at nz
     
     call set_metals

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

              L = L0 + dm*eps_grav ! luminosity at point k
              Lmid = (L0+L)/2

              Pmid = (P+P0)/2

              chiRho_mid = (chiRho0 + chiRho)/2
              chiT_mid = (chiT0 + chiT)/2
              Cp_mid = (Cp0 + Cp)/2
              grada_mid = (grada0 + grada)/2

              do ii = 1, 2
                 Tmid = (T+T0)/2

                 call set_metals
                 
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
                 
                 if (kap_cond .lt. 8 * kap_rad) then
                    call eval_gradT( &
                       s, zbar, x, y, xa, rho_mid, mmid, mstar, rmid, Tmid, log_cr(Tmid), Lmid, Pmid, &
                       chiRho_mid, chiT_mid, Cp_mid, opacity, grada_mid, &
                       lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                       gradT, ierr )
                    if (ierr /= 0) return
                    write(*,*) gradT, res(i_eta)
                    !gradT = 0
                 else
                    gradT = 0.235
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

     call normalize_dqs(nz, dq, ierr)
     if (ierr /= 0) then
        if (s% report_ierr) write(*,*) 'set_qs failed in pre ms model'
        return
     end if
     call set_qs(nz, q, dq, ierr)
     if (ierr /= 0) then
        if (s% report_ierr) write(*,*) 'set_qs failed in pre ms model'
        return
     end if

     contains

    subroutine set_metals ! for opacities
      xc = xa(s% net_iso(ic12))
      xn = xa(s% net_iso(in14))
      xo = xa(s% net_iso(io16))
      xne = xa(s% net_iso(ine20))
      Zbase = s% Zbase
    end subroutine set_metals

       

     subroutine unpack_eos_results
        chiRho = res(i_chiRho)
        chiT = res(i_chiT)
        Cp = res(i_cp)
        grada = res(i_grad_ad)
        lnfree_e = res(i_lnfree_e)
        d_lnfree_e_dlnRho = d_eos_dlnd(i_lnfree_e)
        d_lnfree_e_dlnT = d_eos_dlnT(i_lnfree_e)
     end subroutine unpack_eos_results


  end subroutine build1_pre_ms_model


  subroutine eval_gradT( &
        s, zbar, x, y, xa, rho, m, mstar, r, T, lnT, L, P, &
        chiRho, chiT, Cp, opacity, grada, &
        lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
        gradT, ierr )
     use chem_def, only: ih1
     use my_mlt_info, only: do1_mlt_eval
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

     call do1_mlt_eval( &
        s, 0, cgrav, m, mstar, r, L, x, T, rho, P, &
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

  subroutine normalize_dqs(nz, dq, ierr)
         ! rescale dq's so that add to 1.000
         ! work in from boundaries to meet at largest dq
         integer, intent(in) :: nz
         real(dp), intent(inout) :: dq(:) ! (nz)
         integer, intent(out) :: ierr
         integer :: k, midq
         real(dp) :: dqsum1, dqsum2, dq_min
         include 'formats'
         k = minloc(dq(1:nz),dim=1)
         dq_min = dq(k)
         if (dq_min <= 0d0) then
            write(*,2) 'bad dq', k, dq(k)
            ierr = -1
            !stop
            return
         end if
         midq = maxloc(dq(1:nz),dim=1)
         ! surface inward
         dqsum1 = 0
         do k=1, midq
            dqsum1 = dqsum1 + dq(k)
            if (dq(k) <= 0) then
               ierr = -1
               return
            end if
         end do
         ! center outward
         dqsum2 = 0
         do k=nz, midq+1, -1
            dqsum2 = dqsum2 + dq(k)
            if (dq(k) <= 0) then
               ierr = -1
               return
            end if
         end do
         do k=1,nz
            dq(k) = dq(k)/(dqsum1 + dqsum2)
         end do
      end subroutine normalize_dqs


      subroutine set_qs(nz, q, dq, ierr) ! set q's using normalized dq's
         integer, intent(in) :: nz
         real(dp), intent(inout) :: dq(:) ! (nz)
         real(dp), intent(inout) :: q(:) ! (nz)
         integer, intent(out) :: ierr
         integer :: k, midq
         real(dp) :: dqsum1, dqsum2
         logical :: okay
         include 'formats'
         ierr = 0
         call normalize_dqs(nz, dq, ierr)
         if (ierr /= 0) return
         q(1) = 1d0
         okay = .true.
         do k=2,nz
            q(k) = q(k-1) - dq(k-1)
            if (q(k) < 0d0 .or. q(k) > 1d0) then
               okay = .false.
               exit
            end if
         end do
         if (okay) return
         midq = maxloc(dq(1:nz),dim=1)
         ! surface inward
         dqsum1 = 0
         do k=1, midq
            q(k) = 1d0 - dqsum1
            dqsum1 = dqsum1 + dq(k)
         end do
         ! center outward
         dqsum2 = 0
         do k=nz, midq+1, -1
            dqsum2 = dqsum2 + dq(k)
            q(k) = dqsum2
         end do         
      end subroutine set_qs



      
end module model_builder
