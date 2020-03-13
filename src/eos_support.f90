! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

module my_eos_support

  ! Uses

  use const_def
  use star_def
  use utils_lib, only : is_bad

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: MAX_ITER_FOR_SOLVE = 100

  ! Access specifiers

  private

  public :: get_eos
  public :: solve_eos_given_DE
  public :: solve_eos_given_DP
  public :: solve_eos_given_DS
  public :: solve_eos_given_PT
  public :: solve_eos_given_PgasT
  public :: solve_eos_given_PgasT_auto

  public :: get_peos
  public :: solve_peos_given_DT
  public :: solve_peos_given_DT_auto
  
  ! Procedures

contains

  ! Get eos results data given density & temperature

  subroutine get_eos( &
       s, k, z, xh, abar, zbar, xa, &
       Rho, logRho, T, logT, &
       res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, ierr)

    use eos_lib, only: eosDT_get, eos_gamma_DT_get, &
         eosDT_ideal_gas_get, eosDT_HELMEOS_get
    use eos_def, only: num_eos_basic_results, num_helm_results

    type (star_info), pointer :: s
    integer, intent(in) :: k ! 0 means not being called for a particular cell
    real(dp), intent(in) :: &
         z, xh, abar, zbar, xa(:), Rho, logRho, T, logT
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar
    integer, intent(out) :: ierr

    real(dp) :: eos_z, eos_x, Pgas, Prad, energy, entropy
    real(dp) :: helm_res(num_helm_results)
    integer :: j
    logical :: off_table

    include 'formats'

    ierr = 0
    
    if (s% use_fixed_XZ_for_eos) then
       eos_x = s% fixed_X_for_eos
       eos_z = s% fixed_Z_for_eos
    else
       eos_x = xh
       eos_z = z
    end if

    if (s% use_other_eos) then
       call s% other_eosDT_get( &
            s% id, k, s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            Rho, logRho, T, logT, &
            res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, ierr)
    else if (s% gamma_law_hydro > 0d0) then
       call eos_gamma_DT_get( &
            s% eos_handle, abar, rho, logRho, T, logT, s% gamma_law_hydro, &
            res, dres_dlnRho, dres_dlnT, &
            Pgas, Prad, energy, entropy, ierr)
       if (ierr /= 0) then
          !$OMP critical (get_eos_critical)
          write(*,2) 'eos_gamma_DT_get failed in get_eos', k
          write(*,1) 'abar', abar
          write(*,1) 'rho', rho
          write(*,1) 'log10Rho', logRho
          write(*,1) 'T', T
          write(*,1) 'log10T', logT
          write(*,1) 's% gamma_law_hydro', s% gamma_law_hydro
          !$OMP end critical (get_eos_critical)
       end if
    else if (s% use_eosDT_ideal_gas) then
       call eosDT_ideal_gas_get( &
            s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            Rho, logRho, T, logT, &
            res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, ierr)
    else if (s% use_eosDT_HELMEOS) then
       call eosDT_HELMEOS_get( &
            s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            Rho, logRho, T, logT, &
            s% eosDT_HELMEOS_include_radiation, &
            s% eosDT_HELMEOS_always_skip_elec_pos, &
            s% eosDT_HELMEOS_always_include_elec_pos, &
            s% eos_rq% logT_ion_HELM, s% eos_rq% logT_neutral_HELM, &
            res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, helm_res, off_table, ierr)
    else
       call eosDT_get( &
            s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            Rho, logRho, T, logT, &
            res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, ierr)
    end if

    if (ierr /= 0) then
       if (s% report_ierr) then
          !$OMP critical (get_eos_critical)
          write(*,*) 'get_eos ierr', ierr
          write(*,2) 'k', k
          write(*,1) 'z', z
          write(*,1) 'xh', xh
          write(*,1) 'abar', abar
          write(*,1) 'zbar', zbar
          write(*,1) 'log10Rho', logRho
          write(*,1) 'log10T', logT
          if (s% stop_for_bad_nums .and. &
               is_bad(logRho+logT)) stop 'do_eos_for_cell'
          !$OMP end critical (get_eos_critical)
       end if
       return
    end if

  end subroutine get_eos

  !****

  ! Get Pgas-tabulated eos results data given gas pressure &
  ! temperature

  subroutine get_peos ( &
       s, k, z, xh, abar, zbar, xa, &
       Pgas, logPgas, T, logT, &
       Rho, logRho, dlnRho_dlnPgas, dlnRho_dlnT, &
       res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, ierr)

    use eos_lib, only: eosPT_get, eos_gamma_PT_get
    use eos_def, only: num_eos_basic_results

    type (star_info), pointer :: s
    integer, intent(in) :: k ! 0 means not being called for a particular cell
    real(dp), intent(in) :: &
         z, xh, abar, zbar, xa(:), Pgas, logPgas, T, logT

    real(dp), intent(out) :: &
         Rho, logRho, dlnRho_dlnPgas, dlnRho_dlnT
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar
    integer, intent(out) :: ierr

    real(dp) :: eos_z, eos_x

    include 'formats'

    ierr = 0

    if (s% use_fixed_XZ_for_eos) then
       eos_x = s% fixed_X_for_eos
       eos_z = s% fixed_Z_for_eos
    else
       eos_x = xh
       eos_z = z
    end if

    if (s% use_other_eos) then
       call s% other_eosPT_get( &
            s% id, k, s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            Pgas, logPgas, T, logT, &
            Rho, logRho, dlnRho_dlnPgas, dlnRho_dlnT, &
            res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, ierr)
    else if (s% gamma_law_hydro > 0d0) then
       call eos_gamma_PT_get( &
            s% eos_handle, abar, Pgas, logPgas, T, logT, s% gamma_law_hydro, &
            Rho, logRho, res, dres_dlnRho, dres_dlnT, ierr)
       dlnRho_dlnPgas = 0._dp
       dlnRho_dlnT = 0._dp
    else
       call eosPT_get( &
            s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            Pgas, logPgas, T, logT, &
            Rho, logRho, dlnRho_dlnPgas, dlnRho_dlnT, &
            res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, ierr)
    end if

    if (ierr /= 0) then
       if (s% report_ierr) then
          !$OMP critical (get_peos_critical)
          write(*,*) 'eval_eosPT ierr', ierr
          write(*,2) 'k', k
          write(*,1) 'z', z
          write(*,1) 'xh', xh
          write(*,1) 'abar', abar
          write(*,1) 'zbar', zbar
          write(*,1) 'log10Pgas', logPgas
          write(*,1) 'log10T', logT
          !$OMP end critical (get_peos_critical)
       end if
       return
    end if

  end subroutine get_peos

  !****

  ! Solve for temperature & eos results data given density & energy

  subroutine solve_eos_given_DE( &
       s, k, z, xh, abar, zbar, xa, &
       logRho, logE, logT_guess, logT_tol, logE_tol, &
       logT, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
       ierr)

    use eos_def
    use eos_lib, only: eosDT_get_T

    type (star_info), pointer :: s
    integer, intent(in) :: k ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         z, xh, abar, zbar, xa(:), logRho, logE, &
         logT_guess, logT_tol, logE_tol
    real(dp), intent(out) :: logT
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar
    integer, intent(out) :: ierr

    integer :: eos_calls
    real(dp) :: eos_x, eos_z

    include 'formats'

    ierr = 0

    if (s% use_fixed_XZ_for_eos) then
       eos_x = s% fixed_X_for_eos
       eos_z = s% fixed_Z_for_eos
    else
       eos_x = xh
       eos_z = z
    end if

    if (s% use_other_eos) then
       call s% other_eosDT_get_T( &
            s% id, k, s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            logRho, i_lnE, logE*ln10, &
            logT_tol, logE_tol*ln10, MAX_ITER_FOR_SOLVE, logT_guess, &
            arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
            logT, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
            eos_calls, ierr)
    else if (s% gamma_law_hydro > 0d0) then
       write(*,*) 'cannot call solve_eos_given_DE with gamma_law_hydro > 0'
       ierr = -1
       return
    else if (s% use_eosDT_ideal_gas) then
       write(*,*) 'cannot call solve_eos_given_DE with use_eosDT_ideal_gas set'
       ierr = -1
       return
    else if (s% use_eosDT_HELMEOS) then
       write(*,*) 'cannot call solve_eos_given_DE with use_eosDT_HELMEOS set'
       ierr = -1
       return
    else
       call eosDT_get_T( &
            s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            logRho, i_lnE, logE*ln10, &
            logT_tol, logE_tol*ln10, MAX_ITER_FOR_SOLVE, logT_guess,  &
            arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
            logT, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
            eos_calls, ierr)
    end if

  end subroutine solve_eos_given_DE

  !****

  ! Solve for temperature & eos results data given density & pressure

  subroutine solve_eos_given_DP( &
       s, k, z, xh, abar, zbar, xa, &
       logRho, logP, logT_guess, logT_tol, logP_tol, &
       logT, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
       ierr)

    use eos_def
    use eos_lib, only: eosDT_get_T_given_Ptotal

    type (star_info), pointer :: s
    integer, intent(in) :: k ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         z, xh, abar, zbar, xa(:), logRho, logP, &
         logT_guess, logT_tol, logP_tol
    real(dp), intent(out) :: logT
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar
    integer, intent(out) :: ierr

    integer :: eos_calls
    real(dp) :: eos_x, eos_z

    include 'formats'

    ierr = 0

    if (s% use_fixed_XZ_for_eos) then
       eos_x = s% fixed_X_for_eos
       eos_z = s% fixed_Z_for_eos
    else
       eos_x = xh
       eos_z = z
    end if

    if (s% use_other_eos) then
       write(*,*) 'cannot call solve_eos_given_DP with use_other_eos set'
       ierr = -1
       return
    else if (s% gamma_law_hydro > 0d0) then
       write(*,*) 'cannot call solve_eos_given_DP with gamma_law_hydro > 0'
       ierr = -1
       return
    else if (s% use_eosDT_ideal_gas) then
       write(*,*) 'cannot call solve_eos_given_DP with use_eosDT_ideal_gas set'
       ierr = -1
       return
    else if (s% use_eosDT_HELMEOS) then
       write(*,*) 'cannot call solve_eos_given_DP with use_eosDT_HELMEOS set'
       ierr = -1
       return
    else
       call eosDT_get_T_given_Ptotal( &
            s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            logRho, logP, logT_tol, logP_tol, MAX_ITER_FOR_SOLVE, logT_guess, &
            arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
            logT, res, dres_dlnRho, dres_dlnT, &
            dres_dabar, dres_dzbar, eos_calls, ierr)
    end if
          
  end subroutine solve_eos_given_DP

  !****

  ! Solve for temperature & eos results data for a given density &
  ! entropy

  subroutine solve_eos_given_DS( &
       s, k, z, xh, abar, zbar, xa, &
       logRho, logS, logT_guess, logT_tol, logS_tol, &
       logT, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
       ierr)

    use eos_def
    use eos_lib, only: eosDT_get_T

    type (star_info), pointer :: s
    integer, intent(in) :: k ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         z, xh, abar, zbar, xa(:), logRho, logS, &
         logT_guess, logT_tol, logS_tol
    real(dp), intent(out) :: logT
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar
    integer, intent(out) :: ierr

    integer :: eos_calls
    real(dp) :: eos_x, eos_z

    include 'formats'

    ierr = 0
    
    if (s% use_fixed_XZ_for_eos) then
       eos_x = s% fixed_X_for_eos
       eos_z = s% fixed_Z_for_eos
    else
       eos_x = xh
       eos_z = z
    end if

    if (s% use_other_eos) then
       call s% other_eosDT_get_T( &
            s% id, k, s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            logRho, i_lnS, logS*ln10, &
            logT_tol, logS_tol*ln10, MAX_ITER_FOR_SOLVE, logT_guess, &
            arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
            logT, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
            eos_calls, ierr)
    else if (s% gamma_law_hydro > 0d0) then
       write(*,*) 'cannot call solve_eos_given_DS with gamma_law_hydro > 0'
       ierr = -1
       return
    else if (s% use_eosDT_ideal_gas) then
       write(*,*) 'cannot call solve_eos_given_DS with use_eosDT_ideal_gas set'
       ierr = -1
       return
    else if (s% use_eosDT_HELMEOS) then
       write(*,*) 'cannot call solve_eos_given_DS with use_eosDT_HELMEOS set'
       ierr = -1
       return
    else
       call eosDT_get_T( &
            s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            logRho, i_lnS, logS*ln10, &
            logT_tol, logS_tol*ln10, MAX_ITER_FOR_SOLVE, logT_guess,  &
            arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
            logT, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
            eos_calls, ierr)
    end if

  end subroutine solve_eos_given_DS

  !****

  ! Solve for density & eos results data given pressure & temperature

  subroutine solve_eos_given_PT( &
       s, k, z, xh, abar, zbar, xa, &
       logT, logP, logRho_guess, logRho_tol, logP_tol, &
       logRho, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
       ierr)

    use eos_def
    use eos_lib, only: eosDT_get_Rho_given_Ptotal

    type (star_info), pointer :: s
    integer, intent(in) :: k ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         z, xh, abar, zbar, xa(:), logT, logP, &
         logRho_guess, logRho_tol, logP_tol
    real(dp), intent(out) :: logRho
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar
    integer, intent(out) :: ierr

    integer :: eos_calls
    real(dp) :: eos_x, eos_z

    include 'formats'

    ierr = 0

    if (s% use_fixed_XZ_for_eos) then
       eos_x = s% fixed_X_for_eos
       eos_z = s% fixed_Z_for_eos
    else
       eos_x = xh
       eos_z = z
    end if

    if (s% use_other_eos) then
       write(*,*) 'cannot call solve_eos_given_PT with use_other_eos set'
       ierr = -1
       return
    else if (s% gamma_law_hydro > 0d0) then
       write(*,*) 'cannot call solve_eos_given_PT with gamma_law_hydro > 0'
       ierr = -1
       return
    else if (s% use_eosDT_ideal_gas) then
       write(*,*) 'cannot call solve_eos_given_PT with use_eosDT_ideal_gas set'
       ierr = -1
       return
    else if (s% use_eosDT_HELMEOS) then
       write(*,*) 'cannot call solve_eos_given_PT with use_eosDT_HELMEOS set'
       ierr = -1
       return
    else
       call eosDT_get_Rho_given_Ptotal( &
            s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            logT, logP, logRho_tol, logP_tol, MAX_ITER_FOR_SOLVE, logRho_guess, &
            arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
            logRho, res, dres_dlnRho, dres_dlnT, &
            dres_dabar, dres_dzbar, eos_calls, ierr)
    end if

  end subroutine solve_eos_given_PT

  !****

  ! Solve for density & eos results data given gas pressure &
  ! temperature

  subroutine solve_eos_given_PgasT( &
       s, k, z, xh, abar, zbar, xa, &
       logT, logPgas, logRho_guess, logRho_tol, logPgas_tol, &
       logRho, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
       ierr)

    use eos_def
    use eos_lib, only: eosDT_get_Rho

    type (star_info), pointer :: s
    integer, intent(in) :: k ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         z, xh, abar, zbar, xa(:), logT, logPgas, &
         logRho_guess, logRho_tol, logPgas_tol
    real(dp), intent(out) :: logRho
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar
    integer, intent(out) :: ierr

    integer :: eos_calls
    real(dp) :: eos_x, eos_z

    include 'formats'

    ierr = 0

    if (s% use_fixed_XZ_for_eos) then
       eos_x = s% fixed_X_for_eos
       eos_z = s% fixed_Z_for_eos
    else
       eos_x = xh
       eos_z = z
    end if

    if (s% use_other_eos) then
       call s% other_eosDT_get_Rho( &
            s% id, k, s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            logT, i_lnPgas, logPgas*ln10, &
            logRho_tol, logPgas_tol*ln10, MAX_ITER_FOR_SOLVE, logRho_guess, &
            arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
            logRho, res, dres_dlnRho, dres_dlnT, &
            dres_dabar, dres_dzbar, eos_calls, ierr)
    else if (s% gamma_law_hydro > 0d0) then
       write(*,*) 'cannot call solve_eos_given_PgasT with gamma_law_hydro > 0'
       ierr = -1
       return
    else if (s% use_eosDT_ideal_gas) then
       write(*,*) 'cannot call solve_eos_given_PgasT with use_eosDT_ideal_gas set'
       ierr = -1
       return
    else if (s% use_eosDT_HELMEOS) then
       write(*,*) 'cannot call solve_eos_given_PgasT with use_eosDT_HELMEOS set'
       ierr = -1
       return
    else
       call eosDT_get_Rho( &
            s% eos_handle, eos_z, eos_x, abar, zbar,  &
            s% species, s% chem_id, s% net_iso, xa, &
            logT, i_lnPgas, logPgas*ln10, &
            logRho_tol, logPgas_tol*ln10, MAX_ITER_FOR_SOLVE, logRho_guess, &
            arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
            logRho, res, dres_dlnRho, dres_dlnT, &
            dres_dabar, dres_dzbar, eos_calls, ierr)
    end if

  end subroutine solve_eos_given_PgasT

  !****

  ! Solve for density & eos results data given gas pressure &
  ! temperature, with logRho_guess calculated automatically via an
  ! initial call to get_peos

  subroutine solve_eos_given_PgasT_auto( &
       s, k, z, xh, abar, zbar, xa, &
       logT, logPgas, logRho_tol, logPgas_tol, &
       logRho, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
       ierr)

    use eos_def

    type (star_info), pointer :: s
    integer, intent(in) :: k ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         z, xh, abar, zbar, xa(:), logT, logPgas, &
         logRho_tol, logPgas_tol
    real(dp), intent(out) :: logRho
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar
    integer, intent(out) :: ierr

    real(dp) :: rho_guess, logRho_guess, &
         dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas
    
    call get_peos ( &
       s, k, z, xh, abar, zbar, xa, &
       arg_not_provided, logPgas, arg_not_provided, logT, &
       rho_guess, logRho_guess, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
       res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, ierr)

    call solve_eos_given_PgasT( &
       s, k, z, xh, abar, zbar, xa, &
       logT, logPgas, logRho_guess, logRho_tol, logPgas_tol, &
       logRho, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
       ierr)

  end subroutine solve_eos_given_PgasT_auto
         
  !****

  ! Solve for Pgas-tabulated gas pressure & eos results data given
  ! density & temperature

  subroutine solve_peos_given_DT( &
       s, k, z, xh, abar, zbar, xa, &
       logT, logRho, logPgas_guess, logPgas_tol, logRho_tol, &
       logPgas, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
       ierr)

    use eos_def
    use eos_lib, only: eosPT_get_Pgas_for_Rho

    type (star_info), pointer :: s
    integer, intent(in) :: k ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         z, xh, abar, zbar, xa(:), logT, logRho, &
         logPgas_guess, logPgas_tol, logRho_tol
    real(dp), intent(out) :: logPgas
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar
    integer, intent(out) :: ierr

    integer :: eos_calls
    real(dp) :: eos_x, eos_z, rho, logRho_result, &
         dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas

    include 'formats'

    ierr = 0

    if (s% use_fixed_XZ_for_eos) then
       eos_x = s% fixed_X_for_eos
       eos_z = s% fixed_Z_for_eos
    else
       eos_x = xh
       eos_z = z
    end if

    if (s% use_other_eos) then
       call s% other_eosPT_get_Pgas_for_Rho( &
            s% id, k, s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            logT, logRho, logPgas_tol, logRho_tol, MAX_ITER_FOR_SOLVE, logPgas_guess, &
            arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
            logPgas, rho, logRho_result, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
            res, dres_dlnRho, dres_dlnT, &
            dres_dabar, dres_dzbar, eos_calls, ierr)
    else if (s% gamma_law_hydro > 0d0) then
       write(*,*) 'cannot call solve_peos_given_DT with gamma_law_hydro > 0'
       ierr = -1
       return
    else
       call eosPT_get_Pgas_for_Rho( &
            s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            logT, logRho, logPgas_tol, logRho_tol, MAX_ITER_FOR_SOLVE, logPgas_guess, &
            arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
            logPgas, rho, logRho_result, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
            res, dres_dlnRho, dres_dlnT, &
            dres_dabar, dres_dzbar, eos_calls, ierr)
    end if

  end subroutine solve_peos_given_DT

  !****

  ! Solve for Pgas-tabulated gas pressure & eos results data given
  ! density & temperature, with logPgas_guess calculated automatically
  ! via an initial call to get_eos

  subroutine solve_peos_given_DT_auto( &
       s, k, z, xh, abar, zbar, xa, &
       logT, logRho, logPgas_tol, logRho_tol, &
       logPgas, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
       ierr)

    use eos_def

    type (star_info), pointer :: s
    integer, intent(in) :: k ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         z, xh, abar, zbar, xa(:), logT, logRho, &
         logPgas_tol, logRho_tol
    real(dp), intent(out) :: logPgas
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar
    integer, intent(out) :: ierr

    call get_eos ( &
       s, k, z, xh, abar, zbar, xa, &
       arg_not_provided, logRho, arg_not_provided, logT, &
       res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, ierr)

    call solve_peos_given_DT( &
       s, k, z, xh, abar, zbar, xa, &
       logT, logRho, res(i_lnPgas)/ln10, logPgas_tol, logRho_tol, &
       logPgas, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
       ierr)

  end subroutine solve_peos_given_DT_auto

end module my_eos_support
