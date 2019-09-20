! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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

      module star_utils

      use star_def
      use const_def
      use num_lib
      use utils_lib

      implicit none


      contains


      subroutine foreach_cell(s,nzlo,nzhi,use_omp,do1,ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         logical, intent(in) :: use_omp
         interface
            subroutine do1(s,k,ierr)
               use star_def
               type (star_info), pointer :: s
               integer, intent(in) :: k
               integer, intent(out) :: ierr
            end subroutine do1
         end interface
         integer, intent(out) :: ierr

         integer :: k, op_err
         logical :: okay
         ierr = 0

         if (nzlo == nzhi) then
            call do1(s,nzlo,ierr)
            return
         end if

         if (use_omp) then
            okay = .true.
!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(guided)
            do k = nzlo, nzhi
               if (.not. okay) cycle
               op_err = 0
               call do1(s,k,op_err)
               if (op_err /= 0) okay = .false. ! cannot just exit from a parallel loop
            end do
!$OMP END PARALLEL DO
            if (.not. okay) ierr = -1
         else
            do k = nzlo, nzhi
               call do1(s,k,ierr)
               if (ierr /= 0) exit
            end do
         end if

      end subroutine foreach_cell


      subroutine get_average_Y_and_Z(s, nzlo, nzhi, y_avg, z_avg, ierr)
         use chem_def
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         real(dp), intent(out) :: y_avg, z_avg
         integer, intent(out) :: ierr
         integer :: k, nz,  h1, h2, he3, he4
         real(dp) :: total_mass_h, total_mass_he, total_mass_z, &
            cell_mass, total_mass
         ierr = 0
         nz = s% nz
         h1 = s% net_iso(ih1)
         h2 = s% net_iso(ih2)
         he3 = s% net_iso(ihe3)
         he4 = s% net_iso(ihe4)
         total_mass=0; total_mass_h=0; total_mass_he=0; total_mass_z=0
         do k=nzlo, nzhi
            cell_mass = s% dm(k)
            total_mass = total_mass + cell_mass
            total_mass_h = total_mass_h + cell_mass*s% xa(h1, k)
            if (h2 /= 0) total_mass_h = total_mass_h + cell_mass*s% xa(h2, k)
            total_mass_he = total_mass_he + cell_mass*s% xa(he4, k)
            if (he3 /= 0) total_mass_he = total_mass_he + cell_mass*s% xa(he3, k)
         end do
         total_mass_z = total_mass - (total_mass_h + total_mass_he)
         z_avg = total_mass_z / total_mass
         y_avg = total_mass_he / total_mass
      end subroutine get_average_Y_and_Z


      real(dp) function eval_current_y(s, nzlo, nzhi, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr
         real(dp) :: y_avg, z_avg
         call get_average_Y_and_Z(s, nzlo, nzhi, y_avg, z_avg, ierr)
         eval_current_y = y_avg
      end function eval_current_y


      real(dp) function eval_current_z(s, nzlo, nzhi, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr
         real(dp) :: y_avg, z_avg
         call get_average_Y_and_Z(s, nzlo, nzhi, y_avg, z_avg, ierr)
         eval_current_z = z_avg
      end function eval_current_z


      real(dp) function eval_current_abundance(s, j, nzlo, nzhi, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: j, nzlo, nzhi
         integer, intent(out) :: ierr
         integer :: k, nz
         real(dp) :: cell_mass, jmass, total_mass

         ierr = 0

         if (j == 0) then
            eval_current_abundance = 0
            return
         end if

         nz = s% nz
         total_mass=0; jmass=0
         do k=nzlo, nzhi
            cell_mass = s% dm(k)
            total_mass = total_mass + cell_mass
            jmass = jmass + cell_mass*s% xa(j, k)
         end do
         eval_current_abundance = jmass / total_mass

      end function eval_current_abundance


      subroutine smooth_abundances(s, cnt, nzlo, nzhi, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: cnt ! make this many passes
         integer, intent(in) :: nzlo, nzhi ! only smooth zones nzlo to nzhi inclusive
         integer, intent(out) :: ierr
         integer :: k, j, nz
         ierr = 0
         nz = s% nz
         do j = 1, cnt
            do k = max(nzlo,2), min(nzhi, nz)
               s% xa(:,k) = (s% xa(:,k-1) + s% xa(:,k) + s% xa(:,k+1))/3
            end do
            if (nzhi == nz) s% xa(:,nz) = (s% xa(:,nz-1) + s% xa(:,nz) + s% xa(:,nz))/3
            if (nzlo == 1) s% xa(:,1) = (s% xa(:,2) + s% xa(:,1) + s% xa(:,1))/3
         end do
      end subroutine smooth_abundances


      integer function k_for_q(s, q)
         ! return k s.t. q(k) >= q > q(k)-dq(k)
         type (star_info), pointer :: s
         real(dp), intent(in) :: q
         integer :: k, nz
         nz = s% nz
         if (q >= 1) then
            k_for_q = 1; return
         else if (q <= s% q(nz)) then
            k_for_q = nz; return
         end if
         do k = 1, nz-1
            if (q > s% q(k+1)) then
               k_for_q = k; return
            end if
         end do
         k_for_q = nz
      end function k_for_q


      subroutine get_name_for_restart_file(n, num_digits, num)
         integer, intent(in) :: n, num_digits
         character (len=*), intent(out) :: num
         call get_string_for_model_number('x', n, num_digits, num)
      end subroutine get_name_for_restart_file


      subroutine get_string_for_model_number(prefix, n, num_digits, num)
         character (len=*), intent(in) :: prefix
         integer, intent(in) :: n, num_digits
         character (len=*), intent(out) :: num
         integer :: val
         character (len=32) :: fstring
         include 'formats'
         val = mod(n, 10**num_digits) ! wrap around
         if (val == 0) then
            write(num,*) n
            num = adjustl(num)
            return
         end if
        write(fstring,'( "(a,i",i2.2,".",i2.2,")" )') num_digits, num_digits
        write(num,fstring) trim(prefix), val
      end subroutine get_string_for_model_number


      subroutine report_xa_bad_nums(s,ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: k, j
         ierr = 0
         do k=1,s% nz
            do j=1,s% species
               if (is_bad(s% xa(j,k))) then
                  ierr = -1
                  write(*,*) j, k, s% xa(j,k)
               end if
            end do
         end do
      end subroutine report_xa_bad_nums


      real(dp) function eval_csound(s,k,ierr) result(cs)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         real(dp) :: cs2
         include 'formats'
         ierr = 0
         if (s% use_sr_sound_speed) then
            cs2 = s% gamma1(k)/(1d0 + (s% energy(k) + clight*clight)*s% rho(k)/s% P(k))
            if (cs2 < 0d0) then
               cs = 0d0
               ierr = -1
               return
            end if
            cs = clight*sqrt(cs2)
         else
            cs2 = s% gamma1(k)*s% P(k)/s% rho(k)
            if (cs2 < 0d0) then
               cs = 0d0
               ierr = -1
               return
            end if
            cs = sqrt(cs2)
         end if
      end function eval_csound


      subroutine set_m_grav_and_grav(s) ! using mass_corrections
         type (star_info), pointer :: s
         integer :: k, nz
         real(dp) :: twoGmrc2
         include 'formats'
         nz = s% nz
         if (.not. s% use_mass_corrections) then
            do k=1,nz
               s% m_grav(k) = s% m(k)
            end do
         else
            s% m_grav(nz) = &
               s% M_center + s% dm(nz)*s% mass_correction(nz)
            do k=nz-1,1,-1
               s% m_grav(k) = &
                  s% m_grav(k+1) + s% dm(k)*s% mass_correction(k)
            end do
         end if
         do k=1,nz
            s% grav(k) = s% cgrav(k)*s% m_grav(k)/(s% r(k)*s% r(k))
            if (s% use_gr_factors) then ! GR gravity factor = 1/sqrt(1-2Gm/(rc^2))
               twoGmrc2 = 2*s% cgrav(k)*s% m_grav(k)/(s% r(k)*clight*clight)
               s% grav(k) = s% grav(k)/sqrt(1d0 - twoGmrc2)
            end if
         end do
      end subroutine set_m_grav_and_grav


      subroutine use_xh_to_set_rho_to_dm_div_dV(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: k, nz, i_lnR, i_lnd
         real(dp) :: rL, rR, dm, dV, rho, old_lnd, new_lnd
         include 'formats'
         ierr = 0
         i_lnR = s% i_lnR
         i_lnd = s% i_lnd
         if (i_lnd == 0 .or. i_lnR == 0) return
         nz = s% nz
         rR = s% R_center
         do k = nz, 1, -1
            rL = rR
            rR = exp_cr(s% xh(i_lnR,k))
            dm = s% dm(k)
            dV = 4d0*pi/3d0*(rR*rR*rR - rL*rL*rL)
            rho = dm/dV
            if (rho <= 0d0) then
               write(*,3) 'set_rho_to_dm_div_dV: rho <= 0', &
                  k, nz, rho, dm, dV, rR, rL
               ierr = -1
               return
            end if
            new_lnd = log_cr(rho)
            s% xh(i_lnd,k) = new_lnd
         end do
      end subroutine use_xh_to_set_rho_to_dm_div_dV


      subroutine set_m_and_dm(s)
         type (star_info), pointer :: s
         integer :: k
         include 'formats'
         do k = 1, s% nz
            s% m(k) = s% M_center + s% q(k)*s% xmstar
            s% dm(k) = s% dq(k)*s% xmstar
            if (s% dm(k) <= 0d0 .or. is_bad(s% m(k) + s% dm(k))) then
               write(*,2) 'dm m dq q M_center', k, &
                  s% dm(k), s% m(k), s% dq(k), s% q(k), s% M_center
               if (.not. s% mass_flag) then
                  stop 'set_m_and_dm'
               end if
            end if
         end do
      end subroutine set_m_and_dm


      subroutine set_dm_bar(s, nz, dm, dm_bar)
         type (star_info), pointer :: s
         integer, intent(in) :: nz
         real(dp), intent(in) :: dm(:) ! (nz)
         real(dp), intent(inout) :: dm_bar(:) ! (nz)
         integer :: k
         do k=2,nz-1
            dm_bar(k) = 0.5d0*(dm(k-1) + dm(k))
         end do
         dm_bar(1) = 0.5d0*dm(1)
         if (s% rsp_flag) then ! rsp uses this definition
            dm_bar(nz) = 0.5d0*(dm(nz-1) + dm(nz))
         else
            dm_bar(nz) = 0.5d0*dm(nz-1) + dm(nz)
         end if
      end subroutine set_dm_bar


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
            stop
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


      subroutine set_xqs(nz, xq, dq, ierr) ! set xq's using dq's
         integer, intent(in) :: nz
         real(dp), intent(inout) :: dq(:) ! (nz)
         real(dp), intent(inout) :: xq(:) ! (nz)
         integer, intent(out) :: ierr
         integer :: k
         include 'formats'
         ierr = 0
         xq(1) = 0
         do k=2,nz-1
            xq(k) = xq(k-1) + dq(k-1)
         end do
         xq(nz) = 1 - dq(nz)
         if (xq(nz) < xq(nz-1)) then
            xq(nz) = xq(nz-1) + dq(nz-1)
            dq(nz) = 1 - xq(nz)
            if (dq(nz) <= 0) then
               ierr = -1
               return
            end if
         end if
      end subroutine set_xqs

      subroutine set_mass_eqn_edge_dqs(s)
         type (star_info), pointer :: s
         integer :: k, num_avg
         include 'formats'

         !TODO: describe

         !validate q limits
         if (s% mass_eqn_lower_q_full_off > s% mass_eqn_upper_q_full_off) then
            write(*,*) "bad values for set_mass_eqn_edge_dqs", &
               s% mass_eqn_lower_q_full_off, &
               s% mass_eqn_upper_q_full_off
            stop 'set_mass_eqn_edge_dqs'
         end if

         !remember that dq(k) = q(k)- q(k+1)

         !find upper edge
         !find index k for the cell just above the upper boundary that is rigid,
         !the first dq in the rigid region then corresponds to dq(k-1)
         if (s% q(s% nz) > s% mass_eqn_upper_q_full_off) then
            s% mass_eqn_upper_edge_dq = s% dq(s% nz)
            s% mass_eqn_upper_edge_k = s% nz
         else
            do k = 2, s% nz-1
               if (s% q(k+1) < s% mass_eqn_upper_q_full_off) then
                  s% mass_eqn_upper_edge_dq = s% dq(k-1)
                  s% mass_eqn_upper_edge_k = k
                  exit
               end if
            end do
         end if
         !compute average
         if (s% mass_eqn_upper_edge_k > 2) then
            num_avg = 0
            s% mass_eqn_upper_edge_dq = 0d0
            do k=s% mass_eqn_upper_edge_k, s% mass_eqn_upper_edge_k-s% mass_eqn_cells_for_edge_avg, -1
               if (k<2) exit
               num_avg = num_avg + 1
               s% mass_eqn_upper_edge_dq = s% mass_eqn_upper_edge_dq + s% dq(k-1)
            end do
            s% mass_eqn_upper_edge_dq = s% mass_eqn_upper_edge_dq/num_avg
         end if
         !write(*,*) "check edge", s% mass_eqn_upper_edge_k, &
         !   s%q(s% mass_eqn_upper_edge_k), s% q(s% mass_eqn_upper_edge_k+1), &
         !   s% mass_eqn_upper_q_full_off, num_avg, s% mass_eqn_upper_edge_dq

         !find lower edge
         !find index k for the cell just below the lower boundary that is rigid,
         !the first dq in the rigid region then corresponds to dq(k)
         if (s% q(1) < s% mass_eqn_lower_q_full_off) then
            s% mass_eqn_lower_edge_dq = s% dq(1)
            s% mass_eqn_lower_edge_k = 1
         else
            do k = s% nz, 2, -1
               if (s% q(k-1) > s% mass_eqn_lower_q_full_off) then
                  s% mass_eqn_lower_edge_dq = s% dq(k)
                  s% mass_eqn_lower_edge_k = k
                  exit
               end if
            end do
         end if
         !compute average
         if (s% mass_eqn_lower_edge_k <= s% nz) then
            num_avg = 0
            s% mass_eqn_lower_edge_dq = 0d0
            do k=s% mass_eqn_lower_edge_k, s% mass_eqn_lower_edge_k+s% mass_eqn_cells_for_edge_avg
               if (k>s% nz) exit
               num_avg = num_avg + 1
               s% mass_eqn_lower_edge_dq = s% mass_eqn_lower_edge_dq + s% dq(k)
            end do
            s% mass_eqn_lower_edge_dq = s% mass_eqn_lower_edge_dq/num_avg
         end if
         !write(*,*) "check edge", s% mass_eqn_lower_edge_k, &
         !   s%q(s% mass_eqn_lower_edge_k), s% q(s% mass_eqn_lower_edge_k-1), &
         !   s% mass_eqn_lower_q_full_off, num_avg, s% mass_eqn_lower_edge_dq

         !write(*,*) "check edge dqs", s% mass_eqn_lower_edge_dq, s% mass_eqn_upper_edge_dq, &
         !   s% mass_eqn_upper_edge_k, s% mass_eqn_lower_edge_k
      end subroutine set_mass_eqn_edge_dqs

      subroutine get_gr_gravity_factor(s, k, gr_factor, d_gr_factor_dlnR)
         ! note: this uses gravitational mass, m_grav
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: gr_factor, d_gr_factor_dlnR
         real(dp) :: twoGmrc2, invfac, dtwoGmrc2_dlnR, d_invfac_dlnR
         twoGmrc2 = 2*s% cgrav(k)*s% m_grav(k)/(s% r(k)*clight*clight)
         invfac = sqrt(1d0 - twoGmrc2)
         gr_factor = 1d0/invfac
         dtwoGmrc2_dlnR = -twoGmrc2
         d_invfac_dlnR = -dtwoGmrc2_dlnR/(2*invfac)
         d_gr_factor_dlnR = -d_invfac_dlnR/(invfac*invfac)
      end subroutine get_gr_gravity_factor


      real(dp) function interp_val_to_pt(v,k,sz,dq,str)
         use interp_1d_lib, only: interp_4_to_1
         integer, intent(in) :: k, sz
         real(dp), intent(in) :: v(:), dq(:)
         character (len=*), intent(in) :: str
         integer :: ierr
         include 'formats'
         if (k == 1) then
            interp_val_to_pt = v(k)
            return
         end if
         if (k > 2 .and. k < sz) then
            ierr = 0
            call interp_4_to_1( &
               0.5d0*(dq(k-2)+dq(k-1)), &
               0.5d0*(dq(k-1)+dq(k)), &
               0.5d0*(dq(k)+dq(k+1)), &
               0.5d0*dq(k-2)+dq(k-1), &
               v(k-2), v(k-1), v(k), v(k+1), &
               interp_val_to_pt, str, ierr)
            if (ierr == 0) return
            write(*,1) '0.5d0*(dq(k-2)+dq(k-1))', 0.5d0*(dq(k-2)+dq(k-1))
            write(*,1) '0.5d0*(dq(k-1)+dq(k))', 0.5d0*(dq(k-1)+dq(k))
            write(*,1) '0.5d0*(dq(k)+dq(k+1))', 0.5d0*(dq(k)+dq(k+1))
            write(*,2) 'dq(k-2)', k-2, dq(k-2)
            write(*,2) 'dq(k-1)', k-1, dq(k-1)
            write(*,2) 'dq(k)', k, dq(k)
            write(*,2) 'dq(k+1)', k+1, dq(k+1)

            stop 'interp_val_to_pt'
         endif
         interp_val_to_pt = (v(k)*dq(k-1) + v(k-1)*dq(k))/(dq(k-1) + dq(k))
      end function interp_val_to_pt


      real(dp) function interp_xa_to_pt(xa,j,k,sz,dq,str)
         use interp_1d_lib, only: interp_4_to_1
         real(dp), intent(in) :: xa(:,:), dq(:)
         character (len=*), intent(in) :: str
         integer, intent(in) :: j, k, sz
         integer :: ierr
         include 'formats'
         if (j == 0) then
            interp_xa_to_pt = 0
            return
         end if
         if (k == 1) then
            interp_xa_to_pt = xa(j,k)
            return
         end if
         if (k > 2 .and. k < sz) then
            ierr = 0
            call interp_4_to_1( &
               0.5d0*(dq(k-2)+dq(k-1)), &
               0.5d0*(dq(k-1)+dq(k)), &
               0.5d0*(dq(k)+dq(k+1)), &
               0.5d0*dq(k-2)+dq(k-1), &
               xa(j,k-2), xa(j,k-1), xa(j,k), xa(j,k+1), &
               interp_xa_to_pt, str, ierr)
            interp_xa_to_pt = min(1d0,max(0d0,interp_xa_to_pt))
            if (ierr == 0) return
         endif
         interp_xa_to_pt = (xa(j,k)*dq(k-1) + xa(j,k-1)*dq(k))/(dq(k-1) + dq(k))
         interp_xa_to_pt = min(1d0,max(0d0,interp_xa_to_pt))
      end function interp_xa_to_pt


      real(dp) function get_dtau1(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: k
         real(dp) :: kap, kap_min
         include 'formats'
         ierr = 0
         k = 1
         kap = s% opacity(1)
         if (s% tau_use_kap_floor) then
            if (s% kap_phot_factor_for_kap_floor >= 0d0) then
               kap_min = s% min_kap_floor*s% kap_phot_factor_for_kap_floor
            else
               if (s% Z(k) <= 0.02d0) then
                  kap_min = s% kap_min_Z_0pt02
               else if (s% Z(k) > 0.99d0) then
                  kap_min = s% kap_min_Z_1pt0
               else
                  kap_min = s% kap_min_Z_0pt02 + &
                     (s%Z (k) - 0.02d0)/(1d0 - 0.02d0)* &
                        (s% kap_min_Z_1pt0 - s% kap_min_Z_0pt02)
               end if
            end if
            kap = max(kap_min, kap)
         end if
         get_dtau1 = s% dm(1)*kap/(4*pi*s% rmid(1)*s% rmid(1))
         if (is_bad(get_dtau1)) then
            ierr = -1
            if (.not. s% report_ierr) return
            k = 1
            write(*,2) 'get_dtau1', k, get_dtau1
            write(*,2) 's% dm(1)', k, s% dm(k)
            write(*,2) 's% opacity(1)', k, s% opacity(k)
            write(*,2) 's% rmid(1)', k, s% rmid(k)
            write(*,2) 's% r(1)', k, s% r(k)
            write(*,2) 's% r(2)', 2, s% r(2)
            stop 'get_dtau1'
         end if
      end function get_dtau1


      subroutine get_tau(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         ! tau(k) is optical depth at outer boundary of cell k
         real(dp) :: dtau, dr, Z, kap_min, kap, dm_sum, L_sum
         integer :: k
         logical, parameter :: dbg = .true.
         include 'formats'
         ierr = 0
         dtau = get_dtau1(s, ierr)
         if (ierr /= 0) return
         s% tau(1) = s% tau_factor*s% tau_base
         s% lntau(1) = safe_log_cr(s% tau(1))
         !write(*,1) 'tau(1) from get_tau', tau(1), s% tau_factor, s% tau_base
         !write(*,1) 's% tau_base', s% tau_base
         !write(*,1) 's% tau_factor', s% tau_factor
         !if (tau(1) == 0d0) stop 'get_tau'
         s% tau_start(1) = s% tau(1)
         dm_sum = 0
         L_sum = 0
         do k = 2, s% nz
            s% tau(k) = s% tau(k-1) + dtau
            s% lntau(k) = log_cr(s% tau(k))
            if (s% tau_for_L_BB > 0d0 .and. s% L_for_BB_outer_BC <= 0d0) then
               dm_sum = dm_sum + s% dm(k-1)
               L_sum = L_sum + s% L(k-1)*s% dm(k-1)
               if (s% tau(k) >= s% tau_for_L_BB) then
                  s% L_for_BB_outer_BC = L_sum/dm_sum
                  !write(*,2) 's% L_for_BB_outer_BC', k, s% L_for_BB_outer_BC
               end if
            end if
            if (s% tau_start(k) < 0) s% tau_start(k) = s% tau(k)
            kap = s% opacity(k)
            dtau = s% dm(k)*kap/(4*pi*s% rmid(k)*s% rmid(k))
            if (is_bad(dtau)) then
               ierr = -1
               if (.not. s% report_ierr) return
               write(*,2) 'dtau', k, dtau
               write(*,2) 's% dm(k)', k, s% dm(k)
               write(*,2) 's% opacity(k)', k, s% opacity(k)
               write(*,2) 's% rmid(k)', k, s% rmid(k)
               stop 'get_tau'
            end if
            !write(*,*) 'dtau, dlogtau', k, tau(k) - tau(k-1), &
            !   log10(tau(k)/tau(k-1))
         end do
         if (s% tau_for_L_BB > 0d0 .and. s% L_for_BB_outer_BC < 0d0) then
            write(*,1) 'failed to set s% L_for_BB_outer_BC', s% L_for_BB_outer_BC
            stop 'get_tau'
         end if
      end subroutine get_tau


      integer function find_cell_for_mass(s, m)
         type (star_info), pointer :: s
         real(dp), intent(in) :: m
         integer :: k
         find_cell_for_mass = s% nz
         do k = 1, s% nz-1
            if (s% m(k) >= m .and. m > s% m(k+1)) then
               find_cell_for_mass = k
               return
            end if
         end do
      end function find_cell_for_mass


      subroutine get_delta_Pg(s, nu_max, delta_Pg)
         type (star_info), pointer :: s
         real(dp), intent(in) :: nu_max ! microHz
         real(dp), intent(out) :: delta_Pg ! seconds
         ! g-mode period spacing for l=1
         real(dp) :: integral, N2, omega2, kr2, L2, el, &
            dr, r, r2, cs2, sl2, I_integral, I_integral_limit
         integer :: k, k_sl2
         logical, parameter :: dbg = .false.
         include 'formats'
         if (dbg) then
            write(*,2) 'nu_max', s% model_number, nu_max
            write(*,2) 's% star_mass', s% model_number, s% star_mass
            write(*,2) 's% photosphere_r', s% model_number, s% photosphere_r
            write(*,2) 's% Teff', s% model_number, s% Teff
         end if
         delta_Pg = 0
         integral = 0
         I_integral = 0
         I_integral_limit = 0.5d0
         omega2 = pow2(2*pi*nu_max/1d6)
         if (dbg) write(*,1) 'log omega2', log10_cr(omega2)
         el = 1
         L2 = el*(el+1)
         k_sl2 = 0
         do k = 2, s% nz
            N2 = s% brunt_N2(k)
            r = s% r(k)
            r2 = r*r
            cs2 = s% csound_face(k)*s% csound_face(k)
            sl2 = L2*cs2/r2
            dr = s% rmid(k-1) - s% rmid(k)
            if (omega2 >= sl2) then
               cycle
            end if
            if (k_sl2 == 0) then
               k_sl2 = k
               if (dbg) write(*,2) 'k_sl2', k
            end if
            if (N2 > omega2) then ! in g-cavity
               if (dbg .and. integral == 0) write(*,2) 'enter g-cavity', k
               integral = integral + sqrt(N2)*dr/r
            else ! in decay region
               if (integral == 0) cycle ! ! haven't been in g-cavity yet
               if (dbg .and. I_integral == 0) write(*,2) 'enter decay', k
               ! in decay region below g-cavity; I_integral estimates decay
               kr2 = (1 - n2/omega2)*(1 - Sl2/omega2)*omega2/cs2
               I_integral = I_integral + sqrt(-kr2)*dr
               if (I_integral > I_integral_limit) exit
            end if
         end do

         if (dbg) write(*,2) 'omega2 nu_max integral I_integral', &
            s% model_number, omega2, nu_max, integral, I_integral

         if (integral == 0) return
         delta_Pg = sqrt(2d0)*pi*pi/integral
         if (is_bad(delta_Pg)) delta_Pg = 0

         if (dbg) write(*,2) 'delta_Pg', s% model_number, delta_Pg

      end subroutine get_delta_Pg


      subroutine set_rmid(s, nzlo, nzhi, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr
         integer :: k, nz
         real(dp) :: r003, rp13, rmid, rmid2
         logical :: dbg
         include 'formats'
         ierr = 0
         nz = s% nz
         dbg = .false.
         if (s% RSP_flag) then
            !dbg = s% model_number >= s% max_model_number - 1
            do k=nzlo, nzhi
               if (k < nz) then
                  rmid = 0.5d0*(s% r(k) + s% r(k+1))
               else
                  rmid = 0.5d0*(s% r(k) + s% R_center)
               end if
               s% rmid(k) = rmid
               if (dbg) write(*,3) 'set_rmid s% r(k)', k, s% model_number, s% r(k)
               if (dbg) write(*,3) 'set_rmid s% rmid(k)', k, s% model_number, s% rmid(k)
               if (s% rmid_start(k) < 0) s% rmid_start(k) = s% rmid(k)
               rmid2 = rmid*rmid
               s% drmid_dlnR00(k) = 0.5d0*s% r(k)
               s% drmid2_dlnR00(k) = 2d0*rmid*s% drmid_dlnR00(k)
               if (k < nz) then
                  s% drmid_dlnRp1(k) = 0.5d0*s% r(k+1)
                  s% drmid2_dlnRp1(k) = 2d0*rmid*s% drmid_dlnRp1(k)
               else
                  s% drmid_dlnRp1(k) = 0d0
                  s% drmid2_dlnRp1(k) = 0d0
               end if
            end do
            return
         end if
         do k=nzlo, nzhi
            r003 = s% r(k)*s% r(k)*s% r(k)
            if (k < nz) then
               rp13 = s% r(k+1)*s% r(k+1)*s% r(k+1)
            else
               rp13 = s% R_center*s% R_center*s% R_center
            end if
            rmid = pow_cr(0.5d0*(r003 + rp13),1d0/3d0)
            s% rmid(k) = rmid
            if (s% rmid_start(k) < 0) s% rmid_start(k) = s% rmid(k)
            rmid2 = rmid*rmid
            s% drmid_dlnR00(k) = 0.5d0*r003/rmid2
            s% drmid2_dlnR00(k) = r003/rmid
            if (k < nz) then
               s% drmid_dlnRp1(k) = 0.5d0*rp13/rmid2
               s% drmid2_dlnRp1(k) = rp13/rmid
            else
               s% drmid_dlnRp1(k) = 0d0
               s% drmid2_dlnRp1(k) = 0d0
            end if
         end do
      end subroutine set_rmid


      real(dp) function get_tau_at_r(s, r, ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: r
         integer, intent(out) :: ierr
         real(dp) :: dtau, dr, tau_m1, tau_00
         integer :: k
         logical, parameter :: dbg = .false.
         include 'formats'
         ierr = 0
         dtau = get_dtau1(s, ierr)
         if (ierr /= 0) return
         tau_00 = s% tau_factor*s% tau_base
         get_tau_at_r = tau_00
         if (r >= s% r(1)) return
         do k = 2, s% nz
            tau_m1 = tau_00
            tau_00 = tau_m1 + dtau
            if (r < s% r(k-1) .and. r >= s% r(k)) then
               get_tau_at_r = &
                  (tau_00*(s% r(k-1)-r) + tau_m1*(r-s% r(k)))/(s% r(k-1)-s% r(k))
               return
            end if
            dtau = s% dm(k)*s% opacity(k)/(4*pi*s% rmid(k)*s% rmid(k))
         end do
      end function get_tau_at_r


      integer function find_tau_phot(s, tau00, taup1, ierr)
         ! return k for the cell containing optical depth = tau_base
         type (star_info), pointer :: s
         real(dp), intent(out) :: tau00, taup1
         integer, intent(out) :: ierr
         integer :: k
         real(dp) :: dtau, tau_phot

         include 'formats'
         ierr = 0
         tau00 = 0
         taup1 = 0
         find_tau_phot = 1
         if (s% tau_factor >= 1) return
         tau_phot = s% tau_base
         tau00 = s% tau_factor*tau_phot
         do k = 1, s% nz
            dtau = s% dm(k)*s% opacity(k)/(4*pi*s% rmid(k)*s% rmid(k))
            taup1 = tau00 + dtau
            if (taup1 >= tau_phot) then
               find_tau_phot = k
               return
            end if
            tau00 = taup1
         end do
         ierr = -1
      end function find_tau_phot
      

      
      
      real(dp) function get_r_phot(s)
         type (star_info), pointer :: s  
         real(dp) :: r, m, v, L, T_phot, cs, kap, logg, ysum
         integer :: k_phot
         call get_phot_info(s,r,m,v,L,T_phot,cs,kap,logg,ysum,k_phot)
         get_r_phot = r
      end function get_r_phot


      subroutine get_phot_info(s,r,m,v,L,T_phot,cs,kap,logg,ysum,k_phot)
         type (star_info), pointer :: s
         real(dp), intent(out) :: r, m, v, L, T_phot, cs, kap, logg, ysum
         integer, intent(out) :: k_phot

         integer :: k
         real(dp) :: tau00, taup1, dtau, r003, rp13, r3, tau_phot, &
            Tface_0, Tface_1

         include 'formats'

         tau00 = 0
         taup1 = 0
         ysum = 0
         r = s% r(1)
         m = s% m(1)
         if (s% u_flag) then
            v = s% u(1)
         else if (s% v_flag) then
            v = s% v(1)
         else
            v = 0d0
         end if
         L = max(1d0, s% L(1)) ! don't use negative L(1)
         T_phot = s% T(1)
         cs = s% csound(1)
         kap = s% opacity(1)
         logg = log10_cr(s% cgrav(1)*m/(r*r))
         k_phot = 1
         tau_phot = s% tau_base
         tau00 = s% tau_factor*tau_phot
         if (tau00 >= tau_phot) return
         do k = 1, s% nz-1
            dtau = s% dm(k)*s% opacity(k)/(4*pi*s% rmid(k)*s% rmid(k))
            taup1 = tau00 + dtau
            ysum = ysum + s% rho(k)*(s% r(k) - s% r(k+1))
            if (taup1 >= tau_phot .and. dtau > 0d0) then
               if (k == 1) then
                  Tface_0 = s% T(k)
               else
                  Tface_0 = 0.5d0*(s% T(k) + s% T(k-1))
               end if
               Tface_1 = 0.5d0*(s% T(k) + s% T(k+1))
               T_phot = Tface_0 + (Tface_1 - Tface_0)*(tau_phot - tau00)/dtau
               r003 = s% r(k)*s% r(k)*s% r(k)
               rp13 = s% r(k+1)*s% r(k+1)*s% r(k+1)
               r3 = r003 + (rp13 - r003)*(tau_phot - tau00)/dtau
               r = pow_cr(r3,1d0/3d0)
               m = s% m(k) - s% dm(k)*(tau_phot - tau00)/dtau
               if (s% u_flag) then
                  v = s% u_face(k) + &
                     (s% u_face(k+1) - s% u_face(k))*(tau_phot - tau00)/dtau
               else if (s% v_flag) then
                  v = s% v(k) + (s% v(k+1) - s% v(k))*(tau_phot - tau00)/dtau
               end if
               L = s% L(k) + (s% L(k+1) - s% L(k))*(tau_phot - tau00)/dtau
               k_phot = k
               cs = s% csound(k_phot)
               kap = s% opacity(k_phot)
               logg = log10_cr(s% cgrav(k_phot)*m/(r*r))
               return
            end if
            tau00 = taup1
         end do
         !write(*,*) 'get_phot_info failed to find photosphere'
         k_phot = s% nz
         r = s% R_center
         m = s% m_center
         v = s% v_center
         T_phot = s% T(k_phot)
         L = max(1d0, s% L_center)
         cs = s% csound(k_phot)
         kap = s% opacity(k_phot)
         logg = log10_cr(s% cgrav(k_phot)*m/(r*r))
      end subroutine get_phot_info


      real(dp) function center_value(s, p)
         type (star_info), pointer :: s
         real(dp), intent(in) :: p(:)
         real(dp) :: sum_x, sum_dq, dx, dq
         integer :: k
         sum_x = 0
         sum_dq = 0
         do k = s% nz, 1, -1
            dq = s% dq(k)
            dx = p(k)*dq
            if (sum_dq+dq >= s% center_avg_value_dq) then
               sum_x = sum_x + dx*(s% center_avg_value_dq - sum_dq)/dq
               sum_dq = s% center_avg_value_dq
               exit
            end if
            sum_x = sum_x + dx
            sum_dq = sum_dq + dq
         end do
         center_value = sum_x/sum_dq
      end function center_value


      subroutine interp_q( &
            nz2, nvar_hydro, species, qval, xh, xa, q, dq, struct, comp, ierr)
         use num_lib, only: binary_search
         integer, intent(in) :: nz2, nvar_hydro, species
         real(dp), intent(in) :: qval
         real(dp), intent(in) :: xh(:,:), xa(:,:), q(:), dq(:)
         real(dp), intent(inout) :: struct(:), comp(:)
         integer, intent(out) :: ierr
         integer :: k
         real(dp) :: alfa
         ierr = 0
         if (qval <= q(nz2)) then
            if (nvar_hydro > 0) &
               struct(1:nvar_hydro) = xh(1:nvar_hydro,nz2)
            if (species > 0) &
               comp(1:species) = xa(1:species,nz2)
            return
         end if
         k = binary_search(nz2, q, 0, qval)
         if (k < 1 .or. k >= nz2) then
            ierr = -1
            return
         end if
         if (qval <= q(k) .and. qval > q(k+1)) then
            alfa = (qval - q(k+1)) / dq(k)
            if (nvar_hydro > 0) &
               struct(1:nvar_hydro) = &
                  alfa*xh(1:nvar_hydro,k) + (1-alfa)*xh(1:nvar_hydro,k+1)
            if (species > 0) &
               comp(1:species) = alfa*xa(1:species,k) + (1-alfa)*xa(1:species,k+1)
            return
         end if
         ierr = -1
      end subroutine interp_q


      subroutine std_write_internals_to_file(id, num)
         use utils_lib, only : mkdir
         integer, intent(in) :: num, id
         character (len=strlen) :: fname
         integer :: ierr
         ierr = 0
         call mkdir('plot_data')
         write(fname, '(a, i1, a)') 'plot_data/internals', mod(abs(num), 10), '.data'
         write(*,*) 'call write_internals_to_file ' // trim(fname)
         call write_internals_to_file(id, fname, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in write_internals_to_file ' // trim(fname)
         end if
      end subroutine std_write_internals_to_file


      subroutine write_internals_to_file(id, filename, ierr)
         use utils_lib
         character (len=*), intent(in) :: filename
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         integer :: iounit
         ierr = 0
         iounit = alloc_iounit(ierr); if (ierr /= 0) return
         open(iounit, file=trim(filename), action='write', status='replace', iostat=ierr)
         if (ierr == 0) then
            call write_internals(id, iounit, ierr)
            close(iounit)
         else
            write(*, *) 'failed to open internals file ' // trim(filename)
         end if
         call free_iounit(iounit)
      end subroutine write_internals_to_file


      subroutine write_internals(id, iounit, ierr)
         use chem_def
         integer, intent(in) :: iounit, id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'write_internals: get_star_ptr ierr', ierr
            return
         end if
         call write_model_info(s, iounit, ierr)
      end subroutine write_internals


      subroutine write_model_info(s, iounit, ierr)
         use chem_def
         type (star_info), pointer :: s
         integer, intent(in) :: iounit
         integer, intent(out) :: ierr
         integer, pointer :: chem_id(:)
         integer :: k, i, nz, species
         integer :: he4

         ierr = 0
         
         nz = s% nz
         species = s% species
         chem_id => s% chem_id
         he4 = s% net_iso(ihe4)

         write(iounit,'(a)') '            mass         initial_z       n_shells'
         write(iounit,'(2x,2f15.4,i15)') s% star_mass, s% initial_z, nz
         write(iounit,fmt='(i5)',advance='no') 1
         do i=2,88
            write(iounit,fmt='(i12,15x)',advance='no') i
         end do
         write(iounit,*)
         write(iounit,fmt='(a5,1x,99(a26,1x))',advance='no') &
               'grid', 'r', 'm', 'log_dq', &
               'log10d', 'log10T', 'log10m', 'log10r', 'L', 'r_div_rstar', &
               'log10P', 'log10Pgas', 'chiT', 'chiRho', &
               'dlnRho_dlnPgas_const_T', 'dlnRho_dlnT_const_Pgas', &
               'extra5', 'extra6', 'extra7', 'extra8'
         ! composition info
         do i=1, species
            write(iounit, fmt='(a26, 1x)', advance='no') trim(chem_isos% name(chem_id(i)))
         end do
         do i=1, species
            write(iounit, fmt='(a26, 1x)', advance='no') 'lg_' // trim(chem_isos% name(chem_id(i)))
         end do
         write(iounit,fmt=*)

         do k=1, nz
            write(iounit,'(i5,1x,99(1pe26.16,1x))',advance='no') k,  &
               s% r(k)/Rsun, s% m(k)/Msun, safe_log10_cr(s% dq(k)), &
               s% lnd(k)/ln10, s% lnT(k)/ln10, log10_cr(s% m(k)),  &
               s% lnR(k)/ln10, s% L(k)/Lsun, s% r(k)/s% r(1), &
               s% lnP(k)/ln10, s% lnPgas(k)/ln10, s% chiT(k), s% chiRho(k), &
               s% dlnRho_dlnPgas_const_T(k), s% dlnRho_dlnT_const_Pgas(k), &
               s% profile_extra(k,5), s% profile_extra(k,6), &
               s% profile_extra(k,7), s% profile_extra(k,8)
            do i=1, species
               write(iounit, fmt='(1pe26.16, 1x)', advance='no') s% xa(i, k)
            end do
            do i=1, species
               write(iounit, fmt='(1pe26.16, 1x)', advance='no') safe_log10_cr(s% xa(i, k))
            end do
            write(iounit,*)
         end do

      end subroutine write_model_info


      subroutine dump_model_info_for_ndiff(s, iounit, ierr)
         use chem_def
         type (star_info), pointer :: s
         integer, intent(in) :: iounit
         integer, intent(out) :: ierr
         integer, pointer :: chem_id(:)
         integer :: k, j, nz, species
         include 'formats'
         ierr = 0
         nz = s% nz
         species = s% species
         chem_id => s% chem_id
         write(iounit,*) 'nz', nz
         write(iounit,1) 'star_mass', s% star_mass
         write(iounit,1) 'initial_z', s% initial_z
         do k=1, nz
            do j=1, s% nvar_hydro
               write(iounit,2) trim(s% nameofvar(j)), k, s% xh(j,k)
            end do
            do j=1,species
               write(iounit,2) trim(chem_isos% name(chem_id(j))), k, clip(s% xa(j,k))
            end do
         end do

         contains

         real(dp) function clip(x)
            real(dp), intent(in) :: x
            if (.true. .or. x > 1d-30) then
               clip = x
            else
               clip = 0d0
            end if
         end function clip

      end subroutine dump_model_info_for_ndiff

      
      subroutine set_abs_du_div_cs(s)
         type (star_info), pointer :: s
         
         integer :: k, nz, j
         real(dp) :: abs_du, cs
         include 'formats'
         nz = s% nz

         if (s% u_flag) then
            do k=2,nz-1
               abs_du = &
                  max(abs(s% u_start(k) - s% u_start(k+1)), &
                      abs(s% u_start(k) - s% u_start(k-1)))
               cs = maxval(s% csound(max(1,k-5):min(nz,k+5)))
               s% abs_du_plus_cs(k) = abs_du + cs
               s% abs_du_div_cs(k) = abs_du/cs
            end do
            k = nz
            s% abs_du_plus_cs(k) = &
               abs(s% u_start(k) - s% u_start(k-1)) + s% csound_start(k)
            s% abs_du_div_cs(k) = &
               abs(s% u_start(k) - s% u_start(k-1))/s% csound_start(k)
            k = 2
            s% abs_du_plus_cs(k) = &
               abs(s% u_start(k) - s% u_start(k+1)) + s% csound_start(k)
            s% abs_du_div_cs(k) = &
               abs(s% u_start(k) - s% u_start(k+1))/s% csound_start(k)
            k = 1
            s% abs_du_plus_cs(k) = s% abs_du_plus_cs(k+1)
            s% abs_du_div_cs(k) = s% abs_du_div_cs(k+1)
            do j = 1,3
               do k=2,nz-1
                  s% abs_du_div_cs(k) = sum(s% abs_du_div_cs(k-1:k+1))/3d0
               end do
            end do
         else
            do k=1,nz
               s% abs_du_plus_cs(k) = 1d99
               s% abs_du_div_cs(k) = 1d99
            end do
         end if
      
      end subroutine set_abs_du_div_cs


      real(dp) function rsi_div_rsimelt(s,k,species) 
         ! rsi = ion density parameter for cell k
         ! rsimelt = ion density parameter of quantum melting
         ! rsi < rsimelt => liquid, independent of T
         use chem_def, only: chem_isos
         type (star_info), pointer :: s
         integer, intent(in) :: k, species
         
         integer :: IX, j
         real(dp), dimension(species) :: AZion, ACMI, AY
         real(dp) :: Y, CMImean, Z73, RS, RSI
         real(dp), parameter :: RSIMELT=140d0, TINY=1d-7, &
            AUM=1822.888d0 ! a.m.u./m_e
         
         include 'formats'
         
         ! details from eos/private/pc_eos.f
         
         AZion(1:species) = chem_isos% Z(s% chem_id(1:species))
         ACMI(1:species) = chem_isos% W(s% chem_id(1:species))
         do j=1,species
            if (s% xa(j,k) < s% job% mass_fraction_limit_for_PC) then
               AY(j) = 0
            else
               AY(j) = s% xa(j,k)/ACMI(j)
            end if
         end do
         
         Y=0.d0
         do IX=1,species
            Y=Y+AY(IX)
         end do
         if (dabs(Y-1.d0).gt.TINY) then
           do IX=1,species
              AY(IX)=AY(IX)/Y
           end do
         end if

         CMImean=0.d0
         Z73=0.d0
         do IX=1,species
            if (AY(IX) < TINY) cycle
            Z73 = Z73 + AY(IX)*pow_cr(AZion(IX),7d0/3d0)
            CMImean = CMImean + AY(IX)*ACMI(IX)
         end do

         RS=pow_cr(0.75d0/PI/s% rho(k),1d0/3d0)
         RSI=RS*CMImean*Z73*AUM
         
         if (is_bad(RSI)) then
            write(*,2) 'RSI', k, RSI
            write(*,2) 'Z73', k, Z73
            write(*,2) 'CMImean', k, CMImean
            write(*,2) 'RS', k, RS
            write(*,2) 's% rho(k)', k, s% rho(k)
            !write(*,2) '', k, 
            !write(*,2) '', k, 
            stop 'rsi_div_rsimelt'
         end if
         
         rsi_div_rsimelt = RSI/RSIMELT
         
      end function rsi_div_rsimelt


      subroutine get_shock_info(s)
         type (star_info), pointer :: s
         integer :: k, nz, kk, kmin, k0, k1, ierr
         real(dp) :: v_div_cs_00, v_div_cs_m1, rmax, vmax, alfa
         real(dp), pointer :: v(:)

         include 'formats'
         
         if (s% u_flag) then
            v => s% u ! may not have u_face
         else if (s% v_flag) then
            v => s% v
         else
            return
         end if

         nz = s% nz
         kmin = nz
         do k=1,nz-1
            if (s% q(k) <= s% max_q_for_outer_mach1_location) then
               kmin = k
               exit
            end if
         end do

         s% shock_velocity = 0
         s% shock_csound = 0
         s% shock_lgT = 0
         s% shock_lgRho = 0
         s% shock_lgP = 0
         s% shock_mass = 0
         s% shock_q = 0
         s% shock_radius = 0
         s% shock_gamma1 = 0
         s% shock_entropy = 0
         s% shock_tau = 0
         s% shock_k = 0
         s% shock_pre_lgRho = 0
         
         if (kmin < nz) then ! search inward for 1st shock moving outward
            v_div_cs_00 = v(kmin)/s% csound_face(kmin)
            do k = kmin+1,nz
               v_div_cs_m1 = v_div_cs_00
               v_div_cs_00 = v(k)/s% csound_face(k)
               if (v_div_cs_00 > 0 .and. v_div_cs_m1 > 0 .and. &
                     v_div_cs_00 >= 1d0 .and. v_div_cs_m1 < 1d0) then
                  do kk = k+1, nz ! search inward for local max v
                     if (v(kk) <= v(kk-1)) then
                        s% shock_k = kk - 1
                        exit
                     end if
                  end do
                  exit
               end if
            end do
         end if

         k = s% shock_k
         if (k < nz .and. k > 1) then
            if (v(k) > max(v(k-1),v(k+1))) then
               rmax = s% r(k)
               vmax = v(k)
               s% shock_velocity = vmax
               s% shock_radius = rmax/Rsun
               if (rmax <= s% r(k)) then
                  k0 = k-1
                  k1 = k
               else if (rmax <= s% r(k+1)) then
                  k0 = k
                  k1 = k+1
               else
                  k0 = k+1
                  k1 = k+2
               end if
               alfa = (s% r(k0) - rmax)/(s% r(k0) - s% r(k1))
               s% shock_q = s% q(k0) - alfa*s% dq(k0)
               s% shock_mass = (s% m(k0) - alfa*s% dm(k0))/Msun
               s% shock_csound = s% csound(k0)
               s% shock_lgT = s% lnT(k0)/ln10
               s% shock_lgRho = s% lnd(k0)/ln10
               s% shock_lgP = s% lnP(k0)/ln10
               s% shock_gamma1 = s% gamma1(k0)
               s% shock_entropy = s% entropy(k0)
               s% shock_tau = s% tau(k0)
               s% shock_k = k0
               do kk=k0-1,1,-1
                  if (v(kk) < 0.1d0*s% csound_face(kk)) then
                     s% shock_pre_lgRho = s% lnd(kk)/ln10
                     exit
                  end if
               end do
            end if
         end if

      end subroutine get_shock_info


      real(dp) function min_dr_div_cs(s,min_k) ! seconds
         use crlibm_lib, only: log10_cr
         type (star_info), pointer :: s
         integer, intent(out) :: min_k
         integer :: k, nz, j, k_min
         real(dp) :: dr, dt, D, abs_du, cs, min_q, max_q, &
            min_abs_du_div_cs, r00, rp1, dr_div_cs
         include 'formats'
         nz = s% nz
         min_k = nz
         min_dr_div_cs = 1d99
         if (s% v_flag) then
            r00 = s% R_center
            do k=nz,1,-1
               rp1 = r00
               r00 = s% r(k)
               dr_div_cs = (r00 - rp1)/s% csound(k)
               if (dr_div_cs < min_dr_div_cs) then
                  min_dr_div_cs = dr_div_cs
                  min_k = k
               end if
            end do
            return
         end if
         if (.not. s% u_flag) return
         min_abs_du_div_cs = &
            s% min_abs_du_div_cs_for_dt_div_min_dr_div_cs_limit
         min_q = s% min_q_for_dt_div_min_dr_div_cs_limit
         max_q = s% max_q_for_dt_div_min_dr_div_cs_limit
         k_min = max(1, s% min_k_for_dt_div_min_dr_div_cs_limit)
         do k = k_min, nz
            if (s% q(k) > max_q) cycle
            if (s% q(k) < min_q) exit
            if (s% abs_du_div_cs(k) < min_abs_du_div_cs) cycle
            dr = s% r(k) - s% r(k+1)
            dt = dr/s% abs_du_plus_cs(k)
            if (dt < min_dr_div_cs) then
               min_dr_div_cs = dt
               min_k = k
            end if
         end do
      end function min_dr_div_cs


      subroutine reset_starting_vectors(s)
         type (star_info), pointer :: s
         integer :: k, nz
         nz = s% nz
         do k=1,s% nz
            s% T_start(k) = -1d99
            s% r_start(k) = -1d99
            s% rmid_start(k) = -1d99
            s% v_start(k) = -1d99
            s% u_start(k) = -1d99
            s% lnd_start(k) = -1d99
            s% lnT_start(k) = -1d99
            s% csound_start(k) = -1d99
            s% eta_visc_start(k) = -1d99
            s% rho_start(k) = -1d99
            s% tau_start(k) = -1d99
            s% erad_start(k) = -1d99
            s% alpha_RTI_start(k) = -1d99
            s% w_start(k) = -1d99
            s% dPdr_dRhodr_info(k) = -1d99
            if (.not. s% op_split_burn) cycle
            s% eps_nuc(k) = 0d0
            s% d_epsnuc_dlnd(k) = 0d0
            s% d_epsnuc_dlnT(k) = 0d0
            s% d_epsnuc_dx(:,k) = 0d0
            s% eps_nuc_categories(:,k) = 0d0
            s% dxdt_nuc(:,k) =  0d0
            s% d_dxdt_nuc_dRho(:,k) =  0d0
            s% d_dxdt_nuc_dT(:,k) =  0d0
            s% d_dxdt_nuc_dx(:,:,k) =  0d0
            s% eps_nuc_neu_total(k) = 0d0
         end do
      end subroutine reset_starting_vectors


      ! largest k s.t. for all k' < k, cell k' has Cp(k')*T(k')*mstar_dot < L(k).
      subroutine set_k_CpTMdot_lt_L(s)
         type (star_info), pointer :: s
         integer :: k, nz
         if (s% mstar_dot <= 0d0 .or. s% gamma_law_hydro > 0d0) then
            s% k_CpTMdot_lt_L = 1
            return
         end if
         nz = s% nz
         do k = 2, nz
            if (s% Cp(k)*s% T(k)*s% mstar_dot >= max(1d-99,s% L(k))) then
               s% k_CpTMdot_lt_L = k-1
               return
            end if
         end do
         s% k_CpTMdot_lt_L = nz
      end subroutine set_k_CpTMdot_lt_L


      subroutine set_scale_height(s)
         type (star_info), pointer :: s
         real(dp) :: Hp, alt_Hp, alfa, beta, rho_face, P_face
         integer :: k
         include 'formats'
         !if (s% gamma_law_hydro > 0d0) return
         do k=1,s% nz
            if (s% cgrav(k) == 0) then
               s% scale_height(k) = 0
               cycle
            end if
            if (k == 1) then
               alfa = 1
            else
               alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            end if
            beta = 1 - alfa
            if (alfa == 1) then
               rho_face = s% rho(k)
               P_face = s% P(k)
            else
               rho_face = alfa*s% rho(k) + beta*s% rho(k-1)
               P_face = alfa*s% P(k) + beta*s% P(k-1)
            end if
            Hp = P_face/(rho_face*s% grav(k))
            if (s% cgrav(k) <= 0) then
               alt_Hp = s% r(k)
            else
               alt_Hp = sqrt(P_face / s% cgrav(k)) / rho_face
            end if
            s% scale_height(k) = min(Hp, alt_Hp)
         end do
      end subroutine set_scale_height


      real(dp) function tau_eff(s,k)
         ! tau_eff = tau that gives the local P == P_atm if this location at surface
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: P, g, Pextra_factor
         if (k == 1) then
            tau_eff = s% tau(1)
            return
         end if
         if (s% cgrav(k) <= 0d0) then
            tau_eff = 0d0
            return
         end if
         P = (s% dq(k-1)*s% P(k) + s% dq(k)*s% P(k-1))/(s% dq(k-1) + s% dq(k))
         g = s% cgrav(k)*s% m_grav(k)/(s% r(k)*s% r(k))
         if (s% Pextra_factor < 0) then ! old form
            Pextra_factor = 1.6d-4
         else
            Pextra_factor = s% Pextra_factor
         end if
         tau_eff = s% opacity(k)*(P/g - &
               Pextra_factor*(s% L(k)/s% m_grav(k))/(6d0*pi*clight*s% cgrav(k)))
      end function tau_eff


      real(dp) function eval_Ledd(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         real(dp) :: dtau1, dtau, dr, tau, dqsum, Ledd_sum
         integer :: k
         logical, parameter :: dbg = .false.
         include 'formats'
         ierr = 0
         eval_Ledd = 0d0
         if (s% cgrav(1) <= 0d0) return
         dtau1 = get_dtau1(s, ierr)
         if (ierr /= 0) return
         dtau = dtau1
         tau = s% tau_factor*s% tau_base
         dqsum = s% dq(1)
         Ledd_sum = s% dq(1)*4*pi*clight*s% cgrav(1)*s% m_grav(1)/s% opacity(1)
         do k = 2, s% nz
            tau = tau + dtau
            if (tau > s% surf_avg_tau) exit
            dtau = s% dm(k)*s% opacity(k)/(4*pi*s% rmid(k)*s% rmid(k))
            dqsum = dqsum + s% dq(k)
            Ledd_sum = Ledd_sum + &
               s% dq(k)*4*pi*clight*s% cgrav(1)*s% m_grav(1)/s% opacity(k)
         end do
         eval_Ledd = Ledd_sum/dqsum
      end function eval_Ledd


      real(dp) function eval_min_cell_collapse_time(s,k_lo,k_hi,min_collapse_k,ierr) &
            result(min_collapse_time)
         type (star_info), pointer :: s
         integer, intent(in) :: k_lo, k_hi
         integer, intent(out) :: min_collapse_k, ierr
         real(dp) :: rp1, vp1, r00, v00, time
         integer :: k
         logical, parameter :: dbg = .false.
         include 'formats'
         ierr = 0
         rp1 = s% R_center
         vp1 = s% v_center
         min_collapse_time = 1d99
         min_collapse_k = -1
         if (.not. s% v_flag) return
         do k = k_hi, k_lo, -1
            v00 = s% v(k)
            r00 = s% r(k)
            if (r00 <= rp1) then
               !ierr = -1
               min_collapse_time = -1
               min_collapse_k = -1
               return
               write(*,2) 'bad radii', k, r00, rp1
               stop 'eval_min_cell_collapse_time'
            end if
            if (vp1 > v00) then
               time = (r00 - rp1)/(vp1 - v00)
               if (time < min_collapse_time) then
                  min_collapse_time = time
                  min_collapse_k = k
               end if
            end if
            rp1 = r00
            vp1 = v00
         end do
      end function eval_min_cell_collapse_time


      real(dp) function total_angular_momentum(s) result(J)
         type (star_info), pointer :: s
         include 'formats'
         if (.not. s% rotation_flag) then
            J = 0
         else
            J = dot_product(s% dm_bar(1:s% nz), s% j_rot(1:s% nz))
         end if
      end function total_angular_momentum


      real(dp) function eval_irradiation_heat(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: irradiation_dq, xq, eps
         eval_irradiation_heat = 0
         if (s% irradiation_flux /= 0) then
            irradiation_dq = 4*pi*s% r(1)*s% r(1)*s% column_depth_for_irradiation/s% xmstar
            xq = 1 - s% q(k)
            if (irradiation_dq > xq) then ! add irradiation heat for cell k
               eps = 0.25d0 * s% irradiation_flux / s% column_depth_for_irradiation
               if (irradiation_dq < xq + s% dq(k)) then ! only part of cell gets heated
                  eval_irradiation_heat = eps*(irradiation_dq - xq)/s% dq(k)
               else ! all of cell gets heated
                  eval_irradiation_heat = eps
               end if
            end if
         end if
      end function eval_irradiation_heat


      subroutine start_time(s, time0, total_all_before)
         type (star_info), pointer :: s
         integer, intent(out) :: time0
         real(dp), intent(out) :: total_all_before
         integer :: clock_rate
         if (.not. s% doing_timing) return
         total_all_before = total_times(s)
         call system_clock(time0,clock_rate)
      end subroutine start_time


      subroutine update_time(s, time0, total_all_before, total)
         type (star_info), pointer :: s
         integer, intent(in) :: time0
         real(dp), intent(in) :: total_all_before
         real(dp), intent(inout) :: total
         real(dp) :: total_all_after, other_stuff
         integer :: time1, clock_rate
         if (.not. s% doing_timing) return
         call system_clock(time1,clock_rate)
         total_all_after = total_times(s)
         other_stuff = total_all_after - total_all_before
            ! don't double count any other stuff
         total = total + (dble(time1-time0)/clock_rate - other_stuff)
      end subroutine update_time


      real(dp) function total_times(s)
         type (star_info), pointer :: s
         total_times = &
            s% time_evolve_step + &
            s% time_remesh + &
            s% time_adjust_mass + &
            s% time_element_diffusion + &
            s% time_struct_burn_mix + &
            s% time_newton_matrix + &
            s% time_solve_mix + &
            s% time_solve_burn + &
            s% time_solve_omega_mix + &
            s% time_eos + &
            s% time_neu_kap + &
            s% time_nonburn_net + &
            s% time_mlt + &
            s% time_set_hydro_vars
      end function total_times



      subroutine smooth(dc, sz)
         real(dp), intent(inout) :: dc(:)
         integer, intent(in) :: sz
         integer :: k
         k = 1
         dc(k) = (3*dc(k) + dc(k+1))/4
         do k=2,sz-1
            dc(k) = (dc(k-1) + 2*dc(k) + dc(k+1))/4
         end do
         k = sz
         dc(k) = (dc(k-1) + 3*dc(k))/4
      end subroutine smooth


      subroutine do_boxcar_mixing( &
            s, min_mass, max_mass, boxcar_nominal_mass, number_iterations, ierr)
         ! code based on SNEC
         type (star_info), pointer :: s
         real(dp), intent(in) :: min_mass, max_mass, boxcar_nominal_mass
         integer, intent(in) :: number_iterations
         integer, intent(out) :: ierr
         integer :: nz, iter, k_inner, k_outer, j, k
         real(dp) :: dm_sum, target_mass, actual_mass, xa, min_m, max_m
         include 'formats'
         ierr = 0
         nz = s% nz
         target_mass = boxcar_nominal_mass*Msun
         min_m = min_mass*Msun
         max_m = max_mass*Msun
         write(*,2) 'iters min max box', number_iterations, min_mass, max_mass, boxcar_nominal_mass
         do iter = 1, number_iterations
            do k_inner = nz, 1, -1
               if (s% m(k_inner) < min_m) cycle
               dm_sum = 0
               k_outer = 1
               do k = k_inner, 1, -1
                  dm_sum = dm_sum + s% dm(k)
                  if (dm_sum > target_mass) then
                     k_outer = k
                     exit
                  end if
               end do
               if (s% m(k_outer) > max_m) cycle
               actual_mass = dm_sum
               do j=1,s% species
                  xa = 0d0
                  do k=k_outer,k_inner
                     xa = xa + s% xa(j,k)*s% dm(k)
                  end do
                  do k=k_outer,k_inner
                     s% xa(j,k) = xa/dm_sum
                  end do
               end do
               if (actual_mass < target_mass) exit
            end do
         end do
      end subroutine do_boxcar_mixing


      subroutine get_XYZ(s, xa, X, Y, Z)
         use chem_def, only: ih1, ih2, ihe3, ihe4
         type (star_info), pointer :: s
         real(dp), intent(in) :: xa(:)
         real(dp), intent(out) :: X, Y, Z
         X = 0d0
         if (s% net_iso(ih1) /= 0) X = X + xa(s% net_iso(ih1))
         if (s% net_iso(ih2) /= 0) X = X + xa(s% net_iso(ih2))
         X = min(1d0, max(0d0, X))
         Y = 0d0
         if (s% net_iso(ihe3) /= 0) Y = Y + xa(s% net_iso(ihe3))
         if (s% net_iso(ihe4) /= 0) Y = Y + xa(s% net_iso(ihe4))
         Y = min(1d0, max(0d0, Y))
         Z = min(1d0, max(0d0, 1d0 - (X + Y)))
      end subroutine get_XYZ


      subroutine get_face_values(s, v_mid, v_face, ierr)
         ! simple interpolation by mass
         type (star_info), pointer :: s
         real(dp), intent(in) :: v_mid(:)
         real(dp), intent(inout) :: v_face(:)
         integer, intent(out) :: ierr
         integer :: k
         real(dp) :: dq_sum
         ierr = 0
         v_face(1) = v_mid(1)
         do k=2, s% nz
            dq_sum = s% dq(k-1) + s% dq(k)
            v_face(k) = (v_mid(k)*s% dq(k-1) + v_mid(k-1)*s% dq(k))/dq_sum
         end do
      end subroutine get_face_values


      real(dp) function get_Ledd(s,k) result(Ledd)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: kap_face
         integer :: j
         if (k == 1) then
            j = 2
         else
            j = k
         end if
         kap_face = interp_val_to_pt(s% opacity,j,s% nz,s% dq,'get_Ledd')
         Ledd = pi4*clight*s% cgrav(j)*s% m_grav(j)/kap_face
      end function get_Ledd


      real(dp) function get_Lrad(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         if (k == 1) then
            get_Lrad = s% L(k)
            return
         end if
         get_Lrad = s% L(k) - s% L_conv(k) ! L_conv set by last call on mlt
      end function get_Lrad


      real(dp) function get_Ladv(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: T, Erad, v, r
         T = s% T(k)
         Erad = crad*T*T*T*T
         if (s% u_flag) then
            v = s% u(k)
         else if (s% v_flag) then
            v = s% v(k)
         else
            v = s% r(1)*s% dlnR_dt(1)
         end if
         r = s% rmid(k)
         get_Ladv = 4*pi*r*r*v*Erad
      end function get_Ladv


      real(dp) function get_Lrad_div_Ledd(s,k) result(L_rad_div_Ledd)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer :: j
         real(dp) :: del_m, del_T4, area
         if (s% cgrav(k) <= 0) then
            L_rad_div_Ledd = 0d0
            return
         end if
         if (k == 1) then
            j = 2
         else
            j = k
         end if
         del_m = s% dm_bar(j)
         del_T4 = pow4(s% T(j-1)) - pow4(s% T(j))
         area = 4*pi*s% r(j)*s% r(j)
         L_rad_div_Ledd = &
            -(area*area*crad*(del_T4/del_m)/3)/(pi4*s% cgrav(j)*s% m_grav(j))
      end function get_Lrad_div_Ledd


      subroutine set_work_outward_at_surface(s)
         type (star_info), pointer :: s
         include 'formats'
         s% work_inward_at_center = 4*pi*s% r_center*s% r_center*s% P(s% nz)*s% v_center
         if (s% v_flag) then
            s% work_outward_at_surface = 4*pi*s% r(1)*s% r(1)*s% P_surf*s% v(1)
         else if (s% u_flag) then
            s% work_outward_at_surface = 0
         else
            s% work_outward_at_surface = 0
         end if
         if (s% use_eps_mdot) &
            s% work_outward_at_surface = &
               s% work_outward_at_surface - s% mdot_acoustic_surface + s% mdot_adiabatic_surface
      end subroutine set_work_outward_at_surface


      real(dp) function eval_rms_dvdt_div_v(s, klo, khi)
         type (star_info), pointer :: s
         integer, intent(in) :: klo, khi ! sum from klo to khi
         integer :: k
         real(dp) :: term, sum
         if (khi <= klo) then
            eval_rms_dvdt_div_v = 0d0
            return
         end if
         sum = 0
         do k=klo, khi
            term = s% dv_dt(k)/max(1d-50,abs(s% v(k)))
            sum = sum + term*term
         end do
         eval_rms_dvdt_div_v = sqrt(sum/(khi - klo + 1))
      end function eval_rms_dvdt_div_v
      
      
      real(dp) function cell_specific_KE(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: v
         if (s% u_flag) then
            cell_specific_KE = 0.5d0*s% u(k)*s% u(k)
         else if (s% v_flag) then
            if (k < s% nz) then
               v = 0.5d0*(s% v(k) + s% v(k+1))
            else
               v = 0.5d0*(s% v(k) + s% v_center)
            end if
            cell_specific_KE = 0.5d0*v*v
         else ! ignore kinetic energy if no velocity variables
            cell_specific_KE = 0d0
         end if
      end function cell_specific_KE
      
      
      real(dp) function cell_specific_PE(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: rL, rR, rC, m_cntr, rp13, r003
         if (s% u_flag .or. s% rsp_flag) then
            if (k == s% nz) then
               rL = s% R_center
            else
               rL = s% r(k+1)
            end if
            rR = s% r(k)
            rC = 0.5d0*(rR + rL)
         else ! do not use s% rmid(k) since call this from mesh_adjust
            r003 = s% r(k)*s% r(k)*s% r(k)
            if (k < s% nz) then
               rp13 = s% r(k+1)*s% r(k+1)*s% r(k+1)
            else
               rp13 = s% R_center*s% R_center*s% R_center
            end if
            rC = pow_cr(0.5d0*(r003 + rp13),1d0/3d0)
         end if
         m_cntr = s% m(k) - 0.5d0*s% dm(k)
         cell_specific_PE = -s% cgrav(k)*m_cntr/rC
      end function cell_specific_PE
      
      
      real(dp) function cell_specific_rotational_energy(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: e_00, e_p1
         e_00 = s% i_rot(k)*s% omega(k)*s% omega(k)
         if (k < s% nz) then
            e_p1 = s% i_rot(k+1)*s% omega(k+1)*s% omega(k+1)
         else
            e_p1 = 0
         end if
         cell_specific_rotational_energy = 0.5d0*(e_p1 + e_00)
      end function cell_specific_rotational_energy
      
      
      real(dp) function eval_deltaM_total_from_profile( &
            deltaM, premesh_dm, profile)
         real(dp), intent(in) :: deltaM
         real(dp), intent(in) :: premesh_dm(:), profile(:)
         real(dp) :: sum_dm, dm, total, f
         integer :: k
         include 'formats'
         total = 0
         sum_dm = 0
         do k=1,size(premesh_dm,dim=1)
            if (sum_dm >= deltaM) exit
            dm = premesh_dm(k)
            if (sum_dm + dm > deltaM) then
               f = (deltaM - sum_dm)/dm
               total = total + f*profile(k)
               exit
            end if
            total = total + profile(k)
            sum_dm = sum_dm + dm
         end do
         eval_deltaM_total_from_profile = total
      end function eval_deltaM_total_from_profile
      

      subroutine eval_deltaM_total_energy_integrals( &
            s, klo, khi, deltaM, save_profiles, &
            total_energy_profile, &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total)
         type (star_info), pointer :: s
         integer, intent(in) :: klo, khi ! sum from klo to khi
         real(dp), intent(in) :: deltaM
         logical, intent(in) :: save_profiles
         real(dp), intent(out), dimension(:) :: total_energy_profile
         real(dp), intent(out) :: &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total
         integer :: k
         real(dp) :: rR, rL, rC, dVR, dV, dm, m_cntr, rho, egas, v, sum_dm, &
            cell_total, cell1
         include 'formats'
         
         total_internal_energy = 0d0
         total_gravitational_energy = 0d0
         total_radial_kinetic_energy = 0d0
         total_rotational_kinetic_energy = 0d0
         total_turbulent_energy = 0d0
         sum_total = 0d0

         if (klo < 1 .or. khi > s% nz .or. klo > khi) return
         
         sum_dm = 0
         do k=klo,khi
            if (sum_dm >= deltaM) exit
            cell_total = 0
            dm = s% dm(k)
            if (sum_dm + dm > deltaM) dm = deltaM - sum_dm
            cell1 = dm*s% energy(k)
            cell_total = cell_total + cell1
            total_internal_energy = total_internal_energy + cell1
            if (s% v_flag .or. s% u_flag) then
               cell1 = dm*cell_specific_KE(s,k)
               cell_total = cell_total + cell1
               total_radial_kinetic_energy = total_radial_kinetic_energy + cell1
            end if
            if (.not. s% zero_gravity) then
               cell1 = dm*cell_specific_PE(s,k)
               cell_total = cell_total + cell1
               total_gravitational_energy = total_gravitational_energy + cell1
            end if
            if (s% rotation_flag) then
               cell1 = dm*cell_specific_rotational_energy(s,k)
               total_rotational_kinetic_energy = total_rotational_kinetic_energy + cell1
               if (s% include_rotation_in_total_energy) &
                  cell_total = cell_total + cell1
            end if
            if (s% rsp_flag) then
               cell1 = dm*s% Et(k)
               cell_total = cell_total + cell1
               total_turbulent_energy = total_turbulent_energy + cell1
            end if
            if (save_profiles) then
               total_energy_profile(k) = cell_total
            end if
         end do

         sum_total = total_internal_energy + total_gravitational_energy + &
            total_radial_kinetic_energy + total_turbulent_energy
            
         if (s% include_rotation_in_total_energy) &
            sum_total = sum_total + total_rotational_kinetic_energy

      end subroutine eval_deltaM_total_energy_integrals
      
      
      real(dp) function eval_cell_section_total_energy( &
            s, klo, khi) result(sum_total)
         type (star_info), pointer :: s
         integer, intent(in) :: klo, khi ! sum from klo to khi
         real(dp) :: &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy
         real(dp), pointer, dimension(:) :: total_energy_profile
         nullify(total_energy_profile)
         call eval_deltaM_total_energy_integrals( &
            s, klo, khi, s% mstar, .false., &
            total_energy_profile, &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total)
      end function eval_cell_section_total_energy


      subroutine eval_total_energy_integrals(s, &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total)
         type (star_info), pointer :: s
         real(dp), intent(out) :: &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total
         real(dp), pointer, dimension(:) :: total_energy_profile
         nullify(total_energy_profile)
         call eval_deltaM_total_energy_integrals( &
            s, 1, s% nz, s% mstar, .false., &
            total_energy_profile, &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total)
      end subroutine eval_total_energy_integrals


      real(dp) function get_total_energy_integral(s,k) result(sum_total)
         ! from surface down to k
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy
         real(dp), pointer, dimension(:) :: total_energy_profile
         nullify(total_energy_profile)
         call eval_deltaM_total_energy_integrals( &
            s, 1, k, s% mstar, .false., &
            total_energy_profile, &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total)
      end function get_total_energy_integral


      real(dp) function get_total_energy_integral_outward(s,k) result(sum_total)
         ! from surface down to k
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy
         real(dp), pointer, dimension(:) :: total_energy_profile
         nullify(total_energy_profile)
         call eval_deltaM_total_energy_integrals( &
            s, k, s% nz, s% mstar, .false., &
            total_energy_profile, &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total)
      end function get_total_energy_integral_outward


      real(dp) function get_log_concentration(s,j,k) result(log_c)
         use chem_def, only: chem_isos
         type (star_info), pointer :: s
         integer, intent(in) :: j, k
         ! concentration = number density / number density of electrons
         !  Ci = (Xi/Ai) / sum(Zi*Xi/Ai)   [see Thoul et al, ApJ 421:828-842, 1994]
         integer :: i, cid, species
         real(dp) :: tmp, c
         log_c = -1d99
         if (s% chem_id(j) == 0) return
         species = s% species
         tmp = 0d0
         do i=1,species
            cid = s% chem_id(i)
            tmp = tmp + chem_isos% Z(cid)*s% xa(i,k)/chem_isos% Z_plus_N(cid)
         end do
         cid = s% chem_id(j)
         c = (s% xa(j,k)/chem_isos% Z_plus_N(cid))/tmp
         log_c = safe_log10_cr(c)
      end function get_log_concentration


      real(dp) function get_phi_Joss(s,k) result(phi)
         use eos_def, only: i_lnPgas
         ! Joss, Salpeter, Ostriker, 1973. density inversion when Lrad/Ledd > phi.
         type (star_info), pointer :: s
         integer, intent(in) :: k
         phi = 1d0/(1d0 + (s% Pgas(k)/(4* s% Prad(k)))*s% d_eos_dlnT(i_lnPgas,k))
      end function get_phi_Joss


      logical function after_He_burn(s, he4_limit)
         use chem_def
         type (star_info), pointer :: s
         real(dp), intent(in) :: he4_limit
         integer :: nz, h1, he4
         real(dp) :: small = 1d-4
         after_He_burn = .false.
         nz = s% nz
         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         if (h1 == 0 .or. he4 == 0) return
         if (s% xa(h1,nz) > small .or. s% xa(he4,nz) > he4_limit) return
         after_He_burn = .true.
      end function after_He_burn


      logical function after_C_burn(s, c12_limit)
         use chem_def
         type (star_info), pointer :: s
         real(dp), intent(in) :: c12_limit
         integer :: nz, h1, he4, c12
         real(dp) :: small = 1d-4
         after_C_burn = .false.
         nz = s% nz
         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         c12 = s% net_iso(ic12)
         if (h1 == 0 .or. he4 == 0 .or. c12 == 0) return
         if (s% xa(h1,nz) > small .or. s% xa(he4,nz) > small .or. &
             s% xa(c12,nz) > c12_limit) return
         after_C_burn = .true.
      end function after_C_burn


      subroutine median_smoothing(dd, n, ns, dmed)
         use num_lib, only: qsort
         real(dp), intent(inout) :: dd(:) ! (n)
         integer, intent(in) :: n, ns
         real(dp), intent(inout) :: dmed(:) ! (n) work array

         real(dp) :: x(2*ns+1)
         integer :: i, j, k, nmed, index(2*ns+1)

         nmed = 2*ns+1

         do i=1,n
            if ((i > 1+ns) .and. (i < n-ns)) then
               k = 1
               do j = i-ns, i+ns
                  x(k) = dd(j)
                  k = k+1
               end do
               call qsort(index,nmed,x)
               dmed(i) = x(index(ns+1))
            else
               dmed(i) = dd(i)
            end if
         end do

         do i=1,n
            if (dmed(i) /= 0) dd(i) = dmed(i)
         end do

      end subroutine median_smoothing


      subroutine weighed_smoothing(dd, n, ns, preserve_sign, ddold)
      !     based on routine written by S.-C. Yoon, 18 Sept. 2002
      !     for smoothing  any variable (dd) with size n over 2*ns+1 cells.
         real(dp), intent(inout) :: dd(:) ! (n)
         integer, intent(in) :: n, ns
         logical, intent(in) :: preserve_sign
         real(dp), intent(inout) :: ddold(:) ! (n) work array

         integer :: nweight, mweight, i, j, k
         real(dp) :: weight(2*ns+1), sweight, v0

         include 'formats'

         do i = 1,n
           ddold(i) = dd(i)
         end do

         !--preparation for smoothing --------
         nweight = ns
         mweight = 2*nweight+1
         do i = 1,mweight
            weight(i) = 0d0
         end do
         weight(1) = 1d0
         do i = 1,mweight-1
            do j = i+1,2,-1
               weight(j) = weight(j) + weight(j-1)
            end do
         end do

         !--smoothing ------------------------
         do i=2,n-1
            sweight=0d0
            dd(i)=0d0
            v0 = ddold(i)
            do j = i, max(1,i-nweight), -1
               k=j-i+nweight+1
               if (preserve_sign .and. v0*ddold(j) <= 0) exit
               sweight = sweight+weight(k)
               dd(i) = dd(i)+ddold(j)*weight(k)
            end do
            do j = i+1, min(n,i+nweight)
               k=j-i+nweight+1
               if (preserve_sign .and. v0*ddold(j) <= 0) exit
               sweight = sweight+weight(k)
               dd(i) = dd(i)+ddold(j)*weight(k)
            end do
            if (sweight > 0) then
               sweight = 1d0/sweight
               dd(i) = dd(i)*sweight
            end if
         end do

      end subroutine weighed_smoothing

      
      subroutine threshold_smoothing (dd, dd_thresh, n, ns, preserve_sign, ddold)

        ! Same as weighed_smoothing, but only smooth contiguous regions where |dd| >= dd_thresh

        real(dp), intent(inout) :: dd(:)    ! (n)
        real(dp), intent(in)    :: dd_thresh
        integer, intent(in)     :: n
        integer, intent(in)     :: ns
        logical, intent(in)     :: preserve_sign
        real(dp), intent(inout) :: ddold(:) ! (n) work array

        logical :: in_region
        integer :: i
        integer :: i_a
        integer :: i_b
        
        include 'formats'

        ! Process regions

        in_region = .FALSE.

        i_a = 1
        do i = 1, n

           if (in_region) then

              if (ABS(dd(i)) < dd_thresh) then
                 i_b = i-1
                 if (i_b > i_a) call weighed_smoothing(dd(i_a:i_b), i_b-i_a+1, ns, preserve_sign, ddold(i_a:i_b))
                 in_region = .FALSE.
              endif

           else
              if (ABS(dd(i)) >= dd_thresh) then
                 i_a = i
                 in_region = .TRUE.
              endif

           end if

        end do

        ! Handle the final region

        if (in_region) then

           i_b = n
           if (i_b > i_a) call weighed_smoothing(dd(i_a:i_b), i_b-i_a+1, ns, preserve_sign, ddold(i_a:i_b))

        endif

        ! Finish

        return

      end subroutine threshold_smoothing
      

      real(dp) function eval_kh_timescale(G,M,R,L) result(kh)
         real(dp), intent(in) :: G,M,R,L
         if (L <= 0) then
            kh = 0d0
         else
            kh = 0.75d0*G*M*M/(R*L) ! 0.75 is based on sun.  Hansen & Kawaler eqn 1.30
         end if
      end function eval_kh_timescale



      real(dp) function yrs_for_init_timestep(s)
         type (star_info), pointer :: s
         if (s% initial_mass <= 1) then
            yrs_for_init_timestep = 1d5
         else
            yrs_for_init_timestep = 1d5 / pow_cr(s% initial_mass,2.5d0)
         end if
      end function yrs_for_init_timestep


b
      end module star_utils
