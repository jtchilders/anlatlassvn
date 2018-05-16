module     p12_cbbar_hepneg_matrix
   use msamurai, only: initsamurai, exitsamurai
   use p12_cbbar_hepneg_util, only: square
   use p12_cbbar_hepneg_config, only: ki, &
     & include_helicity_avg_factor, include_color_avg_factor, &
     & debug_lo_diagrams, debug_nlo_diagrams, &
     & include_symmetry_factor, &
     & PSP_check, PSP_verbosity, PSP_rescue, PSP_chk_threshold1, &
     & PSP_chk_threshold2, PSP_chk_kfactor, reduction_interoperation, &
     & convert_to_cdr, &
     & samurai_verbosity, samurai_test, samurai_scalar
   use p12_cbbar_hepneg_kinematics, only: &
       in_helicities, symmetry_factor, num_legs, &
       lo_qcd_couplings, corrections_are_qcd, num_light_quarks, num_gluons
   use p12_cbbar_hepneg_model, only: Nf, NC, sqrt2, init_functions
   use p12_cbbar_hepneg_color, only: TR, CA, CF, numcs, &
     & incolors, init_color
   use p12_cbbar_hepneg_diagramsh0l0, only: amplitude0l0 => amplitude
   use p12_cbbar_hepneg_amplitudeh0, only: samplitudeh0l1 => samplitude, &
        &   finite_renormalisation0 => finite_renormalisation
   use p12_cbbar_hepneg_dipoles, only: insertion_operator

   implicit none
   save

   private

   integer :: banner_ch = 6

   public :: initgolem, exitgolem, samplitude
   public :: samplitudel0, samplitudel1
   public :: ir_subtraction, color_correlated_lo2, spin_correlated_lo2
contains
   !---#[ subroutine banner:
   subroutine     banner()
      implicit none

      character(len=72) :: frame = "+" // repeat("-", 70) // "+"

      if (banner_ch .le. 0) return

      write(banner_ch,'(A72)') frame
      write(banner_ch,'(A72)') "|   __   __   ___   __   __  __                  GoSam                 |"
      write(banner_ch,'(A72)') "|  / _) /  \ / __) (  ) (  \/  )         An Automated One-Loop         |"
      write(banner_ch,'(A72)') "| ( (/\( () )\__ \ /__\  )    (         Matrix Element Generator       |"
      write(banner_ch,'(A72)') "|  \__/ \__/ (___/(_)(_)(_/\/\_)          Version 1.0 Rev: 228         |"
      write(banner_ch,'(A72)') "|                                                                      |"
      write(banner_ch,'(A72)') "|                                (c) The GoSam Collaboration 2011-2013 |"
      write(banner_ch,'(A72)') "|                                                                      |"
      write(banner_ch,'(A72)') "|                AUTHORS:                                              |"
      write(banner_ch,'(A72)') "|                * Gavin Cullen         <gavin.cullen@desy.de>         |"
      write(banner_ch,'(A72)') "|                * Nicolas Greiner      <greiner@mpp.mpg.de>           |"
      write(banner_ch,'(A72)') "|                * Gudrun Heinrich      <gudrun@mpp.mpg.de>            |"
      write(banner_ch,'(A72)') "|                * Gionata Luisoni      <luisonig@mpp.mpg.de>          |"
      write(banner_ch,'(A72)') "|                * Pierpaolo Mastrolia  <pierpaolo.mastrolia@cern.ch>  |"
      write(banner_ch,'(A72)') "|                * Giovanni Ossola      <gossola@citytech.cuny.edu>    |"
      write(banner_ch,'(A72)') "|                * Thomas Reiter        <reiterth@mpp.mpg.de>          |"
      write(banner_ch,'(A72)') "|                * Francesco Tramontano <francesco.tramontano@cern.ch> |"
      write(banner_ch,'(A72)') "|                                                                      |"
      write(banner_ch,'(A72)') "|                                                                      |"
      write(banner_ch,'(A72)') "|    This program is free software: you can redistribute it and/or modify|"
      write(banner_ch,'(A72)') "|    it under the terms of the GNU General Public License as published by|"
      write(banner_ch,'(A72)') "|    the Free Software Foundation, either version 3 of the License, or |"
      write(banner_ch,'(A72)') "|    (at your option) any later version.                               |"
      write(banner_ch,'(A72)') "|                                                                      |"
      write(banner_ch,'(A72)') "|    This program is distributed in the hope that it will be useful,   |"
      write(banner_ch,'(A72)') "|    but WITHOUT ANY WARRANTY; without even the implied warranty of    |"
      write(banner_ch,'(A72)') "|    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     |"
      write(banner_ch,'(A72)') "|    GNU General Public License for more details.                      |"
      write(banner_ch,'(A72)') "|                                                                      |"
      write(banner_ch,'(A72)') "|    You should have received a copy of the GNU General Public License |"
      write(banner_ch,'(A72)') "|    along with this program.  If not, see <http://www.gnu.org/licenses/>.|"
      write(banner_ch,'(A72)') "|                                                                      |"
      write(banner_ch,'(A72)') "|    Scientific publications prepared using the present version of     |"
      write(banner_ch,'(A72)') "|    GoSam or any modified version of it or any code linking to        |"
      write(banner_ch,'(A72)') "|    GoSam or parts of it should make a clear reference to the publication:|"
      write(banner_ch,'(A72)') "|                                                                      |"
      write(banner_ch,'(A72)') "|        G. Cullen et al.,                                             |"
      write(banner_ch,'(A72)') "|        ``Automated One-Loop Calculations with GoSam'',               |"
      write(banner_ch,'(A72)') "|        arXiv:1111.2034 [hep-ph]                                      |"
      write(banner_ch,'(A72)') frame

      banner_ch = 0
   end subroutine banner
   !---#] subroutine banner:

   !---#[ subroutine initgolem :
   subroutine     initgolem(is_first,stage,rndseed)
      implicit none
      logical, optional, intent(in) :: is_first
      integer, optional, intent(in) :: stage
      integer, optional, intent(in) :: rndseed
      logical :: init_third_party
      logical :: file_exists, dir_exists
      integer i, j
      character(len=50) :: file_name
      character(len=9)  :: dir_name = "BadPoints"
      character(len=6)  :: file_numb
      character(len=9)  :: file_pre = "gs_badpts"
      character(len=3)  :: file_ext = "log"
      character(len=1)  :: cstage
      character(len=4)  :: crndseed
      i = 1
      file_exists =.true.

      if(present(is_first)) then
         init_third_party = is_first
      else
         init_third_party = .true.
      end if
      if (init_third_party) then
         call initsamurai('diag',samurai_scalar,&
         &                samurai_verbosity,samurai_test)
      ! call our banner
      call banner()
      if(PSP_check .and. PSP_rescue .and. PSP_verbosity .ge. 1) then
         inquire(file=dir_name, exist=dir_exists)
         if(.not. dir_exists) then
            call system('mkdir BadPoints')
         end if
         if(present(stage)) then
            write(cstage,'(i1)') stage
            write(crndseed,'(i4)') rndseed
            do j=1,4
               if(crndseed(j:j).eq.' ') crndseed(j:j)='0'
            enddo
            file_name = dir_name//"/"//file_pre//"-"//cstage//"-"//crndseed//"."//file_ext
            open(unit=42, file=file_name, status='replace', action='write')
            write(42,'(A22)') "<?xml version='1.0' ?>"
            write(42,'(A5)')  "<run>"
         else
            do while(file_exists)
               write(file_numb, '(I6.1)') i
               file_name = dir_name//"/"//file_pre//trim(adjustl(file_numb))//"."//file_ext
               inquire(file=file_name, exist=file_exists)
               if(file_exists) then
                  write(*,*) "File ", file_name, " already exists!"
                  i = i+1
               else
                  write(*,*) "Bad points stored in file: ", file_name
                  open(unit=42, file=file_name, status='unknown', action='write')
                  write(42,'(A22)') "<?xml version='1.0' ?>"
                  write(42,'(A5)')  "<run>"
               end if
            enddo
         end if
      end if
      end if

      call init_functions()
      call init_color()

   end subroutine initgolem
   !---#] subroutine initgolem :
   !---#[ subroutine exitgolem :
   subroutine     exitgolem(is_last)
      use p12_cbbar_hepneg_groups, only: tear_down_golem95
      implicit none
      logical, optional, intent(in) :: is_last

      logical :: exit_third_party

      if(present(is_last)) then
         exit_third_party = is_last
      else
         exit_third_party = .true.
      end if
      if (exit_third_party) then
         call exitsamurai()
         call tear_down_golem95()
         if(PSP_check .and. PSP_rescue) then
            write(42,'(A6)')  "</run>"
            close(unit=42)
         endif
      end if
   end subroutine exitgolem
   !---#] subroutine exitgolem :

   !---#[ subroutine samplitude :
   subroutine     samplitude(vecs, scale2, amp, ok, h)
      implicit none
      real(ki), dimension(6, 4), intent(in) :: vecs
      real(ki), intent(in) :: scale2
      real(ki), dimension(1:4), intent(out) :: amp
      logical, intent(out), optional :: ok
      integer, intent(in), optional :: h
      real(ki) :: rat2, sam_amp2, sam_amp3, kfac, zero
      integer spprec1, fpprec1, i
      real(ki), dimension(2:3) :: irp
      integer :: tmp_red_int, spprec2, fpprec2 
      call samplitudel01(vecs, scale2, amp, rat2, ok, h)
      if(PSP_check) then
      tmp_red_int=reduction_interoperation
      call ir_subtraction(vecs, scale2, irp)
      if((amp(3)-irp(2)) .ne. 0.0_ki) then
         spprec1 = -int(log10(abs((amp(3)-irp(2))/irp(2))))
      else
         spprec1 = 16
      endif
      if((amp(2)-rat2) .ne. 0.0_ki) then
         fpprec1 = spprec1 + int(log10(abs(amp(2)/(amp(2)-rat2))))
      else
         fpprec1 = -10
      endif
      if(amp(1) .ne. 0.0_ki) then
         kfac = amp(2)/amp(1)
      else
         kfac = 0.0_ki
      endif
      if(spprec1 .lt. PSP_chk_threshold1 .and. spprec1 .gt. -10000 .or. &
           (abs(kfac) > PSP_chk_kfactor .and. PSP_chk_kfactor > 0)) then
         if(PSP_verbosity .eq. 2) write(*,*) "UNSTABLE PHASE SPACE POINT !!"
         if(PSP_rescue) then
            reduction_interoperation = 1
            sam_amp2 = amp(2)
            sam_amp3 = amp(3)
            call samplitudel01(vecs, scale2, amp, rat2, ok, h)
            if((amp(3)-irp(2)) .ne. 0.0_ki) then
               spprec2 = -int(log10(abs((amp(3)-irp(2))/irp(2))))
            else
               spprec2 = 16
            endif
            if((amp(2)-rat2) .ne. 0.0_ki) then
               fpprec2 = spprec2 + int(log10(abs(amp(2)/(amp(2)-rat2))))
            else
               fpprec2 = -10
            endif
            if(amp(1) .ne. 0.0_ki) then
               kfac = amp(2)/amp(1)
            else
               kfac = 0.0_ki
            endif
            if(spprec2 .le. PSP_chk_threshold2 .and. spprec2 .gt. -10000 .or. &
               (abs(kfac) > PSP_chk_kfactor .and. PSP_chk_kfactor > 0)) then
               if(PSP_verbosity .ge. 1) then
                  write(42,'(2x,A7)')"<event>"
                  if(spprec2 .le. PSP_chk_threshold2 .and. spprec2 .gt. -10000) then
                     write(42,'(4x,A15,A16,A18,A3)') "<process name='", &
                          &   "p12_cbbar_hepneg", "' problem='sinpole","'/>"
                  else if(abs(kfac) > PSP_chk_kfactor) then
                     write(42,'(4x,A15,A16,A18,A3)') "<process name='", &
                          &   "p12_cbbar_hepneg", "' problem='kfactor","'/>"
                  end if
                  write(42,'(4x,A27,I2.1,A14,I2.1,A3)') "<pspThresholds threshold1='", &
                       &   PSP_chk_threshold1, "' threshold2='", PSP_chk_threshold2, "'/>"
                  write(42,'(4x,A17,I2.1,A10,I2.1,A3)') "<precSam spprec='", &
                       &   spprec1, "' fpprec='", fpprec1, "'/>"
                  write(42,'(4x,A17,I2.1,A10,I2.1,A3)') "<precGol spprec='", &
                       &   spprec2, "' fpprec='", fpprec2, "'/>"
                  write(42,'(4x,A18,D23.16,A7,D23.16,A6,D23.16,A3)') "<singlePoles sam='", sam_amp3, &
                       &   "' gol='", amp(3), "' ir='", irp(2),"'/>"
                  write(42,'(4x,A17,D23.16,A8,D23.16,2(A7,D23.16),A3)') "<amplitude born='", amp(1), &
                       &   "' rat2='", rat2, "' sam='", sam_amp2, "' gol='", amp(2), "'/>"
                  write(42,'(4x,A9)') "<momenta>"
                  do i=1,6
                     write(42,'(8x,A8,3(D23.16,A6),D23.16,A3)') "<mom e='", vecs(i,1), "' px='", vecs(i,2), &
                          &     "' py='", vecs(i,3), "' pz='", vecs(i,4), "'/>"
                  enddo
                  write(42,'(4x,A10)')"</momenta>"
                  write(42,'(2x,A8)')"</event>"
               endif
               ! Give back a Nan so that point is discarded
               zero = log(1.0_ki)
               amp(2)= 1.0_ki/zero
            else
               if(PSP_verbosity .eq. 2) write(*,*) "POINT SAVED !!"
               if(PSP_verbosity .ge. 2) write(*,*)
            end if
            reduction_interoperation = tmp_red_int
         end if
      end if
   end if
   end subroutine samplitude
   !---#] subroutine samplitude :

   !---#[ subroutine samplitudel01 :
   subroutine     samplitudel01(vecs, scale2, amp, rat2, ok, h)
      use p12_cbbar_hepneg_config, only: &
         & debug_lo_diagrams, debug_nlo_diagrams, logfile, deltaOS, &
         & renormalisation, renorm_beta, renorm_mqwf, renorm_decoupling, &
         & renorm_logs, renorm_mqse, nlo_prefactors
      use p12_cbbar_hepneg_kinematics, only: &
         & inspect_kinematics, init_event
      use p12_cbbar_hepneg_model
      use p12_cbbar_hepneg_dipoles, only: pi
      implicit none
      real(ki), dimension(6, 4), intent(in) :: vecs
      real(ki), intent(in) :: scale2
      real(ki), dimension(4), intent(out) :: amp
      real(ki), intent(out) :: rat2
      logical, intent(out), optional :: ok
      integer, intent(in), optional :: h
      real(ki) :: nlo_coupling

      complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)

      ! Number of heavy quark flavours in loops.
      real(ki), parameter :: NFh = 0.0_ki

      logical :: my_ok

      ! used for m=0 QCD renormalisation
      real(ki) :: beta0

      if(corrections_are_qcd) then
         nlo_coupling = 1.0_ki
      else
         nlo_coupling = 1.0_ki
      end if

      call init_event(vecs)

      if(debug_lo_diagrams .or. debug_nlo_diagrams) then
         write(logfile,'(A7)') "<event>"
         call inspect_kinematics(logfile)
      end if

      
      if (present(h)) then
         amp(1) = samplitudel0(vecs, h)
      else
         amp(1)   = samplitudel0(vecs)
      end if
      select case (renormalisation)
      case (0)
         ! no renormalisation
         deltaOS = 0.0_ki
      case (1)
         ! fully renormalized
         if(renorm_mqse) then
            deltaOS = 1.0_ki
         else
            deltaOS = 0.0_ki
         end if
      case (2)
         ! massive quark counterterms only
         deltaOS = 1.0_ki
      case default
         ! not implemented
         print*, "In p12_cbbar_hepneg_matrix:"
         print*, "  invalid value for renormalisation=", renormalisation
         stop
      end select

      if (present(h)) then
         amp((/4,3,2/)) = samplitudel1(vecs, scale2, my_ok, rat2, h)/nlo_coupling
      else
         amp((/4,3,2/)) = samplitudel1(vecs, scale2, my_ok, rat2)/nlo_coupling
      end if
      select case (renormalisation)
      case (0)
         ! no renormalisation
      case (1)
         ! fully renormalized
         if(corrections_are_qcd) then
            if (renorm_beta) then
               beta0 = (11.0_ki * CA - 4.0_ki * TR * (NF + NFh)) / 6.0_ki
               amp(3) = amp(3) - lo_qcd_couplings * beta0 * amp(1)
               amp(2) = amp(2) + lo_qcd_couplings * CA / 6.0_ki * amp(1)
            end if
            if (renorm_mqwf) then
            end if
         end if
      case (2)
         ! massive quark counterterms only
      case default
         ! not implemented
         print*, "In p12_cbbar_hepneg_matrix:"
         print*, "  invalid value for renormalisation=", renormalisation
         stop
      end select
      if (convert_to_cdr) then
         ! Scheme conversion for infrared structure
         ! Reference:
         ! S. Catani, M. H. Seymour, Z. Trocsanyi,
         ! ``Regularisation scheme independence and unitarity
         !   in QCD cross-sections,''
         ! Phys.Rev. D 55 (1997) 6819
         ! arXiv:hep-ph/9610553
         amp(2) = amp(2) - amp(1) * (&
           &          num_light_quarks * 0.5_ki * CF &
           &        + num_gluons * 1.0_ki/6.0_ki * CA)
      end if
      if (present(ok)) ok = my_ok

      if(debug_lo_diagrams .or. debug_nlo_diagrams) then
         write(logfile,'(A25,E24.16,A3)') &
            & "<result kind='lo' value='", amp(1), "'/>"
         write(logfile,'(A33,E24.16,A3)') &
            & "<result kind='nlo-finite' value='", amp(2), "'/>"
         write(logfile,'(A33,E24.16,A3)') &
            & "<result kind='nlo-single' value='", amp(3), "'/>"
         write(logfile,'(A33,E24.16,A3)') &
            & "<result kind='nlo-double' value='", amp(4), "'/>"
         if(my_ok) then
            write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
         else
            write(logfile,'(A29)') "<flag name='ok' status='no'/>"
         end if
         write(logfile,'(A8)') "</event>"
      end if
      select case(nlo_prefactors)
      case(0)
         ! The result is already in its desired form
      case(1)
         amp(2:4) = amp(2:4) * nlo_coupling
      case(2)
         amp(2:4) = amp(2:4) * nlo_coupling / 8.0_ki / pi / pi
      end select
   end subroutine samplitudel01
   !---#] subroutine samplitudel01 :
   !---#[ function samplitudel0 :
   function     samplitudel0(vecs, h) result(amp)
      use p12_cbbar_hepneg_config, only: logfile
      use p12_cbbar_hepneg_kinematics, only: init_event
      implicit none
      real(ki), dimension(6, 4), intent(in) :: vecs
      integer, optional, intent(in) :: h
      real(ki) :: amp, heli_amp
      complex(ki), dimension(numcs) :: color_vector
      logical, dimension(0:31) :: eval_heli
      real(ki), dimension(6, 4) :: pvecs

      if (present(h)) then
         eval_heli(:) = .false.
         eval_heli(h) = .true.
      else
         eval_heli(:) = .true.
      end if

      amp = 0.0_ki
      if (eval_heli(0)) then
         if (debug_lo_diagrams) then
            write(logfile,*) "<helicity index='0' >"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(3,:)
         pvecs(4,:) = vecs(4,:)
         pvecs(5,:) = vecs(5,:)
         pvecs(6,:) = vecs(6,:)
         call init_event(pvecs, -1)
         !---#] reinitialize kinematics:
         color_vector = amplitude0l0()
         heli_amp = square(color_vector)
         if (debug_lo_diagrams) then
            write(logfile,'(A25,E24.16,A3)') &
                & "<result kind='lo' value='", heli_amp, "'/>"
            write(logfile,*) "</helicity>"
         end if
         amp = amp + heli_amp
      end if
      if (eval_heli(1)) then
         if (debug_lo_diagrams) then
            write(logfile,*) "<helicity index='1' >"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(3,:)
         pvecs(4,:) = vecs(4,:)
         pvecs(5,:) = vecs(5,:)
         pvecs(6,:) = vecs(6,:)
         call init_event(pvecs, +1)
         !---#] reinitialize kinematics:
         color_vector = amplitude0l0()
         heli_amp = square(color_vector)
         if (debug_lo_diagrams) then
            write(logfile,'(A25,E24.16,A3)') &
                & "<result kind='lo' value='", heli_amp, "'/>"
            write(logfile,*) "</helicity>"
         end if
         amp = amp + heli_amp
      end if
      if (include_helicity_avg_factor) then
         amp = amp / real(in_helicities, ki)
      end if
      if (include_color_avg_factor) then
         amp = amp / incolors
      end if
      if (include_symmetry_factor) then
         amp = amp / real(symmetry_factor, ki)
      end if
   end function samplitudel0
   !---#] function samplitudel0 :
   !---#[ function samplitudel1 :
   function     samplitudel1(vecs,scale2,ok,rat2,h) result(amp)
      use p12_cbbar_hepneg_config, only: &
         & debug_nlo_diagrams, logfile, renorm_gamma5
      use p12_cbbar_hepneg_kinematics, only: init_event
      implicit none
      real(ki), dimension(6, 4), intent(in) :: vecs
      logical, intent(out) :: ok
      real(ki), intent(in) :: scale2
      real(ki), intent(out) :: rat2
      integer, optional, intent(in) :: h
      real(ki), dimension(6, 4) :: pvecs
      real(ki), dimension(-2:0) :: amp, heli_amp
      logical :: my_ok
      logical, dimension(0:31) :: eval_heli
      real(ki) :: fr, rational2

      if (present(h)) then
         eval_heli(:) = .false.
         eval_heli(h) = .true.
      else
         eval_heli(:) = .true.
      end if

      amp(:) = 0.0_ki
      rat2 = 0.0_ki
      ok = .true.
      if (eval_heli(0)) then
         if(debug_nlo_diagrams) then
            write(logfile,*) "<helicity index='0'>"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(3,:)
         pvecs(4,:) = vecs(4,:)
         pvecs(5,:) = vecs(5,:)
         pvecs(6,:) = vecs(6,:)
         call init_event(pvecs, -1)
         !---#] reinitialize kinematics:
         heli_amp = samplitudeh0l1(real(scale2,ki),my_ok,rational2)
         if (corrections_are_qcd .and. renorm_gamma5) then
            !---#[ reinitialize kinematics:
            pvecs(1,:) = vecs(1,:)
            pvecs(2,:) = vecs(2,:)
            pvecs(3,:) = vecs(3,:)
            pvecs(4,:) = vecs(4,:)
            pvecs(5,:) = vecs(5,:)
            pvecs(6,:) = vecs(6,:)
            call init_event(pvecs, -1)
            !---#] reinitialize kinematics:
            fr = finite_renormalisation0(real(scale2,ki))
            heli_amp(0) = heli_amp(0) + fr
         end if
         ok = ok .and. my_ok
         amp = amp + heli_amp
         rat2 = rat2 + rational2

         if(debug_nlo_diagrams) then
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-finite' value='", heli_amp(0), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-single' value='", heli_amp(-1), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-double' value='", heli_amp(-2), "'/>"
            if (corrections_are_qcd .and. renorm_gamma5) then
               write(logfile,'(A30,E24.16,A3)') &
                   & "<result kind='fin-ren' value='", fr, "'/>"
            end if
            if(my_ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</helicity>"
         end if
      end if
      if (eval_heli(1)) then
         if(debug_nlo_diagrams) then
            write(logfile,*) "<helicity index='1'>"
         end if
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(3,:)
         pvecs(4,:) = vecs(4,:)
         pvecs(5,:) = vecs(5,:)
         pvecs(6,:) = vecs(6,:)
         call init_event(pvecs, +1)
         !---#] reinitialize kinematics:
         heli_amp = samplitudeh0l1(real(scale2,ki),my_ok,rational2)
         if (corrections_are_qcd .and. renorm_gamma5) then
            !---#[ reinitialize kinematics:
            pvecs(1,:) = vecs(1,:)
            pvecs(2,:) = vecs(2,:)
            pvecs(3,:) = vecs(3,:)
            pvecs(4,:) = vecs(4,:)
            pvecs(5,:) = vecs(5,:)
            pvecs(6,:) = vecs(6,:)
            call init_event(pvecs, +1)
            !---#] reinitialize kinematics:
            fr = finite_renormalisation0(real(scale2,ki))
            heli_amp(0) = heli_amp(0) + fr
         end if
         ok = ok .and. my_ok
         amp = amp + heli_amp
         rat2 = rat2 + rational2

         if(debug_nlo_diagrams) then
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-finite' value='", heli_amp(0), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-single' value='", heli_amp(-1), "'/>"
            write(logfile,'(A33,E24.16,A3)') &
                & "<result kind='nlo-double' value='", heli_amp(-2), "'/>"
            if (corrections_are_qcd .and. renorm_gamma5) then
               write(logfile,'(A30,E24.16,A3)') &
                   & "<result kind='fin-ren' value='", fr, "'/>"
            end if
            if(my_ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</helicity>"
         end if
      end if
      if (include_helicity_avg_factor) then
         amp = amp / real(in_helicities, ki)
      end if
      if (include_color_avg_factor) then
         amp = amp / incolors
      end if
      if (include_symmetry_factor) then
         amp = amp / real(symmetry_factor, ki)
      end if
   end function samplitudel1
   !---#] function samplitudel1 :
   !---#[ subroutine ir_subtraction :
   subroutine     ir_subtraction(vecs,scale2,amp,h)
      use p12_cbbar_hepneg_config, only: &
         & nlo_prefactors
      use p12_cbbar_hepneg_dipoles, only: pi
      use p12_cbbar_hepneg_kinematics, only: &
         & init_event, corrections_are_qcd
      use p12_cbbar_hepneg_model
      implicit none
      real(ki), dimension(6, 4), intent(in) :: vecs
      real(ki), intent(in) :: scale2
      integer, optional, intent(in) :: h
      real(ki), dimension(2), intent(out) :: amp
      real(ki), dimension(2) :: heli_amp
      real(ki), dimension(6, 4) :: pvecs
      complex(ki), dimension(numcs,numcs,2) :: oper
      complex(ki), dimension(numcs) :: color_vectorl0, pcolor
      logical, dimension(0:31) :: eval_heli
      real(ki) :: nlo_coupling

      if (present(h)) then
         eval_heli(:) = .false.
         eval_heli(h) = .true.
      else
         eval_heli(:) = .true.
      end if

      call init_event(vecs)

      if(corrections_are_qcd) then
         nlo_coupling = 1.0_ki
      else
         nlo_coupling = 1.0_ki
      end if

      oper = insertion_operator(real(scale2,ki), vecs)
      amp(:) = 0.0_ki
      if (eval_heli(0)) then
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(3,:)
         pvecs(4,:) = vecs(4,:)
         pvecs(5,:) = vecs(5,:)
         pvecs(6,:) = vecs(6,:)
         call init_event(pvecs, -1)
         !---#] reinitialize kinematics:
         pcolor = amplitude0l0()
         color_vectorl0(1) = pcolor(1)
         heli_amp(1) = square(color_vectorl0, oper(:,:,1))
         heli_amp(2) = square(color_vectorl0, oper(:,:,2))
         amp = amp + heli_amp
      endif
      if (eval_heli(1)) then
         !---#[ reinitialize kinematics:
         pvecs(1,:) = vecs(1,:)
         pvecs(2,:) = vecs(2,:)
         pvecs(3,:) = vecs(3,:)
         pvecs(4,:) = vecs(4,:)
         pvecs(5,:) = vecs(5,:)
         pvecs(6,:) = vecs(6,:)
         call init_event(pvecs, +1)
         !---#] reinitialize kinematics:
         pcolor = amplitude0l0()
         color_vectorl0(1) = pcolor(1)
         heli_amp(1) = square(color_vectorl0, oper(:,:,1))
         heli_amp(2) = square(color_vectorl0, oper(:,:,2))
         amp = amp + heli_amp
      endif
      if (include_helicity_avg_factor) then
         amp = amp / real(in_helicities, ki)
      end if
      if (include_color_avg_factor) then
         amp = amp / incolors
      end if
      if (include_symmetry_factor) then
         amp = amp / real(symmetry_factor, ki)
      end if
      select case(nlo_prefactors)
      case(0)
         ! The result is already in its desired form
      case(1)
         amp(:) = amp(:) * nlo_coupling
      case(2)
         amp(:) = amp(:) * nlo_coupling / 8.0_ki / pi / pi
      end select
   end subroutine ir_subtraction
   !---#] subroutine ir_subtraction :
   !---#[ color correlated ME :
   pure subroutine color_correlated_lo(color_vector,res)
      use p12_cbbar_hepneg_color, only: T1T1, &
      & T1T2, &
      & T1T6, &
      & T2T2, &
      & T2T6, &
      & T6T6
      implicit none
      complex(ki), dimension(numcs), intent(in) :: color_vector
      real(ki), dimension(num_legs,num_legs), intent(out) :: res
      res(:,:)=0.0_ki
      res(1,1) = square(color_vector,T1T1)
      res(1,1) = res(1,1)
      res(1,2) = square(color_vector,T1T2)
      res(2,1) = res(1,2)
      res(1,6) = square(color_vector,T1T6)
      res(6,1) = res(1,6)
      res(2,2) = square(color_vector,T2T2)
      res(2,2) = res(2,2)
      res(2,6) = square(color_vector,T2T6)
      res(6,2) = res(2,6)
      res(6,6) = square(color_vector,T6T6)
      res(6,6) = res(6,6)
   end subroutine color_correlated_lo

   subroutine     color_correlated_lo2(vecs,borncc)
      use p12_cbbar_hepneg_kinematics, only: init_event
      implicit none
      real(ki), dimension(num_legs, 4), intent(in) :: vecs
      real(ki), dimension(num_legs,num_legs), intent(out) :: borncc
      real(ki), dimension(num_legs,num_legs) :: borncc_heli
      real(ki), dimension(num_legs, 4) :: pvecs
      complex(ki), dimension(numcs) :: color_vector

      borncc(:,:) = 0.0_ki
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(3,:)
      pvecs(4,:) = vecs(4,:)
      pvecs(5,:) = vecs(5,:)
      pvecs(6,:) = vecs(6,:)
      call init_event(pvecs, -1)
      !---#] reinitialize kinematics:
      color_vector = amplitude0l0()
      call color_correlated_lo(color_vector,borncc_heli)
      ! The minus is part in the definition according to PowHEG Box.
      ! Since they use it we include it:
      borncc(:,:) = borncc(:,:) - borncc_heli(:,:)
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(3,:)
      pvecs(4,:) = vecs(4,:)
      pvecs(5,:) = vecs(5,:)
      pvecs(6,:) = vecs(6,:)
      call init_event(pvecs, +1)
      !---#] reinitialize kinematics:
      color_vector = amplitude0l0()
      call color_correlated_lo(color_vector,borncc_heli)
      borncc(:,:) = borncc(:,:) - borncc_heli(:,:)
      if (include_helicity_avg_factor) then
         borncc = borncc / real(in_helicities, ki)
      end if
      if (include_color_avg_factor) then
         borncc = borncc / incolors
      end if
      if (include_symmetry_factor) then
         borncc = borncc / real(symmetry_factor, ki)
      end if
   end subroutine color_correlated_lo2
   !---#] color correlated ME :
   !---#[ spin correlated ME :
   subroutine spin_correlated_lo2(vecs, bornsc)
      use p12_cbbar_hepneg_kinematics
      implicit none
      real(ki), dimension(num_legs, 4), intent(in) :: vecs
      real(ki), dimension(num_legs,4,4) :: bornsc
      real(ki), dimension(num_legs, 4) :: pvecs
      complex(ki), dimension(4,4) :: tens
      complex(ki) :: pp, pm, mp, mm
      complex(ki), dimension(numcs) :: heli_amp0
      complex(ki), dimension(numcs) :: heli_amp1
      complex(ki), dimension(4) :: eps6

      bornsc(:,:,:) = 0.0_ki
      !---#[ Initialize helicity amplitudes :
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(3,:)
      pvecs(4,:) = vecs(4,:)
      pvecs(5,:) = vecs(5,:)
      pvecs(6,:) = vecs(6,:)
      call init_event(pvecs, -1)
      !---#] reinitialize kinematics:
      heli_amp0 = amplitude0l0()
      !---#[ reinitialize kinematics:
      pvecs(1,:) = vecs(1,:)
      pvecs(2,:) = vecs(2,:)
      pvecs(3,:) = vecs(3,:)
      pvecs(4,:) = vecs(4,:)
      pvecs(5,:) = vecs(5,:)
      pvecs(6,:) = vecs(6,:)
      call init_event(pvecs, +1)
      !---#] reinitialize kinematics:
      heli_amp1 = amplitude0l0()
      !---#] Initialize helicity amplitudes :
      !---#[ Initialize polarization vectors :
      eps6 = conjg(spvak2k6/Spaa(k2,k6)/sqrt2)
      !---#] Initialize polarization vectors :
      ! Note: By omitting the imaginary parts we lose a term:
      !   Imag(B_j(mu,nu)) = i_ * e_(k_j, mu, q_j, nu) * |Born|^2
      ! where q_j is the reference momentum chosen for the paticle
      ! of momentum k_j. This term should, however not be phenomenologically
      ! relevant.
      !---#[ particle 6 :
      pp  = 0.0_ki &
      &          + square_0l_0l_sc(heli_amp1,heli_amp1)
      pm  = 0.0_ki &
      &          + square_0l_0l_sc(heli_amp1,heli_amp0)
      mp  = 0.0_ki &
      &          + square_0l_0l_sc(heli_amp0,heli_amp1)
      mm  = 0.0_ki &
      &          + square_0l_0l_sc(heli_amp0,heli_amp0)

      call construct_polarization_tensor(conjg(eps6),eps6,tens)
      bornsc(6,:,:) = bornsc(6,:,:) + real(tens(:,:) * pp, ki)
      call construct_polarization_tensor(conjg(eps6),conjg(eps6),tens)
      bornsc(6,:,:) = bornsc(6,:,:) + real(tens(:,:) * pm, ki)
      call construct_polarization_tensor(eps6,eps6,tens)
      bornsc(6,:,:) = bornsc(6,:,:) + real(tens(:,:) * mp, ki)
      call construct_polarization_tensor(eps6,conjg(eps6),tens)
      bornsc(6,:,:) = bornsc(6,:,:) + real(tens(:,:) * mm, ki)
      !---#] particle 6 :
      if (include_helicity_avg_factor) then
         bornsc = bornsc / real(in_helicities, ki)
      end if
      if (include_color_avg_factor) then
         bornsc = bornsc / incolors
      end if
      if (include_symmetry_factor) then
         bornsc = bornsc / real(symmetry_factor, ki)
      end if
   end subroutine spin_correlated_lo2
   !---#] spin correlated ME :
   !---#[ construct polarisation tensor :
   pure subroutine construct_polarization_tensor(eps1, eps2, tens)
      implicit none
      complex(ki), dimension(0:3), intent(in) :: eps1, eps2
      complex(ki), dimension(0:3,0:3), intent(out) :: tens

      integer :: mu, nu

      do mu = 0,3
         do nu = 0, 3
            tens(mu,nu) = eps1(mu) * eps2(nu)
         end do
      end do
   end  subroutine construct_polarization_tensor
   !---#] construct polarisation tensor :
   pure function square_0l_0l_sc(color_vector1, color_vector2) result(amp)
      use p12_cbbar_hepneg_color, only: cmat => CC
      implicit none
      complex(ki), dimension(numcs), intent(in) :: color_vector1, color_vector2
      complex(ki) :: amp
      complex(ki), dimension(numcs) :: v1, v2

      v1 = matmul(cmat, color_vector2)
      v2 = conjg(color_vector1)
      amp = sum(v1(:) * v2(:))
   end function  square_0l_0l_sc

end module p12_cbbar_hepneg_matrix
