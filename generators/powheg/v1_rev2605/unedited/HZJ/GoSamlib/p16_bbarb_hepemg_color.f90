module     p16_bbarb_hepemg_color
   ! file:      /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/Go
   ! Sam_POWHEG/Virtual/p16_bbarb_hepemg/common/color.f90
   ! generator: haggies (1.1)
   use p16_bbarb_hepemg_config, only: ki
   use p16_bbarb_hepemg_model, only: NC, Nf
   implicit none
   save

   private :: ki, NC, Nf

   real(ki), parameter :: TR = 0.5_ki

   complex(ki), parameter, private :: i_ = (0.0_ki, 1.0_ki)
   real(ki), parameter, private :: pi = &
   & 3.1415926535897932384626433832795028841971693993751058209749445920_ki
   real(ki), parameter, private :: pi6 = pi*pi/6.0_ki

   integer, parameter :: numcs = 1
   complex(ki), dimension(numcs, numcs) :: CC
   complex(ki), dimension(numcs, numcs) :: T1T1
   complex(ki), dimension(numcs, numcs) :: T1T2
   complex(ki), dimension(numcs, numcs) :: T1T6
   complex(ki), dimension(numcs, numcs) :: T2T2
   complex(ki), dimension(numcs, numcs) :: T2T6
   complex(ki), dimension(numcs, numcs) :: T6T6
   real(ki) :: incolors

   real(ki) :: CA, CF, KA, KF, gammaA, gammaF

   ! Basis vectors
   real(ki), dimension(numcs), parameter :: c1 = &
      & (/1.0_ki/)
contains
   subroutine     init_color()
      implicit none
      real(ki) :: NA
      real(ki) :: t1
      real(ki) :: t2
      real(ki) :: t3
      real(ki) :: t4
      
      t1 = NC*NC
      t2 = t1-1.0_ki
      NA = t2
      incolors = t1
      CC(1, 1) = (t2*TR)
      CC(1, 1) = (t2*TR)
      t3 = TR*TR
      t4 = (1.0_ki/NC+(t1-2.0_ki)*NC)*t3
      T1T1(1, 1) = t4
      T1T1(1, 1) = t4
      T1T2(1, 1) = ((NC-1.0_ki/NC)*t3)
      T1T2(1, 1) = ((NC-1.0_ki/NC)*t3)
      t1 = (1.0_ki-t1)*t3*NC
      T1T6(1, 1) = t1
      T1T6(1, 1) = t1
      T2T2(1, 1) = t4
      T2T2(1, 1) = t4
      T2T6(1, 1) = t1
      T2T6(1, 1) = t1
      T6T6(1, 1) = (t2*NC*TR)
      T6T6(1, 1) = (t2*NC*TR)

      CA = NC
      CF = TR * NA / NC
      ! KA = Kg in (C.11) [Catani,Seymour]
      KA = (67.0_ki/18.0_ki - pi6) * CA &
         & - 10.0_ki/9.0_ki * TR * Nf
      ! KF = Kq in (C.11) [Catani,Seymour]
      KF = (3.5_ki - pi6) * CF
      ! gammaA = \gamma_g in (C.11) [Catani,Seymour]
      gammaA = 11.0_ki/6.0_ki * CA - 2.0_ki/3.0_ki * TR * Nf
      ! gammaF = \gamma_q in (C.11) [Catani,Seymour]
      gammaF = 1.5_ki * CF
   end subroutine init_color
   subroutine     inspect_color(unit)
      implicit none
      integer, intent(in) :: unit
      integer :: i, j
      character :: ch1, ch2, ch3

      ch3 = ","
      write (unit,'(A13)') "golem_color=["
      do i=1,numcs
         do j=1,numcs
            if (j==1) then
               ch1 = "["
            else
               ch1 = " "
            endif

            if (j == numcs) then
               ch2 = "]"
               if (i == numcs) then
                  ch3 = "]"
               end if
            else
               ch2 = ","
            end if

            if (j == numcs) then
               write (unit,'(3x,A1,A8,G23.16,A1,G23.16,A1,A1,A1)') &
               & ch1, "complex(", real(CC(i,j)), ",", aimag(CC(i,j)), ")", &
               & ch2, ch3
            else
               write (unit,'(3x,A1,A8,G23.16,A1,G23.16,A1,A1)') &
               & ch1, "complex(", real(CC(i,j)), ",", aimag(CC(i,j)), ")", &
               & ch2
            end if
         enddo
      enddo
   end subroutine inspect_color
end module p16_bbarb_hepemg_color
