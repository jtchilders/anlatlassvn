module     p12_sbars_hepemg_model
   ! Model parameters for the model: /home/gionata/Documenti/Lavoro/GoSamPowheg/
   ! POWHEG-BOX/HZJ_tmp/GoSam_POWHEG/Virtual/model/model
   use p12_sbars_hepemg_config, only: ki, &
   & samurai_scalar, samurai_verbosity, samurai_test, &
   & samurai_group_numerators, samurai_istop, &
   & renormalisation, reduction_interoperation, deltaOS, &
   & nlo_prefactors
   implicit none

   private :: ki
   private :: samurai_scalar, samurai_verbosity, samurai_test
   private :: samurai_group_numerators, samurai_istop
   private :: renormalisation, reduction_interoperation, deltaOS
   private :: nlo_prefactors

   real(ki), parameter :: sqrt2 = &
      &1.414213562373095048801688724209698078&
      &5696718753769480731766797379_ki
   real(ki), parameter :: sqrt3 = &
      &1.732050807568877293527446341505872366&
      &9428052538103806280558069795_ki
   
   complex(ki) :: CVBC = (       0.041098933260000_ki,       -0.000000000000000_&
   &ki)
   complex(ki) :: CVBT = (       0.999155438842445_ki,       -0.000000000000000_&
   &ki)
   complex(ki) :: CVBU = (       0.001246715590975_ki,        0.003222990675929_&
   &ki)
   complex(ki) :: CVDC = (      -0.224561465788765_ki,        0.000132461478688_&
   &ki)
   complex(ki) :: CVDT = (       0.007988214712546_ki,        0.003222990675929_&
   &ki)
   complex(ki) :: CVDU = (       0.974436298851474_ki,       -0.000000000000000_&
   &ki)
   complex(ki) :: CVSC = (       0.973591737693919_ki,       -0.000000000000000_&
   &ki)
   complex(ki) :: CVST = (      -0.040341525833692_ki,        0.000724206004881_&
   &ki)
   complex(ki) :: CVSU = (       0.224700000000000_ki,       -0.000000000000000_&
   &ki)
   complex(ki) :: gauge6z = (       0.000000000000000_ki,        0.0000000000000&
   &00_ki)
   real(ki) :: GF =        0.000011663900000_ki
   real(ki), parameter :: gs =        1.000000000000000_ki
   real(ki), parameter :: mB =        0.000000000000000_ki
   real(ki), parameter :: mBMS =        0.000000000000000_ki
   real(ki), parameter :: mC =        0.000000000000000_ki
   real(ki), parameter :: mD =        0.000000000000000_ki
   real(ki), parameter :: me =        0.000000000000000_ki
   real(ki) :: mH =      114.400000000000006_ki
   real(ki), parameter :: mmu =        0.000000000000000_ki
   real(ki), parameter :: mS =        0.000000000000000_ki
   real(ki) :: mT =      171.199999999999989_ki
   real(ki), parameter :: mtau =        0.000000000000000_ki
   real(ki), parameter :: mU =        0.000000000000000_ki
   real(ki) :: mW =       80.418999999999997_ki
   real(ki) :: mZ =       91.188000000000002_ki
   real(ki) :: NC =        3.000000000000000_ki
   real(ki) :: Nf =        5.000000000000000_ki
   real(ki) :: Nfgen =        5.000000000000000_ki
   complex(ki) :: VCB = (       0.041098933260000_ki,        0.000000000000000_k&
   &i)
   complex(ki) :: VCD = (      -0.224561465788765_ki,       -0.000132461478688_k&
   &i)
   complex(ki) :: VCS = (       0.973591737693919_ki,        0.000000000000000_k&
   &i)
   complex(ki) :: VTB = (       0.999155438842445_ki,        0.000000000000000_k&
   &i)
   complex(ki) :: VTD = (       0.007988214712546_ki,       -0.003222990675929_k&
   &i)
   complex(ki) :: VTS = (      -0.040341525833692_ki,       -0.000724206004881_k&
   &i)
   complex(ki) :: VUB = (       0.001246715590975_ki,       -0.003222990675929_k&
   &i)
   complex(ki) :: VUD = (       0.974436298851474_ki,        0.000000000000000_k&
   &i)
   complex(ki) :: VUS = (       0.224700000000000_ki,        0.000000000000000_k&
   &i)
   real(ki), parameter :: wB =        0.000000000000000_ki
   real(ki) :: wchi =        0.000000000000000_ki
   real(ki) :: wghWm =        0.000000000000000_ki
   real(ki) :: wghWp =        0.000000000000000_ki
   real(ki) :: wghZ =        0.000000000000000_ki
   real(ki) :: wH =        0.000000000000000_ki
   real(ki) :: wphi =        0.000000000000000_ki
   real(ki), parameter :: wT =        0.000000000000000_ki
   real(ki) :: wtau =        0.000000000000000_ki
   real(ki) :: wW =        2.047600000000000_ki
   real(ki) :: wZ =        2.495200000000000_ki
   real(ki) :: gUa
   real(ki) :: gWWZZ
   real(ki) :: gBa
   real(ki) :: gtauv
   real(ki) :: gWWZ
   real(ki) :: gtaur
   real(ki) :: gUl
   real(ki) :: gBv
   real(ki) :: gUr
   real(ki) :: gBr
   real(ki) :: gUv
   real(ki) :: gHZZ
   real(ki) :: gSa
   real(ki) :: gtaua
   real(ki) :: gW
   real(ki) :: gHHH
   real(ki) :: gntaur
   real(ki) :: gntauv
   real(ki) :: gAPP
   real(ki) :: gntaul
   real(ki) :: gntaua
   complex(ki) :: gGWX
   real(ki) :: gBl
   real(ki) :: gH
   real(ki) :: gCa
   real(ki) :: ger
   real(ki) :: gev
   real(ki) :: gHC
   real(ki) :: gCl
   real(ki) :: gGWH
   real(ki) :: gea
   real(ki) :: gCr
   real(ki) :: gCv
   real(ki) :: gel
   real(ki) :: gtaul
   real(ki) :: gPWZ
   real(ki) :: gDl
   real(ki) :: gZZPP
   real(ki) :: gHHHH
   real(ki) :: gZZHH
   real(ki) :: gGZWP
   real(ki) :: gHXX
   real(ki) :: gWWPP
   real(ki) :: gDr
   real(ki) :: gDv
   real(ki) :: gGZH
   real(ki) :: gDa
   complex(ki) :: gWAPX
   real(ki) :: gHWW
   real(ki) :: gWWXX
   real(ki) :: gTa
   real(ki) :: gTl
   real(ki) :: gTr
   real(ki) :: gTv
   real(ki) :: gPtau
   real(ki) :: gHPP
   real(ki) :: gWWAZ
   real(ki) :: gWWAA
   real(ki) :: gnmul
   real(ki) :: gAAPP
   real(ki) :: gHtau
   real(ki) :: gZ
   real(ki) :: gXmu
   real(ki) :: gHHPP
   real(ki) :: Nfrat
   real(ki) :: cw
   real(ki) :: gGWZP
   real(ki) :: gnmua
   real(ki) :: NA
   real(ki) :: gner
   real(ki) :: gSv
   real(ki) :: gnev
   real(ki) :: gnmuv
   complex(ki) :: gZXH
   real(ki) :: gnmur
   real(ki) :: gnel
   real(ki) :: gHmu
   real(ki) :: gPPPP
   real(ki) :: gnea
   real(ki) :: gXtau
   real(ki) :: gWWHH
   real(ki) :: gWWWW
   real(ki) :: gWZPH
   real(ki) :: gHe
   complex(ki) :: gWZPX
   real(ki) :: gHHXX
   real(ki) :: gXD
   real(ki) :: gmur
   real(ki) :: gSr
   real(ki) :: gmuv
   real(ki) :: gXB
   real(ki) :: gXC
   real(ki) :: gXT
   real(ki) :: gXU
   real(ki) :: gmua
   real(ki) :: gWPH
   real(ki) :: gXS
   real(ki) :: gSl
   real(ki) :: gmul
   real(ki) :: gHT
   real(ki) :: gHU
   real(ki) :: gHS
   real(ki) :: gPWA
   real(ki) :: gHD
   real(ki) :: gZZXX
   real(ki) :: gHB
   real(ki) :: gPmu
   real(ki) :: gZPP
   real(ki) :: gXe
   real(ki) :: gPe
   real(ki) :: gWAPH
   real(ki) :: gXXXX
   real(ki) :: e
   real(ki) :: gPD
   real(ki) :: sw
   real(ki) :: gPB
   real(ki) :: gPC
   real(ki) :: gXXPP
   real(ki) :: gAZPP
   real(ki) :: gPT
   real(ki) :: gPU
   complex(ki) :: gWPX
   real(ki) :: gPS

   integer, parameter, private :: line_length = 80
   integer, parameter, private :: name_length = max(7,24)
   character(len=name_length), dimension(36) :: names = (/& 
      & "CVBC   ", &
      & "CVBT   ", &
      & "CVBU   ", &
      & "CVDC   ", &
      & "CVDT   ", &
      & "CVDU   ", &
      & "CVSC   ", &
      & "CVST   ", &
      & "CVSU   ", &
      & "gauge6z", &
      & "GF     ", &
      & "mH     ", &
      & "mT     ", &
      & "mW     ", &
      & "mZ     ", &
      & "NC     ", &
      & "Nf     ", &
      & "Nfgen  ", &
      & "VCB    ", &
      & "VCD    ", &
      & "VCS    ", &
      & "VTB    ", &
      & "VTD    ", &
      & "VTS    ", &
      & "VUB    ", &
      & "VUD    ", &
      & "VUS    ", &
      & "wchi   ", &
      & "wghWm  ", &
      & "wghWp  ", &
      & "wghZ   ", &
      & "wH     ", &
      & "wphi   ", &
      & "wtau   ", &
      & "wW     ", &
      & "wZ     "/)
   character(len=1), dimension(3) :: cc = (/'#', '!', ';'/)

   private :: digit, parsereal, names, cc

contains

   function     digit(ch, lnr) result(d)
      implicit none
      character(len=1), intent(in) :: ch
      integer, intent(in) :: lnr
      integer :: d
      d = -1
      select case(ch)
         case('0')
            d = 0
         case('1')
            d = 1
         case('2')
            d = 2
         case('3')
            d = 3
         case('4')
            d = 4
         case('5')
            d = 5
         case('6')
            d = 6
         case('7')
            d = 7
         case('8')
            d = 8
         case('9')
            d = 9
         case default
            write(*,'(A21,1x,I5)') 'invalid digit in line', lnr
         end select
   end function digit

   function     parsereal(str, ierr, lnr) result(num)
      implicit none
      character(len=*), intent(in) :: str
      integer, intent(out) :: ierr
      integer, intent(in) :: lnr
      integer, dimension(0:3,0:4), parameter :: DFA = &
      & reshape( (/1,  1,  2, -1,   & ! state = 0
      &            1, -1,  2,  3,   & ! state = 1
      &            2, -1, -1,  3,   & ! state = 2
      &            4,  4, -1, -1,   & ! state = 3
      &            4, -1, -1, -1/), (/4, 5/))
      real(ki) :: num
      integer :: i, expo, ofs, state, cclass, d, s1, s2
      num = 0.0_ki
      expo = 0
      state = 0
      ofs = 0
      s1 = 1
      s2 = 1
      d = -1
      cclass = -1
      do i=1,len(str)
         select case(str(i:i))
         case('_', '''', ' ')
            cycle
         case('+', '-')
            cclass = 1
         case('.')
            cclass = 2
         case('E', 'e')
            cclass = 3
         case default
            cclass = 0
            d = digit(str(i:i), lnr)
            if (d .eq. -1) then
               ierr = 1
               return
            end if
         end select
         if (cclass .eq. 0) then
            select case(state)
            case(0, 1)
               num = 10.0_ki * num + d
            case(2)
               num = 10.0_ki * num + d
               ofs = ofs - 1
            case(4)
               expo = 10 * expo + d
            end select
         elseif ((cclass .eq. 1) .and. (str(i:i) .eq. '-')) then
            if (state .eq. 0) then
               s1 = -1
            else
               s2 = -1
            endif
         end if
         ! Advance in the DFA
         state = DFA(cclass, state)
         if (state < 0) then
            write(*,'(A21,1x,A1,1x,A7,I5)') 'invalid position for', &
            & str(i:i), 'in line', lnr
            ierr = 1
            return
         end if
      end do
      num = s1 * num * 10.0_ki**(ofs + s2 * expo)
      ierr = 0
   end function parsereal

   subroutine     parseline(line,stat,line_number)
      implicit none
      character(len=*), intent(in) :: line
      integer, intent(out) :: stat
      integer, intent(in), optional :: line_number

      character(len=line_length) :: rvalue, ivalue, value
      character(len=name_length) :: name
      real(ki) :: re, im
      integer :: idx, icomma, idx1, idx2, lnr, nidx, ierr, pdg

      if(present(line_number)) then
         lnr = line_number
      else
         lnr = 0
      end if

      idx = scan(line, '=', .false.)
      if (idx .eq. 0) then
         if(present(line_number)) then
            write(*,'(A13,1x,I5)') 'error at line', line_number
         else
            write(*,'(A18)') 'error in parseline'
         end if
         stat = 1
         return
      end if
      name = adjustl(line(1:idx-1))
      value = adjustl(line(idx+1:len(line)))
      idx = scan(value, ',', .false.)

      if (name .eq. "renormalisation") then
         re = parsereal(value, ierr, lnr)
         if (ierr .ne. 0) then
            stat = 1
            return
         end if
         renormalisation = int(re)
      elseif (name .eq. "nlo_prefactors") then
         re = parsereal(value, ierr, lnr)
         if (ierr .ne. 0) then
            stat = 1
            return
         end if
         nlo_prefactors = int(re)
      elseif (name .eq. "deltaOS") then
         re = parsereal(value, ierr, lnr)
         if (ierr .ne. 0) then
            stat = 1
            return
         end if
         deltaOS = int(re)
      elseif (name .eq. "reduction_interoperation") then
         re = parsereal(value, ierr, lnr)
         if (ierr .ne. 0) then
            stat = 1
            return
         end if
         reduction_interoperation = int(re)
      elseif (name .eq. "samurai_scalar") then
         re = parsereal(value, ierr, lnr)
         if (ierr .ne. 0) then
            stat = 1
            return
         end if
         samurai_scalar = int(re)
      elseif (name .eq. "samurai_verbosity") then
         re = parsereal(value, ierr, lnr)
         if (ierr .ne. 0) then
            stat = 1
            return
         end if
         samurai_verbosity = int(re)
      elseif (name .eq. "samurai_test") then
         re = parsereal(value, ierr, lnr)
         if (ierr .ne. 0) then
            stat = 1
            return
         end if
         samurai_test = int(re)
      elseif (name .eq. "samurai_istop") then
         re = parsereal(value, ierr, lnr)
         if (ierr .ne. 0) then
            stat = 1
            return
         end if
         samurai_istop = int(re)
      elseif (name .eq. "samurai_group_numerators") then
         re = parsereal(value, ierr, lnr)
         if (ierr .ne. 0) then
            stat = 1
            return
         end if
         samurai_group_numerators = (int(re).ne.0)
      elseif (any(names .eq. name)) then
         do nidx=1,size(names)
            if (names(nidx) .eq. name) exit
         end do
         if (idx .gt. 0) then
            rvalue = value(1:idx-1)
            ivalue = value(idx+1:len(value))
            re = parsereal(rvalue, ierr, lnr)
            if (ierr .ne. 0) then
               stat = 1
               return
            end if
            im = parsereal(ivalue, ierr, lnr)
            if (ierr .ne. 0) then
               stat = 1
               return
            end if
         else
            re = parsereal(value, ierr, lnr)
            if (ierr .ne. 0) then
               stat = 1
               return
            end if
            im = 0.0_ki
         end if
         select case (nidx)
         case(1)
            CVBC = cmplx(re, im, ki)
         case(2)
            CVBT = cmplx(re, im, ki)
         case(3)
            CVBU = cmplx(re, im, ki)
         case(4)
            CVDC = cmplx(re, im, ki)
         case(5)
            CVDT = cmplx(re, im, ki)
         case(6)
            CVDU = cmplx(re, im, ki)
         case(7)
            CVSC = cmplx(re, im, ki)
         case(8)
            CVST = cmplx(re, im, ki)
         case(9)
            CVSU = cmplx(re, im, ki)
         case(10)
            gauge6z = cmplx(re, im, ki)
         case(11)
            GF = re
         case(12)
            mH = re
         case(13)
            mT = re
         case(14)
            mW = re
         case(15)
            mZ = re
         case(16)
            NC = re
         case(17)
            Nf = re
         case(18)
            Nfgen = re
         case(19)
            VCB = cmplx(re, im, ki)
         case(20)
            VCD = cmplx(re, im, ki)
         case(21)
            VCS = cmplx(re, im, ki)
         case(22)
            VTB = cmplx(re, im, ki)
         case(23)
            VTD = cmplx(re, im, ki)
         case(24)
            VTS = cmplx(re, im, ki)
         case(25)
            VUB = cmplx(re, im, ki)
         case(26)
            VUD = cmplx(re, im, ki)
         case(27)
            VUS = cmplx(re, im, ki)
         case(28)
            wchi = re
         case(29)
            wghWm = re
         case(30)
            wghWp = re
         case(31)
            wghZ = re
         case(32)
            wH = re
         case(33)
            wphi = re
         case(34)
            wtau = re
         case(35)
            wW = re
         case(36)
            wZ = re
         end select
      elseif (name(1:7).eq."masses(") then
         idx = scan(name, ')', .false.)
         if (idx.eq.0) then
            if(present(line_number)) then
               write(*,'(A13,1x,I5)') 'error at line', line_number
            else
               write(*,'(A18)') 'error in parseline'
            end if
            stat = 1
            return
         endif
         read(name(8:idx-1),*, iostat=ierr) pdg
         if (ierr.ne.0) then
            write(*,*) "Not an integer:", name(8:idx-1)
            if(present(line_number)) then
               write(*,'(A13,1x,I5)') 'error at line', line_number
            else
               write(*,'(A18)') 'error in parseline'
            end if
            stat = 1
            return
         end if
         select case(pdg)
            case(25)
               mH = parsereal(value, ierr, lnr)
            case(6)
               mT = parsereal(value, ierr, lnr)
            case default
               write(*,'(A20,1x,I10)') "Cannot set masses for code:", pdg
               stat = 1
               return
         end select
      elseif (name(1:6).eq."decay(") then
         idx = scan(name, ')', .false.)
         if (idx.eq.0) then
            if(present(line_number)) then
               write(*,'(A13,1x,I5)') 'error at line', line_number
            else
               write(*,'(A18)') 'error in parseline'
            end if
            stat = 1
            return
         endif
         read(name(7:idx-1),*, iostat=ierr) pdg
         if (ierr.ne.0) then
            write(*,*) "Not an integer:", name(7:idx-1)
            if(present(line_number)) then
               write(*,'(A13,1x,I5)') 'error at line', line_number
            else
               write(*,'(A18)') 'error in parseline'
            end if
            stat = 1
            return
         end if
         select case(pdg)
            case(9000003)
               wghZ = parsereal(value, ierr, lnr)
            case(9000004)
               wghWp = parsereal(value, ierr, lnr)
            case(9000005)
               wghWm = parsereal(value, ierr, lnr)
            case(15)
               wtau = parsereal(value, ierr, lnr)
            case(23)
               wZ = parsereal(value, ierr, lnr)
            case(24)
               wW = parsereal(value, ierr, lnr)
            case(25)
               wH = parsereal(value, ierr, lnr)
            case(250)
               wchi = parsereal(value, ierr, lnr)
            case(251)
               wphi = parsereal(value, ierr, lnr)
            case default
               write(*,'(A20,1x,I10)') "Cannot set decay for code:", pdg
               stat = 1
               return
         end select
      elseif (name(1:2).eq."m(" .or. name(1:2).eq."w(") then
         idx = scan(name, ')', .false.)
         if (idx.eq.0) then
            if(present(line_number)) then
               write(*,'(A13,1x,I5)') 'error at line', line_number
            else
               write(*,'(A18)') 'error in parseline'
            end if
            stat = 1
            return
         endif
         read(name(3:idx-1),*, iostat=ierr) pdg
         if (ierr.ne.0) then
            write(*,*) "pdg is not an integer:", name(3:idx-1)
            if(present(line_number)) then
               write(*,'(A13,1x,I5)') 'error at line', line_number
            else
               write(*,'(A18)') 'error in parseline'
            end if
            stat = 1
            return
         end if
         if (name(1:1).eq."m") then
            ! set mass according to PDG code
            select case(pdg)
            case(25)
               mH = parsereal(value, ierr, lnr)
            case(6)
               mT = parsereal(value, ierr, lnr)
            case default
               write(*,'(A20,1x,I10)') "Cannot set mass for PDG code:", pdg
               stat = 1
               return
            end select
         else
            ! set width according to PDG code
            select case(pdg)
            case(9000003)
               wghZ = parsereal(value, ierr, lnr)
            case(9000004)
               wghWp = parsereal(value, ierr, lnr)
            case(9000005)
               wghWm = parsereal(value, ierr, lnr)
            case(15)
               wtau = parsereal(value, ierr, lnr)
            case(23)
               wZ = parsereal(value, ierr, lnr)
            case(24)
               wW = parsereal(value, ierr, lnr)
            case(25)
               wH = parsereal(value, ierr, lnr)
            case(250)
               wchi = parsereal(value, ierr, lnr)
            case(251)
               wphi = parsereal(value, ierr, lnr)
            case default
               write(*,'(A20,1x,I10)') "Cannot set width for PDG code:", pdg
               stat = 1
               return
            end select
         endif
      else
         write(*,'(A20,1x,A20)') 'Unrecognized option:', name
         stat = 1
         return
      end if
      stat = 0
   end subroutine parseline

   subroutine     parse(aunit)
      implicit none
      integer, intent(in) :: aunit
      character(len=line_length) :: line
      integer :: ios, lnr
      lnr = 0
      loop1: do
         read(unit=aunit,fmt='(A80)',iostat=ios) line
         if(ios .ne. 0) exit
         lnr = lnr + 1
         line = adjustl(line)
         if (line .eq. '') cycle loop1
         if (any(cc .eq. line(1:1))) cycle loop1

         call parseline(line,ios,lnr)
         if(ios .ne. 0) then
            write(*,'(A44,I2,A1)') &
            & 'Error while reading parameter file in parse(', aunit, ')'
         end if
      end do loop1
   end subroutine parse
!---#[ SLHA READER:
   subroutine     read_slha(ch, ierr)
      implicit none
      integer, intent(in) :: ch
      integer, intent(out), optional :: ierr

      integer :: lnr, i, l, ofs, ios
      character(len=255) :: line

      integer :: block

      ofs = iachar('A') - iachar('a')

      lnr = 0
      loop1: do
         read(unit=ch,fmt='(A80)',iostat=ios) line
         if(ios .ne. 0) exit
         lnr = lnr + 1

         i = scan(line, '#', .false.)
         if (i .eq. 0) then
            l = len_trim(line)
         else
            l = i - 1
         end if

         if (l .eq. 0) cycle loop1

         ucase: do i = 1, l
            if (line(i:i) >= 'a' .and. line(i:i) <= 'z') then
               line(i:i) = achar(iachar(line(i:i))+ofs)
            end if
         end do ucase

         if (line(1:1) .eq. 'B') then
            if (line(1:5) .eq. 'BLOCK') then
               line = adjustl(line(6:l))
               do i=1,l
                 if (line(i:i) <= ' ') exit
               end do
               l = i
               if ("MASSES" .eq. line(1:l)) then
                  block = 0
               elseif ("DECAY" .eq. line(1:l)) then
                  block = 1
               else
                  block = -1
               end if
            else
               write(*,'(A37,I5)') "Illegal statement in SLHA file, line ", lnr
               if (present(ierr)) ierr = 1
               return
            end if
         elseif (line(1:1) .eq. 'D') then
            if (line(1:5) .eq. 'DECAY') then
               line = adjustl(line(6:l))
               call read_slha_line_decay(line, i)
               block = 2
            else
               write(*,'(A37,I5)') "Illegal statement in SLHA file, line ", lnr
               if (present(ierr)) ierr = 1
               return
            end if
         else
            ! read a parameter line
            select case(block)
            case(0)
               call read_slha_block_masses(line(1:l), i)
               if (i .ne. 0) then
                  if (present(ierr)) ierr = 1
                  write(*,'(A44,I5)') &
                  & "Unrecognized line format in SLHA file, line ", lnr
                  return
               end if
            case(1)
               call read_slha_block_decay(line(1:l), i)
               if (i .ne. 0) then
                  if (present(ierr)) ierr = 1
                  write(*,'(A44,I5)') &
                  & "Unrecognized line format in SLHA file, line ", lnr
                  return
               end if
            case default
               cycle loop1
            end select
         end if
      end do loop1
      if (present(ierr)) ierr = 0
   end subroutine read_slha
   subroutine read_slha_block_masses(line, ierr)
      implicit none
      character(len=*), intent(in) :: line
      integer, intent(out), optional :: ierr
      integer :: idx1,ioerr
      real(ki) :: value

      read(line,*,iostat=ioerr) idx1, value
      if (ioerr .ne. 0) then
         if (present(ierr)) ierr = 1
         return
      end if
      select case(idx1)
      case(25)
         mH = value
      case(6)
         mT = value
      end select
      if (present(ierr)) ierr = 0
   end subroutine read_slha_block_masses
   subroutine read_slha_block_decay(line, ierr)
   !  This subroutine reads the 'branching ratios' of
   !  the slha file: these are just thrown away
      implicit none
      character(len=*), intent(in) :: line
      integer, intent(out), optional :: ierr
      integer :: idx1,idx2,ioerr,nda
      real(ki) :: value
      read(line,*,iostat=ioerr) value, nda, idx1, idx2
      if (ioerr .ne. 0) then
         if (present(ierr)) ierr = 1
         return
      end if
      if (present(ierr)) ierr = 0
   end subroutine read_slha_block_decay
   subroutine read_slha_line_decay(line, ierr)
      implicit none
      character(len=*), intent(in) :: line
      integer, intent(out), optional :: ierr
      integer :: idx1,ioerr
      real(ki) :: value

      read(line,*,iostat=ioerr) idx1, value
      if (ioerr .ne. 0) then
         if (present(ierr)) ierr = 1
         return
      end if
      select case(idx1)
      case(9000003)
         wghZ = value
      case(9000004)
         wghWp = value
      case(9000005)
         wghWm = value
      case(15)
         wtau = value
      case(23)
         wZ = value
      case(24)
         wW = value
      case(25)
         wH = value
      case(250)
         wchi = value
      case(251)
         wphi = value
      end select
      if (present(ierr)) ierr = 0
   end subroutine read_slha_line_decay
!---#] SLHA READER:
!---#[ subroutine init_functions:
   subroutine     init_functions()
      implicit none
      complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
      real(ki), parameter :: pi = 3.14159265358979323846264&
     &3383279502884197169399375105820974944592307816406286209_ki
      complex(ki) :: reg1
      real(ki) :: reg2
      real(ki) :: reg3
      real(ki) :: reg4
      real(ki) :: reg5
      real(ki) :: reg6
      real(ki) :: reg7
      real(ki) :: reg8
      real(ki) :: reg9
      real(ki) :: reg10
      real(ki) :: reg11
      real(ki) :: reg12
      real(ki) :: reg13
      real(ki) :: reg14
      real(ki) :: reg15
      real(ki) :: reg16
      real(ki) :: reg17
      real(ki) :: reg18
      real(ki) :: reg19
!     I've added the following two lines. Carlo Oleari
      real(ki) :: gerpwg,gelpwg
      common/Zlepcoupl/gerpwg,gelpwg
      gAPP = (-1.0_ki)
      gWWAA = (-1.0_ki)
      gAAPP = (2.0_ki)
      Nfrat = (ifpos(Nfgen, (Nf/Nfgen), (1.0_ki)))
      cw = (mW/mZ)
      NA = (NC*NC-(1.0_ki))
      gWZPH = (-1.0_ki/2.0_ki)/cw
      gWZPX = (-1.0_ki/2.0_ki)/cw*i_
      gPWA = -(mW)
      sw = (sqrt(((1.0_ki)-mW*mW/(mZ*mZ))))
      reg2 = cw*cw
      reg3 = sw*sw
      gWWZZ = (-1.0_ki)*reg2/reg3
      reg4 = cw/sw
      gWWZ = -(reg4)
      reg5 = mW/(reg2*sw)
      gHZZ = reg5
      gW = (1.0_ki)/(sqrt2*sw)
      reg6 = mH*mH
      reg7 = mW*sw
      reg8 = reg6/reg7
      gHHH = (-3.0_ki/2.0_ki)*reg8
      gGWX = (-1.0_ki/2.0_ki)*i_*mW/sw
      gH = ((1.0_ki/24.0_ki)/(reg7*pi*pi))
      reg9 = mW*sqrt2*sw
      gPmu = (mmu/reg9)
      reg10 = mW/sw
      gGWH = (-1.0_ki/2.0_ki)*reg10
      reg11 = (1.0_ki/2.0_ki)/reg3
      gWWXX = reg11
      reg12 = cw*sw
      reg13 = (reg2-reg3)/reg12
      gZZPP = reg13
      reg6 = reg6/(reg7*reg7)
      reg14 = (1.0_ki/4.0_ki)*reg6
      gHHPP = -(reg14)
      reg15 = (1.0_ki/2.0_ki)/(reg12*reg12)
      gZZHH = reg15
      gGZWP = ((1.0_ki/2.0_ki)/reg12*mW)
      reg8 = (1.0_ki/2.0_ki)*reg8
      gHXX = -(reg8)
      gWWPP = reg11
      gGZH = (-1.0_ki/2.0_ki)*reg5
      reg1 = (1.0_ki/2.0_ki)*i_/sw
      gWAPX = -(reg1)
      gHWW = reg10
      gPtau = (mtau/reg9)
      gHPP = -(reg8)
      reg5 = (1.0_ki/2.0_ki)/reg7*mBMS
      gHB = -(reg5)
      gWWAZ = reg4
      reg4 = (1.0_ki/2.0_ki)/reg7*mtau
      gHtau = -(reg4)
      gZ = (1.0_ki)/reg12
      reg8 = (1.0_ki/4.0_ki)*gZ
      gUa = reg8
      gBa = -(reg8)
      reg10 = ((1.0_ki/4.0_ki)-reg3)*gZ
      gtauv = -(reg10)
      reg16 = ((1.0_ki/4.0_ki)-(2.0_ki/3.0_ki)*reg3)*gZ
      gUv = reg16
      gUl = (gUa+gUv)
      gUr = (gUv-gUa)
      gtaua = -(reg8)
      gtaur = (gtauv-gtaua)
      gtaul = (gtaua+gtauv)
      gntauv = reg8
      gntaua = reg8
      gntaur = (gntauv-gntaua)
      gntaul = (gntaua+gntauv)
      gCa = reg8
      gev = -(reg10)
      gea = -(reg8)
!     I've commented the following line and added the next one. Carlo Oleari
!      ger = (gev-gea)
      ger = gerpwg
      gCv = reg16
      gCl = (gCa+gCv)
      gCr = (gCv-gCa)
!     I've commented the following line and added the next one. Carlo Oleari
!      gel = (gea+gev)
      gel = gelpwg
      reg17 = ((1.0_ki/4.0_ki)-(1.0_ki/3.0_ki)*reg3)*gZ
      gBv = -(reg17)
      gBr = (gBv-gBa)
      gBl = (gBa+gBv)
      gDv = -(reg17)
      gDa = -(reg8)
      gDr = (gDv-gDa)
      gDl = (gDa+gDv)
      gTa = reg8
      gTv = reg16
      gTl = (gTa+gTv)
      gTr = (gTv-gTa)
      reg16 = (1.0_ki/2.0_ki)/reg7*mmu
      gXmu = reg16
      reg18 = (3.0_ki/4.0_ki)*reg6
      gHHHH = -(reg18)
      reg2 = -(reg3+reg2)
      gGWZP = ((1.0_ki/2.0_ki)/reg12*reg2*mW)
      gnmua = reg8
      reg19 = (1.0_ki/2.0_ki)/reg7*mC
      gXC = -(reg19)
      gnev = reg8
      gnmuv = reg8
      gnmul = (gnmua+gnmuv)
      gZXH = (-1.0_ki/2.0_ki)/reg12*i_
      gnmur = (gnmuv-gnmua)
      gHmu = -(reg16)
      reg16 = (1.0_ki/2.0_ki)/reg7*mT
      gXT = -(reg16)
      gPPPP = (-1.0_ki/2.0_ki)*reg6
      gnea = reg8
      gner = (gnev-gnea)
      gnel = (gnea+gnev)
      gXtau = reg4
      gWWHH = reg11
      gWWWW = (1.0_ki)/reg3
      reg3 = (1.0_ki/2.0_ki)/reg7*mU
      gXU = -(reg3)
      gHHXX = -(reg14)
      reg4 = (1.0_ki/2.0_ki)/reg7*mD
      gXD = reg4
      gmuv = -(reg10)
      gXB = reg5
      gSv = -(reg17)
      gSa = -(reg8)
      gSr = (gSv-gSa)
      reg5 = (1.0_ki/2.0_ki)/reg7*me
      gHe = -(reg5)
      gmua = -(reg8)
      gmur = (gmuv-gmua)
      gPB = (mBMS/reg9)
      reg6 = (1.0_ki/2.0_ki)/reg7*mS
      gXS = reg6
      gSl = (gSa+gSv)
      gmul = (gmua+gmuv)
      gHT = -(reg16)
      gHU = -(reg3)
      gHS = -(reg6)
      gHD = -(reg4)
      gZZXX = reg15
      gPWZ = (-1.0_ki)*reg7/cw
      gHC = -(reg19)
      gZPP = ((1.0_ki/2.0_ki)*reg13)
      gXe = reg5
      gPe = (me/reg9)
      reg3 = (1.0_ki/2.0_ki)/sw
      gWAPH = -(reg3)
      gXXXX = -(reg18)
      e = (sqrt(((8.0_ki)/(sqrt((2.0_ki)))*GF))*mW*sw)
      gPD = (mD/reg9)
      gWPH = -(reg3)
      gPC = (mC/reg9)
      gXXPP = -(reg14)
      gAZPP = (reg2/reg12)
      gPT = (mT/reg9)
      gPU = (mU/reg9)
      gWPX = -(reg1)
      gPS = (mS/reg9)
   end subroutine init_functions
!---#] subroutine init_functions:
!---#[ utility functions for model initialization:
   pure function ifpos(x0, x1, x2)
      implicit none
      real(ki), intent(in) :: x0, x1, x2
      real(ki) :: ifpos

      if (x0 > 0.0_ki) then
         ifpos = x1
      else
         ifpos = x2
      endif
   end  function ifpos

   pure function sort4(m1, m2, m3, m4, n)
      implicit none
      real(ki), intent(in) :: m1, m2, m3, m4
      integer, intent(in) :: n
      real(ki) :: sort4

      real(ki), dimension(4) :: m
      logical :: f
      integer :: i
      real(ki) :: tmp

      m(1) = m1
      m(2) = m2
      m(3) = m3
      m(4) = m4

      ! Bubble Sort
      do
         f = .false.

         do i=1,3
            if (abs(m(i)) .gt. abs(m(i+1))) then
               tmp = m(i)
               m(i) = m(i+1)
               m(i+1) = tmp
               f = .true.
            end if
         end do

         if (.not. f) exit
      end do

      sort4 = m(n)
   end  function sort4
!---#] utility functions for model initialization:
end module p12_sbars_hepemg_model
