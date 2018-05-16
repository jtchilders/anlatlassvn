      subroutine setreal(p,rflav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'
      include 'cvecbos.h' 
      integer nleg
      parameter (nleg=nlegreal)
      real * 8 p(0:3,nleg),p0(0:3,nleg)
      integer rflav(nleg),rflav0(nleg)
      real * 8 amp2


      if(idvecbos.eq.24) then
         rflav0=rflav
         p0=p
      else
c Apply CP to the kinematics
         rflav0=-rflav
         p0=p
         p0(1,:)=-p(1,:)
      endif

      call compreal(p0,rflav0,amp2)

      end


      subroutine compreal(pin,rflav,realamp2)
      use consts_MCFM; use dpinitialization; use define_ampl  
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h' ! for alphas 
      include 'qcdcouple.f' 
      include 'facscale.f' 
      include 'PhysPars.h'
      integer nlegs
      parameter (nlegs=nlegreal)
      real * 8 pin(0:3,nlegs)
      real * 8 p(12,1:4) ! 12 = mxpart in MCFM 
      real * 8 p1(12,1:4) ! 12 = mxpart in MCFM 
      integer rflav(nlegs), rflav1(nlegs) 
      real * 8 realamp2 
      integer i, ud, cs,iq
      double precision msq(-5:5,-5:5)
      logical :: identical, set  
      logical, save :: firsttime = .true. 
      character(len=3) :: chn 
      include 'cvecbos.h'

      if (firsttime) then 
         call dpinitialize(7) 
c-----setting MCFMconsants
	 MCFMwmass  = ph_wmass
	 MCFMwwidth = ph_wwidth
	 MCFMzmass  = ph_zmass
	 MCFMzwidth = ph_zwidth
	 MCFMgw     = ph_unit_e/ph_sthw
         firsttime = .false. 
      endif
      Npoint = 7 
      MCFMgsq    = st_alpha*4d0*pi 
      ason2pi    = MCFMgsq/8d0/pi**2

      do i=1,nlegreal
         p(i,4)   = pin(0,i) 
         p(i,1:3) = pin(1:3,i) 
      enddo
      p(1,:) = -p(1,:) 
      p(2,:) = -p(2,:) 
      rflav1 = rflav

C     get incoming channel 
      if (rflav(1) > 0 .and. rflav(2) < 0) then 
         chn = 'qqb' 
      elseif (rflav(1) < 0 .and. rflav(2) > 0) then 
         chn = 'qbq' 
      elseif (rflav(1) > 0 .and. rflav(2) > 0) then 
         chn = 'qqq' 
      elseif (rflav(1) < 0 .and. rflav(2) < 0) then 
         chn = 'qbb' 
      elseif (rflav(1) > 0 .and. rflav(2) == 0) then 
         chn = 'qgl' 
      elseif (rflav(1) == 0 .and. rflav(2) > 0) then 
         chn = 'glq' 
      elseif (rflav(1) < 0 .and. rflav(2) == 0) then 
         chn = 'qbg' 
      elseif (rflav(1) == 0 .and. rflav(2) < 0) then 
         chn = 'gqb' 
      else
         write(*,*) 'rflav', rflav
         stop 'setreal: incoming flavour not OK' 
      endif

C     count number of families involved (whether quark pair are identical) 
      ud = 0; cs = 0 
      do i = 1,9 
         if (i < 3) then ! skip leptons in counting 
            if (rflav(i) == -1 .or. rflav(i) == 2) then 
               ud = ud +1
            elseif (rflav(i) == -3 .or. rflav(i) == 4) then 
               cs =cs +1
            endif
         elseif (i > 6) then 
            if (rflav(i) == 1 .or. rflav(i) == -2) then 
               ud = ud +1
            elseif (rflav(i) == 3 .or. rflav(i) == -4) then 
               cs =cs +1
            endif
         endif
      enddo
      if (ud == 2 .and. cs == 2) then 
         identical = .false. 
      elseif ((ud == 4 .and. cs == 0).or. (cs == 4 .and. ud == 0)) then 
         identical = .true. 
      else
         write(*,*) 'ud, cs', ud, cs 
         stop 'setreal: number of quark pairs not ok'
      endif

      p1 = p 

C     now make sure that routines are called with momenta in the right order 
      if (chn == 'qqb' .or. chn == 'qbq') then 
C     want 0 -> qb q l l l l qb q g 
         do i=7,9
            if (rflav(i) == 0) then 
               p1(9,:) = p(i,:) 
               rflav1(9) = rflav(i) 
            elseif (rflav(i) > 0) then 
               p1(8,:) = p(i,:) 
               rflav1(8) = rflav(i) 
            elseif (rflav(i) < 0) then 
               p1(7,:) = p(i,:) 
               rflav1(7) = rflav(i) 
            endif
         enddo

      elseif (chn == 'qqq' .or. chn == 'qbb') then 
C     want 0 -> q q' l l l l qb qb' g (1&7 same family, 2&8 same family)  
         iq = 0 
         do i=7,9
            if (rflav(i) == 0) then 
               p1(9,:) = p(i,:) 
               rflav1(9) = rflav(i) 
            elseif (rflav(i)+1 == rflav(1)) then 
               p1(7+iq,:) = p(i,:) 
               rflav1(7+iq) = rflav(i) 
               iq=iq+1
            else
               p1(8,:) = p(i,:) 
               rflav1(8) = rflav(i) 
            endif
         enddo


      elseif (chn == 'glq' .or. chn == 'qgl') then 
         if (chn == 'glq') then 
            iq = 2
         else
            iq = 1
         endif
C     want 0 -> q g l l l l qb' q' qb (1&7 same family, 2&8 same family)  
         set = .false. 
         do i=7,9
            if (rflav(i)+1 == rflav(iq) .and. .not. set) then 
               p1(9,:) = p(i,:) 
               rflav1(9) = rflav(i) 
               set = .true. 
            elseif (rflav(i) < 0 ) then 
               p1(7,:) = p(i,:) 
               rflav1(7) = rflav(i) 
            else
               p1(8,:) = p(i,:) 
               rflav1(8) = rflav(i) 
            endif
         enddo

      elseif (chn == 'gqb' .or. chn == 'qbg') then 
         if (chn == 'gqb') then 
            iq = 2
         else
            iq = 1
         endif
C     want 0 -> qb g l l l l qb' q' q (1&7 same family, 2&8 same family)  
         set = .false. 
         do i=7,9
            if (rflav(i)+1 == rflav(iq).and. .not.set) then 
               p1(9,:) = p(i,:) 
               rflav1(9) = rflav(i) 
               set = .true. 
            elseif (rflav(i) > 0 ) then 
               p1(8,:) = p(i,:) 
               rflav1(8) = rflav(i) 
            else
               p1(7,:) = p(i,:) 
               rflav1(7) = rflav(i) 
            endif
         enddo

      endif

!      if (abs(sum(p1)-sum(p)) > 1d-6) then 
!         write(*,*) 'rflav', rflav 
!         do i=1,9
!            write(*,*) 'p', p(i,:) 
!         enddo
!         do i=1,9
!            write(*,*) 'p1', p1(i,:) 
!         enddo
!         write(*,*) sum(p1), sum(p) 
!         stop 'setreal: p1 not correctly set' 
!      endif

C     now compute msq 
      call qqb_wpwp_qqb_g(p1,msq,chn,identical) 
      
      realamp2 = msq(rflav(1),rflav(2))
      realamp2 = realamp2 * vsymfact 

C     divide out ason2pi as this in included by MCFM file
C     but is included later on in PowHeg 
      realamp2 = realamp2/ason2pi  

      end
