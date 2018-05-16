c     returns 2 Re(M_B * M_V)/(as/(2pi)), 
c     where M_B is the Born amplitude and 
c     M_V is the finite part of the virtual amplitude
c     The as/(2pi) factor is attached at a later point
      subroutine setvirtual(pin,vflav,virtual)
      use consts_MCFM; use dpinitialization; use define_ampl  
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'epinv.f' 
      include 'epinv2.f' 
      include 'scale.f' 
      include 'qcdcouple.f' 
      include 'facscale.f' 
      include 'PhysPars.h'
      real * 8 pin(0:3,nlegborn),pin0(0:3,nlegborn)
      integer vflav(nlegborn),vflav0(nlegborn)
      real * 8 virtual
      real * 8 born
      real *8 dotp, powheginput
      external dotp, powheginput 
C     --------------------------
      character chn*3
      integer mxpart,i,ud, cs
      parameter(mxpart=12) 
      double precision p(mxpart,4),msqB(-5:5,-5:5),msq(-5:5,-5:5)
      double precision p1(mxpart,4)
      logical :: identical, iswap  
      logical, save :: firsttime = .true. 
      logical, save :: polesonly = .false. 
      integer fakevirt
      integer, save :: countampl  = 0 
      save fakevirt 
      include 'cvecbos.h'


      if (firsttime) then 
         call dpinitialize(6)
c-----setting MCFMconstants
	 MCFMwmass  = ph_wmass
	 MCFMwwidth = ph_wwidth
	 MCFMzmass  = ph_zmass
	 MCFMzwidth = ph_zwidth
	 MCFMgw     = ph_unit_e/ph_sthw
      
c Paolo
         fakevirt=powheginput("#fakevirt")
         if (fakevirt == 1) write(*,*) 'WARNING: Using fakevirt !'
         if (polesonly) then 
            write(*,*) 'Scale dependent part only' 
         else
            write(*,*) 'Full virtual' 
         endif
         firsttime  = .false. 
      endif

      if(idvecbos.eq.24) then
         vflav0=vflav
         pin0=pin
      else
c Apply CP to the kinematics
         vflav0=-vflav
         pin0=pin
         pin0(1,:)=-pin(1,:)
      endif
      MCFMgsq    = st_alpha*4d0*pi 
      ason2pi    = MCFMgsq/8d0/pi**2
      facscale=sqrt(st_mufact2)
      scale=sqrt(st_muren2)

      Npoint = 6 
      epinv = 0d0
      epinv2 = 0d0 

      do i=1,nlegborn
         p(i,4)   = pin0(0,i) 
         p(i,1:3) = pin0(1:3,i) 
      enddo
      p(1,:) = -p(1,:)
      p(2,:) = -p(2,:)

      if (vflav0(1) < 0 .and. vflav0(2) < 0) then 
         chn = 'qbb' 
      elseif (vflav0(1) < 0 .and. vflav0(2) > 0) then 
         chn = 'qbq' 
      elseif (vflav0(1) > 0 .and. vflav0(2) < 0) then 
         chn = 'qqb' 
      elseif (vflav0(1) > 0 .and. vflav0(2) > 0) then 
         chn = 'qqq' 
      else
         write(*,*) 'vflav0', vflav0 
         stop 'setvirtual: undefined channel' 
      endif

C     count number of families involved (whether quark pairs are identical) 
      ud = 0;  cs = 0 
      do i = 1,8 
         if (i < 3) then ! skip leptons in counting 
            if (vflav0(i) == -1 .or. vflav0(i) == 2) then 
               ud = ud +1
            elseif (vflav0(i) == -3 .or. vflav0(i) == 4) then 
               cs =cs +1
            endif
         elseif (i > 6) then 
            if (vflav0(i) == 1 .or. vflav0(i) == -2) then 
               ud = ud +1
            elseif (vflav0(i) == 3 .or. vflav0(i) == -4) then 
               cs =cs +1
            endif
         endif
      enddo
      if (ud == 2 .and. cs == 2) then 
         identical = .false. 
      elseif (ud == 4 .or. cs == 4) then 
         identical = .true. 
      else
         write(*,*) 'ud, cs', ud, cs 
         stop 'setborn: number of quark pairs not ok'
      endif
      p1 = p 

C     now make sure that routines are called with momenta in the right order 
      iswap = .false. 
      if (chn == 'qbb') then 
         if (vflav0(1) - 1 .ne. vflav0(7)) then 
            p1(8,:) = p(7,:) 
            p1(7,:) = p(8,:) 
            iswap = .true. 
         endif
      elseif (chn == 'qqq') then 
         if (vflav0(1) - 1 .ne. vflav0(7)) then 
            p1(8,:) = p(7,:) 
            p1(7,:) = p(8,:) 
            iswap = .true. 
         endif
      endif
      call qqb_wpwp_qqb(p1,msqB,chn,identical)
      born = msqB(vflav0(1),vflav0(2))

      if(fakevirt.ne.1) then
C     -- now compute virtual 
         call qqb_wpwp_qqb_v(p1,msq,chn,polesonly,identical)
C     divide out ason2pi as this in included by MCFM file, 
C     but is included later on in PowHeg 
         virtual = msq(vflav0(1),vflav0(2))/ason2pi  
C     scheme change from DRED and coupling constsnt fix to go to MSbar 
C     see p.10 of 1002.2581
         virtual = virtual + born*(Nc/6d0*2d0-4d0*(cf/2d0))
      else
         virtual = 0.2d0*born     ! -- FAKE constant K factor 
      endif
      virtual = virtual * vsymfact 

      countampl = countampl +1 
c      if ((countampl/100000)*100000 == countampl) 
c     1  write(*,*) 'Done', counampl, 'virtual squared amplitudes'  
      end


