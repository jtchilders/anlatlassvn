c     returns 2 Re(M_B * M_V)/(as/(2pi)), 
c     where M_B is the Born amplitude and 
c     M_V is the finite part of the virtual amplitude
c     The as/(2pi) factor is attached at a later point
c     Use MCFM subroutines
      subroutine setvirtual(p,vflav,virtual)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'
c     MCFM include
      include 'MCFM_include/scale.f'
      include 'MCFM_include/qcdcouple.f'
      include 'MCFM_include/constants.f'
      include 'MCFM_include/stopscales.f'
      include 'MCFM_include/masses.f'
cccccccccccccccccccccccccccccccccccccc
c     only to check with madloop paper
      include 'MCFM_include/ewcharge.f'
      include 'MCFM_include/ewcouple.f'
      double precision aemmz
      double precision LOqgmadloop,c0qgmadloop
cccccccccccccccccccccccccccccccccccccc
      double precision p(0:3,nlegborn)
      integer vflav(nlegborn)
      double precision virtual
      double precision pmcfm(mxpart,1:4),pold(0:3,nlegborn)
      double precision virt(-nf:nf, -nf:nf)
      save virt,pold
      logical ini
      save ini
      data ini/.true./
      logical equalp
      integer ileg,imu
      double precision kb_mad(0:3,nlegborn),amp2mad
      double precision LOqg,LOgq,LOqbarg,LOgqbar
      double precision CKM(1:6,1:6)
      common/cckm/CKM
      double precision ckm_light,ckm_tb,ckmfact

      logical checkborn
      parameter (checkborn=.false.)
      logical madloopcheck
      parameter (madloopcheck=.false.)
        
      if (ini) then
         call virtual_initialize_MCFM
         do imu=0,3
            do ileg=1,nlegborn
               pold(imu,ileg)=0d0
            enddo
         enddo
         ini=.false.
      endif
ccccccccccccccccccccccccccccccccccccc
c     !ER:
c     relevant for st 4 flavour
      scale = sqrt(st_muren2)
      renscale_H=sqrt(st_muren2)
      renscale_L=sqrt(st_muren2)
      facscale_H=sqrt(st_mufact2)
      facscale_L=sqrt(st_mufact2)

      as_H = st_alpha
      as_L = st_alpha
      mt=topmass_pow
      mb=bmass_pow
cccccccccccccccccccccccccccccccccccccc

c     MCFM fills an array with all processes.
ccccccccccccccccccccccccccccccccccccccccccc
c     Carlo trick to avoid recomputation of amplitudes already computed
      equalp=.true.
      do imu=0,3
         do ileg=1,nlegborn
            equalp=equalp.and.(p(imu,ileg).eq.pold(imu,ileg))
            pold(imu,ileg)=p(imu,ileg)
         enddo
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c     !ER: uncomment here to recompute always virtual amplitudes
c$$$      equalp=.false.
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(.not.equalp) then
*     q(p1) + g(p2) -> t(p3) + b(p4) + q'(p5)  
         if(ttype*vflav(3).ne.6.or.ttype*vflav(4).ne.-5) then
            write(*,*) 'flavour problem in setvirtual'
            call exit(-1)
         endif
         call mom_to_MCFM(p,pmcfm)
cccccccccccccccccccccccccccccccccccccccccccc
         if(madloopcheck) then
c     Kinematic point in 1103.0621, app. A.2.8
c     To check it, look lines with the string '!A28' in this file
c         pu = (250, 0, 0, 250)
c         pg = (250, 0, 0, -250)
c         pt = (255.3729192644455, 29.17335589243201, 159.7715722928748, -91.96084974966891)
c         pb = (177.5248259329844, -66.11648748945143, -111.8173550700313, 120.9144450003231)
c         pd = (67.10225480257022, 36.94313159701945, -47.95421722284348, -28.95359525065417)

            pmcfm(1,4)= -250.d0
            pmcfm(1,1)= -0.d0
            pmcfm(1,2)= -0.d0
            pmcfm(1,3)= -250.d0

            pmcfm(2,4)= -250.d0
            pmcfm(2,1)= -0.d0  
            pmcfm(2,2)= -0.d0   
            pmcfm(2,3)= +250.d0 
            
            pmcfm(3,4)= 255.3729192644455d0 
            pmcfm(3,1)= 29.17335589243201d0 
            pmcfm(3,2)= 159.7715722928748d0 
            pmcfm(3,3)= -91.96084974966891d0 
            
            pmcfm(4,4)= 177.5248259329844d0 
            pmcfm(4,1)= -66.11648748945143d0 
            pmcfm(4,2)= -111.8173550700313d0 
            pmcfm(4,3)= 120.9144450003231d0 
            
            pmcfm(5,4)= 67.10225480257022d0 
            pmcfm(5,1)= 36.94313159701945d0 
            pmcfm(5,2)= -47.95421722284348d0 
            pmcfm(5,3)= -28.95359525065417d0 
            
            wmass= 80.419d0
            wwidth=2.0476d0     ! I need this width to check with madloop paper (!A28)
            xw = 0.2222465331d0 ! sin^2(theta_W)
            aemmz = 1./132.50698d0 ! alphaem 
            
            mt= 174.3d0
            mb= 4.5d0


            st_muren2=(91.188d0)**2
            scale = sqrt(st_muren2)
            renscale_H=scale
            renscale_L=scale
            facscale_H=scale
            facscale_L=scale
            
c--- Now set up the other derived parameters
            gwsq=fourpi*aemmz/xw
            esq=gwsq*xw
         endif
cccccccccccccccccccccccccccccccccccccccccccc

         call qg_tbq_v(pmcfm,virt,LOqg,LOgq,LOqbarg,LOgqbar)
         if(madloopcheck) then
            LOqgmadloop=7.79629086614075984E-007
            c0qgmadloop=2.31517097632403642E-007
            write(*,*) 'Remember to change alfas to 0.118d0'
            write(*,*) 'Madloop point, Born ratio          : ',
     $           LOqg / LOqgmadloop
            write(*,*) 'Madloop point, virtual finite ratio: ',
     $           virt(2,0) / c0qgmadloop
            stop
         endif
      endif
cccccccccccccccccccccccccccccccccccccccccc
c     CKM
      ckm_tb=CKM(abs(vflav(3)),abs(vflav(4)))**2
      if(vflav(1).eq.0.and.vflav(2).ne.0) then
         ckm_light=CKM(abs(vflav(2)),abs(vflav(5)))**2
      elseif(vflav(1).ne.0.and.vflav(2).eq.0) then
         ckm_light=CKM(abs(vflav(1)),abs(vflav(5)))**2
      else
         write(*,*) 'setvirtual: flavour not valid ',vflav(1),vflav(2)
         call exit(-1)
      endif
      ckmfact=ckm_tb*ckm_light
      virtual=virt(vflav(1),vflav(2))   *   ckmfact
cccccccccccccccccccccccccccc
c     check MCFM Born against Madgraph
      if(checkborn) then
         do ileg=1,nlegborn
            do imu=0,3
               kb_mad(imu,ileg)=p(imu,ileg)
            enddo
         enddo
c     to avoid bugs in HELAS, restore exact masslessness of  incoming partons 
         kb_mad(0,1)=dabs(kb_mad(3,1))
         kb_mad(0,2)=dabs(kb_mad(3,2))
         if(ttype.eq.1) then
            call compborn(kb_mad,vflav,amp2mad)
         elseif(ttype.eq.-1) then
            call compborn_tb(kb_mad,vflav,amp2mad)
         else
            write(*,*) 'wrong ttype in virtual.f, while checking Born'
            call exit(-1)
         endif

         if(vflav(1).eq.0.and.vflav(2).gt.0) then
            ckm_light=CKM(abs(vflav(2)),abs(vflav(5)))**2
            write(*,*) 'GQ ',LOgq    /amp2mad *ckm_tb*ckm_light
         elseif(vflav(1).eq.0.and.vflav(2).lt.0) then
            ckm_light=CKM(abs(vflav(2)),abs(vflav(5)))**2
            write(*,*) 'GQB',LOgqbar /amp2mad *ckm_tb*ckm_light
         elseif(vflav(2).eq.0.and.vflav(1).gt.0) then
            ckm_light=CKM(abs(vflav(1)),abs(vflav(5)))**2
            write(*,*) 'QG ',LOqg    /amp2mad *ckm_tb*ckm_light
         elseif(vflav(2).eq.0.and.vflav(1).lt.0) then
            ckm_light=CKM(abs(vflav(1)),abs(vflav(5)))**2
            write(*,*) 'QBG',LOqbarg /amp2mad *ckm_tb*ckm_light
         else
            write(*,*) 'setvirtual: error while checking LO amplitudes'
            call exit(-1)
         endif
      endif
cccccccccccccccccccccccccccccccccc
     
c     Now divide result by as/(2pi)
      virtual=virtual/(st_alpha/(2*pi))
      end




      subroutine virtual_initialize_MCFM
      implicit none
      include 'MCFM_include/constants.f'
      include 'MCFM_include/ewcharge.f'
      include 'MCFM_include/ewcouple.f'
      include 'MCFM_include/masses.f'
      include 'MCFM_include/scale.f'
      include 'MCFM_include/qcdcouple.f'
c      include 'MCFM_include/zcouple.f'
      include 'MCFM_include/epinv.f'
      include 'MCFM_include/epinv2.f'
      include 'MCFM_include/nflav.f'
      include 'MCFM_include/nwz.f'
      include 'PhysPars.h'
      include 'MCFM_include/b0.f'
c      block data wsalam1
      data Q(-5)/+0.333333333333333d0/
      data Q(-4)/-0.666666666666667d0/
      data Q(-3)/+0.333333333333333d0/
      data Q(-2)/-0.666666666666667d0/
      data Q(-1)/+0.333333333333333d0/
      data Q(0)/+0d0/
      data Q(+1)/-0.333333333333333d0/
      data Q(+2)/+0.666666666666667d0/
      data Q(+3)/-0.333333333333333d0/
      data Q(+4)/+0.666666666666667d0/
      data Q(+5)/-0.333333333333333d0/
      data tau/1d0,-1d0,1d0,-1d0,1d0,0d0,-1d0,1d0,-1d0,1d0,-1d0/
      real * 8 aemmz
      logical gmuscheme

cccccccccccccccccccccccccccccc
c     Only 4 flavours
      data nflav/4/
cccccccccccccccccccccccccccccc


      write(*,*) '*****************************'
      write(*,*) 'Initializing MCFM couplings'
      write(*,*) '*****************************'

      print*, '************************************'
      print*, '************************************'
      print*, 'CHECK SETTINGS MCFM '
      print*, 'CHECK SETTINGS MCFM '
c     !ER: should select t and not tbar
      if(ttype.eq.1) then
         nwz=1      
         print*, 'nwz is set to 1 !!!!!!!'
      elseif(ttype.eq.-1) then
         nwz=-1      
         print*, 'nwz is set to -1 !!!!!!!'
      else
         write(*,*) 'wrong ttype in virtual.f'
         call exit(-1)
      endif
c     !ER: should correspond to b0 (used to renormalize)
c     Check also nflav=4
      b0= (xn*11d0-2d0*nflav)/6d0       
      print*, 'nflav is set to ',nflav
      print*, 'b0 is set to ',b0
c     !ER: also sck should be checked
      print*, 'sck to be checked'
      print*, '************************************'
      print*, '************************************'

      epinv = 0d0
      epinv2 = 0d0
      zmass = ph_Zmass          ! Z mass
      zwidth = ph_Zwidth        ! Z width

      wmass = ph_Wmass          ! W mass
      wwidth = ph_Wwidth        ! W width
      
      gmuscheme = .false.

      if (gmuscheme) then
         Gf = 1.16639d-05
         wmass  = 80.419d0
c--   derived
         xw  = One-(wmass/zmass)**2     ! sin^2(theta_W)
         aemmz  = Rt2*Gf*wmass**2*xw/pi ! alpha_em
      else
         xw = ph_sthw2          ! sin^2(theta_W)
         aemmz = ph_alphaem     ! alpha_em
      endif      

c--- Now set up the other derived parameters
      gwsq=fourpi*aemmz/xw
      esq=gwsq*xw


      end


      subroutine mom_to_MCFM(cmpborn,pmcfm)
      implicit none
      include 'nlegborn.h'
c      include 'pwhg_flst.h'
c      include 'pwhg_kn.h'
      include 'MCFM_include/constants.f'
      double precision cmpborn(0:3,5)
      double precision pmcfm(mxpart,1:4)
      integer ileg,mu
      do ileg=1,2
         pmcfm(ileg,4)=-cmpborn(0,ileg)
         do mu=1,3
            pmcfm(ileg,mu)=-cmpborn(mu,ileg)
         enddo
      enddo
      do ileg=3,5
         pmcfm(ileg,4)=cmpborn(0,ileg)
         do mu=1,3
            pmcfm(ileg,mu)=cmpborn(mu,ileg)
         enddo    
      enddo
c$$$      do ileg=1,5
c$$$         write(*,*) (pmcfm(ileg,mu), mu=1,4)
c$$$      enddo
      end


      subroutine qg_tbq_v(p,msq,LOqg,LOgq,LOqbarg,LOgqbar)
************************************************************************
*     Virtual t-channel single top, with explicit b-quark              *
*                                                                      *
*     q(p1) + g(p2) -> t(p3) + b(p4) + q'(p5)                          *      
*                                                                      *
*     Originally: R. Frederix and F. Tramontano, February 2008         *
*        Adapted: J. Campbell, February 27, 2008                       *
*                                                                      *
************************************************************************
c     u + g  ->  c + s + d  (t-channel single-charm)
      implicit none
      include 'MCFM_include/constants.f'
      include 'MCFM_include/ewcouple.f'
      include 'MCFM_include/masses.f'
      include 'MCFM_include/scheme.f'
      include 'MCFM_include/ckm.f'
      include 'MCFM_include/nwz.f'
      include 'MCFM_include/stopscales.f'
      include 'MCFM_include/sck.f'
      double precision p(mxpart,4),fac
c--- needed for pole check
c     . ,xs,xsn,xsd
      double precision msq(-nf:nf,-nf:nf),msq_qg,msq_gq,msq_qbarg,
     . msq_gqbar,dot,Wprop15,Wprop25,xsqV,xsqR,mq,ma,gsq_H
      double complex 
     . LOamps_qg(2,2,2),Virtamps_qg(2,2,2),
     . LOamps_qbarg(2,2,2),Virtamps_qbarg(2,2,2),
     . LOamps_gq(2,2,2),Virtamps_gq(2,2,2),
     . LOamps_gqbar(2,2,2),Virtamps_gqbar(2,2,2)
c--- needed for pole check
c     . ,lp,lRc1,lRs1
      integer hg,hc,hs,j,k,i3,i4


cccccccccccccccccccc
c     !ER: check born amplitudes
      double precision LOqg,LOgq,LOqbarg,LOgqbar,sqWpropWidth15
      LOqg=0.
      LOgq=0.
      LOqbarg=0.
      LOgqbar=0.
cccccccccccccccccccc

      sck=1d0

      scheme='dred'

c---initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

c--- DEBUG: to check alpha-dependence      
c      return

c--- set mass of quark and antiquark according to nwz
      if (nwz .eq. +1) then
c     !ER:         print*, 'nwz=1'
        mq=mt
	ma=mb
	i3=3
	i4=4
      else
c     !ER:         print*, 'nwz=-1'
        mq=mb
	ma=mt
	i3=4
	i4=3
      endif

c--- variables to pass renormalization scale in to virtual routines    
      xsqV=renscale_H**2
      xsqR=renscale_H**2

c--- note: factor of ason2pi moved inside virtwrap routine compared
c---       to earlier versions, to allow for different scales
      gsq_H=fourpi*as_H
      fac=aveqg*2d0*xn**2*Cf*gsq_H*gwsq**2

c--- propagator for qg and qbarg amplitudes
      Wprop15=1d0/(2d0*dot(p,1,5)-wmass**2)
      sqWpropWidth15=1d0/((2d0*dot(p,1,5)-wmass**2)**2+wmass**2*wwidth**2) !A28

      call virtwrap(p,1,2,i3,i4,5,
     .              mq,ma,xsqV,xsqR,LOamps_qg,Virtamps_qg)
      call virtwrap(p,5,2,i3,i4,1,
     .              mq,ma,xsqV,xsqR,LOamps_qbarg,Virtamps_qbarg)

c--- propagator for gq and gqbar amplitudes
      Wprop25=1d0/(2d0*dot(p,2,5)-wmass**2)

      call virtwrap(p,2,1,i3,i4,5,
     .              mq,ma,xsqV,xsqR,LOamps_gq,Virtamps_gq)
      call virtwrap(p,5,1,i3,i4,2,
     .              mq,ma,xsqV,xsqR,LOamps_gqbar,Virtamps_gqbar)

      msq_qg=0d0
      msq_qbarg=0d0
      msq_gq=0d0
      msq_gqbar=0d0
      do hg=1,2
      do hc=1,2
      do hs=1,2
c      write(6,*) hg,hc,hs,LOamps_gq(hg,hc,hs),Virtamps_gqbar(hg,hc,hs)
      msq_qg=msq_qg+Wprop15**2*dble(
     .       Virtamps_qg(hg,hc,hs)*dconjg(LOamps_qg(hg,hc,hs)))
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     !A28: need propagator with wwidth to check with madloop point
c     comment 2 lines above and uncomment the following 2
c$$$      msq_qg=msq_qg+sqWpropWidth15*dble(
c$$$     .       Virtamps_qg(hg,hc,hs)*dconjg(LOamps_qg(hg,hc,hs)))
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      msq_qbarg=msq_qbarg+Wprop15**2*dble(
     .       Virtamps_qbarg(hg,hc,hs)*dconjg(LOamps_qbarg(hg,hc,hs)))
      msq_gq=msq_gq+Wprop25**2*dble(
     .       Virtamps_gq(hg,hc,hs)*dconjg(LOamps_gq(hg,hc,hs)))
      msq_gqbar=msq_gqbar+Wprop25**2*dble(
     .       Virtamps_gqbar(hg,hc,hs)*dconjg(LOamps_gqbar(hg,hc,hs)))

cccccccccccccccccccccc
c     !ER: check LO squared. This is the sum over helicities
      LOqg    = LOqg     + Wprop15**2
     $     *abs(LOamps_qg(hg,hc,hs))**2
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     !A28: need propagator with wwidth to check with madloop point.
c     comment 2 lines above and uncomment the following 2
c$$$      LOqg    = LOqg     + sqWpropWidth15
c$$$     $     *abs(LOamps_qg(hg,hc,hs))**2
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      LOqbarg = LOqbarg  + Wprop15**2 
     $     *abs(LOamps_qbarg(hg,hc,hs))**2

      LOgq    = LOgq     + Wprop25**2
     $     *abs(LOamps_gq(hg,hc,hs))**2

      LOgqbar = LOgqbar  + Wprop25**2 
     $     *abs(LOamps_gqbar(hg,hc,hs))**2
cccccccccccccccccccccc

      enddo
      enddo
      enddo

c$$$c--- fill matrix elements
c$$$      do j=1,4
c$$$	msq(+j,0)=fac*Vsum(+j)*msq_qg
c$$$	msq(-j,0)=fac*Vsum(-j)*msq_qbarg
c$$$	msq(0,+j)=fac*Vsum(+j)*msq_gq
c$$$	msq(0,-j)=fac*Vsum(-j)*msq_gqbar
c$$$      enddo

c     !ER: taken away Vsum(), since we include the CKM
c     matrix element at a later point, not here.
c--- fill matrix elements
      do j=1,4
	msq(+j,0)=fac*msq_qg
	msq(-j,0)=fac*msq_qbarg
	msq(0,+j)=fac*msq_gq
	msq(0,-j)=fac*msq_gqbar
      enddo

      LOqg    = LOqg    *fac
      LOgq    = LOgq    *fac
      LOqbarg = LOqbarg *fac
      LOgqbar = LOgqbar *fac
c     Notice that LO squared amplitudes at this point
c     match exactly madgraph ones. Only CKM has to be
c     supplied, all the rest is OK.

      return
      end


c--- wrapper to virtual and LO amplitude routines that allows the
c---  momenta to be permuted according to i1,i2,i5
      subroutine virtwrap(p,i1,i2,i3,i4,i5,
     .                    mh,ml,xsqV,xsqR,LOamps,Virtamps)
      implicit none
      include 'MCFM_include/constants.f'
      include 'MCFM_include/b0.f'
      include 'MCFM_include/epinv.f'
      include 'MCFM_include/scale.f'
      include 'MCFM_include/zprods_com.f'
      include 'MCFM_include/stopscales.f'
      include 'MCFM_include/sck.f'
      integer i1,i2,i3,i4,i5,j
      double precision p(mxpart,4),q(mxpart,4),dot,mh,ml,xsqV,xsqR,
     . colA,colB,virt_massless,eta,ren,ason2pi_H,ason2pi_L
      double complex LOamps(2,2,2),Virtamps(2,2,2),
     . Appp,Appm,Apmp,Apmm,Ampp,Ampm,Ammp,Ammm,
     . Bppp,Bppm,Bpmp,Bpmm,Bmpp,Bmpm,Bmmp,Bmmm,
     . xl15,lnrat
      sck=1d0

c--- factors of ason2pi now included in this routine
      ason2pi_H=as_H/twopi
      ason2pi_L=as_L/twopi
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     !ER: to check values here
c$$$      print*, as_H,as_L,scale,mh,ml
c$$$      print*, renscale_H,facscale_H,renscale_L
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do j=1,4
	q(1,j)=p(i1,j)
	q(2,j)=p(i2,j)
        q(3,j)=p(i3,j)-mh**2/2d0/dot(p,i2,i3)*p(i2,j)
        q(4,j)=p(i4,j)-ml**2/2d0/dot(p,i2,i4)*p(i2,j)
        q(5,j)=p(i5,j)
      enddo

c--- set up spinor products
      call spinoru(5,q,za,zb)

      colA=ca/2d0
      colB=(2d0*cf-ca)/2d0

      eta=0d0  !dred scheme
      ren=(-b0*epinv
     .     -cf*(3d0/2d0*epinv + (4d0+1d0-eta)/2d0 + 3d0*log(scale/mh))
     .     -cf*(1d0/2d0*epinv +
     .             sck*(epinv + (4d0+1d0-eta)/2d0 + 3d0*log(scale/ml)))) 


c--- include finite counterterm to go from DR to MSbar scheme
c--- alphas(DR) = alphas(MSbar) * (1+ (Nc / 6) * alphas(MSbar) / (2*pi))
      ren=ren+xn/6d0

c--- include finite renormalization per FT's message of 20/12/2008
c---  -2*cf*LeadingOrder*alphas/2*pi
      ren=ren-2d0*cf

c--- go from OS scheme for a massive quark to MSbar for massless one
c--- (strong-coupling correction)
      ren=ren-2d0/3d0*log(renscale_H/ml) 
     $     *0d0 !ER: non need for this if we are in the 4f

c--- go from OS scheme for a massive quark to MSbar for massless one
c--- (gluon PDF correction)
      ren=ren-2d0/3d0*log(ml/facscale_H) 
     $     *0d0 !ER: non need for this if we are in the 4f
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     !ER: gluon/quarks gammatilde from DR to CDR
c     for the heavy current
      ren=ren -ca/6. -2.*cf/2.
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      

c--- filling common block
      call stop_def(xsqV,xsqR,q,mh,ml)
c---- calling amps(hg,hc,hs)
      call bornampsN(q,mh,ml,LOamps)

      call Aamp_ppp(q,mh,ml,Appp)
      call Aamp_ppm(q,mh,ml,Appm)
      call Aamp_pmp(q,mh,ml,Apmp)
      call Aamp_pmm(q,mh,ml,Apmm)
      call Aamp_mpp(q,mh,ml,Ampp)
      call Aamp_mpm(q,mh,ml,Ampm)
      call Aamp_mmp(q,mh,ml,Ammp)
      call Aamp_mmm(q,mh,ml,Ammm)
      call Bamp_ppp(q,mh,ml,Bppp)
      call Bamp_ppm(q,mh,ml,Bppm)
      call Bamp_pmp(q,mh,ml,Bpmp)
      call Bamp_pmm(q,mh,ml,Bpmm)
      call Bamp_mpp(q,mh,ml,Bmpp)
      call Bamp_mpm(q,mh,ml,Bmpm)
      call Bamp_mmp(q,mh,ml,Bmmp)
      call Bamp_mmm(q,mh,ml,Bmmm)
 
      xl15=lnrat(-2d0*dot(p,i1,i5),renscale_L**2)
c--- correction to the massless line (cf. cv0 in qqb_tbb_v.f)
      virt_massless=-2d0*epinv*(epinv-dble(xl15))-dble(xl15**2)
     .              -3d0*(epinv-dble(xl15))-7d0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     !ER: I need this change (-7 -> -8) to go to tHV also for
c     the massless current.
c     Look also virtual.f in the ST_tch folder (5-flavour scheme)
c     and the checks I did with MCFM.
c     This is clearly the same as doing 
c     ren=ren -2.*cf/2
c     i.e. adding quarks gammatilde to go from DR to CDR for
c     the massless current.
c     I obtain the value in the madloop paper only if I add
c     the gammatilde only for the heavy current and not for the
c     light one, but this looks to me inconsistent, unless the MCFM routine
c     was built such that at the end dred was used only for the heavy current.
c      virt_massless=-2d0*epinv*(epinv-dble(xl15))-dble(xl15**2)
c     .              -3d0*(epinv-dble(xl15))-7d0
c     $     -1d0 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      virt_massless=virt_massless*cf

c--- apply factors of ason2pi now
      colA=colA*ason2pi_H
      colB=colB*ason2pi_H    ! corrections on heavy line
      virt_massless=virt_massless*ason2pi_L ! on light line
      ren=ren*ason2pi_H ! renormalization is like LO
 
      Virtamps(2,2,2)=Appp*colA+Bppp*colB
     .               +(ren+virt_massless)*LOamps(2,2,2)
      Virtamps(2,2,1)=Appm*colA+Bppm*colB
     .               +(ren+virt_massless)*LOamps(2,2,1)
      Virtamps(2,1,2)=Apmp*colA+Bpmp*colB
     .               +(ren+virt_massless)*LOamps(2,1,2)
      Virtamps(2,1,1)=Apmm*colA+Bpmm*colB
     .               +(ren+virt_massless)*LOamps(2,1,1)
      Virtamps(1,2,2)=Ampp*colA+Bmpp*colB
     .               +(ren+virt_massless)*LOamps(1,2,2)
      Virtamps(1,2,1)=Ampm*colA+Bmpm*colB
     .               +(ren+virt_massless)*LOamps(1,2,1)
      Virtamps(1,1,2)=Ammp*colA+Bmmp*colB
     .               +(ren+virt_massless)*LOamps(1,1,2)
      Virtamps(1,1,1)=Ammm*colA+Bmmm*colB
     .               +(ren+virt_massless)*LOamps(1,1,1)

      return



************************************************************************
*   CODE BELOW HERE IS FOR CHECKING POLES ONLY                         *
************************************************************************


c      call checkint(p,xsqV,xsqR,ml**2,mh**2)

c--- needed for pole check

c      xsn=(1d0-dsqrt(1d0-4d0*ml*mh/(2*dot(p,3,4)+2*ml*mh)))
c      xsd=(1d0+dsqrt(1d0-4d0*ml*mh/(2*dot(p,3,4)+2*ml*mh)))
c      xs=-xsn/xsd
c--- these are to check the qg terms
c      lRc1=lnrat(mh*sqrt(xsqR),-2*dot(p,2,3))
c      lRs1=lnrat(ml*sqrt(xsqR),-2*dot(p,2,4))
c      lp=lnrat(xsn,-xsd)

c--- these are to check the gq terms
c      lRc1=lnrat(mh*sqrt(xsqR),-2*dot(p,1,3))
c      lRs1=lnrat(ml*sqrt(xsqR),-2*dot(p,1,4))
c--- pole check
c
c      write(6,*) 'Appp/LOamp(2,2,2)=',Appp/(+LOamps(2,2,2)
c     . *(-2d0*epinv**2 - epinv*(-1d0 + 2d0*lRc1 + 2d0*lRs1)))
c      write(6,*) 'Appm/LOamp(2,2,1)=',Appm/(+LOamps(2,2,1)
c     . *(-2d0*epinv**2 - epinv*(-1d0 + 2d0*lRc1 + 2d0*lRs1)))
c      write(6,*) 'Apmp/LOamp(2,1,2)=',Apmp/(+LOamps(2,1,2)
c     . *(-2d0*epinv**2 - epinv*(-1d0 + 2d0*lRc1 + 2d0*lRs1)))
c      write(6,*) 'Apmm/LOamp(2,1,1)=',Apmm/(+LOamps(2,1,1)
c     . *(-2d0*epinv**2 - epinv*(-1d0 + 2d0*lRc1 + 2d0*lRs1)))
c      write(6,*) 'Ampp/LOamp(1,2,2)=',Ampp/(+LOamps(1,2,2)
c     . *(-2d0*epinv**2 - epinv*(-1d0 + 2d0*lRc1 + 2d0*lRs1)))
c      write(6,*) 'Ampm/LOamp(1,2,1)=',Ampm/(+LOamps(1,2,1)
c     . *(-2d0*epinv**2 - epinv*(-1d0 + 2d0*lRc1 + 2d0*lRs1)))
c      write(6,*) 'Ammp/LOamp(1,1,2)=',Ammp/(+LOamps(1,1,2)
c     . *(-2d0*epinv**2 - epinv*(-1d0 + 2d0*lRc1 + 2d0*lRs1)))
c      write(6,*) 'Ammm/LOamp(1,1,1)=',Ammm/(+LOamps(1,1,1)
c     . *(-2d0*epinv**2 - epinv*(-1d0 + 2d0*lRc1 + 2d0*lRs1)))

c      write(6,*) 'Bppp/LOamp(2,2,2)=',Bppp/(+LOamps(2,2,2)
c     .    *epinv*(-4d0*dot(p,3,4)*lp*xs + mh*ml*(-1d0 + xs**2))/
c     .  (ml*(-1d0 + xs**2))/mh)
c      write(6,*) 'Bppm/LOamp(2,2,1)=',Bppm/(+LOamps(2,2,1)
c     .    *epinv*(-4d0*dot(p,3,4)*lp*xs + mh*ml*(-1d0 + xs**2))/
c     .  (ml*(-1d0 + xs**2))/mh)
c      write(6,*) 'Bpmp/LOamp(2,1,2)=',Bpmp/(+LOamps(2,1,2)
c     .    *epinv*(-4d0*dot(p,3,4)*lp*xs + mh*ml*(-1d0 + xs**2))/
c     .  (ml*(-1d0 + xs**2))/mh)
c      write(6,*) 'Bpmm/LOamp(2,1,1)=',Bpmm/(+LOamps(2,1,1)
c     .    *epinv*(-4d0*dot(p,3,4)*lp*xs + mh*ml*(-1d0 + xs**2))/
c     .  (ml*(-1d0 + xs**2))/mh)
c      write(6,*) 'Bmpp/LOamp(1,2,2)=',Bmpp/(+LOamps(1,2,2)
c     .    *epinv*(-4d0*dot(p,3,4)*lp*xs + mh*ml*(-1d0 + xs**2))/
c     .  (ml*(-1d0 + xs**2))/mh)
c      write(6,*) 'Bmpm/LOamp(1,2,1)=',Bmpm/(+LOamps(1,2,1)
c     .    *epinv*(-4d0*dot(p,3,4)*lp*xs + mh*ml*(-1d0 + xs**2))/
c     .  (ml*(-1d0 + xs**2))/mh)
c      write(6,*) 'Bmmp/LOamp(1,1,2)=',Bmmp/(+LOamps(1,1,2)
c     .    *epinv*(-4d0*dot(p,3,4)*lp*xs + mh*ml*(-1d0 + xs**2))/
c     .  (ml*(-1d0 + xs**2))/mh)
c      write(6,*) 'Bmmm/LOamp(1,1,1)=',Bmmm/(+LOamps(1,1,1)
c     .    *epinv*(-4d0*dot(p,3,4)*lp*xs + mh*ml*(-1d0 + xs**2))/
c     .  (ml*(-1d0 + xs**2))/mh)

c      pause
c      return
      
      end


      subroutine bornampsN(q,mc,ms,amps)
c     u + g  ->  c + s + d  (t-channel single-charm)
      implicit none
      include 'MCFM_include/constants.f'
      include 'MCFM_include/zprods_com.f'
      double precision q(mxpart,4),dot,cDs,gDs,cDg,mc,ms
      double complex trg,trs,trc,trsgc,amps(2,2,2)
      
      cDg=dot(q,3,2)
      gDs=dot(q,4,2)
      cDs=dot(q,3,4)+mc**2*gDs/2d0/cDg
     &              +ms**2*cDg/2d0/gDs

      trg=2d0*za(5,2)*zb(2,1)
      trs=2d0*za(5,4)*zb(4,1)+ms**2*za(5,2)*zb(2,1)/gDs
      trc=2d0*za(5,3)*zb(3,1)+mc**2*za(5,2)*zb(2,1)/cDg
      trsgc=2d0*zb(1,4)*za(4,2)*zb(2,3)*za(3,5)


      amps(2,2,2)=(mc*(cDg*(-ms**2*trg+2d0*gDs*(trg + trs))- 
     -      gDs*trsgc))/(4d0*cDg*gDs)

      amps(2,2,1)=-(mc*ms*(-2d0*cDg*cDs*gDs+gDs**2*mc**2+ 
     -       cDg**2*ms**2)*trg)/(4d0*cDg*gDs**2)

      amps(2,1,2)=(gDs**2*mc**2*trsgc-2d0*cDg*gDs*(gDs*mc**2*trg
     -  +cDs*trsgc)+cDg**2*(4d0*gDs**2*trc-2d0*gDs*trsgc + 
     -       ms**2*trsgc))/(4d0*cDg**2*gDs)

      amps(2,1,1)=-(ms*(-2d0*cDg*cDs*gDs+gDs**2*mc**2+ 
     -   cDg**2*ms**2)*(2d0*cDg*gDs*trc-gDs*mc**2*trg-cDg*trsgc))/
     -  (4d0*cDg**2*gDs**2)

      amps(1,2,2)=-(mc*(-2d0*cDg*cDs*gDs+gDs**2*mc**2+cDg**2*ms**2)*
     -     (-(cDg*ms**2*trg) + 2d0*cDg*gDs*trs - gDs*trsgc)
     -     )/(4d0*cDg**2*gDs**2)

      amps(1,2,1)=(mc*ms*(-2d0*cDg*cDs*gDs + gDs**2*mc**2 + 
     -      cDg**2*ms**2)*trg)/(4d0*cDg**2*gDs)

      amps(1,1,2)=-(-2d0*cDg*gDs*(cDs+gDs)*trsgc+gDs**2*mc**2*trsgc+
     -     cDg**2*(-2d0*gDs*ms**2*trg + 4d0*gDs**2*trs + 
     -        ms**2*trsgc))/(4d0*cDg*gDs**2)

      amps(1,1,1)=(ms*(-gDs*mc**2*trg+cDg*(2d0*gDs*(trc+trg)-trsgc)))/
     -  (4d0*cDg*gDs)


      amps(2,2,2)=amps(2,2,2)/za(2,4)/za(2,3)
      amps(2,2,1)=amps(2,2,1)/za(2,3)**2/zb(4,3)
      amps(2,1,2)=amps(2,1,2)/za(2,4)**2/zb(4,3)
      amps(2,1,1)=amps(2,1,1)/za(2,4)/za(2,3)/zb(4,3)**2
      amps(1,2,2)=amps(1,2,2)/za(4,3)**2/zb(2,4)/zb(2,3)
      amps(1,2,1)=amps(1,2,1)/za(4,3)/zb(2,4)**2
      amps(1,1,2)=amps(1,1,2)/za(4,3)/zb(2,3)**2
      amps(1,1,1)=amps(1,1,1)/zb(2,4)/zb(2,3)

      return
      end
      
      double complex function Lnrat(x,y)
************************************************************************
*     Author: R.K. Ellis                                               *
*     August, 1998.                                                    *
c     Lnrat(x,y)=log(x-i*ep)-log(y-i*ep)                               *
c     this function is hard-wired for sign of epsilon we must adjust   *
c     sign of x and y to get the right sign for epsilon                *
************************************************************************
      implicit none
      include 'MCFM_include/constants.f'
      double precision x,y,htheta
C--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
      htheta(x)=half+half*sign(one,x)
      Lnrat=dcmplx(dlog(abs(x/y)))-impi*dcmplx((htheta(-x)-htheta(-y)))
      return
      end



*
* $Id: dclaus64.F,v 1.2 1996/04/02 16:23:45 mclareni Exp $
*
* $Log: dclaus64.F,v $
* Revision 1.2  1996/04/02 16:23:45  mclareni
* More precise dclaus64 (C326), test added and C344 removed from TESTALL
*
* Revision 1.1.1.1  1996/04/01 15:02:03  mclareni
* Mathlib gen
*
*
      DOUBLE PRECISION FUNCTION DCLAUS(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION A(0:8),B(0:13)
 
      PARAMETER (R1 = 1d0, HF =R1/2d0)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (PI2 = 2d0*PI, PIH = PI/2d0, RPIH = 2d0/PI)
 
      DATA A( 0) / 0.02795 28319 73575 6613D0/
      DATA A( 1) / 0.00017 63088 74389 8116D0/
      DATA A( 2) / 0.00000 12662 74146 1157D0/
      DATA A( 3) / 0.00000 00117 17181 8134D0/
      DATA A( 4) / 0.00000 00001 23006 4129D0/
      DATA A( 5) / 0.00000 00000 01395 2729D0/
      DATA A( 6) / 0.00000 00000 00016 6908D0/
      DATA A( 7) / 0.00000 00000 00000 2076D0/
      DATA A( 8) / 0.00000 00000 00000 0027D0/
 
      DATA B( 0) / 0.63909 70888 57265 341D0/
      DATA B( 1) /-0.05498 05693 01851 716D0/
      DATA B( 2) /-0.00096 12619 45950 606D0/
      DATA B( 3) /-0.00003 20546 86822 550D0/
      DATA B( 4) /-0.00000 13294 61695 426D0/
      DATA B( 5) /-0.00000 00620 93601 824D0/
      DATA B( 6) /-0.00000 00031 29600 656D0/
      DATA B( 7) /-0.00000 00001 66351 954D0/
      DATA B( 8) /-0.00000 00000 09196 527D0/
      DATA B( 9) /-0.00000 00000 00524 004D0/
      DATA B(10) /-0.00000 00000 00030 580D0/
      DATA B(11) /-0.00000 00000 00001 820D0/
      DATA B(12) /-0.00000 00000 00000 110D0/
      DATA B(13) /-0.00000 00000 00000 007D0/
 
      V=MOD(ABS(X),PI2)
      S=SIGN(R1,X)
      IF(V .GT. PI) THEN
       V=PI2-V
       S=-S
      ENDIF
      IF(V .EQ. 0d0 .OR. V .EQ. PI) THEN
       H=0d0
      ELSEIF(V .LT. PIH) THEN
       U=RPIH*V
       H=2d0*U**2-1d0
       ALFA=H+H
       B1=0d0
       B2=0d0
       DO 1 I = 8,0,-1
       B0=A(I)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=V*(1d0-LOG(V)+HF*V**2*(B0-H*B2))
      ELSE
       U=RPIH*V-2d0
       H=2d0*U**2-1d0
       ALFA=H+H
       B1=0d0
       B2=0d0
       DO 2 I = 13,0,-1
       B0=B(I)+ALFA*B1-B2
       B2=B1
    2  B1=B0
       H=(PI-V)*(B0-H*B2)
      ENDIF
      DCLAUS=S*H
      RETURN
      END



      DOUBLE PRECISION FUNCTION DDILOGMCFM(X)

      DOUBLE PRECISION X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO
      DOUBLE PRECISION C(0:18),H,ALFA,B0,B1,B2

      DATA ZERO /0.0D0/, ONE /1.0D0/
      DATA HALF /0.5D0/, MALF /-0.5D0/, MONE /-1.0D0/, MTWO /-2.0D0/
      DATA PI3 /3.28986 81336 96453D0/, PI6 /1.64493 40668 48226D0/

      DATA C( 0) / 0.42996 69356 08137 0D0/
      DATA C( 1) / 0.40975 98753 30771 1D0/
      DATA C( 2) /-0.01858 84366 50146 0D0/
      DATA C( 3) / 0.00145 75108 40622 7D0/
      DATA C( 4) /-0.00014 30418 44423 4D0/
      DATA C( 5) / 0.00001 58841 55418 8D0/
      DATA C( 6) /-0.00000 19078 49593 9D0/
      DATA C( 7) / 0.00000 02419 51808 5D0/
      DATA C( 8) /-0.00000 00319 33412 7D0/
      DATA C( 9) / 0.00000 00043 45450 6D0/
      DATA C(10) /-0.00000 00006 05784 8D0/
      DATA C(11) / 0.00000 00000 86121 0D0/
      DATA C(12) /-0.00000 00000 12443 3D0/
      DATA C(13) / 0.00000 00000 01822 6D0/
      DATA C(14) /-0.00000 00000 00270 1D0/
      DATA C(15) / 0.00000 00000 00040 4D0/
      DATA C(16) /-0.00000 00000 00006 1D0/
      DATA C(17) / 0.00000 00000 00000 9D0/
      DATA C(18) /-0.00000 00000 00000 1D0/

      IF(X .EQ. ONE) THEN
       DDILOGMCFM=PI6
       RETURN
      ELSE IF(X .EQ. MONE) THEN
       DDILOGMCFM=MALF*PI6
       RETURN
      END IF
      T=-X
      IF(T .LE. MTWO) THEN
       Y=MONE/(ONE+T)
       S=ONE
       A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)
      ELSE IF(T .LT. MONE) THEN
       Y=MONE-T
       S=MONE
       A=LOG(-T)
       A=-PI6+A*(A+LOG(ONE+ONE/T))
      ELSE IF(T .LE. MALF) THEN
       Y=(MONE-T)/T
       S=ONE
       A=LOG(-T)
       A=-PI6+A*(MALF*A+LOG(ONE+T))
      ELSE IF(T .LT. ZERO) THEN
       Y=-T/(ONE+T)
       S=MONE
       A=HALF*LOG(ONE+T)**2
      ELSE IF(T .LE. ONE) THEN
       Y=T
       S=ONE
       A=ZERO
      ELSE
       Y=ONE/T
       S=MONE
       A=PI6+HALF*LOG(T)**2
      END IF

      H=Y+Y-ONE
      ALFA=H+H
      B1=ZERO
      B2=ZERO
      DO 1 I = 18,0,-1
      B0=C(I)+ALFA*B1-B2
      B2=B1
    1 B1=B0
      DDILOGMCFM=-(S*(B0-H*B2)+A)
      RETURN
      END


      double precision function dot(p,i,j)
      implicit none
      include 'MCFM_include/constants.f'
      integer i,j
      double precision p(mxpart,4)
      dot=p(i,4)*p(j,4)-p(i,1)*p(j,1)-p(i,2)*p(j,2)-p(i,3)*p(j,3)
      return
      end


      subroutine spinoru(N,p,za,zb)
c---Calculate spinor products
c---extended to deal with negative energies ie with all momenta outgoing
c---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl, 
c---za(i,j)*zb(j,i)=s(i,j)
      implicit none
      include 'MCFM_include/constants.f'
      include 'MCFM_include/zprods_decl.f'
      include 'MCFM_include/sprods_com.f'
      double precision p(mxpart,4),rt(mxpart)
      double complex c23(mxpart),f(mxpart)
      integer i,j,N
      
c---if one of the vectors happens to be zero this routine fails.
      do j=1,N
         za(j,j)=czip
         zb(j,j)=za(j,j)

C-----positive energy case
         if (p(j,4) .gt. 0d0) then
            rt(j)=dsqrt(p(j,4)+p(j,1))
            c23(j)=dcmplx(p(j,3),-p(j,2))
            f(j)=cone
         else
C-----negative energy case
            rt(j)=dsqrt(-p(j,4)-p(j,1))
            c23(j)=dcmplx(-p(j,3),p(j,2))
            f(j)=im
         endif
      enddo
      do i=2,N
         do j=1,i-1
         s(i,j)=two*(p(i,4)*p(j,4)-p(i,1)*p(j,1)
     &              -p(i,2)*p(j,2)-p(i,3)*p(j,3))
         za(i,j)=f(i)*f(j)
     &   *(c23(i)*dcmplx(rt(j)/rt(i))-c23(j)*dcmplx(rt(i)/rt(j)))

         if (abs(s(i,j)).lt.1d-9) then
         zb(i,j)=-(f(i)*f(j))**2*dconjg(za(i,j))
         else
         zb(i,j)=-dcmplx(s(i,j))/za(i,j)
         endif
         za(j,i)=-za(i,j)
         zb(j,i)=-zb(i,j)
         s(j,i)=s(i,j)
         enddo
      enddo

      return
      end


      subroutine stop_def(xsqV,xsqR,q,mc,ms)
c--- Subroutine to be run for each event to calculate log's, bubbles,
c--- triangles and boxes (A0,B0,C0,C00,C001,C002,D00) and fill the
c--- common block in stopf1inc.f
      implicit none
      include 'MCFM_include/constants.f'
      double precision xsqV,xsqR,q(mxpart,4),ms,mc
      double precision muR,ms2,mc2,xsn,xsd,s,t,u,qsq,dot,cDs
cccccccccccccccccccccccccccccccccccccccccccccccccc
c     !ER: changed needed to comply with our include structure
c$$$      include 'MCFM_include/stopf1inc.f'
      double complex lc,ls,lVc,lVs,lp
     &,LsA,LsB1,LsB2,tr1Xfc,tr1Xfs,tr2fu,tr3Xc
     &,tr3c00fs,tr3c001fs,tr3c002fs,tr3Xs,tr3s00ft,tr3s001ft
     &,tr3s002ft,tr4Xc,tr4Xs,tr5Xc,tr5Xs,B0csf,B0cgsf
     &,lRc1,lRc2,lRs1,lRs2,lRcs,BfunX
      common/stopf1inc/ lc,ls,lVc,lVs,lp
     &,LsA,LsB1,LsB2,tr1Xfc,tr1Xfs,tr2fu,tr3Xc
     &,tr3c00fs,tr3c001fs,tr3c002fs,tr3Xs,tr3s00ft,tr3s001ft
     &,tr3s002ft,tr4Xc,tr4Xs,tr5Xc,tr5Xs,B0csf,B0cgsf
     &,lRc1,lRc2,lRs1,lRs2,lRcs,BfunX
cccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex lnrat,LLs1,LLs2,tr1f,tr2f,tr3,tr3c00f,tr3c001f
     .,tr3c002f,tr3s00f,tr3s001f,tr3s002f,tr4,tr5,B0xf,Bfun
      external lnrat,LLs1,LLs2,tr1f,tr2f,tr3,tr3c00f,tr3c001f
     .,tr3c002f,tr3s00f,tr3s001f,tr3s002f,tr4,tr5,B0xf,Bfun

      ms2 = ms**2
      mc2 = mc**2
      cDs = dot(q,3,4)+mc2*dot(q,4,2)/2d0/dot(q,3,2)
     .  +ms2*dot(q,3,2)/2d0/dot(q,4,2)
      qsq = ms2+mc2+2d0*cDs+2d0*dot(q,3,2)+2d0*dot(q,4,2)
      s   = ms2+2d0*dot(q,4,2)
      t   = mc2+2d0*dot(q,3,2)
      u   = mc2+ms2+2d0*cDs
      xsn = (1d0-dsqrt(1d0-4d0*ms*mc/(u-(ms-mc)**2)))
      xsd = (1d0+dsqrt(1d0-4d0*ms*mc/(u-(ms-mc)**2)))
      muR = dsqrt(xsqR)

      lRc1      = lnrat(mc*muR,mc2-t)
      lRs1      = lnrat(ms*muR,ms2-s)
      lRc2      = lRc1**2
      lRs2      = lRs1**2
      lRcs      = lnrat(xsn,-xsd)*lnrat(xsqR,mc*ms)
      lc        = lnrat(mc2-t,mc2)
      ls        = lnrat(ms2-s,ms2)
      lp        = lnrat(xsn,-xsd)
      lVc       = lnrat(xsqV,mc2)
      lVs       = lnrat(xsqV,ms2)
      B0csf     = B0xf(xsqV,u,ms2,mc2)
      B0cgsf    = B0xf(xsqV,qsq,ms2,mc2)
      BfunX     = Bfun(xsqV,qsq,ms2,mc2)
      tr1Xfc    = tr1f(xsqR,t,mc2)
      tr1Xfs    = tr1f(xsqR,s,ms2)
      tr2fu     = tr2f(xsqR,u,ms2,mc2)
      tr3Xc     = tr3(s,qsq,mc2,ms2)
      tr3Xs     = tr3(t,qsq,ms2,mc2)
      tr3s00ft  = tr3s00f(xsqV,t,qsq,ms2,mc2)
      tr3s001ft = tr3s001f(xsqV,t,qsq,ms2,mc2)
      tr3s002ft = tr3s002f(xsqV,t,qsq,ms2,mc2)
      tr3c00fs  = tr3c00f(xsqV,s,qsq,mc2,ms2)
      tr3c001fs = tr3c001f(xsqV,s,qsq,mc2,ms2)
      tr3c002fs = tr3c002f(xsqV,s,qsq,mc2,ms2)
      tr4Xc     = tr4(t,mc2)
      tr4Xs     = tr4(s,ms2)
      tr5Xc     = tr5(u,qsq,mc2,ms2)
      tr5Xs     = tr5(u,qsq,ms2,mc2)
      LsA       = LLs1(xsqR,s,t,qsq,ms2,mc2)
      LsB1      = LLs2(xsqR,s,u,qsq,mc2,ms2)
      LsB2      = LLs2(xsqR,t,u,qsq,ms2,mc2)

      return
      end



