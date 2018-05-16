      subroutine setscalesbtilde
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_flg.h'
      real * 8 pwhg_alphas
      external pwhg_alphas
      real * 8 muf,mur
c signal we will begin by computing Born type contributions
      flg_btildepart='b'
      call set_fac_ren_scales(muf,mur)
      st_mufact2= muf**2*st_facfact**2
      st_muren2 = mur**2*st_renfact**2
      st_alpha  = pwhg_alphas(st_muren2,st_lambda5MSB,st_nlight)
      end

      subroutine setscalesbtlreal
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_flg.h'
      real * 8 pwhg_alphas
      external pwhg_alphas
      real * 8 muf,mur
      logical ini
      data ini/.true./
      save ini
      real * 8 powheginput
      external powheginput
      if(ini) then
         if(powheginput("#btlscalereal").eq.1d0) then
            flg_btlscalereal=.true.
         else
            flg_btlscalereal=.false.
         endif
         ini=.false.
      endif
      if(flg_btlscalereal) then
c if this is active we may compute scales that depends upon
c the real kinematics; the user routine set_fac_ren_scales
c should test the flag flg_btildepart to see if this is the case
         flg_btildepart='r'
         call set_fac_ren_scales(muf,mur)
         st_mufact2= muf**2*st_facfact**2
         st_muren2 = mur**2*st_renfact**2
         st_alpha  = pwhg_alphas(st_muren2,st_lambda5MSB,st_nlight)
      endif
      end

      subroutine setscalesbtlct
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_flg.h'
      real * 8 pwhg_alphas
      external pwhg_alphas
      real * 8 muf,mur
      logical ini
      data ini/.true./
      save ini
      real * 8 powheginput
      external powheginput
      if(ini) then
         if(powheginput("#btlscalereal").eq.1d0) then
            flg_btlscalereal=.true.
         else
            flg_btlscalereal=.false.
         endif
         if(powheginput("#btlscalect").eq.1d0) then
            flg_btlscalect=.true.
         else
            flg_btlscalect=.false.
         endif
      endif
      if(flg_btlscalereal.and.flg_btlscalect) then
c signal we will begin by computing counterterm contributions, in cases
c when it is desirable to have the scales of the counterterm differ from
c those of the real term (the token btlscalect selects this case)
c The user routine should test the flag flg_btildepart to see if
c we are in a counterterm.
         flg_btildepart='c'
         call set_fac_ren_scales(muf,mur)
         st_mufact2= muf**2*st_facfact**2
         st_muren2 = mur**2*st_renfact**2
         st_alpha  = pwhg_alphas(st_muren2,st_lambda5MSB,st_nlight)
      endif
      end

      subroutine set_rad_scales(ptsq)
      implicit none
      real * 8 ptsq
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      include 'pwhg_rad.h'
      real * 8 pwhg_alphas
      integer nf
      external pwhg_alphas
ccccccccccccccccccccccccccccccccc
c     SAER FIX
      real * 8 q2min
      character * 3 whichpdfpk
      external whichpdfpk
      integer ini,mem
      data ini/0/
      save ini,q2min
      
      if (ini.eq.0) then
         if( whichpdfpk().eq.'lha') then    
           mem = 0
           call GetQ2min(mem,q2min)
c     the previous value of q2min is not the value for which pdf is not
c     zero but the minimum value of Q^2 in pdf grids. In the case of
c     heavy quarks involved one should use their masses as minimum value
c     of factorization scale, as we make later on. This works if ptmin
c     is greater or equal to the mass of heavy quark 
        elseif( whichpdfpk().eq.'mlm') then    
c ad hoc value here (mlmpdf does not provide this)
           q2min=2d0
        else
           write(*,*) ' unimplemented pdf package',whichpdfpk()
           stop
        endif 
        ini=1
      endif
      st_mufact2=max(q2min,ptsq) 
cccccccccccccccccccccccccccccccccc
c     In case of final-state radiation, Born and real PDF's
c     should always cancel out in the ratio R/B. If the radiation scale
c     is too low, this cancellation can be spoilt because PDF's can vanish,
c     typically when a heavy flavour is present as initial state.
c     To prevent this, we use a scale higher than the heavy-flavour
c     threshold, so that PDF's are evaluated with a safe value for
c     mufact (50 is an arbitrary choice).
      if(rad_kinreg.ge.2) st_mufact2=50.**2
      st_muren2=ptsq
      st_alpha = pwhg_alphas(st_muren2,st_lambda5MSB,-1)
      if(st_muren2.lt.rad_charmthr2) then
         nf=3
      elseif(st_muren2.lt.rad_bottomthr2) then
         nf=4
      else
         nf=5
      endif
      st_alpha = st_alpha *
     #   (1+st_alpha/(2*pi)*((67d0/18-pi**2/6)*ca-5d0/9*nf))
      end

      subroutine init_rad_lambda
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      real * 8 b0,mu0sq,as,pwhg_alphas
      external pwhg_alphas
      b0=(33-2*5)/(12*pi)
      mu0sq=(2*st_lambda5MSB)**2
c running value of alpha at initial scale (see notes: running_coupling)
      as=pwhg_alphas(mu0sq,st_lambda5MSB,-1)
c for better NLL accuracy (FNO2006, (4.32) and corresponding references)
      as=as*(1+as/(2*pi)*((67d0/18-pi**2/6)*ca-5d0/9*3))
      rad_lamll=sqrt(exp(-1/(b0*as))*mu0sq)
      end

      function pwhg_alphas0(q2,xlam,inf)
      implicit none
      real * 8 pwhg_alphas0,q2,xlam
      integer inf
      real * 8 pi
      parameter (pi=3.141592653589793d0)
      real * 8 b0
      b0=(33-2*inf)/(12*pi)
      pwhg_alphas0=1/(b0*log(q2/xlam**2))
      end
      
C----------------------------------------------------------------------------
C-------------------------------------------------------------------
C------- ALPHA QCD -------------------------------------
c Program to calculate alfa strong with nf flavours,
c as a function of lambda with 5 flavors.
c The value of alfa is matched at the thresholds q = mq.
c When invoked with nf < 0 it chooses nf as the number of
c flavors with mass less then q.
c
      function pwhg_alphas(q2,xlam,inf)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      real * 8 pwhg_alphas,q2,xlam
      integer inf
      real * 8 pi
      parameter (pi=3.141592653589793d0)
      real * 8 olam,b5,bp5,b4,bp4,b3,bp3,xlc,xlb,xllc,xllb,c45,c35,
     # xmc,xmb
      real * 8 q,xlq,xllq
      integer nf
      data olam/0.d0/
      save olam,b5,bp5,b4,bp4,b3,bp3,xlc,xlb,xllc,xllb,c45,c35,xmc,xmb

c      logical ini
c      data ini/.true./
c      save ini
c      pwhg_alphas = 30d0
c      if (ini) then
c         write(*,*) '****************************************'
c         write(*,*) '****************************************'
c         write(*,*) '      RETURN alpha_s = ',pwhg_alphas
c         write(*,*) '****************************************'
c         write(*,*) '****************************************'
c         ini = .false.         
c      endif           
c      return

      if(xlam.ne.olam) then
        olam = xlam
        xmc=sqrt(rad_charmthr2)
        xmb=sqrt(rad_bottomthr2)
        b5  = (33-2*5)/pi/12
        bp5 = (153 - 19*5) / pi / 2 / (33 - 2*5)
        b4  = (33-2*4)/pi/12
        bp4 = (153 - 19*4) / pi / 2 / (33 - 2*4)
        b3  = (33-2*3)/pi/12
        bp3 = (153 - 19*3) / pi / 2 / (33 - 2*3)
        xlc = 2 * log(xmc/xlam)
        xlb = 2 * log(xmb/xlam)
        xllc = log(xlc)
        xllb = log(xlb)
        c45  =  1/( 1/(b5 * xlb) - xllb*bp5/(b5 * xlb)**2 )
     #        - 1/( 1/(b4 * xlb) - xllb*bp4/(b4 * xlb)**2 )
        c35  =  1/( 1/(b4 * xlc) - xllc*bp4/(b4 * xlc)**2 )
     #        - 1/( 1/(b3 * xlc) - xllc*bp3/(b3 * xlc)**2 ) + c45
      endif
      q   = sqrt(q2)
      xlq = 2 * log( q/xlam )
      xllq = log( xlq )
      nf = inf
      if( nf .lt. 0) then
        if( q .gt. xmb ) then
          nf = 5
        elseif( q .gt. xmc ) then
          nf = 4
        else
          nf = 3
        endif
      endif
      if    ( nf .eq. 5 ) then
        pwhg_alphas = 1/(b5 * xlq) -  bp5/(b5 * xlq)**2 * xllq
      elseif( nf .eq. 4 ) then
        pwhg_alphas =
     #    1/( 1/(1/(b4 * xlq) - bp4/(b4 * xlq)**2 * xllq) + c45 )
      elseif( nf .eq. 3 ) then
        pwhg_alphas =
     #    1/( 1/(1/(b3 * xlq) - bp3/(b3 * xlq)**2 * xllq) + c35 )
      else
        print *,'error in alfa: unimplemented # of light flavours',nf
        call exit(1)
      endif
      return
      end
