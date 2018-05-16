      subroutine gg_hg_eval_v(p,res)
!---- Modified for i2mcfm CW August 11
C----Author: R.K. Ellis May 2007
c----Virtual corrections matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
C----f(p1)+f(p2) --> H(b(p3)+b(p4)+g(p5)
      implicit none
      include 'MCFM_Include/constants.f'
!      include 'masses.f'
!     Will take masses from Madgraph
      include 'MCFM_Include/ewcouple.f'
      include 'MCFM_Include/qcdcouple.f'
!     Couplings are filled by pwhg_st.h in driver routine
      include 'coupl.inc'
      include 'MCFM_Include/sprods_com.f'
      include 'MCFM_Include/scheme.f'
      include 'MCFM_Include/interface_settings.f'
C     (Taken from Ravindran, Smith, van Neerven hep-ph/0201114)
C     Modified by overall factors
      integer iglue,j,k
      double precision p(mxpart,4),res,s34,mb,mbsq
      double precision ss,tt,uu,
     . virtgg,virtqa,virtaq,virtqg,virtgq,hdecay,Asq,fac
      logical rem_dec
      parameter(iglue=5)
      double precision oldss,oldtt,olduu
      save oldss,oldtt,olduu,virtgg,virtqa,virtaq,virtqg,virtgq
      data oldss/0d0/
      data oldtt/0d0/
      data olduu/0d0/

      scheme='dred'
!      call check_scheme
      rem_dec=.true. 
      call dotem(iglue,p,s)
      ss=s(1,2)
      tt=s(1,iglue)
      uu=s(2,iglue)

      mb=bmass
      mbsq=mb**2
      Asq=(as/(3d0*pi))**2/vevsq

   !   write(6,*) vevsq,as,gsq 
  !    pause
C   Deal with Higgs decay to b-bbar     
      if(rem_dec.eqv..false.) then 
         s34=s(3,4)+2d0*mb**2
         hdecay=xn*gwsq*mbsq/(4d0*wmass**2)*2d0*(s34-4d0*mb**2)
         hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)
      else
         hdecay=1d0 
      endif

      fac=ason2pi*Asq*gsq*hdecay
      
      if (ss.ne.oldss.and.tt.ne.oldtt.and.uu.ne.olduu) then
         call hjetfill(ss,tt,uu,virtgg,virtqa,virtaq,virtqg,virtgq)
         oldss=ss
         oldtt=tt
         olduu=uu
      endif
     
      res=0d0

!----- identified only by initial two elements (constrains final state) 
!----- however want to ensure that a call of the form -5 0 4 etc is 
!----- not evaluated in this case check ret3 against ret1 and ret2 
!----- easist to check that one of ret1,ret2,ret3 = 0 
      if((ret1.ne.0).and.(ret2.ne.0).and.(ret3.ne.0)) then 
         return
      endif
      if(ret1+ret2-ret3.ne.0) then 
         return 
      endif
      
      if ((ret1.eq.0).and.(ret2.eq.0).and.(ret3.eq.0)) then 
         res=avegg*fac*virtgg
      elseif((ret1.gt.0).and.(ret2.eq.-ret1).and.(ret3.eq.0)) then  
         res=aveqq*fac*virtqa
      elseif ((ret1.lt.0).and.(ret2.eq.-ret1).and.(ret3.eq.0))then 
         res=aveqq*fac*virtaq
!         write(6,*) virtaq 
!         pause
      elseif ((ret1.eq.0).and.(ret2.ne.0).and.(ret3.eq.ret2))then 
         res=aveqg*fac*virtgq
      elseif ((ret1.ne.0).and.(ret2.eq.0).and.(ret3.eq.ret1))then  
         res=aveqg*fac*virtqg
      endif

      return
      end

!      do j=-nf,nf
!      do k=-nf,nf
!      msq(j,k)=0d0
!      if ((j.eq.0).and.(k.eq.0)) msq(j,k)=avegg*fac*virtgg
!      if ((j.gt.0).and.(k.eq.-j)) msq(j,k)=aveqq*fac*virtqa
!      if ((j.lt.0).and.(k.eq.-j)) msq(j,k)=aveqq*fac*virtaq
!      if ((j.eq.0).and.(k.ne.0)) msq(j,k)=aveqg*fac*virtgq
!      if ((j.ne.0).and.(k.eq.0)) msq(j,k)=aveqg*fac*virtqg
!      enddo
!      enddo
