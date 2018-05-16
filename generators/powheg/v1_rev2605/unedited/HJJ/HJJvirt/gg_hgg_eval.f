      subroutine gg_hgg_v_eval(p,res)
c--- Virtual matrix element squared averaged over initial colors and spins
c
c     g(-p1)+g(-p2) -->  H(p3)+g(p_iglue1=5)+g(p_iglue2=6) 
c
c    Calculation is fully analytic
!---- evaulation routine for i2MCFM 
      implicit none
      include 'MCFM_Include/constants.f'
      include 'MCFM_Include/masses.f'
      include 'MCFM_Include/ewcouple.f'
      include 'MCFM_Include/qcdcouple.f'
      include 'MCFM_Include/sprods_com.f'
      include 'MCFM_Include/zprods_com.f'
      include 'MCFM_Include/scheme.f'
      include 'MCFM_Include/nflav.f'
      include 'MCFM_Include/deltar.f'
      include 'MCFM_Include/interface_settings.f'
      include 'MCFM_Include/epinv.f'
      integer j,k,nu
      logical remdec
      integer maxstore,istore,icurr,imax
      parameter (maxstore=10)
      double precision p(mxpart,4),oldp(mxpart,4,maxstore),
     .  oldepinv(maxstore),res,s34
      double precision hdecay,Asq,fac
      double precision qrqr(maxstore),qarb(maxstore),aqbr(maxstore)
     . ,abab(maxstore),qbra(maxstore),bqar(maxstore)
      double precision qaqa(maxstore),aqaq(maxstore),qqqq(maxstore),
     .     aaaa(maxstore)
      double precision qagg(maxstore),aqgg(maxstore),qgqg(maxstore)
     . ,gqqg(maxstore),agag(maxstore),gaag(maxstore),ggqa(maxstore)
      double precision gggg(maxstore)
      double precision Hqarbvsqanal
      double precision Hqaqavsqanal
      double precision HAQggvsqanal
      double precision Hggggvsqanal
      logical CheckEGZ,recalc
      common/CheckEGZ/CheckEGZ
      data istore,imax,icurr/0,0,0/
      save oldp,qrqr,qarb,aqbr,abab,qbra,bqar,qaqa,aqaq,qqqq,aaaa,
     . qagg,aqgg,qgqg,gqqg,agag,gaag,ggqa,gggg,oldepinv,icurr,imax,
     . istore
c This below is a mechanism to avoid recalculating amplitudes.
c A list of 10 previously used moment, epinv values, and amplitude
c results is maintained. If the same p and epinv value is found
c in a call, the store values are used rather than recalculting them.
      do icurr=1,imax
         if(epinv.ne.oldepinv(icurr)) goto 33
         do j=1,6
            do nu=1,4
               if(p(j,nu).ne.oldp(j,nu,icurr)) then
                  goto 33
               endif
            enddo
         enddo
         recalc=.false.
         goto 44
 33      continue
      enddo
      recalc=.true.
      if(imax.lt.maxstore) then
         imax=imax+1
         oldp(:,:,imax)=p
         oldepinv(imax)=epinv
         icurr=imax
         if(imax.eq.maxstore) istore=1
      else
         oldp(:,:,istore)=p
         oldepinv(istore)=epinv
         icurr=istore
         if(istore.lt.maxstore) then
            istore=istore+1
         else
            istore=1
         endif
      endif
 44   continue
      
C*************************************************** 
!      scheme='tH-V'
!     allow users scheme 
      scheme = inscheme 
C*************************************************** 

      if     (scheme .eq. 'dred') then
        deltar=0d0
      elseif (scheme .eq. 'tH-V') then
        deltar=1d0
      else
        write(6,*) 'Invalid scheme in gg_hgg_v.f'
	stop
      endif
      
c--- Set this to true to check squared matrix elements against
c--- hep-ph/0506196 using the point specified in Eq. (51)
      CheckEGZ=.false.
            
c--- Set up spinor products
      if(recalc) call spinoru(6,p,za,zb)

      Asq=(as/(3d0*pi))**2/vevsq

      remdec=.true. 
C   Deal with Higgs decay to b-bbar
!---- Optional removal of decay 
      if(remdec.eqv..false.) then 
         s34=s(3,4)+2d0*mb**2
         hdecay=xn*gwsq*mbsq/(4d0*wmass**2)*2d0*(s34-4d0*mb**2)
         hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)
      else
         hdecay=1d0 
      endif

      fac=ason2pi*Asq*gsq**2*hdecay

c--- for checking EGZ
      if (CheckEGZ) then
        call CheckEGZres
      endif
      
c--- for checking scheme dependence of amplitudes
c      call CheckScheme(1,2,5,6)
      
C--- Note that Hqarbvsqanal(1,2,5,6)=Hqarbvsqanal(6,5,2,1)
C--- and the basic process is q(-k6)+r(-k2)-->q(-k5)+r(-k1)

c--- FOUR-QUARK PROCESSES WITH NON-IDENTICAL QUARKS

C---quark-quark
C     q(1)+r(2)->q(5)+r(6)
      if(recalc) qrqr(icurr)=Hqarbvsqanal(6,2,5,1)
      
C----quark-antiquark annihilation (6-->5-->2-->6) wrt q(1)+r(2)->q(5)+r(6)
c     q(1)+a(2)->r(5)+b(6)
      if(recalc) qarb(icurr)=Hqarbvsqanal(5,6,2,1)

C----antiquark-quark annihilation (1<-->2, 5<-->6) wrt to the above
c     a(1)+q(2)->b(5)+r(6)
c      aqbr=Hqarbvsqanal(6,5,1,2)
      aqbr=qarb
            
C----quark-antiquark scattering (6<-->2) wrt q(1)+r(2)->q(5)+r(6)
c     q(1)+b(2)->r(5)+a(6)
      if(recalc) qbra(icurr)=Hqarbvsqanal(2,6,5,1)

C----antiquark-quark scattering
c     b(1)+q(2)->a(5)+r(6) (1<-->2, 5<-->6) wrt to the above
c      bqar=Hqarbvsqanal(1,5,6,2) 
      bqar(icurr)=qbra(icurr)
      
C---antiquark-antiquark scattering (1<-->5,2<-->6) wrt q(1)+r(2)->q(5)+r(6)
C     a(1)+b(2)->a(5)+b(6)
      if(recalc) abab(icurr)=Hqarbvsqanal(2,6,1,5)

C--- FOUR-QUARK PROCESSES WITH IDENTICAL QUARKS

C     q(1)+q(2)->q(5)+q(6)
      if(recalc) qqqq(icurr)=qrqr(icurr)+Hqarbvsqanal(5,2,6,1)
     .+Hqaqavsqanal(6,2,5,1)

C     a(1)+a(2)->a(5)+a(6) (1<-->5,2<-->6) wrt q(1)+q(2)->q(5)+q(6)
      if(recalc) aaaa(icurr)=abab(icurr)+Hqarbvsqanal(2,5,1,6)
     .+Hqaqavsqanal(2,6,1,5)

C     q(1)+a(2)->q(5)+a(6) (2<-->6) wrt q(1)+q(2)->q(5)+q(6)
      if(recalc) qaqa(icurr)=qbra(icurr)+qarb(icurr)
     .+Hqaqavsqanal(2,6,5,1)

C     a(1)+q(2)->a(5)+q(6) (1<-->2, 5<-->6) wrt the above
C      aqqa=qbra+qarb+Hqaqavsqanal(1,5,6,2)
      aqaq(icurr)=qaqa(icurr)
      
c--- TWO-QUARK, TWO GLUON PROCESSES

C     a(1)+q(2)->g(3)+g(4)
      if(recalc) aqgg(icurr)=+HAQggvsqanal(2,1,5,6)

C     q(1)+g(2)->q(5)+g(6)
      if(recalc) qgqg(icurr)=+HAQggvsqanal(1,5,2,6)

C     g(1)+q(2)->q(5)+g(6)
      if(recalc) gqqg(icurr)=+HAQggvsqanal(2,5,1,6)

C     a(1)+g(2)->a(5)+g(6)
      if(recalc) agag(icurr)=qgqg(icurr)
C      if(recalc) agag(icurr)=+HAQggvsqanal(5,1,2,6)

C     g(1)+a(2)->a(5)+g(6)
      if (recalc) gaag(icurr)=gqqg(icurr)
C      if(recalc) gaag(icurr)=+HAQggvsqanal(5,2,1,6)

C     g(1)+g(2)->q(5)+a(6)
      if(recalc) ggqa(icurr)=+HAQggvsqanal(6,5,1,2)

C     q(1)+a(2)->g(5)+g(6)
      if (recalc) qagg(icurr)=aqgg(icurr)
C      if(recalc) qagg(icurr)=+HAQggvsqanal(1,2,5,6)
      
c--- FOUR GLUON PROCESS
      if(recalc) gggg(icurr)=+Hggggvsqanal(1,2,5,6)
      

C--- DEBUGGING OUTPUT
C      write(6,*) 'qrqr',qrqr
C      write(6,*) 'qarb',qarb
C      write(6,*) 'aqrb',aqrb
C      write(6,*) 'abab',abab
C      write(6,*) 'qbra',qbra
C      write(6,*) 'bqra',bqra

C      write(6,*) 'Identical'
C      write(6,*) 'qaqa',qaqa
C      write(6,*) 'aqqa',aqqa
C      write(6,*) 'qqqq',qqqq
C      write(6,*) 'aaaa',aaaa

      
      if((ret1.eq.0).and.(ret2.eq.0)) then 
         if((ret3.eq.0).and.(ret4.eq.0)) then 
            res=fac*avegg*half*gggg(icurr)
         elseif((ret3.gt.0).and.(ret4.lt.0)) then 
            res=fac*avegg*ggqa(icurr) 
         else 
            res=0d0 
         endif
         return 
      elseif((ret1.gt.0).and.(ret2.gt.0)) then 
         if((ret1.eq.ret2).and.(ret3.eq.ret1).and.(ret4.eq.ret2)) then 
            res=aveqq*fac*half*qqqq(icurr)
        elseif((ret1.ne.ret2).and.(ret3.eq.ret1).and.(ret4.eq.ret2))then 
           res=aveqq*fac*qrqr(icurr)
        else
           res=0d0 
        endif
        return 
      elseif((ret1.lt.0).and.(ret2.lt.0)) then 
         if((ret1.eq.ret2).and.(ret3.eq.ret1).and.(ret4.eq.ret2)) then 
            res=aveqq*fac*half*aaaa(icurr)
        elseif((ret1.ne.ret2).and.(ret3.eq.ret1).and.(ret4.eq.ret2))then 
            res=aveqq*fac*abab(icurr)
         else
           res=0d0 
        endif
        return 
      elseif((ret1.gt.0).and.(ret2.lt.0)) then 
         if((ret1.eq.-ret2).and.((ret3.eq.0).and.(ret4.eq.0)))then 
            res=aveqq*fac*half*qagg(icurr)
         elseif((ret1.eq.-ret2)
     &           .and.((ret3.eq.ret1).and.(ret4.eq.ret2)))then
            res=aveqq*qaqa(icurr)*fac 
         elseif((ret1.eq.-ret2).and.(abs(ret3).ne.abs(ret1))
     &           .and.(ret3.gt.ret4)) then 
            res=aveqq*fac*qarb(icurr) 
         elseif((ret1.ne.-ret2).and.(ret3.gt.ret4)) then 
            res=aveqq*fac*qbra(icurr)
         else
            res=0d0
         endif
         return 
      elseif((ret1.lt.0).and.(ret2.gt.0)) then 
         if((ret1.eq.-ret2).and.((ret3.eq.0).and.(ret4.eq.0)))then 
            res=aveqq*fac*half*aqgg(icurr)
         elseif((ret1.eq.-ret2)
     &           .and.((ret3.eq.ret1).and.(ret4.eq.ret2)))then
            res=aveqq*aqaq(icurr)*fac 
         elseif((ret1.eq.-ret2).and.(abs(ret3).ne.abs(ret1))
     &           .and.(ret3.gt.ret4)) then 
            res=aveqq*fac*aqbr(icurr)
         elseif((ret1.ne.-ret2).and.(ret3.lt.ret4)) then 
            res=aveqq*fac*bqar(icurr) 
         else
            res=0d0
         endif
         return  
      elseif((ret1.eq.0).and.(ret2.gt.0)) then 
         if(ret4.eq.0) then 
            res=aveqg*fac*gqqg(icurr) 
         else
            res=0d0 
         endif
      elseif((ret1.eq.0).and.(ret2.lt.0)) then 
         if(ret4.eq.0) then 
            res=aveqg*fac*gaag(icurr) 
         else
            res=0d0 
         endif
      elseif((ret1.gt.0).and.(ret2.eq.0)) then 
         if(ret4.eq.0) then 
            res=aveqg*fac*qgqg(icurr) 
         else
            res=0d0 
         endif
      elseif((ret1.lt.0).and.(ret2.eq.0)) then 
         if(ret4.eq.0) then 
            res=aveqg*fac*agag(icurr) 
         else
            res=0d0 
         endif
      endif

      return
      end 




  !  do j=-nf,nf
    !  do k=-nf,nf
    !  msq(j,k)=0d0

    !  if ((j.eq.0).and.(k.eq.0)) then
C---gg - all poles cancelled
    !     msq(j,k)=fac*avegg*(half*gggg+dfloat(nflav)*ggqa)

!      elseif ((j.gt.0).and.(k.gt.0)) then
C---qq - all poles cancelled
!         if (j.eq.k) then
!         msq(j,k)=aveqq*fac*half*qqqq
!         else
!         msq(j,k)=aveqq*fac*qrqr
!         endif


!      elseif ((j.lt.0).and.(k.lt.0)) then
C---aa - all poles cancelled
!         if (j.eq.k) then
!         msq(j,k)=aveqq*fac*half*aaaa
!         else
!         msq(j,k)=aveqq*fac*abab
!         endif


!      elseif ((j.gt.0).and.(k.lt.0)) then
C----qa scattering - all poles cancelled
!         if (j.eq.-k) then
!         msq(j,k)=aveqq*fac*(dfloat(nflav-1)*qarb+qaqa+half*qagg)
!             else
!         msq(j,k)=aveqq*fac*qbra
!         endif


!      elseif ((j.lt.0).and.(k.gt.0)) then
C----aq scattering - all poles cancelled
!         if (j.eq.-k) then
!         msq(j,k)=aveqq*fac*(dfloat(nflav-1)*aqbr+aqaq+half*aqgg)
!             else
!         msq(j,k)=aveqq*fac*bqar
!         endif

!     elseif ((j.eq.0).and.(k.gt.0)) then
C----gq scattering - all poles cancelled
!         msq(j,k)=aveqg*fac*gqqg

!      elseif ((j.eq.0).and.(k.lt.0)) then
C----ga scattering - all poles cancelled
!         msq(j,k)=aveqg*fac*gaag

!     elseif ((j.gt.0).and.(k.eq.0)) then
C----qg scattering - all poles cancelled
!         msq(j,k)=aveqg*fac*qgqg

!      elseif ((j.lt.0).and.(k.eq.0)) then
C----ag scattering - all poles cancelled
!         msq(j,k)=aveqg*fac*agag
!      endif

!      enddo
!      enddo
