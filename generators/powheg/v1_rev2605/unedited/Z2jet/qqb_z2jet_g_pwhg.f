      subroutine qqb_z2jet_g_pwhg_nc(p,rflav,jperm,res)
************************************************************************
*     Author: R.K.Ellis                                                *
*     July, 2012.                                                      *
c---  Matrix element squared averaged over initial colors and spins    *
c     q(-p1)+qbar(-p2) -->  Z^0 + f(p5)+f(p6)+g(p7)                    *
c                           |                                          *
c                            --> e^-(p3)+e^+(p4)                       *
c     Routine constructed to work with powheg                          *
c     This routine only works:-                                        *
c     for the (idiosyncratic) ordering in rflav dictated by mcfm       *
************************************************************************
      implicit none 
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'zcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'flags.f'
      include 'lc.f'
      include 'scale.f'
      integer j,k,nquark,rflav(7),nq,i5,nu,i
      integer maxperm,jperm
      parameter (maxperm=6)
      integer, save :: jj(-nf:nf),kk(-nf:nf)
      double precision P(mxpart,4),msq(-nf:nf,-nf:nf),res,fac,ggtemp
      double precision, save:: mmsq_gg(2,2,maxperm),
     & mmsq_qqb(2,2,maxperm),mmsq_qbq(2,2,maxperm),
     & mmsq_qg(2,2,maxperm),mmsq_gq(2,2,maxperm),
     & mmsq_gqb(2,2,maxperm),mmsq_qbg(2,2,maxperm),
     & msqi_qq(2,maxperm),msqi_qbqb(2,maxperm),
     & msqn_qq(2,2,maxperm),msqn_qbqb(2,2,maxperm),
     & msqi_qqb(2,maxperm),msqi_qbq(2,maxperm),
     & msqn_qqb(2,2,maxperm),msqn_qbq(2,2,maxperm),
     & msqi_qqbs(2,maxperm),msqi_qbqs(2,maxperm),
     & msqn_qqbs(2,2,maxperm),msqn_qbqs(2,2,maxperm),
     & msqi_qg(2,maxperm),msqi_qbg(2,maxperm),
     & msqn_qg(2,2,maxperm),msqn_qbg(2,2,maxperm),
     & msqi_gqb(2,maxperm),msqi_gq(2,maxperm),
     & msqn_gqb(2,2,maxperm),msqn_gq(2,2,maxperm)
      double complex, save :: prop
      logical first,qa,aq,qq,aa,qg,ag,gq,ga,gg

      data jj/-1,-2,-1,-2,-1,0,1,2,1,2,1/
      data kk/-1,-2,-1,-2,-1,0,1,2,1,2,1/
      data first/.true./
      save first
C     variables needed to avoid recalculating same stuff
C     variables needed to avoid recalculating same stuff
      integer, parameter :: nchn = 9 
      logical ::  recalc_g(nchn), recalc_q(nchn)
      double precision, save :: opin_g(mxpart,4,nchn,maxperm)
      double precision, save :: opin_q(mxpart,4,nchn,maxperm)
      integer, parameter :: iqa = 1 
      integer, parameter :: iaq = 2 
      integer, parameter :: iqg = 3 
      integer, parameter :: igq = 4 
      integer, parameter :: iag = 5 
      integer, parameter :: iga = 6 
      integer, parameter :: iqq = 7 
      integer, parameter :: iaa = 8 
      integer, parameter :: igg = 9 
      integer :: ichn 

C------first determine if it is a two- or four-quark process 
      nq = 0 
      do j=1,7
         if (j < 3 .or. j > 4) then 
            if (abs(rflav(j)) > 0) then 
               nq = nq+1
            endif
         endif
      enddo
      if (nq == 2) then 
         Gflag = .true. 
         Qflag = .false. 
      elseif (nq == 4) then 
         Qflag = .true. 
         Gflag = .false. 
      else 
         write(*,*) 'nq',nq 
         write(*,*) 'rflav',rflav
         stop 'nq out of range' 
      endif
      colourchoice=0


      j=rflav(1)
      k=rflav(2)
      i5=rflav(5)

      qa=(j > 0) .and. (k < 0)
      aq=(j < 0) .and. (k > 0) 
      qq=(j > 0) .and. (k > 0) 
      aa=(j < 0) .and. (k < 0) 
      qg=(j > 0) .and. (k == 0) 
      ag=(j < 0) .and. (k == 0) 
      gq=(j == 0) .and. (k > 0) 
      ga=(j == 0) .and. (k < 0) 
      gg=(j == 0) .and. (k == 0) 

      if (qa) then
         ichn = iqa
      elseif (aq) then
         ichn = iaq
      elseif (qg) then
         ichn = iqg
      elseif (gq) then
         ichn = igq
      elseif (ag) then
         ichn = iag
      elseif (ga) then
         ichn = iga
      elseif (qq) then
         ichn = iqq
      elseif (aa) then
         ichn = iaa
      elseif (gg) then
         ichn = igg
      else
         stop 'qqb_z2jet_g_pwhg: channel not recognized'
      endif


!     decide if recalc                                                                
      recalc_q = .false.
      recalc_g = .false.


      if (Gflag) then
         do nu=1,4
            do i=1,7
               if(opin_g(i,nu,ichn,jperm).ne.p(i,nu)) then
                  recalc_g(ichn) = .true.
                  goto 10
               endif
            enddo
         enddo
      else
         do nu=1,4
            do i=1,7
               if(opin_q(i,nu,ichn,jperm).ne.p(i,nu)) then
                  recalc_q(ichn) = .true.
                  goto 10
               endif
            enddo
         enddo

      endif

 10       continue

       if (recalc_g(ichn)) then
          opin_g(:,:,ichn,jperm) = p
       elseif (recalc_q(ichn)) then
          opin_q(:,:,ichn,jperm) = p
       endif

      msq(:,:)=0d0



      if (Gflag) then
************************************************************************
*     Calculate contributions from the QQGGG matrix elements            *
************************************************************************

            if (recalc_g(ichn)) then
               call spinoru(7,p,za,zb)
               prop=s(3,4)/dcmplx((s(3,4)-zmass**2),zmass*zwidth)
             if (recalc_g(igg)) then 
c--            matrix elements for gg -> qbq      
               call xzqqggg(5,1,2,7,6,3,4,mmsq_gg(:,:,jperm))
             elseif (recalc_g(iqa)) then 
c--            matrix elements for qqb -> gg      
               call xzqqggg(1,5,6,7,2,3,4,mmsq_qqb(:,:,jperm))
             elseif (recalc_g(iaq)) then 
c--            matrix elements for qbq -> gg      
               call xzqqggg(2,5,6,7,1,3,4,mmsq_qbq(:,:,jperm))
             elseif (recalc_g(iqg)) then 
c--            matrix elements for qg -> qg      
               call xzqqggg(1,2,6,7,5,3,4,mmsq_qg(:,:,jperm))
             elseif (recalc_g(igq)) then 
c--            matrix elements for gq -> gq      
               call xzqqggg(2,1,6,7,5,3,4,mmsq_gq(:,:,jperm))
             elseif (recalc_g(iga)) then 
c--            matrix elements for gqb -> gqb      
               call xzqqggg(5,1,6,7,2,3,4,mmsq_gqb(:,:,jperm))
             elseif (recalc_g(iag)) then 
c--            matrix elements for qbg -> qbg      
               call xzqqggg(5,2,6,7,1,3,4,mmsq_qbg(:,:,jperm))
             endif 
             endif 


      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

            if (gg) then
               ggtemp=0d0
               nquark=rflav(6)
               ggtemp=ggtemp
     &       +abs(Q(nquark)*q1+L(nquark)*l1*prop)**2*mmsq_gg(1,1,jperm)
     &       +abs(Q(nquark)*q1+R(nquark)*r1*prop)**2*mmsq_gg(2,2,jperm)
     &       +abs(Q(nquark)*q1+L(nquark)*r1*prop)**2*mmsq_gg(1,2,jperm)
     &       +abs(Q(nquark)*q1+R(nquark)*l1*prop)**2*mmsq_gg(2,1,jperm)
               msq(j,k)=ggtemp

            elseif (qa) then
               msq(j,k)=abs(Q(j)*q1+L(j)*l1*prop)**2*mmsq_qqb(1,1,jperm)
     &                 +abs(Q(j)*q1+R(j)*r1*prop)**2*mmsq_qqb(2,2,jperm)
     &                 +abs(Q(j)*q1+L(j)*r1*prop)**2*mmsq_qqb(1,2,jperm)
     &                 +abs(Q(j)*q1+R(j)*l1*prop)**2*mmsq_qqb(2,1,jperm)
               msq(j,k)=aveqq/avegg*msq(j,k)/6d0

            elseif (aq) then
               msq(j,k)=abs(Q(k)*q1+L(k)*l1*prop)**2*mmsq_qbq(1,1,jperm)
     &             +abs(Q(k)*q1+R(k)*r1*prop)**2*mmsq_qbq(2,2,jperm)
     &             +abs(Q(k)*q1+L(k)*r1*prop)**2*mmsq_qbq(1,2,jperm)
     &             +abs(Q(k)*q1+R(k)*l1*prop)**2*mmsq_qbq(2,1,jperm)
               msq(j,k)=aveqq/avegg*msq(j,k)/6d0

            elseif (qg) then
               msq(j,k)=+abs(Q(j)*q1+L(j)*l1*prop)**2*mmsq_qg(1,1,jperm)
     &                  +abs(Q(j)*q1+R(j)*r1*prop)**2*mmsq_qg(2,2,jperm)
     &                  +abs(Q(j)*q1+L(j)*r1*prop)**2*mmsq_qg(1,2,jperm)
     &                  +abs(Q(j)*q1+R(j)*l1*prop)**2*mmsq_qg(2,1,jperm)
                msq(j,k)=half*aveqg/avegg*msq(j,k)

            elseif (ag) then
                msq(j,k)=
     &               abs(Q(-j)*q1+L(-j)*l1*prop)**2*mmsq_qbg(1,1,jperm)
     &              +abs(Q(-j)*q1+R(-j)*r1*prop)**2*mmsq_qbg(2,2,jperm)
     &              +abs(Q(-j)*q1+L(-j)*r1*prop)**2*mmsq_qbg(1,2,jperm)
     &              +abs(Q(-j)*q1+R(-j)*l1*prop)**2*mmsq_qbg(2,1,jperm)
                msq(j,k)=half*aveqg/avegg*msq(j,k)

            elseif (gq) then
                msq(j,k)=abs(Q(k)*q1+L(k)*l1*prop)**2*mmsq_gq(1,1,jperm)
     &                  +abs(Q(k)*q1+R(k)*r1*prop)**2*mmsq_gq(2,2,jperm)
     &                  +abs(Q(k)*q1+L(k)*r1*prop)**2*mmsq_gq(1,2,jperm)
     &                  +abs(Q(k)*q1+R(k)*l1*prop)**2*mmsq_gq(2,1,jperm)
                msq(j,k)=half*aveqg/avegg*msq(j,k)

            elseif (ga) then
                msq(j,k)=
     &                abs(Q(-k)*q1+L(-k)*l1*prop)**2*mmsq_gqb(1,1,jperm)
     &               +abs(Q(-k)*q1+R(-k)*r1*prop)**2*mmsq_gqb(2,2,jperm)
     &               +abs(Q(-k)*q1+L(-k)*r1*prop)**2*mmsq_gqb(1,2,jperm)
     &               +abs(Q(-k)*q1+R(-k)*l1*prop)**2*mmsq_gqb(2,1,jperm)
                msq(j,k)=half*aveqg/avegg*msq(j,k)
            endif

   19 continue
      res=msq(j,k)
      return
      endif ! end Gflag


      if (Qflag) then
c--- note the factor of 4d0*xw**2 relative to wbb
      fac=4d0*gsq**3*esq**2*8d0
c--- extra factor of 2**3=8 to compensate for Ta normalization

c--- note: the following two arrays end up being overall 1<->2 symmetric
         if (recalc_q(ichn)) then
            call spinoru(7,p,za,zb)
            prop=s(3,4)/dcmplx((s(3,4)-zmass**2),zmass*zwidth)
            if (recalc_q(iqa)) then
C           qa
               call msq_ZqqQQg(5,1,2,6,7,4,3,msqn_qqbs(:,:,jperm),
     &              msqi_qqbs(:,jperm))
               call msq_ZqqQQg(2,1,5,6,7,4,3,msqn_qqb(:,:,jperm),
     &              msqi_qqb(:,jperm))
            elseif (recalc_q(iaq)) then
C           aq
               call msq_ZqqQQg(1,6,5,2,7,4,3,msqn_qbqs(:,:,jperm),
     &              msqi_qbqs(:,jperm))
               call msq_ZqqQQg(1,2,5,6,7,4,3,msqn_qbq(:,:,jperm),
     &              msqi_qbq(:,jperm))
            elseif (recalc_q(iaa)) then
C           aa
               call msq_ZqqQQg(1,5,2,6,7,4,3,msqn_qbqb(:,:,jperm),
     &              msqi_qbqb(:,jperm))
            elseif (recalc_q(iqq)) then
C           qq
            call msq_ZqqQQg(5,1,6,2,7,4,3,msqn_qq(:,:,jperm),
     &              msqi_qq(:,jperm))
            elseif (recalc_q(iqg)) then
C           qg
               call msq_ZqqQQg(7,1,5,6,2,4,3,msqn_qg(:,:,jperm),
     &              msqi_qg(:,jperm))
            elseif (recalc_q(iga)) then
C           ga
               call msq_ZqqQQg(2,7,5,6,1,4,3,msqn_gqb(:,:,jperm),
     &              msqi_gqb(:,jperm))
            elseif (recalc_q(igq)) then
C           gq
               call msq_ZqqQQg(7,2,5,6,1,4,3,msqn_gq(:,:,jperm),
     &              msqi_gq(:,jperm))
            elseif (recalc_q(iag)) then
C           ag
               call msq_ZqqQQg(1,7,5,6,2,4,3,msqn_qbg(:,:,jperm),
     &              msqi_qbg(:,jperm))
            endif 
         endif                  !end recalc_q
       


         if (qa) then
              if (k.eq.-j) then
                 if (i5 == j) then
                   msq(j,k)=fac*aveqq*msqi_qqbs(jj(j),jperm)
                 else
                   msq(j,k)=fac*aveqq*msqn_qqb(jj(j),jj(i5),jperm)
                 endif
              else
                   msq(j,k)=fac*aveqq*msqn_qqbs(jj(j),-kk(k),jperm)
              endif
         elseif (aq) then
             if (j .eq.-k) then
                if (i5 == k) then
                   msq(j,k)=fac*aveqq*msqi_qbqs(kk(k),jperm)
                else
                   msq(j,k)=fac*aveqq*msqn_qbq(kk(k),kk(i5),jperm)
                endif
             else
                msq(j,k)=fac*aveqq*msqn_qbqs(-jj(j),kk(k),jperm)
             endif
         elseif (qg) then
            if (i5 == j) then  
                msq(j,k)=fac*aveqg*half*msqi_qg(jj(j),jperm)
            else
                msq(j,k)=fac*aveqg*msqn_qg(jj(j),jj(i5),jperm)
            endif
         elseif (ag) then
            if (i5 == -j) then
               msq(j,k)=fac*aveqg*half*msqi_qbg(-jj(j),jperm)
            else 
               msq(j,k)=fac*aveqg*msqn_qbg(-jj(j),jj(i5),jperm)
            endif
         elseif (gq) then
            if (i5 == k) then
               msq(j,k)=fac*aveqg*half*msqi_gq(kk(k),jperm)
            else
               msq(j,k)=fac*aveqg*msqn_gq(kk(k),kk(i5),jperm)
            endif
         elseif (ga) then
            if (i5 == -k) then
               msq(j,k)=fac*aveqg*half*msqi_gqb(-kk(k),jperm)
            else
               msq(j,k)=fac*aveqg*msqn_gqb(-kk(k),kk(i5),jperm)
            endif
         elseif (qq) then
c-qq
            if (j.eq.k) then
               msq(j,k)=half*fac*aveqq*msqi_qq(jj(j),jperm)
            else
               msq(j,k)=fac*aveqq*msqn_qq(jj(j),kk(k),jperm)
            endif
         elseif (aa) then
c-qbqb
            if (j.eq.k) then
               msq(j,k)=msq(j,k)+half*fac*aveqq*msqi_qbqb(-jj(j),jperm)
            else
               msq(j,k)=msq(j,k)
     &              +fac*aveqq*msqn_qbqb(-jj(j),-kk(k),jperm)
            endif
         endif

      res=msq(j,k)

      endif ! end Qflag




      return
      end
