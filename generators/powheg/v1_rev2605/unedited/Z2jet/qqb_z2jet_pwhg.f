      subroutine qqb_z2jet_pwhg(pin,bflav,res)
      implicit none
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) --> Z +g(p5) +g(p6)
c                          |
c                          --> l(p3)+a(p4)
c                            
c--all momenta incoming
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'msq_cs.f'
      include 'flags.f'
      include 'Bmcfm.f'
      include 'qq_cs.f'
      integer i,j,k,pq,pl,nquark,swap(2),swap1(0:2),nup,ndo,
     . j1,j2,j3,icol,nu
      double precision res,msq(-nf:nf,-nf:nf),p(mxpart,4),fac,faclo,
     & pin(mxpart,4),ggtemp(0:2)
      double precision, save ::
     .   BqqbZgg2(4,4,6,2,2),BqgZqg2(4,4,6,2,2),
     .   BqbqZgg2(4,4,6,2,2),BqbgZqbg2(4,4,6,2,2),BgqbZqbg2(4,4,6,2,2),
     .   BgqZqg2(4,4,6,2,2),BggZqbq2(4,4,6,2,2),
     .   qqbZgg2_cs(0:2,2,2),qbqZgg2_cs(0:2,2,2),
     .   qgZqg2_cs(0:2,2,2),gqZqg2_cs(0:2,2,2),
     .   qbgZqbg2_cs(0:2,2,2),gqbZqbg2_cs(0:2,2,2),
     .   ggZqbq2_cs(0:2,2,2)
      double precision
     .   qqbZgg2(2,2),qgZqg2(2,2),
     .   qbqZgg2(2,2),qbgZqbg2(2,2),gqbZqbg2(2,2),
     .   gqZqg2(2,2),ggZqbq2(2,2)
      double precision tup,tdo
      double complex a111,a112,a121,a211,a122,a212,a221,a222
      double complex b111,b112,b121,b211,b122,b212,b221,b222

      double complex, save :: qRb_a(2,2,2),qRb_b(2,2,2)
      double complex, save :: qqb_a(2,2,2),qqb_b(2,2,2)

      double complex, save :: qbq_a(2,2,2),qbq_b(2,2,2)
      double complex, save :: qbR_a(2,2,2),qbR_b(2,2,2)

      double complex, save ::qq_a(2,2,2),qq_b(2,2,2)
      double complex, save :: qR_a(2,2,2),qR_b(2,2,2)

      double complex, save :: qbRb_a(2,2,2),qbRb_b(2,2,2)
      double complex, save :: qbqb_a(2,2,2),qbqb_b(2,2,2)
      double complex, save :: prop

      double complex zamub(mxpart,4,mxpart)
      common/zamub/zamub
      logical ggaq,qgqg,gqqg,agag,gaag,qagg,aqgg,qqqq,aaaa,qaqa,aqqa
      
      data swap/2,1/
      save swap
      data swap1/0,2,1/
      save swap1
      integer bflav(6),nq,nflavmin,nflavmax
      logical same 
      double precision tmp(4) 

C     variables needed to avoid recalculating same stuff 
      logical ::  recalc_g, recalc_q 
      real * 8, save :: opin_g(mxpart,4)
      real * 8, save :: opin_q(mxpart,4)

      res=0d0
C----zero out Bmcfm which may get refilled
      Bmcfm(:,:,:)=0d0

!     Set up momenta 
      p = pin 

      ! now understand if it's a two-quark or four-quark process 
      nq = 0 
      do i=1,6
         if (i < 3 .or. i > 4) then 
            if (abs(bflav(i)) .gt. 0) then 
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
         write(*,*) 'bflav',bflav
         stop 'nq out of range' 
      endif

      recalc_g = .false.
      recalc_q = .false.
      if (Gflag) then
            do i=1,6
               do nu=1,4
                  if(opin_g(i,nu).ne.p(i,nu)) then
                     recalc_g = .true.
                     goto 10
                  endif
               enddo
            enddo
      else 
            do i=1,6
               do nu=1,4
                  if(opin_q(i,nu).ne.p(i,nu)) then
                     recalc_q = .true.
                     goto 10
                  endif
               enddo
            enddo
      endif

 10    continue 
       
       if (recalc_g) then
          opin_g = p
       endif
       if (recalc_q) then
          opin_q = p
       endif  
      

      same = .false. 

!     Gflag case 
      if (bflav(1) == 0 .and. bflav(2) == 0) then 
         if (abs(bflav(5)) == 2 .or. abs(bflav(5)) == 4) then 
            nup = 1
            ndo = 0 
            nflavmin = 2 
            nflavmax = 2
         else
            nup = 0
            ndo = 1
            nflavmin = 1 
            nflavmax = 1
         endif
!     Qflav 
!      qqb -> qqb 
      elseif ((bflav(1).eq.-bflav(2)).and.(bflav(5).eq.-bflav(6))) then 
         
         if (abs(bflav(1)) .eq. abs(bflav(5))) then 
            same = .true. 
         else
            same = .false. 
         endif

            if ((bflav(5)/2)*2 == bflav(5)) then 
            nup = 1
            ndo = 0 
!            same = .false. 
         else
            nup = 0
            ndo = 1
!            same = .false. 
         endif



!     qq -> qq or qbqb -> qbqb 
      elseif ((bflav(1).eq.bflav(2)).and.(bflav(5).eq.bflav(6))) then 
         same = .true. 
      endif
      
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
      

      j=bflav(1)
      k=bflav(2)

      qagg=(j > 0) .and. (k < 0)
      qaqa=(j > 0) .and. (k < 0)

      aqgg=(j < 0) .and. (k > 0) 
      aqqa=(j < 0) .and. (k > 0)

      qgqg=(j > 0) .and. (k == 0)
      agag=(j < 0) .and. (k == 0)

      gqqg=(j == 0) .and. (k > 0)
      gaag=(j == 0) .and. (k < 0)
      ggaq=(j == 0) .and. (k == 0)

      qqqq=(j > 0) .and. (k > 0)
      aaaa=(j < 0) .and. (k < 0)

      fac=v*xn/four*(esq*gsq)**2
c--- calculate 2-quark, 2-gluon amplitudes
      if (Gflag) then
      
      if (recalc_g) then
      call spinoru(6,p,za,zb)
      prop=s(3,4)/dcmplx((s(3,4)-zmass**2),zmass*zwidth)
C     fill <a|\gamma^\mu|b] needed for Bmunu calculation
      call spinormu(6,p,zamub)

        call z2jetsq(1,2,3,4,5,6,za,zb,qqbZgg2)
        call storecsz(qqbZgg2_cs)
        call fillZBmunu(1,2,3,4,5,6,p,BqqbZgg2)
        qqbZgg2_cs(:,:,:) = half*aveqq*qqbZgg2_cs(:,:,:)
        BqqbZgg2=half*aveqq*BqqbZgg2

        call z2jetsq(1,5,3,4,2,6,za,zb,qgZqg2)
        call storecsz(qgZqg2_cs)
        call fillZBmunu(1,5,3,4,2,6,p,BqgZqg2)
        qgZqg2_cs(:,:,:)  = aveqg*qgZqg2_cs(:,:,:)
        BqgZqg2=aveqg*BqgZqg2

        call z2jetsq(2,5,3,4,1,6,za,zb,gqZqg2)
        call storecsz(gqZqg2_cs)
        call fillZBmunu(2,5,3,4,1,6,p,BgqZqg2)
        gqZqg2_cs(:,:,:)  = aveqg*gqZqg2_cs(:,:,:)
        BgqZqg2=aveqg*BgqZqg2

        call z2jetsq(2,1,3,4,5,6,za,zb,qbqZgg2)
        call storecsz(qbqZgg2_cs)
        call fillZBmunu(2,1,3,4,5,6,p,BqbqZgg2)
        qbqZgg2_cs(:,:,:) = half*aveqq*qbqZgg2_cs(:,:,:)
        BqbqZgg2=half*aveqq*BqbqZgg2

        call z2jetsq(5,1,3,4,2,6,za,zb,qbgZqbg2)
        call storecsz(qbgZqbg2_cs)
        call fillZBmunu(5,1,3,4,2,6,p,BqbgZqbg2)
        qbgZqbg2_cs(:,:,:)= aveqg*qbgZqbg2_cs(:,:,:)
        BqbgZqbg2=aveqg*BqbgZqbg2

        call z2jetsq(5,2,3,4,1,6,za,zb,gqbZqbg2)
        call storecsz(gqbZqbg2_cs)
        call fillZBmunu(5,2,3,4,1,6,p,BgqbZqbg2)
        gqbZqbg2_cs(:,:,:)= aveqg*gqbZqbg2_cs(:,:,:)
        BgqbZqbg2=aveqg*BgqbZqbg2

C --NB this is the matrix element for gg->Z qb(5) q(6)
        call z2jetsq(5,6,3,4,1,2,za,zb,ggZqbq2)
        call storecsz(ggZqbq2_cs)        
        call fillZBmunu(5,6,3,4,1,2,p,BggZqbq2)
        ggZqbq2_cs(:,:,:) = avegg*ggZqbq2_cs(:,:,:) 
        BggZqbq2=avegg*BggZqbq2

      endif

      do icol=0,2
      msq_cs(icol,j,k)=zip
      enddo
      
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

      if     (ggaq) then

          do icol=0,2
          ggtemp(icol)=0d0
          do nquark=nflavmin,nflavmax
           ggtemp(icol)=ggtemp(icol)
     .      +abs(Q(nquark)*q1+L(nquark)*l1*prop)**2*ggZqbq2_cs(icol,1,1)
     .      +abs(Q(nquark)*q1+R(nquark)*r1*prop)**2*ggZqbq2_cs(icol,2,2)
     .      +abs(Q(nquark)*q1+L(nquark)*r1*prop)**2*ggZqbq2_cs(icol,1,2)
     .      +abs(Q(nquark)*q1+R(nquark)*l1*prop)**2*ggZqbq2_cs(icol,2,1)
          enddo
          msq_cs(icol,j,k)=fac*ggtemp(icol)
          enddo
          Bmcfm(:,:,:)=0d0
          do nquark=nflavmin,nflavmax
             Bmcfm(:,:,:)=Bmcfm(:,:,:)
     .       +abs(Q(nquark)*q1+L(nquark)*l1*prop)**2*BggZqbq2(:,:,:,1,1)
     .       +abs(Q(nquark)*q1+R(nquark)*r1*prop)**2*BggZqbq2(:,:,:,2,2)
     .       +abs(Q(nquark)*q1+L(nquark)*r1*prop)**2*BggZqbq2(:,:,:,1,2)
     .       +abs(Q(nquark)*q1+R(nquark)*l1*prop)**2*BggZqbq2(:,:,:,2,1)
          enddo
          Bmcfm(:,:,:)=fac*Bmcfm(:,:,:)
       elseif (qagg) then
          do icol=0,2
             msq_cs(icol,j,k)=fac*(
     .       +abs(Q(j)*q1+L(j)*l1*prop)**2*qqbZgg2_cs(icol,1,1)
     .       +abs(Q(j)*q1+R(j)*r1*prop)**2*qqbZgg2_cs(icol,2,2)
     .       +abs(Q(j)*q1+L(j)*r1*prop)**2*qqbZgg2_cs(icol,1,2)
     .       +abs(Q(j)*q1+R(j)*l1*prop)**2*qqbZgg2_cs(icol,2,1))
          enddo
             Bmcfm(:,:,:)=fac*(
     .       +abs(Q(j)*q1+L(j)*l1*prop)**2*BqqbZgg2(:,:,:,1,1)
     .       +abs(Q(j)*q1+R(j)*r1*prop)**2*BqqbZgg2(:,:,:,2,2)
     .       +abs(Q(j)*q1+L(j)*r1*prop)**2*BqqbZgg2(:,:,:,1,2)
     .       +abs(Q(j)*q1+R(j)*l1*prop)**2*BqqbZgg2(:,:,:,2,1))
c---Statistical factor already included above
      elseif (aqgg) then
          do icol=0,2
             msq_cs(icol,j,k)=fac*(
     .       +abs(Q(k)*q1+L(k)*l1*prop)**2*qbqZgg2_cs(icol,1,1)
     .       +abs(Q(k)*q1+R(k)*r1*prop)**2*qbqZgg2_cs(icol,2,2)
     .       +abs(Q(k)*q1+L(k)*r1*prop)**2*qbqZgg2_cs(icol,1,2)
     .       +abs(Q(k)*q1+R(k)*l1*prop)**2*qbqZgg2_cs(icol,2,1))
          enddo
             Bmcfm(:,:,:)=fac*(
     .       +abs(Q(k)*q1+L(k)*l1*prop)**2*BqbqZgg2(:,:,:,1,1)
     .       +abs(Q(k)*q1+R(k)*r1*prop)**2*BqbqZgg2(:,:,:,2,2)
     .       +abs(Q(k)*q1+L(k)*r1*prop)**2*BqbqZgg2(:,:,:,1,2)
     .       +abs(Q(k)*q1+R(k)*l1*prop)**2*BqbqZgg2(:,:,:,2,1))
      elseif (qgqg) then
          do icol=0,2
             msq_cs(icol,j,k)=fac*(
     .       +abs(Q(j)*q1+L(j)*l1*prop)**2*qgZqg2_cs(icol,1,1)
     .       +abs(Q(j)*q1+R(j)*r1*prop)**2*qgZqg2_cs(icol,2,2)
     .       +abs(Q(j)*q1+L(j)*r1*prop)**2*qgZqg2_cs(icol,1,2)
     .       +abs(Q(j)*q1+R(j)*l1*prop)**2*qgZqg2_cs(icol,2,1))
          enddo
             Bmcfm(:,:,:)=fac*(
     .       +abs(Q(j)*q1+L(j)*l1*prop)**2*BqgZqg2(:,:,:,1,1)
     .       +abs(Q(j)*q1+R(j)*r1*prop)**2*BqgZqg2(:,:,:,2,2)
     .       +abs(Q(j)*q1+L(j)*r1*prop)**2*BqgZqg2(:,:,:,1,2)
     .       +abs(Q(j)*q1+R(j)*l1*prop)**2*BqgZqg2(:,:,:,2,1))
      elseif (agag) then
          do icol=0,2
             msq_cs(icol,j,k)=fac*(
     .       +abs(Q(-j)*q1+L(-j)*l1*prop)**2*qbgZqbg2_cs(icol,1,1)
     .       +abs(Q(-j)*q1+R(-j)*r1*prop)**2*qbgZqbg2_cs(icol,2,2)
     .       +abs(Q(-j)*q1+L(-j)*r1*prop)**2*qbgZqbg2_cs(icol,1,2)
     .       +abs(Q(-j)*q1+R(-j)*l1*prop)**2*qbgZqbg2_cs(icol,2,1))
          enddo
             Bmcfm(:,:,:)=fac*(
     .       +abs(Q(-j)*q1+L(-j)*l1*prop)**2*BqbgZqbg2(:,:,:,1,1)
     .       +abs(Q(-j)*q1+R(-j)*r1*prop)**2*BqbgZqbg2(:,:,:,2,2)
     .       +abs(Q(-j)*q1+L(-j)*r1*prop)**2*BqbgZqbg2(:,:,:,1,2)
     .       +abs(Q(-j)*q1+R(-j)*l1*prop)**2*BqbgZqbg2(:,:,:,2,1))
      elseif (gqqg) then
          do icol=0,2
             msq_cs(icol,j,k)=fac*(
     .       +abs(Q(k)*q1+L(k)*l1*prop)**2*gqZqg2_cs(icol,1,1)
     .       +abs(Q(k)*q1+R(k)*r1*prop)**2*gqZqg2_cs(icol,2,2)
     .       +abs(Q(k)*q1+L(k)*r1*prop)**2*gqZqg2_cs(icol,1,2)
     .       +abs(Q(k)*q1+R(k)*l1*prop)**2*gqZqg2_cs(icol,2,1))
          enddo
             Bmcfm(:,:,:)=fac*(
     .       +abs(Q(k)*q1+L(k)*l1*prop)**2*BgqZqg2(:,:,:,1,1)
     .       +abs(Q(k)*q1+R(k)*r1*prop)**2*BgqZqg2(:,:,:,2,2)
     .       +abs(Q(k)*q1+L(k)*r1*prop)**2*BgqZqg2(:,:,:,1,2)
     .       +abs(Q(k)*q1+R(k)*l1*prop)**2*BgqZqg2(:,:,:,2,1))
      elseif (gaag) then
          do icol=0,2
             msq_cs(icol,j,k)=fac*(
     .       +abs(Q(-k)*q1+L(-k)*l1*prop)**2*gqbZqbg2_cs(icol,1,1)
     .       +abs(Q(-k)*q1+R(-k)*r1*prop)**2*gqbZqbg2_cs(icol,2,2)
     .       +abs(Q(-k)*q1+L(-k)*r1*prop)**2*gqbZqbg2_cs(icol,1,2)
     .       +abs(Q(-k)*q1+R(-k)*l1*prop)**2*gqbZqbg2_cs(icol,2,1))
          enddo
             Bmcfm(:,:,:)=fac*(
     .       +abs(Q(-k)*q1+L(-k)*l1*prop)**2*BgqbZqbg2(:,:,:,1,1)
     .       +abs(Q(-k)*q1+R(-k)*r1*prop)**2*BgqbZqbg2(:,:,:,2,2)
     .       +abs(Q(-k)*q1+L(-k)*r1*prop)**2*BgqbZqbg2(:,:,:,1,2)
     .       +abs(Q(-k)*q1+R(-k)*l1*prop)**2*BgqbZqbg2(:,:,:,2,1))
      endif


      msq(j,k)=msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k)

   19 continue
      endif   ! endif Gflag


      faclo=4d0*V*gsq**2*esq**2*aveqq 
      if (Qflag) then
      if (recalc_q) then
         call spinoru(6,p,za,zb)
         prop=s(3,4)/dcmplx((s(3,4)-zmass**2),zmass*zwidth)
c---     qRb->qRb
         call ampqqb_qqb(1,5,2,6,qRb_a,qRb_b)
         call ampqqb_qqb(1,2,5,6,qqb_a,qqb_b)
c---     qR->qR
         call ampqqb_qqb(1,5,6,2,qR_a,qR_b)
         call ampqqb_qqb(1,6,5,2,qq_a,qq_b)
c---     qbR->qbR
         call ampqqb_qqb(6,1,5,2,qbR_a,qbR_b)
         call ampqqb_qqb(2,1,5,6,qbq_a,qbq_b)
c---     qbRb->qbRb
         call ampqqb_qqb(5,1,2,6,qbRb_a,qbRb_b)
          call ampqqb_qqb(6,1,2,5,qbqb_a,qbqb_b)
      endif


      do icol=0,2
      msq_cs(icol,j,k)=zip
      enddo
      




          if (qqqq) then
c----QQ case
            if (j .ne. k) then
            a111=(Q(j)*q1+L(j)*l1*prop)*qR_a(1,1,1)
     .          +(Q(k)*q1+L(k)*l1*prop)*qR_b(1,1,1)
            a121=(Q(j)*q1+L(j)*l1*prop)*qR_a(1,2,1)
     .          +(Q(k)*q1+R(k)*l1*prop)*qR_b(1,2,1)
            a112=(Q(j)*q1+L(j)*r1*prop)*qR_a(1,1,2)
     .          +(Q(k)*q1+L(k)*r1*prop)*qR_b(1,1,2)
            a122=(Q(j)*q1+L(j)*r1*prop)*qR_a(1,2,2)
     .          +(Q(k)*q1+R(k)*r1*prop)*qR_b(1,2,2)
            a211=(Q(j)*q1+R(j)*l1*prop)*qR_a(2,1,1)
     .          +(Q(k)*q1+L(k)*l1*prop)*qR_b(2,1,1)
            a221=(Q(j)*q1+R(j)*l1*prop)*qR_a(2,2,1)
     .          +(Q(k)*q1+R(k)*l1*prop)*qR_b(2,2,1)
            a212=(Q(j)*q1+R(j)*r1*prop)*qR_a(2,1,2)
     .          +(Q(k)*q1+L(k)*r1*prop)*qR_b(2,1,2)
            a222=(Q(j)*q1+R(j)*r1*prop)*qR_a(2,2,2)
     .          +(Q(k)*q1+R(k)*r1*prop)*qR_b(2,2,2)
            msq_cs(0,j,k)=zip
            msq_cs(1,j,k)=
     .      +faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     .             +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
            msq_cs(2,j,k)=zip
            elseif (j .eq. k) then
            a111=(Q(j)*q1+L(j)*l1*prop)*(qR_a(1,1,1)+qR_b(1,1,1))
            b111=(Q(j)*q1+L(j)*l1*prop)*(qq_a(1,1,1)+qq_b(1,1,1))
            a112=(Q(j)*q1+L(j)*r1*prop)*(qR_a(1,1,2)+qR_b(1,1,2))
            b112=(Q(j)*q1+L(j)*r1*prop)*(qq_a(1,1,2)+qq_b(1,1,2))
            a221=(Q(j)*q1+R(j)*l1*prop)*(qR_a(2,2,1)+qR_b(2,2,1))
            b221=(Q(j)*q1+R(j)*l1*prop)*(qq_a(2,2,1)+qq_b(2,2,1))
            a222=(Q(j)*q1+R(j)*r1*prop)*(qR_a(2,2,2)+qR_b(2,2,2))
            b222=(Q(j)*q1+R(j)*r1*prop)*(qq_a(2,2,2)+qq_b(2,2,2))

            a121=(Q(j)*q1+L(j)*l1*prop)*qR_a(1,2,1)
     .          +(Q(k)*q1+R(k)*l1*prop)*qR_b(1,2,1)
            b121=(Q(j)*q1+L(j)*l1*prop)*qq_a(1,2,1)
     .          +(Q(k)*q1+R(k)*l1*prop)*qq_b(1,2,1)
            a122=(Q(j)*q1+L(j)*r1*prop)*qR_a(1,2,2)
     .          +(Q(k)*q1+R(k)*r1*prop)*qR_b(1,2,2)
            b122=(Q(j)*q1+L(j)*r1*prop)*qq_a(1,2,2)
     .          +(Q(k)*q1+R(k)*r1*prop)*qq_b(1,2,2)
            a211=(Q(j)*q1+R(j)*l1*prop)*qR_a(2,1,1)
     .          +(Q(k)*q1+L(k)*l1*prop)*qR_b(2,1,1)
            b211=(Q(j)*q1+R(j)*l1*prop)*qq_a(2,1,1)
     .          +(Q(k)*q1+L(k)*l1*prop)*qq_b(2,1,1)
            a212=(Q(j)*q1+R(j)*r1*prop)*qR_a(2,1,2)
     .          +(Q(k)*q1+L(k)*r1*prop)*qR_b(2,1,2)
            b212=(Q(j)*q1+R(j)*r1*prop)*qq_a(2,1,2)
     .          +(Q(k)*q1+L(k)*r1*prop)*qq_b(2,1,2)

            msq_cs(0,j,k)=half*faclo*(
     .      +Dble(a111*Dconjg(b111))+Dble(a112*Dconjg(b112))
     .      +Dble(a221*Dconjg(b221))+Dble(a222*Dconjg(b222)))*two/xn
            msq_cs(1,j,k)=half*faclo*
     .      (abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     .      +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
            msq_cs(2,j,k)=half*faclo*(
     .      +abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     .      +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
            endif
          elseif (aaaa) then
c----QbQb case
            if (j .ne. k) then
            a111=(Q(-j)*q1+L(-j)*l1*prop)*qbRb_a(1,1,1)
     .          +(Q(-k)*q1+L(-k)*l1*prop)*qbRb_b(1,1,1)
            a121=(Q(-j)*q1+L(-j)*l1*prop)*qbRb_a(1,2,1)
     .          +(Q(-k)*q1+R(-k)*l1*prop)*qbRb_b(1,2,1)

            a112=(Q(-j)*q1+L(-j)*r1*prop)*qbRb_a(1,1,2)
     .          +(Q(-k)*q1+L(-k)*r1*prop)*qbRb_b(1,1,2)
            a122=(Q(-j)*q1+L(-j)*r1*prop)*qbRb_a(1,2,2)
     .          +(Q(-k)*q1+R(-k)*r1*prop)*qbRb_b(1,2,2)

            a211=(Q(-j)*q1+R(-j)*l1*prop)*qbRb_a(2,1,1)
     .          +(Q(-k)*q1+L(-k)*l1*prop)*qbRb_b(2,1,1)
            a221=(Q(-j)*q1+R(-j)*l1*prop)*qbRb_a(2,2,1)
     .          +(Q(-k)*q1+R(-k)*l1*prop)*qbRb_b(2,2,1)

            a212=(Q(-j)*q1+R(-j)*r1*prop)*qbRb_a(2,1,2)
     .          +(Q(-k)*q1+L(-k)*r1*prop)*qbRb_b(2,1,2)
            a222=(Q(-j)*q1+R(-j)*r1*prop)*qbRb_a(2,2,2)
     .          +(Q(-k)*q1+R(-k)*r1*prop)*qbRb_b(2,2,2)
            msq_cs(0,j,k)=zip
            msq_cs(1,j,k)=
     .      +faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     .             +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
            msq_cs(2,j,k)=zip
            elseif (j .eq. k) then

            a111=(Q(-j)*q1+L(-j)*l1*prop)*(qbRb_a(1,1,1)+qbRb_b(1,1,1))
            b111=(Q(-j)*q1+L(-j)*l1*prop)*(qbqb_a(1,1,1)+qbqb_b(1,1,1))
            a112=(Q(-j)*q1+L(-j)*r1*prop)*(qbRb_a(1,1,2)+qbRb_b(1,1,2))
            b112=(Q(-j)*q1+L(-j)*r1*prop)*(qbqb_a(1,1,2)+qbqb_b(1,1,2))
            a221=(Q(-j)*q1+R(-j)*l1*prop)*(qbRb_a(2,2,1)+qbRb_b(2,2,1))
            b221=(Q(-j)*q1+R(-j)*l1*prop)*(qbqb_a(2,2,1)+qbqb_b(2,2,1))
            a222=(Q(-j)*q1+R(-j)*r1*prop)*(qbRb_a(2,2,2)+qbRb_b(2,2,2))
            b222=(Q(-j)*q1+R(-j)*r1*prop)*(qbqb_a(2,2,2)+qbqb_b(2,2,2))


            a121=(Q(-j)*q1+L(-j)*l1*prop)*qbRb_a(1,2,1)
     .          +(Q(-k)*q1+R(-k)*l1*prop)*qbRb_b(1,2,1)
            a122=(Q(-j)*q1+L(-j)*r1*prop)*qbRb_a(1,2,2)
     .          +(Q(-k)*q1+R(-k)*r1*prop)*qbRb_b(1,2,2)
            a211=(Q(-j)*q1+R(-j)*l1*prop)*qbRb_a(2,1,1)
     .          +(Q(-k)*q1+L(-k)*l1*prop)*qbRb_b(2,1,1)
            a212=(Q(-j)*q1+R(-j)*r1*prop)*qbRb_a(2,1,2)
     .          +(Q(-k)*q1+L(-k)*r1*prop)*qbRb_b(2,1,2)

            b121=(Q(-j)*q1+L(-j)*l1*prop)*qbqb_a(1,2,1)
     .          +(Q(-k)*q1+R(-k)*l1*prop)*qbqb_b(1,2,1)
            b122=(Q(-j)*q1+L(-j)*r1*prop)*qbqb_a(1,2,2)
     .          +(Q(-k)*q1+R(-k)*r1*prop)*qbqb_b(1,2,2)
            b211=(Q(-j)*q1+R(-j)*l1*prop)*qbqb_a(2,1,1)
     .          +(Q(-k)*q1+L(-k)*l1*prop)*qbqb_b(2,1,1)
            b212=(Q(-j)*q1+R(-j)*r1*prop)*qbqb_a(2,1,2)
     .          +(Q(-k)*q1+L(-k)*r1*prop)*qbqb_b(2,1,2)


            msq_cs(0,j,k)=half*faclo*(
     .      +Dble(a111*Dconjg(b111))+Dble(a112*Dconjg(b112))
     .      +Dble(a221*Dconjg(b221))+Dble(a222*Dconjg(b222)))*two/xn
            msq_cs(1,j,k)=half*faclo*
     .      (abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     .      +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
            msq_cs(2,j,k)=half*faclo*(
     .      +abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     .      +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
            endif
C---q-qb case
         elseif (qaqa) then
             if (j .ne. -k) then 
            a111=(Q(+j)*q1+L(+j)*l1*prop)*qRb_a(1,1,1)
     .          +(Q(-k)*q1+L(-k)*l1*prop)*qRb_b(1,1,1)
            a112=(Q(+j)*q1+L(+j)*r1*prop)*qRb_a(1,1,2)
     .          +(Q(-k)*q1+L(-k)*r1*prop)*qRb_b(1,1,2)
            a221=(Q(+j)*q1+R(+j)*l1*prop)*qRb_a(2,2,1)
     .          +(Q(-k)*q1+R(-k)*l1*prop)*qRb_b(2,2,1)
            a222=(Q(+j)*q1+R(+j)*r1*prop)*qRb_a(2,2,2)
     .          +(Q(-k)*q1+R(-k)*r1*prop)*qRb_b(2,2,2)

            a121=(Q(+j)*q1+L(+j)*l1*prop)*qRb_a(1,2,1)
     .          +(Q(-k)*q1+R(-k)*l1*prop)*qRb_b(1,2,1)
            a122=(Q(+j)*q1+L(+j)*r1*prop)*qRb_a(1,2,2)
     .          +(Q(-k)*q1+R(-k)*r1*prop)*qRb_b(1,2,2)
            a211=(Q(+j)*q1+R(+j)*l1*prop)*qRb_a(2,1,1)
     .          +(Q(-k)*q1+L(-k)*l1*prop)*qRb_b(2,1,1)
            a212=(Q(+j)*q1+R(+j)*r1*prop)*qRb_a(2,1,2)
     .          +(Q(-k)*q1+L(-k)*r1*prop)*qRb_b(2,1,2)
            msq_cs(0,j,k)=zip
            msq_cs(1,j,k)=zip
            msq_cs(2,j,k)=
     .      +faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     .             +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)

            elseif (j .eq. -k) then
c--case where final state from annihilation diagrams is the same quark
            if (same) then 
            a111=(Q(j)*q1+L(j)*l1*prop)*(qRb_a(1,1,1)+qRb_b(1,1,1))
            b111=(Q(j)*q1+L(j)*l1*prop)*(qqb_a(1,1,1)+qqb_b(1,1,1))

            a112=(Q(j)*q1+L(j)*r1*prop)*(qRb_a(1,1,2)+qRb_b(1,1,2))
            b112=(Q(j)*q1+L(j)*r1*prop)*(qqb_a(1,1,2)+qqb_b(1,1,2))

            a221=(Q(j)*q1+R(j)*l1*prop)*(qRb_a(2,2,1)+qRb_b(2,2,1))
            b221=(Q(j)*q1+R(j)*l1*prop)*(qqb_a(2,2,1)+qqb_b(2,2,1))

            a222=(Q(j)*q1+R(j)*r1*prop)*(qRb_a(2,2,2)+qRb_b(2,2,2))
            b222=(Q(j)*q1+R(j)*r1*prop)*(qqb_a(2,2,2)+qqb_b(2,2,2))

            a121=(Q(+j)*q1+L(+j)*l1*prop)*qRb_a(1,2,1)
     .          +(Q(-k)*q1+R(-k)*l1*prop)*qRb_b(1,2,1)
            a122=(Q(+j)*q1+L(+j)*r1*prop)*qRb_a(1,2,2)
     .          +(Q(-k)*q1+R(-k)*r1*prop)*qRb_b(1,2,2)
            a211=(Q(+j)*q1+R(+j)*l1*prop)*qRb_a(2,1,1)
     .          +(Q(-k)*q1+L(-k)*l1*prop)*qRb_b(2,1,1)
            a212=(Q(+j)*q1+R(+j)*r1*prop)*qRb_a(2,1,2)
     .          +(Q(-k)*q1+L(-k)*r1*prop)*qRb_b(2,1,2)

            b121=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,2,1)
     .          +(Q(-k)*q1+R(-k)*l1*prop)*qqb_b(1,2,1)
            b122=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,2,2)
     .          +(Q(-k)*q1+R(-k)*r1*prop)*qqb_b(1,2,2)
            b211=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,1,1)
     .          +(Q(-k)*q1+L(-k)*l1*prop)*qqb_b(2,1,1)
            b212=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,1,2)
     .          +(Q(-k)*q1+L(-k)*r1*prop)*qqb_b(2,1,2)

            msq_cs(0,j,k)=faclo*(
     .      +Dble(a111*Dconjg(b111))+Dble(a112*Dconjg(b112))
     .      +Dble(a221*Dconjg(b221))+Dble(a222*Dconjg(b222)))*two/xn
            msq_cs(1,j,k)=faclo*(
     .      +abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     .      +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
            msq_cs(2,j,k)=faclo*
     .      (abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     .      +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
         else
!       if ((j.eq.1).or.(j.eq.3).or.(j.eq.5)) then
!           nup=2
!           ndo=nf-3
!       else
!           nup=1
!           ndo=nf-2
!       endif
            b111=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,1,1)
     .          +(Q(+1)*q1+L(+1)*l1*prop)*qqb_b(1,1,1)
            b112=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,1,2)
     .          +(Q(+1)*q1+L(+1)*r1*prop)*qqb_b(1,1,2)
            b221=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,2,1)
     .          +(Q(+1)*q1+R(+1)*l1*prop)*qqb_b(2,2,1)
            b222=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,2,2)
     .          +(Q(+1)*q1+R(+1)*r1*prop)*qqb_b(2,2,2)
            b121=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,2,1)
     .          +(Q(+1)*q1+R(+1)*l1*prop)*qqb_b(1,2,1)
            b122=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,2,2)
     .          +(Q(+1)*q1+R(+1)*r1*prop)*qqb_b(1,2,2)
            b211=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,1,1)
     .          +(Q(+1)*q1+L(+1)*l1*prop)*qqb_b(2,1,1)
            b212=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,1,2)
     .          +(Q(+1)*q1+L(+1)*r1*prop)*qqb_b(2,1,2)
            
      tdo=faclo*(abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     .             +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)

            b111=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,1,1)
     .          +(Q(+2)*q1+L(+2)*l1*prop)*qqb_b(1,1,1)
            b112=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,1,2)
     .          +(Q(+2)*q1+L(+2)*r1*prop)*qqb_b(1,1,2)
            b221=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,2,1)
     .          +(Q(+2)*q1+R(+2)*l1*prop)*qqb_b(2,2,1)
            b222=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,2,2)
     .          +(Q(+2)*q1+R(+2)*r1*prop)*qqb_b(2,2,2)
            b121=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,2,1)
     .          +(Q(+2)*q1+R(+2)*l1*prop)*qqb_b(1,2,1)
            b122=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,2,2)
     .          +(Q(+2)*q1+R(+2)*r1*prop)*qqb_b(1,2,2)
            b211=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,1,1)
     .          +(Q(+2)*q1+L(+2)*l1*prop)*qqb_b(2,1,1)
            b212=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,1,2)
     .          +(Q(+2)*q1+L(+2)*r1*prop)*qqb_b(2,1,2)
            
      tup=faclo*(abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     .          +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)

      msq_cs(1,j,k)=msq_cs(1,j,k)+dfloat(nup)*tup+dfloat(ndo)*tdo
      endif ! if same 
      endif
      elseif (aqqa) then
C---Qb-q case
            if (j .ne. -k) then
            a111=(Q(-j)*q1+L(-j)*l1*prop)*qbR_a(1,1,1)
     .          +(Q(+k)*q1+L(+k)*l1*prop)*qbR_b(1,1,1)
            a121=(Q(-j)*q1+L(-j)*l1*prop)*qbR_a(1,2,1)
     .          +(Q(+k)*q1+R(+k)*l1*prop)*qbR_b(1,2,1)
            a112=(Q(-j)*q1+L(-j)*r1*prop)*qbR_a(1,1,2)
     .          +(Q(+k)*q1+L(+k)*r1*prop)*qbR_b(1,1,2)
            a122=(Q(-j)*q1+L(-j)*r1*prop)*qbR_a(1,2,2)
     .          +(Q(+k)*q1+R(+k)*r1*prop)*qbR_b(1,2,2)
            a211=(Q(-j)*q1+R(-j)*l1*prop)*qbR_a(2,1,1)
     .          +(Q(+k)*q1+L(+k)*l1*prop)*qbR_b(2,1,1)
            a221=(Q(-j)*q1+R(-j)*l1*prop)*qbR_a(2,2,1)
     .          +(Q(+k)*q1+R(+k)*l1*prop)*qbR_b(2,2,1)
            a212=(Q(-j)*q1+R(-j)*r1*prop)*qbR_a(2,1,2)
     .          +(Q(+k)*q1+L(+k)*r1*prop)*qbR_b(2,1,2)
            a222=(Q(-j)*q1+R(-j)*r1*prop)*qbR_a(2,2,2)
     .          +(Q(+k)*q1+R(+k)*r1*prop)*qbR_b(2,2,2)

            msq_cs(0,j,k)=zip
            msq_cs(1,j,k)=zip
            msq_cs(2,j,k)=
     .      +faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     .             +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
            elseif (j .eq. -k) then
            if (same) then 
            a111=(Q(-j)*q1+L(-j)*l1*prop)*(qbR_a(1,1,1)+qbR_b(1,1,1))
            b111=(Q(-j)*q1+L(-j)*l1*prop)*(qbq_a(1,1,1)+qbq_b(1,1,1))
            a112=(Q(-j)*q1+L(-j)*r1*prop)*(qbR_a(1,1,2)+qbR_b(1,1,2))
            b112=(Q(-j)*q1+L(-j)*r1*prop)*(qbq_a(1,1,2)+qbq_b(1,1,2))
            a221=(Q(-j)*q1+R(-j)*l1*prop)*(qbR_a(2,2,1)+qbR_b(2,2,1))
            b221=(Q(-j)*q1+R(-j)*l1*prop)*(qbq_a(2,2,1)+qbq_b(2,2,1))
            a222=(Q(-j)*q1+R(-j)*r1*prop)*(qbR_a(2,2,2)+qbR_b(2,2,2))
            b222=(Q(-j)*q1+R(-j)*r1*prop)*(qbq_a(2,2,2)+qbq_b(2,2,2))

            a121=(Q(-j)*q1+L(-j)*l1*prop)*qbR_a(1,2,1)
     .          +(Q(+k)*q1+R(+k)*l1*prop)*qbR_b(1,2,1)
            a122=(Q(-j)*q1+L(-j)*r1*prop)*qbR_a(1,2,2)
     .          +(Q(+k)*q1+R(+k)*r1*prop)*qbR_b(1,2,2)
            a211=(Q(-j)*q1+R(-j)*l1*prop)*qbR_a(2,1,1)
     .          +(Q(+k)*q1+L(+k)*l1*prop)*qbR_b(2,1,1)
            a212=(Q(-j)*q1+R(-j)*r1*prop)*qbR_a(2,1,2)
     .          +(Q(+k)*q1+L(+k)*r1*prop)*qbR_b(2,1,2)

            b121=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,2,1)
     .          +(Q(+k)*q1+R(+k)*l1*prop)*qbq_b(1,2,1)
            b122=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,2,2)
     .          +(Q(+k)*q1+R(+k)*r1*prop)*qbq_b(1,2,2)
            b211=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,1,1)
     .          +(Q(+k)*q1+L(+k)*l1*prop)*qbq_b(2,1,1)
            b212=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,1,2)
     .          +(Q(+k)*q1+L(+k)*r1*prop)*qbq_b(2,1,2)

            msq_cs(0,j,k)=faclo*(
     .      +Dble(a111*Dconjg(b111))+Dble(a112*Dconjg(b112))
     .      +Dble(a221*Dconjg(b221))+Dble(a222*Dconjg(b222)))*two/xn
            msq_cs(2,j,k)=faclo*
     .      (abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     .      +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
            msq_cs(1,j,k)=faclo*(
     .      +abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     .      +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
         else

c--Here we must also add the contribution of other final state quarks
c  unequal to initial annihilating quarks
!       if ((k.eq.1).or.(k.eq.3).or.(k.eq.5)) then
!           nup=2
!           ndo=nf-3
!       else
!           nup=1
!           ndo=nf-2
!       endif
            b111=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,1,1)
     .          +(Q(+3)*q1+L(+3)*l1*prop)*qbq_b(1,1,1)
            b112=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,1,2)
     .          +(Q(+3)*q1+L(+3)*r1*prop)*qbq_b(1,1,2)
            b221=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,2,1)
     .          +(Q(+3)*q1+R(+3)*l1*prop)*qbq_b(2,2,1)
            b222=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,2,2)
     .          +(Q(+3)*q1+R(+3)*r1*prop)*qbq_b(2,2,2)
            b121=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,2,1)
     .          +(Q(+3)*q1+R(+3)*l1*prop)*qbq_b(1,2,1)
            b122=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,2,2)
     .          +(Q(+3)*q1+R(+3)*r1*prop)*qbq_b(1,2,2)
            b211=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,1,1)
     .          +(Q(+3)*q1+L(+3)*l1*prop)*qbq_b(2,1,1)
            b212=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,1,2)
     .          +(Q(+3)*q1+L(+3)*r1*prop)*qbq_b(2,1,2)
      tdo=faclo*(abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     .          +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)

            b111=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,1,1)
     .          +(Q(+2)*q1+L(+2)*l1*prop)*qbq_b(1,1,1)
            b112=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,1,2)
     .          +(Q(+2)*q1+L(+2)*r1*prop)*qbq_b(1,1,2)
            b221=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,2,1)
     .          +(Q(+2)*q1+R(+2)*l1*prop)*qbq_b(2,2,1)
            b222=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,2,2)
     .          +(Q(+2)*q1+R(+2)*r1*prop)*qbq_b(2,2,2)
            b121=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,2,1)
     .          +(Q(+2)*q1+R(+2)*l1*prop)*qbq_b(1,2,1)
            b122=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,2,2)
     .          +(Q(+2)*q1+R(+2)*r1*prop)*qbq_b(1,2,2)
            b211=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,1,1)
     .          +(Q(+2)*q1+L(+2)*l1*prop)*qbq_b(2,1,1)
            b212=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,1,2)
     .          +(Q(+2)*q1+L(+2)*r1*prop)*qbq_b(2,1,2)
      tup=faclo*(abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     .          +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)

      msq_cs(1,j,k)=msq_cs(1,j,k)+dfloat(nup)*tup+dfloat(ndo)*tdo

      endif                     ! if same 

          endif
          endif
      msq(j,k)=msq(j,k)+msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k)
      endif

C--for safety
      if (msq(j,k) .eq. 0) then
      write(*,*) 'msq=0,msq(j,k),j,k', msq(j,k),j,k
      endif

c--- fill common block qq_cs
      qq_cs(:)=msq_cs(:,j,k)
      res = msq(j,k) 
      return
      end
          
    

      
     
