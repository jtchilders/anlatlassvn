      subroutine qqb_w2jet_pwhg(p,bflav,res)
************************************************************************
*     Author: R.K. Ellis, June 2012                                    * 
*     Routine derived from the MCFM routine qqq_w2jet.f                *
*     to suit the needs of POWHEG                                      *
*    matrix element squared and averaged over initial colours and spins*
*    q(-p1) + aar(-p2) --> W + f(p5) + f(p6)                           *
*                           |                                          *
*                            --> nu(p3) + e^+(p4)                      *
*                     or     --> e^-(p3) + nu~(p4)                     *
*     This routine only work for ckmallowed diagonal                   *
*    all momenta are incoming and in MCFM notation                     *
************************************************************************
      implicit none
      include 'constants.f'
      include 'nwz.f'
      include 'flags.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'ckmallowed.f'
      include 'mmsq_cs.f'
      include 'msq_cs.f'
      include 'Bmcfm.f'
      include 'qq_cs.f'


      integer bflav(6),nq,i,j,i5,i6,nu
      double precision p(mxpart,4),res,facgg,facqq,prop
      double precision, save :: 
     & qaWgg2,aqWgg2,qgWqg2,agWag2,gaWag2,gqWqg2,ggWaq2,
     & qa_ijkk(0:2),qa_ijii(0:2),qa_ijjj(0:2),qa_ijkj(0:2),
     & qa_ijik(0:2),qa_ijkl(0:2),qa_iiij(0:2),qa_iiji(0:2),
     & aq_ijkk(0:2),aq_ijii(0:2),aq_ijjj(0:2),aq_ijkj(0:2),
     & aq_ijik(0:2),aq_ijkl(0:2),aq_iiij(0:2),aq_iiji(0:2),
     & qq_iiji(0:2),qq_ijkj(0:2),qq_ijik(0:2),
     & qq_ijjj(0:2),qq_ijii(0:2),
     & aa_iiji(0:2),aa_ijkj(0:2),aa_ijik(0:2),
     & aa_ijjj(0:2),aa_ijii(0:2),
     & ggBmcfm(4,4,6),qgBmcfm(4,4,6),agBmcfm(4,4,6),gqBmcfm(4,4,6),
     & gaBmcfm(4,4,6),qaBmcfm(4,4,6),aqBmcfm(4,4,6),
     & ggmmsq_cs(0:2),qgmmsq_cs(0:2),agmmsq_cs(0:2),gqmmsq_cs(0:2),
     & gammsq_cs(0:2),qammsq_cs(0:2),aqmmsq_cs(0:2)
      double complex qa1(3),qa2(3),qa3(3),qa4(3),
     &               qq1(4),qq2(4),qq3(4),qq4(4),
     &               aq1(3),aq2(3),aq3(3),aq4(3),
     &               aa1(4),aa2(4),aa3(4),aa4(4)
      logical first,qa,aq,qq,aa,gq,ga,qg,ag,gg,qagg,aqgg,diagonal

      double complex zamub(mxpart,4,mxpart)
      common/zamub/zamub
      data first/.true./
      save first
C     variables needed to avoid recalculating same stuff
      logical ::  recalc_g, recalc_q
      real * 8, save :: opin_g(mxpart,4)
      real * 8, save :: opin_q(mxpart,4)

      if (first) then
      write(6,*) 'nwz',nwz
C----setup allowed values of ckm (diagonal)
            diagonal=.true.
            call setupckmallowed(nwz,diagonal)
      first=.false.
      endif      
 
C----zero out Bmcfm which may get refilled
      Bmcfm(:,:,:)=0d0

C------first determine if it is a two- or four-quark process 
      nq = 0 
      do i=1,6
         if (i < 3 .or. i > 4) then 
            if (abs(bflav(i)) > 0) then 
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


C     decide if recalcaculation is needed 
      recalc_q = .false. 
      recalc_g = .false. 
      

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
      

c--- set up spinors
      call spinoru(6,p,za,zb)
      prop=s(3,4)**2/((s(3,4)-wmass**2)**2+wmass**2*wwidth**2)
      facqq=4d0*V*gsq**2*(gwsq/2d0)**2*aveqq*prop
      facgg=V*xn/four*(gwsq/2d0)**2*gsq**2*prop
      i=bflav(1)
      j=bflav(2)
      i5=bflav(5)
      i6=bflav(6)



c--- calculate 2-quark, 2-gluon amplitudes
      if (Gflag) then
      qagg=(i > 0) .and. (j < 0) .and. ckmallowed(i,j)
      aqgg=(i < 0) .and. (j > 0) .and. ckmallowed(i,j)
      qg=(i > 0) .and. (j == 0) .and. ckmallowed(i,-i5)
      ag=(i < 0) .and. (j == 0) .and. ckmallowed(i,-i5)
      gq=(i == 0) .and. (j > 0) .and. ckmallowed(j,-i5) 
      ga=(i == 0) .and. (j < 0) .and. ckmallowed(j,-i5) 
      gg=(i == 0) .and. (j == 0) .and. ckmallowed(-i5,-i6)  
        
C     fill <a|\gamma^\mu|b] needed for Bmunu calculation
      call spinormu(6,p,zamub)

      if (recalc_g) then
            call w2jetsq(5,6,3,4,1,2,za,zb,ggWaq2)
            call fillBmunu(5,6,3,4,1,2,p,ggBmcfm)
	    ggmmsq_cs(:)=mmsq_cs(:,1,1)

            call   w2jetsq(2,5,3,4,1,6,za,zb,gqWqg2)
            call fillBmunu(2,5,3,4,1,6,p,gqBmcfm)
	    gqmmsq_cs(:)=mmsq_cs(:,1,1)

            call   w2jetsq(5,2,3,4,1,6,za,zb,gaWag2)
            call fillBmunu(5,2,3,4,1,6,p,gaBmcfm)
	    gammsq_cs(:)=mmsq_cs(:,1,1)

            call   w2jetsq(1,5,3,4,2,6,za,zb,qgWqg2)
            call fillBmunu(1,5,3,4,2,6,p,qgBmcfm)
	    qgmmsq_cs(:)=mmsq_cs(:,1,1)

            call   w2jetsq(5,1,3,4,2,6,za,zb,agWag2)
            call fillBmunu(5,1,3,4,2,6,p,agBmcfm)
	    agmmsq_cs(:)=mmsq_cs(:,1,1)

            call   w2jetsq(1,2,3,4,5,6,za,zb,qaWgg2)
            call fillBmunu(1,2,3,4,5,6,p,qaBmcfm)
	    qammsq_cs(:)=mmsq_cs(:,1,1)

            call   w2jetsq(2,1,3,4,5,6,za,zb,aqWgg2)
            call fillBmunu(2,1,3,4,5,6,p,aqBmcfm)
	    aqmmsq_cs(:)=mmsq_cs(:,1,1)
      endif

      if (gg) then
            res=avegg*facgg*ggWaq2
            Bmcfm(:,:,:)=8d0*avegg*facgg*ggBmcfm(:,:,:)
            msq_cs(:,i,j)=avegg*facgg*ggmmsq_cs(:)
      elseif (gq) then
            res=aveqg*facgg*gqWqg2
            Bmcfm(:,:,:)=8d0*aveqg*facgg*gqBmcfm(:,:,:)
            msq_cs(:,i,j)=aveqg*facgg*gqmmsq_cs(:)
      elseif (ga) then
            res=aveqg*facgg*gaWag2
            Bmcfm(:,:,:)=8d0*aveqg*facgg*gaBmcfm(:,:,:)
            msq_cs(:,i,j)=aveqg*facgg*gammsq_cs(:)
      elseif (qg) then
            res=aveqg*facgg*qgWqg2
            Bmcfm(:,:,:)=8d0*aveqg*facgg*qgBmcfm(:,:,:)
            msq_cs(:,i,j)=aveqg*facgg*qgmmsq_cs(:)
      elseif (ag) then
            res=aveqg*facgg*agWag2
            Bmcfm(:,:,:)=8d0*aveqg*facgg*agBmcfm(:,:,:)
            msq_cs(:,i,j)=aveqg*facgg*agmmsq_cs(:)
      elseif (qagg) then
            res=half*aveqq*facgg*qaWgg2
            Bmcfm(:,:,:)=8d0*half*aveqq*facgg*qaBmcfm(:,:,:)
            msq_cs(:,i,j)=half*aveqq*facgg*qammsq_cs(:)
      elseif (aqgg) then
            res=half*aveqq*facgg*aqWgg2
            Bmcfm(:,:,:)=8d0*half*aveqq*facgg*aqBmcfm(:,:,:)
            msq_cs(:,i,j)=half*aveqq*facgg*aqmmsq_cs(:)
      endif

      endif ! end GFLAG

c--- calculate four-quark amplitudes
      if (Qflag) then

      qa=(i>0).and.(j<0).and.(i5>0).and.(i6<0)
      aq=(j>0).and.(i<0).and.(i5>0).and.(i6<0)
      qq=(i>0).and.(j>0).and.(i5>0).and.(i6>0)
      aa=(i<0).and.(j<0).and.(i5<0).and.(i6<0)

c--- basic amplitudes - q a --> W + g* (--> q a) (amps 1 and 3)
c---                and q a --> g* --> q (--> W q) a (amps 2 and 4)           
c--- we label the amplitudes by helicity (qa1 ... qa4)
c--- and by type of contribution qa(1) ... qa(n)

      if (recalc_q) then
      
        call amp_q_QbQ_qb(1,2,5,6,qa1(1),qa2(1),qa3(1),qa4(1))         
c--- crossed - q a --> q a ( --> W a)
        call amp_q_QbQ_qb(1,5,2,6,qa1(2),qa2(2),qa3(2),qa4(2))         
c--- crossed - q a --> q ( --> W q) a 
        call amp_q_QbQ_qb(6,2,5,1,qa1(3),qa2(3),qa3(3),qa4(3))   

c--- now the a q amplitudes
        call amp_q_QbQ_qb(2,1,5,6,aq1(1),aq2(1),aq3(1),aq4(1))         
c--- crossed - a q --> a q ( --> W q)
        call amp_q_QbQ_qb(2,5,1,6,aq1(2),aq2(2),aq3(2),aq4(2))         
c--- crossed - a q --> a ( --> W a) q 
        call amp_q_QbQ_qb(6,1,5,2,aq1(3),aq2(3),aq3(3),aq4(3))   

c--- crossed q q --> q ( --> W q) q
        call amp_q_QbQ_qb(1,5,6,2,qq1(1),qq2(1),qq3(1),qq4(1))   
c--- crossed q q --> q q ( --> W q)
        call amp_q_QbQ_qb(2,5,6,1,qq1(2),qq2(2),qq3(2),qq4(2))   
c--- crossed q q --> q q ( --> W q)
        call amp_q_QbQ_qb(2,6,5,1,qq1(3),qq2(3),qq3(3),qq4(3))   
c--- crossed q q --> q ( --> W q) q
        call amp_q_QbQ_qb(1,6,5,2,qq1(4),qq2(4),qq3(4),qq4(4))   

c--- crossed a a --> a ( --> W a) a
        call amp_q_QbQ_qb(5,1,2,6,aa1(1),aa2(1),aa3(1),aa4(1))   
c--- crossed a a --> a a ( --> W a)
        call amp_q_QbQ_qb(5,2,1,6,aa1(2),aa2(2),aa3(2),aa4(2))   
c--- crossed a a --> a ( --> W a) a
        call amp_q_QbQ_qb(6,2,1,5,aa1(3),aa2(3),aa3(3),aa4(3))   
c--- crossed a a --> a a ( --> W a)
        call amp_q_QbQ_qb(6,1,2,5,aa1(4),aa2(4),aa3(4),aa4(4))   

c--- now square these amplitudes separating into color structures
c   1) Amplitude
c   2) Amplitude with (5<-->6)
c   0) Interference between above
c
        
c--- q(i) a(i) --> g* --> q(j) (--> W q(i)) a(j)
      qa_iiij(1)=abs(qa2(1))**2+abs(qa4(1))**2            
      qa_iiij(2)=abs(qa1(3))**2+abs(qa3(3))**2 
      qa_iiij(0)=2d0/xn*dble(qa2(1)*Dconjg(qa1(3)))
c--- q(i) a(i) --> g* --> q(j) a(j) (--> W a(i))
      qa_iiji(1)=abs(qa2(1))**2+abs(qa4(1))**2            
      qa_iiji(2)=abs(qa1(2))**2+abs(qa3(2))**2 
      qa_iiji(0)=2d0/xn*dble(qa2(1)*Dconjg(qa1(2)))
c--- q(i) a(j) --> g* --> q(l) (--> W q(k)) a(l)
      qa_ijkl(1)=abs(qa2(1))**2+abs(qa4(1))**2
      qa_ijkl(2)=zip
      qa_ijkl(0)=zip
c--- q(i) a(j) --> W + g* (--> q(i) a(i)) i.e. k = i
      qa_ijii(1)=abs(qa1(1))**2+abs(qa3(1))**2
      qa_ijii(2)=abs(qa2(2))**2+abs(qa4(2))**2
      qa_ijii(0)=+2d0/xn*dble(qa1(1)*Dconjg(qa2(2)))
c--- q(i) a(j) --> W + g* (--> q(j) a(j)) i.e. k = j
      qa_ijjj(1)=abs(qa1(1))**2+abs(qa3(1))**2
      qa_ijjj(2)=abs(qa2(3))**2+abs(qa4(3))**2
      qa_ijjj(0)=2d0/xn*dble(qa1(1)*Dconjg(qa2(3)))
c--- q (i) a(j) --> q(i) a(j) ( --> W a(k)) with k != i,j
      qa_ijik(1)=zip
      qa_ijik(2)=abs(qa2(2))**2+abs(qa4(2))**2 
      qa_ijik(0)=zip
c--- q (i) a(j) --> q(i) ( --> W q(k)) a(j) with k != i,j
      qa_ijkj(1)=zip
      qa_ijkj(2)=abs(qa2(3))**2+abs(qa4(3))**2 
      qa_ijkj(0)=zip
c--- q(i) a(j) --> W + g* (--> q(k) a(k)) with k != i,j
      qa_ijkk(1)=abs(qa1(1))**2+abs(qa3(1))**2 
      qa_ijkk(2)=zip
      qa_ijkk(0)=zip

c--- q(i) a(i) --> g* --> q(j) (--> W q(i)) a(j)
      aq_iiij(1)=abs(aq2(1))**2+abs(aq4(1))**2            
      aq_iiij(2)=abs(aq1(3))**2+abs(aq3(3))**2 
      aq_iiij(0)=2d0/xn*dble(aq2(1)*Dconjg(aq1(3)))  
c--- q(i) a(i) --> g* --> q(j) a(j) (--> W a(i))
      aq_iiji(1)=abs(aq2(1))**2+abs(aq4(1))**2            
      aq_iiji(2)=abs(aq1(2))**2+abs(aq3(2))**2 
      aq_iiji(0)=2d0/xn*dble(aq2(1)*Dconjg(aq1(2)))  
c--- q(i) a(j) --> g* --> q(l) (--> W q(k)) a(l)
      aq_ijkl(1)=abs(aq2(1))**2+abs(aq4(1))**2
      aq_ijkl(2)=zip
      aq_ijkl(0)=zip
c--- q(i) a(j) --> W + g* (--> q(k) a(k)) with k != i,j
      aq_ijkk(1)=abs(aq1(1))**2+abs(aq3(1))**2 
      aq_ijkk(2)=zip
      aq_ijkk(0)=zip
c--- q(i) a(j) --> W + g* (--> q(i) a(i)) i.e. k = i
      aq_ijii(1)=abs(aq1(1))**2+abs(aq3(1))**2
      aq_ijii(2)=abs(aq2(2))**2+abs(aq4(2))**2
      aq_ijii(0)=+2d0/xn*dble(aq1(1)*Dconjg(aq2(2)))
c--- q(i) a(j) --> W + g* (--> q(j) a(j)) i.e. k = j
      aq_ijjj(1)=abs(aq1(1))**2+abs(aq3(1))**2
      aq_ijjj(2)=abs(aq2(3))**2+abs(aq4(3))**2
      aq_ijjj(0)=2d0/xn*dble(aq1(1)*Dconjg(aq2(3)))
c--- q (i) a(j) --> q(i) a(j) ( --> W a(k)) with k != i,j
      aq_ijik(2)=abs(aq2(2))**2+abs(aq4(2))**2 
      aq_ijik(1)=zip
      aq_ijik(0)=zip
c--- q (i) a(j) --> q(i) ( --> W q(k)) a(j) with k != i,j
      aq_ijkj(2)=abs(aq2(3))**2+abs(aq4(3))**2 
      aq_ijkj(1)=zip
      aq_ijkj(0)=zip

c--- q(i) q(i) --> q(i) ( --> W q(j) ) q(i)
      qq_iiji(1)=abs(qq1(1))**2+abs(qq3(1))**2
      qq_iiji(2)=abs(qq1(2))**2+abs(qq3(2))**2
      qq_iiji(0)=2d0/xn*dble(qq1(1)*Dconjg(qq1(2)))
c--- q(i) q(j) --> q(i) ( --> W q(k) ) q(j)
      qq_ijkj(1)=abs(qq1(1))**2+abs(qq3(1))**2
      qq_ijkj(2)=zip
      qq_ijkj(0)=zip
c--- q(i) q(j) --> q(i) q(j) ( --> W q(k) )
      qq_ijik(1)=abs(qq1(3))**2+abs(qq3(3))**2
      qq_ijik(2)=zip
      qq_ijik(0)=zip
c--- q(i) q(j) --> q(i) ( --> W q(j) ) q(j)
      qq_ijjj(1)=abs(qq1(1))**2+abs(qq3(1))**2
      qq_ijjj(2)=abs(qq1(4))**2+abs(qq3(4))**2
      qq_ijjj(0)=2d0/xn*dble(qq1(1)*Dconjg(qq1(4)))  
c--- q(i) q(j) --> q(i) q(j) ( --> W q(i) ) 
      qq_ijii(1)=abs(qq1(3))**2+abs(qq3(3))**2
      qq_ijii(2)=abs(qq1(2))**2+abs(qq3(2))**2
      qq_ijii(0)=2d0/xn*dble(qq1(3)*Dconjg(qq1(2)))  

c--- a(i) a(i) --> a(i) ( --> W a(j) ) a(i)
      aa_iiji(1)=abs(aa1(1))**2+abs(aa3(1))**2
      aa_iiji(2)=abs(aa1(2))**2+abs(aa3(2))**2
      aa_iiji(0)=2d0/xn*dble(aa1(1)*Dconjg(aa1(2)))
c--- a(i) a(j) --> a(i) ( --> W a(k) ) a(j)
      aa_ijkj(1)=abs(aa1(1))**2+abs(aa3(1))**2
      aa_ijkj(2)=zip
      aa_ijkj(0)=zip
c--- a(i) a(j) --> a(i) a(j) ( --> W a(k) )
      aa_ijik(1)=abs(aa1(3))**2+abs(aa3(3))**2
      aa_ijik(2)=zip
      aa_ijik(0)=zip
c--- a(i) a(j) --> a(i) ( --> W a(j) ) a(j)
      aa_ijjj(1)=abs(aa1(1))**2+abs(aa3(1))**2
      aa_ijjj(2)=abs(aa1(4))**2+abs(aa3(4))**2
      aa_ijjj(0)=2d0/xn*dble(aa1(1)*Dconjg(aa1(4)))  
c--- a(i) a(j) --> a(i) a(j) ( --> W a(i) ) 
      aa_ijii(2)=abs(aa1(2))**2+abs(aa3(2))**2
      aa_ijii(1)=abs(aa1(3))**2+abs(aa3(3))**2
      aa_ijii(0)=2d0/xn*dble(aa1(2)*Dconjg(aa1(3)))  
      
      endif	! end of recalc_q


      
      if (qa) then
      if (i==-j) then
      if (ckmallowed(j,-i6) .and.(i5==i)) then
c--- q(i) a(i) --> g* --> q(j) (--> W q(i)) a(j)
      msq_cs(0:2,i,j)=facqq*qa_iiij(0:2)
      elseif (ckmallowed(i,-i5).and.(i6==j)) then
c--- q(i) a(i) --> g* --> q(j) a(j) (--> W a(i))
      msq_cs(0:2,i,j)=facqq*qa_iiji(0:2)
      elseif (ckmallowed(-i5,-i6) .and.(i5/=i).and.(i6/=j)) then
c--- q(i) a(j) --> g* --> q(l) (--> W q(k)) a(l)
      msq_cs(0:2,i,j)=facqq*qa_ijkl(0:2)
      endif
      else ! (i/=-j)
      if (ckmallowed(j,-i6) .and.(i5==i).and.(i6==-i5)) then
c--- q(i) a(j) --> W + g* (--> q(i) a(i)) i.e. k = i
      msq_cs(0:2,i,j)=facqq*qa_ijii(0:2)
      elseif (ckmallowed(i,-i5) .and.(i5==-j).and.(i6==-i5)) then
c--- q(i) a(j) --> W + g* (--> q(j) a(j)) i.e. k = j
      msq_cs(0:2,i,j)=facqq*qa_ijjj(0:2)
      elseif (ckmallowed(j,-i6) .and.(i5==i).and.(i6/=-i)) then
c--- q (i) a(j) --> q(i) a(j) ( --> W a(k)) with k != i,j
      msq_cs(0:2,i,j)=facqq*qa_ijik(0:2)
      elseif (ckmallowed(i,-i5) .and.(i5/=j).and.(i6==j)) then
c--- q (i) a(j) --> q(i) ( --> W q(k)) a(j) with k != i,j
      msq_cs(0:2,i,j)=facqq*qa_ijkj(0:2)
      elseif (ckmallowed(i,j) .and.(i5/=i).and.(i6/=j)) then
c--- q(i) a(j) --> W + g* (--> q(k) a(k)) with k != i,j
      msq_cs(0:2,i,j)=facqq*qa_ijkk(0:2)
      endif
      endif
      endif


      if (aq) then
      if (j==-i) then
      if (ckmallowed(i,-i6) .and.(i5==j)) then
c--- q(i) a(i) --> g* --> q(j) (--> W q(i)) a(j)
      msq_cs(0:2,i,j)=facqq*aq_iiij(0:2)
      elseif (ckmallowed(j,-i5) .and.(i6==i)) then
c--- q(i) a(i) --> g* --> q(j) a(j) (--> W a(i))
      msq_cs(0:2,i,j)=facqq*aq_iiji(0:2)
      elseif (ckmallowed(-i5,-i6) .and.(i5/=j).and.(i6/=i)) then
c--- q(i) a(j) --> g* --> q(l) (--> W q(k)) a(l)
      msq_cs(0:2,i,j)=facqq*aq_ijkl(0:2)
      endif
      else ! (j/=-i)
      if (ckmallowed(i,j) .and.(i5/=j).and.(i6/=i)) then
c--- q(i) a(j) --> W + g* (--> q(k) a(k)) with k != i,j
      msq_cs(0:2,i,j)=facqq*aq_ijkk(0:2)
      elseif (ckmallowed(i,-i6) .and.(i5==j).and.(i6==-i5)) then
c--- q(i) a(j) --> W + g* (--> q(i) a(i)) i.e. k = i
      msq_cs(0:2,i,j)=facqq*aq_ijii(0:2)
      elseif (ckmallowed(j,-i5) .and.(i5==-i).and.(i6==-i5)) then
c--- q(i) a(j) --> W + g* (--> q(j) a(j)) i.e. k = j
      msq_cs(0:2,i,j)=facqq*aq_ijjj(0:2)
      elseif (ckmallowed(i,-i6).and.(i5==j).and.(i6/=-j)) then
c--- q (i) a(j) --> q(i) a(j) ( --> W a(k)) with k != i,j
      msq_cs(0:2,i,j)=facqq*aq_ijik(0:2)
      elseif (ckmallowed(j,-i5) .and.(i5/=i).and.(i6==i)) then
c--- q (i) a(j) --> q(i) ( --> W q(k)) a(j) with k != i,j
      msq_cs(0:2,i,j)=facqq*aq_ijkj(0:2)
      endif
      endif
      endif      

      if (qq) then
      if (i==j) then
      if (ckmallowed(i,-i5).and.(i6==j)) then
c--- q(i) q(i) --> q(i) ( --> W q(j) ) q(i)
      msq_cs(0:2,i,j)=facqq*qq_iiji(0:2)
      endif
      else ! (i /= j)
      if (ckmallowed(i,-i5) .and.(i6==j) .and. (i5 /= j)) then
c--- q(i) q(j) --> q(i) ( --> W q(k) ) q(j)
      msq_cs(0:2,i,j)=facqq*qq_ijkj(0:2)
      elseif (ckmallowed(j,-i6) .and.(i5==i).and.(i6/=i)) then
c--- q(i) q(j) --> q(i) q(j) ( --> W q(k) )
      msq_cs(0:2,i,j)=facqq*qq_ijik(0:2)
      elseif (ckmallowed(i,-i5) .and.(i6==j).and.(i5==i6)) then
c--- q(i) q(j) --> q(i) ( --> W q(j) ) q(j)
      msq_cs(0:2,i,j)=half*facqq*qq_ijjj(0:2)
      elseif (ckmallowed(j,-i6) .and.(i5==i).and.(i5==i6)) then
c--- q(i) q(j) --> q(i) q(j) ( --> W q(i) ) 
      msq_cs(0:2,i,j)=half*facqq*qq_ijii(0:2)
      endif
      endif
      endif

      if (aa) then
      if (i==j) then
      if (ckmallowed(i,-i5) .and.(i6==j)) then
c--- a(i) a(i) --> a(i) ( --> W a(j) ) a(i)
      msq_cs(0:2,i,j)=facqq*aa_iiji(0:2)
      endif 
      else ! (i /= j)
      if (ckmallowed(i,-i5) .and.(i6==j).and.(i5/=j)) then
c--- a(i) a(j) --> a(i) ( --> W a(k) ) a(j)
      msq_cs(0:2,i,j)=facqq*aa_ijkj(0:2)
      elseif (ckmallowed(j,-i6) .and.(i5==i).and.(i6/=i)) then
c--- a(i) a(j) --> a(i) a(j) ( --> W a(k) )
      msq_cs(0:2,i,j)=facqq*aa_ijik(0:2)
      elseif (ckmallowed(i,-i5) .and.(i5==j).and.(i5==i6)) then
c--- a(i) a(j) --> a(i) ( --> W a(j) ) a(j)
      msq_cs(0:2,i,j)=half*facqq*aa_ijjj(0:2)
      elseif (ckmallowed(j,-i6) .and.(i5==i).and.(i5==i6)) then
c--- a(i) a(j) --> a(i) a(j) ( --> W a(i) ) 
      msq_cs(0:2,i,j)=half*facqq*aa_ijii(0:2)
      endif
      endif
      endif
      
      res=msq_cs(0,i,j)+msq_cs(1,i,j)+msq_cs(2,i,j)

      endif ! end of if(QFLAG)

c--- fill common block qq_cs
      qq_cs(:)=msq_cs(:,i,j)
c      write(6,*) 'i,j',i,j
c      write(6,*) 'qqb_w2jet_pwhg: qq_cs',qq_cs
      
      return
      end
     
      
