      subroutine qqb_w2jet_v_pwhg(p,vflav,res)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     July 2012.                                                       *
*     Routine derived from the MCFM routine qqb_w2jet_v.f              *
*     to suit the needs of POWHEG                                      *
*                                                                      *
*     Calculate the virtual matrix element squared and                 *
*     subtraction terms for the process                                *
*                                                                      *
*     q(-p1) + qbar(-p2) --> W + j(p5) + j(p6)                         *
*                            |                                         *
*                            --> nu(p3) + e^+(p4)                      *
*                       or   --> e^-(p3) + nu~(p4)                     *
*                                                                      *
*     where the partons are either q(p5) and qbar(p6) [Qflag = .true.] *
*                               or g(p5) and g(p6)    [Gflag = .true.] *
*                                                                      *
*     correction is applied to put alpha_s into the MS bar scheme      *
*     but otherwise the result is in the four-dimensional scheme       *
*                                                                      *
************************************************************************

      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'flags.f'
      include 'nwz.f'
      include 'ckmallowed.f'
      include 'scale.f' 
      integer nu,i,j,i5,i6,vflav(6),nq
      double precision fac,res,subuv,
     & p(mxpart,4),q(mxpart,4),pswap(mxpart,4),res0

      double precision, save :: 
     & mmsq_qqb,mmsq_qbq,mmsq_gq,mmsq_gqb,mmsq_qg,mmsq_qbg,mmsq_gg,
     & qqb_ijkk,qqb_iikl,qqb_ijkj,qqb_ijik,
     & qqb_ijii,qqb_ijjj,qqb_iiij,qqb_iiji,
     & qbq_ijkk,qbq_iikl,qbq_ijkj,qbq_ijik,
     & qbq_ijii,qbq_ijjj,qbq_iiij,qbq_iiji,
     & qq_ijkk,qq_iikl,qq_ijkj,qq_ijik,
     & qq_ijii,qq_ijjj,qq_iiij,qq_iiji,
     & qbqb_ijkk,qbqb_iikl,qbqb_ijkj,qbqb_ijik,
     & qbqb_ijii,qbqb_ijjj,qbqb_iiij,qbqb_iiji
      double complex prop
      logical first,qqb,qbq,qq,qbqb,gq,gqb,qg,qbg,gg,qqbgg,qbqgg,
     & diagonal


      data first/.true./
      save first
C     variables needed to avoid recalculating same stuff 
      logical ::  recalc_g, recalc_q 
      real * 8, save :: opin_g(mxpart,4),scale_g 
      real * 8, save :: opin_q(mxpart,4),scale_q 



      if (first) then
      diagonal=.true.
      call setupckmallowed(nwz,diagonal)
      first=.false.
      endif

      scheme = 'msbr' 
      epinv=0d0

C------first determine if it is a two- or four-quark process 
      nq = 0 
      do j=1,6
         if (j < 3 .or. j > 4) then 
            if (abs(vflav(j)) > 0) then 
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
         write(*,*) 'vflav',vflav
         stop 'nq out of range' 
      endif

C     decide if recalculation is needed 
      recalc_q = .false. 
      recalc_g = .false. 
      

      if (Gflag) then
         if (scale_g .ne. scale) then
            recalc_g = .true.
         else
            do i=1,6
               do nu=1,4
                  if(opin_g(i,nu).ne.p(i,nu)) then
                     recalc_g = .true.
                     goto 10
                  endif
               enddo
            enddo
         endif
      else 
         if (scale_q .ne. scale) then
            recalc_q = .true.
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
      endif

 10    continue 
       
C----Store momenta and scale
       if (recalc_g) then
          opin_g = p
          scale_g = scale
       endif
       if (recalc_q) then
          opin_q = p
          scale_q = scale
       endif  
      

      res=0d0

      i=vflav(1)
      j=vflav(2)
      i5=vflav(5)
      i6=vflav(6)

c--- calculate the lowest order matrix element and fill the
c--- common block twopii with s_{ij}
      call qqb_w2jet_pwhg(p,vflav,res0)
      prop=s(3,4)/dcmplx(s(3,4)-wmass**2,wmass*wwidth)

      fac=V*xn*gw**4*gsq**2*ason2pi
      
************************************************************************
*     Contributions from QQGG matrix elements                          *
************************************************************************            
      if (Gflag) then

c----UV counterterm contains the finite renormalization to arrive
c----at MS bar scheme for alpha_s 
      subuv=2d0*xn*(epinv*(11d0-2d0*dble(nf)/xn)-1d0)/6d0

      qqbgg=(i > 0) .and. (j < 0) .and. ckmallowed(i,j)
      qbqgg=(i < 0) .and. (j > 0) .and. ckmallowed(i,j)
      qg=(i > 0) .and. (j == 0) .and. ckmallowed(i,-i5)
      qbg=(i < 0) .and. (j == 0) .and. ckmallowed(i,-i5)
      gq=(i == 0) .and. (j > 0) .and. ckmallowed(j,-i5) 
      gqb=(i == 0) .and. (j < 0) .and. ckmallowed(j,-i5) 
      gg=(i == 0) .and. (j == 0) .and. ckmallowed(-i5,-i6)  

      if (recalc_g) then 

c---  calculate the qqb terms
CALL    0--> q(p2)+g(p5)+g(p6)+qb(p1)+l(p3)+lbar(p4) 
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
      do nu=1,4
      pswap(1,nu)=p(2,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(5,nu)
      pswap(4,nu)=p(1,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_qqb)

c---  calculate the qbq terms
CALL    0--> q(p1)+g(p5)+g(p6)+qb(p2)+l(p3)+lbar(p4) 
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l+(p6)
      do nu=1,4
      pswap(1,nu)=p(1,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(5,nu)
      pswap(4,nu)=p(2,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_qbq)

c---  calculate the gq terms
CALL    0--> q(p5)+g(p1)+g(p6)+qb(p2)+l(p3)+lbar(p4) 
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l+(p6)
      do nu=1,4
      pswap(1,nu)=p(5,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(1,nu)
      pswap(4,nu)=p(2,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_gq)

c---  calculate the qg terms
CALL    0--> q(p5)+g(p2)+g(p6)+qb(p1)+l(p3)+lbar(p4) 
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l+(p6)
      do nu=1,4
      pswap(1,nu)=p(5,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(1,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_qg)

c---  calculate the gqb terms
CALL    0--> q(p2)+g(p1)+g(p6)+qb(p5)+l(p3)+lbar(p4) 
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l+(p6)
      do nu=1,4
      pswap(1,nu)=p(2,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(1,nu)
      pswap(4,nu)=p(5,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_gqb)

c---  calculate the qbg terms
CALL    0--> q(p1)+g(p2)+g(p6)+qb(p5)+l(p3)+lbar(p4) 
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l+(p6)
      do nu=1,4
      pswap(1,nu)=p(1,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(5,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_qbg)

c--- calculate the gg terms
      do nu=1,4
      pswap(1,nu)=p(5,nu)
      pswap(2,nu)=p(1,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(6,nu)
      pswap(5,nu)=p(3,nu)
      pswap(6,nu)=p(4,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_gg)
      endif ! recalc 
      
************************************************************************
*     Contributions from QQGG matrix elements                          *
************************************************************************            


      if     (qqbgg) then
        res=mmsq_qqb*cdabs(prop)**2
     &           *half*(aveqq/avegg)*(gwsq**2/4d0/esq**2)
      elseif (qbqgg) then
        res=mmsq_qbq*cdabs(prop)**2
     &           *half*(aveqq/avegg)*(gwsq**2/4d0/esq**2)
      elseif (qg) then
        res=mmsq_qg*cdabs(prop)**2
     &           *(aveqg/avegg)*(gwsq**2/4d0/esq**2)
      elseif (qbg) then
        res=mmsq_qbg*cdabs(prop)**2
     &           *(aveqg/avegg)*(gwsq**2/4d0/esq**2)
      elseif (gq) then
        res=mmsq_gq*cdabs(prop)**2
     &           *(aveqg/avegg)*(gwsq**2/4d0/esq**2)
      elseif (gqb) then
        res=mmsq_gqb*cdabs(prop)**2
     &           *(aveqg/avegg)*(gwsq**2/4d0/esq**2)
      elseif (gg) then
        res=mmsq_gg*cdabs(prop)**2*(gwsq**2/4d0/esq**2)
      endif

      endif ! end GFLAG
      
************************************************************************
*     Contributions from QQQQ matrix elements                          *
************************************************************************            
      if (Qflag) then
c--- UV counter-term is already included in a6routine.f
      subuv=0d0
      
      qqb=(i>0).and.(j<0).and.(i5>0).and.(i6<0)
      qbq=(j>0).and.(i<0).and.(i5>0).and.(i6<0)
      qq=(i>0).and.(j>0).and.(i5>0).and.(i6>0)
      qbqb=(i<0).and.(j<0).and.(i5<0).and.(i6<0)

      if (recalc_q) then 

c---  Now transform momenta into a notation 
c---  suitable for calling the BDKW function with notation which is 
c---  q-(-p4)+Q+(-p2)+l-(-p5) ---> q+(p1)+Q-(p3)+l+(p6)
      do nu=1,4
      q(1,nu)=p(2,nu)
      q(2,nu)=p(6,nu)
      q(3,nu)=p(5,nu)
      q(4,nu)=p(1,nu)
      q(5,nu)=p(4,nu)
      q(6,nu)=p(3,nu)
      enddo      

      call spinoru(6,q,za,zb)

c--- set-up qqb matrix elements

!      if (qqb) then
      call qqbw2j_loop(1,2,3,4,5,6,qqb_ijkk,qqb_iikl,qqb_ijkj,qqb_ijik,
     &                             qqb_ijii,qqb_ijjj,qqb_iiij,qqb_iiji)
!      elseif (qbq) then
c--- qbq
      call qqbw2j_loop(4,2,3,1,5,6,qbq_ijkk,qbq_iikl,qbq_ijkj,qbq_ijik,
     &                             qbq_ijii,qbq_ijjj,qbq_iiij,qbq_iiji)
!      elseif (qq) then
c--- qq (note that roles of iiij and ijii are reversed)
      call qqbw2j_loop(2,1,3,4,5,6,qq_ijkk,qq_iikl,qq_ijkj,qq_ijik,
     &                             qq_iiij,qq_ijjj,qq_ijii,qq_iiji)
!      elseif (qbqb) then
c--- qbqb (note that roles of ijjj and iiji are reversed)
      call qqbw2j_loop(1,2,4,3,5,6,qbqb_ijkk,qbqb_iikl,qbqb_ijkj,
     &         qbqb_ijik,qbqb_ijii,qbqb_iiji,qbqb_iiij,qbqb_ijjj)


!      endif
      endif ! recalc_q


      if (qqb) then
      if (i==-j) then
      if (ckmallowed(j,-i6).and.(i5==i)) then
c--- q(i) qb(i) --> g* --> q(j) (--> W q(i)) qb(j) iiij
      res=fac*aveqq*qqb_iiij
      elseif (ckmallowed(i,-i5) .and.(i6==j)) then
c--- q(i) qb(i) --> g* --> q(j) qb(j) (--> W qb(i)) iiji
      res=fac*aveqq*qqb_iiji
      elseif (ckmallowed(-i5,-i6).and.(i5/=i).and.(i6/=j)) then
c--- q(i) qb(j) --> g* --> q(l) (--> W q(k)) qb(l) ijkl
      res=fac*aveqq*qqb_iikl
      endif
      else
      if (ckmallowed(j,-i6) .and.(i5==i).and.(i6==-i5)) then
c--- q(i) qb(j) --> W + g* (--> q(i) qb(i)) i.e. k = i ijii
      res=fac*aveqq*qqb_ijii
      elseif (ckmallowed(i,-i5) .and.(i5==-j).and.(i6==-i5)) then
c--- q(i) qb(j) --> W + g* (--> q(j) qb(j)) i.e. k = j ijjj
      res=fac*aveqq*qqb_ijjj
      elseif (ckmallowed(j,-i6) .and.(i5==i).and.(i6/=-i)) then
c--- q (i) qb(j) --> q(i) qb(j) ( --> W qb(k)) with k != i,j ijik
      res=fac*aveqq*qqb_ijik
      elseif (ckmallowed(i,-i5) .and.(i5/=j).and.(j==i6)) then
c--- q (i) qb(j) --> q(i) ( --> W q(k)) qb(j) with k != i,j ijkj
      res=fac*aveqq*qqb_ijkj
      elseif (ckmallowed(i,j) .and.(i5/=i).and.(i6/=j)) then
c--- q(i) qb(j) --> W + g* (--> q(k) qb(k)) with k != i,j ijkk
      res=fac*aveqq*qqb_ijkk
      endif
      endif
      endif


      if (qbq) then
      if (i==-j) then
      if (ckmallowed(i,-i6) .and.(i5==j)) then
c--- qb(i) q(i) --> g* --> q(j) (--> W q(i)) qb(j) iiij
      res=fac*aveqq*qbq_iiij
      elseif (ckmallowed(j,-i5) .and.(i5/=-i) .and.(i==i6)) then
c--- qb(i) q(i) --> g* --> q(j) qb(j) (--> W qb(i)) iiji
      res=fac*aveqq*qbq_iiji
      elseif (ckmallowed(-i5,-i6) .and.(i5/=j).and.(i6/=i)) then
c--- qb(j) q(i) --> g* --> q(l) (--> W q(k)) qb(l) ijkl
      res=fac*aveqq*qbq_iikl
      endif
      else
      if (ckmallowed(i,j) .and.(i5/=j).and.(i6/=i)) then
c--- qb(j) q(i) --> W + g* (--> q(k) qb(k)) with k != i,j ijkk
      res=fac*aveqq*qbq_ijkk
      elseif (ckmallowed(i,-i6) .and.(i5==j).and.(i6==-i5)) then
c--- qb(j) q(i) --> W + g* (--> q(i) qb(i)) i.e. k = i ijii
      res=fac*aveqq*qbq_ijii
      elseif (ckmallowed(j,-i5) .and.(i5==-i).and.(i6==-i5)) then
c--- qb(j) q(i) --> W + g* (--> q(j) qb(j)) i.e. k = j ijjj
      res=fac*aveqq*qbq_ijjj
      elseif (ckmallowed(i,-i6) .and.(i5==j).and.(i6/=-j)) then
c--- qb (j) q(i) --> q(i) qb(j) ( --> W qb(k)) with k != i,j ijik
      res=fac*aveqq*qbq_ijik
      elseif (ckmallowed(j,-i5) .and.(i5/=i).and.(i==i6)) then
c--- qb(j) q(i) --> q(i) ( --> W q(k)) qb(j) with k != i,j ijkj
      res=fac*aveqq*qbq_ijkj
      endif
      endif
      endif
      
  
      if (qq) then
      if (i==j) then
      if (ckmallowed(i,-i5) .and.(j==i6)) then
c--- q(i) q(i) --> q(i) ( --> W q(j) ) q(i) iiji
      res=fac*aveqq*qq_iiji
      endif
      else
      if (ckmallowed(i,-i5) .and.(j==i6).and.(i5/=j)) then
c--- q(i) q(j) --> q(i) ( --> W q(k) ) q(j) ijkj
      res=fac*aveqq*qq_ijkj
      elseif (ckmallowed(j,-i6) .and.(i==i5).and.(i6/=i)) then
c--- q(i) q(j) --> q(i) q(j) ( --> W q(k) ) ijik
      res=fac*aveqq*qq_ijik
      elseif (ckmallowed(i,-i5) .and.(j==i6).and.(i5==i6)) then
c--- q(i) q(j) --> q(i) ( --> W q(j) ) q(j) ijjj
      res=fac*aveqq*half*qq_ijjj
      elseif (ckmallowed(j,-i6) .and.(i==i5).and.(i5==i6)) then
c--- q(i) q(j) --> q(i) q(j) ( --> W q(i) )  ijii
      res=fac*aveqq*half*qq_ijii
      endif
      endif
      endif


      if (qbqb) then
      if (i==j) then
      if (ckmallowed(i,-i5) .and.(i6==j)) then
c--- qb(i) qb(i) --> qb(i) ( --> W qb(j) ) qb(i) iiji
      res=fac*aveqq*qbqb_iiji
      endif
      else
      if (ckmallowed(i,-i5) .and.(j==i6).and.(i5/=j)) then
c--- qb(i) qb(j) --> qb(i) ( --> W qb(k) ) qb(j) ijkj
      res=fac*aveqq*qbqb_ijkj
      elseif (ckmallowed(j,-i6) .and.(i==i5).and.(i6/=i)) then
c--- qb(i) qb(j) --> qb(i) qb(j) ( --> W qb(k) ) ijik
      res=fac*aveqq*qbqb_ijik
      elseif (ckmallowed(i,-i5) .and.(j==i5).and.(i5==i6)) then
c--- qb(i) qb(j) --> qb(i) ( --> W qb(j) ) qb(j) ijjj
      res=fac*aveqq*half*qbqb_ijjj
      elseif (ckmallowed(j,-i6) .and.(i==i5).and.(i5==i6)) then 
c--- qb(i) qb(j) --> qb(i) qb(j) ( --> W qb(i) ) ijii
      res=fac*aveqq*half*qbqb_ijii
      endif
      endif
      endif

      endif ! End qflag



************************************************************************
*     UV contributions for GFLAG contribitions are included here       *
************************************************************************
      res=res-ason2pi*subuv*res0

C     change from DRED + add finite coupling constant ren to go to MSbar 
      if (scheme == 'msbr') then
      if (Gflag) then 
         res = res + ason2pi*res0*(-2d0*(cf/2d0)- 2d0*Nc/6d0)
      else
         res = res + ason2pi*res0*(-4d0*(cf/2d0)- 0d0*Nc/6d0)
      endif
      endif

C---- divide by ason2pi as wanted by powheg
      res=res/ason2pi        
      return
      end
     
     
 
