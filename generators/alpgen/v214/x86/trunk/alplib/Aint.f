C**************************************************************************
      subroutine matrix0(flvmlm,posi,impul,hel,result)
C**************************************************************************
C
C Input parameters: IMPUL(4,8) particles four momenta,
C                   HEL array of particle helicities
C Output parameters:   RESULT squared matrix element (modulus) 
C
C
C At present:
C
C   IMPUL(J,NPART) ===   J=1 Energy of the NPART-th particle
C   IMPUL(J,NPART) ===   J=2,3,4 x,y,z components  of the NPART-th particle 
C                                                         three momentum
C
      implicit none
C
      integer nmax,npmax        !maximum number of external particles, processes 
c pippoq
      parameter (nmax=10,npmax=200)    
      integer nlb
      parameter (nlb=4)            
      real*8 impul(4,nmax)      !the contribution to the integral
      complex*16 result
      integer lbls,flvmlm(nmax),posi(2),hel(nmax)
      integer aux(nmax),j2,j3,flvstr(nmax,npmax)
C     
      integer indexoutgoing(nmax) !reporting the order of storage of
C     outgoing momenta
      integer indexingoing(nmax) !reporting the order of storage of  
C     ingoing momenta
      integer ningoing
      integer noutgoing         !number of outgoing particles
      real*8  massoutgoing(nmax) !containing the masses of outcoming 
C     particles
      integer j1
      complex*16 elmat
      integer flaginit,nparticle(2),ns
      save flaginit,ningoing,ns
      common/integra/massoutgoing,
     >     noutgoing,indexoutgoing,indexingoing
C     
C     
      integer optmulti 
      common/options/optmulti        
C     
c     real*8 js
C     
      integer rep(nmax),nqrk,nprt,nglu,nlep,ngb,nphot,nh
      common/process/rep,nqrk,nprt,nglu,nlep,ngb,nphot,nh
C     
      data ningoing/2/,flaginit/0/,ns/0/
C     
      optmulti=1
      if (noutgoing.gt.8) then
         write(6,*)'too many particles'
         stop
      endif
      nparticle(1)=nprt
      nparticle(2)=2
      do j1=1,nparticle(1)
         aux(j1)=flvmlm(j1)
      enddo
      j1=1
      j3=0
      do while(j1.le.ns.and.j3.eq.0)
         j2=1
         do while(j2.le.nparticle(1).and.j3.eq.0)
            if(aux(j2).ne.flvstr(j2,j1)) j3=1
            j2=j2+1
         enddo
         if(j3.eq.0)then
            j3=j1
         else
            j3=0
         endif
         j1=j1+1
      enddo
C     
      if (j3.eq.0) then
         ns=ns+1
         lbls=ns
         do j1=1,nparticle(1)
            flvstr(j1,ns)=aux(j1)
         enddo
      else
         lbls=j3
      endif
C     
      if (ns.gt.npmax) then
         write(6,*)'number of different processes ',ns
         write(6,*)'exceeds',npmax
      endif
C     
      if(flaginit.eq.0) then
         flaginit=1
         if (optmulti.eq.0) call processo
      endif
C     
c     optmulti=0
      if (optmulti.eq.1) call processo_h(flvmlm,posi
     >     ,nparticle)       
      call fillmom_lhc(impul,ningoing) !setting the 
      call spincoulor_h(hel)    !setting spin and coulors
      
C     
c     call ottspin(js)
c     call randa(js)
c     js =1.d0
      call itera(elmat)
C     
c     result=elmat*js
      result=elmat
C     
      return
      end
C*********************************************************************
      subroutine genprm(nmb,nx,pr)
C*********************************************************************
      implicit none
C     
c     Given NMB object (NMB<10) it returns in PR a permutation
c     according to the value of NX (0 < NX < NMB ! +1)
C     
      integer nmx
      parameter (nmx=18)
      integer nmb,nx
      integer pr(nmx)
      integer j1,j2,fac(10)
      data fac/1,2,6,24,120,720,5040,40320,362880,3628800/
      integer aux(nmx),aux1(nmx) 
C     
      do j1=1,nmb
         aux(j1)=j1
      enddo
C     
      nx=nx-1
      do j1=nmb,2,-1
         aux1(j1)=nx/fac(j1-1)+1
         nx=mod(nx,fac(j1-1))
      enddo
      aux1(1)=1
C     
      do j1=1,nmb
         aux(j1)=j1
      enddo
      do j1=nmb,1,-1
         pr(j1)=aux(aux1(j1))
         do j2=aux1(j1),j1-1
            aux(j2)=aux(j2+1)
         enddo
      enddo
C     
      return
      end
C**************************************************************************
         subroutine gendual0q(proc,nglu,col,nx,colfac)
C**************************************************************************
C
C Given: a string PROC(nmax) containing the color representation of the 
C particles (2= adjoint 1,-1 3,3bar, 0 singlet), the number of gluons NGLU,
C the SU3 coulor string COL(2*nmax), a random integers NX and an integer
C flag DLMD it returns: the coulor factor COLFAC and it stores in COLOR
C the selected dual amplitudes. 
C If DLMD=0 a coulor string is selected and if this dual amplitude
C does not contribute to the given SU(3) configuration COL, COLFAC
C is returned 0.d0  1 <= NX <= (NGLU !) * NGLU * 2 -1
C If DLMD=1  a coulor string is selected among the the leading 
C (in 1/N) ones only 1 <= NX <= (NGLU !) * NGLU - 1
C
         implicit none
C
         integer nmax
         parameter (nmax=10)
         integer proc(2*nmax),nx,col(2*nmax),nglu
         real*8 colfac,caux
         integer pr(nmax),ant(2*nmax),auxg(nmax),naux
         integer n1,n2,j1
         integer color(2*nmax)    !coulor string
         common/colore/color
C
         nx=nx+1
         call genprm(nglu,nx,pr)
         naux=0
C
         do j1=1,nmax
          if(proc(j1).eq.2) then
           naux=naux+1
           auxg(naux)=j1
          endif
         enddo
         do j1=1,nglu
          pr(j1)=auxg(pr(j1))
         enddo
C
         do j1=1,2*nmax
          color(j1)=0
         enddo
C
         do j1=1,nglu
          ant(2*j1-1)=col(2*pr(j1)-1)
          color(2*pr(j1)-1)=j1
          ant(2*j1)=col(2*pr(j1))
          color(2*pr(j1))=j1-1
         enddo
         color(2*pr(1))=nglu
C
         call trcgl(ant,nglu,colfac)
C
         return
         end

C**************************************************************************
         subroutine gendual2q(proc,antns,nglu,col,nx,colfac)
C**************************************************************************
C
C Given: a string PROC(nmax) containing the color representation of the 
C particles (2= adjoint 1,-1 3,3bar, 0 singlet), a string ANTNS containing
C the labels of the two principal antennas, the number of gluons NGLU,
C the SU3 coulor string COL(2*nmax), a random integers NX and an integer
C flag DLMD it returns: the coulor factor COLFAC and it stores in COLOR
C the selected dual amplitudes. 
C If DLMD=0 a coulor string is selected and if this dual amplitude
C does not contribute to the given SU(3) configuration COL, COLFAC
C is returned 0.d0  1 <= NX <= (NGLU !) * NGLU * 2 -1
C If DLMD=1  a coulor string is selected among the the leading 
C (in 1/N) ones only 1 <= NX <= (NGLU !) * NGLU - 1
C
         implicit none
C
         integer nmax
         parameter (nmax=10)
         integer proc(2*nmax),antns(2),nx,col(2*nmax),nglu
         real*8 colfac,caux
         integer pr(nmax),ant(2*nmax),auxg(nmax),naux
         integer n1,n2,j1
         integer color(2*nmax)    !coulor string
         common/colore/color
C
         nx=nx+1
         call genprm(nglu,nx,pr)
         naux=0
C
         do j1=1,nmax
          if(proc(j1).eq.2) then
           naux=naux+1
           auxg(naux)=j1
          endif
         enddo
         do j1=1,nglu
          pr(j1)=auxg(pr(j1))
         enddo
C
         do j1=1,2*nmax
          color(j1)=0
         enddo
C
         ant(1)=col(2*antns(1)-1)
         color(2*antns(1)-1)=1
         ant(2)=col(2*antns(1))
         color(2*antns(1))=0
         do j1=2,nglu+1
          ant(2*j1-1)=col(2*pr(j1-1)-1)
          color(2*pr(j1-1)-1)=j1
          ant(2*j1)=col(2*pr(j1-1))
          color(2*pr(j1-1))=j1-1
         enddo
         ant(2*nglu+3)=col(2*antns(2)-1)
         ant(2*nglu+4)=col(2*antns(2))
         color(2*antns(2)-1)=0
         color(2*antns(2))=nglu+1
C
         call qbgluq(ant,nglu,colfac)
C
         return
         end
C**************************************************************************
      subroutine gendual4q(proc,antns,nglu,col,nx,colfac)
C**************************************************************************
C     
C     Given: a string PROC(nmax) containing the color representation of
c     the 
C     particles (2= adjoint 1,-1 3,3bar, 0 singlet), a string ANTNS
c     containing
C     the labels of the two principal antennas, the number of gluons
c     NGLU,
C     the SU3 coulor string COL(2*nmax), a random integers NX and an
c     integer
C     flag DLMD it returns: the coulor factor COLFAC and it stores in
c     COLOR
C     the selected dual amplitudes. 
C     If DLMD=0 a coulor string is selected and if this dual amplitude
C     does not contribute to the given SU(3) configuration COL, COLFAC
C     is returned 0.d0  1 <= NX <= (NGLU !) * NGLU * 2 -1
C     If DLMD=1  a coulor string is selected among the the leading 
C     (in 1/N) ones only 1 <= NX <= (NGLU !) * NGLU - 1
C     
      implicit none
C     
      integer nmax
      parameter (nmax=10)
      integer proc(2*nmax),antns(4),nx,col(2*nmax),nglu
      real*8 colfac,caux
      integer pr(nmax),ant(2*nmax),auxg(nmax),naux
      integer n1,n2,j1
      integer color(2*nmax)     !coulor string
      common/colore/color
C     
      n1=mod(nx,nglu+1)
      n2=nglu-n1
      nx=nx/(nglu+1)+1
      call genprm(nglu,nx,pr)
      naux=0
C     
      do j1=1,nmax
         if(proc(j1).eq.2) then
            naux=naux+1
            auxg(naux)=j1
         endif
      enddo
      do j1=1,nglu
         pr(j1)=auxg(pr(j1))
      enddo
C     
      do j1=1,2*nmax
         color(j1)=0
      enddo
C     
      ant(1)=col(2*antns(1)-1)
      color(2*antns(1)-1)=1
      ant(2)=col(2*antns(1))
      color(2*antns(1))=0
      do j1=2,n1+1
         ant(2*j1-1)=col(2*pr(j1-1)-1)
         color(2*pr(j1-1)-1)=j1
         ant(2*j1)=col(2*pr(j1-1))
         color(2*pr(j1-1))=j1-1
      enddo
      ant(2*n1+3)=col(2*antns(2)-1)
      ant(2*n1+4)=col(2*antns(2))
      color(2*antns(2)-1)=0
      color(2*antns(2))=n1+1
C     
      call qbgluq(ant,n1,colfac)
      if(abs(colfac).lt.1.d-20) return
      caux=colfac
C     
      ant(1)=col(2*antns(3)-1)
      color(2*antns(3)-1)=n1+2
      ant(2)=col(2*antns(3))
      color(2*antns(3))=0
      do j1=n1+4,n1+n2+3
         ant(2*(j1-n1-2)-1)=col(2*pr(j1-3)-1)
         color(2*pr(j1-3)-1)=j1-1
         ant(2*(j1-n1-2))=col(2*pr(j1-3))
         color(2*pr(j1-3))=j1-2
      enddo
      ant(2*n2+3)=col(2*antns(4)-1)
      color(2*antns(4)-1)=0
      ant(2*n2+4)=col(2*antns(4))
      color(2*antns(4))=n1+n2+2
C     
      call qbgluq(ant,n2,colfac)
      colfac=colfac*caux
C     
      return
      end
C**************************************************************************
      subroutine smfx3(ir,nn,ntot)
C**************************************************************************
C
C given a random integer IR it returns in NN(j) (j=1,2,3) three integers
C whose sum is NTOT. 1=<IR<=(NTOT+1)(NTOT+2)/2 ; 0 =< NN(j )=< NTOT
C
      implicit none
C     
      integer ir,nn(3),ntot
      integer j1
C     
      if(ir.gt.(ntot+1)*(ntot+2)/2) then
         write(6,*)'wrong IR in SMFX3'
      endif
      nn(1)=-1
      j1=ntot+1
      do while(ir.gt.0)
         ir=ir-j1
         j1=j1-1
         nn(1)=nn(1)+1
      enddo
      nn(2)=-ir
      nn(3)=ntot-nn(1)-nn(2)
C     
      return
      end
C**************************************************************************
      subroutine gendual6q(proc,antns,nglu,col,nx,colfac)
C**************************************************************************
C     
C Given: a string PROC(nmax) containing the color representation of the 
C particles (2= adjoint 1,-1 3,3bar, 0 singlet), a string ANTNS containing
C the labels of the three principal antennas, the number of gluons NGLU,
C the SU3 coulor string COL(2*nmax), a random integers NX and an integer
C flag DLMD it returns: the coulor factor COLFAC and it stores in COLOR
C the selected dual amplitudes. 
C If DLMD=0 a coulor string is selected and if this dual amplitude
C does not contribute to the given SU(3) configuration COL, COLFAC
C is returned 0.d0   0 <= NX <= (NGLU !) * (NGLU+1) (NGLU+2)/2 * 6 -1
C If DLMD=1  a coulor string is selected among the the leading 
C (in 1/N) ones only  0 <= NX <= (NGLU !) * (NGLU+1) (NGLU+2)/2 - 1
C
      implicit none
C     
      integer nmax
      parameter (nmax=10)
      integer proc(2*nmax),antns(6),nx,col(2*nmax),nglu
      real*8 colfac,caux
      integer pr(nmax),ant(2*nmax),auxg(nmax),naux
      integer nn(3),j1
      integer color(2*nmax)     !coulor string
      common/colore/color
C     
      call smfx3(mod(nx,(nglu+1)*(nglu+2)/2)+1,nn,nglu)
      nx=2*nx/(nglu+1)/(nglu+2)+1
      call genprm(nglu,nx,pr)
      naux=0
C     
      do j1=1,nmax
         if(proc(j1).eq.2) then
            naux=naux+1
            auxg(naux)=j1
         endif
      enddo
      do j1=1,nglu
         pr(j1)=auxg(pr(j1))
      enddo
C     
      do j1=1,2*nmax
         color(j1)=0
      enddo
C     
      ant(1)=col(2*antns(1)-1)
      color(2*antns(1)-1)=1
      ant(2)=col(2*antns(1))
      color(2*antns(1))=0
      do j1=2,nn(1)+1
         ant(2*j1-1)=col(2*pr(j1-1)-1)
         color(2*pr(j1-1)-1)=j1
         ant(2*j1)=col(2*pr(j1-1))
         color(2*pr(j1-1))=j1-1
      enddo
      ant(2*nn(1)+3)=col(2*antns(2)-1)
      ant(2*nn(1)+4)=col(2*antns(2))
      color(2*antns(2)-1)=0
      color(2*antns(2))=nn(1)+1
C     
      call qbgluq(ant,nn(1),colfac)
      if(abs(colfac).lt.1.d-20) return
      caux=colfac
C     
      ant(1)=col(2*antns(3)-1)
      color(2*antns(3)-1)=nn(1)+2
      ant(2)=col(2*antns(3))
      color(2*antns(3))=0
      do j1=nn(1)+4,nn(1)+nn(2)+3
         ant(2*(j1-nn(1)-2)-1)=col(2*pr(j1-3)-1)
         color(2*pr(j1-3)-1)=j1-1
         ant(2*(j1-nn(1)-2))=col(2*pr(j1-3))
         color(2*pr(j1-3))=j1-2
      enddo
      ant(2*nn(2)+3)=col(2*antns(4)-1)
      color(2*antns(4)-1)=0
      ant(2*nn(2)+4)=col(2*antns(4))
      color(2*antns(4))=nn(1)+nn(2)+2
C     
      call qbgluq(ant,nn(2),colfac)
      if(abs(colfac).lt.1.d-20) return
      caux=colfac*caux
C     
      ant(1)=col(2*antns(5)-1)
      color(2*antns(5)-1)=nn(1)+nn(2)+3
      ant(2)=col(2*antns(5))
      color(2*antns(5))=0
      do j1=nn(1)+nn(2)+1,nn(1)+nn(2)+nn(3)
         ant(2*(j1-nn(1)-nn(2))+1)=col(2*pr(j1)-1)
         color(2*pr(j1)-1)=j1+3
         ant(2*(j1-nn(1)-nn(2))+2)=col(2*pr(j1))
         color(2*pr(j1))=j1+2
      enddo
      ant(2*nn(3)+3)=col(2*antns(6)-1)
      color(2*antns(6)-1)=0
      ant(2*nn(3)+4)=col(2*antns(6))
      color(2*antns(6))=nn(1)+nn(2)+nn(3)+3
C     
      call qbgluq(ant,nn(3),colfac)
      colfac=colfac*caux
C     
      return
      end
C****************************************************************************
      subroutine selcol(flvmlm,posi,rnd1,rnd2,allowed)
c     select colours
C****************************************************************************
C
      implicit none
C     
      integer nmax
      parameter (nmax=10)       !maximum number of external particles. 
      integer flvmlm(nmax),posi(2)
      integer color(2*nmax)     !coulor string
      common/colore/color
      integer nlb
      parameter (nlb=4)            
      integer nrep
      parameter (nrep=100)
      integer allowed           !0 if coulor flow possible 1 otherwise
      integer j1,m1,m2,count,col(nmax),cola(nmax)
      integer gell(2,8),spin
      integer spin1(8),spin2(8),nqq
      integer rep(nmax),nqrk,nprt,nglu,nlep,ngb,nphot,nh
      real*8 rnd1,rnd2          ! random numbers to select coulor
      real*8 n1,n2,qq(3,3)
      common/process/rep,nqrk,nprt,nglu,nlep,ngb,nphot,nh
      data  gell/1,1,1,2,1,3,2,1,2,2,2,3,3,1,3,2/
      data  spin1/0,1,2,-1,0,1,-2,-1/,spin2/0,-1,0,1,0,1,0,-1/
      data qq/1,4,7,2,5,8,3,6,1/
c     data repa/-1,0,-1,1,0,1,2,2,2,2,   !d nu b ub e- bb g g g g   incoming
c     >         -1,0,-1,-1,1,0,1,1,2,2,  !d nu c b ub e- cb bb g g  incoming
c     >         -1,-1,0,-1,1,1,0,1,2,2   !d d nu b ub db e- bb g g  incoming
c     >             ,970*10  /
C
         integer flvax(nmax)
         common/flv/flvax
C
         do j1=1,nmax
           flvax(j1)=flvmlm(j1)
         enddo      
C
      call flvtocol(flvmlm,posi)
C     
      n2=8**nglu
      n1=3**nqrk
      nqq=nqrk/2
      if (nqrk.eq.0.) then
      elseif (nqrk.eq.2.) then
      elseif (nqrk.eq.4) then
      elseif (nqrk.eq.6) then
      else
         write(6,*)'wrong flavour assignment',flvmlm
         stop
      endif
C     
      count=0
      m1=int(rnd1*n1)+1
      do j1=1,nqrk
         count=count+1
         cola(count)=mod(m1,3)+1
         m1=m1/3
      enddo
      m2=int(rnd2*n2)+1
      do j1=nqrk+1,nqrk+nglu
         count=count+1
         cola(count)=mod(m2,8)+1
         m2=m2/8
      enddo
      do j1=1,nqq
         col(j1)=qq(cola(j1),cola(j1+nqq))
      enddo
      do j1=nqq+1,nqq+nglu
         col(j1)=cola(nqq+j1)
      enddo
C     
      allowed=0
      spin=0
      do j1=1,nqq+nglu
         spin=spin+spin1(col(j1))
      enddo
      if(spin.ne.0) then
         allowed=1
         return
      endif
      spin=0
      do j1=1,nqq+nglu
         spin=spin+spin2(col(j1))
      enddo
      if(spin.ne.0) then
         allowed=1
         return
      endif
C     
      count=0
      do j1=1,nprt
         if (rep(j1).eq.0) then
            color(2*j1-1)=0
            color(2*j1)=0
         elseif (rep(j1).eq.-1) then
            count=count+1
            color(2*j1-1)=cola(count)
            color(2*j1)=0
         elseif (rep(j1).eq.1) then
            count=count+1
            color(2*j1-1)=0
            color(2*j1)=cola(count)
         elseif (rep(j1).eq.2) then
            count=count+1
            color(2*j1-1)=gell(1,cola(count))
            color(2*j1)=gell(2,cola(count))
         else
            write(6,*)'wrong coulor representation in INCOLBBW'
            stop
         endif
      enddo
C     
      return
      end                                                                 
C********************************************************
      subroutine qbgluq(ant,nglu,colfac)
C********************************************************
C     
      implicit none
C     
      integer nmax
      parameter (nmax=10)
      integer ant(2*nmax),nglu
      real*8 colfac
      integer flg,j1,j2,j4
      real*8 gel(3,3,3,3),aux(3),res(3)
      data gel/81*0./
      data flg/0/
      save flg,gel
C     
      if (flg.eq.0) then
         flg=1
         gel(1,1,1,1)=1./sqrt(2.d0)
         gel(1,1,3,3)=-1./sqrt(2.d0)
         gel(1,2,1,2)=1.
         gel(1,3,1,3)=1.
         gel(2,1,2,1)=1.
         gel(2,2,1,1)=-1./sqrt(6.d0)
         gel(2,2,2,2)=sqrt(2.d0/3.d0)
         gel(2,2,3,3)=-1./sqrt(6.d0)
         gel(2,3,2,3)=1.
         gel(3,1,3,1)=1.
         gel(3,2,3,2)=1.
      endif
C     
      do j1=1,3
         res(j1)=0
         aux(j1)=0
      enddo
      aux(ant(1))=1.
      do j4=2,nglu+1
         do j1=1,3
            res(j1)=0.
            do j2=1,3
               res(j1)=res(j1)+aux(j2)*
     >              gel(ant(2*j4-1),ant(2*j4),j1,j2)
            enddo
         enddo
         do j1=1,3
            aux(j1)=res(j1)
         enddo
      enddo
C     
      if(nglu.eq.0) res(ant(1))=1.
C     
      colfac=res(ant(2*nglu+4))
C     
      return
      end
C********************************************************
      subroutine trcgl(ant,nglu,colfac)
C********************************************************
C     
      implicit none
C     
      integer nmax
      parameter (nmax=10)
      integer ant(2*nmax),nglu
      real*8 colfac
      integer flg,j1,j2,j3,j4
      real*8 gel(3,3,3,3),aux(3,3),res(3,3)
      data gel/81*0./
      data flg/0/
      save flg,gel
C     
      if (flg.eq.0) then
         flg=1
         gel(1,1,1,1)=1./sqrt(2.d0)
         gel(1,1,3,3)=-1./sqrt(2.d0)
         gel(1,2,1,2)=1.
         gel(1,3,1,3)=1.
         gel(2,1,2,1)=1.
         gel(2,2,1,1)=-1./sqrt(6.d0)
         gel(2,2,2,2)=sqrt(2.d0/3.d0)
         gel(2,2,3,3)=-1./sqrt(6.d0)
         gel(2,3,2,3)=1.
         gel(3,1,3,1)=1.
         gel(3,2,3,2)=1.
      endif
C     
      do j1=1,3
         do j2=1,3
            aux(j1,j2)=gel(ant(1),ant(2),j1,j2)
         enddo
      enddo
      do j4=2,nglu
         do j1=1,3
            do j2=1,3
               res(j1,j2)=0.
               do j3=1,3
                  res(j1,j2)=res(j1,j2)+aux(j1,j3)*
     >                 gel(ant(2*j4-1),ant(2*j4),j3,j2)
               enddo
            enddo
         enddo
         do j1=1,3
            do j2=1,3
               aux(j1,j2)=res(j1,j2)
            enddo
         enddo
      enddo
C     
      colfac=0.
      do j1=1,3
         colfac=colfac+res(j1,j1)
      enddo
C     
      return
      end
C********************************************************
      subroutine prm(col,rnd,all,nglu,colfac)      
C********************************************************
C     
      implicit none
C     
      integer nmax
      parameter (nmax=10)
C     
      integer col(2*nmax),ant(2*nmax),all,nglu
      real*8 rnd(nmax-1),colfac,aux1
      integer j1,j2,j3,aux(nmax),cmb(2)
C     
      integer color(2*nmax)     !coulor string
      common/colore/color
C     
      do j1=1,nglu
         aux(j1)=j1
      enddo
      do j1=nglu-1,1,-1
         do j2=2,j1
            if (rnd(j2).lt.rnd(j2-1)) then
               aux1=rnd(j2-1)
               rnd(j2-1)=rnd(j2)
               rnd(j2)=aux1
               aux1=aux(j2-1)
               aux(j2-1)=aux(j2)
               aux(j2)=aux1
            endif
         enddo
      enddo
C     
      do j1=1,nglu
         j2=2*j1-1
         j3=2*aux(j1)-1
         ant(j2)=col(j3)
         ant(j2+1)=col(j3+1)
      enddo
C     
      all=1
      j1=1
      cmb(1)=ant(1)
      cmb(2)=ant(2)
      do while(all.eq.1.and.j1.lt.2*nglu-1)
         j1=j1+2
         call prdglu(cmb,ant(j1))
         if(cmb(1).eq.0) all=0
      enddo
C     
      if(all.eq.1) then
         call trcgl(ant,nglu,colfac)
         if(abs(colfac).lt.1.d-10) then
            all=0
            return
         else 
            do j1=1,nglu
               color(2*aux(j1)-1)=j1
               color(2*aux(j1))=j1+1
               if(j1.eq.nglu) color(2*j1)= 1
            enddo
         endif
      endif
C     
      return
      end
C***************************************************************
      subroutine prdglu(cmb,cmb1)
C***************************************************************
C     
      implicit none
C     
      integer cmb(2),cmb1(2)
C     
      if(cmb(1).eq.cmb(2)) then
         cmb(1)=cmb1(1)
         cmb(2)=cmb1(2)
         return
      elseif(cmb1(1).eq.cmb1(2)) then
         return
      elseif(cmb(2).eq.cmb1(1)) then
         cmb(2)=cmb1(2)
      else
         cmb(1)=0
         cmb(2)=0
      endif
C     
      return
      end
C**********************************************************************
      subroutine flvtocol(flvmlm,posi)
C**********************************************************************
C     
      implicit none
C     
      integer nmax
      parameter (nmax=10)
      integer flvmlm(nmax),posi(2)
      integer rep(nmax),nprt,nglu,nqrk,nlep,ngb,nphot,nh
      integer flv(nmax),alp(nmax),j1,j2,ist1,ist2
      integer tabmlm(-100:100,2),tmp,dummy
      data tabmlm/70*99,
     >     6*99,0,3*99,
     >     8*99,0,0,
     >     4*99,-1,-1,-1,-1,-1,-1,
     >     2,
     >     1,1,1,1,1,1,4*99,
     >     0,0,8*99,
     >     99,4*0,5*99,      ! (*,1)
     >     70*99,                
     >     70*0,                
     >     6*0,28,3*0,
     >     8*0,3,4,
     >     4*0,9,10,5,6,1,2,
     >     30,
     >     14,13,18,17,22,21,4*0,
     >     16,15,8*0,
     >     0,29,26,27,25,5*0,
     >     70*0/                ! (*,2)
C     
      save tabmlm
C     
      common/process/rep,nqrk,nprt,nglu,nlep,ngb,nphot,nh
C     
      nglu=0
      nqrk=0
      nprt=0
      nphot=0
      nlep=0
      ngb=0
      nh=0
 1    nprt= nprt+1
      flv(nprt)=flvmlm(nprt)
      if(flv(nprt).ne.1001) then
         if(abs(tabmlm(flvmlm(nprt),1)).eq.2) then
            nglu=nglu+1
         elseif(abs(tabmlm(flvmlm(nprt),1)).eq.1) then
            nqrk=nqrk+1
         elseif(abs(tabmlm(flvmlm(nprt),1)).eq.0) then
            dummy=abs(tabmlm(flvmlm(nprt),2))
            if(dummy.eq.29) then
               nphot=nphot+1
            elseif(dummy.gt.25) then
               ngb=ngb+1
            elseif(dummy.eq.25) then
               nh=nh+1
            elseif(dummy.gt.0) then
               nlep=nlep+1
            else
               write(6,*)'wrong TABMLM in FLVTOCOL'
            endif
         else
            write(6,*)'WRONG FLAVOUR ASSIGNMENT IN FLVTOCOL'
            write(6,*)tabmlm(flvmlm(nprt),1),flvmlm(nprt),nprt
            stop
         endif
         if(nprt.ne.nmax) goto 1
      else
         nprt=nprt-1         
      endif
      do j1=1,2
         flv(posi(j1))=-flv(posi(j1)) 
      enddo
      do j1=1,nprt
         alp(j1)=tabmlm(flv(j1),2)
      enddo
      ist1=flv(posi(1))
      ist2=flv(posi(2))
C     
      do j1=nprt-1,1,-1
         do j2=1,j1
            if (alp(j2).gt.alp(j2+1))then
               tmp=alp(j2+1)
               alp(j2+1)=alp(j2)
               alp(j2)=tmp
               tmp=flv(j2+1)
               flv(j2+1)=flv(j2)
               flv(j2)=tmp
               tmp=flvmlm(j2+1)
               flvmlm(j2+1)=flvmlm(j2)
               flvmlm(j2)=tmp
            endif
         enddo
      enddo
C     
      do j1=1,nprt
         rep(j1)=tabmlm(flv(j1),1)
         if(rep(j1).eq.9) then
            write(6,*)'something wrong in FLVTOCOL'
            stop
         endif
      enddo
C     
      j1=1
      do while(ist1.ne.flv(j1).and.ist2.ne.flv(j1))
         j1=1+j1
      enddo
      posi(1)=j1
      if(flv(j1).eq.ist1) then
         ist1=ist2
      elseif(flv(j1).eq.ist2) then
      else
         write(6,*)'failure in FLVTOCOL'
      endif
      j1=j1+1
      do while(ist1.ne.flv(j1))
         j1=1+j1
      enddo
      posi(2)=j1
C     
      return
      end
      

C**************************************************************************
         subroutine selspin(rspin,hel)
C**************************************************************************
C
C Input parameters:   RSPIN random number to select spin configuration
C Output parameters:  HEL helicity configuration 
C
         implicit none
C
         integer nmax
         parameter (nmax=10)
C
         real*8 rspin
         integer hel(nmax)
C
         integer rep(nmax),nqrk,nprt,nglu,nlep,ngb,nphot,nh
         common/process/rep,nqrk,nprt,nglu,nlep,ngb,nphot,nh
         integer nspmax
         parameter (nspmax=1024)
         integer j1,j2
C
         nh=nprt-nglu-nqrk-ngb-nlep-nphot
C
         j2=rspin*2**(nglu+nqrk+nphot+nlep)*3**ngb
         do j1=1,nqrk+nlep
           hel(j1)=mod(j2,2)
           if(hel(j1).eq.0)hel(j1)=-1
           j2=j2/2       
         enddo
         do j1=nqrk+nlep+1,nqrk+nlep+nh
            hel(j1)=0
         enddo
         do j1=nqrk+nlep+nh+1,nqrk+nlep+nh+ngb
           hel(j1)=mod(j2,3)-1
           j2=j2/3       
         enddo
         do j1=nqrk+nlep+nh+ngb+1,nqrk+nlep+nh+ngb+nglu+nphot
            hel(j1)=mod(j2,2)
            if(hel(j1).eq.0)hel(j1)=-1
            j2=j2/2       
         enddo
C
         return
         end



C*********************************************************************
        subroutine fltspn(flvmlm,hel,posi,allowed)
C*********************************************************************
C
        implicit none
C
        integer nmax
        parameter (nmax=10)
        integer flvmlm(nmax), hel(nmax),allowed,posi(2)
C
        integer j1,hs(-12:12),j2,nf(-12:12),flv(nmax)
C
        do j1=1,nmax
         flv(j1)=flvmlm(j1)
        enddo
        flv(posi(1))=-flv(posi(1))
        flv(posi(2))=-flv(posi(2))
C
        do j1=-12,12
         hs(j1)=0
         nf(j1)=0
        enddo
        j1=1
        do j1=1,nmax
           j2=flv(j1)
           if(j2.eq.1001) goto 111
           if(abs(j2).le.12) then
              hs(j2)=hs(j2)+hel(j1)
              nf(j2)=nf(j2)+1
           endif
           
        enddo
C
 111    allowed=1
        j2=abs(abs(nf(1)-nf(-1))+abs(nf(2)-nf(-2)))
        if(mod(j2,2).ne.0) then
         write(*,*)'wrong flavour assignment in fltspn'
         stop
        endif
        j2=j2/2
        if(nf(1).gt.nf(-1)) then
         if(hs(1)-hs(-1).ne.j2) return
        else
         if(hs(-1)-hs(1).ne.j2) return
        endif
        if(nf(2).gt.nf(-2)) then
         if(hs(2)-hs(-2).ne.j2) return
        else
         if(hs(-2)-hs(2).ne.j2) return
        endif
c
        j2=abs(abs(nf(3)-nf(-3))+abs(nf(4)-nf(-4)))
        if(mod(j2,2).ne.0) then
         write(*,*)'wrong flavour assignment in fltspn'
         stop
        endif
        j2=j2/2
        if(nf(3).gt.nf(-3)) then
         if(hs(3)-hs(-3).ne.j2) return
        else
         if(hs(-3)-hs(3).ne.j2) return
        endif
        if(nf(4).gt.nf(-4)) then
         if(hs(4)-hs(-4).ne.j2) return
        else
         if(hs(-4)-hs(4).ne.j2) return
        endif
c
        j2=abs(abs(nf(7)-nf(-7))+abs(nf(8)-nf(-8)))
        if(mod(j2,2).ne.0) then
         write(*,*)'wrong flavour assignment in fltspn'
         stop
        endif
        j2=j2/2
        if(nf(7).gt.nf(-7)) then
         if(hs(7)-hs(-7).ne.j2) return
        else
         if(hs(-7)-hs(7).ne.j2) return
        endif
        if(nf(8).gt.nf(-8)) then
         if(hs(8)-hs(-8).ne.j2) return
        else
         if(hs(-8)-hs(8).ne.j2) return
        endif
c
        j2=abs(abs(nf(9)-nf(-9))+abs(nf(10)-nf(-10)))
        if(mod(j2,2).ne.0) then
         write(*,*)'wrong flavour assignment in fltspn'
         stop
        endif
        j2=j2/2
        if(nf(9).gt.nf(-9)) then
         if(hs(9)-hs(-9).ne.j2) return
        else
         if(hs(-9)-hs(9).ne.j2) return
        endif
        if(nf(10).gt.nf(-10)) then
         if(hs(10)-hs(-10).ne.j2) return
        else
         if(hs(-10)-hs(10).ne.j2) return
        endif
c
        j2=abs(abs(nf(11)-nf(-11))+abs(nf(12)-nf(-12)))
        if(mod(j2,2).ne.0) then
         write(*,*)'wrong flavour assignment in fltspn'
         stop
        endif
        j2=j2/2
        if(nf(12).ne.hs(12)) return
        if(nf(-12).ne.hs(-12)) return
        if(nf(11).gt.nf(-11)) then
         if(hs(11)-hs(-11).ne.j2) return
        else
         if(hs(-11)-hs(11).ne.j2) return
        endif
C
        allowed=0
C
        return
        end
c***************************************************************************
      subroutine vtop(ptop,mt,mb,mw,vt)
c***************************************************************************
c
C  Input  PTOP(4) : top four momentum (E_t, p_x,p_y,p_z)
C         MT,MB,MW  top, bottom and W mass
c  Output VT(4) top-bar anti-spinor
c
      implicit none
      real*8 ptop(4),mt,mb,mw
      complex*16 vt(4)
C
      complex*16  gam(4,4,5),vb(4),ubf(4),vfb(4),pl(4,4),gamp(4,4,4) 
      complex*16  gams(4,4,5),pg(4,4) 
      complex*16  eps(4),tmp,epsl(4,4),ptsl(4,4),tmpu(4),iden(4,4)
      real*8 pf(4),pfb(4),pb(4),r1,sp,mtlnz(4),analitic,ew,pw
      real*8 gvt,gat,cosvma
      common/topcoup/cosvma
      integer flag,j1,j2,j3,j4
C
      common/pdec/pf,pfb,pb
C
      data gam/(1.,0.),(0.,0),(0.,0.),(0.,0.),(0.,0.)
     >     ,(1.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     >  (-1.,0.),(0.,0.),(0.,0.),(0.,0.), (0.,0.),(-1.,0.),          !gamma0
     >  (0.,0.),(0.,0.),(0.,0.),(1.,0.),(0.,0.),
     >  (0.,0.),(1.,0.),(0.,0.),(0.,0.),(-1.,0.),
     >  (0.,0.),(0.,0.),(-1.,0.),(0.,0.),(0.,0.),(0.,0.),             !gamma1
     >  (0.,0.),(0.,0.),(0.,0.),(0.,-1.),(0.,0.), 
     >  (0.,0.),(0.,1.),(0.,0.),(0.,0.),(0.,1.), 
     >  (0.,0.),(0.,0.),(0.,-1.),(0.,0.),(0.,0.), (0.,0.),           !gamma2
     >  (0.,0.),(0.,0.),(1.,0.),(0.,0.),(0.,0.),
     >  (0.,0.),(0.,0.),(-1.,0.),(-1.,0.),(0.,0.),
     >  (0.,0.),(0.,0.),(0.,0.),(1.,0.),(0.,0.),(0.,0.),              !gamma3
     >  (0.,0.),(0.,0.),(1.,0.),(0.,0.),(0.,0.), 
     >  (0.,0.),(0.,0.),(1.,0.),(1.,0.),(0.,0.),
     >   (0.,0.),(0.,0.),(0.,0.),(1.,0.),(0.,0.),(0.,0.)/           !gamma5
        data iden/(1.,0.), (0.,0.), (0.,0.), (0.,0.),
     *   (0.,0.), (1.,0.),(0.,0.), (0.,0.), (0.,0.), (0.,0.),
     *   (1.,0.), (0.,0.),(0.,0.), (0.,0.), (0.,0.), (1.,0.)/
        data flag/0/,mtlnz/1.d0,-1.d0,-1.d0,-1.d0/
C
      save gamp,flag,mtlnz,gam,gams
C
      if(flag.eq.0) then
* top-W couplings
        gvt= cosvma+sqrt(1.d0-cosvma**2)
        gat=-cosvma+sqrt(1.d0-cosvma**2)
c     repeat the loop twice, to get rid of g77 3.1 compilation bug on
C     linux PC's
        do flag=1,2
          do j1=1,4
            do j2=1,4
              pl(j1,j2)=iden(j1,j2)-gam(j2,j1,5)
              pg(j1,j2)=gvt*iden(j1,j2)+gat*gam(j2,j1,5)
            enddo
          enddo
          do j1=1,4
            do j2=1,4
              do j3=1,4
                gamp(j1,j2,j3)=(0.d0,0.d0)
                gams(j1,j2,j3)=(0.d0,0.d0)
                do j4=1,4
                  gamp(j1,j2,j3)=gamp(j1,j2,j3)+gam(j4,j1,j3)*pl(j4,j2)
                  gams(j1,j2,j3)=gams(j1,j2,j3)+gam(j4,j1,j3)*pg(j4,j2)
                enddo
              enddo
            enddo
          enddo
        enddo
        flag=1
      endif
C
      call tdeckin(ptop,mt,mb,mw,pf,pfb,pb)
C
c debug
c      r1=(pf(1)+pfb(1))**2
c      do j1=2,4
c        r1=r1 - (pf(j1)+pfb(j1))**2
c      enddo
c      write(*,*) 'Aint dilepton mass: ',sqrt(r1)
c      r1=(ptop(1))**2
c      do j1=2,4
c        r1=r1 - (ptop(j1))**2
c      enddo
c      write(*,*) 'Aint top mass: ',sqrt(r1)
c end debug

      call randa(r1)
      if(r1.gt.0.5d0) then
        sp=0.51d0
      else
        sp=-0.51d0
      endif      
C
      call fermionsources_d(sp,vb,pb,mb)
      sp=0.49d0
      call fermionbarsources_d(sp,ubf,pf,0.d0)
      sp=0.51d0
      call fermionsources_d(sp,vfb,pfb,0.d0)
C
      do j1=1,4
        eps(j1)=(0.d0,0.d0)
        do j2=1,4
          tmp=0.d0
          do j3=1,4
            tmp=tmp+gamp(j2,j3,j1)*vfb(j3)
          enddo
          eps(j1)=eps(j1)+ubf(j2)*tmp
        enddo
      enddo
C
      do j1=1,4
        do j2=1,4
          epsl(j1,j2)=(0.d0,0.d0)
          do j3=1,4
            epsl(j1,j2)=epsl(j1,j2)+mtlnz(j3)*eps(j3)*gams(j1,j2,j3)
          enddo
        enddo
      enddo
C
      do j1=1,4
        tmpu(j1)=(0.d0,0.d0)
        do j2=1,4
          tmpu(j1)=tmpu(j1)+vb(j2)*epsl(j1,j2)
        enddo
      enddo
C
      do j1=1,4
        do j2=1,4
          ptsl(j1,j2)=(0.d0,0.d0)
          do j3=1,4
            ptsl(j1,j2)=ptsl(j1,j2)+mtlnz(j3)*ptop(j3)*gam(j2,j1,j3)
          enddo
          ptsl(j1,j2)=ptsl(j1,j2)-mt*iden(j1,j2)
        enddo
      enddo
C
      do j1=1,4
        vt(j1)=(0.d0,0.d0)
        do j2=1,4
          vt(j1)=vt(j1)+tmpu(j2)*ptsl(j1,j2)
        enddo
      enddo
C
      ew=0.5d0*(mt**2+mw**2-mb**2)/mt
      pw=sqrt(ew**2-mw**2)
      analitic= (128*mt**2*((gAt**2 + gVt**2)*(3.*ew**2 - pw**2)*mt - 
     -      3.*(ew*(gAt**2 + gVt**2) + (-gAt**2 + gVt**2)*mb)*mw**2
     -      ))/6.
      analitic=sqrt(4.d0*mt/analitic)
      do j1=1,4
        vt(j1)=analitic*vt(j1)
      enddo
C
      return
      end
c***************************************************************************
      subroutine utopb(ptop,mt,mb,mw,utb)
c***************************************************************************
c
C  Input  PTOP(4) : top four momentum (E_t, p_x,p_y,p_z)
C         MT,MB,MW  top, bottom and W mass
c  Output UT(4) top spinor
c
      implicit none
      real*8 ptop(4),mt,mb,mw
      complex*16 utb(4)
C
      complex*16  gam(4,4,5),ubb(4),ubf(4),vfb(4),pl(4,4),gamp(4,4,4) 
      complex*16  gams(4,4,5),pg(4,4) 
      complex*16  eps(4),tmp,epsl(4,4),ptsl(4,4),tmpu(4),iden(4,4)
      real*8 pf(4),pfb(4),pb(4),r1,sp,mtlnz(4),analitic,ew,pw
      real*8 gvt,gat,cosvma
      common/topcoup/cosvma
      integer flag,j1,j2,j3,j4
C
      common/pdec/pf,pfb,pb
C
      data gam/(1.,0.),(0.,0),(0.,0.),(0.,0.),(0.,0.)
     >     ,(1.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     >  (-1.,0.),(0.,0.),(0.,0.),(0.,0.), (0.,0.),(-1.,0.),          !gamma0
     >  (0.,0.),(0.,0.),(0.,0.),(1.,0.),(0.,0.),
     >  (0.,0.),(1.,0.),(0.,0.),(0.,0.),(-1.,0.),
     >  (0.,0.),(0.,0.),(-1.,0.),(0.,0.),(0.,0.),(0.,0.),             !gamma1
     >  (0.,0.),(0.,0.),(0.,0.),(0.,-1.),(0.,0.), 
     >  (0.,0.),(0.,1.),(0.,0.),(0.,0.),(0.,1.), 
     >  (0.,0.),(0.,0.),(0.,-1.),(0.,0.),(0.,0.), (0.,0.),           !gamma2
     >  (0.,0.),(0.,0.),(1.,0.),(0.,0.),(0.,0.),
     >  (0.,0.),(0.,0.),(-1.,0.),(-1.,0.),(0.,0.),
     >  (0.,0.),(0.,0.),(0.,0.),(1.,0.),(0.,0.),(0.,0.),              !gamma3
     >  (0.,0.),(0.,0.),(1.,0.),(0.,0.),(0.,0.), 
     >  (0.,0.),(0.,0.),(1.,0.),(1.,0.),(0.,0.),
     >   (0.,0.),(0.,0.),(0.,0.),(1.,0.),(0.,0.),(0.,0.)/           !gamma5
        data iden/(1.,0.), (0.,0.), (0.,0.), (0.,0.),
     *   (0.,0.), (1.,0.),(0.,0.), (0.,0.), (0.,0.), (0.,0.),
     *   (1.,0.), (0.,0.),(0.,0.), (0.,0.), (0.,0.), (1.,0.)/
        data flag/0/,mtlnz/1.d0,-1.d0,-1.d0,-1.d0/
C
      save gamp,flag,mtlnz,gam,gams
C
      if(flag.eq.0) then
* top-W couplings
        gvt= cosvma+sqrt(1.d0-cosvma**2)
        gat=-cosvma+sqrt(1.d0-cosvma**2)
c     repeat the loop twice, to get rid of g77 3.1 compilation bug on
C     linux PC's
        do flag=1,2
          do j1=1,4
            do j2=1,4
              pl(j1,j2)=iden(j1,j2)-gam(j2,j1,5)
              pg(j1,j2)=gvt*iden(j1,j2)+gat*gam(j2,j1,5)
            enddo
          enddo
          do j1=1,4
            do j2=1,4
              do j3=1,4
                gamp(j1,j2,j3)=(0.d0,0.d0)
                gams(j1,j2,j3)=(0.d0,0.d0)
                do j4=1,4
                  gamp(j1,j2,j3)=gamp(j1,j2,j3)+gam(j4,j1,j3)*pl(j4,j2)
                  gams(j1,j2,j3)=gams(j1,j2,j3)+gam(j4,j1,j3)*pg(j4,j2)
                enddo
              enddo
            enddo
          enddo
        enddo
        flag=1
      endif
C
      r1=(ptop(1))**2
      do j1=2,4
        r1=r1 - (ptop(j1))**2
      enddo
c      write(*,*) 'Aint top mass before: ',sqrt(r1)
      call tdeckin(ptop,mt,mb,mw,pf,pfb,pb)
C
c debug
c      r1=(pf(1)+pfb(1))**2
c      do j1=2,4
c        r1=r1 - (pf(j1)+pfb(j1))**2
c      enddo
c      write(*,*) 'Aint dilepton mass: ',sqrt(r1),' mw=',mw
c      r1=(ptop(1))**2
c      do j1=2,4
c        r1=r1 - (ptop(j1))**2
c      enddo
c      write(*,*) 'Aint top mass after: ',sqrt(r1)
c end debug

      call randa(r1)
      if(r1.gt.0.5d0) then
        sp=0.49d0
      else
        sp=-0.49d0
      endif      
C
      call fermionbarsources_d(sp,ubb,pb,mb)
      sp=0.49d0
      call fermionbarsources_d(sp,ubf,pf,0.d0)
      sp=0.51d0
      call fermionsources_d(sp,vfb,pfb,0.d0)
C
      do j1=1,4
        eps(j1)=(0.d0,0.d0)
        do j2=1,4
          tmp=0.d0
          do j3=1,4
            tmp=tmp+gamp(j2,j3,j1)*vfb(j3)
          enddo
          eps(j1)=eps(j1)+ubf(j2)*tmp
        enddo
      enddo
C
      do j1=1,4
        do j2=1,4
          epsl(j1,j2)=(0.d0,0.d0)
          do j3=1,4
            epsl(j1,j2)=epsl(j1,j2)+mtlnz(j3)*eps(j3)*gams(j1,j2,j3)
          enddo
        enddo
      enddo
C
      do j1=1,4
        tmpu(j1)=(0.d0,0.d0)
        do j2=1,4
          tmpu(j1)=tmpu(j1)+ubb(j2)*epsl(j2,j1)
        enddo
      enddo
C
      do j1=1,4
        do j2=1,4
          ptsl(j1,j2)=(0.d0,0.d0)
          do j3=1,4
            ptsl(j1,j2)=ptsl(j1,j2)+mtlnz(j3)*ptop(j3)*gam(j2,j1,j3)
          enddo
          ptsl(j1,j2)=ptsl(j1,j2)+mt*iden(j1,j2)
        enddo
      enddo
C
      do j1=1,4
        utb(j1)=(0.d0,0.d0)
        do j2=1,4
          utb(j1)=utb(j1)+tmpu(j2)*ptsl(j2,j1)
        enddo
      enddo
C
      ew=0.5d0*(mt**2+mw**2-mb**2)/mt
      pw=sqrt(ew**2-mw**2)
      analitic= (128*mt**2*((gAt**2 + gVt**2)*(3.*ew**2 - pw**2)*mt - 
     -      3.*(ew*(gAt**2 + gVt**2) + (-gAt**2 + gVt**2)*mb)*mw**2
     -      ))/6.
      analitic=sqrt(4.d0*mt/analitic)
      do j1=1,4
        utb(j1)=analitic*utb(j1)
      enddo
C
      return
      end
c***************************************************************************
      subroutine tdeckin(ptop,mt,mb,mw,pf,pfb,pb)
c***************************************************************************
C
C  Input  PTOP(4) : top four momentum (E_t, p_x,p_y,p_z)
C         MT,MB,MW  top, bottom and W mass
C  Output PF(4), PFB(4), PB(4) : fermion, antifermion, bottom four-momenta
C                                (En, p_x,p_y,p_z)
C
C
      implicit none
C
      real*8 ptop(4),mt,mb,mw
      real*8 pf(4),pfb(4),pb(4)
C
      real*8 ew,pw,r1,r2,ct,phi,pi,st,qw(4),pfs(4),pbs(4)
      integer j1
      real*8 uw(3),u1(3),u2(3)
C
      pi=acos(-1.d0)
      ew=(mt*mt+mw*mw-mb*mb)/(2.d0*mt)
      pw=sqrt(ew*ew-mw*mw)
C
      call randa(r1)
      call randa(r2)
      ct=2.d0*r1-1.d0
      st=sqrt(1.d0-ct*ct)
      phi=2.d0*pi*r2
      uw(3)=ct
      uw(1)=st*sin(phi)
      uw(2)=st*cos(phi)
      do j1=2,4
       pbs(j1)=pw*uw(j1-1)
      enddo
      pbs(1)=sqrt(mb*mb+pbs(4)*pbs(4)+pbs(3)*pbs(3)+pbs(2)*pbs(2))
      qw(1)=mt-pbs(1)
      do j1=2,4
        qw(j1)=-pbs(j1)
      enddo
C
      if(uw(1).eq.0.d0) then
       u1(1)=1.d0
       u1(2)=0.d0
       u1(3)=0.d0
      else
       ct=1.d0/sqrt(uw(1)**2+uw(2)**2)
       u1(1)=-uw(2)*ct       
       u1(2)=uw(1)*ct
       u1(3)=0
      endif
      u2(1)=uw(2)*u1(3)-uw(3)*u1(2)
      u2(2)=uw(3)*u1(1)-uw(1)*u1(3)
      u2(3)=uw(1)*u1(2)-uw(2)*u1(1)
C
      pw=0.5d0*mw
      call randa(r1)
      call randa(r2)
      ct=2.d0*r1-1.d0
      st=sqrt(1.d0-ct*ct)      
      phi=2.d0*pi*r2
C
      pfs(1)=pw
      do j1=1,3
       pfs(j1+1)=pw*(ct*uw(j1)+st*(cos(phi)*u1(j1)+sin(phi)*u2(j1)))
      enddo
C
      call boost(0,mt*mt,ptop,pbs,pb)
      call boost(0,mw*mw,qw,pfs,pf)
      do j1=1,4
        pfs(j1)=pf(j1)
      enddo
      call boost(0,mt*mt,ptop,pfs,pf)
C
      do j1=1,4
       pfb(j1)=ptop(j1)
      enddo
      do j1=1,4
       pfb(j1)=pfb(j1)-pb(j1)-pf(j1)
      enddo
C
      end
c
C***********************************************************************
          subroutine fermionsources_d(spinsoura,fieldaux ,mom,mass)
C***********************************************************************
C
C This subroutine return the initialized fermion field.
C Elicity heigenstates
C 
          implicit none
C
          real*8  mass         !containing the mass of the external particle
          real*8  mom(4)          !array containing the momenta of the external particle
          real*8 spinsoura    !array containing the spin of the source
          real*8 p,p3p,p3m,mp,coeffp,coeffm
C                               
          complex*16 fieldaux(4)  !array returning the fermion field configuration
          complex*16 im  !immaginary unity in complex representation
          complex*16 p1p
          data im/(0.,1.)/
          integer spinsour
C
C Static variables
C
          save im
C
C  
C
          p=sqrt(mom(2)**2+mom(3)**2+mom(4)**2)
          p3p=p+mom(4)
          p3m=p-mom(4)
          mp=mass+mom(1)
          p1p=mom(2)+mom(3)*im
C
          if(abs(spinsoura).lt.0.5) spinsoura=-spinsoura
          spinsour=nint(100*spinsoura)
          if(abs(p3m).lt.1.d-10.or.abs(p3p).lt.1.d-10)then
           call fs_d(spinsour,fieldaux ,mom,mass)
           return
          endif
C
          coeffp=1./Sqrt(2.*p*mp*p3p)
          coeffm=1./Sqrt(2.*p*mp*p3m)
C
C "spin up" ingoing fermion
C
           if (spinsour.eq.49) then
            fieldaux(1)=coeffp*p3p*mp
            fieldaux(2)=coeffp*p1p*mp
            fieldaux(3)=coeffp*p3p*p
            fieldaux(4)=coeffp*p1p*p
           endif
C
C "spin down" ingoing fermion
C
           if (spinsour.eq.-49) then
            fieldaux(1)=coeffm*p3m*mp
            fieldaux(2)=-coeffm*p1p*mp
            fieldaux(3)=-coeffm*p3m*p
            fieldaux(4)=coeffm*p1p*p
           endif
C
C "spin up" outgoing antifermion
C
           if (spinsour.eq.51) then
            fieldaux(1)=-coeffm*p3m*p
            fieldaux(2)=coeffm*p1p*p
            fieldaux(3)=coeffm*p3m*mp
            fieldaux(4)=-coeffm*p1p*mp
           endif
C
C "spin down" outgoing antifermion
C
           if (spinsour.eq.-51) then
            fieldaux(1)=coeffp*p3p*p 
            fieldaux(2)=coeffp*p1p*p
            fieldaux(3)=coeffp*p3p*mp
            fieldaux(4)=coeffp*p1p*mp
           endif
C
           return
           end

C***********************************************************************
          subroutine fermionbarsources_d(spinsoura,fieldaux ,mom,mass)
C***********************************************************************
C
C This subroutine return the initialized fermionbar field.
C Elicity heigenstates
C 
          implicit none
C
          real*8  mass         !containing the mass of the external particle
          real*8  mom(4)          !array containing the momenta of the external particle
          real*8 spinsoura    !array containing the spin of the source
          real*8 p,p3p,p3m,mp,coeffp,coeffm
C                               
          complex*16 fieldaux(4)     !array returning the fermion field configuration
          complex*16 im  !immaginary unity in complex representation
          complex*16 p1p
          data im/(0.,1.)/
          integer spinsour
C
C Static variables
C
          save im
C
C
          p=sqrt(mom(2)**2+mom(3)**2+mom(4)**2)
          p3p=p+mom(4)
          p3m=p-mom(4)
          mp=mass+mom(1)
          p1p=mom(2)-mom(3)*im
C
          if(abs(spinsoura).lt.0.5) spinsoura=-spinsoura
          spinsour=nint(100*spinsoura)
          if(abs(p3m).lt.1.d-10.or.abs(p3p).lt.1.d-10)then
           call fbs_d(spinsour,fieldaux ,mom,mass)
           return
          endif
C
          coeffp=1./Sqrt(2.*p*mp*p3p)
          coeffm=1./Sqrt(2.*p*mp*p3m)
C
C "spin up" ingoing fermion
C
           if (spinsour.eq.49) then
            fieldaux(1)=coeffp*p3p*mp
            fieldaux(2)=coeffp*p1p*mp
            fieldaux(3)=-coeffp*p3p*p
            fieldaux(4)=-coeffp*p1p*p
           endif
C
C "spin down" ingoing fermion
C
           if (spinsour.eq.-49) then
            fieldaux(1)=coeffm*p3m*mp
            fieldaux(2)=-coeffm*p1p*mp
            fieldaux(3)=coeffm*p3m*p
            fieldaux(4)=-coeffm*p1p*p
           endif
C
C "spin up" outgoing antifermion
C
           if (spinsour.eq.51) then
            fieldaux(1)=-coeffm*p3m*p
            fieldaux(2)=coeffm*p1p*p
            fieldaux(3)=-coeffm*p3m*mp
            fieldaux(4)=coeffm*p1p*mp
           endif
C
C "spin down" outgoing antifermion
C
           if (spinsour.eq.-51) then
            fieldaux(1)=coeffp*p3p*p 
            fieldaux(2)=coeffp*p1p*p
            fieldaux(3)=-coeffp*p3p*mp
            fieldaux(4)=-coeffp*p1p*mp
           endif
C
           return
           end

C***********************************************************************
          subroutine fs_d(spinsour,fieldaux ,mom,mass)
C***********************************************************************
C
C This subroutine return the initialized fermion field.
C 
          implicit none
C
          real*8  mass         !containing the mass of the external particle
          real*8  mom(4)       !array containing the momenta of the external particle
          integer spinsour    !array containing the spin of the source
C                               
          complex*16 fieldaux(4)     !array returning the fermion field configuration
          complex*16 im  !immaginary unity in complex representation
          data im/(0.,1.)/
C
C Static variables
C
          save im
C
C
C "spin up" ingoing fermion
C
           if(mom(4).lt.0) spinsour=-spinsour
           if (spinsour.eq.49) then
            fieldaux(1)=sqrt((mom(1)+mass))
            fieldaux(2)=0. 
            fieldaux(3)=mom(4)/sqrt((mass+mom(1)))         
            fieldaux(4)=mom(2)/sqrt((mass+mom(1)))
     >                         + mom(3)/sqrt((mass+mom(1)))*im          
           endif
C
C "spin down" ingoing fermion
C
           if (spinsour.eq.-49) then
            fieldaux(1)=0.
            fieldaux(2)=sqrt((mom(1)+mass)) 
            fieldaux(3)= mom(2)/sqrt((mass+mom(1)))
     >                         - mom(3)/sqrt((mass+mom(1)))*im         
            fieldaux(4)=-mom(4)/sqrt((mass+mom(1)))
                    endif
C
C "spin up" outgoing antifermion
C
           if (spinsour.eq.51) then
            fieldaux(1)=mom(4)/sqrt((mass+mom(1)))         
            fieldaux(2)=mom(2)/sqrt((mass+mom(1)))
     >                         + mom(3)/sqrt((mass+mom(1)))*im          
            fieldaux(3)=sqrt((mom(1)+mass))
            fieldaux(4)=0. 
           endif
C
C "spin down" outgoing antifermion
C
           if (spinsour.eq.-51) then
            fieldaux(1)= mom(2)/sqrt((mass+mom(1)))
     >                         - mom(3)/sqrt((mass+mom(1)))*im         
            fieldaux(2)=-mom(4)/sqrt((mass+mom(1)))
            fieldaux(3)=0.
            fieldaux(4)=sqrt((mom(1)+mass))          
           endif
C
           return
           end
C***********************************************************************
          subroutine fbs_d(spinsour,fieldaux,mom,mass)
C***********************************************************************
C
C This subroutine return the initialized fermion field.
C 
          implicit none
C
          real*8  mass        !containing the mass of the external particle
          real*8  mom(4)      !array containing the momenta of the external particle
          integer spinsour    !array containing the spin of the source
C                               
          complex*16 fieldaux(4)   !array returning the fermionbar field configuration
          complex*16 im           !immaginary unity in complex representation
          data im/(0.,1.)/
C
C Static variables
C
          save im
C
C
C "spin up" outgoing fermion
C
           if(mom(4).gt.0) spinsour=-spinsour
           if (spinsour.eq.49) then
            fieldaux(1)=sqrt((mom(1)+mass))
            fieldaux(2)=0. 
            fieldaux(3)=-mom(4)/sqrt((mass+mom(1)))         
            fieldaux(4)=-mom(2)/sqrt((mass+mom(1)))
     >                         +mom(3)/sqrt((mass+mom(1)))*im          
           endif
C
C "spin down" outgoing fermion
C
           if (spinsour.eq.-49) then
            fieldaux(1)=0.
            fieldaux(2)=sqrt((mom(1)+mass)) 
            fieldaux(3)= -mom(2)/sqrt((mass+mom(1)))
     >                         - mom(3)/sqrt((mass+mom(1)))*im         
            fieldaux(4)=mom(4)/sqrt((mass+mom(1)))
                    endif
C
C "spin up" ingoing antifermion
C
           if (spinsour.eq.51) then
            fieldaux(1)=-mom(4)/sqrt((mass+mom(1)))         
            fieldaux(2)=-mom(2)/sqrt((mass+mom(1)))
     >                          +mom(3)/sqrt((mass+mom(1)))*im          
            fieldaux(3)=sqrt((mom(1)+mass))
            fieldaux(4)=0. 
           endif
C
C "spin down" ingoing antifermion
C
           if (spinsour.eq.-51) then
            fieldaux(1)= -mom(2)/sqrt((mass+mom(1)))
     >                         - mom(3)/sqrt((mass+mom(1)))*im         
            fieldaux(2)=mom(4)/sqrt((mass+mom(1)))
            fieldaux(3)=0.
            fieldaux(4)=sqrt((mom(1)+mass))          
           endif
C
           return
           end

C*******************************************************************************
      subroutine hvydec(flv)
C*******************************************************************************
C
C Input : FLV  contains the labels of the particles of the given process. 
C         It MUST be equal to ALPHAFL in EVTGEN
C
C Common
C
      implicit none
C
      integer nmax        !maximum number of external particles, change here
      parameter (nmax=10)          
      integer flv(nmax)
C
      real*8 idec(4,3,nmax)
      common/decay/idec     ! idec(1:4,i,j) four momentum of the I-th 
C            decay product of the J-th particle where J is given according
C            to alpha ordering; the four momentum is (E,p_x,p_y,p_z).
C            for the top(topbr) I=1 is b (bbar) I=2 is fermion 
C            I=3 is fermionbar (W-> f fb) 
*
      common/lblt/jtpl,jtmn
*
C
      integer j1,jtpl(2),jtmn(2),ntpl,ntmn    ! at the end of the do loop
C             ntpl contains the number of top particles jtpl(j) the
C             location of the corresponding decay product momenta in
C             the array IDEC (IDEC(*,*,jtpl(1)) are the three four-momenta
C             related to the first top+, ,jtpl(2) to the second and analogously
C             for the top-
C
      jtpl(1)=0
      jtpl(2)=0
      jtmn(1)=0
      jtmn(2)=0
      ntpl=0
      ntmn=0
      do j1=1,nmax
        if(flv(j1).eq.6) then                !top
          ntpl=ntpl+1
          jtpl(ntpl)=j1
        elseif(flv(j1).eq.-6) then           !topbar
          ntmn=ntmn+1
          jtmn(ntmn)=j1
        endif
      enddo
C
      return
      end
C**************************************************************************
      subroutine Vpol(mv,pw,ga,gv,pf,pfb,eps)
C**************************************************************************
C
C Input   MV,PW (vector boson mass and four momentum); GV, GA  vectorial and assial coupling
C Output  PF, PFB  massless fermion/antifermion four momentum (PF = (E,p_x,p_y,p_z))
C         EPS vector boson polarization
C
      implicit none
C
      real*8 mv,gv,ga,pw(4)
      real*8 pf(4), pfb(4)
      complex*16 eps(4) 
C
      complex*16 gam(4,4,5),ubf(4),vfb(4),prj(4,4),iden(4,4)
      real*8 pl,sp,sp0,pfb0(4),pf0(4),ct,phi,jc
      integer j1,j2,j3,j4
C
      data gam/(1.,0.),(0.,0),(0.,0.),(0.,0.),(0.,0.)
     >     ,(1.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     >  (-1.,0.),(0.,0.),(0.,0.),(0.,0.), (0.,0.),(-1.,0.),          !gamma0
     >  (0.,0.),(0.,0.),(0.,0.),(1.,0.),(0.,0.),
     >  (0.,0.),(1.,0.),(0.,0.),(0.,0.),(-1.,0.),
     >  (0.,0.),(0.,0.),(-1.,0.),(0.,0.),(0.,0.),(0.,0.),             !gamma1
     >  (0.,0.),(0.,0.),(0.,0.),(0.,-1.),(0.,0.), 
     >  (0.,0.),(0.,1.),(0.,0.),(0.,0.),(0.,1.), 
     >  (0.,0.),(0.,0.),(0.,-1.),(0.,0.),(0.,0.), (0.,0.),           !gamma2
     >  (0.,0.),(0.,0.),(1.,0.),(0.,0.),(0.,0.),
     >  (0.,0.),(0.,0.),(-1.,0.),(-1.,0.),(0.,0.),
     >  (0.,0.),(0.,0.),(0.,0.),(1.,0.),(0.,0.),(0.,0.),              !gamma3
     >  (0.,0.),(0.,0.),(1.,0.),(0.,0.),(0.,0.), 
     >  (0.,0.),(0.,0.),(1.,0.),(1.,0.),(0.,0.),
     >   (0.,0.),(0.,0.),(0.,0.),(1.,0.),(0.,0.),(0.,0.)/           !gamma5
C
      data iden/(1.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     >          (1.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     >          (1.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(1.,0.)/
               
C
      save gam      
C
      call randa(ct)
      ct=2.d0*ct-1.d0
      call randa(phi)
      phi=2.d0*acos(-1.d0)*phi
C
      pf0(1)=0.5d0*mv
      pf0(4)=pf0(1)*ct
      pf0(2)=pf0(1)*sqrt(1.d0-ct*ct)
      pf0(3)=pf0(2)*sin(phi)
      pf0(2)=pf0(2)*cos(phi)
      pfb0(1)=pf0(1)
      do j1=2,4
        pfb0(j1)=-pf0(j1)
      enddo
C
      call boost(0,mv*mv,pw,pf0,pf)
      call boost(0,mv*mv,pw,pfb0,pfb)
C
c      if(abs(abs(gv)-abs(ga)).gt.1.d-10) then
        call randa(sp)
c        if(sp.lt.0.5d0) then
c          sp0=1.d0
c        else
c          sp0=-1.d0
c        endif
c        jc=2.d0
c      else
c        jc=1.d0
c        sp0=1.d0
c      endif
      pl= (gv-ga)*(gv-ga)/(2.d0*(gv*gv+ga*ga))
      if(sp.lt.pl) then
         sp0= 1.d0
      else 
         sp0= -1.d0
      endif
      sp=0.49d0*sp0
      call fermionbarsources_d(sp,ubf,pf,0.d0)
      sp=0.51d0*sp0
      call fermionsources_d(sp,vfb,pfb,0.d0)
C
      do j1=1,4
        do j2=1,4
          prj(j1,j2)=iden(j1,j2)-sp0*gam(j2,j1,5)
        enddo
      enddo
C
      sp=sqrt(1.d0/(8.d0*mv*mv))
      do j1=1,4
        eps(j1)=(0.d0,0.d0)
        do j2=1,4
          do j3=1,4
            do j4=1,4
              eps(j1)=eps(j1)+ubf(j2)*vfb(j4)*gam(j3,j2,j1)*prj(j3,j4)
            enddo
          enddo
        enddo
        eps(j1)=eps(j1)*sp
      enddo
C
      return  
      end




