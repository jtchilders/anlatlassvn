C********************************************************************
      function Soft(i,j,k,p)
C********************************************************************
c     Calculates the soft limit of the dipol((ij),k) with i beeing
C     the soft parton
C     Automatically created on Wed Feb 16 11:42:30 CET 2005
C     by ./CreateSoft.csh
C---------------------------------------------------------------------------
C
C     This software is not designed or intended for use in on-line control of
C     aircraft, air traffic, aircraft navigation or aircraft communications; 
C     or in the design, construction, operation or maintenance of any nuclear
C     facility. 
C
C---------------------------------------------------------------------------
      implicit none
      double precision Soft
      integer iem,i,j,k
      double precision p(4,-2:4),pn(4,-2:3)
      double precision massq(-2:4)
      double precision C(6,6) ! color correlation matrix
      double precision p1(4),p2(4),p3(4)
      double precision q1(4),q2(4),r1(4),r2(4),sq(4),sqb(4)
      double precision ksoft(4),ki(4),kk(4),nvec(4)
      double precision hein1(4),hein2(4),dtemp
      double precision kq(4),kqb(4),sqrest(4),sqbrest(4)
      double complex AMMM,AMMP,AMPM,AMPP,APMM,APMP,APPM,APPP
      double complex A1MMM,A1MMP,A1MPM,A1MPP,A1PMM,A1PMP,A1PPM,A1PPP
      double complex A2MMM,A2MMP,A2MPM,A2MPP,A2PMM,A2PMP,A2PPM,A2PPP
      double complex A3MMM,A3MMP,A3MPM,A3MPP,A3PMM,A3PMP,A3PPM,A3PPP
      double complex A4MMM,A4MMP,A4MPM,A4MPP,A4PMM,A4PMP,A4PPM,A4PPP
      double complex A5MMM,A5MMP,A5MPM,A5MPP,A5PMM,A5PMP,A5PPM,A5PPP
      double complex A6MMM,A6MMP,A6MPM,A6MPP,A6PMM,A6PMP,A6PPM,A6PPP
      double complex A1MMMc,A1MMPc,A1MPMc,A1MPPc,A1PMMc,A1PMPc,A1PPMc
      double complex A1PPPc,A2MMMc,A2MMPc,A2MPMc,A2MPPc,A2PMMc,A2PMPc
      double complex A2PPMc,A2PPPc,A3MMMc,A3MMPc,A3MPMc,A3MPPc,A3PMMc
      double complex A3PMPc,A3PPMc,A3PPPc,A4MMMc,A4MMPc,A4MPMc,A4MPPc
      double complex A4PMMc,A4PMPc,A4PPMc,A4PPPc,A5MMMc,A5MMPc,A5MPMc
      double complex A5MPPc,A5PMMc,A5PMPc,A5PPMc,A5PPPc,A6MMMc,A6MMPc
      double complex A6MPMc,A6MPPc,A6PMMc,A6PMPc,A6PPMc,A6PPPc     
      double complex dummy
      double precision betaq,betaqb,mqq,miq
      double precision dotp
      integer mu,ispincount,ihelp,err
      data massq/0.d0,0.d0,-1.d0,1.d0,1.d0,0.d0,0.d0/
#include "Soft.dec"
#include "common.F"
#include "constants.F"
      Soft = 0.d0
C     Set the masses:
      do mu =1,4
         nvec(mu) = p(mu,1)
      enddo
      massq(1) = dotp(nvec,nvec)
      massq(2) = massq(1)
C
C     find out the emitter
C 
      if ((j.eq.3).or.(j.eq.4)) then
         iem = 3
      else
         iem = j
      endif
C
C set the n-parton configuration:
C
      do ihelp = -2,2
        do mu = 1,4
          pn(mu,ihelp) = p(mu,ihelp) 
        enddo
      enddo

      if (i.eq.3) then
        do mu = 1,4
          pn(mu,3) = p(mu,4) 
        enddo
      else 
        do mu = 1,4
          pn(mu,3) = p(mu,3) 
        enddo
      endif

      do mu = 1,4
	 ksoft(mu) = p(mu,i)
	 ki(mu) = p(mu,j)
	 kk(mu) = p(mu,k)
      enddo

      do mu = 1,4
        p1(mu) = pn(mu,-1) 
        p2(mu) = pn(mu,-2) 
        p3(mu) = -pn(mu,3) 
	kq(mu) = pn(mu,1) 
	kqb(mu) = pn(mu,2) 
      enddo
C
C set the color matrix, it is bad code but just for testing
C
      C(1,1) = -4711.d0
      if ((iem.eq.-1).and.(k.eq.-2)) 
     >     call cpColor(cg1g2,C)
      if ((iem.eq.-1).and.(k.eq.1)) 
     >     call cpcolor(cg1q,C)
      if ((iem.eq.-1).and.(k.eq.2)) 
     >     call cpcolor(cg1qb,C)
      if ((iem.eq.-1).and.(k.eq.3)) 
     >     call cpcolor(cg1g3,C)
      if ((iem.eq.-1).and.(k.eq.4)) 
     >     call cpcolor(cg1g3,C)
      if ((iem.eq.-2).and.(k.eq.-1)) 
     >     call cpcolor(cg1g2,C)
      if ((iem.eq.-2).and.(k.eq.1)) 
     >     call cpcolor(cg2q,C)
      if ((iem.eq.-2).and.(k.eq.2)) 
     >     call cpcolor(cg2qb,C)
      if ((iem.eq.-2).and.(k.eq.3)) 
     >     call cpcolor(cg2g3,C)
      if ((iem.eq.-2).and.(k.eq.4)) 
     >     call cpcolor(cg2g3,C)
      if ((iem.eq.3).and.(k.eq.-1)) 
     >     call cpcolor(cg1g3,C)
      if ((iem.eq.3).and.(k.eq.-2)) 
     >     call cpcolor(cg2g3,C)
      if ((iem.eq.3).and.(k.eq.1)) 
     >     call cpcolor(cg3q,C)
      if ((iem.eq.3).and.(k.eq.2)) 
     >     call cpcolor(cg3qb,C)
      if ((iem.eq.1).and.(k.eq.-1)) 
     >     call cpcolor(cg1q,C)
      if ((iem.eq.1).and.(k.eq.-2)) 
     >     call cpcolor(cg2q,C)
      if ((iem.eq.1).and.(k.eq.3)) 
     >     call cpcolor(cg3q,C)
      if ((iem.eq.1).and.(k.eq.4)) 
     >     call cpcolor(cg3q,C)
      if ((iem.eq.1).and.(k.eq.2)) 
     >     call cpcolor(cqqb,C)
      if ((iem.eq.2).and.(k.eq.-1)) 
     >     call cpcolor(cg1qb,C)
      if ((iem.eq.2).and.(k.eq.-2)) 
     >     call cpcolor(cg2qb,C)
      if ((iem.eq.2).and.(k.eq.3)) 
     >     call cpcolor(cg3qb,C)
      if ((iem.eq.2).and.(k.eq.4)) 
     >     call cpcolor(cg3qb,C)
      if ((iem.eq.2).and.(k.eq.1)) 
     >     call cpcolor(cqqb,C)
      if (C(1,1).eq.-4711.d0) then
          write(*,*) 'Soft; Colormatrix not assigned'
          stop
      endif
      mqq = dotp(kq,kq)
#ifdef DEBUG
      write(*,*) 'Soft: mq1,mq2 = ',dsqrt(dotp(kq,kq)) 
     >     ,dsqrt(dotp(kqb,kqb)) 
#endif
      betaq = dsqrt( 1.d0 - mqq/kq(4)**2 )
      betaqb = dsqrt( 1.d0 - mqq/kqb(4)**2 )
      sqrest(4) = 0
      sqbrest(4) = 0
      do ihelp=1,3
C     we can take whatever we want here....
         sqrest(ihelp) = kq(ihelp)/betaq/kq(4)
         sqbrest(ihelp) = kqb(ihelp)/betaqb/kqb(4)
      enddo
      call boost(kq,sqrest,sq,err)
      call boost(kqb,sqbrest,sqb,err) 
#ifdef DEBUG
      write(*,*) 'Soft: check the spin vectors for consitency'
      write(*,*) 'Soft: kq.sq = ',dotp(kq,sq) 
      write(*,*) 'Soft: kqb.sqb = ',dotp(kqb,sqb)
      write(*,*) 'Soft: sq.sq = ', dotp(sq,sq)  
      write(*,*) 'Soft: sqb.sqb = ', dotp(sqb,sqb)  
#endif   
      call evalqr(kq,sq,kqb,sqb,q1,q2,r1,r2)
      Soft = 0.d0
#ifdef DEBUG
      write(*,*) 'Soft: some further checks on the kinematics'
      write(*,*) 'Soft: 0 = p1.p1 = ', dotp(p1,p1)
      write(*,*) 'Soft: 0 = p2.p2 = ', dotp(p2,p2)
      write(*,*) 'Soft: 0 = p3.p3 = ', dotp(p3,p3)
      write(*,*) 'Soft: mq^2 = 2*q1.q2 = ', 2.d0 * dotp(q1,q2)
      write(*,*) 'Soft: mq^2 = 2*r1.r2 = ', 2.d0 * dotp(r1,r2)
      do ihelp=1,4
         hein1(ihelp) = p1(ihelp)+p2(ihelp)
	 hein2(ihelp) = -p3(ihelp) + q1(ihelp) + q2(ihelp) 
     .                + r1(ihelp) + r2(ihelp)
         write(*,*) 'Soft:', p1(ihelp)+p2(ihelp)+p3(ihelp), ' = ',
     .        q1(ihelp) + q2(ihelp) + r1(ihelp) + r2(ihelp)
      enddo
      write(*,*) 'Soft: s = ', dotp(hein1,hein1)
      write(*,*) 'Soft: s = ', dotp(hein2,hein2)
#endif 
C     We don't study polarized quarks so we have to sum over all
C     the polarization. The easiest way to do that is buy just
C     caclulating 2 * Soft(q1,q2,r1,r2) + Soft(q1,q2,r2,r1)
C     Note that this works only if there is no term linear in the 
C     spin vector, with other word we assume the structure
C       Soft = (indep.st,stb) +  (something * st.stb) 
C     To add the two configuration we just run two times the following
C     code and exchange at the end r1 <--> r2
      do ispincount=1,2
C
C     There are six different permutations of the partons which we just
C     label 1-6. For every permutation we have in addition 8 different
C     helicity configuration. All together we have thus 48 amplitides
C     to evaluate:
C
      A1MMM = AMMM(p1,p2,p3,q1,q2,r1,r2)
      A1MMMc = conjg(A1MMM)
      A1MMP = AMMP(p1,p2,p3,q1,q2,r1,r2)
      A1MMPc = conjg(A1MMP)
      A1MPM = AMPM(p1,p2,p3,q1,q2,r1,r2)
      A1MPMc = conjg(A1MPM)
      A1MPP = AMPP(p1,p2,p3,q1,q2,r1,r2)
      A1MPPc = conjg(A1MPP)
      A1PMM = APMM(p1,p2,p3,q1,q2,r1,r2)
      A1PMMc = conjg(A1PMM)
      A1PMP = APMP(p1,p2,p3,q1,q2,r1,r2)
      A1PMPc = conjg(A1PMP)
      A1PPM = APPM(p1,p2,p3,q1,q2,r1,r2)
      A1PPMc = conjg(A1PPM)
      A1PPP = APPP(p1,p2,p3,q1,q2,r1,r2)
      A1PPPc = conjg(A1PPP)
      A2MMM = AMMM(p1,p3,p2,q1,q2,r1,r2)
      A2MMMc = conjg(A2MMM)
      A2MMP = AMPM(p1,p3,p2,q1,q2,r1,r2)
      A2MMPc = conjg(A2MMP)
      A2MPM = AMMP(p1,p3,p2,q1,q2,r1,r2)
      A2MPMc = conjg(A2MPM)
      A2MPP = AMPP(p1,p3,p2,q1,q2,r1,r2)
      A2MPPc = conjg(A2MPP)
      A2PMM = APMM(p1,p3,p2,q1,q2,r1,r2)
      A2PMMc = conjg(A2PMM)
      A2PMP = APPM(p1,p3,p2,q1,q2,r1,r2)
      A2PMPc = conjg(A2PMP)
      A2PPM = APMP(p1,p3,p2,q1,q2,r1,r2)
      A2PPMc = conjg(A2PPM)
      A2PPP = APPP(p1,p3,p2,q1,q2,r1,r2)
      A2PPPc = conjg(A2PPP)
      A3MMM = AMMM(p2,p1,p3,q1,q2,r1,r2)
      A3MMMc = conjg(A3MMM)
      A3MMP = AMMP(p2,p1,p3,q1,q2,r1,r2)
      A3MMPc = conjg(A3MMP)
      A3MPM = APMM(p2,p1,p3,q1,q2,r1,r2)
      A3MPMc = conjg(A3MPM)
      A3MPP = APMP(p2,p1,p3,q1,q2,r1,r2)
      A3MPPc = conjg(A3MPP)
      A3PMM = AMPM(p2,p1,p3,q1,q2,r1,r2)
      A3PMMc = conjg(A3PMM)
      A3PMP = AMPP(p2,p1,p3,q1,q2,r1,r2)
      A3PMPc = conjg(A3PMP)
      A3PPM = APPM(p2,p1,p3,q1,q2,r1,r2)
      A3PPMc = conjg(A3PPM)
      A3PPP = APPP(p2,p1,p3,q1,q2,r1,r2)
      A3PPPc = conjg(A3PPP)
      A4MMM = AMMM(p2,p3,p1,q1,q2,r1,r2)
      A4MMMc = conjg(A4MMM)
      A4MMP = AMPM(p2,p3,p1,q1,q2,r1,r2)
      A4MMPc = conjg(A4MMP)
      A4MPM = APMM(p2,p3,p1,q1,q2,r1,r2)
      A4MPMc = conjg(A4MPM)
      A4MPP = APPM(p2,p3,p1,q1,q2,r1,r2)
      A4MPPc = conjg(A4MPP)
      A4PMM = AMMP(p2,p3,p1,q1,q2,r1,r2)
      A4PMMc = conjg(A4PMM)
      A4PMP = AMPP(p2,p3,p1,q1,q2,r1,r2)
      A4PMPc = conjg(A4PMP)
      A4PPM = APMP(p2,p3,p1,q1,q2,r1,r2)
      A4PPMc = conjg(A4PPM)
      A4PPP = APPP(p2,p3,p1,q1,q2,r1,r2)
      A4PPPc = conjg(A4PPP)
      A5MMM = AMMM(p3,p1,p2,q1,q2,r1,r2)
      A5MMMc = conjg(A5MMM)
      A5MMP = APMM(p3,p1,p2,q1,q2,r1,r2)
      A5MMPc = conjg(A5MMP)
      A5MPM = AMMP(p3,p1,p2,q1,q2,r1,r2)
      A5MPMc = conjg(A5MPM)
      A5MPP = APMP(p3,p1,p2,q1,q2,r1,r2)
      A5MPPc = conjg(A5MPP)
      A5PMM = AMPM(p3,p1,p2,q1,q2,r1,r2)
      A5PMMc = conjg(A5PMM)
      A5PMP = APPM(p3,p1,p2,q1,q2,r1,r2)
      A5PMPc = conjg(A5PMP)
      A5PPM = AMPP(p3,p1,p2,q1,q2,r1,r2)
      A5PPMc = conjg(A5PPM)
      A5PPP = APPP(p3,p1,p2,q1,q2,r1,r2)
      A5PPPc = conjg(A5PPP)
      A6MMM = AMMM(p3,p2,p1,q1,q2,r1,r2)
      A6MMMc = conjg(A6MMM)
      A6MMP = APMM(p3,p2,p1,q1,q2,r1,r2)
      A6MMPc = conjg(A6MMP)
      A6MPM = AMPM(p3,p2,p1,q1,q2,r1,r2)
      A6MPMc = conjg(A6MPM)
      A6MPP = APPM(p3,p2,p1,q1,q2,r1,r2)
      A6MPPc = conjg(A6MPP)
      A6PMM = AMMP(p3,p2,p1,q1,q2,r1,r2)
      A6PMMc = conjg(A6PMM)
      A6PMP = APMP(p3,p2,p1,q1,q2,r1,r2)
      A6PMPc = conjg(A6PMP)
      A6PPM = AMPP(p3,p2,p1,q1,q2,r1,r2)
      A6PPMc = conjg(A6PPM)
      A6PPP = APPP(p3,p2,p1,q1,q2,r1,r2)
      A6PPPc = conjg(A6PPP)
      t1 = C(1,3)
      t4 = C(1,5)
      t7 = C(1,4)
      t10 = C(2,3)
      t13 = C(2,4)
      t16 = C(2,5)
      t21 = C(2,2)
      t24 = C(5,6)
      t29 = C(3,5)
      t32 = C(3,6)
      t47 = t1*A1PMM*A3PMMc+t4*A1MMP*A5MMPc+t7*A4PMP*A1PMPc+t10*A2MPP*A3
     #MPPc+t13*A2PMM*A4PMMc+t16*A2MPM*A5MPMc+t1*A3PPM*A1PPMc+t21*A2MMM*A
     #2MMMc+t24*A6PMM*A5PMMc+t1*A1PPM*A3PPMc+t29*A5MMM*A3MMMc+t32*A3PPM*
     #A6PPMc+t7*A4PPM*A1PPMc+t32*A3PMP*A6PMPc+t29*A5PMM*A3PMMc+t21*A2MMP
     #*A2MMPc+t4*A5PPP*A1PPPc+t29*A5MMP*A3MMPc
      t52 = C(3,3)
      t57 = C(1,6)
      t60 = C(4,5)
      t63 = C(4,6)
      t66 = C(2,6)
      t75 = C(4,4)
      t82 = C(5,5)
      t87 = C(6,6)
      t92 = t10*A2MMP*A3MMPc+t1*A1MPM*A3MPMc+t52*A3MMP*A3MMPc+t32*A3MPP*
     #A6MPPc+t57*A1MPP*A6MPPc+t60*A5PPM*A4PPMc+t63*A4PMM*A6PMMc+t66*A2PM
     #M*A6PMMc+t13*A4PMM*A2PMMc+t24*A5PMM*A6PMMc+t24*A6MPP*A5MPPc+t75*A4
     #PMP*A4PMPc+t60*A4PPP*A5PPPc+t29*A3PPM*A5PPMc+t82*A5PMM*A5PMMc+t4*A
     #1PPP*A5PPPc+t87*A6MMP*A6MMPc+t10*A2MPM*A3MPMc
      t94 = C(1,2)
      t97 = C(3,4)
      t132 = t94*A1PMP*A2PMPc+t97*A4PPM*A3PPMc+t66*A6PPP*A2PPPc+t63*A4PM
     #P*A6PMPc+t16*A5PMM*A2PMMc+t75*A4PMM*A4PMMc+t24*A5PPP*A6PPPc+t66*A6
     #MPP*A2MPPc+t82*A5PPM*A5PPMc+t24*A6PPP*A5PPPc+t24*A5MMP*A6MMPc+t29*
     #A5MPM*A3MPMc+t52*A3PMP*A3PMPc+t32*A3PPP*A6PPPc+t29*A3PPP*A5PPPc+t2
     #1*A2MPP*A2MPPc+t60*A5MPM*A4MPMc+t13*A2MMM*A4MMMc
      t169 = t66*A6PMP*A2PMPc+t94*A2MPP*A1MPPc+t97*A3PPM*A4PPMc+t7*A4MMP
     #*A1MMPc+t1*A3PPP*A1PPPc+t94*A1MPP*A2MPPc+t24*A5MPM*A6MPMc+t57*A6MP
     #P*A1MPPc+t63*A4PPM*A6PPMc+t10*A2PPP*A3PPPc+t87*A6MPP*A6MPPc+t10*A3
     #PMP*A2PMPc+t7*A4MPP*A1MPPc+t66*A6PPM*A2PPMc+t60*A5MPP*A4MPPc+t29*A
     #5MPP*A3MPPc+t94*A1PMM*A2PMMc+t97*A4MMP*A3MMPc
      t172 = C(1,1)
      t209 = t172*A1PMM*A1PMMc+t87*A6PPP*A6PPPc+t24*A6PPM*A5PPMc+t94*A1M
     #MP*A2MMPc+t1*A3MMP*A1MMPc+t16*A2PMP*A5PMPc+t10*A2PPM*A3PPMc+t82*A5
     #PMP*A5PMPc+t97*A3MMP*A4MMPc+t97*A4PMM*A3PMMc+t172*A1PMP*A1PMPc+t66
     #*A6MPM*A2MPMc+t7*A4MMM*A1MMMc+t16*A5MPP*A2MPPc+t16*A2PMM*A5PMMc+t2
     #1*A2MPM*A2MPMc+t60*A5PMP*A4PMPc+t172*A1PPM*A1PPMc
      t246 = t4*A1MMM*A5MMMc+t172*A1PPP*A1PPPc+t82*A5MPM*A5MPMc+t97*A3PM
     #P*A4PMPc+t75*A4PPM*A4PPMc+t13*A2PPM*A4PPMc+t75*A4MMP*A4MMPc+t32*A3
     #MMP*A6MMPc+t66*A2MMP*A6MMPc+t7*A4PPP*A1PPPc+t16*A2MMP*A5MMPc+t7*A1
     #MPP*A4MPPc+t29*A5PPP*A3PPPc+t16*A2MPP*A5MPPc+t57*A1MMP*A6MMPc+t94*
     #A1MMM*A2MMMc+t60*A5MMP*A4MMPc+t172*A1MMM*A1MMMc
      t284 = t13*A2MMP*A4MMPc+t10*A3MMM*A2MMMc+t1*A1PMP*A3PMPc+t75*A4MPM
     #*A4MPMc+t57*A6PPM*A1PPMc+t52*A3PPP*A3PPPc+t16*A2PPP*A5PPPc+t57*A1P
     #PP*A6PPPc+t32*A6PMP*A3PMPc+t32*A3MMM*A6MMMc+t32*A6MPP*A3MPPc+t57*A
     #6PMP*A1PMPc+t75*A4PPP*A4PPPc+t32*A3PMM*A6PMMc+t57*A6MPM*A1MPMc+t82
     #*A5MMM*A5MMMc+t4*A5MMP*A1MMPc+t24*A6PMP*A5PMPc
      t321 = t16*A5MMP*A2MMPc+t66*A6PMM*A2PMMc+t24*A6MPM*A5MPMc+t13*A4PM
     #P*A2PMPc+t10*A2MMM*A3MMMc+t13*A4MPM*A2MPMc+t32*A6PPM*A3PPMc+t4*A1M
     #PP*A5MPPc+t94*A2MPM*A1MPMc+t16*A2PPM*A5PPMc+t94*A1MPM*A2MPMc+t10*A
     #3MMP*A2MMPc+t94*A2MMP*A1MMPc+t1*A3PMM*A1PMMc+t13*A4MMP*A2MMPc+t63*
     #A6PPP*A4PPPc+t21*A2PMM*A2PMMc+t97*A3PMM*A4PMMc
      t361 = t60*A4PPM*A5PPMc+t13*A2PPP*A4PPPc+t60*A4MMP*A5MMPc+t63*A4MM
     #P*A6MMPc+t57*A6MMM*A1MMMc+t66*A6MMM*A2MMMc+t32*A6MMM*A3MMMc+t63*A6
     #PPM*A4PPMc+t63*A6MMM*A4MMMc+t7*A1PMP*A4PMPc+t94*A2MMM*A1MMMc+t57*A
     #1PPM*A6PPMc+t87*A6MMM*A6MMMc+t66*A2PMP*A6PMPc+t4*A1PMM*A5PMMc+t172
     #*A1MMP*A1MMPc+t1*A1MMP*A3MMPc+t29*A5PPM*A3PPMc
      t398 = t16*A5PPM*A2PPMc+t4*A5PPM*A1PPMc+t13*A4PPM*A2PPMc+t52*A3PPM
     #*A3PPMc+t10*A3PPM*A2PPMc+t66*A2PPM*A6PPMc+t57*A1PMP*A6PMPc+t29*A3P
     #MM*A5PMMc+t29*A3MMP*A5MMPc+t4*A1PPM*A5PPMc+t13*A4PPP*A2PPPc+t10*A3
     #PPP*A2PPPc+t4*A5PMM*A1PMMc+t32*A6PPP*A3PPPc+t24*A5MMM*A6MMMc+t60*A
     #5MMM*A4MMMc+t16*A5MMM*A2MMMc+t4*A5MMM*A1MMMc
      t436 = t63*A4MMM*A6MMMc+t60*A4MMM*A5MMMc+t97*A4MMM*A3MMMc+t13*A4MM
     #M*A2MMMc+t29*A3MMM*A5MMMc+t97*A3MMM*A4MMMc+t52*A3MMM*A3MMMc+t1*A3M
     #MM*A1MMMc+t66*A2MMM*A6MMMc+t16*A2MMM*A5MMMc+t24*A6MMM*A5MMMc+t57*A
     #1MMM*A6MMMc+t7*A1MMM*A4MMMc+t1*A1MMM*A3MMMc+t7*A1MMP*A4MMPc+t94*A2
     #PMM*A1PMMc+t57*A1PMM*A6PMMc+t87*A6PMM*A6PMMc
      t473 = t63*A6MPP*A4MPPc+t24*A5MPP*A6MPPc+t82*A5MPP*A5MPPc+t4*A5MPP
     #*A1MPPc+t63*A4MPP*A6MPPc+t75*A4MMM*A4MMMc+t60*A4MPP*A5MPPc+t75*A4M
     #PP*A4MPPc+t97*A4MPP*A3MPPc+t13*A4MPP*A2MPPc+t29*A3MPP*A5MPPc+t97*A
     #3MPP*A4MPPc+t52*A3MPP*A3MPPc+t10*A3MPP*A2MPPc+t1*A3MPP*A1MPPc+t66*
     #A2MPP*A6MPPc+t13*A2MPP*A4MPPc+t1*A1MPP*A3MPPc
      t512 = t172*A1MPP*A1MPPc+t87*A6MPM*A6MPMc+t63*A6MPM*A4MPMc+t32*A6M
     #PM*A3MPMc+t16*A5MPM*A2MPMc+t4*A5MPM*A1MPMc+t63*A4MPM*A6MPMc+t60*A4
     #MPM*A5MPMc+t97*A4MPM*A3MPMc+t7*A4MPM*A1MPMc+t32*A3MPM*A6MPMc+t29*A
     #3MPM*A5MPMc+t97*A3MPM*A4MPMc+t52*A3MPM*A3MPMc+t10*A3MPM*A2MPMc+t1*
     #A3MPM*A1MPMc+t66*A2MPM*A6MPMc+t13*A2MPM*A4MPMc
      t549 = t57*A1MPM*A6MPMc+t4*A1MPM*A5MPMc+t7*A1MPM*A4MPMc+t172*A1MPM
     #*A1MPMc+t63*A6MMP*A4MMPc+t32*A6MMP*A3MMPc+t66*A6MMP*A2MMPc+t57*A6M
     #MP*A1MMPc+t82*A5MMP*A5MMPc+t21*A2PPM*A2PPMc+t24*A5PMP*A6PMPc+t94*A
     #2PPM*A1PPMc+t7*A1PPM*A4PPMc+t94*A1PPM*A2PPMc+t87*A6PMP*A6PMPc+t63*
     #A6PMP*A4PMPc+t29*A5PMP*A3PMPc+t16*A5PMP*A2PMPc
      t587 = t4*A5PMP*A1PMPc+t60*A4PMP*A5PMPc+t97*A4PMP*A3PMPc+t29*A3PMP
     #*A5PMPc+t1*A3PMP*A1PMPc+t13*A2PMP*A4PMPc+t10*A2PMP*A3PMPc+t21*A2PM
     #P*A2PMPc+t24*A6MMP*A5MMPc+t94*A2PMP*A1PMPc+t4*A1PMP*A5PMPc+t63*A6P
     #MM*A4PMMc+t32*A6PMM*A3PMMc+t57*A6PMM*A1PMMc+t60*A5PMM*A4PMMc+t60*A
     #4PMM*A5PMMc+t7*A4PMM*A1PMMc+t52*A3PMM*A3PMMc
      t624 = t10*A3PMM*A2PMMc+t10*A2PMM*A3PMMc+t57*A6PPP*A1PPPc+t60*A5PP
     #P*A4PPPc+t82*A5PPP*A5PPPc+t16*A5PPP*A2PPPc+t7*A1PMM*A4PMMc+t63*A4P
     #PP*A6PPPc+t97*A4PPP*A3PPPc+t97*A3PPP*A4PPPc+t66*A2PPP*A6PPPc+t21*A
     #2PPP*A2PPPc+t94*A2PPP*A1PPPc+t7*A1PPP*A4PPPc+t1*A1PPP*A3PPPc+t94*A
     #1PPP*A2PPPc+t87*A6PPM*A6PPMc+t24*A5PPM*A6PPMc
      dummy = t47+t92+t132+t169+t209+t246+t284+t321+t361+t398+t436+t473+
     #t512+t549+t587+t624
      Soft = Soft + dummy
C
C     exchange the vectors r1 <-> r2 it cooresponds to sqb --> -sqb
C
        do ihelp = 1,4
          dtemp = r1(ihelp)
          r1(ihelp) = r2(ihelp)
          r2(ihelp) = dtemp
        enddo
      enddo ! ispincount
      miq = massq(j)
C     multiply by the overall singular factor:      
      Soft = -8.d0 * pi * alphas * 
     >    Soft / dotp(ki,ksoft)
     >    * (dotp(ki,kk)/(dotp(ki,ksoft)+dotp(kk,ksoft))
     >    - miq/2.d0/dotp(ki,ksoft))
*      write(*,*) 'Soft: ', soft
      Soft = Soft * g2**3
      end
