C********************************************************************
      function Dqgq(iem,pn,C,V)
C********************************************************************
C     Automatically created on Wed Feb 16 11:42:23 CET 2005
C     by ./CreateDqgq.csh
C---------------------------------------------------------------------------
C
C     This software is not designed or intended for use in on-line control of
C     aircraft, air traffic, aircraft navigation or aircraft communications; 
C     or in the design, construction, operation or maintenance of any nuclear
C     facility. 
C
C---------------------------------------------------------------------------
      implicit none
      double precision Dqgq,V 
      integer iem
      double precision pn(4,-2:3)
      double precision p1(4),p2(4),p3(4)
      double precision q1(4),q2(4),r1(4),r2(4),sq(4),sqb(4)
      double precision kq(4),kqb(4),sqrest(4),sqbrest(4),pem(4)
      double precision hein1(4),hein2(4)
      double precision C(6,6) ! color correlation matrix
      double precision dtemp
      double complex smmm,smmp,smpm,smpp,spmm,spmp,sppm,sppp
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
      double precision betaq,betaqb,mqq
      double precision dotp
      integer mu,j,ispincount,error
#include "Dqgq.dec"
#include "common.F"
#include "constants.F"
      Dqgq = 0.d0
#ifdef DEBUG
      CALL PRINTCONF(PN)
#endif
      do mu = 1,4
        p1(mu) = pn(mu,-1) 
        p2(mu) = pn(mu,-2) 
        p3(mu) = -pn(mu,3) 
	kq(mu) = pn(mu,1) 
	kqb(mu) = pn(mu,2) 
	pem(mu) = pn(mu,iem) 
      enddo

      mqq = dotp(kq,kq)
#ifdef DEBUG
      write(*,*) 'Dqgq: mq1,mq2 = ',dsqrt(dotp(kq,kq)) 
     >     ,dsqrt(dotp(kqb,kqb)) 
#endif
      betaq = dsqrt( 1.d0 - mqq/kq(4)**2 )
      betaqb = dsqrt( 1.d0 - mqq/kqb(4)**2 )
      sqrest(4) = 0
      sqbrest(4) = 0
      do j=1,3
C     we can take whatever we want here....
         sqrest(j) = kq(j)/betaq/kq(4)
         sqbrest(j) = kqb(j)/betaqb/kqb(4)
      enddo
      call boost(kq,sqrest,sq,error) 
      if (error.lt.0) then
         ierr=error
         Dqgq=0.d0
         return
      endif
      call boost(kqb,sqbrest,sqb,error) 
      if (error.lt.0) then
         ierr=error
         Dqgq=0.d0
	 return
      endif
#ifdef DEBUG
      write(*,*) 'Dqgq: check the spin vectors for consitency'
      write(*,*) 'Dqgq: kq.sq = ',dotp(kq,sq) 
      write(*,*) 'Dqgq: kqb.sqb = ',dotp(kqb,sqb)
      write(*,*) 'Dqgq: sq.sq = ', dotp(sq,sq)  
      write(*,*) 'Dqgq: sqb.sqb = ', dotp(sqb,sqb)  
#endif   
      call evalqr(kq,sq,kqb,sqb,q1,q2,r1,r2)
      Dqgq = 0.d0
#ifdef DEBUG
      write(*,*) 'Dqgq: some further checks on the kinematics'
      write(*,*) 'Dqgq: 0 = p1.p1 = ', dotp(p1,p1)
      write(*,*) 'Dqgq: 0 = p2.p2 = ', dotp(p2,p2)
      write(*,*) 'Dqgq: 0 = p3.p3 = ', dotp(p3,p3)
      write(*,*) 'Dqgq: mq^2 = 2*q1.q2 = ', 2.d0 * dotp(q1,q2)
      write(*,*) 'Dqgq: mq^2 = 2*r1.r2 = ', 2.d0 * dotp(r1,r2)
      do j=1,4
         hein1(j) = p1(j)+p2(j)
	 hein2(j) = -p3(j) + q1(j) + q2(j) + r1(j) + r2(j)
         write(*,*) 'Dqgq:', p1(j)+p2(j)+p3(j), ' = ',
     .        q1(j) + q2(j) + r1(j) + r2(j)
      enddo
      write(*,*) 'Dqgq: s = ', dotp(hein1,hein1)
      write(*,*) 'Dqgq: s = ', dotp(hein2,hein2)
#endif 
C     We don't study polarized quarks so we have to sum over all
C     the polarization. The easiest way to do that is buy just
C     caclulating 2 * Dqgq(q1,q2,r1,r2) + Dqgq(q1,q2,r2,r1)
C     Note that this works only if there is no term linear in the 
C     spin vector, with other word we assume the structure
C       Dqgq = (indep.st,stb) +  (something * st.stb) 
C     To add the two configuration we just run two times the following
C     code and exchange at the end r1 <--> r2
      do ispincount=1,2
C
C     There are six different permutations of the partons which we just
C     label 1-6. For every permutation we have in addition 8 different
C     helicity configuration. All together we have thus 48 amplitudes
C     to evaluate:
C
      call gggtthel(p1,p2,p3,q1,q2,r1,r2,
     . smmm,smmp,smpm,smpp,spmm,spmp,sppm,sppp)
      A1MMM = sMMM
      A1MMMc = conjg(A1MMM)
      A1MMP = sMMP
      A1MMPc = conjg(A1MMP)
      A1MPM = sMPM
      A1MPMc = conjg(A1MPM)
      A1MPP = sMPP
      A1MPPc = conjg(A1MPP)
      A1PMM = sPMM
      A1PMMc = conjg(A1PMM)
      A1PMP = sPMP
      A1PMPc = conjg(A1PMP)
      A1PPM = sPPM
      A1PPMc = conjg(A1PPM)
      A1PPP = sPPP
      A1PPPc = conjg(A1PPP)
      call gggtthel(p1,p3,p2,q1,q2,r1,r2,
     . smmm,smmp,smpm,smpp,spmm,spmp,sppm,sppp)
      A2MMM = sMMM
      A2MMMc = conjg(A2MMM)
      A2MMP = sMPM
      A2MMPc = conjg(A2MMP)
      A2MPM = sMMP
      A2MPMc = conjg(A2MPM)
      A2MPP = sMPP
      A2MPPc = conjg(A2MPP)
      A2PMM = sPMM
      A2PMMc = conjg(A2PMM)
      A2PMP = sPPM
      A2PMPc = conjg(A2PMP)
      A2PPM = sPMP
      A2PPMc = conjg(A2PPM)
      A2PPP = sPPP
      A2PPPc = conjg(A2PPP)
      call gggtthel(p2,p1,p3,q1,q2,r1,r2,
     . smmm,smmp,smpm,smpp,spmm,spmp,sppm,sppp)
      A3MMM = sMMM
      A3MMMc = conjg(A3MMM)
      A3MMP = sMMP
      A3MMPc = conjg(A3MMP)
      A3MPM = sPMM
      A3MPMc = conjg(A3MPM)
      A3MPP = sPMP
      A3MPPc = conjg(A3MPP)
      A3PMM = sMPM
      A3PMMc = conjg(A3PMM)
      A3PMP = sMPP
      A3PMPc = conjg(A3PMP)
      A3PPM = sPPM
      A3PPMc = conjg(A3PPM)
      A3PPP = sPPP
      A3PPPc = conjg(A3PPP)
      call gggtthel(p2,p3,p1,q1,q2,r1,r2,
     . smmm,smmp,smpm,smpp,spmm,spmp,sppm,sppp)
      A4MMM = sMMM
      A4MMMc = conjg(A4MMM)
      A4MMP = sMPM
      A4MMPc = conjg(A4MMP)
      A4MPM = sPMM
      A4MPMc = conjg(A4MPM)
      A4MPP = sPPM
      A4MPPc = conjg(A4MPP)
      A4PMM = sMMP
      A4PMMc = conjg(A4PMM)
      A4PMP = sMPP
      A4PMPc = conjg(A4PMP)
      A4PPM = sPMP
      A4PPMc = conjg(A4PPM)
      A4PPP = sPPP
      A4PPPc = conjg(A4PPP)
      call gggtthel(p3,p1,p2,q1,q2,r1,r2,
     . smmm,smmp,smpm,smpp,spmm,spmp,sppm,sppp)
      A5MMM = sMMM
      A5MMMc = conjg(A5MMM)
      A5MMP = sPMM
      A5MMPc = conjg(A5MMP)
      A5MPM = sMMP
      A5MPMc = conjg(A5MPM)
      A5MPP = sPMP
      A5MPPc = conjg(A5MPP)
      A5PMM = sMPM
      A5PMMc = conjg(A5PMM)
      A5PMP = sPPM
      A5PMPc = conjg(A5PMP)
      A5PPM = sMPP
      A5PPMc = conjg(A5PPM)
      A5PPP = sPPP
      A5PPPc = conjg(A5PPP)
      call gggtthel(p3,p2,p1,q1,q2,r1,r2,
     . smmm,smmp,smpm,smpp,spmm,spmp,sppm,sppp)
      A6MMM = sMMM
      A6MMMc = conjg(A6MMM)
      A6MMP = sPMM
      A6MMPc = conjg(A6MMP)
      A6MPM = sMPM
      A6MPMc = conjg(A6MPM)
      A6MPP = sPPM
      A6MPPc = conjg(A6MPP)
      A6PMM = sMMP
      A6PMMc = conjg(A6PMM)
      A6PMP = sPMP
      A6PMPc = conjg(A6PMP)
      A6PPM = sMPP
      A6PPMc = conjg(A6PPM)
      A6PPP = sPPP
      A6PPPc = conjg(A6PPP)
      t1 = C(1,4)
      t4 = C(3,4)
      t7 = C(1,3)
      t10 = C(2,3)
      t15 = C(2,4)
      t20 = C(1,5)
      t27 = C(4,6)
      t34 = C(1,6)
      t41 = C(6,6)
      t44 = C(2,5)
      t47 = t1*A1MPP*A4MPPc+t4*A4PPM*A3PPMc+t7*A1MMM*A3MMMc+t10*A2PPP*A3
     #PPPc+t10*A2PMP*A3PMPc+t15*A2PMM*A4PMMc+t7*A1MMP*A3MMPc+t20*A1PMP*A
     #5PMPc+t4*A4MMM*A3MMMc+t15*A4PPP*A2PPPc+t27*A4PMM*A6PMMc+t27*A4MMM*
     #A6MMMc+t4*A3MPP*A4MPPc+t34*A1MMM*A6MMMc+t7*A1MPM*A3MPMc+t4*A4PMM*A
     #3PMMc+t41*A6MMP*A6MMPc+t44*A5MPP*A2MPPc
      t48 = C(1,1)
      t55 = C(3,5)
      t64 = C(5,6)
      t67 = C(2,6)
      t74 = C(4,4)
      t79 = C(2,2)
      t90 = t48*A1MPP*A1MPPc+t15*A2MPP*A4MPPc+t41*A6MPM*A6MPMc+t55*A5MPP
     #*A3MPPc+t41*A6PPM*A6PPMc+t27*A6PPP*A4PPPc+t10*A3MPM*A2MPMc+t64*A6P
     #PP*A5PPPc+t67*A6MMM*A2MMMc+t20*A5MMM*A1MMMc+t41*A6PPP*A6PPPc+t74*A
     #4MPM*A4MPMc+t34*A6MPP*A1MPPc+t79*A2PMP*A2PMPc+t10*A3MMP*A2MMPc+t64
     #*A5PPP*A6PPPc+t48*A1PMM*A1PMMc+t20*A5PPM*A1PPMc
      t98 = C(4,5)
      t119 = C(1,2)
      t124 = C(3,3)
      t131 = t1*A4PPM*A1PPMc+t7*A3MPP*A1MPPc+t55*A3PMM*A5PMMc+t98*A4PPM*
     #A5PPMc+t20*A1MMM*A5MMMc+t34*A1MPP*A6MPPc+t98*A4MPP*A5MPPc+t44*A2MM
     #P*A5MMPc+t55*A5MPM*A3MPMc+t79*A2PPM*A2PPMc+t4*A4PMP*A3PMPc+t15*A4M
     #MP*A2MMPc+t27*A4MPM*A6MPMc+t119*A2MPP*A1MPPc+t74*A4PMP*A4PMPc+t124
     #*A3PMM*A3PMMc+t10*A3MPP*A2MPPc+t27*A6MPP*A4MPPc
      t138 = C(3,6)
      t169 = t55*A5MMM*A3MMMc+t27*A6MPM*A4MPMc+t15*A2MMM*A4MMMc+t138*A3P
     #PP*A6PPPc+t119*A2PPM*A1PPMc+t10*A2MMM*A3MMMc+t98*A5PPM*A4PPMc+t98*
     #A5MMM*A4MMMc+t64*A5PPM*A6PPMc+t4*A4MMP*A3MMPc+t119*A2MPM*A1MPMc+t2
     #7*A6MMM*A4MMMc+t55*A5PMM*A3PMMc+t20*A5PMP*A1PMPc+t44*A5PMM*A2PMMc+
     #t138*A3MMM*A6MMMc+t4*A4MPP*A3MPPc+t1*A4MMM*A1MMMc
      t208 = t138*A3MMP*A6MMPc+t98*A4MPM*A5MPMc+t64*A6MPP*A5MPPc+t98*A5P
     #MP*A4PMPc+t34*A6MMP*A1MMPc+t64*A5MPM*A6MPMc+t44*A2MPP*A5MPPc+t44*A
     #5MPM*A2MPMc+t74*A4MPP*A4MPPc+t138*A6MMP*A3MMPc+t55*A5PPP*A3PPPc+t9
     #8*A5PPP*A4PPPc+t79*A2MMP*A2MMPc+t7*A3MMM*A1MMMc+t98*A4MMM*A5MMMc+t
     #64*A6MMP*A5MMPc+t4*A4MPM*A3MPMc+t119*A1MMM*A2MMMc
      t223 = C(5,5)
      t246 = t34*A6MMM*A1MMMc+t1*A4MMP*A1MMPc+t138*A6PMP*A3PMPc+t10*A2MP
     #M*A3MPMc+t41*A6MMM*A6MMMc+t4*A3PPM*A4PPMc+t27*A4MMP*A6MMPc+t223*A5
     #PPP*A5PPPc+t119*A1MPP*A2MPPc+t124*A3PPP*A3PPPc+t34*A6PPP*A1PPPc+t1
     #19*A2MMP*A1MMPc+t10*A3PMM*A2PMMc+t64*A5MMP*A6MMPc+t67*A6PPP*A2PPPc
     #+t41*A6MPP*A6MPPc+t44*A2MMM*A5MMMc+t138*A6PPP*A3PPPc
      t284 = t15*A4PMP*A2PMPc+t1*A1MMM*A4MMMc+t67*A6MPP*A2MPPc+t34*A6PMM
     #*A1PMMc+t27*A6MMP*A4MMPc+t7*A3PPM*A1PPMc+t55*A3MMP*A5MMPc+t34*A6MP
     #M*A1MPMc+t119*A1PPP*A2PPPc+t15*A4MPP*A2MPPc+t10*A3MMM*A2MMMc+t20*A
     #5MPM*A1MPMc+t34*A1MMP*A6MMPc+t7*A1PPP*A3PPPc+t79*A2MPM*A2MPMc+t64*
     #A6PPM*A5PPMc+t15*A2PPM*A4PPMc+t67*A2MPM*A6MPMc
      t321 = t119*A2PPP*A1PPPc+t55*A5PPM*A3PPMc+t7*A1MPP*A3MPPc+t44*A5PP
     #M*A2PPMc+t138*A6MMM*A3MMMc+t64*A5MMM*A6MMMc+t64*A6MPM*A5MPMc+t10*A
     #2PMM*A3PMMc+t67*A2MMM*A6MMMc+t20*A1PMM*A5PMMc+t7*A1PMM*A3PMMc+t1*A
     #4PMP*A1PMPc+t15*A4MPM*A2MPMc+t34*A1MPM*A6MPMc+t27*A4PMP*A6PMPc+t1*
     #A4PPP*A1PPPc+t98*A4PMM*A5PMMc+t1*A1PPP*A4PPPc
      t361 = t55*A3MPP*A5MPPc+t138*A3MPP*A6MPPc+t55*A5PMP*A3PMPc+t20*A1M
     #PM*A5MPMc+t44*A2MPM*A5MPMc+t44*A5MMM*A2MMMc+t223*A5MPM*A5MPMc+t98*
     #A4PMP*A5PMPc+t15*A4PMM*A2PMMc+t67*A6PMP*A2PMPc+t10*A2MPP*A3MPPc+t2
     #7*A4MPP*A6MPPc+t4*A3PMM*A4PMMc+t20*A1PPP*A5PPPc+t67*A2MMP*A6MMPc+t
     #20*A5PMM*A1PMMc+t138*A6MPM*A3MPMc+t4*A3MMP*A4MMPc
      t398 = t1*A4MPP*A1MPPc+t223*A5PMP*A5PMPc+t98*A5MPM*A4MPMc+t1*A4MPM
     #*A1MPMc+t20*A5MPP*A1MPPc+t7*A3MPM*A1MPMc+t34*A1PPP*A6PPPc+t67*A2MP
     #P*A6MPPc+t119*A2PMM*A1PMMc+t223*A5MMP*A5MMPc+t1*A1MMP*A4MMPc+t10*A
     #3PPP*A2PPPc+t4*A3MPM*A4MPMc+t55*A3MMM*A5MMMc+t98*A4MMP*A5MMPc+t124
     #*A3MPP*A3MPPc+t138*A6PPM*A3PPMc+t79*A2MPP*A2MPPc
      t436 = t34*A1PMP*A6PMPc+t74*A4MMM*A4MMMc+t119*A2MMM*A1MMMc+t79*A2P
     #MM*A2PMMc+t15*A2MPM*A4MPMc+t44*A5PPP*A2PPPc+t67*A2PPM*A6PPMc+t64*A
     #5MPP*A6MPPc+t55*A3PPM*A5PPMc+t4*A3MMM*A4MMMc+t67*A6MMP*A2MMPc+t124
     #*A3MMP*A3MMPc+t15*A4MMM*A2MMMc+t67*A6MPM*A2MPMc+t10*A2MMP*A3MMPc+t
     #74*A4MMP*A4MMPc+t48*A1MMP*A1MMPc+t7*A1PMP*A3PMPc
      t473 = t223*A5MMM*A5MMMc+t67*A6PMM*A2PMMc+t55*A3MPM*A5MPMc+t48*A1P
     #MP*A1PMPc+t74*A4PMM*A4PMMc+t10*A3PMP*A2PMPc+t55*A5MMP*A3MMPc+t7*A3
     #PMP*A1PMPc+t67*A2PMP*A6PMPc+t98*A5PMM*A4PMMc+t124*A3MMM*A3MMMc+t20
     #*A5MMP*A1MMPc+t119*A1PMP*A2PMPc+t41*A6PMM*A6PMMc+t64*A6PMM*A5PMMc+
     #t1*A1PMM*A4PMMc+t27*A6PMM*A4PMMc+t44*A2PMP*A5PMPc
      t512 = t138*A6PMM*A3PMMc+t64*A5PMM*A6PMMc+t119*A2PMP*A1PMPc+t15*A2
     #PMP*A4PMPc+t34*A1PMM*A6PMMc+t223*A5PMM*A5PMMc+t138*A3MPM*A6MPMc+t1
     #*A4PMM*A1PMMc+t4*A3PMP*A4PMPc+t15*A4PPM*A2PPMc+t7*A3PPP*A1PPPc+t67
     #*A2PMM*A6PMMc+t64*A5PMP*A6PMPc+t119*A1MMP*A2MMPc+t223*A5PPM*A5PPMc
     #+t44*A5PMP*A2PMPc+t44*A5MMP*A2MMPc+t20*A1MMP*A5MMPc
      t549 = t44*A2PMM*A5PMMc+t48*A1MPM*A1MPMc+t27*A6PMP*A4PMPc+t44*A2PP
     #P*A5PPPc+t48*A1MMM*A1MMMc+t20*A1MPP*A5MPPc+t64*A6MMM*A5MMMc+t7*A1P
     #PM*A3PPMc+t10*A2PPM*A3PPMc+t98*A5MPP*A4MPPc+t34*A1PPM*A6PPMc+t79*A
     #2MMM*A2MMMc+t20*A1PPM*A5PPMc+t67*A6PPM*A2PPMc+t119*A1PPM*A2PPMc+t4
     #8*A1PPM*A1PPMc+t41*A6PMP*A6PMPc+t64*A6PMP*A5PMPc
      t587 = t34*A6PMP*A1PMPc+t138*A3PMP*A6PMPc+t223*A5MPP*A5MPPc+t55*A3
     #PMP*A5PMPc+t1*A1MPM*A4MPMc+t124*A3PMP*A3PMPc+t119*A1MPM*A2MPMc+t74
     #*A4PPP*A4PPPc+t4*A4PPP*A3PPPc+t55*A3PPP*A5PPPc+t4*A3PPP*A4PPPc+t67
     #*A2PPP*A6PPPc+t124*A3MPM*A3MPMc+t15*A2PPP*A4PPPc+t79*A2PPP*A2PPPc+
     #t48*A1PPP*A1PPPc+t27*A6PPM*A4PPMc+t1*A1PPM*A4PPMc
      t624 = t34*A6PPM*A1PPMc+t27*A4PPM*A6PPMc+t74*A4PPM*A4PPMc+t138*A3P
     #PM*A6PPMc+t124*A3PPM*A3PPMc+t10*A3PPM*A2PPMc+t44*A2PPM*A5PPMc+t7*A
     #3PMM*A1PMMc+t7*A3MMP*A1MMPc+t15*A2MMP*A4MMPc+t119*A1PMM*A2PMMc+t20
     #*A5PPP*A1PPPc+t27*A4PPP*A6PPPc+t98*A4PPP*A5PPPc+t138*A6MPP*A3MPPc+
     #t1*A1PMP*A4PMPc+t98*A5MMP*A4MMPc+t138*A3PMM*A6PMMc
      dummy = t47+t90+t131+t169+t208+t246+t284+t321+t361+t398+t436+t473+
     #t512+t549+t587+t624
      Dqgq = Dqgq + dummy
C
C     exchange the vectors r1 <-> r2 it corresponds to sqb --> -sqb
C
        do j = 1,4
          dtemp = r1(j)
          r1(j) = r2(j)
          r2(j) = dtemp
        enddo
      enddo ! ispincount
      Dqgq = Dqgq * V
C     multiply with the coupling 
C     The additional factor 2. comes from the spin sum, r1 <-> r2 takes
C     takes only two of the four possibilities into account.
      Dqgq = 2.d0 * Dqgq * g2**3
      end
