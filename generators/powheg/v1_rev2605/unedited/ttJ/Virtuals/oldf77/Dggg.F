C********************************************************************
      function Dggg(iem,pn,C,a,b,nvec)
C********************************************************************
C     Automatically created on Wed Feb 16 11:41:56 CET 2005
C     by ./CreateDggg.csh
C---------------------------------------------------------------------------
C
C     This software is not designed or intended for use in on-line control of
C     aircraft, air traffic, aircraft navigation or aircraft communications; 
C     or in the design, construction, operation or maintenance of any nuclear
C     facility. 
C
C---------------------------------------------------------------------------
      implicit none
      double precision Dggg
      integer iem
      double precision pn(4,-2:3)
      double precision p1(4),p2(4),p3(4)
      double precision q1(4),q2(4),r1(4),r2(4),sq(4),sqb(4)
      double precision kq(4),kqb(4),sqrest(4),sqbrest(4),pem(4)
      double precision C(6,6) ! color correlation matrix
      double precision a ! coeff. of -g_(mu,nu)
      double precision b, nvec(4) ! spin correlation
      double precision n1(4),n2(4),svec(4), alpha,nvec_nvec
      double precision dtemp
      double precision q(4), phi, theta, rn
      double complex Em, Ep
      double complex s1,B1,B2,B3,B4,B5,B6
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
      double complex spa,spb
      double complex dummy
      double precision betaq,betaqb,mqq,hein1(4),hein2(4)
      double precision dotp,norm3
      integer mu,j,rho,ispincount,error
#include "Dggg.dec"
#include "common.F"
#include "constants.F"
      Dggg = 0.d0
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
      write(*,*) 'Dggg: iem = ',iem
      write(*,*) 'Dggg: mq1,mq2 = ',dsqrt(dotp(kq,kq)) 
     >     ,dsqrt(dotp(kqb,kqb)) 
      write(*,*) 'Dggg: a= ', a
      write(*,*) 'Dggg: b= ', b
      write(*,*) 'Dggg: nvec(mu)= ', nvec
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
         Dggg=0.d0
	 return
      endif
      call boost(kqb,sqbrest,sqb,error)  
      if (error.lt.0) then
         ierr=error
         Dggg=0.d0
	 return
      endif
#ifdef DEBUG
      write(*,*) 'Dggg: kq.sq = ',dotp(kq,sq) 
      write(*,*) 'Dggg: kqb.sqb = ',dotp(kqb,sqb)
      write(*,*) 'Dggg: sqrest.sqrest = ', dotp(sqrest,sqrest)  
      write(*,*) 'Dggg: sq.sq = ', dotp(sq,sq)  
      write(*,*) 'Dggg: sqbrest.sqbrest = ', dotp(sqbrest,sqbrest)  
      write(*,*) 'Dggg: sqb.sqb = ', dotp(sqb,sqb)  
#endif   
      call evalqr(kq,sq,kqb,sqb,q1,q2,r1,r2)
      Dggg = 0.d0
C
C     decompose the vector nvec in two massless vectors:
C
      nvec_nvec = dotp(nvec,nvec)
      if (nvec_nvec.lt.0.d0) then
         alpha = dsqrt(nvec(4)**2/norm3(nvec)**2/(-nvec_nvec))
         do j=1,3
            svec(j) = nvec(j) * alpha
         enddo
         svec(4) = alpha * norm3(nvec)**2/nvec(4)
         alpha = dsqrt(-nvec_nvec)
         do j=1,4
            n1(j) = 0.5d0 * (nvec(j) + alpha * svec(j))
            n2(j) = 0.5d0 * (nvec(j) - alpha * svec(j))
#ifdef DEBUG
         write(*,*) 'Dggg: ',nvec(j),' = ', n1(j)+n2(j)	     
#endif
         enddo
#ifdef DEBUG
         write(*,*) 'Dggg: ',dotp(nvec,nvec)
         write(*,*) 'Dggg: ',dotp(n1,n1)+2.d0 * dotp(n1,n2)
     .        + dotp(n2,n2)
         write(*,*) 'Dggg: 1 =',dotp(svec,svec)
         write(*,*) 'Dggg: 0 =',dotp(svec,nvec)
         write(*,*) 'Dggg: 0 =',dotp(n1,n1)
         write(*,*) 'Dggg: 0 =',dotp(n2,n2)
#endif
C
C       use n2 as reference momenta for helicty
#ifdef DEBUG
         write(*,*) 'ktij.n = ', dotp(pem,nvec)
#endif
        Ep = spa(n2,n1)*spb(n1,pem)/sqrt2/spa(n2,pem) 
#ifdef DEBUG 
         write(*,*) 'EP1 = ',Ep
#endif
C
C     Different choices for the reference momentum are possible:
C         if (iem.eq.3) then
C            Ep = spa(p1,n1)*spb(n1,pem)/sqrt2/spa(p1,pem) 
C      .        + spa(p1,n2)*spb(n2,pem)/sqrt2/spa(p1,pem) 
C         else
C            Ep = spa(p3,n1)*spb(n1,pem)/sqrt2/spa(p3,pem) 
C      .        + spa(p3,n2)*spb(n2,pem)/sqrt2/spa(p3,pem) 
C         endif
C         write(*,*) 'EP2 = ',Ep                
C         theta = pi * rn(0.1233d0)
C         phi = 2.d0 * pi * rn(0.2331d0)
C         q(1) = dsin(theta)*dcos(phi)
C         q(2) = dsin(theta)*dsin(phi)
C         q(3) = dcos(theta)
C         q(4) = 1.
C         Ep = spa(q,n1)*spb(n1,pem)/sqrt2/spa(q,pem) 
C      .        + spa(q,n2)*spb(n2,pem)/sqrt2/spa(q,pem)
C         write(*,*) 'Ep3 = ',Ep 
C       Em = - spb(n2,n1)*spa(n1,pem)/sqrt2/spb(n2,pem)
C     somehow I skrewed up my notation....:
         Em = Ep
         Ep = conjg(Ep)
      else
         write(*,*) 'Dggg: dotp(n,n) =', dotp(nvec,nvec)
         write(*,*) 'Dggg: not yet implemnted'
	 ierr = -2
	 Dggg = 0.d0
	 return
      endif
#ifdef DEBUG
      write(*,*) 'Dggg: some checks on the kinematics'
      write(*,*) 'Dggg: 0 = p1.p1 = ', dotp(p1,p1)
      write(*,*) 'Dggg: 0 = p2.p2 = ', dotp(p2,p2)
      write(*,*) 'Dggg: 0 = p3.p3 = ', dotp(p3,p3)
      write(*,*) 'Dggg: mq^2 = 2*q1.q2 = ', 2.d0 * dotp(q1,q2)
      write(*,*) 'Dggg: mq^2 = 2*r1.r2 = ', 2.d0 * dotp(r1,r2)
      do j=1,4
         hein1(j) = p1(j)+p2(j)
	 hein2(j) = -p3(j) + q1(j) + q2(j) + r1(j) + r2(j)
         write(*,*) 'Dggg:', p1(j)+p2(j)+p3(j), ' = ',
     .        q1(j) + q2(j) + r1(j) + r2(j)
      enddo
      write(*,*) 'Dggg: s = ', dotp(hein1,hein1)
      write(*,*) 'Dggg: s = ', dotp(hein2,hein2)
#endif 
C
C     We don't study polarized quarks so we have to sum over all
C     the polarization. The easiest way to do that is buy just
C     caclulating 2 * Dggg(q1,q2,r1,r2) + Dggg(q1,q2,r2,r1)
C     Note that this works only if there is no term linear in the 
C     spin vector, with other word we assume the structure
C       Dggg = (indep.st,stb) +  (something * st.stb) 
C     To add the two configuration we just run two times the following
C     code and exchange at the end r1 <--> r2
      do ispincount=1,2
C
C     There are six different permutations of the partons which we just
C     label 1-6. For every permutation we have in addition 8 different
C     helicity configuration. All together we have thus 48 amplitides
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
      if (iem.eq.-1) then
      B1 = Em*A1MMM+Ep*A1PMM
      B2 = Em*A2MMM+Ep*A2PMM
      B3 = Em*A3MMM+Ep*A3PMM
      B4 = Em*A4MMM+Ep*A4PMM
      B5 = Em*A5MMM+Ep*A5PMM
      B6 = Em*A6MMM+Ep*A6PMM
      endif
      if (iem.eq.-2) then
      B1 = Em*A1MMM+Ep*A1MPM
      B2 = Em*A2MMM+Ep*A2MPM
      B3 = Em*A3MMM+Ep*A3MPM
      B4 = Em*A4MMM+Ep*A4MPM
      B5 = Em*A5MMM+Ep*A5MPM
      B6 = Em*A6MMM+Ep*A6MPM
      endif
      if (iem.eq.3) then
      B1 = Em*A1MMM+Ep*A1MMP
      B2 = Em*A2MMM+Ep*A2MMP
      B3 = Em*A3MMM+Ep*A3MMP
      B4 = Em*A4MMM+Ep*A4MMP
      B5 = Em*A5MMM+Ep*A5MMP
      B6 = Em*A6MMM+Ep*A6MMP
      endif
      t6 = b*B1
      t7 = conjg(B1)
      t11 = C(1,4)
      t16 = b*B4
      t20 = C(1,2)
      t25 = b*B2
      t34 = conjg(B2)
      t38 = C(3,5)
      t43 = b*B3
      t44 = conjg(B5)
      t53 = conjg(B4)
      t64 = C(3,4)
      t69 = conjg(B3)
      t73 = C(3,6)
      t78 = conjg(B6)
      t86 = b*B5
      t90 = C(1,3)
      t98 = C(2,5)
      t121 = C(2,4)
      t143 = C(4,5)
      s1 = C(1,1)*(a*(A1MMM*A1MMMc+A1MPM*A1MPMc)+t6*t7)+t11*(a*(A4MMM*A1
     #MMMc+A4MPM*A1MPMc)+t16*t7)+t20*(a*(A2MMM*A1MMMc+A2MPM*A1MPMc)+t25*
     #t7)+C(2,2)*(a*(A2MMM*A2MMMc+A2MPM*A2MPMc)+t25*t34)+t38*(a*(A3MMM*A
     #5MMMc+A3MPM*A5MPMc)+t43*t44)+C(4,4)*(a*(A4MMM*A4MMMc+A4MPM*A4MPMc)
     #+t16*t53)+t20*(a*(A1MMM*A2MMMc+A1MPM*A2MPMc)+t6*t34)+t64*(a*(A4MMM
     #*A3MMMc+A4MPM*A3MPMc)+t16*t69)+t73*(a*(A3MMM*A6MMMc+A3MPM*A6MPMc)+
     #t43*t78)
      t151 = s1+t38*(a*(A5MMM*A3MMMc+A5MPM*A3MPMc)+t86*t69)+t90*(a*(A1MM
     #M*A3MMMc+A1MPM*A3MPMc)+t6*t69)+t98*(a*(A5MMM*A2MMMc+A5MPM*A2MPMc)+
     #t86*t34)+t64*(a*(A3MMM*A4MMMc+A3MPM*A4MPMc)+t43*t53)+C(3,3)*(a*(A3
     #MMM*A3MMMc+A3MPM*A3MPMc)+t43*t69)+t121*(a*(A4MMM*A2MMMc+A4MPM*A2MP
     #Mc)+t16*t34)+t90*(a*(A3MMM*A1MMMc+A3MPM*A1MPMc)+t43*t7)+t11*(a*(A1
     #MMM*A4MMMc+A1MPM*A4MPMc)+t6*t53)+t143*(a*(A4MMM*A5MMMc+A4MPM*A5MPM
     #c)+t16*t44)
      t159 = C(1,5)
      t167 = C(4,6)
      t189 = C(2,6)
      t204 = C(5,6)
      t224 = b*B6
      t235 = C(1,6)
      t250 = C(2,3)
      s1 = t98*(a*(A2MMM*A5MMMc+A2MPM*A5MPMc)+t25*t44)+t159*(a*(A5MMM*A1
     #MMMc+A5MPM*A1MPMc)+t86*t7)+t167*(a*(A4MMM*A6MMMc+A4MPM*A6MPMc)+t16
     #*t78)+t143*(a*(A5MMM*A4MMMc+A5MPM*A4MPMc)+t86*t53)+t159*(a*(A1MMM*
     #A5MMMc+A1MPM*A5MPMc)+t6*t44)+t189*(a*(A2MMM*A6MMMc+A2MPM*A6MPMc)+t
     #25*t78)+t121*(a*(A2MMM*A4MMMc+A2MPM*A4MPMc)+t25*t53)+t204*(a*(A5MM
     #M*A6MMMc+A5MPM*A6MPMc)+t86*t78)+C(5,5)*(a*(A5MMM*A5MMMc+A5MPM*A5MP
     #Mc)+t86*t44)
      t287 = s1+t73*(a*(A6MMM*A3MMMc+A6MPM*A3MPMc)+t224*t69)+t189*(a*(A6
     #MMM*A2MMMc+A6MPM*A2MPMc)+t224*t34)+t235*(a*(A6MMM*A1MMMc+A6MPM*A1M
     #PMc)+t224*t7)+t204*(a*(A6MMM*A5MMMc+A6MPM*A5MPMc)+t224*t44)+t250*(
     #a*(A2MMM*A3MMMc+A2MPM*A3MPMc)+t25*t69)+t235*(a*(A1MMM*A6MMMc+A1MPM
     #*A6MPMc)+t6*t78)+C(6,6)*(a*(A6MMM*A6MMMc+A6MPM*A6MPMc)+t224*t78)+t
     #167*(a*(A6MMM*A4MMMc+A6MPM*A4MPMc)+t224*t53)+t250*(a*(A3MMM*A2MMMc
     #+A3MPM*A2MPMc)+t43*t34)
      dummy = t151+t287
      Dggg = Dggg + dummy
      if (iem.eq.-1) then
      B1 = Em*A1MMP+Ep*A1PMP
      B2 = Em*A2MMP+Ep*A2PMP
      B3 = Em*A3MMP+Ep*A3PMP
      B4 = Em*A4MMP+Ep*A4PMP
      B5 = Em*A5MMP+Ep*A5PMP
      B6 = Em*A6MMP+Ep*A6PMP
      endif
      if (iem.eq.-2) then
      B1 = Em*A1MMP+Ep*A1MPP
      B2 = Em*A2MMP+Ep*A2MPP
      B3 = Em*A3MMP+Ep*A3MPP
      B4 = Em*A4MMP+Ep*A4MPP
      B5 = Em*A5MMP+Ep*A5MPP
      B6 = Em*A6MMP+Ep*A6MPP
      endif
      if (iem.eq.3) then
      B1 = Em*A1MPM+Ep*A1MPP
      B2 = Em*A2MPM+Ep*A2MPP
      B3 = Em*A3MPM+Ep*A3MPP
      B4 = Em*A4MPM+Ep*A4MPP
      B5 = Em*A5MPM+Ep*A5MPP
      B6 = Em*A6MPM+Ep*A6MPP
      endif
      t1 = C(1,3)
      t6 = b*B1
      t7 = conjg(B3)
      t11 = C(4,5)
      t16 = b*B4
      t17 = conjg(B5)
      t26 = conjg(B4)
      t30 = C(3,6)
      t35 = b*B3
      t36 = conjg(B6)
      t40 = C(1,4)
      t45 = conjg(B1)
      t49 = C(2,4)
      t54 = conjg(B2)
      t58 = C(4,6)
      t66 = C(2,3)
      t71 = b*B2
      t105 = C(3,5)
      t110 = b*B5
      t135 = C(1,2)
      s1 = t1*(a*(A1MMP*A3MMPc+A1MPP*A3MPPc)+t6*t7)+t11*(a*(A4MMP*A5MMPc
     #+A4MPP*A5MPPc)+t16*t17)+C(4,4)*(a*(A4MMP*A4MMPc+A4MPP*A4MPPc)+t16*
     #t26)+t30*(a*(A3MMP*A6MMPc+A3MPP*A6MPPc)+t35*t36)+t40*(a*(A4MMP*A1M
     #MPc+A4MPP*A1MPPc)+t16*t45)+t49*(a*(A4MMP*A2MMPc+A4MPP*A2MPPc)+t16*
     #t54)+t58*(a*(A4MMP*A6MMPc+A4MPP*A6MPPc)+t16*t36)+t66*(a*(A2MMP*A3M
     #MPc+A2MPP*A3MPPc)+t71*t7)+C(3,3)*(a*(A3MMP*A3MMPc+A3MPP*A3MPPc)+t3
     #5*t7)
      t151 = s1+t66*(a*(A3MMP*A2MMPc+A3MPP*A2MPPc)+t35*t54)+t40*(a*(A1MM
     #P*A4MMPc+A1MPP*A4MPPc)+t6*t26)+C(1,1)*(a*(A1MMP*A1MMPc+A1MPP*A1MPP
     #c)+t6*t45)+t105*(a*(A5MMP*A3MMPc+A5MPP*A3MPPc)+t110*t7)+t1*(a*(A3M
     #MP*A1MMPc+A3MPP*A1MPPc)+t35*t45)+t49*(a*(A2MMP*A4MMPc+A2MPP*A4MPPc
     #)+t71*t26)+t11*(a*(A5MMP*A4MMPc+A5MPP*A4MPPc)+t110*t26)+t135*(a*(A
     #2MMP*A1MMPc+A2MPP*A1MPPc)+t71*t45)+C(2,2)*(a*(A2MMP*A2MMPc+A2MPP*A
     #2MPPc)+t71*t54)
      t152 = C(1,5)
      t160 = C(3,4)
      t182 = C(2,6)
      t190 = C(2,5)
      t209 = b*B6
      t228 = C(5,6)
      t236 = C(1,6)
      s1 = t152*(a*(A5MMP*A1MMPc+A5MPP*A1MPPc)+t110*t45)+t160*(a*(A4MMP*
     #A3MMPc+A4MPP*A3MPPc)+t16*t7)+t160*(a*(A3MMP*A4MMPc+A3MPP*A4MPPc)+t
     #35*t26)+t105*(a*(A3MMP*A5MMPc+A3MPP*A5MPPc)+t35*t17)+t182*(a*(A2MM
     #P*A6MMPc+A2MPP*A6MPPc)+t71*t36)+t190*(a*(A2MMP*A5MMPc+A2MPP*A5MPPc
     #)+t71*t17)+t190*(a*(A5MMP*A2MMPc+A5MPP*A2MPPc)+t110*t54)+t182*(a*(
     #A6MMP*A2MMPc+A6MPP*A2MPPc)+t209*t54)+C(5,5)*(a*(A5MMP*A5MMPc+A5MPP
     #*A5MPPc)+t110*t17)
      t287 = s1+t135*(a*(A1MMP*A2MMPc+A1MPP*A2MPPc)+t6*t54)+t228*(a*(A5M
     #MP*A6MMPc+A5MPP*A6MPPc)+t110*t36)+t236*(a*(A6MMP*A1MMPc+A6MPP*A1MP
     #Pc)+t209*t45)+C(6,6)*(a*(A6MMP*A6MMPc+A6MPP*A6MPPc)+t209*t36)+t228
     #*(a*(A6MMP*A5MMPc+A6MPP*A5MPPc)+t209*t17)+t58*(a*(A6MMP*A4MMPc+A6M
     #PP*A4MPPc)+t209*t26)+t30*(a*(A6MMP*A3MMPc+A6MPP*A3MPPc)+t209*t7)+t
     #152*(a*(A1MMP*A5MMPc+A1MPP*A5MPPc)+t6*t17)+t236*(a*(A1MMP*A6MMPc+A
     #1MPP*A6MPPc)+t6*t36)
      dummy = t151+t287
      Dggg = Dggg + dummy
      if (iem.eq.-1) then
      B1 = Em*A1MPM+Ep*A1PPM
      B2 = Em*A2MPM+Ep*A2PPM
      B3 = Em*A3MPM+Ep*A3PPM
      B4 = Em*A4MPM+Ep*A4PPM
      B5 = Em*A5MPM+Ep*A5PPM
      B6 = Em*A6MPM+Ep*A6PPM
      endif
      if (iem.eq.-2) then
      B1 = Em*A1PMM+Ep*A1PPM
      B2 = Em*A2PMM+Ep*A2PPM
      B3 = Em*A3PMM+Ep*A3PPM
      B4 = Em*A4PMM+Ep*A4PPM
      B5 = Em*A5PMM+Ep*A5PPM
      B6 = Em*A6PMM+Ep*A6PPM
      endif
      if (iem.eq.3) then
      B1 = Em*A1PMM+Ep*A1PMP
      B2 = Em*A2PMM+Ep*A2PMP
      B3 = Em*A3PMM+Ep*A3PMP
      B4 = Em*A4PMM+Ep*A4PMP
      B5 = Em*A5PMM+Ep*A5PMP
      B6 = Em*A6PMM+Ep*A6PMP
      endif
      t1 = C(4,5)
      t6 = b*B4
      t7 = conjg(B5)
      t11 = C(1,6)
      t16 = b*B6
      t17 = conjg(B1)
      t21 = C(1,4)
      t29 = C(5,6)
      t34 = b*B5
      t35 = conjg(B6)
      t44 = conjg(B4)
      t48 = C(3,5)
      t53 = b*B3
      t62 = conjg(B3)
      t71 = b*B1
      t75 = C(3,4)
      t83 = C(1,3)
      t91 = C(2,3)
      t96 = conjg(B2)
      t107 = C(2,5)
      t112 = b*B2
      t116 = C(2,6)
      t132 = C(2,4)
      s1 = t1*(a*(A4PMM*A5PMMc+A4PPM*A5PPMc)+t6*t7)+t11*(a*(A6PMM*A1PMMc
     #+A6PPM*A1PPMc)+t16*t17)+t21*(a*(A4PMM*A1PMMc+A4PPM*A1PPMc)+t6*t17)
     #+t29*(a*(A5PMM*A6PMMc+A5PPM*A6PPMc)+t34*t35)+C(4,4)*(a*(A4PMM*A4PM
     #Mc+A4PPM*A4PPMc)+t6*t44)+t48*(a*(A3PMM*A5PMMc+A3PPM*A5PPMc)+t53*t7
     #)+C(3,3)*(a*(A3PMM*A3PMMc+A3PPM*A3PPMc)+t53*t62)+C(1,1)*(a*(A1PMM*
     #A1PMMc+A1PPM*A1PPMc)+t71*t17)+t75*(a*(A3PMM*A4PMMc+A3PPM*A4PPMc)+t
     #53*t44)
      t154 = s1+t83*(a*(A3PMM*A1PMMc+A3PPM*A1PPMc)+t53*t17)+t91*(a*(A3PM
     #M*A2PMMc+A3PPM*A2PPMc)+t53*t96)+t11*(a*(A1PMM*A6PMMc+A1PPM*A6PPMc)
     #+t71*t35)+t107*(a*(A2PMM*A5PMMc+A2PPM*A5PPMc)+t112*t7)+t116*(a*(A2
     #PMM*A6PMMc+A2PPM*A6PPMc)+t112*t35)+C(2,2)*(a*(A2PMM*A2PMMc+A2PPM*A
     #2PPMc)+t112*t96)+t132*(a*(A2PMM*A4PMMc+A2PPM*A4PPMc)+t112*t44)+t13
     #2*(a*(A4PMM*A2PMMc+A4PPM*A2PPMc)+t6*t96)+t75*(a*(A4PMM*A3PMMc+A4PP
     #M*A3PPMc)+t6*t62)
      t155 = C(4,6)
      t163 = C(3,6)
      t222 = C(1,5)
      t237 = C(1,2)
      s1 = t155*(a*(A4PMM*A6PMMc+A4PPM*A6PPMc)+t6*t35)+t163*(a*(A3PMM*A6
     #PMMc+A3PPM*A6PPMc)+t53*t35)+t29*(a*(A6PMM*A5PMMc+A6PPM*A5PPMc)+t16
     #*t7)+t155*(a*(A6PMM*A4PMMc+A6PPM*A4PPMc)+t16*t44)+t116*(a*(A6PMM*A
     #2PMMc+A6PPM*A2PPMc)+t16*t96)+t163*(a*(A6PMM*A3PMMc+A6PPM*A3PPMc)+t
     #16*t62)+C(6,6)*(a*(A6PMM*A6PMMc+A6PPM*A6PPMc)+t16*t35)+C(5,5)*(a*(
     #A5PMM*A5PMMc+A5PPM*A5PPMc)+t34*t7)+t21*(a*(A1PMM*A4PMMc+A1PPM*A4PP
     #Mc)+t71*t44)
      t287 = s1+t222*(a*(A5PMM*A1PMMc+A5PPM*A1PPMc)+t34*t17)+t91*(a*(A2P
     #MM*A3PMMc+A2PPM*A3PPMc)+t112*t62)+t237*(a*(A2PMM*A1PMMc+A2PPM*A1PP
     #Mc)+t112*t17)+t107*(a*(A5PMM*A2PMMc+A5PPM*A2PPMc)+t34*t96)+t83*(a*
     #(A1PMM*A3PMMc+A1PPM*A3PPMc)+t71*t62)+t222*(a*(A1PMM*A5PMMc+A1PPM*A
     #5PPMc)+t71*t7)+t48*(a*(A5PMM*A3PMMc+A5PPM*A3PPMc)+t34*t62)+t237*(a
     #*(A1PMM*A2PMMc+A1PPM*A2PPMc)+t71*t96)+t1*(a*(A5PMM*A4PMMc+A5PPM*A4
     #PPMc)+t34*t44)
      dummy = t154+t287
      Dggg = Dggg + dummy
      if (iem.eq.-1) then
      B1 = Em*A1MPP+Ep*A1PPP
      B2 = Em*A2MPP+Ep*A2PPP
      B3 = Em*A3MPP+Ep*A3PPP
      B4 = Em*A4MPP+Ep*A4PPP
      B5 = Em*A5MPP+Ep*A5PPP
      B6 = Em*A6MPP+Ep*A6PPP
      endif
      if (iem.eq.-2) then
      B1 = Em*A1PMP+Ep*A1PPP
      B2 = Em*A2PMP+Ep*A2PPP
      B3 = Em*A3PMP+Ep*A3PPP
      B4 = Em*A4PMP+Ep*A4PPP
      B5 = Em*A5PMP+Ep*A5PPP
      B6 = Em*A6PMP+Ep*A6PPP
      endif
      if (iem.eq.3) then
      B1 = Em*A1PPM+Ep*A1PPP
      B2 = Em*A2PPM+Ep*A2PPP
      B3 = Em*A3PPM+Ep*A3PPP
      B4 = Em*A4PPM+Ep*A4PPP
      B5 = Em*A5PPM+Ep*A5PPP
      B6 = Em*A6PPM+Ep*A6PPP
      endif
      t1 = C(1,2)
      t6 = b*B1
      t7 = conjg(B2)
      t16 = conjg(B1)
      t20 = C(2,6)
      t25 = b*B6
      t29 = C(3,6)
      t34 = conjg(B3)
      t38 = C(1,6)
      t51 = b*B3
      t55 = C(1,4)
      t60 = conjg(B4)
      t64 = C(3,4)
      t69 = b*B4
      t73 = C(5,6)
      t78 = b*B5
      t79 = conjg(B6)
      t87 = conjg(B5)
      t91 = C(1,5)
      t115 = C(3,5)
      t123 = C(2,3)
      t128 = b*B2
      t132 = C(2,4)
      t140 = C(2,5)
      t148 = C(1,3)
      s1 = t1*(a*(A1PMP*A2PMPc+A1PPP*A2PPPc)+t6*t7)+C(1,1)*(a*(A1PMP*A1P
     #MPc+A1PPP*A1PPPc)+t6*t16)+t20*(a*(A6PMP*A2PMPc+A6PPP*A2PPPc)+t25*t
     #7)+t29*(a*(A6PMP*A3PMPc+A6PPP*A3PPPc)+t25*t34)+t38*(a*(A6PMP*A1PMP
     #c+A6PPP*A1PPPc)+t25*t16)+C(3,3)*(a*(A3PMP*A3PMPc+A3PPP*A3PPPc)+t51
     #*t34)+t55*(a*(A1PMP*A4PMPc+A1PPP*A4PPPc)+t6*t60)+t64*(a*(A4PMP*A3P
     #MPc+A4PPP*A3PPPc)+t69*t34)+t73*(a*(A5PMP*A6PMPc+A5PPP*A6PPPc)+t78*
     #t79)
      t156 = s1+t73*(a*(A6PMP*A5PMPc+A6PPP*A5PPPc)+t25*t87)+t91*(a*(A5PM
     #P*A1PMPc+A5PPP*A1PPPc)+t78*t16)+C(6,6)*(a*(A6PMP*A6PMPc+A6PPP*A6PP
     #Pc)+t25*t79)+C(5,5)*(a*(A5PMP*A5PMPc+A5PPP*A5PPPc)+t78*t87)+t115*(
     #a*(A5PMP*A3PMPc+A5PPP*A3PPPc)+t78*t34)+t123*(a*(A2PMP*A3PMPc+A2PPP
     #*A3PPPc)+t128*t34)+t132*(a*(A4PMP*A2PMPc+A4PPP*A2PPPc)+t69*t7)+t14
     #0*(a*(A5PMP*A2PMPc+A5PPP*A2PPPc)+t78*t7)+t148*(a*(A1PMP*A3PMPc+A1P
     #PP*A3PPPc)+t6*t34)
      t165 = C(4,5)
      t208 = C(4,6)
      s1 = C(4,4)*(a*(A4PMP*A4PMPc+A4PPP*A4PPPc)+t69*t60)+t165*(a*(A5PMP
     #*A4PMPc+A5PPP*A4PPPc)+t78*t60)+t148*(a*(A3PMP*A1PMPc+A3PPP*A1PPPc)
     #+t51*t16)+t64*(a*(A3PMP*A4PMPc+A3PPP*A4PPPc)+t51*t60)+t115*(a*(A3P
     #MP*A5PMPc+A3PPP*A5PPPc)+t51*t87)+t29*(a*(A3PMP*A6PMPc+A3PPP*A6PPPc
     #)+t51*t79)+t55*(a*(A4PMP*A1PMPc+A4PPP*A1PPPc)+t69*t16)+t208*(a*(A4
     #PMP*A6PMPc+A4PPP*A6PPPc)+t69*t79)+t165*(a*(A4PMP*A5PMPc+A4PPP*A5PP
     #Pc)+t69*t87)
      t287 = s1+t1*(a*(A2PMP*A1PMPc+A2PPP*A1PPPc)+t128*t16)+t38*(a*(A1PM
     #P*A6PMPc+A1PPP*A6PPPc)+t6*t79)+C(2,2)*(a*(A2PMP*A2PMPc+A2PPP*A2PPP
     #c)+t128*t7)+t20*(a*(A2PMP*A6PMPc+A2PPP*A6PPPc)+t128*t79)+t140*(a*(
     #A2PMP*A5PMPc+A2PPP*A5PPPc)+t128*t87)+t132*(a*(A2PMP*A4PMPc+A2PPP*A
     #4PPPc)+t128*t60)+t123*(a*(A3PMP*A2PMPc+A3PPP*A2PPPc)+t51*t7)+t91*(
     #a*(A1PMP*A5PMPc+A1PPP*A5PPPc)+t6*t87)+t208*(a*(A6PMP*A4PMPc+A6PPP*
     #A4PPPc)+t25*t60)
      dummy = t156+t287
      Dggg = Dggg + dummy
#ifdef DEBUG
      write(*,*) 'Dggg: b1,...', b1,b2,b3,b4,b5,b6
#endif
C     exchange the vectors r1 <-> r2 it corresponds to sqb --> -sqb
      do rho = 1,4
        dtemp = r1(rho)
        r1(rho) = r2(rho)
        r2(rho) = dtemp
      enddo
      enddo
C     Multiply with the coupling g^6
C     The additional factor 2. comes from the spin sum, r1 <-> r2 takes
C     takes only two of the four possibilities into account.
      Dggg = 2.d0 * Dggg * g2**3      
      end
