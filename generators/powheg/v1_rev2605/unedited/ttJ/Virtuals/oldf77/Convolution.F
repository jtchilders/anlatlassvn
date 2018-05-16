C     $Modified: Mon Aug  9 16:08:38 2004 by puwer $
C     Finite contribution from factorization+coll.-subtraction
***#define COMPARE 
      function conv(x)
      implicit none
      double precision conv
      double precision conva,convb
      double precision conva1,convb1
      double precision convax,convbx
      double precision x(6)
      double precision z,pa(4),pb(4),kq(4),kqb(4),kg(4)
      integer icont
      integer iset
      double precision setz(2)
      double precision pqqreg,pgqreg,pqgreg,pggreg
      double precision pggnreg
      double precision psqq,psqg,psgq,psgg
      double precision kax(4),kbx(4),kqx(4),kqbx(4),kgx(4)
      double precision ka1(4),kb1(4),kq1(4),kqb1(4),kg1(4)
      double precision squarex,corra1x,corrb1x,corrabx,
     &     corraqx,corraqbx,corrbqx,corrbqbx
      double precision square1,corra11,corrb11,corrab1,
     &     corraq1,corraqb1,corrbq1,corrbqb1
      double precision jacobx,jacob1,ecmsx,ecms1
      double precision sv,mj,mj2,sja
      double precision APa,APb,Kconta,Kcontb
      double precision APax,APbx,Kcontax,Kcontbx
      double precision APa1,APb1,Kconta1,Kcontb1
      double precision g_q,g_g,K_q,K_g,nf
      double precision Tgq, Tqq
      double precision Kgq_q_reg, Kgq_q_p, Kqq_q_reg,Kqq_q_p
      double precision Kqg_q_reg,Kqg_q_p,Kgg_q_reg,Kgg_q_p,kqq_q_pp
      double precision Kgg_g_reg,Kgg_g_p, Kgg_g_pp
      double precision JgQ,JQQ
      double precision mf,muqq,muqqx,muqq1
      double precision Kbargg_reg, Kbargg_p
      double precision suba, subb
      double precision dotp,theta
      double precision Kcal,delta,deltaa,deltab,APdeltaa
      double precision m(3),pset(4,3),eventx(4,-2:3),event1(4,-2:3)
      double precision obs3, obsval
      integer err,mu
      logical linit
      data linit/.true./
#include "common.F"
#include "constants.F"
#include "IDipolsDef.F"
C
      if (linit) then
         write(*,*) 'Check the definition of Nf'
         write(*,*) 'TO DO'
         write(*,*) 'ADD CONTRIBUTION FROM +-PRESCRIPTION'
         write(*,*) 'ADD PHASE SPACE FACTOR IN CORRX,CORR1...'
         linit = .false.
      endif
C      
#ifdef COMPARE
C
C     Two kinematical configuration each for x=1 and x <> 1.
C
      setz(1) = 0.7d0
      setz(2) = 0.3d0
      iset = 1
      z = setz(mod(iset,10))
C
C     z <> 1
C
      call stefan(iset,z,kax,kbx,kqx,kqbx,kgx)
      call ColorCorr(kax,kbx,kqx,kqbx,kgx,squarex,corra1x,
     &     corrb1x,corrabx,corraqx,corraqbx,corrbqx,corrbqbx) 
C
       squarex=0.d0
*       corra1x = 0.d0  ! T_2 T_1
       corrb1x = 0.d0
*       corrabx = 0.d0  ! T_2 T_3
*       corraqx = 0.d0  ! T_2 T_0
*       corraqbx = 0.d0 ! T_2 T_4 
       corrbqx = 0.d0
       corrbqbx =0.d0
C
      write(*,'(A,I1,A)') ' --------- set ',iset,'  ------------'
      write(*,*) '--- color correlated MEs  ---'
      write(*,*) '--------- z <> 1 ------------'
      write(*,*) '  M^2   = ', squarex
      write(*,*) 'T_2 T_0 = ', corraqx
      write(*,*) 'T_2 T_1 = ', corra1x
      write(*,*) 'T_2 T_3 = ', corrabx
      write(*,*) 'T_2 T_4 = ', corraqbx
C
C     Check the color decomposition of the squared matrix element
C
      write(*,*) -1.d0/CA*(corraqx+corraqbx+corrabx+corra1x)
      write(*,*) squarex
C
C     Decompose the squared matrix element according to the color structure
C
      squarex = -1.d0/CA*(corraqx+corraqbx+corrabx+corra1x)
C
C     z = 1
C
      call stefan(iset,1.d0,ka1,kb1,kq1,kqb1,kg1)
      call ColorCorr(ka1,kb1,kq1,kqb1,kg1,square1,corra11,
     &     corrb11,corrab1,corraq1,corraqb1,corrbq1,corrbqb1)
C       square1=0.d0
*       corra11=0.d0  ! T_2 T_1
       corrb11=0.d0
*       corrab1=0.d0  ! T_2 T_3
*       corraq1=1.d0  ! T_2 T_0
*       corraqb1=0.d0 ! T_2 T_4
       corrbq1 =0.d0
       corrbqb1=0.d0
C
      write(*,*) '--- color correlated MEs  ---'
      write(*,*) '--------- z = 1 ------------'
      write(*,*) '  M^2   = ', square1
      write(*,*) 'T_2 T_0 = ', corraq1
      write(*,*) 'T_2 T_1 = ', corra11
      write(*,*) 'T_2 T_3 = ', corrab1
      write(*,*) 'T_2 T_4 = ', corraqb1
C
C     Check the color decomposition of the squared matrix element
C
      write(*,*) -1.d0/CA*(corraq1+corraqb1+corrab1+corra11)
      write(*,*) square1
C
C     Decompose the squared matrix element according to the color structure
C
      square1 = -1.d0/CA*(corraq1+corraqb1+corrab1+corra11)
C
      write(*,*) '----------------------------'
      write(*,*) 'mt       = ', mt,dsqrt(dotp(kqx,kqx))
     &     ,dsqrt(dotp(kqx,kqx))
      write(*,*) 'z       = ', z
      write(*,*) 'mu_fac  = ', mu_fac
C
#else
C     
      z = x(6)
C
C     Check if we have enough energy to produce the ttbar pair in
C     the sub-process, if not we return 0
C
      if ((4.d0*mt**2).gt.z*s) then
         conv = 0.d0
         return
      endif
C
C     Calculate the final state momenta for 
C     z*ka+kb -->  kt + ktb + kg, we work in (z ka,kb) CMS
C    
      ecmsx = sqrt_s * dsqrt(z)

C
      kax(1) = 0.d0
      kax(2) = ecmsx/2.d0
      kax(3) = 0.d0
      kax(4) = ecmsx/2.d0

      kbx(1) = 0.d0
      kbx(2) = -ecmsx/2.d0
      kbx(3) = 0.d0
      kbx(4) = ecmsx/2.d0

      m(1) = mt
      m(2) = mt
      m(3) = 0.d0

      call eventn(3,s*z,m,x,pset,jacobx,err)

      do mu=1,4
         eventx(mu,-1) = kax(mu)
         eventx(mu,-2) = kbx(mu)
         eventx(mu,1) = pset(mu,1)
         eventx(mu,2) = pset(mu,2)
         eventx(mu,3) = pset(mu,3)
         
         kqx(mu) = pset(mu,1)
         kqbx(mu) = pset(mu,2)
         kgx(mu) = pset(mu,3)
      enddo

      if (jacobx.gt.0.d0) then
         jacobx = jacobx * obs3(eventx)
         if (jacobx.ne.0.d0) then
            call ColorCorr(kax,kbx,kqx,kqbx,kgx,squarex,corra1x,
     &           corrb1x,corrabx,corraqx,corraqbx,corrbqx,corrbqbx) 

            squarex = squarex * jacobx 
            corra1x = corra1x * jacobx
            corrb1x = corrb1x * jacobx
            corrabx = corrabx * jacobx
            corraqx = corraqx * jacobx
            corraqbx = corraqbx * jacobx
            corrbqx = corrbqx * jacobx
            corrbqbx = corrbqbx * jacobx
         endif
      endif
C
      ecms1 = sqrt_s 
      ka1(1) = 0.d0
      ka1(2) = ecms1/2.d0
      ka1(3) = 0.d0
      ka1(4) = ecms1/2.d0

      kb1(1) = 0.d0
      kb1(2) = -ecms1/2.d0
      kb1(3) = 0.d0
      kb1(4) = ecms1/2.d0

      m(1) = mt
      m(2) = mt
      m(3) = 0.d0

      call eventn(3,s,m,x,pset,jacob1,err)

      do mu=1,4
         event1(mu,-1) = ka1(mu)
         event1(mu,-2) = kb1(mu)
         event1(mu,1) = pset(mu,1)
         event1(mu,2) = pset(mu,2)
         event1(mu,3) = pset(mu,3)

         kq1(mu) = pset(mu,1)
         kqb1(mu) = pset(mu,2)
         kg1(mu) = pset(mu,3)
      enddo

      

      if (jacob1.gt.0.d0) then
         jacob1 = jacob1 * obs3(event1)
         if (jacob1.ne.0) then
            call ColorCorr(ka1,kb1,kq1,kqb1,kg1,square1,corra11,
     &           corrb11,corrab1,corraq1,corraqb1,corrbq1,corrbqb1)  
            
            square1 = square1 * jacob1
            corra11 = corra11 * jacob1
            corrb11 = corrb11 * jacob1
            corrab1 = corrab1 * jacob1
            corraq1 = corraq1 * jacob1
            corraqb1 = corraqb1 * jacob1
            corrbq1 = corrbq1 * jacob1
            corrbqb1 = corrbqb1 * jacob1
         endif
      endif
#endif
      Apax = 0.d0
      Apbx = 0.d0
      Kcontax = 0.d0
      Kcontbx = 0.d0
      if (jacobx.ne.0.d0) then
C     Calculate first the contribution from z<>1:
C
C     Calculate the contribution from the AP-kernel
C     given that we work in the kax,kb CMS, we can calculate also
C     the second contribution where b is splitting.
C
      APax = Pggreg(z) * (
     &     + corra1x  * dlog(mu_fac**2/2.d0/dotp(kgx,kax))
     &     + corraqx  * dlog(mu_fac**2/2.d0/dotp(kqx,kax))
     &     + corraqbx * dlog(mu_fac**2/2.d0/dotp(kqbx,kax))
     &     + corrabx  * dlog(mu_fac**2/2.d0/dotp(kbx,kax)) ) 
C
      APbx = Pggreg(z) * (
     &     + corrb1x  * dlog(mu_fac**2/2.d0/dotp(kgx,kbx))
     &     + corrbqx  * dlog(mu_fac**2/2.d0/dotp(kqx,kbx))
     &     + corrbqbx * dlog(mu_fac**2/2.d0/dotp(kqbx,kbx))
     &     + corrabx  * dlog(mu_fac**2/2.d0/dotp(kax,kbx)) ) 
C
C     Now the terms involving the +-prescription
C
      APax = APax + Pggnreg(z) * ( 
     &     + corra1x*dlog(mu_fac**2/2.d0/dotp(kgx,kax))
     &     + corraqx  * dlog(mu_fac**2/2.d0/dotp(kqx,kax))
     &     + corraqbx * dlog(mu_fac**2/2.d0/dotp(kqbx,kax)) 
     &     + corrabx * dlog(mu_fac**2/2.d0/dotp(kbx,kax))
     &     )
C
      APbx = APbx + Pggnreg(z) * ( 
     &     + corrb1x*dlog(mu_fac**2/2.d0/dotp(kgx,kbx))
     &     + corrbqx  * dlog(mu_fac**2/2.d0/dotp(kqx,kbx))
     &     + corrbqbx * dlog(mu_fac**2/2.d0/dotp(kqbx,kbx)) 
     &     + corrabx * dlog(mu_fac**2/2.d0/dotp(kax,kbx))
     &     )

C
      APax = alphas/2.d0/pi/CA*APax
      APbx = alphas/2.d0/pi/CA*APbx


C     
C     Calculate the K contribution, first the regular part:
C
      Kcontax = 
C     First term:
     &     Kbargg_reg(z) * squarex
C     First sum:
     &     - Kgg_g_reg(z,2.d0*dotp(kax,kgx)/z,0.d0) * corra1x
     &     - Kgg_q_reg(z,2.d0*dotp(kax,kqx)/z,mt)   * corraqx
     &     - Kgg_q_reg(z,2.d0*dotp(kax,kqbx)/z,mt)  * corraqbx
C     Second sum:
C     -- j = gluon
C     the contribution vanishes due to m_g = 0
C     -- j = q
     &     -1.d0/CA * corraqx * Pggreg(z) * 
     &     dlog((1.d0-z)* 2.d0 *dotp(kqx,kax)/z
     &     /((1.d0-z)* 2.d0 *dotp(kqx,kax)/z+mt**2))
C     -- j = qbar
     &     -1.d0/CA * corraqbx * Pggreg(z) * 
     &     dlog((1.d0-z)*2.d0*dotp(kqbx,kax)/z/
     &     ((1.d0-z)*2.d0*dotp(kqbx,kax)/z+mt**2))
C     -- j = incoming gluon
     &     -1.d0/CA * corrabx * Pggreg(z) * dlog(1.d0-z) 

      Kcontbx = 
C     First term:
     &     Kbargg_reg(z) * squarex
C     First sum:
     &     - Kgg_g_reg(z,2.d0*dotp(kbx,kgx)/z,0.d0) * corrb1x
     &     - Kgg_q_reg(z,2.d0*dotp(kbx,kqx)/z,mt)   * corrbqx
     &     - Kgg_q_reg(z,2.d0*dotp(kbx,kqbx)/z,mt)  * corrbqbx
C     Second sum:
C     -- j = gluon
C     the contribution vanishes due to m_g = 0
C     -- j = q
     &     -1.d0/CA * corrbqx * Pggreg(z) * 
     &     dlog((1.d0-z)* 2.d0 *dotp(kqx,kbx)/z
     &     /((1.d0-z)* 2.d0 *dotp(kqx,kbx)/z+mt**2))
C     -- j = qbar
     &     -1.d0/CA * corrbqbx * Pggreg(z) * 
     &     dlog((1.d0-z)*2.d0*dotp(kqbx,kbx)/z/
     &     ((1.d0-z)*2.d0*dotp(kqbx,kbx)/z+mt**2))
C     -- j = incoming gluon
     &     -1.d0/CA * corrabx * Pggreg(z) * dlog(1.d0-z) 


C
C     Now the simple plus prescriptions
C
      Kcontax = Kcontax 
     &     + Kbargg_p(z) * squarex 
C     Kgg_g
     &     - Kgg_g_p(z,2.d0*dotp(kax,kgx)/z,0.d0) * corra1x
C     Kgg_q
     &     - (
     &     2.d0*dlog(1.d0-z)/(1.d0-z) * (corraqx)
     &     + 2.d0/(1.d0-z) *( 
     &     dlog(((2.d0-z)*2.d0*dotp(kqx,kax)/z)/
     &     ((2.d0-z)*2.d0*dotp(kqx,kax)/z+mt**2))*corraqx
     &     )
     &     )
C     Kgg_qb
     &     - (
     &     2.d0*dlog(1.d0-z)/(1.d0-z) * (corraqbx)
     &     + 2.d0/(1.d0-z) *( 
     &     dlog(((2.d0-z)*2.d0*dotp(kqbx,kax)/z)/
     &     ((2.d0-z)*2.d0*dotp(kqbx,kax)/z+mt**2))*corraqbx
     &     )
     &     )
C     From the inclusion of b
     &     - 2.d0*(dlog(1.d0-z)/(1.d0-z))*(corrabx)

C
      Kcontbx = Kcontbx 
     &     + Kbargg_p(z) * ( squarex)
C     Kgg_g
     &     - Kgg_g_p(z,2.d0*dotp(kbx,kgx)/z,0.d0) * (corrb1x)
C     Kgg_q
     &     - (
     &     2.d0*dlog(1.d0-z)/(1.d0-z) * (corrbqx)
     &     + 2.d0/(1.d0-z) *( 
     &     dlog(((2.d0-z)*2.d0*dotp(kqx,kbx)/z)/
     &     ((2.d0-z)*2.d0*dotp(kqx,kbx)/z+mt**2))*corrbqx
     &     )
     &     )
C     Kgg_qb
     &     - (
     &     2.d0*dlog(1.d0-z)/(1.d0-z) * (corrbqbx)
     &     + 2.d0/(1.d0-z) *( 
     &     dlog(((2.d0-z)*2.d0*dotp(kqbx,kbx)/z)/
     &     ((2.d0-z)*2.d0*dotp(kqbx,kbx)/z+mt**2))*corrbqbx
     &     )
     &     )
C     From the inclusion of a
     &     - 2.d0*(dlog(1.d0-z)/(1.d0-z))*(corrabx)


C
C     Now the more nasty once...
C
C     JgQ - contribution
C
C     Gluon a
      muqqx = mt**2/2.d0/dotp(kqx,kax)*z
      muqq1 = mt**2/2.d0/dotp(kq1,ka1)
      Kcontax = Kcontax 
     &     - (
     &     ((1.d0-z)/2.d0/(1.d0-z+muqqx)**2 
     &     - 2.d0/(1.d0-z)*(1.d0 + dlog(1.d0-z+muqqx)))
     &     * corraqx
     &     + 2.d0/(1.d0 -z)*(dlog(2.d0+muqqx-z)*corraqx)
     &     )
      muqqx = mt**2/2.d0/dotp(kqbx,kax)*z
      muqq1 = mt**2/2.d0/dotp(kqb1,ka1)
      Kcontax = Kcontax 
     &     - ( 
     &     ((1.d0-z)/2.d0/(1.d0-z+muqqx)**2 
     &     - 2.d0/(1.d0-z)*(1.d0 + dlog(1.d0-z+muqqx)))
     &     * corraqbx
     &     + 2.d0/(1.d0 -z)*(dlog(2.d0+muqqx-z)*corraqbx)
     &     )

C     Gluon b
      muqqx = mt**2/2.d0/dotp(kqx,kbx)*z
      muqq1 = mt**2/2.d0/dotp(kq1,kb1)
      Kcontbx = Kcontbx 
     &     - (
     &     ((1.d0-z)/2.d0/(1.d0-z+muqqx)**2 
     &     - 2.d0/(1.d0-z)*(1.d0 + dlog(1.d0-z+muqqx)))
     &     * corrbqx
     &     + 2.d0/(1.d0 -z)*(dlog(2.d0+muqqx-z)*corrbqx)
     &     )
      muqqx = mt**2/2.d0/dotp(kqbx,kbx)*z
      muqq1 = mt**2/2.d0/dotp(kqb1,kb1)
      Kcontbx = Kcontbx 
     &     - ( 
     &     ((1.d0-z)/2.d0/(1.d0-z+muqqx)**2 
     &     - 2.d0/(1.d0-z)*(1.d0 + dlog(1.d0-z+muqqx)))
     &     * corrbqbx
     &     + 2.d0/(1.d0 -z)*(dlog(2.d0+muqqx-z)*corrbqbx)
     &     )

      Kcontax = alphas/2.d0/pi * Kcontax
      Kcontbx = alphas/2.d0/pi * Kcontbx

      endif ! if (jacobx.ne.0.d0) ...
C
C     Endpoint contribution:
C
      Apa1 = 0.d0
      Apb1 = 0.d0
      Kconta1 = 0d0
      Kcontb1 = 0d0
      
      if (jacob1.ne.0.d0) then
      APa1 = Pggnreg(z) * ( 
     &     - corra11*dlog(mu_fac**2/2.d0/dotp(kg1,ka1))
     &     - corraq1  * dlog(mu_fac**2/2.d0/dotp(kq1,ka1))
     &     - corraqb1 * dlog(mu_fac**2/2.d0/dotp(kqb1,ka1)) 
     &     - corrab1 * dlog(mu_fac**2/2.d0/dotp(kb1,ka1)) 
     &     )
C
      APb1 = Pggnreg(z) * ( 
     &     - corrb11*dlog(mu_fac**2/2.d0/dotp(kg1,kb1))
     &     - corrbq1  * dlog(mu_fac**2/2.d0/dotp(kq1,kb1))
     &     - corrbqb1 * dlog(mu_fac**2/2.d0/dotp(kqb1,kb1)) 
     &     - corrab1 * dlog(mu_fac**2/2.d0/dotp(ka1,kb1)) 
     &     )

      APa1 = alphas/2.d0/pi/CA*APa1
      APb1 = alphas/2.d0/pi/CA*APb1
C
C     The lower boundary of the z integration is not 0, we get
C     an additional contribution evaluating the +-prescription.
C
****************      APstuff = Apstuff + ...


C
C     Now the simple plus prescriptions
C
      Kconta1 = 
     &     + Kbargg_p(z) * ( - square1)
C     Kgg_g
     &     - Kgg_g_p(z,2.d0*dotp(kax,kgx)/z,0.d0) * (-corra11)
C     Kgg_q
     &     - (
     &     2.d0*dlog(1.d0-z)/(1.d0-z) * (-corraq1)
     &     + 2.d0/(1.d0-z) *( 
     &     -dlog((2.d0*dotp(kq1,ka1))/
     &     (2.d0*dotp(kq1,ka1)+mt**2))*corraq1)
     &     )
C     Kgg_qb
     &     - (
     &     2.d0*dlog(1.d0-z)/(1.d0-z) * (-corraqb1)
     &     + 2.d0/(1.d0-z) *( 
     &     -dlog((2.d0*dotp(kqb1,ka1))/
     &     (2.d0*dotp(kqb1,ka1)+mt**2))*corraqb1)
     &     )
C     From the inclusion of b
     &     - 2.d0*(dlog(1.d0-z)/(1.d0-z))*(-corrab1)

C
      Kcontb1 = 
     &     + Kbargg_p(z) * ( - square1 )
C     Kgg_g
     &     - Kgg_g_p(z,2.d0*dotp(kbx,kgx)/z,0.d0) * (-corrb11)
C     Kgg_q
     &     - (
     &     2.d0*dlog(1.d0-z)/(1.d0-z) * (-corrbq1)
     &     + 2.d0/(1.d0-z) *( 
     &     -dlog((2.d0*dotp(kq1,kb1))/
     &     (2.d0*dotp(kq1,kb1)+mt**2))*corrbq1)
     &     )
C     Kgg_qb
     &     - (
     &     2.d0*dlog(1.d0-z)/(1.d0-z) * (-corrbqb1)
     &     + 2.d0/(1.d0-z) *( 
     &     -dlog((2.d0*dotp(kqb1,kb1))/
     &     (2.d0*dotp(kqb1,kb1)+mt**2))*corrbqb1)
     &     )
C     From the inclusion of a
     &     - 2.d0*(dlog(1.d0-z)/(1.d0-z))*(-corrab1)


C     Now the more nasty once...
C
C     JgQ - contribution
C
C     Gluon a
      muqqx = mt**2/2.d0/dotp(kqx,kax)*z
      muqq1 = mt**2/2.d0/dotp(kq1,ka1)
      Kconta1 = Kconta1 
     &     - (
     &     - ((1.d0-z)/2.d0/(1.d0-z+muqq1)**2 
     &     - 2.d0/(1.d0-z)*(1.d0 + dlog(1.d0-z+muqq1)))
     &     * corraq1
     &     + 2.d0/(1.d0 -z)*(
     &     - dlog(1.d0+muqq1)*corraq1)
     &     )
      muqqx = mt**2/2.d0/dotp(kqbx,kax)*z
      muqq1 = mt**2/2.d0/dotp(kqb1,ka1)
      Kconta1 = Kconta1 
     &     - ( 
     &     -((1.d0-z)/2.d0/(1.d0-z+muqq1)**2 
     &     - 2.d0/(1.d0-z)*(1.d0 + dlog(1.d0-z+muqq1)))
     &     * corraqb1
     &     + 2.d0/(1.d0 -z)*(
     &     - dlog(1.d0+muqq1)*corraqb1)
     &     )

C     Gluon b
      muqqx = mt**2/2.d0/dotp(kqx,kbx)*z
      muqq1 = mt**2/2.d0/dotp(kq1,kb1)
      Kcontb1 = Kcontb1 
     &     - (
     &     - ((1.d0-z)/2.d0/(1.d0-z+muqq1)**2 
     &     - 2.d0/(1.d0-z)*(1.d0 + dlog(1.d0-z+muqq1)))
     &     * corrbq1
     &     + 2.d0/(1.d0 -z)*(
     &     - dlog(1.d0+muqq1)*corrbq1)
     &     )
      muqqx = mt**2/2.d0/dotp(kqbx,kbx)*z
      muqq1 = mt**2/2.d0/dotp(kqb1,kb1)
      Kcontb1 = Kcontb1 
     &     - ( 
     &     -((1.d0-z)/2.d0/(1.d0-z+muqq1)**2 
     &     - 2.d0/(1.d0-z)*(1.d0 + dlog(1.d0-z+muqq1)))
     &     * corrbqb1
     &     + 2.d0/(1.d0 -z)*(
     &     - dlog(1.d0+muqq1)*corrbqb1)
     &     )


C
C
      Kconta1 = alphas/2.d0/pi * Kconta1
      Kcontb1 = alphas/2.d0/pi * Kcontb1

      endif ! if (jacobx.ne.0.d0) ...

C
C     Put the two contributions together:
C
      convax = APax + Kcontax
      convbx = APbx + Kcontbx
      conva1 = APa1 + Kconta1
      convb1 = APb1 + Kcontb1
      
      conva = convax + conva1
      convb = convbx + convb1

      conv = conva ! + convb

*      write(*,*) Apax, APa1
*      conv =    APax + Apa1   

      conv =   convb
C
C     Now the delta functions
C     First from the AP-kernel
C
C     Gluon a
      deltaa = g_g  * (
     &     + corra11  * dlog(mu_fac**2/2.d0/dotp(kg1,ka1))
     &     + corraq1  * dlog(mu_fac**2/2.d0/dotp(kq1,ka1))
     &     + corraqb1 * dlog(mu_fac**2/2.d0/dotp(kqb1,ka1))
     &     + corrab1  * dlog(mu_fac**2/2.d0/dotp(kb1,ka1)) 
     &     ) / CA

C     Just for the comparison with Stefan
      APdeltaa = alphas/2.d0/pi * deltaa

      deltaa = deltaa
C     Delta term in Kbar:
     &     - square1 * (g_g + K_g - 5.d0/6.d0 * pi**2 * CA)
C     Next term comes from \cal K:
C     Kgg_g:
     &     - corra11 * (-g_g/CA)
C     Kgg_q:
      sja = 2.d0*dotp(kq1,ka1)
      mj2 = mt**2
      deltaa = deltaa - corraq1 * (-g_q/Cf 
     &     + mj2/sja*dlog(mj2/(sja+mj2))
     &     + 1.d0/2.d0 * mj2/(sja+mj2))
C     Kgg_qb:
      sja = 2.d0*dotp(kqb1,ka1)
      mj2 = mt**2
      deltaa = deltaa - corraqb1 * (-g_q/Cf 
     &     + mj2/sja*dlog(mj2/(sja+mj2))
     &     + 1.d0/2.d0 * mj2/(sja+mj2))
C     The contribution from the third sum:
C     j=quark
      mj = mt
      mj2 = mt**2
      sja = 2.d0*dotp(kq1,ka1)
      deltaa = deltaa 
     &     - 1.d0/CA *  corraq1 * g_g *
     &     ( dlog((sja-2.d0*mj*dsqrt(sja+mj2)+2.d0*mj2)/sja)
     &     + 2.d0*mj/(dsqrt(sja+mj2)+mj) )
C     j=anti-quark
      mj = mt
      mj2 = mt**2
      sja = 2.d0*dotp(kqb1,ka1)
      deltaa = deltaa 
     &     - 1.d0/CA *  corraqb1 * g_g *
     &     ( dlog((sja-2.d0*mj*dsqrt(sja+mj2)+2.d0*mj2)/sja)
     &     + 2.d0*mj/(dsqrt(sja+mj2)+mj) )
C     
C     And to the end, the term which appears explicitly in eq. 6.68
      deltaa = deltaa + corrab1 * pi**2/3.d0
C     
C     Gluon b
      deltab = g_g  * (
     &     + corrb11  * dlog(mu_fac**2/2.d0/dotp(kg1,kb1))
     &     + corrbq1  * dlog(mu_fac**2/2.d0/dotp(kq1,kb1))
     &     + corrbqb1 * dlog(mu_fac**2/2.d0/dotp(kqb1,kb1))
     &     + corrab1  * dlog(mu_fac**2/2.d0/dotp(kb1,ka1)) 
     &     ) / CA

      deltab = deltab
C     Delta term in Kbar:
     &     - square1 * (g_g + K_g - 5.d0/6.d0 * pi**2 * CA)
C     Next term comes from \cal K:
C     Kgg_g:
     &     - corrb11 * (-g_g/CA)
C     Kgg_q:
      sja = 2.d0*dotp(kq1,kb1)
      mj2 = mt**2
      deltab = deltab - corrbq1 * (-g_q/Cf 
     &     + mj2/sja*dlog(mj2/(sja+mj2))
     &     + 1.d0/2.d0 * mj2/(sja+mj2))
C     Kgg_qb:
      sja = 2.d0*dotp(kqb1,kb1)
      mj2 = mt**2
      deltab = deltab - corrbqb1 * (-g_q/Cf 
     &     + mj2/sja*dlog(mj2/(sja+mj2))
     &     + 1.d0/2.d0 * mj2/(sja+mj2))
C     The contribution from the third sum:
C     j=quark
      mj = mt
      mj2 = mt**2
      sja = 2.d0*dotp(kq1,kb1)
      deltab = deltab 
     &     - 1.d0/CA *  corrbq1 * g_g *
     &     ( dlog((sja-2.d0*mj*dsqrt(sja+mj2)+2.d0*mj2)/sja)
     &     + 2.d0*mj/(dsqrt(sja+mj2)+mj) )
C     j=anti-quark
      mj = mt
      mj2 = mt**2
      sja = 2.d0*dotp(kqb1,kb1)
      deltab = deltab 
     &     - 1.d0/CA *  corrbqb1 * g_g *
     &     ( dlog((sja-2.d0*mj*dsqrt(sja+mj2)+2.d0*mj2)/sja)
     &     + 2.d0*mj/(dsqrt(sja+mj2)+mj) )
C     
C     And to the end, the term which appears explicitly in eq. 6.68
      deltab = deltab + corrab1 * pi**2/3.d0

C
      deltaa = alphas/2.d0/pi * deltaa
      deltab = alphas/2.d0/pi * deltab



#ifdef COMPARE
      write(*,*) '----------------------------------------------'
      write(*,*) '------------- Stefans Notation ---------------'
      write(*,*) '---------- the subtracted term is ------------'
      write(*,*) '--------- included in the endpoint -----------'
      write(*,*) '----------------------------------------------'
      Write(*,*) '----- when comparing with Stefan we should ---'
      write(*,*) '----- the color ordered MEs to 1 -------------'
      write(*,*) '--------- The P contrib.: --------------------'
      write(*,*) 'APa = ', (Apax)/(alphas/2.d0/pi)
      write(*,*) 'APb = ', (Apbx)/(alphas/2.d0/pi)
      write(*,*) '--------- The P endpoint contrib.: -----------'
      write(*,*) 'APa = ', (APa1 + APdeltaa)/(alphas/2.d0/pi)
*      write(*,*) 'APb = ', (APb1 + APdeltab)/(alphas/2.d0/pi)
      write(*,*) '--------- The K contrib.: --------------------'
      write(*,*) (Kcontax)/(alphas/2.d0/pi)
      write(*,*) '--------- The K endpoint contrib.: -----------'
      write(*,*) (Kconta1+deltaa-APdeltaa)/(alphas/2.d0/pi)
      write(*,*) '--------- The sum of all: --------------------'
      write(*,*) '-- Here the color ordered MEs are included: --'
      write(*,*) '----- Endpoint:-------------------------------'
      write(*,*) (conva1+deltaa)/64.d0
      write(*,*) '----- z<>1 -----------------------------------'
      write(*,*) (convax)/64.d0
      write(*,*) '----------------------------------------------'
      stop
#endif
      end

#ifdef COMPARE
C
C     Just for comparison with Stefan....
C
      function Kcal(x,sij,mj2)
      implicit none
      double precision kcal, x,sij,mj2
      Kcal = -( 2.d0/(1.d0-x)*dlog(1.d0-x)
     &     -2.d0/(1.d0-x)*dlog(2.d0-x)
     &     +( (1.d0-x)/2.d0/(1.d0-x+mj2/sij)**2 
     &     - 2.d0/(1.d0-x)*(1.d0+dlog(1.d0-x+mj2/sij))
     &     +2.d0/(1.d0-x)*dlog(2.d0+mj2/sij-x))
     &     + 2.d0/(1.d0-x)*dlog((2.d0-x)*sij/((2.d0-x)*sij+mj2))
     &     )
     &     -2.d0*mj2/x/sij*log(mj2/((1.d0-x)*sij+mj2))
     &     -2.d0*((1.d0-x)/x-1.d0+x*(1.d0-x)) 
     &     * dlog((1.d0-x)*sij/((1.d0-x)*sij+mj2));
      end

      subroutine stefan(iset,z,ka,kb,kq,kqb,kg)
      implicit none
      integer iset
      integer mu
      double precision z,R(3,3),phi
      double precision ka(4),kb(4),kq(4),kqb(4),kg(4)
C code produced from Stefans input by Converter.C
      if (iset.eq.1) then
      if (z.eq.0.7d0) then
         kq(1) = -91.12817618641557d0
         kq(2) = -224.4144988454333d0
         kq(3) = -93.98674376768055d0
         kq(4) = 312.6910932392566d0
         kg(1) = 284.3721344703247d0
         kg(2) = 569.0124665251295d0
         kg(3) = 598.1478105631707d0
         kg(4) = 873.16865564777d0
         ka(1) = 0d0
         ka(2) = 0d0
         ka(3) = 686d0
         ka(4) = -686d0
         kb(1) = 0d0
         kb(2) = 0d0
         kb(3) = -980d0
         kb(4) = -980d0
         kqb(1) = -193.2439582839092d0
         kqb(2) = -344.5979676796965d0
         kqb(3) = -210.1610667954903d0
         kqb(4) = 480.1402511129734d0
      else
         kq(1) = -110.5466594152423d0
         kq(2) = -272.2349355588941d0
         kq(3) = -183.8387820821969d0
         kq(4) = 387.8209920106065d0
         kg(1) = 344.968930708861d0
         kg(2) = 690.2632091671036d0
         kg(3) = 547.2726520216188d0
         kg(4) = 946.0307694440274d0
         ka(1) = 0d0
         ka(2) = 0d0
         ka(3) = 980d0
         ka(4) = -980d0
         kb(1) = 0d0
         kb(2) = 0d0
         kb(3) = -980d0
         kb(4) = -980d0
         kqb(1) = -234.4222712936187d0
         kqb(2) = -418.0282736082098d0
         kqb(3) = -363.4338699394223d0
         kqb(4) = 626.148238545366d0
      endif
      endif
      if (iset.eq.2) then
      if (z.eq.0.3d0) then
         kq(1) = -55.75700837376461d0
         kq(2) = -137.3085868164699d0
         kq(3) = 47.57342474484788d0
         kq(4) = 233.4560402694306d0
         kg(1) = 173.9938199848536d0
         kg(2) = 348.1517373497873d0
         kg(3) = 632.481058359739d0
         kg(4) = 742.6410780427193d0
         ka(1) = 0d0
         ka(2) = 0d0
         ka(3) = 294d0
         ka(4) = -294d0
         kb(1) = 0d0
         kb(2) = 0d0
         kb(3) = -980d0
         kb(4) = -980d0
         kqb(1) = -118.236811611089d0
         kqb(2) = -210.8431505333176d0
         kqb(3) = 5.945516895412899d0
         kqb(4) = 297.90288168785d0
      else
         kq(1) = -110.5466594152423d0
         kq(2) = -272.2349355588941d0
         kq(3) = -183.8387820821969d0
         kq(4) = 387.8209920106065d0
         kg(1) = 344.968930708861d0
         kg(2) = 690.2632091671036d0
         kg(3) = 547.2726520216188d0
         kg(4) = 946.0307694440274d0
         ka(1) = 0d0
         ka(2) = 0d0
         ka(3) = 980d0
         ka(4) = -980d0
         kb(1) = 0d0
         kb(2) = 0d0
         kb(3) = -980d0
         kb(4) = -980d0
         kqb(1) = -234.4222712936187d0
         kqb(2) = -418.0282736082098d0
         kqb(3) = -363.4338699394223d0
         kqb(4) = 626.148238545366d0
      endif
      endif

      if (iset.eq.11) then
      if (z.eq.0.7d0) then
      kq(1) = -193.2439582839092d0
      kq(2) = -344.5979676796965d0
      kq(3) = -210.1610667954903d0
      kq(4) = 480.1402511129734d0
      kg(1) = 284.3721344703247d0
      kg(2) = 569.0124665251295d0
      kg(3) = 598.1478105631707d0
      kg(4) = 873.16865564777d0
      ka(1) = 0d0
      ka(2) = 0d0
      ka(3) = 686d0
      ka(4) = -686d0
      kb(1) = 0d0
      kb(2) = 0d0
      kb(3) = -980d0
      kb(4) = -980d0
      kqb(1) = -91.12817618641557d0
      kqb(2) = -224.4144988454333d0
      kqb(3) = -93.98674376768055d0
      kqb(4) = 312.6910932392566d0
      else
      kq(1) = -234.4222712936187d0
      kq(2) = -418.0282736082098d0
      kq(3) = -363.4338699394223d0
      kq(4) = 626.148238545366d0
      kg(1) = 344.968930708861d0
      kg(2) = 690.2632091671036d0
      kg(3) = 547.2726520216188d0
      kg(4) = 946.0307694440274d0
      ka(1) = 0d0
      ka(2) = 0d0
      ka(3) = 980d0
      ka(4) = -980d0
      kb(1) = 0d0
      kb(2) = 0d0
      kb(3) = -980d0
      kb(4) = -980d0
      kqb(1) = -110.5466594152423d0
      kqb(2) = -272.2349355588941d0
      kqb(3) = -183.8387820821969d0
      kqb(4) = 387.8209920106065d0
      endif
      endif

      do mu=1,4
         ka(mu)=-ka(mu)
         kb(mu)=-kb(mu)
      enddo
      phi = 0.8d0
      R(1,1) = dcos(phi)
      R(3,3) = dcos(phi)
      R(1,3) = dsin(phi)
      R(3,1) = -dsin(phi)
      R(2,2) = 1.d0
      R(1,2) = 0.d0
      R(2,1) = 0.d0
      R(2,3) = 0.d0
      R(3,2) = 0.d0
      
      write(*,*) kq
      call rotate(R,ka)
      call rotate(R,kb)
      call rotate(R,kg)
      call rotate(R,kq)
      call rotate(R,kqb)  

      end
#endif
