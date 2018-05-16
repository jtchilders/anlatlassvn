C     $Modified: Wed Jan  3 14:38:17 2007 by puwer $
**#define VERBOSE
**#define DEBUG
      function subtest(p1,p2,kq,kqb,p3,p4)
      implicit none
      double precision subtest,p(4,-2:4)
      double precision p1(4),p2(4),kq(4),kqb(4),p3(4),p4(4)
      integer mu
      double precision subtractions,D
      do mu = 1,4
          p(mu,-2) = p2(mu)
          p(mu,-1) = p1(mu)
          p(mu,1) = kq(mu)
          p(mu,2) = kqb(mu)
          p(mu,3) = p3(mu)
          p(mu,4) = p4(mu)
      enddo
**      subtest =  D(1,3,2,p)+D(1,3,4,p)+ D(1,4,2,p)+ D(1,4,3,p)
      subtest = subtractions(p)
      end

      function subtractions(p)
      implicit none
      double precision subtractions
      double precision p(4,-2:4)
      double precision D
      double precision t1,t2,t3
#include "common.F"
      subtractions = 0.d0
#ifdef DEBUG
      write(*,*) 'Subtractions: D1...'
#endif
C     Final state emitter, final state spectator (F-F)
      subtractions = 
     >       D(1,3,2,p) + D(1,3,4,p) + D(1,4,2,p) + D(1,4,3,p)
     >     + D(2,3,1,p) + D(2,3,4,p) + D(2,4,1,p) + D(2,4,3,p)
     >     + D(3,4,1,p) + D(3,4,2,p)
C     Final state emitter, initial state spectator (F-I)      
#ifdef DEBUG
      write(*,*) 'Subtractions: D2...'
#endif
      subtractions = subtractions 
     >     + D(1,3,-2,p) + D(1,3,-1,p) + D(1,4,-2,p) + D(1,4,-1,p)
     >     + D(2,3,-2,p) + D(2,3,-1,p) + D(2,4,-2,p) + D(2,4,-1,p)
     >     + D(3,4,-2,p) + D(3,4,-1,p)
C     Initial state emitter, final state spectator (I-F) 
#ifdef DEBUG
      write(*,*) 'Subtractions: D3...'
#endif
      subtractions = subtractions     
     >     + D(-2,3,1,p) + D(-2,3,2,p) + D(-2,3,4,p) + D(-2,4,1,p)
     >     + D(-2,4,2,p) + D(-2,4,3,p) + D(-1,3,1,p) + D(-1,3,2,p)
     >     + D(-1,3,4,p) + D(-1,4,1,p) + D(-1,4,2,p) + D(-1,4,3,p)
C     Initial state emitter, initial state spectator (I-I) 
#ifdef DEBUG
      write(*,*) 'Subtractions: D4...'
#endif
      subtractions = subtractions + 
     >     D(-2,3,-1,p) + D(-2,4,-1,p) + D(-1,3,-2,p) + D(-1,4,-2,p)
      subtractions = subtractions * (-1.d0)
      end
C*************************************************************************
      function D(i,j,k,p)
      implicit none
      double precision D
      integer i,j,k
      double precision p(4,-2:4)
      double precision  Dij_k, Dij_a, Dai_j,Dai_b
#ifdef DEBUG
      if ((i.lt.0) .and.(j.lt.0)) then
         call error('D: error both coll. part. in the in. state')
         stop
      endif
#endif
      if ( i*j .gt. 0 ) then
         if (k.gt.0) then 
C--------------------------------------------------------------------
C     Final-state emitter and final-state spectator:
C--------------------------------------------------------------------
            D = Dij_k(i,j,k,p)
#ifdef DEBUG
            write(*,*) 'D:', i,j,k, D
#endif
            return 
         else
C--------------------------------------------------------------------
C     Final-state emitter and initial-state spectator:
C--------------------------------------------------------------------
            D = Dij_a(i,j,k,p)
#ifdef DEBUG
            write(*,*) 'D:', i,j,k, D
#endif
            return 
         endif
      else
         if (k.gt.0) then
C--------------------------------------------------------------------
C     Initial-state emitter and final-state spectator:
C--------------------------------------------------------------------
            D = Dai_j(i,j,k,p)
#ifdef DEBUG
            write(*,*) 'D:', i,j,k, D
#endif
            return 
         else
C--------------------------------------------------------------------
C     Initial-state emitter and initial-state spectator:
C--------------------------------------------------------------------
            D = Dai_b(i,j,k,p)
#ifdef DEBUG
            write(*,*) 'D:', i,j,k, D
#endif
            return 
         endif
      endif
      end
C*************************************************************************
C*************************************************************************
      function Dai_j(i,j,k,p)
      implicit none 
      integer i,j,k
      double precision Dai_j, p(4,-2:4),ka(4),ki(4),kj(4),kk(4)
      double precision pn(4,-2:3),ktai(4),ktj(4)
      double precision nvec(4),n1(4),n2(4), a, b
      double precision xij_a,zti,ztj,Tijq
      integer ihelp 
      double precision obs3
      double precision Dggg
#include "subtractions.h"
#include "common.F"
#include "constants.F"
      integer plist(-2:4),temitter, tspectator,rho
      data plist /gluon,gluon,-4711,gluon,gluon,quark,antiquark/
      integer findEmitter
      double precision dotp
      Dai_j = 0.d0
C
C     ka || ki, kj is the final state spectator
C
      if (i.gt.j) then
         write(*,*) 'Dai_j: the first parton should be the',
     .        ' initial state parton'
         stop
      endif
      do rho=1,4
         ka(rho) = p(rho,i)
         ki(rho) = p(rho,j)
         kj(rho) = p(rho,k)
      enddo

      xij_a = ( dotp(ka,ki) + dotp(ka,kj)- dotp(ki,kj))
     .     / ( dotp(ka,ki) + dotp(ka,kj) )
      zti = dotp(ka,ki) / (dotp(ka,ki) + dotp(ka,kj))
      ztj = 1.d0 - zti
      do rho = 1,4
         ktai(rho) = xij_a * ka(rho)
         ktj(rho) = ki(rho)+ kj(rho) - (1.d0 - xij_a) * ka(rho) 
      enddo

      if (i.eq.-1) then
         do rho=1,4
            pn(rho,-1) = ktai(rho) 
            pn(rho,-2) = p(rho,-2)
         enddo
      else
         do rho=1,4
            pn(rho,-1) = p(rho,-1)
            pn(rho,-2) = ktai(rho) 
         enddo
      endif
C
C     temitter = findEmitter(4,i,j,plist)
C     for the moment we only study the case that the emitter is 
C     a gluon
      temitter = gluon

      if (temitter.eq.gluon) then
C     
C     g --> gg splitting
C     
         Tijq = CA
         a =  16.d0 * pi * alphas * Ca * (1.d0/(2.d0-xij_a-ztj)
     .        -1.d0 + xij_a* (1.d0 - xij_a))
         b = 16.d0 * pi* alphas * CA *
     .        (1.d0-xij_a)/xij_a*zti*ztj/dotp(ki,kj)

         do rho = 1,4
            nvec(rho) = ki(rho)/zti - kj(rho)/ztj
         enddo
         if ((k.eq.3).or.(k.eq.4)) then
            do rho=1,4
               pn(rho,1) = p(rho,1)
               pn(rho,2) = p(rho,2)
               pn(rho,3) = ktj(rho) 
            enddo
            Dai_j = obs3(pn)
            if (Dai_j.eq.0.d0) return
            if (i.eq.-1) 
     .           Dai_j = Dai_j * Dggg(i,pn,Cg1g3,a,b,nvec)
            if (i.eq.-2) 
     .           Dai_j = Dai_j * Dggg(i,pn,Cg2g3,a,b,nvec)
         else
            if (j.eq.3) then
               ihelp = 4
            else
               ihelp = 3
            endif
            if (k.eq.1) then
               do rho=1,4
                  pn(rho,1) = ktj(rho) 
                  pn(rho,2) = p(rho,2)
                  pn(rho,3) = p(rho,ihelp)
               enddo
               Dai_j = obs3(pn)
               if (Dai_j.eq.0.d0) return
               if (i.eq.-1) 
     .              Dai_j = Dai_j * Dggg(i,pn,Cg1q,a,b,nvec)
               if (i.eq.-2) 
     .              Dai_j = Dai_j * Dggg(i,pn,Cg2q,a,b,nvec)
            else
               if (k.eq.2) then 
                  do rho=1,4
                     pn(rho,1) = p(rho,1)
                     pn(rho,2) = ktj(rho) 
                     pn(rho,3) = p(rho,ihelp)
                  enddo
                  Dai_j = obs3(pn)
                  if (Dai_j.eq.0.d0) return
                  if (i.eq.-1) 
     .                 Dai_j = Dai_j * Dggg(i,pn,Cg1qb,a,b,nvec)
                  if (i.eq.-2) 
     .                 Dai_j = Dai_j * Dggg(i,pn,Cg2qb,a,b,nvec)
               else
                  write(*,*) 'Dai_j: inconsistency'
                  stop
               endif
            endif
         endif
      endif

      if (temitter.eq.quark) then
         stop
      endif
      if (temitter.eq.antiquark) then
         stop
      endif
      Dai_j = -Dai_j/2.d0/dotp(ka,ki)/xij_a/Tijq
#ifdef VERBOSE
*      write(*,*) 'Prefactor = ',-1.d0/2.d0/dotp(ka,ki)/xij_a/Tijq
      write(*,'(A6,I3,I3,I3,A4,$)') 'Dai_j(',i,j,k,') = '
      write(*,'(E25.15)') Dai_j
#endif
      end
C*************************************************************************
C*************************************************************************
      function Dai_b(i,j,k,p)
      implicit none 
      double precision Dai_b, p(4,-2:4),ka(4),ki(4),kb(4)
      double precision ktai(4),kv(4),ktv(4),pn(4,-2:3),nvec(4)
      double precision k1(4),k2(4),k3(4)
      double precision Dggg, obs3
      integer i,j,k
      integer rho
      double precision xi_ab,a,b
      double precision kvq,ktvq,ksumq
#include "subtractions.h"
#include "common.F"
#include "constants.F"
      integer plist(-2:4),temitter, tspectator
      data plist /gluon,gluon,-4711,gluon,gluon,quark,antiquark/
      integer findEmitter
      double precision dotp

      Dai_b = 0.d0
      
      if (i.gt.0) then 
         write(*,*) 'Dai_b: i>j, i should be in the initial state,'
         stop
      endif

      do rho=1,4
         ka(rho) = p(rho,i)
         ki(rho) = p(rho,j)
         kb(rho) = p(rho,k)
         k1(rho) = p(rho,1)
         k2(rho) = p(rho,2)
      enddo

      if (j.eq.4) then
         do rho=1,4
            k3(rho) = p(rho,3)
         enddo
      else
         do rho=1,4
            k3(rho) = p(rho,4)
         enddo
      endif

      xi_ab = (dotp(ka,kb)-dotp(ki,ka)-dotp(ki,kb))/dotp(ka,kb)
      do rho =1,4
         ktai(rho) = xi_ab * ka(rho)
         kv(rho) = ka(rho)+kb(rho) - ki(rho)
         ktv(rho) = ktai(rho)+kb(rho) 
      enddo

      kvq = dotp(kv,kv)
      ktvq = dotp(ktv,ktv)
      ksumq = kvq + ktvq + 2.d0 * dotp(kv,ktv)

      do rho =1,4
         pn(rho,1) = k1(rho)
     .        - 2.d0 * (dotp(k1,kv)+dotp(k1,ktv))/ksumq 
     .        * (kv(rho) + ktv(rho))
     .        + 2.d0 * dotp(k1,kv)/kvq*ktv(rho)
         pn(rho,2) = k2(rho)
     .        - 2.d0 * (dotp(k2,kv)+dotp(k2,ktv))/ksumq 
     .        * (kv(rho) + ktv(rho))
     .        + 2.d0 * dotp(k2,kv)/kvq*ktv(rho)
         pn(rho,3) = k3(rho)
     .        - 2.d0 * (dotp(k3,kv)+dotp(k3,ktv))/ksumq 
     .        * (kv(rho) + ktv(rho))
     .        + 2.d0 * dotp(k3,kv)/kvq*ktv(rho)
      enddo

      if (i.eq.-1) then
         do rho =1,4
            pn(rho,-1) = ktai(rho)
            pn(rho,-2) = p(rho,-2)
         enddo
      else
         do rho =1,4
            pn(rho,-1) = p(rho,-1)
            pn(rho,-2) = ktai(rho)
         enddo
      endif

      a = 16.d0 * pi * alphas * CA * ( xi_ab / (1.d0 -xi_ab)
     .     + xi_ab * (1.d0 - xi_ab) )
      b = 16.d0 * pi * alphas * CA * (1.d0 -xi_ab)/xi_ab
     .     * dotp(ka,kb)/dotp(ki,ka)/dotp(ki,kb)
      do rho = 1,4
         nvec(rho) = ki(rho)-dotp(ki,ka)/dotp(kb,ka)*kb(rho)
      enddo
      Dai_b = obs3(pn)
      if (Dai_b.eq.0.d0) return
      Dai_b = Dai_b * Dggg(i,pn,Cg1g2,a,b,nvec)
      Dai_b = -1.d0/2.d0/dotp(ka,ki)/xi_ab/CA*Dai_b
#ifdef VERBOSE
*      write(*,*) 'Prefactor = ',-1.d0/2.d0/dotp(ka,ki)/xi_ab/CA
      write(*,'(A6,I3,I3,I3,A4,$)') 'Dai_b(',i,j,k,') = '
      write(*,'(E25.15)') Dai_b
#endif
      end
C*************************************************************************
C*************************************************************************
      function Dij_a(i,j,k,p)
      implicit none  
      double precision massq(-2:4)
      double precision Dij_a, p(4,-2:4),Q(4),ki(4),kj(4),ka(4)
      double precision pn(4,-2:3),kta(4),ktij(4)
      double precision help
      integer i,j,k,rho
      double precision mui,muj,muij,miq, mjq, mijq,xij_a,zti,ztj
      double precision mqq,ztim,ztjm
      double precision a,b,nvec(4),n1(4),n2(4)
      double precision Tijq
      data massq/0.d0,0.d0,-1.d0,1.d0,1.d0,0.d0,0.d0/
#include "subtractions.h"
#include "common.F"
#include "constants.F"
      integer plist(-2:4),temitter, tspectator
      data plist /gluon,gluon,-4711,quark,antiquark,gluon,gluon/
      integer findEmitter
      double precision V
      double precision dotp,dggg,dqgq,Dgqqb ,obs3   

#ifdef VERBOSE
      write(*,*) 'Dij_a:'
#endif 
C
      Dij_a = 0.d0

C     Set the masses:
      do rho =1,4
         nvec(rho) = p(rho,1)
      enddo
      
      massq(1) = dotp(nvec,nvec)
      massq(2) = massq(1)

      miq = massq(i)
      mjq = massq(j)

      do rho = 1,4
         ki(rho) = p(rho,i)
         kj(rho) = p(rho,j)
         ka(rho) = p(rho,k)
      enddo
      
      do rho = 1,4
         Q(rho) = ki(rho) + kj(rho) - ka(rho)
      enddo

      temitter = findEmitter(4,i,j,plist)
      tspectator = plist(k)
#ifdef VERBOSE
      write(*,*) 'Dij_a: i = ',i 
      write(*,*) 'Dij_a: j = ',j 
      write(*,*) 'Dij_a: k = ',k 
      if (temitter.eq.gluon) write(*,*) 'Dij_a: emit. = gluon'
      if (temitter.eq.quark) write(*,*) 'Dij_a: emit. = quark'
      if (temitter.eq.antiquark) write(*,*) 'Dij_a: emit. = aquark'
      if (tspectator.eq.gluon) write(*,*) 'Dij_a: spec. = gluon'
      if (tspectator.eq.quark) write(*,*) 'Dij_a: spec. = quark'
      if (tspectator.eq.antiquark) write(*,*) 'Dij_a: spec. = aquark'
#endif
      if (temitter.eq.gluon) then
         mijq = 0.d0
         Tijq = CA
      else
C     
C     if the emitter is a Quark or an anti-quark one of the collinear
C     partons is massless
         mijq = massq(i)+massq(j)
         Tijq = CF
      endif
      
      xij_a = ( dotp(ka,ki) + dotp(ka,kj) - dotp(ki,kj)
     .     + 0.5d0 * (mijq - miq - mjq ) )/( dotp(ka,ki)+dotp(ka,kj) )
      zti = dotp(ka,ki)/(dotp(ka,ki)+dotp(ka,kj))
      ztj = 1.d0 - zti
      
      do rho = 1,4
         kta(rho) = xij_a * ka(rho)
         ktij(rho) = ki(rho) + kj(rho)
     .        - ( 1.d0 - xij_a ) * ka(rho)
      enddo

#ifdef VERBOSE      
      write(*,*) 'Dij_a: ktij.ktij = ', dotp(ktij,ktij)
      write(*,*) 'Dij_a: kta.kta = ',dotp(kta,kta)
      write(*,*) 'Dij_a: Q(mu) in old and new variables:'
      do rho=1,4
         write(*,*) Q(rho), ' = ', ktij(rho) - kta(rho)
      enddo
#endif
      help = 2.d0 * dotp(ktij,ka)
      mui  = dsqrt(miq) / dsqrt(help)
      muj  = dsqrt(mjq) / dsqrt(help)
      muij = dsqrt(mijq) / dsqrt(help)
      
      if (k.eq.-1) then
         do rho =1,4
            pn(rho,-1) = kta(rho)
            pn(rho,-2) = p(rho,-2)
         enddo
      else
         do rho =1,4
            pn(rho,-1) = p(rho,-1)
            pn(rho,-2) = kta(rho)
         enddo
      endif
      
      if ( temitter.eq.gluon) then 
         if ((i.eq.quark).or.(j.eq.quark)) then
C     g--> Q Qb
            stop
            if (massq(i) .ne.0.d0) then
C     massive quarks:
               a = 8.d0 * pi * alphas * TR 
               b = - 32.d0 * alphas * pi *TR / 
     .              (massq(i) + 2.d0*dotp(ki,kj) + massq(j))
               do rho = 1,4
                  n1(rho) = zti * ki(rho) - ztj * kj(rho)
               enddo
            else
C     massless quarks:
               a = 8.d0 * pi * alphas * TR 
               b = - 16.d0 * alphas * pi *TR / dotp(ki,kj)
               do rho = 1,4
                  n1(rho) = zti * ki(rho) - ztj * kj(rho)
               enddo
               stop
            endif
            stop
            Dij_a = obs3(pn)
            if (Dij_a.eq.0.d0) then
               return
            endif
            Dij_a = Dgqqb(2,pn,Cxy,a,b,nvec)
         else
            if ((plist(i).eq.gluon).and.(plist(j).eq.gluon)) then
C     g--> gg
               a = 16.d0 * pi * alphas * CA * 
     .              ( 1.d0/(1.d0-zti+(1.d0-xij_a))
     .              + 1.d0/(1.d0-ztj+(1.d0-xij_a)) 
     .              - 2.d0 )
               b = 16.d0 * pi * alphas * CA / dotp(ki,kj)
               do rho = 1,4 
                  nvec(rho) = zti * ki(rho) - ztj * kj(rho)
               enddo
               do rho=1,4 
                  pn(rho,1) = p(rho,1)
                  pn(rho,2) = p(rho,2)
                  pn(rho,3) = ktij(rho)
               enddo
            endif
            Dij_a = obs3(pn)
            if (Dij_a.eq.0.d0) return
            if (k.eq.-1) then
               Dij_a = Dij_a * Dggg(3,pn,Cg1g3,a,b,nvec)
            else
               Dij_a = Dij_a * Dggg(3,pn,Cg2g3,a,b,nvec)
            endif
         endif
      endif
      if (( temitter.eq.quark).or.( temitter.eq.antiquark)) then
#ifdef VERBOSE
      write(*,*) 'Dij_a: q --> q g, qb --> qb g'
#endif
          mqq = mijq
          if (plist(i).eq.gluon) then
             V = 8.d0 * pi * alphas * Cf * (2.d0/(2.d0 - xij_a - ztj)
     .            - 1.d0 - ztj - mqq/dotp(ki,kj))
          else
             V = 8.d0 * pi * alphas * Cf * (2.d0/(2.d0 - xij_a - zti)
     .            - 1.d0 - zti - mqq/dotp(ki,kj))
          endif
          if ((i.eq.3).or.(j.eq.3)) then
             do rho=1,4 
                pn(rho,3) = p(rho,4)
             enddo
          else   
             do rho=1,4 
                pn(rho,3) = p(rho,3)
             enddo
          endif
          if (temitter.eq.quark) then
             do rho=1,4 
                pn(rho,1) = ktij(rho)
                pn(rho,2) = p(rho,2)
             enddo
             Dij_a = obs3(pn)
             if (Dij_a.eq.0.d0) return

             if (k.eq.-1) then
                Dij_a = Dij_a * Dqgq(1,pn,Cg1q,V)
             else
                Dij_a = Dij_a * Dqgq(1,pn,Cg2q,V)
             endif
          else
             do rho=1,4 
                pn(rho,1) = p(rho,1)
                pn(rho,2) = ktij(rho)
             enddo
             Dij_a = obs3(pn)
             if (Dij_a.eq.0.d0) return

             if (k.eq.-1) then
                Dij_a = Dij_a * Dqgq(2,pn,Cg1qb,V)
             else
                Dij_a = Dij_a * Dqgq(2,pn,Cg2qb,V)
             endif
          endif
       endif
      Dij_a = -Dij_a / (miq+2.d0*dotp(ki,kj)+mjq-mijq)/xij_a/Tijq 
#ifdef VERBOSE
*      write(*,*) 'Prefactor = ',
*     &     -1.d0 / (miq+2.d0*dotp(ki,kj)+mjq-mijq)/xij_a/Tijq 
      write(*,'(A6,I3,I3,I3,A4,$)') 'Dij_a(',i,j,k,') = '
      write(*,'(E25.15)') Dij_a
#endif
      end
C*************************************************************************
C
C     Final state emitter, final state spectator
C
C*************************************************************************
      function Dij_k(i,j,k,p)
      implicit none 
      double precision Dij_k, p(4,-2:4),Q(4),ki(4),kj(4),kk(4)
      double precision massq(-2:4)
      integer i,j,k,ihelp,mu
      double precision ktij(4),ktk(4),nvec(4),n1(4),n2(4), pn(4,-2:3)
      double precision a,b
      double precision mui,muj,muij,muk
      double precision miq,mjq,mijq,mkq,mquark2
      double precision vtij_k, yij_k, zti,ztj,vij_k
      double precision vij_i,ym,yp,zp,zm,ztim,ztjm
      double precision VgQ
      double precision QQ,sqrtQ
      double precision lambda,x,y,z,rho 
      double precision Dggg, Dqgq, Dgqqb,obs3
      double precision Tijq
      data massq/0.d0,0.d0,-1.d0,1.d0,1.d0,0.d0,0.d0/
#include "subtractions.h"
#include "common.F"
#include "constants.F"
      integer plist(-2:4),temitter, tspectator
      data plist /gluon,gluon,-4711,quark,antiquark,gluon,gluon/
      integer findEmitter
      double precision dotp      
      lambda(x,y,z) = x**2+y**2+z**2-2.d0*x*y-2.d0*x*z-2.d0*y*z
      
      Dij_k = 0.d0

C     Set the masses:
      do mu =1,4
         nvec(mu) = p(mu,1)
      enddo
      massq(1) = dotp(nvec,nvec)
      massq(2) = massq(1)

C     check if the partons are all in the final state:
#ifdef VERBOSE
      if ((i.le.0).or.(j.le.0).or.(k.le.0)) then
         call 
     .    error('Dij_k: one of the partons is in the initial state')
         stop
      endif
#endif
      do mu = 1,4
         ki(mu) = p(mu,i)
         kj(mu) = p(mu,j)
         kk(mu) = p(mu,k)
      enddo

      do rho=1,4
         Q(rho) = ki(rho) + kj(rho) + kk(rho)
      enddo

      QQ = dotp(Q,Q)
      sqrtQ = dsqrt(QQ)
C     set the type of the emitter and spectator
      temitter = findEmitter(4,i,j,plist)
      tspectator = plist(k)
#ifdef VERBOSE
      if (temitter.eq.gluon) write(*,*) 'Dij_k: emit. = gluon'
      if (temitter.eq.quark) write(*,*) 'Dij_k: emit. = quark'
      if (temitter.eq.antiquark) write(*,*) 'Dij_k: emit. = aquark'
      if (tspectator.eq.gluon) write(*,*) 'Dij_k: spec. = gluon'
      if (tspectator.eq.quark) write(*,*) 'Dij_k: spec. = quark'
      if (tspectator.eq.antiquark) write(*,*) 'Dij_k: spec. = aquark'
#endif

      miq = massq(i)
      mjq = massq(j)
      mkq = massq(k)
      
#ifdef VERBOSE
      write(*,*) 'Dij_k: miq = ',miq
      write(*,*) 'Dij_k: mjq = ',mjq
      write(*,*) 'Dij_k: mkq = ',mkq
      if ((miq.lt.0).or.(mjq.lt.0).or.(mkq.lt.0)) then
         write(*,*) 'Dij_k: negative mass square'
         stop
      endif
#endif
      
      mui = dsqrt(miq) / sqrtQ
      muj = dsqrt(mjq) / sqrtQ
      muk = dsqrt(mkq) / sqrtQ


#ifdef VERBOSE
      write(*,*) 'Dij_k: mui = ',mui
      write(*,*) 'Dij_k: muj = ',muj
      write(*,*) 'Dij_k: muk = ',muk
#endif
C
C     set the emitter mass and the corrsponding casimir
C
      if (temitter .eq. gluon) then
         muij = 0.d0
         Tijq = CA
      endif
C
      if ((temitter .eq. quark).or.(temitter .eq. antiquark)) then
C     we assume QCD here: if the emitter is a quark i or j is
C     a massless gluon
         muij = dsqrt((massq(i)+massq(j))/QQ)
         Tijq = CF
      endif
      mijq = muij**2 * QQ
C
C     calculate the kinematics
C
      vtij_k = dsqrt( lambda(1.d0,muij**2,muk**2) )/
     .     (1.d0 - muij**2 - muk**2) 
C
      yij_k = dotp(ki,kj)/(dotp(ki,kj)+dotp(ki,kk)+dotp(kj,kk))
C
      zti = dotp(ki,kk)/(dotp(ki,kk)+dotp(kj,kk))
C
      ztj = 1.d0 - zti
C
      vij_k = dsqrt((2.d0*muk**2+(1.d0-mui**2-muj**2-muk**2)
     .     *(1.d0-yij_k))**2-4.d0*muk**2)/(1.d0-mui**2-muj**2-muk**2)/
     .     (1.d0-yij_k)
C
      vij_i = dsqrt((1.d0 -mui**2-muj**2-muk**2)**2*yij_k**2
     .     -4.d0*mui**2*muj**2)
     .     / ( (1.d0-mui**2-muj**2-muk**2)*yij_k + 2.d0 * mui**2)
C
      ym = 2.d0*mui*muj/(1.d0-mui**2-muj**2-muk**2)
C
      yp = 1.d0 - 2.d0*muk*(1.d0-muk)/(1.d0-mui**2-muj**2-muk**2)
C
      zp = (2.d0 * mui**2+(1.d0-mui**2-muj**2-muk**2)*yij_k)/
     .     2.d0/(mui**2+muj**2+(1.d0-mui**2-muj**2-muk**2)*yij_k)
C
      zm = zp * (1.d0 - vij_i*vij_k)
C
      zp = zp * (1.d0 + vij_i*vij_k)
C
      ztim = zti - 0.5d0 * (1.d0-vij_k)
C
      ztjm = ztj - 0.5d0 * (1.d0-vij_k)

      do rho=1,4
         ktk(rho) = dsqrt(lambda(QQ,mijq,mkq))/
     .        dsqrt(
     .        lambda(QQ,miq+2.d0*dotp(ki,kj)+mjq,mkq) )
     .        * ( kk(rho) - dotp(Q,kk)/QQ*Q(rho) )
C     .        * ( kk(rho) * (1.d0 - dotp(Q,kk)/QQ) 
C     .        - dotp(Q,kk)/QQ *( ki(rho) + kj(rho) ) )
     .        + (QQ+mkq-mijq)/2.d0/QQ*Q(rho)
         ktij(rho) = Q(rho) - ktk(rho)
      enddo
#ifdef VERBOSE
      write(*,*) 'Dij_k: dotp(ktk,ktk) = ',dotp(ktk,ktk)
      write(*,*) 'Dij_k: dotp(ktij,ktij) = ',dotp(ktij,ktij)
#endif
C
      if (temitter .eq. gluon) then
         if (plist(i).eq.gluon) then
#ifdef VERBOSE
      write(*,*) 'Dij_k:  g --> gg dipol'
#endif
C     g --> gg
            a = 16.d0 * pi * alphas * CA *
     .           ( 1.d0/(1.d0-zti*(1.d0-yij_k))
     .           + 1.d0/(1.d0-ztj*(1.d0-yij_k))
     .           - (2.d0-kappa*zp*zm)/vij_k )
            b = 16.d0 * pi * alphas * CA / vij_k /dotp(ki,kj)
            do rho = 1,4
               nvec(rho) = ztim * ki(rho) - ztjm * kj(rho)
            enddo
#ifdef VERBOSE
            write(*,*) 'Dij_k:  a,b', a,b
            write(*,*) 'Dij_k: dotp(n,n) = ', dotp(nvec,nvec) 
            write(*,*) 'Dij_k: dotp(n,ktij) = ',dotp(nvec,ktij) 
            write(*,*) 'Dij_k: dotp(ktij,ktij) = ',dotp(ktij,ktij) 
#endif
            if (tspectator.eq.quark) then
               do rho=1,4
                  pn(rho,-1) = p(rho,-1)
                  pn(rho,-2) = p(rho,-2)
                  pn(rho,1) = ktk(rho) 
                  pn(rho,2) = p(rho,2)
                  pn(rho,3) = ktij(rho)
               enddo
               Dij_k = obs3(pn)
               if (Dij_k.eq.0.d0) return
               Dij_k = Dij_k * Dggg(3,pn,Cg3q,a,b,nvec)
            else
               if (tspectator.eq.antiquark) then
                  do rho=1,4
                     pn(rho,-1) = p(rho,-1)
                     pn(rho,-2) = p(rho,-2)
                     pn(rho,1)  = p(rho,1)
                     pn(rho,2)  = ktk(rho)
                     pn(rho,3)  = ktij(rho)
                  enddo
                  Dij_k = obs3(pn)
                  if (Dij_k.eq.0.d0) return
                  Dij_k = Dij_k * Dggg(3,pn,Cg3qb,a,b,nvec)
               endif
            endif
         endif
         if ((plist(i).eq.quark).or.(plist(i).eq.quark)) then
C
C     g --> QQb
C     thats for the heavy quarks, although I think we will
C     not use it.
C
      write(*,*) 'Dij_k:  g --> QQb dipol'
      write(*,*) 'Dij_k: not checked so far'
      stop
            mquark2 = dotp(ki,ki)
            a = -8.d0*pi*alphas*TR/vij_k*(
     .           1.d0-2.d0*kappa * (zp*zm 
     .           - mquark2/(dotp(ki,ki)+2.d0*dotp(ki,kj)+dotp(kj,kj))))
            b = -8.d0*pi*alphas*TR / vij_k*4.d0 /
     .           (dotp(ki,ki)+2.d0*dotp(ki,kj)+dotp(kj,kj))
            do rho = 1,4
               nvec(rho) = ztim * ki(rho) 
     .              - ztjm * kj(rho)
            enddo
            if (tspectator.eq.gluon) then
               if (k.eq.3) then
                  ihelp = 4
               else
                  ihelp = 3
               endif
               do rho=1,4
                  pn(rho,-1) = p(rho,-1)
                  pn(rho,-2) = p(rho,-2)
                  pn(rho,1) =  p(rho,ihelp)
                  pn(rho,2) = ktij(rho) 
                  pn(rho,3) = ktk(rho)
               enddo
               Dij_k = obs3(pn)
               if (Dij_k.eq.0.d0) then
                  return
               endif
               Dij_k = Dij_k * Dgqqb(2,pn,Cxy,a,b,nvec)
            endif
         endif
      endif
      if ((temitter .eq. quark).or.(temitter .eq. antiquark)) then
C     
C     q  --> qg
C     qb --> qbg 
C
#ifdef VERBOSE
         write(*,*) 'Dij_k: q  --> qg, qb --> qbg  dipol'
         write(*,*) 'Dij_k: CF = ', CF
#endif
         if (plist(i).eq.gluon) then
            vgQ = 8.d0 * pi * alphas * CF * 
     .           (2.d0/(1.d0-ztj*(1.d0-yij_k))
     .           -vtij_k/vij_k * (1.d0+ztj+mijq/dotp(ki,kj)))
         else
            vgQ = 8.d0 * pi * alphas * CF * 
     .           (2.d0/(1.d0-zti*(1.d0-yij_k))
     .           -vtij_k/vij_k * (1.d0+zti+mijq/dotp(ki,kj)))
         endif
#ifdef VERBOSE
         write(*,*) 'Dij_k: vgQ = ', vgQ
         write(*,*) 'Dij_k: vtij_k = ',vtij_k
         write(*,*) 'Dij_k: vij_k = ',vij_k
#endif
C     
         if (temitter .eq. quark) then
            do rho=1,4
               pn(rho,-1) = p(rho,-1)
               pn(rho,-2) = p(rho,-2)
               pn(rho,1) =  ktij(rho)
            enddo
            if (tspectator.eq.gluon) then
               do rho=1,4
                  pn(rho,2) =  p(rho,2)
                  pn(rho,3) = ktk(rho)
               enddo
               Dij_k = obs3(pn)
               if (Dij_k.eq.0.d0) return
               Dij_k = Dij_k * Dqgq(1,pn,Cg3q,vgQ)
            else  
C     the spectator is the anti-quark, so we have to find out what
C     gluon is the collinear one and which is the remaining one
               if ((i.eq.3).or.(j.eq.3)) then
                  ihelp = 4
               else
                  ihelp = 3
               endif
               do rho=1,4
                  pn(rho,2) = ktk(rho) 
                  pn(rho,3) = p(rho,ihelp)
               enddo 
               Dij_k = obs3(pn)
               if (Dij_k.eq.0.d0) return
               Dij_k = Dij_k * Dqgq(1,pn,Cqqb,vgQ)
            endif
         endif
C
         if (temitter .eq. antiquark) then
C     
            do rho=1,4
               pn(rho,-1) = p(rho,-1)
               pn(rho,-2) = p(rho,-2)
               pn(rho,2) =  ktij(rho)
            enddo
            if (tspectator.eq.gluon) then
#ifdef VERBOSE
               write(*,*) 'Dij_k: calculate D_qg,g dipole '
#endif
               do rho=1,4
                  pn(rho,1) =  p(rho,1)
                  pn(rho,3) =  ktk(rho)
               enddo     
               Dij_k = obs3(pn)
               if (Dij_k.eq.0.d0) return
               Dij_k = Dij_k * Dqgq(2,pn,cg3qb,vgQ)
            else  
C     the spectator is the anti-quark, so we have to find out what
C     gluon is the collinear one and which is the remaining one
               if ((i.eq.3).or.(j.eq.3)) then
                  ihelp = 4
               else
                  ihelp = 3
               endif     
               do rho=1,4          
                  pn(rho,1) = ktk(rho) 
                  pn(rho,3) = p(rho,ihelp)
               enddo    
               Dij_k = obs3(pn)
               if (Dij_k.eq.0.d0) return
               Dij_k = Dij_k * Dqgq(2,pn,Cqqb,vgQ)
            endif
         endif
      endif
      Dij_k = -Dij_k / ( miq + 2.d0 * dotp(ki,kj) + mjq - mijq )
     >     / Tijq
#ifdef VERBOSE
*      write(*,*) 'Prefactor = ',
*     &      - 1.d0/ ( miq + 2.d0 * dotp(ki,kj) + mjq - mijq )
      write(*,'(A6,I3,I3,I3,A4,$)') 'Dij_k(',i,j,k,') = '
      write(*,'(E25.15)') Dij_k
#endif
      end
C*************************************************************************
      function findEmitter(n,i,j,plist)
      implicit none 
      integer findEmitter,n,i,j
      integer plist(-2:n)
#include "subtractions.h"      
C
C     Get the type of the emitter
C 
      if (((plist(i).eq.quark).and.(plist(j).eq.antiquark))
     _    .or.((plist(i).eq.antiquark).and.(plist(j).eq.quark))) 
     _    then
         findEmitter = gluon
         return
      endif

      if ((plist(i).eq.quark).or.(plist(j).eq.quark)) then
         findEmitter = quark
         return
      endif
C
      if ((plist(i).eq.antiquark).or.(plist(j).eq.antiquark))
     .     then
         findEmitter = antiquark
         return
      endif

      if ((plist(i).eq.gluon).and.(plist(j).eq.gluon)) then
         findEmitter = gluon
      endif
      end
C****************************************************************
      subroutine error(string)
      character*(*) string
      write(*,*) string
      end
C****************************************************************
C****************************************************************
      subroutine cpColor(Csource,Cdest)
      implicit none
      integer i,j
      double precision Csource(6,6),Cdest(6,6)
      do i=1,6
         do j= 1,6
            Cdest(i,j) = Csource(i,j)
         enddo
      enddo
      end
C****************************************************************
C     DUMMY FUNCTIONS FOR TESTING:
      function Dgqqb(iem,pn,Cxy,a,b,nvec)
      implicit none
      integer iem
      double precision Dgqqb,pn(-2:4),a,b,nvec(4),Cxy(6,6)
      Dgqqb = 0.d0
      end
      function Dqbgqb(iem,pn,Cxy,a,b,nvec)
      implicit none
      integer iem
      double precision Dqbgqb,pn(-2:4),a,b,nvec(4),Cxy(6,6)
      Dqbgqb = 0.d0
      end
