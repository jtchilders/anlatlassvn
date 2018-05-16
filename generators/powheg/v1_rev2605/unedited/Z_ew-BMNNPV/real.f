      subroutine setreal(p,fermion_flav,amp2)
      implicit none
      include 'nlegborn.h'
*
      real * 8 p(0:3,nlegreal)
      real * 8 pp(0:3,nlegreal)
      integer fermion_flav(nlegreal)
      real * 8 amp2

      integer tmp
      integer fflav(nlegreal)
      real * 8 const

*

      if (     fermion_flav(1).ne.22
     &    .and.fermion_flav(2).ne.22
     &    .and.fermion_flav(nlegreal).ne.22) then


          call setreal_st(p,fermion_flav,amp2)


      else

          pp = p
          fflav = fermion_flav

          const = 1d0

*

          if (fermion_flav(1).eq.22) then
                  
              pp(:,1) = -p(:,5)
              pp(:,5) = -p(:,1)

              fflav(5) = fermion_flav(1)
              fflav(1) = -fermion_flav(5)

              ! the 3 comes from the colour
              const = -3d0

          elseif (fermion_flav(2).eq.22) then
                  
              pp(:,2) = -p(:,5)
              pp(:,5) = -p(:,2)

              fflav(5) = fermion_flav(2)
              fflav(2) = -fermion_flav(5)

              ! the 3 comes from the colour
              const = -3d0

          endif

          call setreal_ew(pp,fflav,amp2)

          amp2 = const*amp2

      endif

*
      end

*

      subroutine setreal_st(p,rflav,amp2real)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      real * 8 p(0:3,nlegreal)
      integer rflav(nlegreal)
      real * 8 amp2real

      real * 8 charge(-5:5)
      data charge
     #     / 0.33333333333333333333d0, !   1d0/3
     #      -0.66666666666666666667d0, !  -2d0/3
     #       0.33333333333333333333d0, !   1d0/3 
     #      -0.66666666666666666667d0, !   -2d0/3
     #       0.33333333333333333333d0, !   1d0/3 
     #       0d0,                      !   0d0   
     #      -0.33333333333333333333d0, !   -1d0/3
     #       0.66666666666666666667d0, !   2d0/3   
     #      -0.33333333333333333333d0, !   -1d0/3
     #       0.66666666666666666667d0, !   2d0/3 
     #      -0.33333333333333333333d0/ !   -1d0/3

      real *8 amp2,ploc(0:3,nlegreal)

      integer i,mu,iii(nlegreal)
      integer ferm_type(nlegreal)
      real * 8 ferm_charge(nlegreal)

c     vector boson id and decay
      integer vdecaymode
      common/cvecbos/vdecaymode

      logical check,cond
      parameter (check=.false.)

      do i=1,nlegreal
         iii(i)=rflav(i)
         do mu=0,3
            ploc(mu,i)=p(mu,i)
         enddo
      enddo

      if(check) then
         do mu=0,3
            if(
     $dabs(ploc(mu,1)+ploc(mu,2)-ploc(mu,3)-ploc(mu,4)-ploc(mu,5))
     $.gt.1d-6) then
               write(*,*) 'momentum non conservation in setreal'
               stop
            endif
         enddo
         if(iii(3).le.0.or.iii(4).ge.0) then
            write(*,*) 'Error in setreal 1'
            stop
         endif
      endif
      ferm_type(3) = +1
      ferm_type(4) = -1

      if(mod(vdecaymode,2).eq.1) then
         ferm_charge(3) = -1d0
         ferm_charge(4) = +1d0
      elseif(mod(vdecaymode,2).eq.0) then
         ferm_charge(3) = 0d0
         ferm_charge(4) = 0d0
      else
         write(*,*) 'Error in vdecaymode in setreal'
         stop
      endif

      if((iii(1)+iii(2).eq.0).and.(iii(1)*iii(2).ne.0)) then
c     q q -> e- e+ g
         cond=.false.
         if((iii(1).eq.1).or.(iii(1).eq.3).or.(iii(1).eq.5)) then
c     d d~ -> Z g
            ferm_charge(1) = charge(1)
            ferm_charge(2) = charge(-1)
            ferm_charge(5) = 0d0
            ferm_type(1) = 1
            ferm_type(2) = -1
            ferm_type(5) = 0
            cond=.true.
         elseif((iii(1).eq.-1).or.(iii(1).eq.-3).or.(iii(1).eq.-5)) then
c     d~ d -> Z g
            ferm_charge(1) = charge(-1)
            ferm_charge(2) = charge(1)
            ferm_charge(5) = 0d0
            ferm_type(1) = -1
            ferm_type(2) = 1
            ferm_type(5) = 0
            cond=.true.
         elseif((iii(1).eq.2).or.(iii(1).eq.4)) then
c     u u~ -> Z g
            ferm_charge(1) = charge(2)
            ferm_charge(2) = charge(-2)
            ferm_charge(5) = 0d0
            ferm_type(1) = 1
            ferm_type(2) = -1
            ferm_type(5) = 0
            cond=.true.
         elseif((iii(1).eq.-2).or.(iii(1).eq.-4)) then
c     u~ u -> Z g
            ferm_charge(1) = charge(-2)
            ferm_charge(2) = charge(2)
            ferm_charge(5) = 0d0
            ferm_type(1) = -1
            ferm_type(2) = 1
            ferm_type(5) = 0
            cond=.true.
         endif
         
         if(cond) then
            call q_aq_to_l_al_g(ploc,ferm_type,ferm_charge,amp2)         
         else
            write(*,*) 'Error in setreal 2'
            stop
         endif

      elseif((iii(1).eq.0).and.(iii(2).eq.iii(5))) then
c     g q -> e- e+ q
         cond=.false.
         if((iii(2).eq.1).or.(iii(2).eq.3).or.(iii(2).eq.5)) then
c     g d -> Z d
            ferm_charge(1) = 0d0
            ferm_charge(2) = charge(1)
            ferm_charge(5) = charge(1)
            ferm_type(1) = 0
            ferm_type(2) = 1
            ferm_type(5) = 1
            cond=.true.
         elseif((iii(2).eq.-1).or.(iii(2).eq.-3).or.(iii(2).eq.-5)) then
c     g d~ -> Z d~
            ferm_charge(1) = 0d0
            ferm_charge(2) = charge(-1)
            ferm_charge(5) = charge(-1)
            ferm_type(1) = 0
            ferm_type(2) = -1
            ferm_type(5) = -1   
            cond=.true.
         elseif((iii(2).eq.2).or.(iii(2).eq.4)) then
c     g u -> Z u
            ferm_charge(1) = 0d0
            ferm_charge(2) = charge(2)
            ferm_charge(5) = charge(2)
            ferm_type(1) = 0
            ferm_type(2) = 1
            ferm_type(5) = 1
            cond=.true.
         elseif((iii(2).eq.-2).or.(iii(2).eq.-4)) then
c     g u~ -> Z u~
            ferm_charge(1) = 0d0
            ferm_charge(2) = charge(-2)
            ferm_charge(5) = charge(-2)
            ferm_type(1) = 0
            ferm_type(2) = -1
            ferm_type(5) = -1            
            cond=.true.
         endif
         if(cond) then
            call g_q_to_l_al_q(ploc,ferm_type,ferm_charge,amp2)
         else
            write(*,*) 'Error in setreal 3'
            stop
         endif
         
      elseif((iii(2).eq.0).and.(iii(1).eq.iii(5))) then
c     q g -> e- e+ q
         cond=.false.
         if((iii(1).eq.1).or.(iii(1).eq.3).or.(iii(1).eq.5)) then
c     d g -> Z d
            ferm_charge(2) = 0d0
            ferm_charge(1) = charge(1)
            ferm_charge(5) = charge(1)
            ferm_type(2) = 0
            ferm_type(1) = 1
            ferm_type(5) = 1
            cond=.true.
         elseif((iii(1).eq.-1).or.(iii(1).eq.-3).or.(iii(1).eq.-5)) then
c     d~ g -> Z d~
            ferm_charge(2) = 0d0
            ferm_charge(1) = charge(-1)
            ferm_charge(5) = charge(-1)
            ferm_type(2) = 0
            ferm_type(1) = -1
            ferm_type(5) = -1
            cond=.true.
         elseif((iii(1).eq.2).or.(iii(1).eq.4)) then
c     u g  -> Z u
            ferm_charge(2) = 0d0
            ferm_charge(1) = charge(2)
            ferm_charge(5) = charge(2)
            ferm_type(2) = 0
            ferm_type(1) = 1
            ferm_type(5) = 1
            cond=.true.
         elseif((iii(1).eq.-2).or.(iii(1).eq.-4)) then
c     u~ g -> Z u~
            ferm_charge(2) = 0d0
            ferm_charge(1) = charge(-2)
            ferm_charge(5) = charge(-2)
            ferm_type(2) = 0
            ferm_type(1) = -1
            ferm_type(5) = -1   
            cond=.true.
         endif
         if(cond) then
            call q_g_to_l_al_q(ploc,ferm_type,ferm_charge,amp2)
         else
            write(*,*) 'Error in setreal 4'
            stop
         endif

      else
         write(*,*) 'Error in setreal 5'
         stop
      endif

      amp2real=amp2

c     cancel as/(2pi) associated with amp2. It will be put back by real_ampsq
      amp2real = amp2real/(st_alpha/(2*pi))

      end




c this subroutine compute the Born amplitude for the process
c q(p1) aq(p2) -> Z(p3+p4) g(p5)  con Z -> f(p3) af(p4) 
c NUMERICALLY, with the bra/ket formalism, not by squaring the analytic 
c amplitude
c It gets the matrix p with all the momenta and gives   
c the amplitude squared (amp2) averaged over initial 
c polarization     
c
c         q  --->----------->------ g
c                     |
c                     |            f
c                     |          /  
c         aq ---<-----/\/\/\/\/\/
c                       Z/gamma \
c                                \ af
c     ferm_type = 1 fermion
c     ferm_type = -1 antifermion
c     fermion_charge = +2/3, -1/3, -2/3, +1/3

      subroutine q_aq_to_l_al_g(pphy,fermion_type,fermion_charge,
     #     amp2)
   
      implicit none
      include 'PhysPars.h'
      include 'pwhg_math.h'  
      include 'mathx.h'  
      include 'pwhg_st.h'  
      include 'nlegborn.h'  
      include 'pwhg_physpar.h'  
c the nleg 4-momentum vectors
c p(i,1) is the i-th component of vector p1...   
      integer fermion_type(nlegreal)
      real * 8 fermion_charge(nlegreal)
      real * 8 pphy(0:3,nlegreal)
      real * 8 amp2
*
      real*8 mlep2
      common/leptmass/mlep2
*
      real * 8 dotp
      external dotp
*
      real * 8 p1(0:3),p2(0:3),k1(0:3),k2(0:3),k(0:3)
      integer nu,i
      real * 8 p(0:3,nlegreal)
      real * 8 ferm_charge(nlegreal)
      integer ferm_type(nlegreal)
      integer utype_q,utype_f
      real * 8 q_q,v_q,a_q,q_f,v_f,a_f
      real*8 qq2,ql2
      real*8 ml2,ml4
      complex*16 glm,glmd,glp,glpd
      complex*16 gqm,gqmd,gqp,gqpd
      complex*16 glm2,glp2
      complex*16 gqm2,gqp2
      real*8 tmp
      real*8 s,sp,t,u
      real*8 s2,u2
      real*8 k1k,k2k,p1k,p2k
      real*8 p1k2,p2k2
      complex*16 densp,denspd,densp2,denspg

      integer ifirst
      data ifirst /0/
      save ifirst
      save glm,glp,glmd,glpd,glm2,glp2

      if (ifirst.eq.0) then
          ifirst = 1

          if (ql.ne.0d0) then
              glm = gl(0)
              glp = gl(1)
          else
              glm = gn(0)
              glp = gn(1)
          endif
          glmd = conjg(glm)
          glpd = conjg(glp)
          glm2 = glm*glmd
          glp2 = glp*glpd

          ml2 = mlep2
          ml4 = ml2*ml2
      endif


      if ((fermion_type(3).ne.1).and.(fermion_type(4).ne.-1)) then
         write(*,*) 'ERROR: this subroutine deals only with Z decay'
         stop
      endif
     

c  copy of local variables
      do i=1,nlegreal
         do nu=0,3
            p(nu,i) = pphy(nu,i)
         enddo
         ferm_charge(i) = fermion_charge(i)
         ferm_type(i) = fermion_type(i)
      enddo

c exchance particle 1 and 2
      if (ferm_type(1).eq.-1) then
        if (ferm_type(2).eq.1) then
            call exchange_momenta(p(0,1),p(0,2))
            tmp = ferm_charge(1)
            ferm_charge(1)=-ferm_charge(2)
            ferm_charge(2)=-tmp
            tmp = ferm_type(1)
            ferm_type(1)=ferm_type(2)
            ferm_type(2)=tmp
         else
            write(*,*) 'Error in the type of the quark 1-2'
            stop
         endif
      endif

c utype = +1 if this is an up-type quark (u,c,ubar,cbar)
c utype = -1 otherwise
      if (abs(ferm_charge(1)).eq.2d0/3) then
         utype_q = +1
         q_q = 2d0/3
      elseif (abs(ferm_charge(1)).eq.1d0/3) then
         utype_q = -1
         q_q = -1d0/3
      else
         write(*,*) 'What charge is this??', ferm_charge(1)
         stop
      endif

                 
      if (abs(ferm_charge(3)).eq.1d0) then
         utype_f = -1
         q_f = -1d0
      elseif (abs(ferm_charge(3)).eq.0d0) then
         utype_f = +1
         q_f = 0d0
      elseif (abs(ferm_charge(3)).eq.1d0/3) then
         utype_f = -1
         q_f = -1d0/3
      elseif (abs(ferm_charge(3)).eq.2d0/3) then
         utype_f = +1
         q_f = 2d0/3
      else
         write(*,*) 'What charge is this??',ferm_charge(3)
         stop
      endif

      ql = q_f
      qq = q_q
      
      if (abs(qq).eq.1d0/3d0) then
          gqm = gd(0)
          gqp = gd(1)
          qq = -1d0/3d0
      else
          gqm = gu(0)
          gqp = gu(1)
          qq = 2d0/3d0
      endif

      qq2 = qq*qq
      ql2 = ql*ql
      gqmd = conjg(gqm)
      gqpd = conjg(gqp)

      gqm2 = gqm*gqmd
      gqp2 = gqp*gqpd

      p1 = p(:,1)
      p2 = p(:,2)

      k1=p(:,3)
      k2=p(:,4)
      k=p(:,5)

      p1k=dotp(p1,k)
      p2k=dotp(p2,k)
      k1k=dotp(k1,k)
      k2k=dotp(k2,k)

      p1k2 = p1k**2
      p2k2 = p2k**2

      s =  2d0*dotp(p1,p2)
      t = -2d0*dotp(p1,k1) + ml2
      u = -2d0*dotp(p2,k1) + ml2
      sp = 2d0*dotp(k1,k2) + 2d0*ml2

      s2 = s**2
      u2 = u**2

      densp = 1d0/(sp  - ph_Zmass2 + ii*ph_ZmZw)
      denspd = dconjg(densp)
      densp2 = densp*denspd
      denspg = 1d0/sp


      amp2 =   (4*(-(densp2*((glmd*glp + glm*glpd)*(gqm2 + gqp2)*ml2*
     -            (-2*(p1k2 + p2k2) + 2*(p1k + p2k)*s - s2) + 
     -           glm2*(gqp2*(ml4*s + 
     -                t*(2*p1k2 + 2*p2k2 + 2*p1k*(p2k - s) + s2 + 
     -                    p2k*(-3*s + t)) + 
     -             (p1k*(-2*p2k + s - t) + (-p2k + s)*t)*u + p1k*u2 - 
     -                 ml2*(2*p1k2 + 2*p2k2 + p2k*(-3*s + t - u) - 
     -                    p1k*(s + t - u) + s*(s + t + u)) + 
     -                k1k*(4*p2k2 + 2*ml2*(p1k + p2k - s) - 2*p1k*t - 
     -                    2*p2k*(2*s + u) + s*(s + t + u))) + 
     -              gqm2*(ml4*s + p2k*t*(-2*p1k + s + t) + 
     -                 (2*p1k2 + 2*p2k2 + p1k*(2*p2k - 3*s - t) + 
     -                    s*(s + t) - p2k*(2*s + t))*u + p1k*u2 + 
     -                 k1k*(4*p1k2 + 2*ml2*(p1k + p2k - s) - 
     -                    2*p1k*(2*s + t) - 2*p2k*u + s*(s + t + u)) - 
     -                 ml2*(2*p1k2 + 2*p2k2 + p1k*(-3*s - t + u) - 
     -                    p2k*(s - t + u) + s*(s + t + u)))) + 
     -           glp2*(gqm2*(ml4*s + 
     -                t*(2*p1k2 + 2*p2k2 + 2*p1k*(p2k - s) + s2 + 
     -                    p2k*(-3*s + t)) + 
     -             (p1k*(-2*p2k + s - t) + (-p2k + s)*t)*u + p1k*u2 - 
     -                 ml2*(2*p1k2 + 2*p2k2 + p2k*(-3*s + t - u) - 
     -                    p1k*(s + t - u) + s*(s + t + u)) + 
     -                k1k*(4*p2k2 + 2*ml2*(p1k + p2k - s) - 2*p1k*t - 
     -                    2*p2k*(2*s + u) + s*(s + t + u))) + 
     -              gqp2*(ml4*s + p2k*t*(-2*p1k + s + t) + 
     -                 (2*p1k2 + 2*p2k2 + p1k*(2*p2k - 3*s - t) + 
     -                    s*(s + t) - p2k*(2*s + t))*u + p1k*u2 + 
     -                 k1k*(4*p1k2 + 2*ml2*(p1k + p2k - s) - 
     -                    2*p1k*(2*s + t) - 2*p2k*u + s*(s + t + u)) - 
     -                 ml2*(2*p1k2 + 2*p2k2 + p1k*(-3*s - t + u) - 
     -                    p2k*(s - t + u) + s*(s + t + u)))))) - 
     -      denspg*qq*(2*denspg*qq*
     -          (2*ml4*s + t*
     -        (2*p1k2 + 2*p2k2 - 2*p1k*s + s2 + 2*p2k*(-s + t)) + 
     -           (2*p1k2 + 2*p2k2 - 2*p1k*(s + t) - 2*p2k*(s + t) + 
     -               s*(s + 2*t))*u + 2*p1k*u2 + 
     -            2*k1k*(2*p1k2 + 2*p2k2 + 2*ml2*(p1k + p2k - s) - 
     -          2*p2k*s + s2 + s*t - 2*p1k*(s + t) - 2*p2k*u + s*u) - 
     -            2*ml2*(4*p1k2 + 4*p2k2 + p2k*(-4*s + t - u) + 
     -               p1k*(-4*s - t + u) + s*(2*s + t + u))) + 
     -         densp*ql*(glp*(gqp*
     -                (ml4*s + 
     -                  (2*p1k - s - t)*(2*k1k*p1k - k1k*s - p2k*t) + 
     -                  (2*p1k2 - 2*k1k*p2k + 2*p1k*p2k + 2*p2k2 + 
     -                     k1k*s - 3*p1k*s - 2*p2k*s + s2 - 
     -                     (p1k + p2k - s)*t)*u + p1k*u2) + 
     -               gqm*(ml4*s + 2*p1k2*t + 
     -                  t*(2*p2k2 + p2k*(-3*s + t - u) + s*(s + u)) + 
     -                p1k*(2*p2k*(t - u) + s*(-2*t + u) + u*(-t + u)) + 
     -                  k1k*(4*p2k2 - 2*p1k*t - 2*p2k*(2*s + u) + 
     -                     s*(s + t + u)))) + 
     -            glm*(gqm*(ml4*s + 
     -                  (2*p1k - s - t)*(2*k1k*p1k - k1k*s - p2k*t) + 
     -                  (2*p1k2 - 2*k1k*p2k + 2*p1k*p2k + 2*p2k2 + 
     -                     k1k*s - 3*p1k*s - 2*p2k*s + s2 - 
     -                     (p1k + p2k - s)*t)*u + p1k*u2) + 
     -               gqp*(ml4*s + 2*p1k2*t + 
     -                  t*(2*p2k2 + p2k*(-3*s + t - u) + s*(s + u)) + 
     -                p1k*(2*p2k*(t - u) + s*(-2*t + u) + u*(-t + u)) + 
     -                  k1k*(4*p2k2 - 2*p1k*t - 2*p2k*(2*s + u) + 
     -                     s*(s + t + u))))) + 
     -         denspd*ql*(glpd*
     -             (gqpd*(ml4*s + 
     -                  (2*p1k - s - t)*(2*k1k*p1k - k1k*s - p2k*t) + 
     -                  (2*p1k2 - 2*k1k*p2k + 2*p1k*p2k + 2*p2k2 + 
     -                     k1k*s - 3*p1k*s - 2*p2k*s + s2 - 
     -                     (p1k + p2k - s)*t)*u + p1k*u2) + 
     -               gqmd*(ml4*s + 2*p1k2*t + 
     -                  t*(2*p2k2 + p2k*(-3*s + t - u) + s*(s + u)) + 
     -                p1k*(2*p2k*(t - u) + s*(-2*t + u) + u*(-t + u)) + 
     -                  k1k*(4*p2k2 - 2*p1k*t - 2*p2k*(2*s + u) + 
     -                     s*(s + t + u)))) + 
     -            glmd*(gqmd*(ml4*s + 
     -                  (2*p1k - s - t)*(2*k1k*p1k - k1k*s - p2k*t) + 
     -                  (2*p1k2 - 2*k1k*p2k + 2*p1k*p2k + 2*p2k2 + 
     -                     k1k*s - 3*p1k*s - 2*p2k*s + s2 - 
     -                     (p1k + p2k - s)*t)*u + p1k*u2) + 
     -               gqpd*(ml4*s + 2*p1k2*t + 
     -                  t*(2*p2k2 + p2k*(-3*s + t - u) + s*(s + u)) + 
     -                p1k*(2*p2k*(t - u) + s*(-2*t + u) + u*(-t + u)) + 
     -                  k1k*(4*p2k2 - 2*p1k*t - 2*p2k*(2*s + u) + 
     -                     s*(s + t + u))))))))/(p1k*p2k)

      if (amp2.lt.0d0) amp2 = -amp2
      amp2 = amp2* el2_scheme*dconjg(el2_scheme) * (4*pi*st_alpha)
c     1/4 from average over initial-state polarization
c     1/nc^2 * nc = 1/nc from average over initial-state colors and sum over 
c     quark colors 
      amp2=  amp2*CF/4d0/nc

      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C AMPLITUDES RELATES BY CROSSING:
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c this subroutine compute the born amplitude for the process
c q(p1) g(p2) -> Z(p3+p4) q(p5)  con Z -> l-(p3) l+(p4) 
c NUMERICALLY, with the bra/ket formalism, not by squaring the analytic 
c amplitude. This amplitude is obtained by crossing the formula for 
c q_aq_to_l_al_g
      
      subroutine q_g_to_l_al_q(pphy,fermion_type,fermion_charge,
     #     amp2)
      implicit none
      integer nleg
      parameter (nleg=5)
      integer fermion_type(nleg)
      real * 8 fermion_charge(nleg)
      real * 8 pphy(0:3,nleg)
      real * 8 pp(0:3,nleg), ferm_charge(nleg)
      integer ferm_type(nleg)
      real * 8 amp2
      integer mu,i


      do i = 1,nleg
         do mu=0,3
            pp(mu,i) = pphy(mu,i)
         enddo
         ferm_charge(i) = fermion_charge(i)
         ferm_type(i) = fermion_type(i)
      enddo

      do mu=0,3
c     exchange initial gluon <-> final quark
         pp(mu,5) = -pphy(mu,2)
         pp(mu,2) = -pphy(mu,5)
      enddo
      
c no useful information is in ferm_type(2) or ferm_charge(2), 
c since it's the gluon, BEFORE the following exchange
      ferm_type(2) = -ferm_type(5)
c NOTE the MINUS sign     !!!
      ferm_charge(2) = -ferm_charge(5)

c     if the following two lines are missing 
      ferm_type(5)=0
      ferm_charge(5)=0d0 
c     ferm_type(5) and ferm_charge(5) don't contain
c     their correct values, but this does not affect 
c     the correct call of
      
      call q_aq_to_l_al_g(pp,ferm_type,ferm_charge,
     #     amp2)

c     correct for color average
      amp2 = amp2 * 3d0/8d0

      end


ccccccccccccccccccccccccccccccccccccccccccccc
c this subroutine compute the Born amplitude for the process
c g(p1) q(p2) -> Z(p3+p4) q(p5)   con Z -> l-(p3) l+(p4) 
c NUMERICALLY, with the bra/ket formalism, not by squaring the analytic 
c amplitude. This amplitude is obtained by crossing the formula for 
c q_aq_to_l_al_g


      subroutine g_q_to_l_al_q(pphy,fermion_type,fermion_charge,
     #     amp2)
      
      implicit none
      integer nleg
      parameter (nleg=5)
      integer fermion_type(nleg)
      real * 8 fermion_charge(nleg)
      real * 8 pphy(0:3,nleg)
      real * 8 pp(0:3,nleg), ferm_charge(nleg)
      integer ferm_type(nleg)
      real * 8 amp2
      integer mu,i


      do i = 1,nleg
         do mu=0,3
            pp(mu,i) = pphy(mu,i)
         enddo
         ferm_charge(i) = fermion_charge(i)
         ferm_type(i) = fermion_type(i)
      enddo

      do mu=0,3
c     exchange initial gluon <-> final quark
         pp(mu,5) = -pphy(mu,1)
         pp(mu,1) = -pphy(mu,5)
      enddo
      
c no useful information is in ferm_type(1) or ferm_charge(1), 
c since it's the gluon, BEFORE the following exchange
      ferm_type(1) = -ferm_type(5)
c NOTE the MINUS sign!!!
      ferm_charge(1) = -ferm_charge(5)

          
      call q_aq_to_l_al_g(pp,ferm_type,ferm_charge,
     #     amp2)
      
c     correct for color average
      amp2 = amp2 *3d0/8d0
     
      end
*
**
*

      subroutine setreal_ew(p,fermion_flav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'pwhg_physpar.h'
      include 'pwhg_math.h'
      include 'mathx.h'
*
      real * 8 p(0:3,nlegreal)
      integer fermion_flav(nlegreal)
      real * 8 amp2
*

      real*8 mlep2
      common/leptmass/mlep2

*
      real*8 dotp
      external dotp
      
      real*8 chargeofparticle
      external chargeofparticle

      complex*16 e_
      external e_
*
      real*8 p1(0:3),p2(0:3),k1(0:3),k2(0:3),k(0:3)
      real*8 qq2,ql2
      complex*16 glm,glmd,glp,glpd
      complex*16 gqm,gqmd,gqp,gqpd
      complex*16 glm2,glp2
      complex*16 gqm2,gqp2
      real*8 ml2,ml4
      complex*16 mz2m1,mz2dm1,mz2m2
      real*8 s,sp,t,u
      real*8 s2,u2,t2
      real*8 k1k,k2k,p1k,p2k
      real*8 p1k2,p2k2
      real*8 k1k2,k2k2
      complex*16 dens,densd,dens2,densg
      complex*16 densp,denspd,densp2,denspg

      integer ifirst
*
      data ifirst /0/
      save ifirst
      save glm,glp,glmd,glpd,glm2,glp2
      save mz2m1,mz2dm1,mz2m2

      ml2 = mlep2
      ml4 = ml2*ml2

      qq = chargeofparticle(fermion_flav(1))
      ql = chargeofparticle(fermion_flav(3))

      qq2 = qq**2
      ql2 = ql**2

      if (ifirst.eq.0) then
          ifirst = 1

          if (ql.ne.0d0) then
              glm = gl(0)
              glp = gl(1)
          else
              glm = gn(0)
              glp = gn(1)
          endif
          glmd = conjg(glm)
          glpd = conjg(glp)
          glm2 = glm*glmd
          glp2 = glp*glpd

          mz2m1  = 1d0/mz2 
          mz2dm1 = conjg(mz2m1)
          mz2m2  = mz2m1*mz2dm1
      endif

      if (qq.lt.0d0) then
          gqm = gd(0)
          gqp = gd(1)
      else
          gqm = gu(0)
          gqp = gu(1)
      endif

      gqmd = conjg(gqm)
      gqpd = conjg(gqp)

      gqm2 = gqm*gqmd
      gqp2 = gqp*gqpd

*
      if (fermion_flav(1).lt.0d0) then
          p2 = p(:,1)
          p1 = p(:,2)
      else
          p1 = p(:,1)
          p2 = p(:,2)
      endif
      k1 = p(:,3)
      k2 = p(:,4)
      k  = p(:,5)

      p1k=dotp(p1,k)
      p2k=dotp(p2,k)
      k1k=dotp(k1,k)
      k2k=dotp(k2,k)

      p1k2 = p1k**2
      p2k2 = p2k**2
      k1k2 = k1k**2
      k2k2 = k2k**2

      s =  2d0*dotp(p1,p2)
      t = -2d0*dotp(p1,k1) + ml2
      u = -2d0*dotp(p2,k1) + ml2
      sp = 2d0*dotp(k1,k2) + 2d0*ml2

      s2 = s**2
      t2 = t**2
      u2 = u**2

      dens = 1d0/(s - ph_Zmass2 + ii*ph_ZmZw)
      densd = dconjg(dens)
      dens2 = dens*densd
      densg = 1d0/s

      densp = 1d0/(sp - ph_Zmass2 + ii*ph_ZmZw)
      denspd = dconjg(densp)
      densp2 = densp*denspd
      denspg = 1d0/sp

*
      amp2 =
     -  (-4*densg**2*qq2*(k1k*
     -       (t*(2*p2k*t + s*(2*k1k - 2*p1k + s + t)) + 
     -         (s2 - 2*p1k*t + 4*s*t + 2*t2 - 2*p2k*(s + t) + 
     -            2*k1k*(s + 2*t))*u + (2*p1k + s + 2*t)*u2) - 
     -      4*k2k2*(p2k*(s + t) + p1k*(-4*p2k + s + u)) - 
     -      2*k2k*(4*k1k2*s + t*(-2*p2k2 + 2*p1k*(p2k - s) + s*t) - 
     -         2*(p1k2 - p1k*p2k + p2k*s)*u + s*u2 + 
     -      2*k1k*(-(p1k*s) - p2k*s + 2*p2k*t + s*t + 2*p1k*u + s*u))))/
     -  (k1k*k2k2) +
     -  (-2*dens2*(glp2*(gqp2*(4*k2k2*p2k*(2*p1k - s - t) - 
     -            4*k2k*(k1k - p2k)*(k1k*s + p2k*t) - 
     -         4*k2k*(p1k*(k1k + p2k) + (k1k - p2k)*s)*u - 2*k2k*s*u2 + 
     -         k1k*u*(2*p1k*u + (s + t)*(2*k1k - 2*p2k + s + t + u))) + 
     -         gqm2*(4*k2k2*p1k*(2*p2k - s - u) + 
     -          k1k*t*(2*p2k*t + (s + u)*(2*k1k - 2*p1k + s + t + u)) - 
     -            2*k2k*(2*k1k2*s + t*(2*p1k*(p2k - s) + s*t) - 
     -             2*p1k2*u + 2*k1k*((p2k + s)*t + p1k*(-s + u))))) + 
     -      glm2*(gqm2*(4*k2k2*p2k*(2*p1k - s - t) - 
     -            4*k2k*(k1k - p2k)*(k1k*s + p2k*t) - 
     -         4*k2k*(p1k*(k1k + p2k) + (k1k - p2k)*s)*u - 2*k2k*s*u2 + 
     -         k1k*u*(2*p1k*u + (s + t)*(2*k1k - 2*p2k + s + t + u))) + 
     -         gqp2*(4*k2k2*p1k*(2*p2k - s - u) + 
     -         k1k*t*(2*p2k*t + (s + u)*(2*k1k - 2*p1k + s + t + u)) - 
     -            2*k2k*(2*k1k2*s + t*(2*p1k*(p2k - s) + s*t) - 
     -             2*p1k2*u + 2*k1k*((p2k + s)*t + p1k*(-s + u)))))))/
     -  (k1k*k2k2) +
     -         (2*densd*densg*qq*(glpd*
     -       (gqpd*(4*k2k2*p2k*(2*p1k - s - t) - 
     -            4*k2k*(k1k - p2k)*(k1k*s + p2k*t) - 
     -         4*k2k*(p1k*(k1k + p2k) + (k1k - p2k)*s)*u - 2*k2k*s*u2 + 
     -         k1k*u*(2*p1k*u + (s + t)*(2*k1k - 2*p2k + s + t + u))) + 
     -         gqmd*(4*k2k2*p1k*(2*p2k - s - u) + 
     -         k1k*t*(2*p2k*t + (s + u)*(2*k1k - 2*p1k + s + t + u)) - 
     -            2*k2k*(2*k1k2*s + t*(2*p1k*(p2k - s) + s*t) - 
     -               2*p1k2*u + 2*k1k*((p2k + s)*t + p1k*(-s + u))))) + 
     -      glmd*(gqmd*(4*k2k2*p2k*(2*p1k - s - t) - 
     -            4*k2k*(k1k - p2k)*(k1k*s + p2k*t) - 
     -         4*k2k*(p1k*(k1k + p2k) + (k1k - p2k)*s)*u - 2*k2k*s*u2 + 
     -         k1k*u*(2*p1k*u + (s + t)*(2*k1k - 2*p2k + s + t + u))) + 
     -         gqpd*(4*k2k2*p1k*(2*p2k - s - u) + 
     -          k1k*t*(2*p2k*t + (s + u)*(2*k1k - 2*p1k + s + t + u)) - 
     -            2*k2k*(2*k1k2*s + t*(2*p1k*(p2k - s) + s*t) - 
     -               2*p1k2*u + 2*k1k*((p2k + s)*t + p1k*(-s + u)))))))/
     -  (k1k*k2k2) +
     -         (2*dens*densg*qq*(glp*(gqp*
     -          (4*k2k2*p2k*(2*p1k - s - t) - 
     -            4*k2k*(k1k - p2k)*(k1k*s + p2k*t) - 
     -         4*k2k*(p1k*(k1k + p2k) + (k1k - p2k)*s)*u - 2*k2k*s*u2 + 
     -         k1k*u*(2*p1k*u + (s + t)*(2*k1k - 2*p2k + s + t + u))) + 
     -         gqm*(4*k2k2*p1k*(2*p2k - s - u) + 
     -          k1k*t*(2*p2k*t + (s + u)*(2*k1k - 2*p1k + s + t + u)) - 
     -            2*k2k*(2*k1k2*s + t*(2*p1k*(p2k - s) + s*t) - 
     -             2*p1k2*u + 2*k1k*((p2k + s)*t + p1k*(-s + u))))) + 
     -      glm*(gqm*(4*k2k2*p2k*(2*p1k - s - t) - 
     -            4*k2k*(k1k - p2k)*(k1k*s + p2k*t) - 
     -         4*k2k*(p1k*(k1k + p2k) + (k1k - p2k)*s)*u - 2*k2k*s*u2 + 
     -         k1k*u*(2*p1k*u + (s + t)*(2*k1k - 2*p2k + s + t + u))) + 
     -         gqp*(4*k2k2*p1k*(2*p2k - s - u) + 
     -          k1k*t*(2*p2k*t + (s + u)*(2*k1k - 2*p1k + s + t + u)) - 
     -            2*k2k*(2*k1k2*s + t*(2*p1k*(p2k - s) + s*t) - 
     -             2*p1k2*u + 2*k1k*((p2k + s)*t + p1k*(-s + u)))))))/
     -  (k1k*k2k2) +
     -       (-8*densg*denspg*qq**3*
     -    (2*k1k**3*(p1k - p2k)*s + 
     -      2*k1k2*(((p1k - p2k)*p2k - p1k*s)*t + 
     -         (p1k2 - p1k*(p2k + 2*t) + p2k*(s + 2*t))*u) - 
     -      k2k*(p2k*t - p1k*u)*
     -       (-(s*t) - (s + 2*t)*u + p2k*(s + 3*t + u) + 
     -         p1k*(-4*p2k + s + t + 3*u)) + 
     -      k1k*(p1k*(2*p1k*s - (-2*p2k + s)**2)*t + 
     -         (2*p1k2*(2*p2k + t) - 2*p1k*s*(2*p2k + t) + 
     -            p2k*(s2 + 2*s*t + 2*t2 - 2*p2k*(s + t)))*u - 
     -         2*p1k*t*u2 - k2k*(p1k - p2k)*
     -       (8*p1k*p2k - s2 - 4*p2k*t + s*t + (-4*p1k + s + 4*t)*u))))/
     -  (k1k*k2k*p1k*p2k) +
     -         (2*densd*denspg*qq2*(glpd*
     -       (gqpd*(2*k1k**3*p1k*s + 
     -            2*k1k2*(p1k*p2k*(-s + t) + 
     -               (p1k2 - p1k*t + p2k*(s + t))*u) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (-((s + t)*u) + 2*p1k*(-p2k + u) + p2k*(s + t + u)) - 
     -            k1k*(k2k*p2k*(4*p1k - s)*(2*p1k - s - t) - 
     -               2*p1k2*p2k*u + p2k*(2*p2k - s - t)*(s + t)*u - 
     -               k2k*(4*p1k2 - 2*p1k*(s + t) + p2k*(s + 2*t))*u + 
     -               p1k*(2*p2k2*t + 2*p2k*(2*s + t)*u + t*u2))) + 
     -         gqmd*(-2*k1k**3*p2k*s + 
     -            2*k1k2*(p1k*s*(p2k - t) - p1k*(p2k + t)*u + 
     -               p2k*t*(-p2k + u)) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (t*(2*p2k - s - u) + p1k*(-2*p2k + s + t + u)) + 
     -            k1k*(p2k*t2*u + 2*p1k2*(s*t + (p2k + t)*u) - 
     -               p1k*t*(2*p2k2 + (s + u)**2 - 2*p2k*(2*s + u)) + 
     -               k2k*(2*p2k*t*(-2*p2k + s + u) + 
     -                  p1k*(8*p2k2 - 6*p2k*s + s2 - s*t - 4*p2k*u + 
     -                     s*u - 2*t*u))))) + 
     -      glmd*(gqmd*(2*k1k**3*p1k*s + 
     -            2*k1k2*(p1k*p2k*(-s + t) + 
     -               (p1k2 - p1k*t + p2k*(s + t))*u) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (-((s + t)*u) + 2*p1k*(-p2k + u) + p2k*(s + t + u)) - 
     -            k1k*(k2k*p2k*(4*p1k - s)*(2*p1k - s - t) - 
     -               2*p1k2*p2k*u + p2k*(2*p2k - s - t)*(s + t)*u - 
     -               k2k*(4*p1k2 - 2*p1k*(s + t) + p2k*(s + 2*t))*u + 
     -               p1k*(2*p2k2*t + 2*p2k*(2*s + t)*u + t*u2))) + 
     -         gqpd*(-2*k1k**3*p2k*s + 
     -            2*k1k2*(p1k*s*(p2k - t) - p1k*(p2k + t)*u + 
     -               p2k*t*(-p2k + u)) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (t*(2*p2k - s - u) + p1k*(-2*p2k + s + t + u)) + 
     -            k1k*(p2k*t2*u + 2*p1k2*(s*t + (p2k + t)*u) - 
     -               p1k*t*(2*p2k2 + (s + u)**2 - 2*p2k*(2*s + u)) + 
     -               k2k*(2*p2k*t*(-2*p2k + s + u) + 
     -                  p1k*(8*p2k2 - 6*p2k*s + s2 - s*t - 4*p2k*u + 
     -                     s*u - 2*t*u)))))))/(k1k*k2k*p1k*p2k) +
     -         (2*dens*denspg*qq2*(glp*
     -       (gqp*(2*k1k**3*p1k*s + 
     -            2*k1k2*(p1k*p2k*(-s + t) + 
     -               (p1k2 - p1k*t + p2k*(s + t))*u) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (-((s + t)*u) + 2*p1k*(-p2k + u) + p2k*(s + t + u)) - 
     -            k1k*(k2k*p2k*(4*p1k - s)*(2*p1k - s - t) - 
     -               2*p1k2*p2k*u + p2k*(2*p2k - s - t)*(s + t)*u - 
     -               k2k*(4*p1k2 - 2*p1k*(s + t) + p2k*(s + 2*t))*u + 
     -               p1k*(2*p2k2*t + 2*p2k*(2*s + t)*u + t*u2))) + 
     -         gqm*(-2*k1k**3*p2k*s + 
     -            2*k1k2*(p1k*s*(p2k - t) - p1k*(p2k + t)*u + 
     -               p2k*t*(-p2k + u)) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (t*(2*p2k - s - u) + p1k*(-2*p2k + s + t + u)) + 
     -            k1k*(p2k*t2*u + 2*p1k2*(s*t + (p2k + t)*u) - 
     -               p1k*t*(2*p2k2 + (s + u)**2 - 2*p2k*(2*s + u)) + 
     -               k2k*(2*p2k*t*(-2*p2k + s + u) + 
     -                  p1k*(8*p2k2 - 6*p2k*s + s2 - s*t - 4*p2k*u + 
     -                     s*u - 2*t*u))))) + 
     -      glm*(gqm*(2*k1k**3*p1k*s + 
     -            2*k1k2*(p1k*p2k*(-s + t) + 
     -               (p1k2 - p1k*t + p2k*(s + t))*u) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (-((s + t)*u) + 2*p1k*(-p2k + u) + p2k*(s + t + u)) - 
     -            k1k*(k2k*p2k*(4*p1k - s)*(2*p1k - s - t) - 
     -               2*p1k2*p2k*u + p2k*(2*p2k - s - t)*(s + t)*u - 
     -               k2k*(4*p1k2 - 2*p1k*(s + t) + p2k*(s + 2*t))*u + 
     -               p1k*(2*p2k2*t + 2*p2k*(2*s + t)*u + t*u2))) + 
     -         gqp*(-2*k1k**3*p2k*s + 
     -            2*k1k2*(p1k*s*(p2k - t) - p1k*(p2k + t)*u + 
     -               p2k*t*(-p2k + u)) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (t*(2*p2k - s - u) + p1k*(-2*p2k + s + t + u)) + 
     -            k1k*(p2k*t2*u + 2*p1k2*(s*t + (p2k + t)*u) - 
     -               p1k*t*(2*p2k2 + (s + u)**2 - 2*p2k*(2*s + u)) + 
     -               k2k*(2*p2k*t*(-2*p2k + s + u) + 
     -                  p1k*(8*p2k2 - 6*p2k*s + s2 - s*t - 4*p2k*u + 
     -                     s*u - 2*t*u)))))))/(k1k*k2k*p1k*p2k) +
     -        (-8*denspg**2*qq**4*(t*(2*p2k2 + s2 + 2*p2k*(-s + t)) + 
     -      (2*p2k2 - 2*p2k*(s + t) + s*(s + 2*t))*u + 2*p1k2*(t + u) - 
     -      2*p1k*((t - u)*u + s*(t + u)) + 
     -      2*k1k*(2*p1k2 + 2*p2k2 - 2*p1k*(s + t) - 2*p2k*(s + u) + 
     -         s*(s + t + u))))/(p1k*p2k) +
     -         (-4*densp2*qq2*(glp2*
     -       (gqp2*((2*p1k - s - t)*(2*k1k*p1k - k1k*s - p2k*t) + 
     -            (2*p1k2 - 2*k1k*p2k + 2*p1k*p2k + 2*p2k2 + k1k*s - 
     -               3*p1k*s - 2*p2k*s + s2 - (p1k + p2k - s)*t)*u + 
     -            p1k*u2) + gqm2*
     -          (2*p1k2*t + p1k*
     -             (-2*s*t + 2*p2k*(t - u) + s*u - t*u + u2) + 
     -            t*(2*p2k2 + p2k*(-3*s + t - u) + s*(s + u)) + 
     -         k1k*(4*p2k2 - 2*p1k*t - 2*p2k*(2*s + u) + s*(s + t + u)))
     -     ) + glm2*(gqm2*((2*p1k - s - t)*(2*k1k*p1k - k1k*s - p2k*t) + 
     -            (2*p1k2 - 2*k1k*p2k + 2*p1k*p2k + 2*p2k2 + k1k*s - 
     -               3*p1k*s - 2*p2k*s + s2 - (p1k + p2k - s)*t)*u + 
     -            p1k*u2) + gqp2*
     -          (2*p1k2*t + p1k*
     -             (-2*s*t + 2*p2k*(t - u) + s*u - t*u + u2) + 
     -            t*(2*p2k2 + p2k*(-3*s + t - u) + s*(s + u)) + 
     -         k1k*(4*p2k2 - 2*p1k*t - 2*p2k*(2*s + u) + s*(s + t + u)))
     -         )))/(p1k*p2k)+
     -         (2*densg*denspd*qq2*(glmd*
     -       (-(gqpd*(2*k1k**3*p2k*s + 
     -              2*k1k2*(p1k*s*(t + p2k*(-1 + mz2dm1*t)) + 
     -               p2k*t*(p2k - u) + p1k*(p2k + t + mz2dm1*p2k*t)*u) + 
     -              k2k*(p2k*t - p1k*u)*
     -               (t*(2*p2k - s - u) + p1k*(-2*p2k + s + t + u)) + 
     -              k1k*(-(p2k*t2*u) - 
     -                 2*p1k2*
     -               (s*(t + mz2dm1*p2k*t) + (p2k + t + mz2dm1*p2k*t)*u)
     -                  + k2k*(2*p2k*t*(2*p2k - s - u) + 
     -                 p1k*(-8*p2k2 + 6*p2k*s - s2 + s*t + 4*p2k*u - 
     -                       s*u + 2*t*u)) + 
     -              p1k*t*((s + u)**2 - 2*p2k2*(-1 + mz2dm1*(s + u)) + 
     -                    p2k*(mz2dm1*s2 + 
     -                       s*(-4 + mz2dm1*t + 2*mz2dm1*u) + 
     -                       u*(-2 + mz2dm1*(t + u))))))) + 
     -         gqmd*(2*k1k**3*p1k*s + 
     -            2*k1k2*(p1k*p2k*(-s + t) + 
     -          (p1k2 + p2k*(s + t) - p1k*(t + mz2dm1*p2k*(s + t)))*u)
     -              - k2k*(p2k*t - p1k*u)*
     -             (-((s + t)*u) + 2*p1k*(-p2k + u) + p2k*(s + t + u)) - 
     -            k1k*(k2k*p2k*(4*p1k - s)*(2*p1k - s - t) + 
     -               p2k*(2*p2k - s - t)*(s + t)*u - 
     -               2*p1k2*p2k*(1 + mz2dm1*(s + t))*u - 
     -               k2k*(4*p1k2 - 2*p1k*(s + t) + p2k*(s + 2*t))*u + 
     -               p1k*(t*u2 + 2*p2k2*(t - mz2dm1*(s + t)*u) + 
     -              p2k*u*(mz2dm1*s2 + s*(4 + 2*mz2dm1*t + mz2dm1*u) + 
     -                     t*(2 + mz2dm1*(t + u))))))) + 
     -      glpd*(-(gqmd*(2*k1k**3*p2k*s + 
     -              2*k1k2*(p1k*s*(t + p2k*(-1 + mz2dm1*t)) + 
     -              p2k*t*(p2k - u) + p1k*(p2k + t + mz2dm1*p2k*t)*u) + 
     -              k2k*(p2k*t - p1k*u)*
     -               (t*(2*p2k - s - u) + p1k*(-2*p2k + s + t + u)) + 
     -              k1k*(-(p2k*t2*u) - 
     -                 2*p1k2*
     -               (s*(t + mz2dm1*p2k*t) + (p2k + t + mz2dm1*p2k*t)*u)
     -                  + k2k*(2*p2k*t*(2*p2k - s - u) + 
     -                 p1k*(-8*p2k2 + 6*p2k*s - s2 + s*t + 4*p2k*u - 
     -                       s*u + 2*t*u)) + 
     -              p1k*t*((s + u)**2 - 2*p2k2*(-1 + mz2dm1*(s + u)) + 
     -                    p2k*(mz2dm1*s2 + 
     -                       s*(-4 + mz2dm1*t + 2*mz2dm1*u) + 
     -                       u*(-2 + mz2dm1*(t + u))))))) + 
     -         gqpd*(2*k1k**3*p1k*s + 
     -            2*k1k2*(p1k*p2k*(-s + t) + 
     -            (p1k2 + p2k*(s + t) - p1k*(t + mz2dm1*p2k*(s + t)))*u)
     -              - k2k*(p2k*t - p1k*u)*
     -             (-((s + t)*u) + 2*p1k*(-p2k + u) + p2k*(s + t + u)) - 
     -            k1k*(k2k*p2k*(4*p1k - s)*(2*p1k - s - t) + 
     -               p2k*(2*p2k - s - t)*(s + t)*u - 
     -               2*p1k2*p2k*(1 + mz2dm1*(s + t))*u - 
     -               k2k*(4*p1k2 - 2*p1k*(s + t) + p2k*(s + 2*t))*u + 
     -               p1k*(t*u2 + 2*p2k2*(t - mz2dm1*(s + t)*u) + 
     -              p2k*u*(mz2dm1*s2 + s*(4 + 2*mz2dm1*t + mz2dm1*u) + 
     -                 t*(2 + mz2dm1*(t + u)))))))))/(k1k*k2k*p1k*p2k) +
     -     (-2*dens*denspd*qq*(glp2*
     -       (gqp2*(2*k1k**3*p1k*s + 
     -            2*k1k2*(p1k*p2k*(-s + t) + 
     -               (p1k2 - p1k*t + p2k*(s + t))*u) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (-((s + t)*u) + 2*p1k*(-p2k + u) + p2k*(s + t + u)) - 
     -            k1k*(k2k*p2k*(4*p1k - s)*(2*p1k - s - t) - 
     -               2*p1k2*p2k*u + p2k*(2*p2k - s - t)*(s + t)*u - 
     -               k2k*(4*p1k2 - 2*p1k*(s + t) + p2k*(s + 2*t))*u + 
     -               p1k*(2*p2k2*t + 2*p2k*(2*s + t)*u + t*u2))) + 
     -         gqm2*(-2*k1k**3*p2k*s + 
     -            2*k1k2*(p1k*s*(p2k - t) - p1k*(p2k + t)*u + 
     -               p2k*t*(-p2k + u)) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (t*(2*p2k - s - u) + p1k*(-2*p2k + s + t + u)) + 
     -            k1k*(p2k*t2*u + 2*p1k2*(s*t + (p2k + t)*u) - 
     -               p1k*t*(2*p2k2 + (s + u)**2 - 2*p2k*(2*s + u)) + 
     -               k2k*(2*p2k*t*(-2*p2k + s + u) + 
     -                  p1k*(8*p2k2 - 6*p2k*s + s2 - s*t - 4*p2k*u + 
     -                     s*u - 2*t*u))))) + 
     -      glm2*(gqm2*(2*k1k**3*p1k*s + 
     -            2*k1k2*(p1k*p2k*(-s + t) + 
     -               (p1k2 - p1k*t + p2k*(s + t))*u) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (-((s + t)*u) + 2*p1k*(-p2k + u) + p2k*(s + t + u)) - 
     -            k1k*(k2k*p2k*(4*p1k - s)*(2*p1k - s - t) - 
     -               2*p1k2*p2k*u + p2k*(2*p2k - s - t)*(s + t)*u - 
     -               k2k*(4*p1k2 - 2*p1k*(s + t) + p2k*(s + 2*t))*u + 
     -               p1k*(2*p2k2*t + 2*p2k*(2*s + t)*u + t*u2))) + 
     -         gqp2*(-2*k1k**3*p2k*s + 
     -            2*k1k2*(p1k*s*(p2k - t) - p1k*(p2k + t)*u + 
     -               p2k*t*(-p2k + u)) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (t*(2*p2k - s - u) + p1k*(-2*p2k + s + t + u)) + 
     -            k1k*(p2k*t2*u + 2*p1k2*(s*t + (p2k + t)*u) - 
     -               p1k*t*(2*p2k2 + (s + u)**2 - 2*p2k*(2*s + u)) + 
     -               k2k*(2*p2k*t*(-2*p2k + s + u) + 
     -                  p1k*(8*p2k2 - 6*p2k*s + s2 - s*t - 4*p2k*u + 
     -                     s*u - 2*t*u)))))))/(k1k*k2k*p1k*p2k)+
     -         (4*denspd*denspg*qq**3*
     -    (glpd*(gqpd*((2*p1k - s - t)*(2*k1k*p1k - k1k*s - p2k*t) + 
     -            (2*p1k2 - 2*k1k*p2k + 2*p1k*p2k + 2*p2k2 + k1k*s - 
     -               3*p1k*s - 2*p2k*s + s2 - (p1k + p2k - s)*t)*u + 
     -            p1k*u2) + gqmd*
     -          (2*p1k2*t + p1k*
     -             (-2*s*t + 2*p2k*(t - u) + s*u - t*u + u2) + 
     -            t*(2*p2k2 + p2k*(-3*s + t - u) + s*(s + u)) + 
     -         k1k*(4*p2k2 - 2*p1k*t - 2*p2k*(2*s + u) + s*(s + t + u)))
     -     ) + glmd*(gqmd*((2*p1k - s - t)*(2*k1k*p1k - k1k*s - p2k*t) + 
     -            (2*p1k2 - 2*k1k*p2k + 2*p1k*p2k + 2*p2k2 + k1k*s - 
     -               3*p1k*s - 2*p2k*s + s2 - (p1k + p2k - s)*t)*u + 
     -            p1k*u2) + gqpd*
     -          (2*p1k2*t + p1k*
     -             (-2*s*t + 2*p2k*(t - u) + s*u - t*u + u2) + 
     -            t*(2*p2k2 + p2k*(-3*s + t - u) + s*(s + u)) + 
     -         k1k*(4*p2k2 - 2*p1k*t - 2*p2k*(2*s + u) + s*(s + t + u)))
     -         )))/(p1k*p2k)+
     -        (2*densg*densp*qq2*(glp*
     -       (gqp*(2*k1k**3*p1k*s + 
     -            2*k1k2*(p1k*p2k*(-s + t) + 
     -               (p1k2 - p1k*t + p2k*(s + t))*u) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (-((s + t)*u) + 2*p1k*(-p2k + u) + p2k*(s + t + u)) - 
     -            k1k*(k2k*p2k*(4*p1k - s)*(2*p1k - s - t) - 
     -               2*p1k2*p2k*u + p2k*(2*p2k - s - t)*(s + t)*u - 
     -               k2k*(4*p1k2 - 2*p1k*(s + t) + p2k*(s + 2*t))*u + 
     -               p1k*(2*p2k2*t + 2*p2k*(2*s + t)*u + t*u2))) + 
     -         gqm*(-2*k1k**3*p2k*s + 
     -            2*k1k2*(p1k*s*(p2k - t) - p1k*(p2k + t)*u + 
     -               p2k*t*(-p2k + u)) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (t*(2*p2k - s - u) + p1k*(-2*p2k + s + t + u)) + 
     -            k1k*(p2k*t2*u + 2*p1k2*(s*t + (p2k + t)*u) - 
     -               p1k*t*(2*p2k2 + (s + u)**2 - 2*p2k*(2*s + u)) + 
     -               k2k*(2*p2k*t*(-2*p2k + s + u) + 
     -                  p1k*(8*p2k2 - 6*p2k*s + s2 - s*t - 4*p2k*u + 
     -                     s*u - 2*t*u))))) + 
     -      glm*(gqm*(2*k1k**3*p1k*s + 
     -            2*k1k2*(p1k*p2k*(-s + t) + 
     -               (p1k2 - p1k*t + p2k*(s + t))*u) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (-((s + t)*u) + 2*p1k*(-p2k + u) + p2k*(s + t + u)) - 
     -            k1k*(k2k*p2k*(4*p1k - s)*(2*p1k - s - t) - 
     -               2*p1k2*p2k*u + p2k*(2*p2k - s - t)*(s + t)*u - 
     -               k2k*(4*p1k2 - 2*p1k*(s + t) + p2k*(s + 2*t))*u + 
     -               p1k*(2*p2k2*t + 2*p2k*(2*s + t)*u + t*u2))) + 
     -         gqp*(-2*k1k**3*p2k*s + 
     -            2*k1k2*(p1k*s*(p2k - t) - p1k*(p2k + t)*u + 
     -               p2k*t*(-p2k + u)) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (t*(2*p2k - s - u) + p1k*(-2*p2k + s + t + u)) + 
     -            k1k*(p2k*t2*u + 2*p1k2*(s*t + (p2k + t)*u) - 
     -               p1k*t*(2*p2k2 + (s + u)**2 - 2*p2k*(2*s + u)) + 
     -               k2k*(2*p2k*t*(-2*p2k + s + u) + 
     -                  p1k*(8*p2k2 - 6*p2k*s + s2 - s*t - 4*p2k*u + 
     -                     s*u - 2*t*u)))))))/(k1k*k2k*p1k*p2k) +
     -         (-2*densd*densp*qq*(glp2*
     -       (gqp2*(2*k1k**3*p1k*s + 
     -            2*k1k2*(p1k*p2k*(-s + t) + 
     -               (p1k2 - p1k*t + p2k*(s + t))*u) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (-((s + t)*u) + 2*p1k*(-p2k + u) + p2k*(s + t + u)) - 
     -            k1k*(k2k*p2k*(4*p1k - s)*(2*p1k - s - t) - 
     -               2*p1k2*p2k*u + p2k*(2*p2k - s - t)*(s + t)*u - 
     -               k2k*(4*p1k2 - 2*p1k*(s + t) + p2k*(s + 2*t))*u + 
     -               p1k*(2*p2k2*t + 2*p2k*(2*s + t)*u + t*u2))) + 
     -         gqm2*(-2*k1k**3*p2k*s + 
     -            2*k1k2*(p1k*s*(p2k - t) - p1k*(p2k + t)*u + 
     -               p2k*t*(-p2k + u)) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (t*(2*p2k - s - u) + p1k*(-2*p2k + s + t + u)) + 
     -            k1k*(p2k*t2*u + 2*p1k2*(s*t + (p2k + t)*u) - 
     -               p1k*t*(2*p2k2 + (s + u)**2 - 2*p2k*(2*s + u)) + 
     -               k2k*(2*p2k*t*(-2*p2k + s + u) + 
     -                  p1k*(8*p2k2 - 6*p2k*s + s2 - s*t - 4*p2k*u + 
     -                     s*u - 2*t*u))))) + 
     -      glm2*(gqm2*(2*k1k**3*p1k*s + 
     -            2*k1k2*(p1k*p2k*(-s + t) + 
     -               (p1k2 - p1k*t + p2k*(s + t))*u) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (-((s + t)*u) + 2*p1k*(-p2k + u) + p2k*(s + t + u)) - 
     -            k1k*(k2k*p2k*(4*p1k - s)*(2*p1k - s - t) - 
     -               2*p1k2*p2k*u + p2k*(2*p2k - s - t)*(s + t)*u - 
     -               k2k*(4*p1k2 - 2*p1k*(s + t) + p2k*(s + 2*t))*u + 
     -               p1k*(2*p2k2*t + 2*p2k*(2*s + t)*u + t*u2))) + 
     -         gqp2*(-2*k1k**3*p2k*s + 
     -            2*k1k2*(p1k*s*(p2k - t) - p1k*(p2k + t)*u + 
     -               p2k*t*(-p2k + u)) - 
     -            k2k*(p2k*t - p1k*u)*
     -             (t*(2*p2k - s - u) + p1k*(-2*p2k + s + t + u)) + 
     -            k1k*(p2k*t2*u + 2*p1k2*(s*t + (p2k + t)*u) - 
     -               p1k*t*(2*p2k2 + (s + u)**2 - 2*p2k*(2*s + u)) + 
     -               k2k*(2*p2k*t*(-2*p2k + s + u) + 
     -                  p1k*(8*p2k2 - 6*p2k*s + s2 - s*t - 4*p2k*u + 
     -                     s*u - 2*t*u)))))))/(k1k*k2k*p1k*p2k)+
     -         (4*densp*denspg*qq**3*(glp*
     -       (gqp*((2*p1k - s - t)*(2*k1k*p1k - k1k*s - p2k*t) + 
     -            (2*p1k2 - 2*k1k*p2k + 2*p1k*p2k + 2*p2k2 + k1k*s - 
     -               3*p1k*s - 2*p2k*s + s2 - (p1k + p2k - s)*t)*u + 
     -            p1k*u2) + gqm*
     -          (2*p1k2*t + p1k*
     -             (-2*s*t + 2*p2k*(t - u) + s*u - t*u + u2) + 
     -            t*(2*p2k2 + p2k*(-3*s + t - u) + s*(s + u)) + 
     -         k1k*(4*p2k2 - 2*p1k*t - 2*p2k*(2*s + u) + s*(s + t + u)))
     -       ) + glm*(gqm*((2*p1k - s - t)*(2*k1k*p1k - k1k*s - p2k*t) + 
     -            (2*p1k2 - 2*k1k*p2k + 2*p1k*p2k + 2*p2k2 + k1k*s - 
     -               3*p1k*s - 2*p2k*s + s2 - (p1k + p2k - s)*t)*u + 
     -            p1k*u2) + gqp*
     -          (2*p1k2*t + p1k*
     -             (-2*s*t + 2*p2k*(t - u) + s*u - t*u + u2) + 
     -            t*(2*p2k2 + p2k*(-3*s + t - u) + s*(s + u)) + 
     -         k1k*(4*p2k2 - 2*p1k*t - 2*p2k*(2*s + u) + s*(s + t + u)))
     -         )))/(p1k*p2k) +
     -         (2*ml4*(-(dens2*p1k*p2k*
     -         (2*(glmd*glp + glm*glpd)*(gqm2 + gqp2)*(k1k2 + k2k2)*s + 
     -           glp2*(gqp2*(2*k1k**3 + 
     -                 2*k2k2*(-2*p1k - 2*p2k + s + t + u) + 
     -                 2*k1k*k2k*(-2*p1k - 2*p2k + s + 2*(t + u)) + 
     -                 k1k2*(2*p1k - 2*p2k + 5*(s + t + u))) + 
     -           gqm2*(2*k1k**3 + 2*k2k2*(-2*p1k - 2*p2k + s + t + u) + 
     -                 2*k1k*k2k*(-2*p1k - 2*p2k + s + 2*(t + u)) + 
     -                 k1k2*(-2*p1k + 2*p2k + 5*(s + t + u)))) + 
     -           glm2*(gqm2*(2*k1k**3 + 
     -                 2*k2k2*(-2*p1k - 2*p2k + s + t + u) + 
     -                 2*k1k*k2k*(-2*p1k - 2*p2k + s + 2*(t + u)) + 
     -                 k1k2*(2*p1k - 2*p2k + 5*(s + t + u))) + 
     -           gqp2*(2*k1k**3 + 2*k2k2*(-2*p1k - 2*p2k + s + t + u) + 
     -                 2*k1k*k2k*(-2*p1k - 2*p2k + s + 2*(t + u)) + 
     -                 k1k2*(-2*p1k + 2*p2k + 5*(s + t + u)))))) + 
     -      qq*(densd*(densg*p1k*p2k*
     -             (glpd*(gqpd*
     -                   (2*k1k**3 + 
     -                     2*k2k2*(-2*p1k - 2*p2k + 2*s + t + u) + 
     -                     2*k1k*k2k*(-2*p1k - 2*p2k + s + 2*(t + u)) + 
     -                     k1k2*(2*p1k - 2*p2k + 7*s + 5*(t + u))) + 
     -                  gqmd*(2*k1k**3 + 
     -                     2*k2k2*(-2*p1k - 2*p2k + 2*s + t + u) + 
     -                     2*k1k*k2k*(-2*p1k - 2*p2k + s + 2*(t + u)) + 
     -                     k1k2*(-2*p1k + 2*p2k + 7*s + 5*(t + u)))) + 
     -               glmd*(gqmd*
     -                   (2*k1k**3 + 
     -                     2*k2k2*(-2*p1k - 2*p2k + 2*s + t + u) + 
     -                     2*k1k*k2k*(-2*p1k - 2*p2k + s + 2*(t + u)) + 
     -                     k1k2*(2*p1k - 2*p2k + 7*s + 5*(t + u))) + 
     -                  gqpd*(2*k1k**3 + 
     -                     2*k2k2*(-2*p1k - 2*p2k + 2*s + t + u) + 
     -                     2*k1k*k2k*(-2*p1k - 2*p2k + s + 2*(t + u)) + 
     -                     k1k2*(-2*p1k + 2*p2k + 7*s + 5*(t + u))))) + 
     -            denspg*k1k*k2k*qq*
     -             (glmd*(gqpd*
     -                   (2*k1k2*(-p1k + p2k) + 
     -                     k1k*
     -                      (3*p1k2 + 2*k2k*(-p1k + p2k) + 
     -                        p1k*(2*p2k - 3*s - t - 2*u) + 
     -                        p2k*(-p2k + s + 2*t + u)) + 
     -                     k2k*
     -                      (3*p1k2 + p2k*(-3*p2k + 2*(s + t) + u) - 
     -                        p1k*(2*s + t + 2*u))) + 
     -                  gqmd*(2*k1k2*(-p1k + p2k) + 
     -                     k1k*
     -                      (p1k2 + 2*k2k*(-p1k + p2k) + 
     -                        p2k*(-3*p2k + 3*s + 2*t + u) - 
     -                        p1k*(2*p2k + s + t + 2*u)) + 
     -                     k2k*
     -                      (3*p1k2 + p2k*(-3*p2k + 2*(s + t) + u) - 
     -                        p1k*(2*s + t + 2*u)))) + 
     -               glpd*(gqmd*
     -                   (2*k1k2*(-p1k + p2k) + 
     -                     k1k*
     -                      (3*p1k2 + 2*k2k*(-p1k + p2k) + 
     -                        p1k*(2*p2k - 3*s - t - 2*u) + 
     -                        p2k*(-p2k + s + 2*t + u)) + 
     -                     k2k*
     -                      (3*p1k2 + p2k*(-3*p2k + 2*(s + t) + u) - 
     -                        p1k*(2*s + t + 2*u))) + 
     -                  gqpd*(2*k1k2*(-p1k + p2k) + 
     -                     k1k*
     -                      (p1k2 + 2*k2k*(-p1k + p2k) + 
     -                        p2k*(-3*p2k + 3*s + 2*t + u) - 
     -                        p1k*(2*p2k + s + t + 2*u)) + 
     -                     k2k*
     -                      (3*p1k2 + p2k*(-3*p2k + 2*(s + t) + u) - 
     -                        p1k*(2*s + t + 2*u))))) + 
     -            densp*k1k*k2k*
     -             (-((glmd*glp + glm*glpd)*(gqm2 + gqp2)*(k1k + k2k)*
     -                  (p1k - p2k)*(p1k + p2k - s)) + 
     -               glp2*(gqp2*
     -                   (2*k1k2*(p1k - p2k) + 
     -                     k1k*
     -                      (2*k2k*(p1k - p2k) + 
     -                        p2k*(2*p2k - 2*(s + t) - u) + 
     -                        p1k*(2*p2k + t + 2*u)) + 
     -                     k2k*
     -                      (-2*p1k2 + 2*p2k2 - p2k*(s + 2*t + u) + 
     -                        p1k*(s + t + 2*u))) + 
     -                  gqm2*(2*k1k2*(p1k - p2k) + 
     -                     k2k*
     -                      (-2*p1k2 + 2*p2k2 - p2k*(s + 2*t + u) + 
     -                        p1k*(s + t + 2*u)) + 
     -                     k1k*
     -                   (-2*p1k2 + 2*k2k*(p1k - p2k) - p2k*(2*t + u) + 
     -                        p1k*(-2*p2k + 2*s + t + 2*u)))) + 
     -               glm2*(gqm2*
     -                   (2*k1k2*(p1k - p2k) + 
     -                     k1k*
     -                      (2*k2k*(p1k - p2k) + 
     -                        p2k*(2*p2k - 2*(s + t) - u) + 
     -                        p1k*(2*p2k + t + 2*u)) + 
     -                     k2k*
     -                      (-2*p1k2 + 2*p2k2 - p2k*(s + 2*t + u) + 
     -                        p1k*(s + t + 2*u))) + 
     -                  gqp2*(2*k1k2*(p1k - p2k) + 
     -                     k2k*
     -                      (-2*p1k2 + 2*p2k2 - p2k*(s + 2*t + u) + 
     -                        p1k*(s + t + 2*u)) + 
     -                     k1k*
     -                   (-2*p1k2 + 2*k2k*(p1k - p2k) - p2k*(2*t + u) + 
     -                        p1k*(-2*p2k + 2*s + t + 2*u)))))) + 
     -         dens*(densg*p1k*p2k*
     -             (glp*(gqp*(2*k1k**3 + 
     -                     2*k2k2*(-2*p1k - 2*p2k + 2*s + t + u) + 
     -                     2*k1k*k2k*(-2*p1k - 2*p2k + s + 2*(t + u)) + 
     -                     k1k2*(2*p1k - 2*p2k + 7*s + 5*(t + u))) + 
     -                  gqm*(2*k1k**3 + 
     -                     2*k2k2*(-2*p1k - 2*p2k + 2*s + t + u) + 
     -                     2*k1k*k2k*(-2*p1k - 2*p2k + s + 2*(t + u)) + 
     -                     k1k2*(-2*p1k + 2*p2k + 7*s + 5*(t + u)))) + 
     -               glm*(gqm*(2*k1k**3 + 
     -                     2*k2k2*(-2*p1k - 2*p2k + 2*s + t + u) + 
     -                     2*k1k*k2k*(-2*p1k - 2*p2k + s + 2*(t + u)) + 
     -                     k1k2*(2*p1k - 2*p2k + 7*s + 5*(t + u))) + 
     -                  gqp*(2*k1k**3 + 
     -                     2*k2k2*(-2*p1k - 2*p2k + 2*s + t + u) + 
     -                     2*k1k*k2k*(-2*p1k - 2*p2k + s + 2*(t + u)) + 
     -                     k1k2*(-2*p1k + 2*p2k + 7*s + 5*(t + u))))) + 
     -            denspg*k1k*k2k*qq*
     -             (glm*(gqp*(2*k1k2*(-p1k + p2k) + 
     -                     k1k*
     -                      (3*p1k2 + 2*k2k*(-p1k + p2k) + 
     -                        p1k*(2*p2k - 3*s - t - 2*u) + 
     -                        p2k*(-p2k + s + 2*t + u)) + 
     -                     k2k*
     -                      (3*p1k2 + p2k*(-3*p2k + 2*(s + t) + u) - 
     -                        p1k*(2*s + t + 2*u))) + 
     -                  gqm*(2*k1k2*(-p1k + p2k) + 
     -                     k1k*
     -                      (p1k2 + 2*k2k*(-p1k + p2k) + 
     -                        p2k*(-3*p2k + 3*s + 2*t + u) - 
     -                        p1k*(2*p2k + s + t + 2*u)) + 
     -                     k2k*
     -                      (3*p1k2 + p2k*(-3*p2k + 2*(s + t) + u) - 
     -                        p1k*(2*s + t + 2*u)))) + 
     -               glp*(gqm*(2*k1k2*(-p1k + p2k) + 
     -                     k1k*
     -                      (3*p1k2 + 2*k2k*(-p1k + p2k) + 
     -                        p1k*(2*p2k - 3*s - t - 2*u) + 
     -                        p2k*(-p2k + s + 2*t + u)) + 
     -                     k2k*
     -                      (3*p1k2 + p2k*(-3*p2k + 2*(s + t) + u) - 
     -                        p1k*(2*s + t + 2*u))) + 
     -                  gqp*(2*k1k2*(-p1k + p2k) + 
     -                     k1k*
     -                      (p1k2 + 2*k2k*(-p1k + p2k) + 
     -                        p2k*(-3*p2k + 3*s + 2*t + u) - 
     -                        p1k*(2*p2k + s + t + 2*u)) + 
     -                     k2k*
     -                      (3*p1k2 + p2k*(-3*p2k + 2*(s + t) + u) - 
     -                        p1k*(2*s + t + 2*u))))) + 
     -            denspd*k1k*k2k*
     -             (-((glmd*glp + glm*glpd)*(gqm2 + gqp2)*(k1k + k2k)*
     -                  (p1k - p2k)*(p1k + p2k - s)) + 
     -               glp2*(gqp2*
     -                   (2*k1k2*(p1k - p2k) + 
     -                     k1k*
     -                      (2*k2k*(p1k - p2k) + 
     -                        p2k*(2*p2k - 2*(s + t) - u) + 
     -                        p1k*(2*p2k + t + 2*u)) + 
     -                     k2k*
     -                      (-2*p1k2 + 2*p2k2 - p2k*(s + 2*t + u) + 
     -                        p1k*(s + t + 2*u))) + 
     -                  gqm2*(2*k1k2*(p1k - p2k) + 
     -                     k2k*
     -                      (-2*p1k2 + 2*p2k2 - p2k*(s + 2*t + u) + 
     -                        p1k*(s + t + 2*u)) + 
     -                     k1k*
     -                   (-2*p1k2 + 2*k2k*(p1k - p2k) - p2k*(2*t + u) + 
     -                        p1k*(-2*p2k + 2*s + t + 2*u)))) + 
     -               glm2*(gqm2*
     -                   (2*k1k2*(p1k - p2k) + 
     -                     k1k*
     -                      (2*k2k*(p1k - p2k) + 
     -                        p2k*(2*p2k - 2*(s + t) - u) + 
     -                        p1k*(2*p2k + t + 2*u)) + 
     -                     k2k*
     -                      (-2*p1k2 + 2*p2k2 - p2k*(s + 2*t + u) + 
     -                        p1k*(s + t + 2*u))) + 
     -                  gqp2*(2*k1k2*(p1k - p2k) + 
     -                     k2k*
     -                      (-2*p1k2 + 2*p2k2 - p2k*(s + 2*t + u) + 
     -                        p1k*(s + t + 2*u)) + 
     -                     k1k*
     -                   (-2*p1k2 + 2*k2k*(p1k - p2k) - p2k*(2*t + u) + 
     -                        p1k*(-2*p2k + 2*s + t + 2*u)))))) - 
     -         qq*(-2*k1k2*k2k2*
     -             (-(densp2*(glm2 + glp2)*(gqm2 + gqp2)) + 
     -               denspg*qq*
     -                (densp*(glm + glp)*(gqm + gqp) + 
     -              denspd*(glmd + glpd)*(gqmd + gqpd) - 4*denspg*qq))*s
     -             + 4*densg**2*p1k*p2k*
     -             (2*k1k**3 + 2*k2k2*(-2*p1k - 2*p2k + 2*s + t + u) +
     -               2*k1k*k2k*(-2*p1k - 2*p2k + s + 2*(t + u)) + 
     -               k1k2*(7*s + 5*(t + u))) + 
     -            densg*k1k*k2k*
     -             (8*denspg*qq*
     -                (2*k1k2*(-p1k + p2k) + 
     -                  k2k*(3*p1k2 + p2k*(-3*p2k + 2*(s + t) + u) - 
     -                     p1k*(2*s + t + 2*u)) + 
     -                  k1k*(2*p1k2 + 2*k2k*(-p1k + p2k) + 
     -                p2k*(-2*p2k + 2*(s + t) + u) - p1k*(2*s + t + 2*u)
     -                     )) + 
     -               densp*(glp*
     -                   (gqp*(2*k1k2*(p1k - p2k) + 
     -                        k1k*
     -                         (-p1k2 + 2*k2k*(p1k - p2k) + 
     -                          p2k*(3*p2k - 3*s - 2*t - u) + 
     -                          p1k*(2*p2k + s + t + 2*u)) + 
     -                        k2k*
     -                         (-3*p1k2 + p2k*(3*p2k - 2*(s + t) - u) + 
     -                          p1k*(2*s + t + 2*u))) + 
     -                     gqm*
     -                      (2*k1k2*(p1k - p2k) + 
     -                        k2k*
     -                         (-3*p1k2 + p2k*(3*p2k - 2*(s + t) - u) + 
     -                          p1k*(2*s + t + 2*u)) + 
     -                        k1k*
     -                         (-3*p1k2 + 2*k2k*(p1k - p2k) + 
     -                          p2k*(p2k - s - 2*t - u) + 
     -                          p1k*(-2*p2k + 3*s + t + 2*u)))) + 
     -                  glm*(gqm*
     -                      (2*k1k2*(p1k - p2k) + 
     -                        k1k*
     -                         (-p1k2 + 2*k2k*(p1k - p2k) + 
     -                          p2k*(3*p2k - 3*s - 2*t - u) + 
     -                          p1k*(2*p2k + s + t + 2*u)) + 
     -                        k2k*
     -                         (-3*p1k2 + p2k*(3*p2k - 2*(s + t) - u) + 
     -                          p1k*(2*s + t + 2*u))) + 
     -                     gqp*
     -                      (2*k1k2*(p1k - p2k) + 
     -                        k2k*
     -                         (-3*p1k2 + p2k*(3*p2k - 2*(s + t) - u) + 
     -                          p1k*(2*s + t + 2*u)) + 
     -                        k1k*
     -                         (-3*p1k2 + 2*k2k*(p1k - p2k) + 
     -                          p2k*(p2k - s - 2*t - u) + 
     -                          p1k*(-2*p2k + 3*s + t + 2*u))))) + 
     -               denspd*(glmd*
     -                   (gqpd*
     -                      (2*k1k2*(p1k - p2k + mz2dm1*p1k*p2k) + 
     -                        k2k*
     -                         (-3*p1k2 + p2k*(3*p2k - 2*(s + t) - u) + 
     -                          p1k*(2*s + t + 2*u)) + 
     -                        k1k*
     -                         (2*k2k*(p1k - p2k) - 
     -                          p1k2*(3 + 2*mz2dm1*p2k) + 
     -                          p2k*(p2k - s - 2*t - u) + 
     -                          p1k*
     -                          (3*s + t + 2*u + 
     -                     p2k*(-2 + mz2dm1*(-2*p2k + 5*s + 3*(t + u))))
     -                          )) + 
     -                     gqmd*
     -                      (2*k1k2*(p1k - p2k + mz2dm1*p1k*p2k) + 
     -                        k2k*
     -                         (-3*p1k2 + p2k*(3*p2k - 2*(s + t) - u) + 
     -                          p1k*(2*s + t + 2*u)) + 
     -                        k1k*
     -                         (2*k2k*(p1k - p2k) - 
     -                          p1k2*(1 + 2*mz2dm1*p2k) + 
     -                          p2k*(3*p2k - 3*s - 2*t - u) + 
     -                          p1k*
     -                          (s + t + 2*u + 
     -                    p2k*(2 + mz2dm1*(-2*p2k + 5*s + 3*(t + u)))))
     -                        )) + 
     -                  glpd*(gqmd*
     -                      (2*k1k2*(p1k - p2k + mz2dm1*p1k*p2k) + 
     -                        k2k*
     -                         (-3*p1k2 + p2k*(3*p2k - 2*(s + t) - u) + 
     -                          p1k*(2*s + t + 2*u)) + 
     -                        k1k*
     -                         (2*k2k*(p1k - p2k) - 
     -                          p1k2*(3 + 2*mz2dm1*p2k) + 
     -                          p2k*(p2k - s - 2*t - u) + 
     -                          p1k*
     -                          (3*s + t + 2*u + 
     -                     p2k*(-2 + mz2dm1*(-2*p2k + 5*s + 3*(t + u))))
     -                          )) + 
     -                     gqpd*
     -                      (2*k1k2*(p1k - p2k + mz2dm1*p1k*p2k) + 
     -                        k2k*
     -                         (-3*p1k2 + p2k*(3*p2k - 2*(s + t) - u) + 
     -                          p1k*(2*s + t + 2*u)) + 
     -                        k1k*
     -                         (2*k2k*(p1k - p2k) - 
     -                          p1k2*(1 + 2*mz2dm1*p2k) + 
     -                          p2k*(3*p2k - 3*s - 2*t - u) + 
     -                          p1k*
     -                          (s + t + 2*u + 
     -                     p2k*(2 + mz2dm1*(-2*p2k + 5*s + 3*(t + u)))))
     -                        ))))))))/(k1k2*k2k2*p1k*p2k) +
     -     (8*densg**2*ml2*qq2*(2*k1k**3*(2*s + t + u) - 
     -      2*k1k*k2k*(-p1k2 + 2*p1k*p2k - p2k2 - p1k*s - p2k*s + 
     -      k2k*(p1k + p2k + s) + 2*p2k*t + s*t + (2*p1k + s - 2*t)*u) + 
     -      k2k2*(8*p1k*p2k - 2*p1k*s - 2*p2k*s - 4*p2k*t + s*t + 
     -         (-4*p1k + s + 2*t)*u) + 
     -      k1k2*(-3*p1k*s - 3*p2k*s + s2 - 4*k2k*(p1k + p2k + 2*s) - 
     -       p1k*t + p2k*t + 4*s*t + t2 + (p1k - p2k + 4*s + 6*t)*u + 
     -         u2)))/(k1k2*k2k2) +
     -       (2*dens2*ml2*(-2*(glmd*glp + glm*glpd)*(gqm2 + gqp2)*k1k*
     -       (k2k2*s + k1k*(-k1k + p1k + p2k)*s + 
     -         k2k*(p1k*(4*p2k - s) + s*(2*k1k - p2k + t + u))) + 
     -      glm2*(gqp2*(2*k2k2*(2*p1k - t)*(2*p2k - s - u) + 
     -            2*k1k**3*(s + t + u) + 
     -            4*k1k*k2k*(-(k2k*p1k) + p1k2 + p1k*(p2k - u) + 
     -               t*(-p2k + u)) + 
     -            k1k2*(s2 - 4*k2k*(p1k + p2k + s) + 4*p2k*t + 6*s*t + 
     -               t2 + 2*s*u + 6*t*u + u2 - 2*p1k*(s + t + u))) + 
     -         gqm2*(2*k2k2*(2*p1k - s - t)*(2*p2k - u) + 
     -            2*k1k**3*(s + t + u) + 
     -            4*k1k*k2k*(-(k2k*p2k) + p2k*(p1k + p2k - t) + 
     -               (-p1k + t)*u) + 
     -            k1k2*(-4*k2k*(p1k + p2k + s) + (s + t)**2 + 4*p1k*u + 
     -               6*(s + t)*u + u2 - 2*p2k*(s + t + u)))) + 
     -      glp2*(gqm2*(2*k2k2*(2*p1k - t)*(2*p2k - s - u) + 
     -            2*k1k**3*(s + t + u) + 
     -            4*k1k*k2k*(-(k2k*p1k) + p1k2 + p1k*(p2k - u) + 
     -               t*(-p2k + u)) + 
     -            k1k2*(s2 - 4*k2k*(p1k + p2k + s) + 4*p2k*t + 6*s*t + 
     -               t2 + 2*s*u + 6*t*u + u2 - 2*p1k*(s + t + u))) + 
     -         gqp2*(2*k2k2*(2*p1k - s - t)*(2*p2k - u) + 
     -            2*k1k**3*(s + t + u) + 
     -            4*k1k*k2k*(-(k2k*p2k) + p2k*(p1k + p2k - t) + 
     -               (-p1k + t)*u) + 
     -          k1k2*(-4*k2k*(p1k + p2k + s) + (s + t)**2 + 4*p1k*u + 
     -               6*(s + t)*u + u2 - 2*p2k*(s + t + u))))))/
     -  (k1k2*k2k2) +
     -         (-2*densd*densg*ml2*qq*
     -    (glmd*(gqpd*(2*k2k2*(2*p1k - t)*(2*p2k - s - u) + 
     -            2*k1k**3*(2*s + t + u) - 
     -            2*k1k*k2k*(-2*p1k2 + 2*p1k*p2k - p1k*s - p2k*s + 
     -            k2k*(2*p1k + s) + 2*p2k*t + s*t + (2*p1k + s - 2*t)*u)
     -             + k1k2*(-2*p2k*s + s2 - 4*k2k*(p1k + p2k + 2*s) +
     -               4*p2k*t + 6*s*t + t2 + 2*s*u + 6*t*u + u2 - 
     -               2*p1k*(2*s + t + u))) + 
     -         gqmd*(2*k2k2*(2*p1k - s - t)*(2*p2k - u) + 
     -            2*k1k**3*(2*s + t + u) + 
     -            k1k2*(-2*p1k*s - 4*p2k*s + s2 - 
     -               4*k2k*(p1k + p2k + 2*s) - 2*p2k*t + 2*s*t + t2 + 
     -               4*p1k*u - 2*p2k*u + 6*s*u + 6*t*u + u2) - 
     -            2*k1k*k2k*(k2k*(2*p2k + s) - (2*p2k + s)*(p2k - t) + 
     -               (s - 2*t)*u + p1k*(2*p2k - s + 2*u)))) + 
     -      glpd*(gqmd*(2*k2k2*(2*p1k - t)*(2*p2k - s - u) + 
     -            2*k1k**3*(2*s + t + u) - 
     -            2*k1k*k2k*(-2*p1k2 + 2*p1k*p2k - p1k*s - p2k*s + 
     -            k2k*(2*p1k + s) + 2*p2k*t + s*t + (2*p1k + s - 2*t)*u)
     -             + k1k2*(-2*p2k*s + s2 - 4*k2k*(p1k + p2k + 2*s) + 
     -               4*p2k*t + 6*s*t + t2 + 2*s*u + 6*t*u + u2 - 
     -               2*p1k*(2*s + t + u))) + 
     -         gqpd*(2*k2k2*(2*p1k - s - t)*(2*p2k - u) + 
     -            2*k1k**3*(2*s + t + u) + 
     -            k1k2*(-2*p1k*s - 4*p2k*s + s2 - 
     -               4*k2k*(p1k + p2k + 2*s) - 2*p2k*t + 2*s*t + t2 + 
     -               4*p1k*u - 2*p2k*u + 6*s*u + 6*t*u + u2) - 
     -            2*k1k*k2k*(k2k*(2*p2k + s) - (2*p2k + s)*(p2k - t) + 
     -             (s - 2*t)*u + p1k*(2*p2k - s + 2*u))))))/(k1k2*k2k2)+
     -        (-2*dens*densg*ml2*qq*(glm*
     -       (gqp*(2*k2k2*(2*p1k - t)*(2*p2k - s - u) + 
     -            2*k1k**3*(2*s + t + u) - 
     -            2*k1k*k2k*(-2*p1k2 + 2*p1k*p2k - p1k*s - p2k*s + 
     -            k2k*(2*p1k + s) + 2*p2k*t + s*t + (2*p1k + s - 2*t)*u)
     -             + k1k2*(-2*p2k*s + s2 - 4*k2k*(p1k + p2k + 2*s) +
     -               4*p2k*t + 6*s*t + t2 + 2*s*u + 6*t*u + u2 - 
     -               2*p1k*(2*s + t + u))) + 
     -         gqm*(2*k2k2*(2*p1k - s - t)*(2*p2k - u) + 
     -            2*k1k**3*(2*s + t + u) + 
     -            k1k2*(-2*p1k*s - 4*p2k*s + s2 - 
     -               4*k2k*(p1k + p2k + 2*s) - 2*p2k*t + 2*s*t + t2 + 
     -               4*p1k*u - 2*p2k*u + 6*s*u + 6*t*u + u2) - 
     -            2*k1k*k2k*(k2k*(2*p2k + s) - (2*p2k + s)*(p2k - t) + 
     -               (s - 2*t)*u + p1k*(2*p2k - s + 2*u)))) + 
     -      glp*(gqm*(2*k2k2*(2*p1k - t)*(2*p2k - s - u) + 
     -            2*k1k**3*(2*s + t + u) - 
     -            2*k1k*k2k*(-2*p1k2 + 2*p1k*p2k - p1k*s - p2k*s + 
     -            k2k*(2*p1k + s) + 2*p2k*t + s*t + (2*p1k + s - 2*t)*u)
     -             + k1k2*(-2*p2k*s + s2 - 4*k2k*(p1k + p2k + 2*s) +
     -               4*p2k*t + 6*s*t + t2 + 2*s*u + 6*t*u + u2 - 
     -               2*p1k*(2*s + t + u))) + 
     -         gqp*(2*k2k2*(2*p1k - s - t)*(2*p2k - u) + 
     -            2*k1k**3*(2*s + t + u) + 
     -            k1k2*(-2*p1k*s - 4*p2k*s + s2 - 
     -               4*k2k*(p1k + p2k + 2*s) - 2*p2k*t + 2*s*t + t2 + 
     -               4*p1k*u - 2*p2k*u + 6*s*u + 6*t*u + u2) - 
     -            2*k1k*k2k*(k2k*(2*p2k + s) - (2*p2k + s)*(p2k - t) + 
     -            (s - 2*t)*u + p1k*(2*p2k - s + 2*u))))))/(k1k2*k2k2)+
     -         (8*densg*denspg*ml2*qq**3*
     -    (2*k1k2*(p1k - p2k)*(p1k + p2k - 2*(s + t + u)) + 
     -      k2k*(p1k2*(-12*p2k + s + t + 9*u) + 
     -         p2k*(s*(5*t + u) - p2k*(s + 9*t + u) + 2*t*(t + 2*u)) + 
     -         p1k*(12*p2k2 - 2*u*(2*t + u) - s*(t + 5*u))) + 
     -      k1k*(4*k2k*(p1k - p2k)*(p1k + p2k - s - t - u) + 
     -         2*p1k2*(2*p2k + 3*s + t + 2*u) + 
     -         p2k*(3*s2 + 2*s*(2*t + u) - 2*p2k*(3*s + 2*t + u) + 
     -            2*t*(t + 2*u)) - 
     -         p1k*(3*s2 + 2*
     -         (2*p2k2 + s*t + 2*(s + t)*u + u2 + p2k*(-t + u))))))/
     -  (k1k*k2k*p1k*p2k) + 
     -        (2*densd*denspg*ml2*qq2*
     -    (-(glpd*(gqpd*(2*k1k2 + 5*k2k*s + k1k*(4*k2k - t - 3*u)) + 
     -           gqmd*(-2*k1k2 - 5*k2k*s + k1k*(-4*k2k + 3*t + u)))) - 
     -      glmd*(gqmd*(2*k1k2 + 5*k2k*s + k1k*(4*k2k - t - 3*u)) + 
     -         gqpd*(-2*k1k2 - 5*k2k*s + k1k*(-4*k2k + 3*t + u))) + 
     -      (p2k*(glmd*(2*gqpd*k1k*(k1k + 2*k2k + s) + 
     -              gqpd*(k1k + 5*k2k)*t + 
     -              gqmd*(k2k*(s + 4*t + u) + k1k*(4*s + 3*t + 2*u))) + 
     -       glpd*(2*gqmd*k1k*(k1k + 2*k2k + s) + gqmd*(k1k + 5*k2k)*t + 
     -           gqpd*(k2k*(s + 4*t + u) + k1k*(4*s + 3*t + 2*u)))))/p1k
     -       - (glpd*(2*gqpd*(k1k2*p1k + k2k*p2k*(-5*p1k + p2k) + 
     -               k1k*(2*k2k*p1k - p2k*(p1k + 3*p2k) + p1k*s)) + 
     -            gqpd*(k1k + 5*k2k)*p1k*u + 
     -            gqmd*(k1k*(2*p2k2 + p1k*(6*p2k + 4*s + 2*t + 3*u)) + 
     -               k2k*(10*p2k2 + p1k*(-2*p2k + s + t + 4*u)))) + 
     -         glmd*(2*gqmd*(k1k2*p1k + k2k*p2k*(-5*p1k + p2k) + 
     -               k1k*(2*k2k*p1k - p2k*(p1k + 3*p2k) + p1k*s)) + 
     -            gqmd*(k1k + 5*k2k)*p1k*u + 
     -            gqpd*(k1k*(2*p2k2 + p1k*(6*p2k + 4*s + 2*t + 3*u)) + 
     -            k2k*(10*p2k2 + p1k*(-2*p2k + s + t + 4*u)))))/p2k + 
     -      (glpd*(gqpd*(k1k*s*(k1k + 3*k2k + s) + 2*k1k*(k1k + k2k)*t + 
     -               (2*k1k2 + 2*k1k*k2k + k1k*s + 3*k2k*s + 
     -                  2*(k1k + k2k)*t)*u + (k1k + k2k)*u2) + 
     -            gqmd*(k2k*(s*t + 2*(s + t)*u + u2) + 
     -               k1k2*(3*s + 2*(t + u)) + 
     -           k1k*((s + u)*(2*(s + t) + u) + k2k*(s + 2*(t + u))))) + 
     -      glmd*(gqmd*(k1k*s*(k1k + 3*k2k + s) + 2*k1k*(k1k + k2k)*t + 
     -               (2*k1k2 + 2*k1k*k2k + k1k*s + 3*k2k*s + 
     -                  2*(k1k + k2k)*t)*u + (k1k + k2k)*u2) + 
     -            gqpd*(k2k*(s*t + 2*(s + t)*u + u2) + 
     -               k1k2*(3*s + 2*(t + u)) + 
     -           k1k*((s + u)*(2*(s + t) + u) + k2k*(s + 2*(t + u))))))/
     -       p2k - (glmd*(gqpd*
     -             (k2k*t*(3*s + t + 2*u) + 
     -               k1k*(3*k2k*s + s2 + 2*k2k*t + s*t + t2 + 
     -                  2*(k2k + t)*u) + k1k2*(s + 2*(t + u))) + 
     -            gqmd*(k1k2*(3*s + 2*(t + u)) + 
     -               k2k*(s*(2*t + u) + t*(t + 2*u)) + 
     -          k1k*((s + t)*(2*s + t + 2*u) + k2k*(s + 2*(t + u))))) + 
     -         glpd*(gqmd*(k2k*t*(3*s + t + 2*u) + 
     -               k1k*(3*k2k*s + s2 + 2*k2k*t + s*t + t2 + 
     -                  2*(k2k + t)*u) + k1k2*(s + 2*(t + u))) + 
     -            gqpd*(k1k2*(3*s + 2*(t + u)) + 
     -               k2k*(s*(2*t + u) + t*(t + 2*u)) + 
     -          k1k*((s + t)*(2*s + t + 2*u) + k2k*(s + 2*(t + u))))))/
     -       p1k))/(k1k*k2k) +
     -        (2*dens*denspg*ml2*qq2*
     -   ((p2k*(glm*(2*gqp*k1k*(k1k + 2*k2k + s) + gqp*(k1k + 5*k2k)*t + 
     -              gqm*(k2k*(s + 4*t + u) + k1k*(4*s + 3*t + 2*u))) + 
     -          glp*(2*gqm*k1k*(k1k + 2*k2k + s) + gqm*(k1k + 5*k2k)*t + 
     -          gqp*(k2k*(s + 4*t + u) + k1k*(4*s + 3*t + 2*u)))))/p1k - 
     -      (glp*(gqp*(2*k1k2*(p1k + p2k) + 
     -               k2k*(p2k*(2*p2k + 5*s) + 5*p1k*(-2*p2k + u)) + 
     -               k1k*(4*k2k*(p1k + p2k) + p1k*(-2*p2k + 2*s + u) - 
     -                  p2k*(6*p2k + t + 3*u))) + 
     -            gqm*(-2*k1k2*p2k + 
     -               k1k*(-4*k2k*p2k + p2k*(2*p2k + 3*t + u) + 
     -                  p1k*(6*p2k + 4*s + 2*t + 3*u)) + 
     -            k2k*(5*p2k*(2*p2k - s) + p1k*(-2*p2k + s + t + 4*u))))
     -          + glm*(gqm*(2*k1k2*(p1k + p2k) + 
     -               k2k*(p2k*(2*p2k + 5*s) + 5*p1k*(-2*p2k + u)) + 
     -               k1k*(4*k2k*(p1k + p2k) + p1k*(-2*p2k + 2*s + u) - 
     -                  p2k*(6*p2k + t + 3*u))) + 
     -            gqp*(-2*k1k2*p2k + 
     -               k1k*(-4*k2k*p2k + p2k*(2*p2k + 3*t + u) + 
     -                  p1k*(6*p2k + 4*s + 2*t + 3*u)) + 
     -          k2k*(5*p2k*(2*p2k - s) + p1k*(-2*p2k + s + t + 4*u)))))/
     -       p2k + (glp*(gqp*(k1k*s*(k1k + 3*k2k + s) + 
     -               2*k1k*(k1k + k2k)*t + 
     -               (2*k1k2 + 2*k1k*k2k + k1k*s + 3*k2k*s + 
     -                  2*(k1k + k2k)*t)*u + (k1k + k2k)*u2) + 
     -            gqm*(k2k*(s*t + 2*(s + t)*u + u2) + 
     -               k1k2*(3*s + 2*(t + u)) + 
     -           k1k*((s + u)*(2*(s + t) + u) + k2k*(s + 2*(t + u))))) + 
     -         glm*(gqm*(k1k*s*(k1k + 3*k2k + s) + 2*k1k*(k1k + k2k)*t + 
     -               (2*k1k2 + 2*k1k*k2k + k1k*s + 3*k2k*s + 
     -                  2*(k1k + k2k)*t)*u + (k1k + k2k)*u2) + 
     -            gqp*(k2k*(s*t + 2*(s + t)*u + u2) + 
     -               k1k2*(3*s + 2*(t + u)) + 
     -           k1k*((s + u)*(2*(s + t) + u) + k2k*(s + 2*(t + u))))))/
     -       p2k - (glm*(gqp*(k2k*t*(3*s + t + 2*u) + 
     -               k1k*(3*k2k*s + s2 + 2*k2k*t + s*t + t2 + 
     -                  2*(k2k + t)*u) + k1k2*(s + 2*(t + u))) + 
     -            gqm*(k1k2*(3*s + 2*(t + u)) + 
     -               k2k*(s*(2*t + u) + t*(t + 2*u)) + 
     -           k1k*((s + t)*(2*s + t + 2*u) + k2k*(s + 2*(t + u))))) + 
     -         glp*(gqm*(k2k*t*(3*s + t + 2*u) + 
     -               k1k*(3*k2k*s + s2 + 2*k2k*t + s*t + t2 + 
     -                  2*(k2k + t)*u) + k1k2*(s + 2*(t + u))) + 
     -            gqp*(k1k2*(3*s + 2*(t + u)) + 
     -               k2k*(s*(2*t + u) + t*(t + 2*u)) + 
     -           k1k*((s + t)*(2*s + t + 2*u) + k2k*(s + 2*(t + u))))))/
     -       p1k))/(k1k*k2k)+
     -        (16*denspg**2*ml2*qq**4*
     -    (4*p1k2 + 4*p2k2 - 2*k1k*(p1k + p2k - s) - 4*p1k*s - 4*p2k*s + 
     -      2*s2 - p1k*t + p2k*t + s*t + (p1k - p2k + s)*u))/(p1k*p2k)+
     -         (4*densp2*ml2*qq2*((glmd*glp + glm*glpd)*(gqm2 + gqp2)*
     -       (2*(p1k2 + p2k2) - 2*(p1k + p2k)*s + s2) + 
     -      glm2*(gqp2*(2*p1k2 + 2*p2k2 - 2*k1k*(p1k + p2k - s) - 
     -            p1k*s - 3*p2k*s + s2 - p1k*t + p2k*t + s*t + 
     -            (p1k - p2k + s)*u) + 
     -         gqm2*(2*p1k2 + 2*p2k2 - 2*k1k*(p1k + p2k - s) - 3*p1k*s - 
     -          p2k*s + s2 - p1k*t + p2k*t + s*t + (p1k - p2k + s)*u)) + 
     -      glp2*(gqm2*(2*p1k2 + 2*p2k2 - 2*k1k*(p1k + p2k - s) - 
     -            p1k*s - 3*p2k*s + s2 - p1k*t + p2k*t + s*t + 
     -            (p1k - p2k + s)*u) + 
     -         gqp2*(2*p1k2 + 2*p2k2 - 2*k1k*(p1k + p2k - s) - 3*p1k*s - 
     -         p2k*s + s2 - p1k*t + p2k*t + s*t + (p1k - p2k + s)*u))))/
     -  (p1k*p2k)+
     -         (2*densg*denspd*ml2*qq2*
     -((p2k*(glmd*(2*gqpd*k1k*(k1k + 2*k2k + s) + gqpd*(k1k + 5*k2k)*t + 
     -              gqmd*(k2k*(s + 4*t + u) + k1k*(4*s + 3*t + 2*u))) + 
     -       glpd*(2*gqmd*k1k*(k1k + 2*k2k + s) + gqmd*(k1k + 5*k2k)*t + 
     -           gqpd*(k2k*(s + 4*t + u) + k1k*(4*s + 3*t + 2*u)))))/p1k
     -       - (p1k*(glpd*(2*gqpd*k1k*(k1k + 2*k2k + s) + 
     -              gqpd*(k1k + 5*k2k)*u + 
     -              gqmd*(k1k*(4*s + 2*t + 3*u) + k2k*(s + t + 4*u))) + 
     -       glmd*(2*gqmd*k1k*(k1k + 2*k2k + s) + gqmd*(k1k + 5*k2k)*u + 
     -           gqpd*(k1k*(4*s + 2*t + 3*u) + k2k*(s + t + 4*u)))))/p2k
     -       + glmd*(gqpd*(k2k*s + 2*k1k2*(1 + mz2dm1*(2*s + t + u)) + 
     -            k1k*(4*k2k + 2*mz2dm1*s2 - 3*t - u + 
     -               s*(-4 + 5*mz2dm1*t + 3*mz2dm1*u) + 
     -               mz2dm1*(t2 + 4*t*u + u2))) + 
     -         gqmd*(-(k2k*s) + 2*k1k2*(-1 + mz2dm1*(2*s + t + u)) + 
     -            k1k*(-4*k2k + 2*mz2dm1*s2 + t + 3*u + 
     -               s*(4 + 3*mz2dm1*t + 5*mz2dm1*u) + 
     -               mz2dm1*(t2 + 4*t*u + u2)))) + 
     -      glpd*(gqmd*(k2k*s + 2*k1k2*(1 + mz2dm1*(2*s + t + u)) + 
     -            k1k*(4*k2k + 2*mz2dm1*s2 - 3*t - u + 
     -               s*(-4 + 5*mz2dm1*t + 3*mz2dm1*u) + 
     -               mz2dm1*(t2 + 4*t*u + u2))) + 
     -         gqpd*(-(k2k*s) + 2*k1k2*(-1 + mz2dm1*(2*s + t + u)) + 
     -            k1k*(-4*k2k + 2*mz2dm1*s2 + t + 3*u + 
     -               s*(4 + 3*mz2dm1*t + 5*mz2dm1*u) + 
     -               mz2dm1*(t2 + 4*t*u + u2)))) - 
     -      2*p2k*(glmd*(gqpd*(3*k2k + k1k*(-1 + mz2dm1*(s + t + u))) + 
     -            gqmd*(3*k2k + k1k*(-1 + mz2dm1*(3*s + t + u)))) + 
     -         glpd*(gqmd*(3*k2k + k1k*(-1 + mz2dm1*(s + t + u))) + 
     -            gqpd*(3*k2k + k1k*(-1 + mz2dm1*(3*s + t + u))))) - 
     -      2*p1k*(glpd*(gqpd*(-3*k2k + k1k*(1 + mz2dm1*(s + t + u))) + 
     -            gqmd*(-3*k2k + k1k*(1 + mz2dm1*(3*s + t + u)))) + 
     -         glmd*(gqmd*(-3*k2k + k1k*(1 + mz2dm1*(s + t + u))) + 
     -            gqpd*(-3*k2k + k1k*(1 + mz2dm1*(3*s + t + u))))) + 
     -      (glpd*(gqpd*(k1k*s*(k1k + 3*k2k + s) + 2*k1k*(k1k + k2k)*t + 
     -               (2*k1k2 + 2*k1k*k2k + k1k*s + 3*k2k*s + 
     -                  2*(k1k + k2k)*t)*u + (k1k + k2k)*u2) + 
     -            gqmd*(k2k*(s*t + 2*(s + t)*u + u2) + 
     -               k1k2*(3*s + 2*(t + u)) + 
     -           k1k*((s + u)*(2*(s + t) + u) + k2k*(s + 2*(t + u))))) + 
     -       glmd*(gqmd*(k1k*s*(k1k + 3*k2k + s) + 2*k1k*(k1k + k2k)*t + 
     -               (2*k1k2 + 2*k1k*k2k + k1k*s + 3*k2k*s + 
     -                  2*(k1k + k2k)*t)*u + (k1k + k2k)*u2) + 
     -            gqpd*(k2k*(s*t + 2*(s + t)*u + u2) + 
     -               k1k2*(3*s + 2*(t + u)) + 
     -           k1k*((s + u)*(2*(s + t) + u) + k2k*(s + 2*(t + u))))))/
     -       p2k - (glmd*(gqpd*
     -             (k2k*t*(3*s + t + 2*u) + 
     -               k1k*(3*k2k*s + s2 + 2*k2k*t + s*t + t2 + 
     -                  2*(k2k + t)*u) + k1k2*(s + 2*(t + u))) + 
     -            gqmd*(k1k2*(3*s + 2*(t + u)) + 
     -               k2k*(s*(2*t + u) + t*(t + 2*u)) + 
     -           k1k*((s + t)*(2*s + t + 2*u) + k2k*(s + 2*(t + u))))) + 
     -         glpd*(gqmd*(k2k*t*(3*s + t + 2*u) + 
     -               k1k*(3*k2k*s + s2 + 2*k2k*t + s*t + t2 + 
     -                  2*(k2k + t)*u) + k1k2*(s + 2*(t + u))) + 
     -            gqpd*(k1k2*(3*s + 2*(t + u)) + 
     -               k2k*(s*(2*t + u) + t*(t + 2*u)) + 
     -           k1k*((s + t)*(2*s + t + 2*u) + k2k*(s + 2*(t + u))))))/
     -       p1k))/(k1k*k2k)+
     -
     -
     -        (2*dens*denspd*ml2*qq*(-2*
     -       (((glm2 - 2*glmd*glp - glp2)*gqm2 + 
     -            (-glm2 + glp2 - 2*glm*glpd)*gqp2)*k1k + 
     -         ((3*glm2 + glp2 + 2*glm*glpd)*gqm2 + 
     -            (glm2 + 2*glmd*glp + 3*glp2)*gqp2)*k2k)*p1k + 
     -      2*(((-glm2 + glp2 - 2*glm*glpd)*gqm2 + 
     -            (glm2 - 2*glmd*glp - glp2)*gqp2)*k1k + 
     -         ((glm2 + 2*glmd*glp + 3*glp2)*gqm2 + 
     -            (3*glm2 + glp2 + 2*glm*glpd)*gqp2)*k2k)*p2k + 
     -      glm2*(gqp2*(-2*k1k2 + 2*k1k*(-2*k2k + s + 2*t) + 
     -            k2k*(-3*s + t - u)) + 
     -    gqm2*(2*k1k2 + k1k*(4*k2k - 2*s - 4*u) + k2k*(3*s + t - u)))
     -        + glp2*(gqm2*(-2*k1k2 + 2*k1k*(-2*k2k + s + 2*t) + 
     -            k2k*(-3*s + t - u)) + 
     -    gqp2*(2*k1k2 + k1k*(4*k2k - 2*s - 4*u) + k2k*(3*s + t - u)))
     -        + (k1k + k2k)*(-2*(glmd*glp - glm*glpd)*(gqm2 - gqp2)*s - 
     -         (glmd*glp + glm*glpd)*(gqm2 + gqp2)*t + 
     -         (glmd*glp + glm*glpd)*(gqm2 + gqp2)*u) - 
     -      (2*glp2*gqm2*k1k2*s + glm*glpd*gqm2*k1k2*s + 
     -         2*glm2*gqp2*k1k2*s + glm*glpd*gqp2*k1k2*s + 
     -         2*glm2*gqm2*k1k*k2k*s + glm*glpd*gqm2*k1k*k2k*s + 
     -         2*glp2*gqp2*k1k*k2k*s + glm*glpd*gqp2*k1k*k2k*s + 
     -         glp2*gqm2*k1k*s2 + glm*glpd*gqm2*k1k*s2 + 
     -         glm2*gqp2*k1k*s2 + glm*glpd*gqp2*k1k*s2 + 
     -         2*glm2*gqm2*k1k2*t + 2*glp2*gqm2*k1k2*t + 
     -         2*glm2*gqp2*k1k2*t + 2*glp2*gqp2*k1k2*t + 
     -         2*glm2*gqm2*k1k*k2k*t + 2*glp2*gqm2*k1k*k2k*t + 
     -         2*glm2*gqp2*k1k*k2k*t + 2*glp2*gqp2*k1k*k2k*t + 
     -         2*glp2*gqm2*k1k*s*t + 2*glm2*gqp2*k1k*s*t + 
     -         glp2*gqm2*k2k*s*t + glm2*gqp2*k2k*s*t + 
     -         (glm*glpd*(gqm2 + gqp2)*(k1k + k2k)*s + 
     -            glp2*(2*gqp2*(k1k2 + k1k*(k2k + t) + k2k*(s + t)) + 
     -           gqm2*(2*k1k2 + 2*k1k*(k2k + s + t) + k2k*(s + 2*t)))
     -          + glm2*(2*gqm2*(k1k2 + k1k*(k2k + t) + k2k*(s + t)) + 
     -          gqp2*(2*k1k2 + 2*k1k*(k2k + s + t) + k2k*(s + 2*t))))*
     -          u + (glm2 + glp2)*(gqm2 + gqp2)*(k1k + k2k)*u2 + 
     -    glmd*glp*(gqm2 + gqp2)*s*(k1k*(k1k + k2k + s) + (k1k + k2k)*u)
     -         )/p2k - (p2k*((glmd*glp + glm*glpd)*(gqm2 + gqp2)*
     -            (2*k1k*s + (k1k + k2k)*t) + 
     -           glm2*(2*gqp2*(k1k2 + 2*k1k*k2k + 2*k2k*t) + 
     -              gqm2*(2*k1k*(s + t + u) + k2k*(s + 3*t + u))) + 
     -           glp2*(2*gqm2*(k1k2 + 2*k1k*k2k + 2*k2k*t) + 
     -            gqp2*(2*k1k*(s + t + u) + k2k*(s + 3*t + u)))))/p1k + 
     -      (p1k*((glmd*glp + glm*glpd)*(gqm2 + gqp2)*
     -            (2*k1k*s + (k1k + k2k)*u) + 
     -           glp2*(2*gqp2*(k1k2 + 2*k1k*k2k + 2*k2k*u) + 
     -              gqm2*(2*k1k*(s + t + u) + k2k*(s + t + 3*u))) + 
     -           glm2*(2*gqm2*(k1k2 + 2*k1k*k2k + 2*k2k*u) + 
     -            gqp2*(2*k1k*(s + t + u) + k2k*(s + t + 3*u)))))/p2k + 
     -      (glm*glpd*gqm2*k1k2*s + 2*glp2*gqp2*k1k2*s + 
     -         glm*glpd*gqp2*k1k2*s + 2*glp2*gqm2*k1k*k2k*s + 
     -         glm*glpd*gqm2*k1k*k2k*s + glm*glpd*gqp2*k1k*k2k*s + 
     -         glm*glpd*gqm2*k1k*s2 + glp2*gqp2*k1k*s2 + 
     -         glm*glpd*gqp2*k1k*s2 + 2*glp2*gqm2*k1k2*t + 
     -         2*glp2*gqp2*k1k2*t + 2*glp2*gqm2*k1k*k2k*t + 
     -         2*glp2*gqp2*k1k*k2k*t + glm*glpd*gqm2*k1k*s*t + 
     -         2*glp2*gqp2*k1k*s*t + glm*glpd*gqp2*k1k*s*t + 
     -         2*glp2*gqm2*k2k*s*t + glm*glpd*gqm2*k2k*s*t + 
     -         glp2*gqp2*k2k*s*t + glm*glpd*gqp2*k2k*s*t + 
     -         glp2*gqm2*k1k*t2 + glp2*gqp2*k1k*t2 + 
     -         glp2*gqm2*k2k*t2 + glp2*gqp2*k2k*t2 + 
     -         glmd*glp*(gqm2 + gqp2)*s*
     -          (k1k*(k1k + k2k + s) + (k1k + k2k)*t) + 
     -         glp2*(2*gqm2*(k1k + k2k)*(k1k + t) + 
     -        gqp2*(2*k1k2 + 2*k1k*(k2k + s + t) + k2k*(s + 2*t)))*u + 
     -         glm2*(gqm2*(2*k1k + s + t)*(k2k*t + k1k*(s + t)) + 
     -         gqm2*(2*k1k2 + 2*k1k*(k2k + s + t) + k2k*(s + 2*t))*u + 
     -            gqp2*(2*k1k2*(t + u) + k2k*t*(2*s + t + 2*u) + 
     -         k1k*(2*k2k*(s + t + u) + t*(t + 2*u)))))/p1k))/(k1k*k2k)+
     -
     -
     -         (-4*denspd*denspg*ml2*qq**3*
     -    (glmd*(gqpd*(4*p1k2 + 4*p2k2 - 2*k1k*(p1k + p2k - s) - 
     -            3*p1k*s - 5*p2k*s + 2*s2 - p1k*t + p2k*t + s*t + 
     -            (p1k - p2k + s)*u) + 
     -       gqmd*(4*p1k2 + 4*p2k2 - 2*k1k*(p1k + p2k - s) - 5*p1k*s - 
     -       3*p2k*s + 2*s2 - p1k*t + p2k*t + s*t + (p1k - p2k + s)*u)
     -     ) + glpd*(gqmd*(4*p1k2 + 4*p2k2 - 2*k1k*(p1k + p2k - s) - 
     -            3*p1k*s - 5*p2k*s + 2*s2 - p1k*t + p2k*t + s*t + 
     -            (p1k - p2k + s)*u) + 
     -     gqpd*(4*p1k2 + 4*p2k2 - 2*k1k*(p1k + p2k - s) - 5*p1k*s - 
     -       3*p2k*s + 2*s2 - p1k*t + p2k*t + s*t + (p1k - p2k + s)*u)
     -         )))/(p1k*p2k)+
     -         (-2*densg*densp*ml2*qq2*
     -    (glp*(gqp*(k1k2*(2*p1k2 + p1k*(2*p2k - s - 2*(t + u)) + 
     -               p2k*(3*s + 2*(t + u))) + 
     -            k2k*(p1k2*(-6*p2k + 5*u) + 
     -               p1k*(6*p2k2 + p2k*s - u*(3*s + 2*t + u)) + 
     -           p2k*(2*s*t + t2 + s*u + 2*t*u - p2k*(s + 4*t + u))) + 
     -            k1k*(p1k2*(2*(p2k + s) + u) + 
     -           p2k*((s + t)*(2*s + t + 2*u) - p2k*(4*s + 3*t + 2*u)) - 
     -               p1k*(2*p2k2 + s2 + s*u + 2*t*u + u2 + 
     -                  p2k*(4*s + t + 3*u)) + 
     -               k2k*(4*p1k2 + p1k*(4*p2k - 3*s - 2*(t + u)) + 
     -                  p2k*(s + 2*(t + u))))) + 
     -         gqm*(k2k*(p2k*t*(-5*p2k + 3*s + t + 2*u) + 
     -               p1k2*(-6*p2k + s + t + 4*u) + 
     -             p1k*(6*p2k2 - p2k*s - s*t - 2*(s + t)*u - u2)) + 
     -            k1k2*(p2k*(-2*p2k + s + 2*(t + u)) - 
     -               p1k*(2*p2k + 3*s + 2*(t + u))) + 
     -            k1k*(p1k2*(2*(p2k + 2*s + t) + 3*u) + 
     -               p2k*(s2 + s*t + t2 - p2k*(2*s + t) + 2*t*u) - 
     -               p1k*(2*p2k2 - p2k*(4*s + 3*t + u) + 
     -                  (s + u)*(2*(s + t) + u)) - 
     -               k2k*(p2k*(4*p2k - 3*s - 2*(t + u)) + 
     -                  p1k*(4*p2k + s + 2*(t + u)))))) + 
     -      glm*(gqm*(k1k2*(2*p1k2 + p1k*(2*p2k - s - 2*(t + u)) + 
     -               p2k*(3*s + 2*(t + u))) + 
     -            k2k*(p1k2*(-6*p2k + 5*u) + 
     -               p1k*(6*p2k2 + p2k*s - u*(3*s + 2*t + u)) + 
     -           p2k*(2*s*t + t2 + s*u + 2*t*u - p2k*(s + 4*t + u))) + 
     -            k1k*(p1k2*(2*(p2k + s) + u) + 
     -           p2k*((s + t)*(2*s + t + 2*u) - p2k*(4*s + 3*t + 2*u)) - 
     -               p1k*(2*p2k2 + s2 + s*u + 2*t*u + u2 + 
     -                  p2k*(4*s + t + 3*u)) + 
     -               k2k*(4*p1k2 + p1k*(4*p2k - 3*s - 2*(t + u)) + 
     -                  p2k*(s + 2*(t + u))))) + 
     -         gqp*(k2k*(p2k*t*(-5*p2k + 3*s + t + 2*u) + 
     -               p1k2*(-6*p2k + s + t + 4*u) + 
     -             p1k*(6*p2k2 - p2k*s - s*t - 2*(s + t)*u - u2)) + 
     -            k1k2*(p2k*(-2*p2k + s + 2*(t + u)) - 
     -               p1k*(2*p2k + 3*s + 2*(t + u))) + 
     -            k1k*(p1k2*(2*(p2k + 2*s + t) + 3*u) + 
     -               p2k*(s2 + s*t + t2 - p2k*(2*s + t) + 2*t*u) - 
     -               p1k*(2*p2k2 - p2k*(4*s + 3*t + u) + 
     -                  (s + u)*(2*(s + t) + u)) - 
     -               k2k*(p2k*(4*p2k - 3*s - 2*(t + u)) + 
     -              p1k*(4*p2k + s + 2*(t + u))))))))/(k1k*k2k*p1k*p2k)+
     -
     -
     -
     -         (2*densd*densp*ml2*qq*(4*glm*glpd*gqm2*k1k*p1k2*p2k + 
     -      4*glmd*glp*gqp2*k1k*p1k2*p2k - 
     -      4*glmd*glp*gqm2*k2k*p1k2*p2k - 
     -      4*glm*glpd*gqp2*k2k*p1k2*p2k - 
     -      4*glmd*glp*gqm2*k1k*p1k*p2k2 - 
     -      4*glm*glpd*gqp2*k1k*p1k*p2k2 + 
     -      4*glm*glpd*gqm2*k2k*p1k*p2k2 + 
     -     4*glmd*glp*gqp2*k2k*p1k*p2k2 - glmd*glp*gqm2*k1k2*p1k*s - 
     -      glm*glpd*gqm2*k1k2*p1k*s - glmd*glp*gqp2*k1k2*p1k*s - 
     -      glm*glpd*gqp2*k1k2*p1k*s - glmd*glp*gqm2*k1k*k2k*p1k*s - 
     -      glm*glpd*gqm2*k1k*k2k*p1k*s - glmd*glp*gqp2*k1k*k2k*p1k*s - 
     -      glm*glpd*gqp2*k1k*k2k*p1k*s + 2*glmd*glp*gqm2*k1k*p1k2*s + 
     -     2*glm*glpd*gqm2*k1k*p1k2*s + 2*glmd*glp*gqp2*k1k*p1k2*s + 
     -      2*glm*glpd*gqp2*k1k*p1k2*s + glmd*glp*gqm2*k1k2*p2k*s + 
     -      glm*glpd*gqm2*k1k2*p2k*s + glmd*glp*gqp2*k1k2*p2k*s + 
     -      glm*glpd*gqp2*k1k2*p2k*s + glmd*glp*gqm2*k1k*k2k*p2k*s + 
     -      glm*glpd*gqm2*k1k*k2k*p2k*s + glmd*glp*gqp2*k1k*k2k*p2k*s + 
     -    glm*glpd*gqp2*k1k*k2k*p2k*s + 2*glmd*glp*gqm2*k1k*p1k*p2k*s - 
     -   2*glm*glpd*gqm2*k1k*p1k*p2k*s - 2*glmd*glp*gqp2*k1k*p1k*p2k*s + 
     -   2*glm*glpd*gqp2*k1k*p1k*p2k*s + 2*glmd*glp*gqm2*k2k*p1k*p2k*s - 
     -   2*glm*glpd*gqm2*k2k*p1k*p2k*s - 2*glmd*glp*gqp2*k2k*p1k*p2k*s + 
     -    2*glm*glpd*gqp2*k2k*p1k*p2k*s - 2*glmd*glp*gqm2*k1k*p2k2*s - 
     -     2*glm*glpd*gqm2*k1k*p2k2*s - 2*glmd*glp*gqp2*k1k*p2k2*s - 
     -      2*glm*glpd*gqp2*k1k*p2k2*s - glmd*glp*gqm2*k1k*p1k*s2 - 
     -      glm*glpd*gqm2*k1k*p1k*s2 - glmd*glp*gqp2*k1k*p1k*s2 - 
     -      glm*glpd*gqp2*k1k*p1k*s2 + glmd*glp*gqm2*k1k*p2k*s2 + 
     -      glm*glpd*gqm2*k1k*p2k*s2 + glmd*glp*gqp2*k1k*p2k*s2 + 
     -      glm*glpd*gqp2*k1k*p2k*s2 - glmd*glp*gqm2*k1k*p1k*p2k*t - 
     -      glm*glpd*gqm2*k1k*p1k*p2k*t - glmd*glp*gqp2*k1k*p1k*p2k*t - 
     -      glm*glpd*gqp2*k1k*p1k*p2k*t - glmd*glp*gqm2*k2k*p1k*p2k*t - 
     -      glm*glpd*gqm2*k2k*p1k*p2k*t - glmd*glp*gqp2*k2k*p1k*p2k*t - 
     -      glm*glpd*gqp2*k2k*p1k*p2k*t - glmd*glp*gqm2*k1k*p2k2*t - 
     -      glm*glpd*gqm2*k1k*p2k2*t - glmd*glp*gqp2*k1k*p2k2*t - 
     -      glm*glpd*gqp2*k1k*p2k2*t - glmd*glp*gqm2*k2k*p2k2*t - 
     -      glm*glpd*gqm2*k2k*p2k2*t - glmd*glp*gqp2*k2k*p2k2*t - 
     -      glm*glpd*gqp2*k2k*p2k2*t + glmd*glp*gqm2*k1k*p2k*s*t + 
     -      glm*glpd*gqm2*k1k*p2k*s*t + glmd*glp*gqp2*k1k*p2k*s*t + 
     -      glm*glpd*gqp2*k1k*p2k*s*t + glmd*glp*gqm2*k2k*p2k*s*t + 
     -      glm*glpd*gqm2*k2k*p2k*s*t + glmd*glp*gqp2*k2k*p2k*s*t + 
     -      glm*glpd*gqp2*k2k*p2k*s*t + 
     -      (glmd*glp + glm*glpd)*(gqm2 + gqp2)*(k1k + k2k)*p1k*
     -       (p1k + p2k - s)*u + 
     -      glm2*(gqp2*(-2*k1k2*
     -             (p2k*(p2k - t - u) + p1k*(p2k + s + t + u)) + 
     -            k2k*(p2k*t*(-4*p2k + 2*s + t + 2*u) + 
     -               p1k2*(-2*p2k + s + t + 3*u) + 
     -               p1k*(6*p2k2 + p2k*(-3*s + t - u) - s*(t + u) - 
     -                  u*(2*t + u))) + 
     -            k1k*(2*p1k2*(p2k + s + t + u) + p2k*t*(t + 2*u) + 
     -               2*k2k*(-(p1k*(2*p2k + t + u)) + 
     -                  p2k*(-2*p2k + s + t + u)) + 
     -          p1k*(2*p2k2 + 2*p2k*(s + 2*t) - (s + u)*(s + 2*t + u))
     -               )) + gqm2*
     -        (2*k1k2*(p1k2 + p1k*(p2k - t - u) + p2k*(s + t + u)) + 
     -            k1k*(-2*p1k2*p2k - p2k*(2*p2k - s - t)*(s + t) + 
     -               2*p2k*(-p2k + s + t)*u - 
     -               p1k*(2*p2k*(p2k + s) + 2*(2*p2k + t)*u + u2) + 
     -          2*k2k*(2*p1k2 + p1k*(2*p2k - s - t - u) + p2k*(t + u))
     -               ) + k2k*(p1k2*(-6*p2k + 4*u) + 
     -               p2k*(t*(s + t) + (s + 2*t)*u - p2k*(s + 3*t + u)) + 
     -          p1k*(2*p2k2 + p2k*(3*s + t - u) - u*(2*(s + t) + u))))
     -         ) + glp2*(gqm2*(-2*k1k2*
     -             (p2k*(p2k - t - u) + p1k*(p2k + s + t + u)) + 
     -            k2k*(p2k*t*(-4*p2k + 2*s + t + 2*u) + 
     -               p1k2*(-2*p2k + s + t + 3*u) + 
     -               p1k*(6*p2k2 + p2k*(-3*s + t - u) - s*(t + u) - 
     -                  u*(2*t + u))) + 
     -            k1k*(2*p1k2*(p2k + s + t + u) + p2k*t*(t + 2*u) + 
     -               2*k2k*(-(p1k*(2*p2k + t + u)) + 
     -                  p2k*(-2*p2k + s + t + u)) + 
     -          p1k*(2*p2k2 + 2*p2k*(s + 2*t) - (s + u)*(s + 2*t + u))
     -               )) + gqp2*
     -        (2*k1k2*(p1k2 + p1k*(p2k - t - u) + p2k*(s + t + u)) + 
     -            k1k*(-2*p1k2*p2k - p2k*(2*p2k - s - t)*(s + t) + 
     -               2*p2k*(-p2k + s + t)*u - 
     -               p1k*(2*p2k*(p2k + s) + 2*(2*p2k + t)*u + u2) + 
     -          2*k2k*(2*p1k2 + p1k*(2*p2k - s - t - u) + p2k*(t + u))
     -               ) + k2k*(p1k2*(-6*p2k + 4*u) + 
     -               p2k*(t*(s + t) + (s + 2*t)*u - p2k*(s + 3*t + u)) + 
     -          p1k*(2*p2k2 + p2k*(3*s + t - u) - u*(2*(s + t) + u))))
     -         )))/(k1k*k2k*p1k*p2k)+
     -
     -
     -
     -   (-4*densp*denspg*ml2*qq**3*
     -  (glm*(gqp*(4*p1k2 + 4*p2k2 - 2*k1k*(p1k + p2k - s) - 3*p1k*s - 
     -       5*p2k*s + 2*s2 - p1k*t + p2k*t + s*t + (p1k - p2k + s)*u)
     -           + gqm*(4*p1k2 + 4*p2k2 - 2*k1k*(p1k + p2k - s) - 
     -            5*p1k*s - 3*p2k*s + 2*s2 - p1k*t + p2k*t + s*t + 
     -            (p1k - p2k + s)*u)) + 
     -      glp*(gqm*(4*p1k2 + 4*p2k2 - 2*k1k*(p1k + p2k - s) - 
     -            3*p1k*s - 5*p2k*s + 2*s2 - p1k*t + p2k*t + s*t + 
     -            (p1k - p2k + s)*u) + 
     -     gqp*(4*p1k2 + 4*p2k2 - 2*k1k*(p1k + p2k - s) - 5*p1k*s - 
     -       3*p2k*s + 2*s2 - p1k*t + p2k*t + s*t + (p1k - p2k + s)*u)
     -         )))/(p1k*p2k)+
     -       (2*ml2*ml4*(2*dens2*(glm2 + glp2)*(gqm2 + gqp2)*
     -       (2*k1k2 + 2*k1k*k2k + k2k2)*p1k*p2k - 
     -      qq*(-(densg*qq*(densp*(glm + glp)*(gqm + gqp)*k1k*k2k*
     -               (k1k + k2k)*(p1k - p2k) + 
     -              8*densg*(2*k1k2 + 2*k1k*k2k + k2k2)*p1k*p2k + 
     -              denspd*(glmd + glpd)*(gqmd + gqpd)*k1k*k2k*
     -            ((k1k + k2k)*p1k - (k1k + k2k - 2*k1k*mz2dm1*p1k)*p2k)
     -               - 8*denspg*k1k*k2k*(k1k + k2k)*(p1k - p2k)*qq)) + 
     -     dens*(denspd*(glm2 + glp2)*(gqm2 + gqp2)*k1k*k2k*(k1k + k2k)*
     -             (p1k - p2k) + 
     -            (glm + glp)*(gqm + gqp)*
     -             (2*densg*(2*k1k2 + 2*k1k*k2k + k2k2)*p1k*p2k - 
     -               denspg*k1k*k2k*(k1k + k2k)*(p1k - p2k)*qq)) + 
     -     densd*(densp*(glm2 + glp2)*(gqm2 + gqp2)*k1k*k2k*(k1k + k2k)*
     -             (p1k - p2k) + 
     -            (glmd + glpd)*(gqmd + gqpd)*
     -             (2*densg*(2*k1k2 + 2*k1k*k2k + k2k2)*p1k*p2k - 
     -               denspg*k1k*k2k*(k1k + k2k)*(p1k - p2k)*qq)))))/
     -  (k1k2*k2k2*p1k*p2k) +
     -        (-4*qq*(dens*(denspd*((glmd*glp + glm*glpd)*(gqm2 - gqp2)*
     -             (k1k + k2k)*ml2*(p1k - p2k) + 
     -            glp2*(gqp2*(2*k1k2*p1k - 2*k1k*p1k*(ml2 + p2k - u) + 
     -                  k2k*p2k*(-2*p1k + s + t - u)) + 
     -               gqm2*(2*k1k2*p2k - 2*k1k*p2k*(ml2 + p1k - t) + 
     -                  k2k*p1k*(-2*p2k + s - t + u))) - 
     -            glm2*(gqm2*(2*k1k2*p1k - 2*k1k*p1k*(ml2 + p2k - u) + 
     -                  k2k*p2k*(-2*p1k + s + t - u)) + 
     -               gqp2*(2*k1k2*p2k - 2*k1k*p2k*(ml2 + p1k - t) + 
     -                  k2k*p1k*(-2*p2k + s - t + u)))) + 
     -         denspg*qq*(-(glp*
     -               (gqp*(2*k1k2*p1k + 
     -                    k2k*
     -                   (ml2*(-p1k + p2k) + p2k*(-2*p1k + s + t - u)) + 
     -                    k1k*(ml2*(-3*p1k + p2k) + 2*p1k*(-p2k + u))) + 
     -                 gqm*(2*k1k2*p2k + 
     -                    k1k*(ml2*(p1k - 3*p2k) + 2*p2k*(-p1k + t)) + 
     -               k2k*(ml2*(p1k - p2k) + p1k*(-2*p2k + s - t + u)))))
     -              + glm*(gqm*
     -                (2*k1k2*p1k + 
     -              k2k*(ml2*(-p1k + p2k) + p2k*(-2*p1k + s + t - u)) + 
     -                  k1k*(ml2*(-3*p1k + p2k) + 2*p1k*(-p2k + u))) + 
     -               gqp*(2*k1k2*p2k + 
     -                  k1k*(ml2*(p1k - 3*p2k) + 2*p2k*(-p1k + t)) + 
     -              k2k*(ml2*(p1k - p2k) + p1k*(-2*p2k + s - t + u))))))
     -       + densd*(densp*(-((glmd*glp + glm*glpd)*(gqm2 - gqp2)*
     -               (k1k + k2k)*ml2*(p1k - p2k)) - 
     -            glp2*(gqp2*(2*k1k2*p1k - 2*k1k*p1k*(ml2 + p2k - u) + 
     -                  k2k*p2k*(-2*p1k + s + t - u)) + 
     -               gqm2*(2*k1k2*p2k - 2*k1k*p2k*(ml2 + p1k - t) + 
     -                  k2k*p1k*(-2*p2k + s - t + u))) + 
     -            glm2*(gqm2*(2*k1k2*p1k - 2*k1k*p1k*(ml2 + p2k - u) + 
     -                  k2k*p2k*(-2*p1k + s + t - u)) + 
     -               gqp2*(2*k1k2*p2k - 2*k1k*p2k*(ml2 + p1k - t) + 
     -                  k2k*p1k*(-2*p2k + s - t + u)))) + 
     -         denspg*qq*(glpd*
     -             (gqpd*(2*k1k2*p1k + 
     -              k2k*(ml2*(-p1k + p2k) + p2k*(-2*p1k + s + t - u)) + 
     -                  k1k*(ml2*(-3*p1k + p2k) + 2*p1k*(-p2k + u))) + 
     -               gqmd*(2*k1k2*p2k + 
     -                  k1k*(ml2*(p1k - 3*p2k) + 2*p2k*(-p1k + t)) + 
     -              k2k*(ml2*(p1k - p2k) + p1k*(-2*p2k + s - t + u)))) - 
     -            glmd*(gqmd*(2*k1k2*p1k + 
     -               k2k*(ml2*(-p1k + p2k) + p2k*(-2*p1k + s + t - u)) + 
     -                  k1k*(ml2*(-3*p1k + p2k) + 2*p1k*(-p2k + u))) + 
     -               gqpd*(2*k1k2*p2k + 
     -                  k1k*(ml2*(p1k - 3*p2k) + 2*p2k*(-p1k + t)) + 
     -              k2k*(ml2*(p1k - p2k) + p1k*(-2*p2k + s - t + u))))))
     -       + densg*qq*(densp*
     -          (glp*(gqp*(2*k1k2*p1k + 
     -               k2k*(ml2*(-p1k + p2k) + p2k*(-2*p1k + s + t - u)) + 
     -                  k1k*(ml2*(-3*p1k + p2k) + 2*p1k*(-p2k + u))) + 
     -               gqm*(2*k1k2*p2k + 
     -                  k1k*(ml2*(p1k - 3*p2k) + 2*p2k*(-p1k + t)) + 
     -              k2k*(ml2*(p1k - p2k) + p1k*(-2*p2k + s - t + u)))) - 
     -            glm*(gqm*(2*k1k2*p1k + 
     -               k2k*(ml2*(-p1k + p2k) + p2k*(-2*p1k + s + t - u)) + 
     -                  k1k*(ml2*(-3*p1k + p2k) + 2*p1k*(-p2k + u))) + 
     -               gqp*(2*k1k2*p2k + 
     -                  k1k*(ml2*(p1k - 3*p2k) + 2*p2k*(-p1k + t)) + 
     -              k2k*(ml2*(p1k - p2k) + p1k*(-2*p2k + s - t + u)))))
     -          + denspd*(-(glpd*
     -               (gqpd*(2*k1k2*p1k + 
     -                    k2k*
     -                  (ml2*(-p1k + p2k) + p2k*(-2*p1k + s + t - u)) + 
     -                    k1k*(ml2*(-3*p1k + p2k) + 2*p1k*(-p2k + u))) + 
     -                 gqmd*(2*k1k2*p2k + 
     -                    k1k*(ml2*(p1k - 3*p2k) + 2*p2k*(-p1k + t)) + 
     -               k2k*(ml2*(p1k - p2k) + p1k*(-2*p2k + s - t + u)))))
     -              + glmd*(gqmd*
     -                (2*k1k2*p1k + 
     -              k2k*(ml2*(-p1k + p2k) + p2k*(-2*p1k + s + t - u)) + 
     -                  k1k*(ml2*(-3*p1k + p2k) + 2*p1k*(-p2k + u))) + 
     -               gqpd*(2*k1k2*p2k + 
     -                  k1k*(ml2*(p1k - 3*p2k) + 2*p2k*(-p1k + t)) + 
     -             k2k*(ml2*(p1k - p2k) + p1k*(-2*p2k + s - t + u)))))))
     -     *e_(p1,p2,k,k1))/(k1k*k2k*p1k*p2k)

      amp2 = dble( amp2 * el2_scheme*dconjg(el2_scheme) * el2 
     +             /4/3/(st_alpha/(2*pi)) )

      end
*
**
*
      complex*16 function e_(q1,q2,q3,q4)
      implicit none
      real*8 q1(0:3),q2(0:3),q3(0:3),q4(0:3)

      complex*16 esign

      esign = -(0.d0,1.d0)

      e_ =q1(1)*q2(2)*q3(3)*q4(0)-q1(1)*q2(2)*q3(0)*q4(3)+
     .    q1(1)*q2(3)*q3(0)*q4(2)-q1(1)*q2(3)*q3(2)*q4(0)+
     .    q1(1)*q2(0)*q3(2)*q4(3)-q1(1)*q2(0)*q3(3)*q4(2)+
     .    q1(2)*q2(1)*q3(0)*q4(3)-q1(2)*q2(1)*q3(3)*q4(0)+
     .    q1(2)*q2(3)*q3(1)*q4(0)-q1(2)*q2(3)*q3(0)*q4(1)+
     .    q1(2)*q2(0)*q3(3)*q4(1)-q1(2)*q2(0)*q3(1)*q4(3)+
     .    q1(3)*q2(1)*q3(2)*q4(0)-q1(3)*q2(1)*q3(0)*q4(2)+
     .    q1(3)*q2(2)*q3(0)*q4(1)-q1(3)*q2(2)*q3(1)*q4(0)+
     .    q1(3)*q2(0)*q3(1)*q4(2)-q1(3)*q2(0)*q3(2)*q4(1)+
     .    q1(0)*q2(1)*q3(3)*q4(2)-q1(0)*q2(1)*q3(2)*q4(3)+
     .    q1(0)*q2(2)*q3(1)*q4(3)-q1(0)*q2(2)*q3(3)*q4(1)+
     .    q1(0)*q2(3)*q3(2)*q4(1)-q1(0)*q2(3)*q3(1)*q4(2)

      e_ = esign * e_ 

      return
      end function e_
