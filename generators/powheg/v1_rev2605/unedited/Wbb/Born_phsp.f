      subroutine born_phsp(xborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      real * 8 xborn(ndiminteg-3)
      real * 8 m2,xjac,taumin,tau,y,betaCM,vec(3)
      integer mu,j
      logical ini
      data ini/.true./
      save ini
      logical check
      parameter(check=.false.)
      real * 8 epsilon
      parameter (epsilon=1d-10)
      real*8 rsh
      real*8 alam,alam1,sh,s23,rs23
      real*8 ct,st,ct23,st23,ph23,ph
      real*8 p10,p1b,gp230,gp23b,p230,p23b,p20,p2b,p30,p3b
      real * 8 z
      external alam
      real * 8 mass(3)
      real * 8 compmass
      external compmass
      real * 8 gam,massa,r
      real * 8 random
      external random
      real * 8 xmax,xmin,norm, rmin,rmax,atanxmax, atanxmin


      
      if(ini) then
c     set initial- and final-state masses for Born and real         
         kn_masses(1)=0
         kn_masses(2)=0
         kn_masses(4)=ph_mb
         kn_masses(5)=ph_mb
         kn_masses(nlegreal)=0
         ini=.false.        
      endif
     
      xmax = (ph_Wmass2high-ph_Wmass2)/ph_WmWw
      xmin = (ph_Wmass2low-ph_Wmass2)/ph_WmWw
      atanxmax = atan(xmax)
      atanxmin = atan(xmin)
      norm = 1/(atanxmax-atanxmin)
      rmin = norm* atanxmin
      rmax = norm* atanxmax

c      r = rmin + (rmax-rmin)*xborn(7)
      r = rmin + (rmax-rmin)*random()

      m2 = ph_Wmass2 +  ph_WmWw *tan(r/norm)
           
c 11   m2 = virt2(ph_Wmass,ph_Wwidth,xborn(7))
c      if (m2.lt.ph_Wmass2low.or.m2.gt.ph_Wmass2high) then
c         goto 11
c         write(*,*) m2
c      endif

c     set the W mass to its virtuality!!
      kn_masses(3)=sqrt(m2)

      xjac = 1d0
      kn_minmass=2*ph_mb+kn_masses(3)
      m2=kn_minmass**2

c     d x1 d x2 = d tau d y;
      taumin=m2/kn_sbeams
      tau=exp(log(taumin)*(1-xborn(1)**2))
      xjac=xjac*tau*abs(log(taumin))*2*xborn(1)
      kn_sborn=kn_sbeams*tau
c     ymax=|log(tau)|/2
      y=-(1-2*xborn(2))*log(tau)/2
      xjac=-xjac*log(tau)
c     Build kinematics
      kn_xb1=sqrt(tau)*exp(y)
      kn_xb2=tau/kn_xb1

      rsh=sqrt(kn_sborn)
      sh=kn_sborn
      mass(1)=kn_masses(4)
      mass(2)=kn_masses(5)
      mass(3)=kn_masses(3)
      
      kn_cmpborn(1,1)=0.d0
      kn_cmpborn(2,1)=0.d0
      kn_cmpborn(3,1)=rsh/2d0
      kn_cmpborn(0,1)=rsh/2d0
      kn_cmpborn(1,2)=0.d0
      kn_cmpborn(2,2)=0.d0
      kn_cmpborn(3,2)=-rsh/2d0
      kn_cmpborn(0,2)=rsh/2d0     


c generation of invariants and angles:
      s23=(rsh-mass(1))**2*xborn(3)+(1d0-xborn(3))*(mass(2)+mass(3))**2
      rs23=dsqrt(s23)
      xjac=xjac*((rsh-mass(1))**2-(mass(2)+mass(3))**2)
      ct23=2d0*xborn(4)**2-1d0
      xjac=xjac*4*xborn(4)
c      ct23=2d0*xborn(4)-1d0
      st23=dsqrt(dabs((1d0-ct23)*(1d0+ct23)))
c
c      ct=2d0*xborn(5)-1d0
c      xjac=xjac*2d0
      z=1-2*xborn(5)
      xjac=xjac*2
      ct=1.5d0*(z-z**3/3)
      xjac=xjac*1.5d0*(1-z**2)

      st=dsqrt(dabs((1d0-ct)*(1d0+ct)))
      ph23=2d0*pi*xborn(6)
c      ph=2d0*pi*x(5)
      ph=0d0
      xjac=xjac*(2d0*pi)**2

c decay a+b->p1+p23 in the CMS vec(a+b)=vec(0) 
c with vec(a) is pos. z-axis:
      p10=(sh+mass(1)**2-s23)/2d0/rsh
      p230=(sh+s23-mass(1)**2)/2d0/rsh
      alam1=alam(rsh,mass(1),rs23)
      p1b=alam1/2d0/rsh
      p23b=p1b
      gp230=p230/rs23
      gp23b=p23b/rs23
      kn_cmpborn(1,4)=p1b*st*dcos(ph)
      kn_cmpborn(2,4)=p1b*st*dsin(ph)
      kn_cmpborn(3,4)=p1b*ct
      kn_cmpborn(0,4)=p10

c no flux factor (1/(2*kn_sborn))
      xjac=xjac*alam1/8d0/sh

c decay p23->p2+p3 in the CMS vec(p23)=vec(0) 
c with vec(p23) is pos. z-axis:
      p20=(s23+mass(2)**2-mass(3)**2)/2d0/rs23
      p30=(s23+mass(3)**2-mass(2)**2)/2d0/rs23
      alam1=alam(rs23,mass(2),mass(3))
      p2b=alam1/2d0/rs23
      p3b=p2b
      kn_cmpborn(1,5)=p2b*st23*dcos(ph23)
      kn_cmpborn(2,5)=p2b*st23*dsin(ph23)
      kn_cmpborn(3,5)=p2b*ct23
      kn_cmpborn(0,5)=p20
      kn_cmpborn(1,3)=-p3b*st23*dcos(ph23)
      kn_cmpborn(2,3)=-p3b*st23*dsin(ph23)
      kn_cmpborn(3,3)=-p3b*ct23
      kn_cmpborn(0,3)=p30
      xjac=xjac*alam1/8d0/s23
c
c normalization of the phase space weight:
      xjac=xjac/(2d0*pi)**5

      kn_jacborn=xjac

c boost and rotation, so that pa is positive z axis:
      call boost_z(gp230,gp23b,kn_cmpborn(0,3),kn_cmpborn(3,3))
      call rotation_2(-1d0,st,ct,dsin(ph),dcos(ph),
     $     kn_cmpborn(1,3),kn_cmpborn(2,3),kn_cmpborn(3,3))

      call boost_z(gp230,gp23b,kn_cmpborn(0,5),kn_cmpborn(3,5))
      call rotation_2(-1d0,st,ct,dsin(ph),dcos(ph),
     $     kn_cmpborn(1,5),kn_cmpborn(2,5),kn_cmpborn(3,5))


c     boost back in the lab frame
c     now boost everything along 3rd axis
      betaCM=(kn_xb1-kn_xb2)/(kn_xb1+kn_xb2)
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(nlegborn,vec,betaCM,kn_cmpborn(0,1),kn_pborn(0,1))      
      
      if (check) then
         write(*,*)
         do j=1,nlegborn
            write(*,*) 'mom CM',j,(kn_cmpborn(mu,j),mu=0,3)
         enddo
         do j=1,nlegborn
            write(*,*) 'mom ',j,(kn_pborn(mu,j),mu=0,3)
            write(*,*) 'mass ',j,compmass(kn_pborn(0,j))
         enddo
         
         call checkmomzero(nlegborn,kn_pborn)
      endif
      end





********************************************************
      subroutine boost_z(ecm,pcm,p4,p3)
********************************************************
c this subroutine performs a boost along the z direction
********************************************************
      implicit none
      real*8 ecm,pcm,p4,p3
      real*8 gamma,bgamma,E,Pp
      gamma=ecm
      bgamma=pcm
      E=p4
      Pp=p3
      p4=gamma*E+bgamma*Pp
      p3=gamma*Pp+bgamma*E 
      return
      end
********************************************************
      subroutine rotation_2(vz,st,ct,sp,cp,p1,p2,p3)
********************************************************
c this subroutine performs two rotations
********************************************************
      implicit none
      real*8 st,ct,sp,cp,p1,p2,p3
      real*8 p1s,p2s,p3s,vz
      p1s=p1
      p2s=p2
      p3s=p3
      p1=vz*(cp*(ct*p1s+st*p3s)-sp*p2s)
      p2=vz*(sp*(ct*p1s+st*p3s)+cp*p2s)
      p3=vz*(-st*p1s+ct*p3s)
      return
      end
********************************************************
      real*8 function alam(m,m1,m2)
      implicit none
      real*8 m,m1,m2
      real*8 s,s1,s2
      real*8 aux
      s=m**2
      s1=m1**2
      s2=m2**2
      aux=s**2+s1**2+s2**2-2*s*s1-2*s*s2-2*s1*s2
      if(aux.lt.0.d0)aux=0.d0	
      alam=dsqrt(aux)
      return
      end	
********************************************************





      function compmass(p)
      implicit none
      real * 8 p(0:3),compmass
      compmass = sqrt(abs(p(0)**2-p(1)**2-p(2)**2-p(3)**2))
      end



      subroutine born_suppression(fact)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      logical ini
      data ini/.true./
      real * 8 fact,pt2,pt2supp,powheginput,pt
      save ini,pt2supp,pt     
c CAVEAT!!!  process dependent subroutine
      if (ini) then
         pt = powheginput("#bornsuppfact")         
         if (pt.gt.0d0) then
            write(*,*)
            write(*,*) "**** IMPROVED STATISTICS AT HIGH PT_W  ****"
            write(*,*)
         endif
         ini = .false.
         pt2supp = pt**2
      endif
      if (pt.le.0) then
         fact=1d0
      else         
c     improve statistics at high pt_W with this factor
         pt2=kn_pborn(1,3)**2+kn_pborn(2,3)**2
         fact=pt2/(pt2+pt2supp)         
      endif
      end


      subroutine set_fac_ren_scales(muf,mur)
      implicit none
      include 'PhysPars.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      real * 8 muf,mur
      muf=ph_Wmass+2*ph_mb
      mur=ph_Wmass+2*ph_mb
      end




