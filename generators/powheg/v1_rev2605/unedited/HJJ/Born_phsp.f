      subroutine born_phsp(xborn)
      implicit none
      include 'brinclude.h'
      include 'pwhg_kn.h'
      include 'coupl.inc'
      include 'PhysPars.h'
      real * 8 xborn(ndiminteg-3),jac,tmpvec(0:3),bwcutoff
      integer k
      logical ini,fullphsp
      data ini/.true./
      save ini,fullphsp
      real * 8 powheginput
      if(ini) then
c     set initial- and final-state masses for Born and real
         do k=1,nlegborn
            kn_masses(k)=0
         enddo
         kn_masses(nlegreal)=0
         kn_masses(3)=hmass
         ph_HmHw=hmass*hwidth
         ph_Hmass=hmass
         ph_Hwidth=hwidth
         ph_Hmass2=hmass**2
         bwcutoff=powheginput("bwcutoff")
         ph_Hmass2low=max(hmass-bwcutoff*hwidth,0d0)**2
         ph_Hmass2high=min(hmass+bwcutoff*hwidth,sqrt(kn_sbeams))**2
         if(powheginput("#fullphsp").eq.1) then
            fullphsp = .true.
         else
            fullphsp = .false.
         endif
         ini=.false.
      endif      

      brkn_ktmin=0
      call born_phsp_hj(xborn)
c      kn_masses(3)=sqrt(brkn_cmpborn(0,3)**2-brkn_cmpborn(1,3)**2
c     $     -brkn_cmpborn(2,3)**2-brkn_cmpborn(3,3)**2)

      if(fullphsp) then
         if(xborn(ndiminteg-3).lt.0.5d0) then
            brkn_emitter=4
            call br_real_phsp_fsr(xborn(ndiminteg-6),jac)
            jac=jac * 1/brkn_dijterm(4,5)
     1    /(1/brkn_dijterm(0,4)+1/brkn_dijterm(0,5)+1/brkn_dijterm(4,5))
     2        * brkn_cmpreal(0,4)/(brkn_cmpreal(0,4)+brkn_cmpreal(0,5))
            kn_cmpborn=brkn_cmpreal
            kn_pborn=brkn_preal
            kn_xb1=brkn_x1
            kn_xb2=brkn_x2
            if(xborn(ndiminteg-3).lt.0.25d0) then
               tmpvec=kn_cmpborn(:,4)
               kn_cmpborn(:,4)=kn_cmpborn(:,5)
               kn_cmpborn(:,5)=tmpvec
               tmpvec=kn_pborn(:,4)
               kn_pborn(:,4)=kn_pborn(:,5)
               kn_pborn(:,5)=tmpvec
            endif
         else
            brkn_emitter=0
            call br_real_phsp_isr(xborn(ndiminteg-6),jac)
            jac=jac * 1/brkn_dijterm(0,5)
     1    /(1/brkn_dijterm(0,4)+1/brkn_dijterm(0,5)+1/brkn_dijterm(4,5))
            kn_cmpborn=brkn_cmpreal
            kn_pborn=brkn_preal
            kn_xb1=brkn_x1
            kn_xb2=brkn_x2
            if(xborn(ndiminteg-3).lt.0.75d0) then
               tmpvec=kn_cmpborn(:,4)
               kn_cmpborn(:,4)=kn_cmpborn(:,5)
               kn_cmpborn(:,5)=tmpvec
               tmpvec=kn_pborn(:,4)
               kn_pborn(:,4)=kn_pborn(:,5)
               kn_pborn(:,5)=tmpvec
            endif
         endif
      else
         brkn_emitter=0
         call br_real_phsp_isr(xborn(ndiminteg-6),jac)
         kn_cmpborn=brkn_cmpreal
         kn_pborn=brkn_preal
         kn_xb1=brkn_x1
         kn_xb2=brkn_x2
         jac=jac/4
      endif
c factor of 4 to compensate for the 0.25 range for each
c of the above contributions.
      kn_jacborn=4*brkn_jacborn*jac

c     set the CMS energy 
      kn_sborn=brkn_sreal

c     minimal final state mass 
c     kn_minmass=sqrt(ph_Hmass2low)
c      kn_minmass=kn_ktmin + sqrt(kn_ktmin**2 + kn_masses(3)**2)
      kn_minmass=kn_ktmin + sqrt(kn_ktmin**2 + ph_Hmass2low)
      end

         

      subroutine born_phsp_hj(xborn)
      implicit none
      include 'brinclude.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      real * 8 xborn(ndiminteg-4)
      real * 8 m2,xjac,taumin,tau,y,beta,vec(3),cth,cthdec,phidec,s,
     # z,zhigh,zlow,dir(3),khiggs,cthmax,BW_fixed,BW_running
      real * 8 powheginput
      integer mu,k,j
      logical ini,higgsfixedwidth
      data ini/.true./
      save ini,higgsfixedwidth
      external powheginput
c Phase space:
c d m2/(2 pi) d omega /(8*(2 pi)^2)  (s-m2)/2
c omega: 3d angle in CM system
      if(ini) then
         higgsfixedwidth=powheginput("#higgsfixedwidth").gt.0
         ini=.false.
      endif
      zlow=atan((ph_Hmass2low  - ph_Hmass2)/ph_HmHw)
      zhigh=atan((min(ph_Hmass2high,kn_sbeams)  - ph_Hmass2)/ph_HmHw)
      z=zlow+(zhigh-zlow)*xborn(1)
      xjac=zhigh-zlow
      m2=ph_HmHw*tan(z)+ph_Hmass2
c     The BW integrates to Pi ==> divide by Pi
      xjac=xjac/pi

      if(.not.higgsfixedwidth) then
c     running width
         BW_fixed=ph_HmHw/((m2-ph_Hmass2)**2 + ph_HmHw**2)
         BW_running= (m2*ph_Hwidth/ph_Hmass) /
     $        ((m2-ph_Hmass2)**2+(m2*ph_Hwidth/ph_Hmass)**2)
         xjac = xjac * BW_running/BW_fixed
      endif

c     assign the Higgs boson mass
      kn_masses(3)=sqrt(m2)

c d x1 d x2 = d tau d y;
      taumin=( sqrt(m2+brkn_ktmin**2) + brkn_ktmin )**2/kn_sbeams
      tau=exp(log(taumin)*(1-xborn(2)**2))
      xjac=xjac*tau*abs(log(taumin))*2*xborn(2)
      s=kn_sbeams*tau
      brkn_sborn=s
c compute H momentum in partonic cm
      khiggs=(s-m2)/(2*sqrt(s))
c ymax=|log(tau)|/2
      y=-(1-2*xborn(3))*log(tau)/2
      xjac=-xjac*log(tau)
      cthmax=sqrt(1-(brkn_ktmin/khiggs)**2)
      z=1-2*xborn(4)
      xjac=xjac*2
      cth=1.5d0*(z-z**3/3)
      xjac=xjac*1.5d0*(1-z**2)
      cth=cth*cthmax
c      brkn_born_pt2=(1-cth**2)*khiggs**2
      xjac=xjac*cthmax
c supply 2 pi for azimuthal integration (not performed)
      xjac=xjac*2*pi
      xjac=xjac*(s-m2)/s/(8*(2*pi)**2)
c
      brkn_jacborn=xjac
c Build kinematics
c velocity of H in partonic CM
      brkn_cmpborn(1,3)=sqrt(1-cth**2)*khiggs
      brkn_cmpborn(2,3)=0
      brkn_cmpborn(3,3)=cth*khiggs
      brkn_cmpborn(0,3)=sqrt(m2+khiggs**2)
      brkn_xb1=sqrt(tau)*exp(y)
      brkn_xb2=tau/brkn_xb1
      brkn_cmpborn(1,4)=-brkn_cmpborn(1,3)
      brkn_cmpborn(2,4)=-brkn_cmpborn(2,3)
      brkn_cmpborn(3,4)=-brkn_cmpborn(3,3)
      brkn_cmpborn(0,4)=khiggs
      brkn_cmpborn(0,1)=sqrt(s)/2
      brkn_cmpborn(0,2)=brkn_cmpborn(0,1)
      brkn_cmpborn(3,1)=brkn_cmpborn(0,1)
      brkn_cmpborn(3,2)=-brkn_cmpborn(0,2)
      brkn_cmpborn(1,1)=0
      brkn_cmpborn(1,2)=0
      brkn_cmpborn(2,1)=0
      brkn_cmpborn(2,2)=0      
c now boost everything along 3
      beta=(brkn_xb1-brkn_xb2)/(brkn_xb1+brkn_xb2)
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(nlegborn-3,vec,beta,brkn_cmpborn(0,3),brkn_pborn(0,3))
      do mu=0,3
         brkn_pborn(mu,1)=brkn_xb1*kn_beams(mu,1)
         brkn_pborn(mu,2)=brkn_xb2*kn_beams(mu,2)
      enddo
c      call checkmomzero(nlegborn,brkn_pborn)
c      call checkmass(2,brkn_pborn(0,3))

c minimal final state mass 
c      brkn_minmass=sqrt(ph_Hmass2low)
      call br_compute_csimax_fsr
      end


c Mappings of the underlying born configuration in
c brkn_cmpborn(0:3,br_nlegborn), and the xrad(1:3) variables
c in the unit cube, into brkn_real(0:3,nlegreal).
c The factor jac_over_csi*csi*brkn_csimax, multiplied
c by the Born phase space jacobian, yields the real phase
c space jacobian.
c More explicitly:
c d Phi_n = d^3 xrad jac_over_csi csi csimax d Phi_{n-1}
c Since
c  d Phi_n = d phi d y d csi Jrad d Phi_{n-1}
c (where Jrad is given in FNO2006) we get
c                                  d phi d y d csi
c csimax csi jac_over_csi = Jrad  ----------------
c                                    d^3 xrad
c Notice that using d csi=d csitilde csimax the csimax
c factor cancels, and jac_over_csi is as given in the
c code below (see notes on xscaled.tm).
c br_real_phsp_fsr: provides the mapping for the final state
c radiation, assuming that the emitter is the brkn_emitter-th
c particle, and the emitted particle is the nlegreal-th particle
c br_real_phsp_isr: mapping for the initial state radiation
      subroutine br_real_phsp_fsr(xrad,jac)
      implicit none
      real * 8 xrad(3),jac
      include 'pwhg_math.h'
      include 'brinclude.h'
      real * 8 q0,q2,xjac
c Boost the underlying Born variables to their cm frame
      q0=2*brkn_cmpborn(0,1)
      q2=brkn_sborn
      brkn_csitilde=xrad(1)
      xjac=1
      brkn_y=1-2*xrad(2)
      xjac=xjac*2
c importance sampling for brkn_y
      xjac=xjac*1.5d0*(1-brkn_y**2)
      brkn_y=1.5d0*(brkn_y-brkn_y**3/3)
      brkn_azi=2*pi*xrad(3)
      xjac=xjac*2*pi
      brkn_csimax=brkn_csimax_arr(brkn_emitter)
      brkn_csi=brkn_csitilde*brkn_csimax     
c remember: no csimax in the jacobian factor, we are integrating in csitilde 
      call br_real_phsp_fsr_rad
      jac=xjac*brkn_jacreal*brkn_csimax
      end

c br_real_phsp_fsr_rad: provides the mapping for the final state
c radiation, assuming that we are considering the region rad_kinreg
c and the emitted particle is the nlegreal-th particle,
c for given brkn_csi, brkn_y, brkn_azi. Sets the jacobian
c brkn_jacreal so that brkn_jacreal d brkn_csi d brkn_y d brkn_azi times
c the underlying Born jacobian is the phase space volume
      subroutine br_real_phsp_fsr_rad
      implicit none
      include 'pwhg_math.h'
      include 'brinclude.h'
      real * 8 vec(3),q0,beta
      integer i
      data vec/0d0,0d0,1d0/
      save vec
      q0=2*brkn_cmpborn(0,1)
c remember: no csimax factor, we are integrating in csitilde 
      call barradmap(br_nlegborn-2,brkn_emitter-2,q0,brkn_cmpborn(0,3),
     1    brkn_csi,brkn_y,brkn_azi,brkn_preal(0,3),brkn_jacreal)
      beta=(brkn_xb1-brkn_xb2)/(brkn_xb1+brkn_xb2)
      call mboost(br_nlegreal-2,vec,beta,
     1    brkn_preal(0,3),brkn_preal(0,3))
      do i=0,3
         brkn_preal(i,1)=brkn_pborn(i,1)
         brkn_preal(i,2)=brkn_pborn(i,2)
      enddo
      brkn_x1=brkn_xb1
      brkn_x2=brkn_xb2
      brkn_sreal=brkn_sborn
c      call checkmomzero(br_nlegreal,brkn_preal)
      call br_compcmkin
      call br_compdij
      end

      subroutine br_real_phsp_isr(xrad,jac)
      implicit none
      real * 8 xrad(3),jac
      include 'pwhg_math.h'
      include 'brinclude.h'
      real * 8 xjac
      brkn_csitilde=(3-2*xrad(1))*xrad(1)**2
      xjac=6*(1-xrad(1))*xrad(1)
      brkn_y=1-2*xrad(2)
      xjac=xjac*2
      xjac=xjac*1.5d0*(1-brkn_y**2)
      brkn_y=1.5d0*(brkn_y-brkn_y**3/3)
      brkn_azi=2*pi*xrad(3)
      xjac=xjac*2*pi
      call br_compcsimax
      brkn_csi=brkn_csitilde*brkn_csimax
      call br_real_phsp_isr_rad
      jac=xjac*brkn_jacreal*brkn_csimax
      end

      subroutine br_compcsimax
      implicit none
      include 'brinclude.h'
      real * 8 y,xb1,xb2
      xb1=brkn_xb1
      xb2=brkn_xb2
      y=brkn_y
      brkn_csimax=1-max(2*(1+y)*xb1**2/
     1    (sqrt((1+xb1**2)**2*(1-y)**2+16*y*xb1**2)+(1-y)*(1-xb1**2)),
     1            2*(1-y)*xb2**2/
     1    (sqrt((1+xb2**2)**2*(1+y)**2-16*y*xb2**2)+(1+y)*(1-xb2**2)))
      end

      subroutine br_real_phsp_isr_rad
      implicit none
      include 'pwhg_math.h'
      include 'brinclude.h'
      include 'pwhg_kn.h'
      real * 8 y,xb1,xb2,x1,x2,betal,betat,vecl(3),vect(3),
     1         cth,sth,cph,sph,csi,pt2
      integer i,mu
      real * 8 dotp
      external dotp
c the following call sets brkn_csimax, brkn_csimaxp, brkn_csimaxm
c also when br_real_phsp_isr_rad is called directly
c (i.e. not through br_real_phsp_isr_rad0)
      call br_compcsimax
      y=brkn_y
      xb1=brkn_xb1
      xb2=brkn_xb2
      csi=brkn_csi
      cth=y
      sth=sqrt(1-cth**2)
      cph=cos(brkn_azi)
      sph=sin(brkn_azi)
      x1=xb1/sqrt(1-csi)*sqrt((2-csi*(1-y))/(2-csi*(1+y)))
      x2=xb2/sqrt(1-csi)*sqrt((2-csi*(1+y))/(2-csi*(1-y)))
      brkn_x1=x1
      brkn_x2=x2
      do mu=0,3
         brkn_preal(mu,1)=kn_beams(mu,1)*x1
         brkn_preal(mu,2)=kn_beams(mu,2)*x2
      enddo
      brkn_sreal=brkn_sborn/(1-csi)
c Build k_n+1 in the rest frame of brkn_preal
c      write(*,*) ' br_nlegreal ',br_nlegreal
      brkn_preal(0,br_nlegreal)=sqrt(brkn_sreal)*csi/2
      brkn_preal(1,br_nlegreal)=brkn_preal(0,br_nlegreal)*sth*sph
      brkn_preal(2,br_nlegreal)=brkn_preal(0,br_nlegreal)*sth*cph
      brkn_preal(3,br_nlegreal)=brkn_preal(0,br_nlegreal)*cth
c boost it to the frame of brkn_preal
      do i=1,3
         vecl(i)=(brkn_preal(i,1)+brkn_preal(i,2))
     1          /(brkn_preal(0,1)+brkn_preal(0,2))
      enddo      
      betal=sqrt(vecl(1)**2+vecl(2)**2+vecl(3)**2)
      if(betal.gt.0) then
         do i=1,3
            vecl(i)=vecl(i)/betal
         enddo
      else
         vecl(1)=1
         vecl(2)=0
         vecl(3)=0
      endif
      call mboost(1,vecl,betal,
     1    brkn_preal(0,br_nlegreal),brkn_preal(0,br_nlegreal))
c longitudinal boost of underlying Born to zero rapidity frame
      do i=1,3
         vecl(i)=(brkn_pborn(i,1)+brkn_pborn(i,2))
     1          /(brkn_pborn(0,1)+brkn_pborn(0,2))
      enddo
      betal=sqrt(vecl(1)**2+vecl(2)**2+vecl(3)**2)
      if(betal.gt.0) then
         do i=1,3
            vecl(i)=vecl(i)/betal
         enddo
      else
         vecl(1)=1
         vecl(2)=0
         vecl(3)=0
      endif
      call mboost(br_nlegborn-2,vecl,-betal,
     1 brkn_pborn(0,3),brkn_preal(0,3))
c      call printtot(br_nlegborn,brkn_preal(0,1))
c construct transverse boost velocity
      vect(3)=0
      vect(1)=brkn_preal(1,br_nlegreal)
      vect(2)=brkn_preal(2,br_nlegreal)
      pt2=vect(1)**2+vect(2)**2
c      betat=1/sqrt(1+(brkn_sreal*(1-csi))/pt2)
      betat=sqrt(pt2/(pt2+brkn_sreal*(1-csi)))
      if(pt2.eq.0) then
         vect(1)=1
         vect(2)=0
      else
         vect(1)=vect(1)/sqrt(pt2)
         vect(2)=vect(2)/sqrt(pt2)
      endif
c     write(*,*) ' k+1: ',(brkn_preal(mu,br_nlegreal),mu=0,3)
         call mboost(br_nlegborn-2,vect,-betat,
     1        brkn_preal(0,3),brkn_preal(0,3))
c      call printtot(nlegborn,brkn_preal(0,1))
c longitudinal boost in opposite direction
      call mboost(br_nlegborn-2,vecl,betal,
     1 brkn_preal(0,3),brkn_preal(0,3))
c      call printtot(br_nlegreal,brkn_preal(0,1))
      brkn_jacreal=brkn_sreal/(4*pi)**3*csi/(1-csi)
      call br_compcmkin
      call br_compdij
      end


      subroutine br_compcmkin
      implicit none
      include 'brinclude.h'
      real * 8 vecl(3),betal
      data vecl/0d0,0d0,1d0/
      save vecl
      betal=-(brkn_preal(3,1)+brkn_preal(3,2))/
     1 (brkn_preal(0,1)+brkn_preal(0,2))
      call mboost(br_nlegreal,vecl,betal,brkn_preal,brkn_cmpreal)
      end

      subroutine br_compdij
      implicit none
      include 'brinclude.h'
      integer j,k
      real * 8 y
      real * 8 crossp,dotp
      external crossp,dotp
      do j=4,br_nlegreal
         y=1-dotp(brkn_cmpreal(0,1),brkn_cmpreal(0,j))
     1 /(brkn_cmpreal(0,1)*brkn_cmpreal(0,j))
         brkn_dijterm(0,j)=(brkn_cmpreal(0,j)**2
     1 *(1-y**2))**brpar_diexp
         brkn_dijterm(1,j)=(brkn_cmpreal(0,j)**2
     1 *2*(1-y))**brpar_diexp
         brkn_dijterm(2,j)=(brkn_cmpreal(0,j)**2
     1 *2*(1+y))**brpar_diexp
      enddo
      do j=4,br_nlegreal
         do k=j+1,br_nlegreal
            brkn_dijterm(j,k)=
     1(2*dotp(brkn_cmpreal(0,k),brkn_cmpreal(0,j))*
     1       brkn_cmpreal(0,k)*brkn_cmpreal(0,j)
     2    /  (brkn_cmpreal(0,k)+brkn_cmpreal(0,j))**2)**brpar_dijexp
c     2    /  ((brkn_cmpreal(1,k)+brkn_cmpreal(1,j))**2+
c     3        (brkn_cmpreal(2,k)+brkn_cmpreal(2,j))**2+
c     4        (brkn_cmpreal(3,k)+brkn_cmpreal(3,j))**2))**brpar_dijexp
         enddo
      enddo
      end

      subroutine br_compute_csimax_fsr
      implicit none
c Compute csimax for all possible final state emitters;
c for initial state emitters it is not possible, since
c csimax depends upon y in this case.
      include 'brinclude.h'
      integer j
      real * 8 q0,mrec2
      logical valid_emitter
      external valid_emitter
      j=4
      q0=2*brkn_cmpborn(0,1)
      mrec2=(q0-brkn_cmpborn(0,j))**2
     1     -brkn_cmpborn(1,j)**2-brkn_cmpborn(2,j)**2
     1     -brkn_cmpborn(3,j)**2
      brkn_csimax_arr(j)=1-mrec2/q0**2
      end



      subroutine born_suppression(fact)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      include 'coupl.inc'
      real * 8 fact,ptmin
      real * 8 pt4,pt5,pt45,dotp
      real * 8 powheginput
      logical ini
      data ini/.true./
      save ptmin,ini
      if(ini) then
         ptmin=powheginput('#bornsuppfact')
         ini=.false.
      endif
      if(ptmin.gt.0) then
         pt4=kn_cmpborn(1,4)**2+kn_cmpborn(2,4)**2
         pt5=kn_cmpborn(1,5)**2+kn_cmpborn(2,5)**2
         pt45=2*dotp(kn_cmpborn(0,4),kn_cmpborn(0,5))
     1     *  kn_cmpborn(0,4)*kn_cmpborn(0,5)
     2        /(kn_cmpborn(0,4)**2+kn_cmpborn(0,5)**2)
         fact=((1/ptmin**2)/(1/pt4+1/pt5+1/pt45+1/ptmin**2))**2
      else
         fact=1
      endif
      end



      subroutine set_fac_ren_scales(muf,mur)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'coupl.inc'
      real * 8 muf,mur
      logical ini
      data ini/.true./
      logical runningscales
      real * 8 pt1,pt2,ptHsq,Ht
      real * 8 powheginput
      external powheginput      
      save ini,runningscales


      if (ini) then
         if (powheginput("#runningscales").eq.1) then
            runningscales=.true.
            write(*,*) '****************************************'
            write(*,*) '****************************************'
            write(*,*) '**   mur=Ht  used for Bbar function   **'
            write(*,*) '**   muf=Ht  used for Bbar function   **'
            write(*,*) '****************************************'
            write(*,*) '****************************************'
         else
            runningscales=.false.
         endif
         ini=.false.
      endif
      
      if (runningscales) then
         pt1=sqrt(kn_pborn(1,4)**2+kn_pborn(2,4)**2)
         pt2=sqrt(kn_pborn(1,5)**2+kn_pborn(2,5)**2)
         ptHsq=kn_pborn(1,3)**2+kn_pborn(2,3)**2
         Ht=sqrt(hmass**2+ptHsq)+pt1+pt2
         mur=Ht
         muf=mur
      else
         muf=hmass
         mur=hmass
      endif
      end



