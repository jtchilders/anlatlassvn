      subroutine KMampsqggQQbdkn(n,in,q1,q2,q3,q4,q5,q6,q7,q8,res)
      implicit none
c----notation                                                       *
C     0--> g(q1)+g(q2)+nu(q3)+e+(q4)+b(q5)+bbar(q6)+e-(q7)+nubar(q8)
c---
c--- Computed using the helicity amplitudes from:
c---
c--- \bibitem{Korner:2002hy}
c--- J.~G.~Korner and Z.~Merebashvili,
c--- %``One-loop corrections to four-point functions with two external massive
c--- %fermions and two external massless partons,''
c--- Phys.\ Rev.\  D {\bf 66}, 054023 (2002)
c--- [arXiv:hep-ph/0207054].
c--- %%CITATION = PHRVA,D66,054023;%%
      include 'constants.f'
      include 'masses.f'
      include 'spinorsw.f'
      integer j,nu,in,icol
      double precision q1(4),q2(4),q3(4),q4(4),q5(4),q6(4),q7(4),q8(4),
     & n(4),k1(4),k2(4),k3(4),k4(4),k5(4),k6(4),k7(4),k8(4),kn(4),
     & res(0:2)
      double complex zp1(4),zp2(4),zp3(4),zp4(4),ze1(4),ze2(4),zn(4),
     & tab1,tab2,tab3,Ubm(4),Vm(4),Bs,Bt,Bu,loab,loba,loqed
      
c--- set up common variable to indicate KM routines should be used
c--- to calculate spinor strings
      spinorsw='KM'      

c--- these routines expect energy in the 1st component and, as usual,
c--- we modifiy this solution by exchanging (px <--> pz) and changing 
C--- the sign of py, (to remain in right handed system) [c.f. spinor_w.f]
      do j=1,4
      if (j .eq. 1) nu=4
      if (j .eq. 2) nu=3
      if (j .eq. 3) nu=2
      if (j .eq. 4) nu=1
      k1(j)=q1(nu)
      k2(j)=q2(nu)
      k3(j)=q3(nu)
      k4(j)=q4(nu)
      k5(j)=q5(nu)
      k6(j)=q6(nu)
      k7(j)=q7(nu)
      k8(j)=q8(nu)
      kn(j)=n(nu)
      if (j .eq. 3) then
      k1(j)=-k1(j)
      k2(j)=-k2(j)
      k3(j)=-k3(j)
      k4(j)=-k4(j)
      k5(j)=-k5(j)
      k6(j)=-k6(j)
      k7(j)=-k7(j)
      k8(j)=-k8(j)
      kn(j)=-kn(j)
      endif
      enddo

C-----remember that diagrams are Qbar(k1)+Q(k2)+q(k3)+qbar(k4)
C     0--> qbar(q1)+q(q2)+Q(q345)+Qb(q678)
C     0--> qbar(k4)+q(k3)+Q(k2)+Qb(k1)
      do nu=1,4
      zp1(nu)=dcmplx(k1(nu))
      zp2(nu)=dcmplx(k2(nu))
      zp3(nu)=dcmplx(k3(nu)+k4(nu)+k5(nu))
      zp4(nu)=dcmplx(k6(nu)+k7(nu)+k8(nu))
      zn(nu)=dcmplx(kn(nu))
      enddo

      call ubarWspinor(k3,k4,k5,mt,wmass,Ubm)
      call uWspinor(k6,k7,k8,mt,wmass,Vm)

      do icol=0,2
      res(icol)=zip
      enddo

      do j=-1,1,2

        if     (in .eq. 1) then
	  call pol_mless(zp2,j,ze2)
          tab1=Bs(zp1,zp2,zp3,zp4,zn,ze2,Ubm,Vm)
          tab3=Bt(zp1,zp2,zp3,zp4,zn,ze2,Ubm,Vm)
          tab2=Bu(zp1,zp2,zp3,zp4,zn,ze2,Ubm,Vm)
        elseif (in .eq. 2) then
	  call pol_mless(zp1,j,ze1)
          tab1=Bs(zp1,zp2,zp3,zp4,ze1,zn,Ubm,Vm)
          tab3=Bt(zp1,zp2,zp3,zp4,ze1,zn,Ubm,Vm)
          tab2=Bu(zp1,zp2,zp3,zp4,ze1,zn,Ubm,Vm)
        else
          write(6,*) 'KMampsqggttbn: Unimplemented value of in',in
          stop
        endif

        loab=tab2+tab1
        loba=tab3-tab1
        loqed=tab2+tab3

        res(1)=res(1)+V*xn/4d0*abs(loab)**2
        res(2)=res(2)+V*xn/4d0*abs(loba)**2
        res(0)=res(0)-V/xn/4d0*abs(loqed)**2

      enddo
      
      return
      end
      

