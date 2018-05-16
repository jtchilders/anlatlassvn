      double complex function Bs(zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm)
      implicit none
      include 'constants.f'
c----- \bibitem{Korner:2002hy}
c----- J.~G.~Korner and Z.~Merebashvili,
c----- %``One-loop corrections to four-point functions with two external massive
c----- %fermions and two external massless partons,''
c----- Phys.\ Rev.\  D {\bf 66}, 054023 (2002)
c----- [arXiv:hep-ph/0207054].
c----- %%CITATION = PHRVA,D66,054023;%%

C--- This is an implementation of KM, Eq. (2.6)
c--- Note: need to make sure to switch sign of contributions containing
c--- p1 or p2, due to difference in incoming/outgoing momenta

c--- Note that this is proportional to the colour factor (Ta*Tb-Tb*Ta)
      double complex zp1(4),zp2(4),zp3(4),zp4(4),s
      double complex e1(4),e2(4),Ubm(4),Vm(4),cdot,string1
      s=2d0*cdot(zp1,zp2)
      Bs=2d0*im/s*(
     &  -cdot(e1,e2)*string1(Ubm,zp1,Vm)
     &  -cdot(e1,zp2)*string1(Ubm,e2,Vm)
     &  +cdot(e2,zp1)*string1(Ubm,e1,Vm))
      return
      end
