!!
!! Copyright (C) 2011 Andreas van Hameren. 
!!
!! This file is part of OneLOop-2.2.1.
!!
!! OneLOop-2.2.1 is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! OneLOop-2.2.1 is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with OneLOop-2.2.1.  If not, see <http://www.gnu.org/licenses/>.
!!


module avh_olo_tri
  use avh_olo_kinds
  use avh_olo_units
  use avh_olo_func
  implicit none
  private
  public :: tria0,tria1,tria2,tria3,tria4,trif0,trif1,trif2,trif3 &
           ,trif3HV &
           ,permtable,casetable,base
  complex(kindc2) ,parameter :: oieps=C1P0+CiP0*epsilon(R1P0)**2
  integer ,parameter :: permtable(3,0:7)=reshape((/ &
       1,2,3 &! 0, 0 masses non-zero, no permutation
      ,1,2,3 &! 1, 1 mass non-zero,   no permutation
      ,3,1,2 &! 2, 1 mass non-zero,   1 cyclic permutation
      ,1,2,3 &! 3, 2 masses non-zero, no permutation
      ,2,3,1 &! 4, 1 mass non-zero,   2 cyclic permutations
      ,2,3,1 &! 5, 2 masses non-zero, 2 cyclic permutations
      ,3,1,2 &! 6, 2 masses non-zero, 1 cyclic permutation
      ,1,2,3 &! 7, 3 masses non-zero, no permutation
      /) ,(/3,8/))                     ! 0,1,2,3,4,5,6,7
  integer ,parameter :: casetable(0:7)=(/0,1,1,2,1,2,2,3/)
  integer ,parameter :: base(3)=(/4,2,1/)

contains

   subroutine tria4( rslt ,cpp,cm2,cm3 ,rmu2 )
!*******************************************************************
! calculates
!               C   /             d^(Dim)q
!            ------ | ----------------------------------
!            i*pi^2 / q^2 [(q+k1)^2-m2] [(q+k1+k2)^2-m3]
!
! with  k1^2=m2, k2^2=pp, (k1+k2)^2=m3.
! m2,m3 should NOT be identically 0d0.
!*******************************************************************
   use avh_olo_logc ,only: logc
   use avh_olo_logc2 ,only: logc2
   use avh_olo_li2c2 ,only: li2c2
   complex(kindc2) ,intent(out) :: rslt(0:2)
   complex(kindc2) ,intent(in)  :: cm2,cm3,cpp
   real(kindr2)    ,intent(in)  :: rmu2
   type(qmplx_type) :: q23,qm3,q32
   complex(kindc2) :: sm2,sm3,k23,r23,d23,cc
!
!  write(*,*) 'MESSAGE from OneLOop tria4: you are calling me' !CALLINGME
!
   sm2 = mysqrt(cm2)
   sm3 = mysqrt(cm3)
   k23 = (cm2+cm3-cpp)/(sm2*sm3)
   call rfun( r23,d23, k23 )
   if (r23.eq.-C1P0) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop tria4: ' &
       ,'threshold singularity, returning 0'
     rslt = C0P0
     return
   endif
   q23 = qonv(r23,-1)
   qm3 = qonv(cm3/rmu2,-1)
   q32 = qonv(sm3)/qonv(sm2)
!
   rslt(2) = C0P0
   cc = logc2(q23) * r23/(C1P0+r23)/(sm2*sm3)
   rslt(1) = -cc
   rslt(0) = cc*( logc(qm3) - logc(q23) ) &
           - li2c2(q32*q23,q32/q23) / cm2 &
           + li2c2(q23*q23,qonv(C1P0)) * r23/(sm2*sm3)
   end subroutine


   subroutine tria3( rslt ,cp2,cp3,cm3 ,rmu2 )
!*******************************************************************
! calculates
!               C   /          d^(Dim)q
!            ------ | -----------------------------
!            i*pi^2 / q^2 (q+k1)^2 [(q+k1+k2)^2-m3]
!
! with  p2=k2^2, p3=(k1+k2)^2.
! mm should NOT be identically 0d0,
! and p2 NOR p3 should be identical to mm. 
!*******************************************************************
   use avh_olo_logc ,only: logc
   use avh_olo_logc2 ,only: logc2
   use avh_olo_li2c2 ,only: li2c2
   complex(kindc2) ,intent(out) :: rslt(0:2)
   complex(kindc2) ,intent(in)  :: cp2,cp3,cm3
   real(kindr2)    ,intent(in)  :: rmu2
   type(qmplx_type) :: q13,q23,qm3,x1,x2
   complex(kindc2) :: r13,r23
!
!  write(*,*) 'MESSAGE from OneLOop tria3: you are calling me' !CALLINGME
!
   r13 = cm3-cp3
   r23 = cm3-cp2
   q13 = qonv(r13,-1)
   q23 = qonv(r23,-1)
   qm3 = qonv(cm3,-1)
   x1 = q23/qm3
   x2 = q13/qm3
   rslt(2) = C0P0
   rslt(1) = -logc2( q23/q13 )/r13
   rslt(0) = -li2c2( x1,x2 )/cm3 &
           - rslt(1)*( logc(x1*x2)+logc(qm3/rmu2) )
   end subroutine


   subroutine tria2( rslt ,cp3,cm3 ,rmu2 )
!*******************************************************************
! calculates
!               C   /          d^(Dim)q
!            ------ | -----------------------------
!            i*pi^2 / q^2 (q+k1)^2 [(q+k1+k2)^2-m3]
!
! with  k1^2 = 0 , k2^2 = m3  and  (k1+k2)^2 = p3.
! mm should NOT be identically 0d0,
! and pp should NOT be identical to mm. 
!*******************************************************************
   use avh_olo_logc ,only: logc
   use avh_olo_li2c ,only: li2c
   complex(kindc2) ,intent(out) :: rslt(0:2)
   complex(kindc2) ,intent(in)  :: cp3,cm3
   real(kindr2)    ,intent(in)  :: rmu2
   complex(kindc2) ,parameter :: const=C1P0*TWOPI*TWOPI/96
   complex(kindc2) ,parameter :: half=C1P0/2
   type(qmplx_type) :: q13,qm3,qxx
   complex(kindc2) :: r13,logm,z2,z1,z0,cc
!
!  write(*,*) 'MESSAGE from OneLOop tria2: you are calling me' !CALLINGME
!
   r13 = cm3-cp3
   q13 = qonv(r13,-1)
   qm3 = qonv(cm3,-1)
   logm = logc( qm3/rmu2 )
   qxx = qm3/q13
   z2 = half
   z1 = logc(qxx)
   z0 = const + z1*z1/2 - li2c(qxx)
   cc = -C1P0/r13
   rslt(2) = cc*z2
   rslt(1) = cc*(z1 - z2*logm)
   rslt(0) = cc*(z0 + (z2*logm/2-z1)*logm)
   end subroutine


   subroutine tria1( rslt ,cm3 ,rmu2 )
!*******************************************************************
! calculates
!               C   /          d^(Dim)q
!            ------ | -----------------------------
!            i*pi^2 / q^2 (q+k1)^2 [(q+k1+k2)^2-m3]
!
! with  k1^2 = (k1+k2)^2 = m3.
! mm should NOT be identically 0d0.
!*******************************************************************
   use avh_olo_logc ,only: logc
   complex(kindc2) ,intent(out) :: rslt(0:2)
   complex(kindc2) ,intent(in)  :: cm3
   real(kindr2)    ,intent(in)  :: rmu2
   complex(kindc2) :: zm
!
!  write(*,*) 'MESSAGE from OneLOop tria1: you are calling me' !CALLINGME
!
   zm = C1P0/(2*cm3)
   rslt(2) = C0P0
   rslt(1) = -zm
   rslt(0) = zm*( 2*C1P0 + logc(qonv(cm3/rmu2,-1)) )
   end subroutine


   subroutine tria0( rslt ,cp ,ap ,rmu2 )
!*******************************************************************
! calculates
!               C   /         d^(Dim)q
!            ------ | ------------------------
!            i*pi^2 / q^2 (q+k1)^2 (q+k1+k2)^2
!
! with  Dim = 4-2*eps
!         C = pi^eps * mu^(2*eps) * exp(gamma_Euler*eps)
!
! input:  p1 = k1^2,  p2 = k2^2,  p3 = k3^2
! output: rslt(0) = eps^0   -coefficient
!         rslt(1) = eps^(-1)-coefficient
!         rslt(2) = eps^(-2)-coefficient
!
! If any of these numbers is IDENTICALLY 0d0, the corresponding
! IR-singular case is returned.
!*******************************************************************
   use avh_olo_loga ,only: loga
   use avh_olo_loga2 ,only: loga2
   complex(kindc2) ,intent(out) :: rslt(0:2)
   complex(kindc2) ,intent(in)  :: cp(3)
   real(kindr2)    ,intent(in)  :: ap(3),rmu2
   real(kindr2) :: pp(3),rp1,rp2,rp3
   complex(kindc2) ,parameter :: const=C1P0*TWOPI*TWOPI/48
   complex(kindc2) :: log2,log3
   integer :: icase,i1,i2,i3
!
   pp(1)=real(cp(1))
   pp(2)=real(cp(2))
   pp(3)=real(cp(3))
!
   icase = 0
   if (ap(1).gt.R0P0) icase = icase + base(1)
   if (ap(2).gt.R0P0) icase = icase + base(2)
   if (ap(3).gt.R0P0) icase = icase + base(3)
   rp1 = pp(permtable(1,icase))
   rp2 = pp(permtable(2,icase))
   rp3 = pp(permtable(3,icase))
   icase  = casetable(  icase)
!
   i1=0 ;if (-rp1.lt.R0P0) i1=-1
   i2=0 ;if (-rp2.lt.R0P0) i2=-1
   i3=0 ;if (-rp3.lt.R0P0) i3=-1
!
   if     (icase.eq.0) then
! 0 masses non-zero
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop tria0: ' &
       ,'all external masses equal zero, returning 0'
     rslt = C0P0
   elseif (icase.eq.1) then
! 1 mass non-zero
!  write(*,*) 'MESSAGE from OneLOop tria0 1: you are calling me' !CALLINGME
    log3 = loga( -rp3/rmu2 , i3 )
    rslt(2) = cmplx( R1P0/rp3 ,kind=kindc2 )
    rslt(1) = -log3/rp3
    rslt(0) = ( log3**2/2 - const )/rp3
  elseif (icase.eq.2) then
! 2 masses non-zero
!  write(*,*) 'MESSAGE from OneLOop tria0 2: you are calling me' !CALLINGME
    log2 = loga( -rp2/rmu2 ,i2 )
    log3 = loga( -rp3/rmu2 ,i3 )
    rslt(2) = C0P0
    rslt(1) = loga2( rp3/rp2 ,i3-i2 )/rp2
    rslt(0) = -rslt(1)*(log3+log2)/2
  elseif (icase.eq.3) then
! 3 masses non-zero
    call trif0( rslt ,cp(1),cp(2),cp(3) )
  endif
  end subroutine


   subroutine trif0( rslt ,p1,p2,p3 )
!*******************************************************************
! Finite 1-loop scalar 3-point function with all internal masses
! equal zero. Obtained from the formulas for 4-point functions in
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
! by sending one internal mass to infinity.
!*******************************************************************
   use avh_olo_logc ,only: logc
   use avh_olo_logc2 ,only: logc2
   use avh_olo_li2c2 ,only: li2c2
   complex(kindc2) ,intent(out) :: rslt(0:2)
   complex(kindc2) ,intent(in)  :: p1,p2,p3
   type(qmplx_type) :: q23,q24,q34,qx1,qx2
   complex(kindc2) :: r23,r24,r34,aa,bb,cc,dd,x1,x2
   real(kindr2) :: hh
!
!  write(*,*) 'MESSAGE from OneLOop trif0: you are calling me' !CALLINGME
!
   r23 = -p1
   r24 = -p3
   r34 = -p2
!
   aa = r34*r24
   bb = r24 + r34 - r23
   cc = C1P0
   hh = real(r23)
   dd = mysqrt( bb*bb - 4*aa*cc , -real(aa)*hh )
   call solabc( x1,x2,dd ,aa,bb,cc ,1 )
   x1 = -x1
   x2 = -x2
!
   qx1 = qonv(x1, hh)
   qx2 = qonv(x2,-hh)
   q23 = qonv(r23,-1)
   q24 = qonv(r24,-1)
   q34 = qonv(r34,-1)
!
   rslt = C0P0
!
   rslt(0) = li2c2( qx1*q34 ,qx2*q34 )*r34 &
           + li2c2( qx1*q24 ,qx2*q24 )*r24 &
           - logc2( qx1/qx2 )*logc( qx1*qx2 )/(x2*2) &
           - logc2( qx1/qx2 )*logc( q23 )/x2
!
   rslt(0) = rslt(0)/aa
   end subroutine


   subroutine trif1( rslt ,p1i,p2i,p3i ,m3i )
!*******************************************************************
! Finite 1-loop scalar 3-point function with one internal masses
! non-zero. Obtained from the formulas for 4-point functions in
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
! by sending one internal mass to infinity.
!*******************************************************************
   use avh_olo_logc ,only: logc
   use avh_olo_logc2 ,only: logc2
   use avh_olo_li2c2 ,only: li2c2
   complex(kindc2) ,intent(out) :: rslt(0:2)
   complex(kindc2) ,intent(in)  :: p1i,p2i,p3i ,m3i 
   type(qmplx_type) :: q23,q24,q34,qm4,qx1,qx2,qss
   complex(kindc2) :: p2,p3,p4,p12,p23,m4,sm2,sm3,sm4 &
                     ,aa,bb,cc,dd,x1,x2,r23,r24,r34
   real(kindr2) :: mhh
   complex(kindc2) ,parameter :: oieps=C1P0+CiP0*epsilon(R1P0)**2
!
!  write(*,*) 'MESSAGE from OneLOop trif1: you are calling me' !CALLINGME
!
!   p1 = nul
   p2 = p1i
   p3 = p2i
   p4 = p3i
   p12 = p1i
   p23 = p3i
!   m1 = infinite
!   m2 = m1i = C0P0
!   m3 = m2i = C0P0
   m4 = m3i
!
   sm4 = mysqrt(m4)
   mhh = abs(sm4)
   sm3 = cmplx(mhh,kind=kindc2)
   sm2 = sm3
!
   r24 = C0P0
   r34 = C0P0
                  r23 = (   -p2 *oieps )/(sm2*sm3)
   if (m4.ne.p23) r24 = ( m4-p23*oieps )/(sm2*sm4)
   if (m4.ne.p3 ) r34 = ( m4-p3 *oieps )/(sm3*sm4)     
!
   aa = r34*r24 - r23
!
   if (aa.eq.C0P0) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop trif1: ' &
       ,'threshold singularity, returning 0'
     rslt = C0P0
     return
   endif
!
   bb = r24/sm3 + r34/sm2 - r23/sm4
   cc = C1P0/(sm2*sm3)
!   hh = real(r23)
!   dd = mysqrt( bb*bb - 4*aa*cc , -real(aa)*hh )
   call solabc( x1,x2,dd ,aa,bb,cc ,0 )
   x1 = -x1
   x2 = -x2
!
   qx1 = qonv(x1 ,1) ! x1 SHOULD HAVE im. part
   qx2 = qonv(x2 ,1) ! x2 SHOULD HAVE im. part
   q23 = qonv(r23,-1)
   q24 = qonv(r24,-1)
   q34 = qonv(r34,-1)
   qm4 = qonv(sm4,-1)
!
   rslt = C0P0
!
   rslt(0) = -logc2( qx1/qx2 )*logc( qx1*qx2/(qm4*qm4) )/(x2*2) &
             -li2c2( qx1*qm4 ,qx2*qm4 )*sm4
!
   if (r34.ne.C0P0) then
     qss = q34*mhh
     rslt(0) = rslt(0) + li2c2( qx1*qss ,qx2*qss )*r34*sm3
   endif
!
   if (r24.ne.C0P0) then
     qss = q24*mhh
     rslt(0) = rslt(0) + li2c2( qx1*qss ,qx2*qss )*r24*sm2
   endif
!
   rslt(0) = rslt(0) - logc2( qx1/qx2 )*logc( q23*(mhh*mhh) )/x2
!
   rslt(0) = rslt(0)/(aa*sm2*sm3*sm4)
   end subroutine


   subroutine trif2( rslt ,p1i,p2i,p3i ,m2i,m3i )
!*******************************************************************
! Finite 1-loop scalar 3-point function with two internal masses
! non-zero. Obtained from the formulas for 4-point functions in
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
! by sending one internal mass to infinity.
!*******************************************************************
   use avh_olo_logc ,only: logc
   use avh_olo_logc2 ,only: logc2
   use avh_olo_li2c2 ,only: li2c2
   complex(kindc2) ,intent(out) :: rslt(0:2)
   complex(kindc2) ,intent(in)  :: p1i,p2i,p3i ,m2i,m3i
   type(qmplx_type) :: q23,q34,q24,qm2,qm3,qm4,qx1,qx2,qss,qy1,qy2
   complex(kindc2) :: p2,p3,p23,m2,m4,sm2,sm3,sm4,aa,bb,cc,dd,x1,x2 &
                     ,r23,k24,r34,r24,d24
  complex(kindc2) ,parameter :: oieps=C1P0+CiP0*epsilon(R1P0)**2
!
!  write(*,*) 'MESSAGE from OneLOop trif2: you are calling me' !CALLINGME
!
!   p1 = nul
   p2 = p3i
   p3 = p1i
!   p4 = p2i
!   p12 = p3i
   p23 = p2i
!   m1 = infinite
   m2 = m3i
!   m3 = m1i = C0P0
   m4 = m2i
!
!   sm1 = infinite
   sm2 = mysqrt(m2)
   sm3 = cmplx(abs(sm2),kind=kindc2) !mysqrt(m3)
   sm4 = mysqrt(m4)
!
   r23 = C0P0
   k24 = C0P0
   r34 = C0P0
   if (m2   .ne.p2 ) r23 = (    m2-p2 *oieps )/(sm2*sm3) ! p2
   if (m2+m4.ne.p23) k24 = ( m2+m4-p23*oieps )/(sm2*sm4) ! p2+p3
   if (m4   .ne.p3 ) r34 = (    m4-p3 *oieps )/(sm3*sm4) ! p3
!
   call rfun( r24,d24 ,k24 )
!
   aa = r34/r24 - r23
!
   if (aa.eq.C0P0) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop trif2: ' &
       ,'threshold singularity, returning 0'
     rslt = C0P0
     return
   endif
!
   bb = -d24/sm3 + r34/sm2 - r23/sm4
   cc = (sm4/sm2 - r24)/(sm3*sm4)
!   hh = dreal(r23 - r24*r34)
!   dd = mysqrt( bb*bb - 4*aa*cc , -dreal(aa)*hh )
   call solabc(x1,x2,dd ,aa,bb,cc ,0)
   x1 = -x1
   x2 = -x2
!
   qx1 = qonv(x1 ,1 ) ! x1 SHOULD HAVE im. part
   qx2 = qonv(x2 ,1 ) ! x2 SHOULD HAVE im. part
   q23 = qonv(r23,-1)
   q24 = qonv(r24,-1)
   q34 = qonv(r34,-1)
   qm2 = qonv(sm2,-1)
   qm3 = qonv(sm3,-1)
   qm4 = qonv(sm4,-1)
!
   rslt = C0P0
!
   qy1 = qx1/q24
   qy2 = qx2/q24
!
   rslt(0) = li2c2( qy1*qm2 ,qy2*qm2 )/r24*sm2
!
   if (x2.ne.C0P0) then ! better to put a threshold on cc 
     rslt(0) = rslt(0) + ( logc2( qy1/qy2 )*logc( qy1*qy2/(qm2*qm2) ) &
                          -logc2( qx1/qx2 )*logc( qx1*qx2/(qm4*qm4) ) )/(x2*2)
   endif
!
   rslt(0) = rslt(0) - li2c2( qx1*qm4 ,qx2*qm4 )*sm4
!
   if (r23.ne.C0P0) then
     qss = q23*qm3/q24
     rslt(0) = rslt(0) - li2c2( qx1*qss ,qx2*qss )*r23*sm3/r24
   endif
!
   if (r34.ne.C0P0) then
     qss = q34*qm3
     rslt(0) = rslt(0) + li2c2( qx1*qss ,qx2*qss )*r34*sm3
   endif
!
   rslt(0) = rslt(0)/(aa*sm2*sm3*sm4)
   end subroutine


   subroutine trif3( rslt ,p1i,p2i,p3i ,m1i,m2i,m3i )
!*******************************************************************
! Finite 1-loop scalar 3-point function with all internal masses
! non-zero. Obtained from the formulas for 4-point functions in
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
! by sending one internal mass to infinity.
!*******************************************************************
   use avh_olo_logc ,only: logc
   use avh_olo_logc2 ,only: logc2
   use avh_olo_li2c2 ,only: li2c2
   complex(kindc2) ,intent(out) :: rslt(0:2)
   complex(kindc2) ,intent(in)  :: p1i,p2i,p3i,m1i,m2i,m3i
   type(qmplx_type) :: q12,q13,q23,qm1,qm2,qm3,qx1,qx2,qz1,qz2,qtt
   complex(kindc2) :: p1,p2,p3,m1,m2,m3,sm1,sm2,sm3,aa,bb,cc,dd,x1,x2 &
                     ,k12,k13,k23,r12,r13,r23,d12,d13,d23 
   real(kindr2) :: h1,h2,h3
  complex(kindc2) ,parameter :: oieps=C1P0+CiP0*epsilon(R1P0)**2
!
!  write(*,*) 'MESSAGE from OneLOop trif3: you are calling me' !CALLINGME
!
   h1 = -aimag(m1i)
   h2 = -aimag(m2i)
   h3 = -aimag(m3i)
   if (h2.ge.h1.and.h2.ge.h3) then
     p1=p3i ;p2=p1i ;p3=p2i ;m1=m3i ;m2=m1i ;m3=m2i
   else
     p1=p1i ;p2=p2i ;p3=p3i ;m1=m1i ;m2=m2i ;m3=m3i
   endif
!
   sm1 = mysqrt(m1)
   sm2 = mysqrt(m2)
   sm3 = mysqrt(m3)
!
   k12 = C0P0
   k13 = C0P0
   k23 = C0P0
   if (m1+m2.ne.p1) k12 = ( m1+m2-p1*oieps )/(sm1*sm2) ! p1
   if (m1+m3.ne.p3) k13 = ( m1+m3-p3*oieps )/(sm1*sm3) ! p1+p2 => p12
   if (m2+m3.ne.p2) k23 = ( m2+m3-p2*oieps )/(sm2*sm3) ! p2
!
   call rfun( r12,d12 ,k12 )
   call rfun( r13,d13 ,k13 )
   call rfun( r23,d23 ,k23 )
!
   aa = sm2/sm3 - k23 + r13*(k12 - sm2/sm1)
!
   if (aa.eq.C0P0) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop trif3: ' &
       ,'threshold singularity, returning 0'
     rslt = C0P0
     return
   endif
!
   bb = d13/sm2 + k12/sm3 - k23/sm1
   cc = ( sm1/sm3 - C1P0/r13 )/(sm1*sm2)
!   hh = dreal( (r13-sm1/sm3)/(sm1*sm2) )
!   dd = mysqrt( bb*bb - 4*aa*cc , -dreal(aa)*hh )
   call solabc( x1,x2,dd ,aa,bb,cc ,0 )
   x1 = -x1
   x2 = -x2
!
   qx1 = qonv(x1 ,1) ! x1 SHOULD HAVE im. part
   qx2 = qonv(x2 ,1) ! x2 SHOULD HAVE im. part
   q12 = qonv(r12,-1)
   q13 = qonv(r13,-1)
   q23 = qonv(r23,-1)
   qm1 = qonv(sm1,-1)
   qm2 = qonv(sm2,-1)
   qm3 = qonv(sm3,-1)
!
   rslt = C0P0
!
   qz1 = qx1*qm2
   qz2 = qx2*qm2
   rslt(0) = rslt(0) + ( li2c2( qz1*q12 ,qz2*q12 )*r12 &
                        +li2c2( qz1/q12 ,qz2/q12 )/r12 )*sm2
   qtt = q13*qm2
   qz1 = qx1*qtt
   qz2 = qx2*qtt
   rslt(0) = rslt(0) - ( li2c2( qz1*q23 ,qz2*q23 )*r23 &
                        +li2c2( qz1/q23 ,qz2/q23 )/r23 )*r13*sm2
   qz1 = qx1*q13
   qz2 = qx2*q13
   rslt(0) = rslt(0) + li2c2( qz1*qm3 ,qz2*qm3 )*r13*sm3 &
                     - li2c2( qx1*qm1 ,qx2*qm1 )*sm1
   if (x2.ne.C0P0) then
     rslt(0) = rslt(0) + ( logc2( qz1/qz2 )*logc( qz1*qz2/(qm3*qm3) ) &
                          -logc2( qx1/qx2 )*logc( qx1*qx2/(qm1*qm1) ) )/(x2*2)
   endif
!
   rslt(0) = rslt(0)/(aa*sm1*sm2*sm3)
   end subroutine
   

   subroutine trif3HV( rslt ,pp,mm ,ap ,smax )
!*******************************************************************
! Finite 1-loop scalar 3-point function with all internal masses
! non-zero. Based on the fomula of 't Hooft & Veltman
!*******************************************************************
   use avh_olo_logc ,only: logc
   use avh_olo_li2c ,only: li2c
   use avh_olo_logc2 ,only: logc2
   complex(kindc2) ,intent(out) :: rslt(0:2)
   complex(kindc2) ,intent(in)  :: pp(3),mm(3)
   real(kindr2)    ,intent(in)  :: ap(3),smax
   complex(kindc2) :: p1,p2,p3,m1,m2,m3,slam,yy
   complex(kindc2) :: sm1,sm2,sm3,sp1,sp2
   type(qmplx_type) :: qm1,qm2,qm3
   real(kindr2) :: a12,a23,a31,thrs,a1,a2,a3
   real(kindc2) ,parameter :: eps1=epsilon(R1P0),eps2=eps1*eps1
   real(kindc2) ,parameter :: small=eps1*100
   complex(kindc2) ,parameter :: ieps2=CiP0*eps2
!
!  write(*,*) 'MESSAGE from OneLOop trif3HV: you are calling me' !CALLINGME
!
! Order squared momenta, first one smallest
   if     (ap(1).le.ap(2).and.ap(1).le.ap(3)) then
     if (ap(2).le.ap(3)) then
       a1=ap(1) ;a2=ap(2) ;a3=ap(3)
       p1=pp(1) ;p2=pp(2) ;p3=pp(3)
       m1=mm(1) ;m2=mm(2) ;m3=mm(3)
     else
       a1=ap(1) ;a2=ap(3) ;a3=ap(2)
       p1=pp(1) ;p2=pp(3) ;p3=pp(2)
       m1=mm(2) ;m2=mm(1) ;m3=mm(3)
     endif
   elseif (ap(2).le.ap(3).and.ap(2).le.ap(1)) then
     if (ap(3).le.ap(1)) then
       a1=ap(2) ;a2=ap(3) ;a3=ap(1)
       p1=pp(2) ;p2=pp(3) ;p3=pp(1)
       m1=mm(2) ;m2=mm(3) ;m3=mm(1)
     else
       a1=ap(2) ;a2=ap(1) ;a3=ap(3)
       p1=pp(2) ;p2=pp(1) ;p3=pp(3)
       m1=mm(3) ;m2=mm(2) ;m3=mm(1)
     endif
   else
     if (ap(1).le.ap(2)) then
       a1=ap(3) ;a2=ap(1) ;a3=ap(2)
       p1=pp(3) ;p2=pp(1) ;p3=pp(2)
       m1=mm(3) ;m2=mm(1) ;m3=mm(2)
     else
       a1=ap(3) ;a2=ap(2) ;a3=ap(1)
       p1=pp(3) ;p2=pp(2) ;p3=pp(1)
       m1=mm(1) ;m2=mm(3) ;m3=mm(2)
     endif
   endif
!
! Need to cut out negligible squared momenta
   thrs = smax*small
!
! Add infinitesimal imaginary parts to masses
   m1=m1-abs(real(m1))*ieps2 ;m2=m2-abs(real(m2))*ieps2 ;m3=m3-abs(real(m3))*ieps2
!       
   if (a1.gt.thrs) then ! 3 non-zero squared momenta
     a1=sign(R1P0,real(p1)) ;a2=sign(R1P0,real(p2)) ;a3=sign(R1P0,real(p3))
     if     (a1.ne.a2) then ;slam=(p3-p1-p2)**2-4*p1*p2
     elseif (a2.ne.a3) then ;slam=(p1-p2-p3)**2-4*p2*p3
     elseif (a3.ne.a1) then ;slam=(p2-p3-p1)**2-4*p3*p1
     else
       sp1=mysqrt(p1,1) ;sp2=mysqrt(p2,1)
       slam=(p3-(sp1+sp2)**2)*(p3-(sp1-sp2)**2)
     endif
     if (slam.eq.C0P0) then
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop trif3HV: ' &
         ,'threshold singularity, returning 0'
       rslt = C0P0
       return
     endif
     slam = mysqrt( slam ,1 )
     sm1=mysqrt(m1,-1) ;sm2=mysqrt(m2,-1) ;sm3=mysqrt(m3,-1)
     rslt = C0P0
     rslt(0) = s3fun( p1,sm1,sm2 , (m2-m3)+p2    ,p3-p1-p2 ,p2 ,slam ) &
             - s3fun( p3,sm1,sm3 ,-(m1-m2)+p3-p2 ,p2-p1-p3 ,p1 ,slam ) &
             + s3fun( p2,sm2,sm3 ,-(m1-m2)+p3-p2 ,p1+p2-p3 ,p1 ,slam )
     rslt(0) = -rslt(0)/slam
!
   elseif (a2.gt.thrs) then ! 2 non-zero squared momenta
     if (p2.eq.p3) then
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop trif3HV: ' &
         ,'threshold singularity, returning 0'
       rslt = C0P0
       return
     endif
     sm1=mysqrt(m1,-1) ;sm2=mysqrt(m2,-1) ;sm3=mysqrt(m3,-1)
     yy = ( (m1-m2)-p3+p2 )/( p2-p3 )
     rslt = C0P0
     rslt(0) = s3fun( p3,sm1,sm3 ,yy ) &
             - s3fun( p2,sm2,sm3 ,yy )
     rslt(0) = rslt(0)/(p2-p3)
!
   elseif (a3.gt.thrs) then ! 1 non-zero squared momentum
     sm1=mysqrt(m1,-1) ;sm3=mysqrt(m3,-1)
     rslt = C0P0
     yy = -( (m1-m2)-p3 )/p3
     rslt(0) = s3fun( p3,sm1,sm3 ,yy ) &
             - s2fun( m2-m3 ,m3 ,yy )
     rslt(0) = -rslt(0)/p3
!
   else ! all squared momenta zero
     rslt = C0P0
     a12=abs(m1-m2) ;a23=abs(m2-m3) ;a31=abs(m3-m1)
     if     (a12.ge.a23.and.a12.ge.a31) then
       if (a12.eq.R0P0) then ;rslt(0)=-C1P0/(2*m3) ;else
       qm1=qonv(m1) ;qm2=qonv(m2) ;qm3=qonv(m3)
       rslt(0) = ( logc2(qm3/qm1) - logc2(qm3/qm2) )/(m1-m2)
       endif
     elseif (a23.ge.a12.and.a23.ge.a31) then
       if (a23.eq.R0P0) then ;rslt(0)=-C1P0/(2*m1) ;else
       qm1=qonv(m1) ;qm2=qonv(m2) ;qm3=qonv(m3)
       rslt(0) = ( logc2(qm1/qm2) - logc2(qm1/qm3) )/(m2-m3)
       endif
     else
       if (a31.eq.R0P0) then ;rslt(0)=-C1P0/(2*m2) ;else
       qm1=qonv(m1) ;qm2=qonv(m2) ;qm3=qonv(m3)
       rslt(0) = ( logc2(qm2/qm3) - logc2(qm2/qm1) )/(m3-m1)
       endif
     endif
   endif
!
   contains
!
     function s3fun( aa,s1,s2 ,t1,t2,t3,t4 ) result(rslt)
!***************************************************************
! int( ( ln(a*y^2+b*y+c) - ln(a*y0^2+b*y0+c) )/(y-y0) ,y=0..1 )
! with  b=s1^2-s2^2-aa  and  c=s2^2
! and with  y0  in terms of t1,t2,t3,t4 defined at the "present"
! function below.
! t4  should be  sqrt(lambda(aa,t2,t3))
!***************************************************************
     complex(kindc2)         ,intent(in) :: aa,s1,s2,t1
     complex(kindc2),optional,intent(in) :: t2,t3
     complex(kindc2),optional,intent(inout) :: t4
     complex(kindc2) :: rslt ,cc,bb,dd,y0,y1,y2,zz,hh,alpha
     real(kindr2) :: rez,arez,aimz
     type(qmplx_type) :: q1,q2
     real(kindr2) ,parameter :: lower=eps1/100,upper=eps1*100
!
     bb = (s1+s2)*(s1-s2)-aa
     cc = s2*s2
     dd = (aa-(s1+s2)**2)*(aa-(s1-s2)**2)
     dd = sqrt( dd )!+ sign(abs(dd),real(aa))*ieps2 )
     call solabc( y1,y2 ,dd ,aa,bb,cc ,1 )
!
     if (present(t4)) then
       call solabc( alpha,hh ,t4 ,aa,t2,t3 ,1 )
       y0 = -(t1+bb*alpha)/t4
     else
       y0 = t1
     endif
!
     q1 = qonv(y0-y1)
     q2 = qonv(y0-y2)
     rslt = li2c(qonv(-y1)/q1) - li2c(qonv(C1P0-y1)/q1) &
          + li2c(qonv(-y2)/q2) - li2c(qonv(C1P0-y2)/q2)
! Take some care about the imaginary part of  a*y0^2+b*y0+c=a*(y0-y1)*(y0-y2)
     zz = y0*(aa*y0+bb)
     rez=real(zz)  ;arez=abs(rez) ;aimz=abs(aimag(zz))
     if (arez*lower.le.aimz.and.aimz.le.arez*upper) then
! Here, the value of Imz is just numerical noise due to cancellations.
! Realize that |Imz|~eps2 indicates there were no such cancellations,
! so the lower limit is needed in in the if-statement!
       zz = (rez + cc)/aa
     else
       zz = (zz + cc)/aa
     endif
     hh = eta3(-y1,-y2,cc/aa) - eta3(y0-y1,y0-y2,zz)
     if (real(aa).lt.R0P0.and.aimag(zz).lt.R0P0) hh = hh - TWOPI*CiP0 
     if (hh.ne.C0P0) rslt = rslt + hh*logc(qonv((y0-C1P0)/y0,1))
!     write(*,*) 'y0',y0 !DEBUG
!     write(*,*) 'y1',y1 !DEBUG
!     write(*,*) 'y2',y2 !DEBUG
!     write(*,*) 'zz',(y0-y1)*(y0-y2),aimag(cc/aa) !DEBUG
!     write(*,*) 'OLO E',aimag(-y1),aimag(-y1),aimag(cc/aa),(eta3(-y1,-y2,cc/aa).ne.C0P0) !DEBUG
!     write(*,*) 'OLO F',aimag(y0-y1),aimag(y0-y1),aimag(zz),(eta3(y0-y1,y0-y2,zz).ne.C0P0) !DEBUG
!     write(*,*) 'OLO G',real(aa),aimag(zz),(real(aa).lt.R0P0.and.aimag(zz).lt.R0P0) !DEBUG
!
     end function
!
     function s2fun( aa,bb ,y0 ) result(rslt)
!**************************************************
! int( ( ln(a*y+b) - ln(a*y0+b) )/(y-y0) ,y=0..1 )
!**************************************************
     complex(kindc2) ,intent(in) :: aa,bb,y0
     complex(kindc2) :: rslt ,y1,hh
     type(qmplx_type) :: q1
     y1 = -bb/aa
     q1 = qonv(y0-y1)
     rslt = li2c(qonv(-y1,-1)/q1) - li2c(qonv(C1P0-y1,-1)/q1)
! aa may have imaginary part, so  theta(-aa)*theta(-Im(y0-y1))  is not
! sufficient and need the following:
     hh = eta5( aa ,-y1,bb ,y0-y1,aa*(y0-y1) )
     if (hh.ne.C0P0) rslt = rslt + hh*logc(qonv((y0-C1P0)/y0,1))
!     write(*,*) 'OLO s2fun aa',aa !DEBUG
!     write(*,*) 'OLO s2fun bb',bb !DEBUG
!     write(*,*) 'OLO s2fun y0',y0 !DEBUG
!     write(*,*) 'OLO s2fun y1',y1 !DEBUG
     end function
!
   end subroutine


end module
