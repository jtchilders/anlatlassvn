	double complex function mylog(z)
	implicit none
	double complex z
#include "HPL.h"

        if (abs(dimag(z)).lt.EPS) then
	   if (abs(z).lt.EPS) then
	      print*,'divergenza log(0)'
              stop
	   else if (dreal(z).lt.0d0) then
	      mylog = log(-z) -dcmplx(0d0,PI)
	   else
	      mylog = log(z)
	   end if
	else
	   mylog = log(z)
	end if

	end function

***************************************************************************

	double complex function myli2(z)
	implicit none
	double complex z,mylog,H0m1
	double precision theta, cl2, ddilog,x,y
#include "HPL.h"

	external cl2,ddilog,mylog, H0m1

        myli2 = -H0m1(-z)

c$$$	x = dreal(z)
c$$$	y = dimag(z)
c$$$
c$$$	if (abs(y).lt.EPS) then
c$$$	   if (x.gt.1d0) then
c$$$	      myli2 = ddilog(x)-
c$$$     $               dcmplx(0d0,PI)*log(x)
c$$$	   else if (abs(z-1d0).lt.EPS) then
c$$$	      myli2 = PI**2/6d0
c$$$	   else
c$$$	      myli2 = ddilog(x)
c$$$	   end if
c$$$	else
c$$$	   if (abs(z-1d0).lt.EPS) then
c$$$	      theta = atan(y/x)
c$$$	      myli2 = cl2(theta)
c$$$	   else
c$$$
c$$$	      print*,'non lo so fare'
c$$$	      stop
c$$$	   end if
c$$$	end if
	
	end function

***************************************************************************

	double complex function myli3(z)
	implicit none
	double complex z
	double precision theta, cl3, dtrilog
	double complex H00m1
#include "HPL.h"

	external cl3, dtrilog, h00m1

	myli3 = - H00m1(-z)

	write(*,*) 'myli3', myli3

c$$$	if (abs(dimag(z)).lt.EPS) then
c$$$	   if (dreal(z).gt.1d0) then
c$$$	      myli3 = dtrilog(dreal(z)) -
c$$$     -                   dcmplx(0d0,PI)/2d0*log(z)**2
c$$$	   else if (abs(z-1d0).lt.EPS) then
c$$$	      myli3 = 1.2020569031d0
c$$$	   else
c$$$	      myli3 = dtrilog(dreal(z))
c$$$	   end if
c$$$	else
c$$$	   if (abs(z-1d0).lt.EPS) then
c$$$	      theta = atan(dimag(z)/dreal(z))
c$$$	      myli3 = cl3(theta)
c$$$	   else
c$$$	      print*,'non lo so fare'
c$$$	      stop
c$$$	   end if
c$$$	end if
	
	end function

***************************************************************************

	double complex function myli4(z)
	implicit none
	double complex z, cquadrilog
	external cquadrilog
#include "HPL.h"

	if (abs(dimag(z)).lt.EPS) then
	   if (dreal(z).gt.1d0) then
	      myli4 = dreal(cquadrilog(z)) -
     -                   dcmplx(0d0,PI)/6d0*log(z)**3
	   else if (abs(z-1d0).lt.EPS) then
	      myli4 = 1.0823232337d0
	   else
	      myli4 = cquadrilog(z)
	   end if
	else
	   if (abs(z-1d0).lt.EPS) then
	      myli4 = cquadrilog(z)
	   else
	      print*,'non lo so fare'
	      stop
	   end if
	end if
	
	end function

***************************************************************************


      double complex function H1(index,z)
      implicit none

#include "HPL.h"
      integer index
      double complex z
      double complex Hm1,H0,Hp1,Hm4,Hmr

      external Hm1,H0,Hp1


      if (index.eq.-4) H1 = Hm4(z)
      if (index.eq.-2) H1 = Hmr(z)
      if (index.eq.-1) H1 = Hm1(z)
      if (index.eq.0)  H1 = H0(z)
      if (index.eq.1)  H1 = Hp1(z)

      end 


      double complex function H2(index1,index2,z)
      implicit none
#include "HPL.h"
      integer index1, index2
      double complex z
      double complex Hmrm1,Hm4mr,Hmrmr,H0m1,H00

      external Hmrm1, Hm4mr, Hmrmr, H0m1,H00

*      if ((index1.eq.-4).and.(index2.eq.-2)) H2 = Hm4mr(x)
      if ((index1.eq.-2).and.(index2.eq.-1)) H2 = Hmrm1(z)
      if ((index1.eq.-2).and.(index2.eq.-2)) H2 = Hmrmr(z)
      if ((index1.eq.0).and.(index2.eq.-1)) H2 = H0m1(z)
      if ((index1.eq.0).and.(index2.eq.0)) H2 = H00(z)

      end 


      double complex function H3(in1,in2,in3,z)
      implicit none
#include "HPL.h"

      double complex z
      integer in1,in2,in3
      double complex Hm4mrm1,Hmrmrmr,H00m1,H0mrmr,Hmrmrm1,Hmr0m1
	double complex H000
      external Hm4mrm1,Hmrmrmr,H00m1,H0mrmr,Hmrmrm1,Hmr0m1,H000

      if ((in1.eq.-4).and.(in2.eq.-2).and.(in3.eq.-1)) H3=Hm4mrm1(z)
      if ((in1.eq.-2).and.(in2.eq.-2).and.(in3.eq.-2)) H3=Hmrmrmr(z)
      if ((in1.eq.-2).and.(in2.eq.-2).and.(in3.eq.-1)) H3=Hmrmrm1(z)
      if ((in1.eq.-2).and.(in2.eq.0).and.(in3.eq.-1)) H3 = Hmr0m1(z)
      if ((in1.eq.0).and.(in2.eq.-2).and.(in3.eq.-2)) H3 = H0mrmr(z)
      if ((in1.eq.0).and.(in2.eq.0).and.(in3.eq.-1)) H3 = H00m1(z)
      if ((in1.eq.0).and.(in2.eq.0).and.(in3.eq.0)) H3 = H000(z)

      end 








******************************************************************************

      double complex function Hm1(z)
      implicit none
#include "HPL.h"

      double precision x,y
      double complex z

	x = dreal(z)
        y = dimag(z)

      if ((x.eq.-1d0).and.(abs(y).lt.EPS)) then
         print*,'errore Hm1: singolare' 
	 stop
      end if

      if (abs(y).lt.EPS) then
         Hm1 = log(dcmplx(x+1,-EPS))
      else
	 Hm1 = log(z+1)
      end if
      end function


      double complex function H0(z)
      implicit none
#include "HPL.h"
      double precision x,y
      double complex z

	x = dreal(z)
	y = dimag(z)


      if ((x.eq.0d0).and.(abs(y).lt.EPS)) then
         print*,'errore H0: singolare' 
	 stop
      end if
      if (abs(y).lt.EPS) then
	 H0 = log(dcmplx(x,-EPS))
      else
	 H0 = log(z)
      endif

      end function


      double complex function Hp1(z)
      implicit none
#include "HPL.h"
      double precision x,y
      double complex z

	x = dreal(z)
	y = dimag(z)

      if ((x.eq.1d0).and.(abs(y).lt.EPS)) then
         print*,'errore Hp1: singolare' 
	 stop
      end if


      if (abs(y).lt.EPS) then
         Hp1 = log(dcmplx(x-1,-EPS))
      else
	 Hp1 = log(z-1)
      end if
 
      end function


      double complex function Hm4(z)
      implicit none
#include "HPL.h"

      double precision x,y
      double complex z

	x = dreal(z)
	y = dimag(z)



      if ((x.eq.-4d0).and.(abs(y).lt.EPS)) then
         print*,'errore Hm4: singolare' 
	 stop
      end if

	if (abs(y).lt.EPS) y = -EPS

      Hm4 = log(dcmplx(x+4,y))-log(4d0)

      end function



      double complex function Hmr(z)
      implicit none
#include "HPL.h"

      double precision x,y
      double complex myasinh,myacosh,myasin
      external myasinh,myacosh,myasin
      double complex z,tmp

	x = dreal(z)
	y = dimag(z)


      if (x.ge.0d0) then
         Hmr = 2d0*myasinh(sqrt(z/4d0))
      else if (x.lt.0d0.and.x.gt.-4d0) then
         tmp = 2d0*myasin(sqrt(-z/4d0))
c         Hmr = dcmplx(0d0,tmp)
         Hmr = dcmplx(0d0,-1d0)*tmp
      else if (x.eq.-4d0) then
         Hmr = dcmplx(0d0,-PI)
      else
         Hmr = 2d0*myacosh(sqrt(-z/4d0))+
     -         dcmplx(0d0,-PI)

      end if

      end function



**************************************************************************
	double complex function H00(z)
	implicit none
	double complex z, mylog
	external mylog
#include "HPL.h"

	H00 = 0.5d0*mylog(z)**2

	end function

	double complex function H000(z)
	implicit none
	double complex z, mylog
	external mylog
#include "HPL.h"

	H000 = 1d0/6d0*mylog(z)**3

	end function

	double complex function H0000(z)
	implicit none
	double complex z, mylog
	external mylog
#include "HPL.h"

	H0000 = 1d0/24d0*mylog(z)**4

	end function

***************************************************************************

	double complex function H0m10(z)
	implicit none
	double complex z, mylog,myli2,myli3
	external mylog,myli2,myli3
#include "HPL.h"

        H0m10=  -mylog(z)*myli2(-z) +2d0*myli3(-z);

	end function


	double complex function H0m100(z)
	implicit none
	double complex z, mylog,myli2,myli3,myli4
	external mylog,myli2,myli3,myli4
#include "HPL.h"

        H0m100= -1d0/2d0* mylog(z)**2 *myli2(-z) +
     -           2d0* mylog(z)* myli3(-z)-3d0* myli4(-z)

	end function



	double complex function H00m10(z)
	implicit none
	double complex z, mylog,myli2,myli3,myli4
	external mylog,myli2,myli3,myli4
#include "HPL.h"

        H00m10= -mylog(z) *myli3(-z) + 3d0* myli4(-z)

	end function



	double complex function H100(z)
	implicit none
	double complex z, mylog,myli2,myli3
	external mylog,myli2,myli3
#include "HPL.h"

        H100= -1d0/2d0* mylog(1-z)* mylog(z)**2 - 
     -        mylog(z)* myli2(z) + myli3(z)

	end function



	double complex function H010(z)
	implicit none
	double complex z, mylog,myli2,myli3
	external mylog,myli2,myli3
#include "HPL.h"

        H010= mylog(z)* myli2(z)- 2d0* myli3(z)

	end function



	double complex function H1000(z)
	implicit none
	double complex z, mylog,myli2,myli3,myli4
	external mylog,myli2,myli3,myli4
#include "HPL.h"

        H1000= -1d0/6d0* mylog(1-z)* mylog(z)**3-
     -          1d0/2d0* mylog(z)**2* myli2(z)+
     -          mylog(z)* myli3(z) - myli4(z)

	end function



	double complex function H0100(z)
	implicit none
	double complex z, mylog,myli2,myli3, myli4
	external mylog,myli2,myli3,myli4
#include "HPL.h"

        H0100= 1d0/2d0* mylog(z)**2* myli2(z) -
     -         2d0* mylog(z)* myli3(z) +3d0* myli4(z)

	end function



	double complex function H0001(z)
	implicit none
	double complex z, myli4
	external myli4
#include "HPL.h"

	H0001= myli4(z)

	end function



	double complex function H0010(z)
	implicit none
	double complex z, mylog,myli2,myli3, myli4
	external mylog,myli2,myli3,myli4
#include "HPL.h"

        H0010 = mylog(z)* myli3(z) - 3d0* myli4(z)

	end function


	double complex function H0m1m10(z)
	implicit none
	double complex z
#include "HPL.h"

c        H0m1m10 = -PI**4/72d0 + PI**2/6d0*mylog(z)*mylog(1+z) -
c     -            0.5d0*mylog(-z)*mylog(z)*mylog(1+z)**2 -
c     -            mylog(z)*mylog(1+z)*myli2(-z) -
c     -            0.5d0*myli2(-z)**2 -
c     -            mylog(z)*myli3(1+z)

	H0m1m10 = 0d0
	stop

	end function


	double complex function H10m10(z)
	implicit none
	double complex z
#include "HPL.h"

        H10m10 = 0d0
	stop


	end function





************************************************************************

      double complex function Hmrm1(z)
      implicit none
#include "HPL.h"
      double precision x,iy,ix2,iy2
      integer dim,i
      parameter (dim=128)
      double precision y(dim),w(dim)
      double complex z,iz,sum,sum2,ccost

      double complex myasinh,myacosh,myasin
      external myasinh,myacosh,myasin

      x=dreal(z)
      iy=dimag(z)

      call gauleg(0d0,x,y,w,dim)

      sum =dcmplx(0d0,0d0)
      sum2=dcmplx(0d0,0d0)

      if (x.gt.0d0) then

          do i=1,dim
            iz = dcmplx(y(i),-EPS)
            sum=sum+w(i)*log(1+iz)/sqrt(iz*(4+iz))
          end do
          Hmrm1 = sum
          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy,y,w,dim)
            do i=1,dim
	       iz = dcmplx(x,y(i))
               sum2=sum2+w(i)*log(1+iz)/sqrt(iz*(4+iz))*dcmplx(0d0,1d0)
	    end do
          Hmrm1 = sum*flagreal+sum2*flagcmplx
          endif

      else if (x.eq.0d0) then
         Hmrm1 = 0d0

      else if (x.gt.-1d0.and.x.lt.0d0) then
            ix2 = dreal(myasin(sqrt(-z/4d0)))
            iy2 = dimag(myasin(sqrt(-z/4d0)))

c          call gauleg(0d0,asin(sqrt(-x/4d0)),y,w,dim)
          call gauleg(0d0,ix2,y,w,dim)
          do i=1,dim
            iz = dcmplx(y(i),-EPS)
            sum=sum+w(i)*log(1-4d0*sin(iz)*sin(iz))*dcmplx(0d0,-2d0)
          end do
          Hmrm1 = sum

          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy2,y,w,dim)
            do i=1,dim
	     iz = dcmplx(ix2,y(i))
             sum2=sum2+w(i)*log(1-4d0*sin(iz)*sin(iz))*dcmplx(0d0,-2d0)*
     -                 dcmplx(0d0,1d0)
	    end do
          Hmrm1 = Hmrm1*flagreal+sum2*flagcmplx
          endif

      else if (x.eq.-1d0) then
         Hmrm1 = dcmplx(0d0,0.67662773760643575d0)

      else if (x.gt.-4d0.and.x.lt.-1d0) then
          ccost = dcmplx(0d0,0.67662773760643575d0)
            ix2 = dreal(myasin(sqrt(-z/4d0)))
            iy2 = dimag(myasin(sqrt(-z/4d0)))

          call gauleg(PI/6d0,ix2,y,w,dim)
          do i=1,dim
            iz = dcmplx(y(i),-EPS)
            sum=sum+w(i)*log(4d0*sin(iz)*sin(iz)-1d0)*dcmplx(0d0,-2d0)
          end do

          Hmrm1 = ccost+sum-
     -      2d0*PI*(ix2-PI/6d0)

          if (abs(iy).gt.EPS) then

            call gauleg(0d0,iy2,y,w,dim)
            do i=1,dim
	       iz = dcmplx(ix2,y(i))
               sum2=sum2+w(i)*log(4d0*sin(iz)*sin(iz)-1d0)*
     -                        dcmplx(0d0,-2d0)

	    end do
	    sum2 = sum2+(-2d0*PI)*iy2
	    sum2 = sum2*dcmplx(0d0,1d0)

          Hmrm1 = Hmrm1*flagreal+sum2*flagcmplx
          endif

      else if (x.eq.-4d0) then
         Hmrm1 = -2d0/3d0*PI*PI

      else if (x.lt.-4d0) then
            ix2=dreal(myacosh(sqrt(-z/4d0)))
            iy2=dimag(myacosh(sqrt(-z/4d0)))

          Hmrm1 = -2d0/3d0*PI*PI
          call gauleg(0d0,ix2,y,w,dim)
          do i=1,dim
            iz = dcmplx(y(i),-EPS)
            sum=sum+w(i)*
     -      2d0*(log(4d0*cosh(y(i))*cosh(y(i))-1d0)-
     -           dcmplx(0d0,PI))
          end do
          Hmrm1 = Hmrm1+sum
c	  print*,z,'  sum=',sum
          if (abs(iy).gt.EPS) then

            call gauleg(0d0,iy2,y,w,dim)
            do i=1,dim
	       iz = dcmplx(ix2,y(i))
               sum2=sum2+w(i)*
     -              2d0*(log((cdexp(iz)+cdexp(-iz))**2-1d0)-
     -              dcmplx(0d0,PI))*dcmplx(0d0,1d0)
	    end do
          Hmrm1 = Hmrm1*flagreal+sum2*flagcmplx
          endif
c	  Hmrm1=sum

      end if


      end function



      double complex function Hmrmr(z)
      implicit none
#include "HPL.h"
      double precision x,iy
      integer dim,i
      parameter (dim=32)
      double precision y(dim),w(dim)
      double complex z,iz,sum,sum2
      double complex myasinh,myacosh,myasin
      external myasinh,myacosh,myasin
      double precision phaseover2factor
      external phaseover2factor

      x=dreal(z)
      iy=dimag(z)


      sum=dcmplx(0d0,0d0)
      sum2=dcmplx(0d0,0d0)

      if (x.gt.0d0) then
	 call gauleg(0d0,x,y,w,dim)
	 do i=1,dim
	    iz = dcmplx(y(i),-EPS)
	    sum=sum+w(i)*2d0*myasinh(sqrt(iz/4d0))/sqrt(iz*(4+iz))
	 end do
	 Hmrmr = sum

          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy,y,w,dim)
            do i=1,dim
            iz = dcmplx(x,y(i))
            sum2=sum2+w(i)*2d0*myasinh(sqrt(iz/4d0))/sqrt(iz*(4+iz))*
     -          dcmplx(0d0,1d0)
	    end do
          Hmrmr = sum*flagreal+sum2*flagcmplx
          endif



      else if ((x.eq.0d0).and.(abs(iy).lt.EPS)) then 
         Hmrmr = dcmplx(0d0,0d0)

      else if (x.eq.0d0) then
         
          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy,y,w,dim)
            do i=1,dim
            iz = dcmplx(x,y(i))
            sum2=sum2+w(i)*myasinh(sqrt(iz/4d0))/sqrt(iz*(4+iz))*
     -       dcmplx(0d0,1d0)
	    end do
          Hmrmr = sum2*flagcmplx
          endif

      else if (x.lt.0d0.and.x.gt.-4) then

         Hmrmr = -2d0*asin(sqrt(-x)/2d0)*asin(sqrt(-x)/2d0)

          if (abs(iy).gt.EPS) then
c$$$            call gauleg(0d0,iy,y,w,dim)
c$$$            do i=1,dim
c$$$            iz = dcmplx(x,y(i))
c$$$            sum2=sum2+
c$$$     -           w(i)*(-2d0)*myasin(sqrt(-iz/4d0))**2*
c$$$     -           dcmplx(0d0,1d0)
c$$$
c$$$	    end do
c$$$          Hmrmr = Hmrmr*flagreal+sum2*flagcmplx
         Hmrmr = -2d0*myasin(sqrt(-z)/2d0)*myasin(sqrt(-z)/2d0)

          endif


      else if ((x.eq.-4).and.(abs(iy).lt.EPS)) then 
         Hmrmr = -PI*PI/2d0

      else if ((x.eq.-4).and.(abs(iy).gt.EPS)) then 
c         Hmrmr = -PI*PI/2d0
         Hmrmr = -2d0*myasin(sqrt(-z)/2d0)*myasin(sqrt(-z)/2d0)

c$$$          if (abs(iy).gt.EPS) then
c$$$            call gauleg(0d0,iy,y,w,dim)
c$$$            do i=1,dim
c$$$            iz = dcmplx(x,y(i))
c$$$            sum2=sum2+w(i)*myasin(sqrt(-iz/4d0))**2*
c$$$     -                dcmplx(0d0,1d0)
c$$$	    end do
c$$$          Hmrmr = Hmrmr*flagreal+sum2*flagcmplx
c$$$          endif




      else if (x.lt.-4) then

         Hmrmr = -PI*PI/2d0 -
     -    2d0*PI*
     -    dcmplx(0d0,1d0)* myacosh(sqrt(-z/4)) + 
     -    2d0* myacosh(sqrt(-z/4))*myacosh(sqrt(-z/4))


      end if

      end function



      double complex function H0m1(z)
      implicit none
#include "HPL.h"
	integer i,dim
        parameter (dim=32)
        double precision yi(dim),w(dim)
      double precision x,y,ddilog
      double complex z,cli2,iz,mylog,check,sum
      external cli2,mylog,ddilog

      x = dreal(z)
      y = dimag(z)


      if (abs(y).lt.EPS) then

	 if (x.ge.-1d0) then
	    iz = dcmplx(-x,0d0)
	    H0m1 = -cli2(iz)
	 else 
	    iz = dcmplx(-x,1d-10)
	    H0m1 = -cli2(iz)
	 end if
	

      else 

	 if (x.ge.-1d0) then
	    iz = dcmplx(-x,0d0)
	    H0m1 = -cli2(iz)

	    sum = dcmplx(0d0,0d0)
	    call gauleg(0d0,y,yi,w,dim)
	    do i=1,dim

	       iz = dcmplx(x,yi(i))
	       sum = sum+w(i)*log(iz+1d0)/(iz)
	    end do

	    H0m1 = H0m1+sum*dcmplx(0d0,1d0)

	 else

	    H0m1 = PI**2*5d0/6d0 -
     -                dcmplx(0d0,1d0)*PI*log(z)+
     -                log(-z)*log(-z-1)+
     -                cli2(dcmplx(x+1,0d0))

	    sum = dcmplx(0d0,0d0)
	    call gauleg(0d0,-y,yi,w,dim)
	    do i=1,dim
	       iz = dcmplx(-x,yi(i))
	       sum = sum+w(i)*(-1d0)*log(iz)/(iz-1d0)
	    end do


	    H0m1 = H0m1+sum*dcmplx(0d0,1d0)

	   endif


      endif

      end function H0m1



      double complex function Hm10(z)
      implicit none
#include "HPL.h"
      double precision x,y
      double complex z,cli2,mylog,iz
      external cli2,mylog

	x = dreal(z)
	y = dimag(z)

      if (x.gt.0d0) then
c         iz = dcmplx(x,0d0)
*         print*,'z1=',z
         Hm10 = cli2(-z) + log(z)*log(1+z)
      else if ((x.eq.0d0).and.(abs(y).lt.EPS)) then
         Hm10 = dcmplx(0d0,0d0)
      else if ((x.eq.0d0).and.(abs(y).gt.EPS)) then
         Hm10 = cli2(-z) + log(z)*log(1+z)
      else if (x.gt.-1d0.and.x.lt.0d0) then
c         iz = dcmplx(x,-1d-10)
*         print*,'z2=',z
         Hm10 = cli2(-z)+mylog(z)*mylog(1+z)
      else if ((x.eq.-1d0).and.(y.eq.0)) then
	 print*,'Hm10 singular'
      else 
         Hm10 = cli2(-z)+mylog(z)*mylog(1+z)

      end if


      end function Hm10


      double complex function H10(z)
      implicit none
#include "HPL.h"
      double precision x,y
      double complex z,cli2,iz
      external cli2

	x = dreal(z)
	y = dimag(z)

        if (abs(y).lt.EPS) then
	   iz = dcmplx(x,-EPS)
	else
	   iz =z 
	end if

	H10 = -cli2(1-iz)


      end function H10

*****************************************************************************

      double complex function H00m1(z)
      implicit none
#include "HPL.h"
      double precision x,iy, dtrilog
      integer dim,i
      parameter (dim=64)
      double precision y(dim),w(dim)
      double complex z,sum,H0m1,iz,sum2,cli2
      external dtrilog,H0m1,cli2

	x = dreal(z)
	iy = dimag(z)
	sum = dcmplx(0d0,0d0)
	sum2= dcmplx(0d0,0d0)

      if (x.gt.-1d0) then
          H00m1 = -dtrilog(-x)

         if (abs(iy).gt.EPS) then
           call gauleg(0d0,iy,y,w,dim)
           do i=1,dim
           iz = dcmplx(x,y(i))
           sum2=sum2+w(i)*cli2(-iz)/(-iz)*dcmplx(0d0,1d0)
	   end do
          H00m1 = H00m1*flagreal+sum2*flagcmplx
          endif

      else

          call gauleg(0d0,x,y,w,dim)

          do i=1,dim
            sum=sum+w(i)*H0m1(dcmplx(y(i),0d0))/y(i)
          end do
c	  print*,'z=',z
c          print*,'sum=',sum
c          print*,'cfr=',-dtrilog(-x)
c	  print*,'  '
          H00m1 = dcmplx(-dtrilog(-x),dimag(sum))

          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy,y,w,dim)
            do i=1,dim
            iz = dcmplx(x,y(i))
            sum2=sum2+w(i)*cli2(-iz)/(-iz)*dcmplx(0d0,1d0)
	    end do
          H00m1 = H00m1*flagreal+sum2*flagcmplx
          endif


      end if


      end function H00m1

      double complex function Hmrmrmr(z)
      implicit none
#include "HPL.h"
      double precision x,iy,iy2,ix2
      integer dim,i
      parameter (dim=32)
      double precision y(dim),w(dim),t
      double complex z,sum,Hmrmr,fm4,sum2,iz,ct
      double complex myasinh,myacosh,mycosh,myasin,mysinh
      external myasinh,myacosh,mycosh,myasin,mysinh
      double precision phaseover2factor
      external phaseover2factor

      external Hmrmr

      fm4 = dcmplx(0d0,5.16771d0)

	x = dreal(z)
	iy = dimag(z)
	sum = dcmplx(0d0,0d0)
	sum2= dcmplx(0d0,0d0)


      if (x.gt.0d0) then
            ix2=dreal(myasinh(sqrt(z/4d0)))
            iy2=dimag(myasinh(sqrt(z/4d0)))


          call gauleg(0d0,ix2,y,w,dim)

          do i=1,dim
	     ct = 4d0*mysinh(dcmplx(y(i),0d0))**2
            sum=sum+w(i)*Hmrmr(ct)*2d0
c	     t = dcmplx(y(i),0d0)
c            sum=sum+w(i)*Hmrmr(t)/sqrt(y(i)*(4d0+y(i)))
          end do
          Hmrmrmr = sum

          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy2,y,w,dim)
            do i=1,dim
	     ct = 4d0*mysinh(dcmplx(ix2,y(i)))**2

            sum2=sum2+w(i)*Hmrmr(ct)*2d0*
     -              dcmplx(0d0,1d0)

c            iz = dcmplx(x,y(i))
c            sum2=sum2+w(i)*Hmrmr(iz)/sqrt(iz*(4d0+iz))*
c     -                     dcmplx(0d0,1d0)
	    end do
            Hmrmrmr = Hmrmrmr*flagreal+sum2*flagcmplx
          endif


      else if (x.eq.0d0) then

          Hmrmrmr = dcmplx(0d0,0d0)

          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy,y,w,dim)
            do i=1,dim
            iz = dcmplx(x,y(i))
            sum2=sum2+w(i)*Hmrmr(iz)/sqrt(iz*(4d0+iz))*
     -                dcmplx(0d0,1d0)
	    end do
            Hmrmrmr = Hmrmrmr*flagreal+sum2*flagcmplx
          endif


      else if (x.lt.0d0.and.x.gt.-4) then
            ix2=dreal(myasin(sqrt(-z/4d0)))
            iy2=dimag(myasin(sqrt(-z/4d0)))

          call gauleg(0d0,ix2,y,w,dim)

          do i=1,dim
	     iz=dcmplx(-4d0*sin(y(i))**2,0d0)
            sum=sum+w(i)*Hmrmr(iz)*
     -                   dcmplx(0d0,-2d0)
          end do
          Hmrmrmr = sum

          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy2,y,w,dim)
	    sum2=dcmplx(0d0,0d0)
            do i=1,dim
            iz = dcmplx(ix2,y(i))

            sum2=sum2+w(i)*Hmrmr(-4d0*sin(iz)**2)*
     -                dcmplx(0d0,-2d0)*dcmplx(0d0,1d0)
	    end do

            Hmrmrmr = Hmrmrmr*flagreal+sum2*flagcmplx
          endif

      else if ((x.eq.-4d0).and.(abs(iy).lt.EPS)) then
         Hmrmrmr = fm4

      else if ((x.eq.-4d0).and.(abs(iy).gt.EPS)) then
         Hmrmrmr = fm4
          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy,y,w,dim)
            do i=1,dim
            iz = dcmplx(x,y(i))
            sum2=sum2+w(i)*Hmrmr(iz)/sqrt(iz*(4d0+iz))*
     -                dcmplx(0d0,1d0)
	    end do
            Hmrmrmr = Hmrmrmr*flagreal+sum2*flagcmplx
          endif

      else if (x.lt.-4d0) then
            ix2=dreal(myacosh(sqrt(-z/4d0)))
            iy2=dimag(myacosh(sqrt(-z/4d0)))

          call gauleg(0d0,ix2,y,w,dim)

          do i=1,dim
	     t = -4d0*cosh(y(i))**2
            sum=sum+w(i)*Hmrmr(dcmplx(t,0d0))*2d0
          end do
          Hmrmrmr = sum+fm4

          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy2,y,w,dim)
            do i=1,dim
            iz = dcmplx(ix2,y(i))
            sum2=sum2+w(i)*Hmrmr(-4d0*mycosh(iz)**2)*2d0*
     -                dcmplx(0d0,1d0)
	    end do
            Hmrmrmr = Hmrmrmr*flagreal+sum2*flagcmplx
          endif


      end if

      end function Hmrmrmr


      double complex function Hmrmrm1(z)
      implicit none
#include "HPL.h"
      double precision x,iy,ix2,iy2
      integer dim,i
      parameter (dim=64)
      double precision y(dim),w(dim)
      double complex z,sum,Hmrm1,fm4,iz,sum2,t
      double complex myasinh,myacosh,mysinh,mycosh,myasin
      external myasinh,myacosh,mysinh,mycosh,myasin
      double precision phaseover2factor
      external phaseover2factor

      external Hmrm1

      fm4 = dcmplx(1.86987d0, 6.89028d0);

	x  = dreal(z)
	iy = dimag(z)
	sum = dcmplx(0d0,0d0)
	sum2= dcmplx(0d0,0d0)

      if (x.gt.0d0) then
            ix2=dreal(myasinh(sqrt(z/4d0)))
            iy2=dimag(myasinh(sqrt(z/4d0)))

          call gauleg(0d0,ix2,y,w,dim)
          sum = dcmplx(0d0,0d0)
          do i=1,dim
            sum=sum+w(i)*Hmrm1(4d0*mysinh(dcmplx(y(i),0d0))**2)*2d0
          end do
          Hmrmrm1 = sum

          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy2,y,w,dim)
            do i=1,dim
            iz = dcmplx(ix2,y(i))
            sum2=sum2+w(i)*Hmrm1(4d0*mysinh(iz)**2)*2d0*
     -                     dcmplx(0d0,1d0)
	    end do
            Hmrmrm1 = Hmrmrm1*flagreal+sum2*flagcmplx
          endif


      else if ((x.eq.0d0).and.(abs(iy).lt.EPS)) then
          Hmrmrm1 = dcmplx(0d0,0d0)

      else if ((x.eq.0d0).and.(abs(iy).gt.EPS)) then
          Hmrmrm1 = dcmplx(0d0,0d0)
            ix2=dreal(myasinh(sqrt(z/4d0)))
            iy2=dimag(myasinh(sqrt(z/4d0)))
          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy2,y,w,dim)
            do i=1,dim
            iz = dcmplx(ix2,y(i))
            sum2=sum2+w(i)*Hmrm1(4d0*mysinh(iz)**2)*2d0*
     -                     dcmplx(0d0,1d0)
	    end do
            Hmrmrm1 = Hmrmrm1*flagreal+sum2*flagcmplx
          endif


      else if (x.lt.0d0.and.x.gt.-4) then
            ix2=dreal(myasin(sqrt(-z/4d0)))
            iy2=dimag(myasin(sqrt(-z/4d0)))

          call gauleg(0d0,ix2,y,w,dim)

          do i=1,dim
            iz = -4d0*sin(dcmplx(y(i),0d0))**2
            sum=sum+w(i)*Hmrm1(iz)*dcmplx(0d0,-2d0)
          end do
          Hmrmrm1 = sum
          if (abs(iy).gt.EPS) then
	     sum2=dcmplx(0d0,0d0)
            call gauleg(0d0,iy2,y,w,dim)
            do i=1,dim
            iz = dcmplx(ix2,y(i))
            sum2=sum2+w(i)*Hmrm1(-4d0*sin(iz)**2)*
     -                dcmplx(0d0,-2d0)*
     -                dcmplx(0d0,1d0)
	    end do
            Hmrmrm1 = Hmrmrm1*flagreal+sum2*flagcmplx
          endif


      else if ((x.eq.-4d0).and.(abs(iy).lt.EPS)) then
         Hmrmrm1 = fm4

      else if ((x.eq.-4d0).and.(abs(iy).gt.EPS)) then
         Hmrmrm1 = fm4
            ix2=dreal(myasin(sqrt(-z/4d0)))
            iy2=dimag(myasin(sqrt(-z/4d0)))

          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy2,y,w,dim)
            do i=1,dim
            iz = dcmplx(ix2,y(i))
            sum2=sum2+w(i)*Hmrm1(-4d0*sin(iz)**2)*
     -                dcmplx(0d0,-2d0)*
     -                dcmplx(0d0,1d0)

	    end do
            Hmrmrm1 = Hmrmrm1*flagreal+sum2*flagcmplx
          endif


      else if (x.lt.-4d0) then
            ix2=dreal(myacosh(sqrt(-z/4d0)))
            iy2=dimag(myacosh(sqrt(-z/4d0)))

          call gauleg(0d0,ix2,y,w,dim)
          sum = dcmplx(0d0,0d0)
          do i=1,dim
	     t = dcmplx(y(i),0d0)
            sum=sum+w(i)*Hmrm1(-4d0*mycosh(t)**2)*2d0
          end do
          Hmrmrm1 = sum+fm4


          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy2,y,w,dim)
            do i=1,dim
            iz = dcmplx(ix2,y(i))
            sum2=sum2+w(i)*Hmrm1(-4d0*mycosh(iz)**2)*2d0*
     -                dcmplx(0d0,1d0)

	    end do
            Hmrmrm1 = Hmrmrm1*flagreal+sum2*flagcmplx
          endif

      end if



      end function Hmrmrm1


      double complex function Hm4mrm1(z)
      implicit none
#include "HPL.h"
      double precision x,iy,ix2,iy2
      integer dim,i
      parameter (dim=2048)
      double precision y(dim),w(dim),t
      double complex z,sum,Hmrm1,ieps,sum2,iz,sum3,ccost
      double complex myasinh,myacosh,mycosh,myasin
      external myasinh,myacosh,mycosh,myasin

c        ieps = dcmplx(0d0,1d-9)
        ieps = dimag(z)
	x = dreal(z)
	iy = dimag(z)
	sum = dcmplx(0d0,0d0)
	sum2= dcmplx(0d0,0d0)
	sum3= dcmplx(0d0,0d0)

      if (x.gt.-4d0) then
c$$$          call gauleg(0d0,x,y,w,dim)
c$$$
c$$$          do i=1,dim
c$$$            sum=sum+w(i)*Hmrm1(dcmplx(y(i),0d0))/(4d0+y(i))
c$$$          end do
c$$$          Hm4mrm1 = sum
c$$$
c$$$          if (abs(iy).gt.EPS) then
c$$$            call gauleg(0d0,iy,y,w,dim)
c$$$            do i=1,dim
c$$$            iz = dcmplx(x,y(i))
c$$$            sum3=sum3+w(i)*Hmrm1(iz)/(4d0+iz)*
c$$$     -                     dcmplx(0d0,1d0)
c$$$	    end do
c$$$          Hm4mrm1 = Hm4mrm1*flagreal+sum3*flagcmplx
c$$$          endif
	 
          ccost = (log(z+4d0)-log(4d0))*Hmrm1(z)

	  if (x.gt.-1d0) then
	     ix2=dreal(myasin(sqrt(-z/4d0)))
	     iy2=dimag(myasin(sqrt(-z/4d0)))
	     call gauleg(0d0,ix2,y,w,dim)
	     do i=1,dim
		sum=sum+w(i)*log(1d0-4d0*sin(y(i))**2)*
     -		  (log(cos(y(i))**2))
	     end do
	     sum = sum*dcmplx(0d0,-2d0)

	     if (abs(iy).gt.EPS) then
		call gauleg(0d0,iy2,y,w,dim)
		do i=1,dim
		   iz = dcmplx(ix2,y(i))
		   sum2=sum2+w(i)*(log(1d0-4d0*sin(iz)**2))*
     -		  (log(cos(iz)**2))
		end do
		sum2=sum2*dcmplx(0d0,1d0)*dcmplx(0d0,-2d0)

	     end if


	  else if (x.eq.-1d0) then
	     sum = sum + 0.0665077d0*dcmplx(0d0,-2d0)
	  else

	     sum = sum + 0.0665077d0
	     ix2=dreal(myasin(sqrt(-z/4d0)))
	     iy2=dimag(myasin(sqrt(-z/4d0)))

	     call gauleg(0.523598775598d0,ix2,y,w,dim)
	     do i=1,dim
		sum=sum+
     -            w(i)*(log(4d0*sin(y(i))**2-1d0)-
     -                  PI*dcmplx(0d0,1d0) )*
     -		  (log(cos(y(i))**2))
	     end do

	     sum = sum*dcmplx(0d0,-2d0)

	     if (abs(iy).gt.EPS) then
		call gauleg(0d0,iy2,y,w,dim)
		do i=1,dim
		   iz = dcmplx(ix2,y(i))
		   sum2=sum2+w(i)*(log(4d0*sin(iz)**2-1d0)-
     -                  PI*dcmplx(0d0,1d0) )*
     -		  (log(cos(iz)**2))
		end do
		sum2=sum2*dcmplx(0d0,1d0)*dcmplx(0d0,-2d0)

	     end if


	  end if

c	  print*,ccost,sum,sum2
          Hm4mrm1 = ccost-sum-sum2


      else if ((x.eq.-4d0).and.(abs(iy).lt.EPS)) then
          ccost = (log(-ieps)-log(4d0))*Hmrm1(z)
          call gauleg(0d0,1.57079632d0,y,w,dim)

          do i=1,dim
            sum=sum+w(i)*log(1d0-4d0*sin(y(i))**2-ieps)*
     -                   log(cos(y(i))**2-ieps)
          end do
            sum = sum*dcmplx(0d0,-2d0)
          Hm4mrm1 = ccost-sum

      else if ((x.eq.-4d0).and.(abs(iy).gt.EPS)) then
          ccost = (log(iy)-log(4d0))*Hmrm1(z)
          call gauleg(0d0,1.57079632d0,y,w,dim)

          do i=1,dim
            sum=sum+w(i)*log(1d0-4d0*sin(y(i))**2-ieps)*
     -                   log(cos(y(i))**2-ieps)*
     -                   dcmplx(0d0,-2d0)
          end do
          Hm4mrm1 = ccost-sum



            call gauleg(0d0,iy,y,w,dim)
            do i=1,dim
            iz = dcmplx(x,y(i))
            sum3=sum3+w(i)*log(1d0-4d0*sin(iz)**2)*
     -                     log(cos(iz)**2)*
     -                     dcmplx(0d0,-2d0)*
     -                     dcmplx(0d0,1d0)
	    end do
          Hm4mrm1 = Hm4mrm1*flagreal+sum3*flagcmplx



      else if (x.lt.-4d0) then

c$$$	 if (abs(iy).lt.EPS) then
c$$$	    ieps = dcmplx(0d0,-1d-9)
c$$$	 else
c$$$	    ieps = dcmplx(0d0,iy)
c$$$	 end if
c$$$
c$$$	 Hm4mrm1 = (log(4d0+x+ieps)-log(4d0))*Hmrm1(z)
c$$$
c$$$
c$$$	  ieps = dcmplx(0d0,1d-9)
c$$$
c$$$          call gauleg(0d0,iy,y,w,dim)
c$$$          do i=1,dim
c$$$c             iz = -4d0*sin(dcmplx(0d0,y(i)))**2
c$$$             iz = dcmplx(0d0,y(i))
c$$$            sum=sum+w(i)/sqrt(iz*(4d0+iz))*log(1d0+iz)*
c$$$     -                   (log(4d0+iz)-log(4d0))*dcmplx(0d0,1d0)
c$$$           
c$$$          end do
c$$$
c$$$
c$$$          call gauleg(0d0,x,y,w,dim)
c$$$          do i=1,dim
c$$$c             iz = -4d0*sin(dcmplx(y(i),iy))**2
c$$$             iz = dcmplx(y(i),iy)
c$$$            sum2=sum2+w(i)/sqrt(iz*(4d0+iz))*log(1d0+iz)*
c$$$     -                   (log(4d0+iz)-log(4d0))
c$$$          end do
c$$$	  print*, z,Hm4mrm1 ,sum ,sum2
c$$$
c$$$	  Hm4mrm1 = Hm4mrm1 - sum - sum2
c$$$



c  H(-4,-r,-1) = H(-4) H(-r,-1) - H(-r,-4,-1) - H(-r,-1,-4)
c  H(-r,-4,-1)+H(-r,-1,-4) = \int_0^x dt (1/sqrt(t (t+4)) H(-1) H(-4)

	 if (abs(iy).lt.EPS) then
	    ieps = dcmplx(0d0,-1d-9)
	 else
	    ieps = dcmplx(0d0,iy)
	 end if

          Hm4mrm1 = (log(4d0+x+ieps)-log(4d0))*Hmrm1(z)

	  ieps = dcmplx(0d0,1d-9)

          call gauleg(0d0,PI/2d0,y,w,dim)
          do i=1,dim
             t = -4d0*sin(y(i))**2
            sum=sum+w(i)*log(1d0+t-ieps)*
     -                   (log(4d0+t-ieps)-log(4d0))*dcmplx(0d0,-2d0)
          end do

          ccost = sum

          ix2=dreal(myacosh(sqrt(-z/4d0)))
          iy2=dimag(myacosh(sqrt(-z/4d0)))
          call gauleg(0d0,ix2,y,w,dim)
          sum2 = dcmplx(0d0,0d0)
          do i=1,dim
             t = -4d0*cosh(y(i))**2
            sum2=sum2+w(i)*log(1d0+t-ieps)*
     -                   (log(4d0+t-ieps)-log(4d0))*(2d0)
          end do

          Hm4mrm1 = Hm4mrm1-ccost-sum2

          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy2,y,w,dim)
            do i=1,dim
            iz = dcmplx(ix2,y(i))
            iz = -4d0*mycosh( dcmplx(ix2,y(i)) )**2

            sum3=sum3+w(i)*(-log(1d0+iz-ieps)*
     -                      (log(4d0+iz-ieps)-log(4d0))*(2d0))*
     -                dcmplx(0d0,1d0)

	    end do

          Hm4mrm1 = Hm4mrm1*flagreal+sum3*flagcmplx
          endif



      end  if

      end function Hm4mrm1




      double complex function H0mrmr(z)
      implicit none
#include "HPL.h"
      double precision x,iy
      integer dim,i
      parameter (dim=128)
      double precision y(dim),w(dim)
      double complex z,sum,Hmrmr,sum2,iz
      double complex myasinh,myacosh
      external myasinh,myacosh

	x = dreal(z)
	iy = dimag(z)
	sum = dcmplx(0d0,0d0)
	sum2= dcmplx(0d0,0d0)

      if (x.gt.0d0) then

          call gauleg(0d0,x,y,w,dim)

          do i=1,dim
            sum=sum+w(i)*Hmrmr(dcmplx(y(i),0d0))/y(i)
          end do
          H0mrmr = sum

          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy,y,w,dim)
            do i=1,dim
            iz = dcmplx(x,y(i))
            sum2=sum2+w(i)*Hmrmr(iz)/iz*dcmplx(0d0,1d0)
	    end do
          H0mrmr = H0mrmr*flagreal+sum2*flagcmplx
          endif

      else if ((x.eq.0d0).and.(abs(iy).lt.EPS)) then
         H0mrmr = 0d0

      else if ((x.eq.0d0).and.(abs(iy).gt.EPS)) then
         H0mrmr = 0d0
          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy,y,w,dim)
            do i=1,dim
            iz = dcmplx(x,y(i))
            sum2=sum2+w(i)*Hmrmr(iz)/iz*dcmplx(0d0,1d0)
	    end do
          H0mrmr = H0mrmr*flagreal+sum2*flagcmplx
          endif

      else if (x.lt.0d0) then

          call gauleg(0d0,x,y,w,dim)

          do i=1,dim
            sum=sum+w(i)*Hmrmr(dcmplx(y(i),0d0))/y(i)
          end do
          H0mrmr = sum

          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy,y,w,dim)
            do i=1,dim
            iz = dcmplx(x,y(i))
            sum2=sum2+w(i)*Hmrmr(iz)/iz*dcmplx(0d0,1d0)
	    end do
          H0mrmr = H0mrmr*flagreal+sum2*flagcmplx
          endif

      end if

      end function H0mrmr





      double complex function Hmr0m1(z)
      implicit none
#include "HPL.h"
      double precision x,iy,ix2,iy2
      integer dim,i
      parameter (dim=128)
      double precision y(dim),w(dim)
      double complex z,sum,H0m1,ccost,sum2,iz,iz2,t
      double complex myasinh,myacosh,mycosh,myasin,mysinh
      external myasinh,myacosh,mycosh,myasin,H0m1,mysinh

	x = dreal(z)
	iy = dimag(z)
	sum = dcmplx(0d0,0d0)
	sum2= dcmplx(0d0,0d0)

      if (x.gt.0d0) then
            ix2=dreal(myasinh(sqrt(z/4d0)))
            iy2=dimag(myasinh(sqrt(z/4d0)))

          call gauleg(0d0,ix2,y,w,dim)

          do i=1,dim
            t=4d0*mysinh(dcmplx(y(i),0d0))**2
            sum=sum+w(i)*H0m1(t)*2d0
          end do
          Hmr0m1 = sum
      
          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy2,y,w,dim)
            do i=1,dim
            iz = 4d0*mysinh(dcmplx(ix2,y(i)))**2
            sum2=sum2+w(i)*H0m1(iz)*2d0*dcmplx(0d0,1d0)
	    end do
          Hmr0m1 = Hmr0m1*flagreal+sum2*flagcmplx
          endif

      else if ((x.eq.0d0).and.(abs(iy).lt.EPS)) then
         Hmr0m1 = 0d0

      else if ((x.eq.0d0).and.(abs(iy).gt.EPS)) then
	 ix2=dreal(myasinh(sqrt(z/4d0)))
	 iy2=dimag(myasinh(sqrt(z/4d0)))
         Hmr0m1 = 0d0

	 call gauleg(0d0,iy2,y,w,dim)
	 do i=1,dim
	    iz = 4d0*mysinh(dcmplx(ix2,y(i)))**2
	    sum2=sum2+w(i)*H0m1(iz)*dcmplx(0d0,1d0)
	 end do
	 Hmr0m1 = Hmr0m1*flagreal+sum2*flagcmplx*2d0


      else if (x.lt.0d0.and.x.gt.-4d0) then
            ix2=dreal(myasin(sqrt(-z/4d0)))
            iy2=dimag(myasin(sqrt(-z/4d0)))

          call gauleg(0d0,ix2,y,w,dim)

          do i=1,dim
            t=-4d0*sin(dcmplx(y(i),0d0))**2
            sum=sum+w(i)*H0m1(t)*dcmplx(0d0,-2d0)
          end do
          Hmr0m1 = sum


          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy2,y,w,dim)
            do i=1,dim
            iz = dcmplx(ix2,y(i))
            iz2 = -4d0*sin(iz)**2
            sum2=sum2+w(i)*H0m1(iz2)*dcmplx(0d0,-2d0)*
     -                dcmplx(0d0,1d0)
	    end do
          Hmr0m1 = Hmr0m1*flagreal+sum2*flagcmplx
          endif


       else if (x.le.-4d0) then

          call gauleg(0d0,PI/2d0,y,w,dim)
	  sum = dcmplx(0d0,0d0)
          do i=1,dim
            t = -4d0*sin(dcmplx(y(i),0d0))**2
            sum=sum+w(i)*H0m1(t)*dcmplx(0d0,-2d0)
          end do
          ccost = sum


	  ix2=dreal(myacosh(sqrt(-z/4d0)))
	  iy2=dimag(myacosh(sqrt(-z/4d0)))

          call gauleg(0d0,ix2,y,w,dim)
	  sum = dcmplx(0d0,0d0)
          do i=1,dim
            t=-4d0*mycosh(dcmplx(y(i),0d0))**2
            sum=sum+w(i)*H0m1(t)*2d0
          end do
          Hmr0m1 = ccost+sum


          if (abs(iy).gt.EPS) then
            call gauleg(0d0,iy2,y,w,dim)
            do i=1,dim
            iz = dcmplx(ix2,y(i))
            iz2 = -4d0*mycosh(iz)**2
            sum2=sum2+w(i)*H0m1(iz2)*2d0
	    end do
	    sum2 = (sum2+0d0*ccost*iy2)*dcmplx(0d0,1d0)

	    Hmr0m1 = Hmr0m1*flagreal+sum2*flagcmplx
          endif


      end if

      end function Hmr0m1



**************************************************************************


      subroutine gauleg(x1,x2,x,w,n)
      implicit none
      integer n
      double precision x1,x2,x(n),w(n)
      double precision gEPS
      parameter (gEPS=3.d-14)

      integer i,j,m
      double precision p1,p2,p3,pp,xl,xm,z,z1

      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)


      do i=1,m
         z=cos(3.141592654d0*(i-0.25d0)/(n+0.5d0))
 1       continue
         p1=1.d0
         p2=0.d0
         do j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
         end do
         pp=n*(z*p1-p2)/(z*z-1.d0)
         z1=z
         z=z1-p1/pp
         if (abs(z-z1).gt.gEPS) goto 1
         x(i)=xm-xl*z
         x(n+1-i)=xm+xl*z
         w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
         w(n+1-i)=w(i)
      end do
      return
      end


******************************************************************************

C     REAL DILOGA
c$$$      DOUBLE PRECISION FUNCTION DDILOG(X)
c$$$ 
c$$$      DOUBLE PRECISION X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO
c$$$      DOUBLE PRECISION C(0:18),H,ALFA,B0,B1,B2
c$$$ 
c$$$      DATA ZERO /0.0D0/, ONE /1.0D0/
c$$$      DATA HALF /0.5D0/, MALF /-0.5D0/, MONE /-1.0D0/, MTWO /-2.0D0/
c$$$      DATA PI3 /3.28986 81336 96453D0/, PI6 /1.64493 40668 48226D0/
c$$$ 
c$$$      DATA C( 0) / 0.42996 69356 08137 0D0/
c$$$      DATA C( 1) / 0.40975 98753 30771 1D0/
c$$$      DATA C( 2) /-0.01858 84366 50146 0D0/
c$$$      DATA C( 3) / 0.00145 75108 40622 7D0/
c$$$      DATA C( 4) /-0.00014 30418 44423 4D0/
c$$$      DATA C( 5) / 0.00001 58841 55418 8D0/
c$$$      DATA C( 6) /-0.00000 19078 49593 9D0/
c$$$      DATA C( 7) / 0.00000 02419 51808 5D0/
c$$$      DATA C( 8) /-0.00000 00319 33412 7D0/
c$$$      DATA C( 9) / 0.00000 00043 45450 6D0/
c$$$      DATA C(10) /-0.00000 00006 05784 8D0/
c$$$      DATA C(11) / 0.00000 00000 86121 0D0/
c$$$      DATA C(12) /-0.00000 00000 12443 3D0/
c$$$      DATA C(13) / 0.00000 00000 01822 6D0/
c$$$      DATA C(14) /-0.00000 00000 00270 1D0/
c$$$      DATA C(15) / 0.00000 00000 00040 4D0/
c$$$      DATA C(16) /-0.00000 00000 00006 1D0/
c$$$      DATA C(17) / 0.00000 00000 00000 9D0/
c$$$      DATA C(18) /-0.00000 00000 00000 1D0/
c$$$ 
c$$$      IF(X .EQ. ONE) THEN
c$$$       DDILOG=PI6
c$$$       RETURN
c$$$      ELSE IF(X .EQ. MONE) THEN
c$$$       DDILOG=MALF*PI6
c$$$       RETURN
c$$$      END IF
c$$$      T=-X
c$$$      IF(T .LE. MTWO) THEN
c$$$       Y=MONE/(ONE+T)
c$$$       S=ONE
c$$$       A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)
c$$$      ELSE IF(T .LT. MONE) THEN
c$$$       Y=MONE-T
c$$$       S=MONE
c$$$       A=LOG(-T)
c$$$       A=-PI6+A*(A+LOG(ONE+ONE/T))
c$$$      ELSE IF(T .LE. MALF) THEN
c$$$       Y=(MONE-T)/T
c$$$       S=ONE
c$$$       A=LOG(-T)
c$$$       A=-PI6+A*(MALF*A+LOG(ONE+T))
c$$$      ELSE IF(T .LT. ZERO) THEN
c$$$       Y=-T/(ONE+T)
c$$$       S=MONE
c$$$       A=HALF*LOG(ONE+T)**2
c$$$      ELSE IF(T .LE. ONE) THEN
c$$$       Y=T
c$$$       S=ONE
c$$$       A=ZERO
c$$$      ELSE
c$$$       Y=ONE/T
c$$$       S=MONE
c$$$       A=PI6+HALF*LOG(T)**2
c$$$      END IF
c$$$ 
c$$$      H=Y+Y-ONE
c$$$      ALFA=H+H
c$$$      B1=ZERO
c$$$      B2=ZERO
c$$$      DO 1 I = 18,0,-1
c$$$      B0=C(I)+ALFA*B1-B2
c$$$      B2=B1
c$$$    1 B1=B0
c$$$      DDILOG=-(S*(B0-H*B2)+A)
c$$$      RETURN
c$$$      END
c$$$ 
c$$$
C     COMPLEX DILOG 
      double COMPLEX FUNCTION CLI2(Z)
      double COMPLEX Z,Z1,CLI20,MOD2
      REAL MOD,PI3,PI6
      DATA PI6/1.6449341/,PI3/3.2898681/

      IF(Z.EQ.(1.,0.)) THEN
        CLI2=CMPLX(PI6,0.)
        RETURN
      ENDIF
      LAB=0
      Z1=Z
      CLI20=(0.,0.)
      MOD=SQRT(REAL(Z)**2+AIMAG(Z)**2)
      IF(MOD.GE.1.)  Z1=1./Z
      IF(REAL(Z1).GE.0.5)  THEN
          LAB=1
          Z1=1.-Z1
      ENDIF
      DO 20 ITE=1,100
      MOD2=Z1**ITE/(ITE**2)
      CLI20=MOD2+CLI20
      IF(SQRT(REAL(MOD2)**2+AIMAG(MOD2)**2).LE.10E-07) GO TO 30
20    CONTINUE
30    IF(LAB.EQ.1)  CLI20=-CLI20+PI6-LOG(Z1)*LOG(1.-Z1)
      IF(MOD.LE.1) THEN
        CLI2=CLI20
        RETURN
      ELSE
        CLI2=-PI6-CLI20-CDLOG(-Z)**2/2.
      ENDIF
      RETURN
      END


      FUNCTION DTRILOG(X)
C     *******************
C CALCULATES THE TRILOGARITHM FOR REAL ARGUMENTS X </= 1.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 WGPLG
      DTRILOG=DREAL(WGPLG(2,1,X))
      RETURN
C***********************************************************************
      END

      COMPLEX*16 FUNCTION WGPLG (N,P,X)
 
      INTEGER P,P1,NC(10),INDEX(31)
      DOUBLE PRECISION FCT(0:4),SGN(0:4),U(0:4),S1(4,4),C(4,4)
      DOUBLE PRECISION A(0:30,10)
      DOUBLE PRECISION X,X1,H,ALFA,R,Q,C1,C2,B0,B1,B2,ZERO,HALF
 
      COMPLEX*16 V(0:5),SK,SM
 
      DATA FCT /1.0D0,1.0D0,2.0D0,6.0D0,24.0D0/
      DATA SGN /1.0D0,-1.0D0,1.0D0,-1.0D0,1.0D0/
      DATA ZERO /0.0D0/, HALF /0.5D0/
      DATA C1 /1.33333 33333 333D0/, C2 /0.33333 33333 3333D0/
 
      DATA S1(1,1) /1.64493 40668 482D0/
      DATA S1(1,2) /1.20205 69031 596D0/
      DATA S1(1,3) /1.08232 32337 111D0/
      DATA S1(1,4) /1.03692 77551 434D0/
      DATA S1(2,1) /1.20205 69031 596D0/
      DATA S1(2,2) /2.70580 80842 778D-1/
      DATA S1(2,3) /9.65511 59989 444D-2/
      DATA S1(3,1) /1.08232 32337 111D0/
      DATA S1(3,2) /9.65511 59989 444D-2/
      DATA S1(4,1) /1.03692 77551 434D0/
 
      DATA C(1,1) / 1.64493 40668 482D0/
      DATA C(1,2) / 1.20205 69031 596D0/
      DATA C(1,3) / 1.08232 32337 111D0/
      DATA C(1,4) / 1.03692 77551 434D0/
      DATA C(2,1) / 0.00000 00000 000D0/
      DATA C(2,2) /-1.89406 56589 945D0/
      DATA C(2,3) /-3.01423 21054 407D0/
      DATA C(3,1) / 1.89406 56589 945D0/
      DATA C(3,2) / 3.01423 21054 407D0/
      DATA C(4,1) / 0.00000 00000 000D0/
 
      DATA INDEX /1,2,3,4,6*0,5,6,7,7*0,8,9,8*0,10/
 
      DATA NC /24,26,28,30,22,24,26,19,22,17/
 
      DATA A( 0,1) / .96753 21504 3498D0/
      DATA A( 1,1) / .16607 30329 2785D0/
      DATA A( 2,1) / .02487 93229 2423D0/
      DATA A( 3,1) / .00468 63619 5945D0/
      DATA A( 4,1) / .00100 16274 9616D0/
      DATA A( 5,1) / .00023 20021 9609D0/
      DATA A( 6,1) / .00005 68178 2272D0/
      DATA A( 7,1) / .00001 44963 0056D0/
      DATA A( 8,1) / .00000 38163 2946D0/
      DATA A( 9,1) / .00000 10299 0426D0/
      DATA A(10,1) / .00000 02835 7538D0/
      DATA A(11,1) / .00000 00793 8705D0/
      DATA A(12,1) / .00000 00225 3670D0/
      DATA A(13,1) / .00000 00064 7434D0/
      DATA A(14,1) / .00000 00018 7912D0/
      DATA A(15,1) / .00000 00005 5029D0/
      DATA A(16,1) / .00000 00001 6242D0/
      DATA A(17,1) / .00000 00000 4827D0/
      DATA A(18,1) / .00000 00000 1444D0/
      DATA A(19,1) / .00000 00000 0434D0/
      DATA A(20,1) / .00000 00000 0131D0/
      DATA A(21,1) / .00000 00000 0040D0/
      DATA A(22,1) / .00000 00000 0012D0/
      DATA A(23,1) / .00000 00000 0004D0/
      DATA A(24,1) / .00000 00000 0001D0/
 
      DATA A( 0,2) / .95180 88912 7832D0/
      DATA A( 1,2) / .43131 13184 6532D0/
      DATA A( 2,2) / .10002 25071 4905D0/
      DATA A( 3,2) / .02442 41559 5220D0/
      DATA A( 4,2) / .00622 51246 3724D0/
      DATA A( 5,2) / .00164 07883 1235D0/
      DATA A( 6,2) / .00044 40792 0265D0/
      DATA A( 7,2) / .00012 27749 4168D0/
      DATA A( 8,2) / .00003 45398 1284D0/
      DATA A( 9,2) / .00000 98586 9565D0/
      DATA A(10,2) / .00000 28485 6995D0/
      DATA A(11,2) / .00000 08317 0847D0/
      DATA A(12,2) / .00000 02450 3950D0/
      DATA A(13,2) / .00000 00727 6496D0/
      DATA A(14,2) / .00000 00217 5802D0/
      DATA A(15,2) / .00000 00065 4616D0/
      DATA A(16,2) / .00000 00019 8033D0/
      DATA A(17,2) / .00000 00006 0204D0/
      DATA A(18,2) / .00000 00001 8385D0/
      DATA A(19,2) / .00000 00000 5637D0/
      DATA A(20,2) / .00000 00000 1735D0/
      DATA A(21,2) / .00000 00000 0536D0/
      DATA A(22,2) / .00000 00000 0166D0/
      DATA A(23,2) / .00000 00000 0052D0/
      DATA A(24,2) / .00000 00000 0016D0/
      DATA A(25,2) / .00000 00000 0005D0/
      DATA A(26,2) / .00000 00000 0002D0/
 
      DATA A( 0,3) / .98161 02799 1365D0/
      DATA A( 1,3) / .72926 80632 0726D0/
      DATA A( 2,3) / .22774 71490 9321D0/
      DATA A( 3,3) / .06809 08329 6197D0/
      DATA A( 4,3) / .02013 70118 3064D0/
      DATA A( 5,3) / .00595 47848 0197D0/
      DATA A( 6,3) / .00176 76901 3959D0/
      DATA A( 7,3) / .00052 74821 8502D0/
      DATA A( 8,3) / .00015 82746 1460D0/
      DATA A( 9,3) / .00004 77492 2076D0/
      DATA A(10,3) / .00001 44792 0408D0/
      DATA A(11,3) / .00000 44115 4886D0/
      DATA A(12,3) / .00000 13500 3870D0/
      DATA A(13,3) / .00000 04148 1779D0/
      DATA A(14,3) / .00000 01279 3307D0/
      DATA A(15,3) / .00000 00395 9070D0/
      DATA A(16,3) / .00000 00122 9055D0/
      DATA A(17,3) / .00000 00038 2658D0/
      DATA A(18,3) / .00000 00011 9459D0/
      DATA A(19,3) / .00000 00003 7386D0/
      DATA A(20,3) / .00000 00001 1727D0/
      DATA A(21,3) / .00000 00000 3687D0/
      DATA A(22,3) / .00000 00000 1161D0/
      DATA A(23,3) / .00000 00000 0366D0/
      DATA A(24,3) / .00000 00000 0116D0/
      DATA A(25,3) / .00000 00000 0037D0/
      DATA A(26,3) / .00000 00000 0012D0/
      DATA A(27,3) / .00000 00000 0004D0/
      DATA A(28,3) / .00000 00000 0001D0/
 
      DATA A( 0,4) /1.06405 21184 614 D0/
      DATA A( 1,4) /1.06917 20744 981 D0/
      DATA A( 2,4) / .41527 19325 1768D0/
      DATA A( 3,4) / .14610 33293 6222D0/
      DATA A( 4,4) / .04904 73264 8784D0/
      DATA A( 5,4) / .01606 34086 0396D0/
      DATA A( 6,4) / .00518 88935 0790D0/
      DATA A( 7,4) / .00166 29871 7324D0/
      DATA A( 8,4) / .00053 05827 9969D0/
      DATA A( 9,4) / .00016 88702 9251D0/
      DATA A(10,4) / .00005 36832 8059D0/
      DATA A(11,4) / .00001 70592 3313D0/
      DATA A(12,4) / .00000 54217 4374D0/
      DATA A(13,4) / .00000 17239 4082D0/
      DATA A(14,4) / .00000 05485 3275D0/
      DATA A(15,4) / .00000 01746 7795D0/
      DATA A(16,4) / .00000 00556 7550D0/
      DATA A(17,4) / .00000 00177 6234D0/
      DATA A(18,4) / .00000 00056 7224D0/
      DATA A(19,4) / .00000 00018 1313D0/
      DATA A(20,4) / .00000 00005 8012D0/
      DATA A(21,4) / .00000 00001 8579D0/
      DATA A(22,4) / .00000 00000 5955D0/
      DATA A(23,4) / .00000 00000 1911D0/
      DATA A(24,4) / .00000 00000 0614D0/
      DATA A(25,4) / .00000 00000 0197D0/
      DATA A(26,4) / .00000 00000 0063D0/
      DATA A(27,4) / .00000 00000 0020D0/
      DATA A(28,4) / .00000 00000 0007D0/
      DATA A(29,4) / .00000 00000 0002D0/
      DATA A(30,4) / .00000 00000 0001D0/
 
      DATA A( 0,5) / .97920 86066 9175D0/
      DATA A( 1,5) / .08518 81314 8683D0/
      DATA A( 2,5) / .00855 98522 2013D0/
      DATA A( 3,5) / .00121 17721 4413D0/
      DATA A( 4,5) / .00020 72276 8531D0/
      DATA A( 5,5) / .00003 99695 8691D0/
      DATA A( 6,5) / .00000 83806 4065D0/
      DATA A( 7,5) / .00000 18684 8945D0/
      DATA A( 8,5) / .00000 04366 6087D0/
      DATA A( 9,5) / .00000 01059 1733D0/
      DATA A(10,5) / .00000 00264 7892D0/
      DATA A(11,5) / .00000 00067 8700D0/
      DATA A(12,5) / .00000 00017 7654D0/
      DATA A(13,5) / .00000 00004 7342D0/
      DATA A(14,5) / .00000 00001 2812D0/
      DATA A(15,5) / .00000 00000 3514D0/
      DATA A(16,5) / .00000 00000 0975D0/
      DATA A(17,5) / .00000 00000 0274D0/
      DATA A(18,5) / .00000 00000 0077D0/
      DATA A(19,5) / .00000 00000 0022D0/
      DATA A(20,5) / .00000 00000 0006D0/
      DATA A(21,5) / .00000 00000 0002D0/
      DATA A(22,5) / .00000 00000 0001D0/
 
      DATA A( 0,6) / .95021 85196 3952D0/
      DATA A( 1,6) / .29052 52916 1433D0/
      DATA A( 2,6) / .05081 77406 1716D0/
      DATA A( 3,6) / .00995 54376 7280D0/
      DATA A( 4,6) / .00211 73389 5031D0/
      DATA A( 5,6) / .00047 85947 0550D0/
      DATA A( 6,6) / .00011 33432 1308D0/
      DATA A( 7,6) / .00002 78473 3104D0/
      DATA A( 8,6) / .00000 70478 8108D0/
      DATA A( 9,6) / .00000 18278 8740D0/
      DATA A(10,6) / .00000 04838 7492D0/
      DATA A(11,6) / .00000 01303 3842D0/
      DATA A(12,6) / .00000 00356 3769D0/
      DATA A(13,6) / .00000 00098 7174D0/
      DATA A(14,6) / .00000 00027 6586D0/
      DATA A(15,6) / .00000 00007 8279D0/
      DATA A(16,6) / .00000 00002 2354D0/
      DATA A(17,6) / .00000 00000 6435D0/
      DATA A(18,6) / .00000 00000 1866D0/
      DATA A(19,6) / .00000 00000 0545D0/
      DATA A(20,6) / .00000 00000 0160D0/
      DATA A(21,6) / .00000 00000 0047D0/
      DATA A(22,6) / .00000 00000 0014D0/
      DATA A(23,6) / .00000 00000 0004D0/
      DATA A(24,6) / .00000 00000 0001D0/
 
      DATA A( 0,7) / .95064 03218 6777D0/
      DATA A( 1,7) / .54138 28546 5171D0/
      DATA A( 2,7) / .13649 97959 0321D0/
      DATA A( 3,7) / .03417 94232 8207D0/
      DATA A( 4,7) / .00869 02788 3583D0/
      DATA A( 5,7) / .00225 28408 4155D0/
      DATA A( 6,7) / .00059 51608 9806D0/
      DATA A( 7,7) / .00015 99561 7766D0/
      DATA A( 8,7) / .00004 36521 3096D0/
      DATA A( 9,7) / .00001 20747 4688D0/
      DATA A(10,7) / .00000 33801 8176D0/
      DATA A(11,7) / .00000 09563 2476D0/
      DATA A(12,7) / .00000 02731 3129D0/
      DATA A(13,7) / .00000 00786 6968D0/
      DATA A(14,7) / .00000 00228 3195D0/
      DATA A(15,7) / .00000 00066 7205D0/
      DATA A(16,7) / .00000 00019 6191D0/
      DATA A(17,7) / .00000 00005 8018D0/
      DATA A(18,7) / .00000 00001 7246D0/
      DATA A(19,7) / .00000 00000 5151D0/
      DATA A(20,7) / .00000 00000 1545D0/
      DATA A(21,7) / .00000 00000 0465D0/
      DATA A(22,7) / .00000 00000 0141D0/
      DATA A(23,7) / .00000 00000 0043D0/
      DATA A(24,7) / .00000 00000 0013D0/
      DATA A(25,7) / .00000 00000 0004D0/
      DATA A(26,7) / .00000 00000 0001D0/
 
      DATA A( 0,8) / .98800 01167 2229D0/
      DATA A( 1,8) / .04364 06760 9601D0/
      DATA A( 2,8) / .00295 09117 8278D0/
      DATA A( 3,8) / .00031 47780 9720D0/
      DATA A( 4,8) / .00004 31484 6029D0/
      DATA A( 5,8) / .00000 69381 8230D0/
      DATA A( 6,8) / .00000 12464 0350D0/
      DATA A( 7,8) / .00000 02429 3628D0/
      DATA A( 8,8) / .00000 00504 0827D0/
      DATA A( 9,8) / .00000 00109 9075D0/
      DATA A(10,8) / .00000 00024 9467D0/
      DATA A(11,8) / .00000 00005 8540D0/
      DATA A(12,8) / .00000 00001 4127D0/
      DATA A(13,8) / .00000 00000 3492D0/
      DATA A(14,8) / .00000 00000 0881D0/
      DATA A(15,8) / .00000 00000 0226D0/
      DATA A(16,8) / .00000 00000 0059D0/
      DATA A(17,8) / .00000 00000 0016D0/
      DATA A(18,8) / .00000 00000 0004D0/
      DATA A(19,8) / .00000 00000 0001D0/
 
      DATA A( 0,9) / .95768 50654 6350D0/
      DATA A( 1,9) / .19725 24967 9534D0/
      DATA A( 2,9) / .02603 37031 3918D0/
      DATA A( 3,9) / .00409 38216 8261D0/
      DATA A( 4,9) / .00072 68170 7110D0/
      DATA A( 5,9) / .00014 09187 9261D0/
      DATA A( 6,9) / .00002 92045 8914D0/
      DATA A( 7,9) / .00000 63763 1144D0/
      DATA A( 8,9) / .00000 14516 7850D0/
      DATA A( 9,9) / .00000 03420 5281D0/
      DATA A(10,9) / .00000 00829 4302D0/
      DATA A(11,9) / .00000 00206 0784D0/
      DATA A(12,9) / .00000 00052 2823D0/
      DATA A(13,9) / .00000 00013 5066D0/
      DATA A(14,9) / .00000 00003 5451D0/
      DATA A(15,9) / .00000 00000 9436D0/
      DATA A(16,9) / .00000 00000 2543D0/
      DATA A(17,9) / .00000 00000 0693D0/
      DATA A(18,9) / .00000 00000 0191D0/
      DATA A(19,9) / .00000 00000 0053D0/
      DATA A(20,9) / .00000 00000 0015D0/
      DATA A(21,9) / .00000 00000 0004D0/
      DATA A(22,9) / .00000 00000 0001D0/
 
      DATA A( 0,10) / .99343 65167 1347D0/
      DATA A( 1,10) / .02225 77012 6826D0/
      DATA A( 2,10) / .00101 47557 4703D0/
      DATA A( 3,10) / .00008 17515 6250D0/
      DATA A( 4,10) / .00000 89997 3547D0/
      DATA A( 5,10) / .00000 12082 3987D0/
      DATA A( 6,10) / .00000 01861 6913D0/
      DATA A( 7,10) / .00000 00317 4723D0/
      DATA A( 8,10) / .00000 00058 5215D0/
      DATA A( 9,10) / .00000 00011 4739D0/
      DATA A(10,10) / .00000 00002 3652D0/
      DATA A(11,10) / .00000 00000 5082D0/
      DATA A(12,10) / .00000 00000 1131D0/
      DATA A(13,10) / .00000 00000 0259D0/
      DATA A(14,10) / .00000 00000 0061D0/
      DATA A(15,10) / .00000 00000 0015D0/
      DATA A(16,10) / .00000 00000 0004D0/
      DATA A(17,10) / .00000 00000 0001D0/
 
      IF(N .LT. 1 .OR. N .GT. 4 .OR. P .LT. 1 .OR. P .GT. 4 .OR.
     1   N+P .GT. 5) THEN
       WGPLG=ZERO
       PRINT 1000, N,P
       RETURN
      END IF
      IF(X .EQ. SGN(0)) THEN
       WGPLG=S1(N,P)
       RETURN
      END IF
 
      IF(X .GT. FCT(2) .OR. X .LT. SGN(1)) THEN
       X1=SGN(0)/X
       H=C1*X1+C2
       ALFA=H+H
       V(0)=SGN(0)
       V(1)=LOG(DCMPLX(-X,ZERO))
       DO 33 L = 2,N+P
   33  V(L)=V(1)*V(L-1)/L
       SK=ZERO
       DO 34 K = 0,P-1
       P1=P-K
       R=X1**P1/(FCT(P1)*FCT(N-1))
       SM=ZERO
       DO 35 M = 0,K
       N1=N+K-M
       L=INDEX(10*N1+P1-10)
       B1=ZERO
       B2=ZERO
       DO 31 I = NC(L),0,-1
       B0=A(I,L)+ALFA*B1-B2
       B2=B1
   31  B1=B0
       Q=(FCT(N1-1)/FCT(K-M))*(B0-H*B2)*R/P1**N1
   35  SM=SM+V(M)*Q
   34  SK=SK+SGN(K)*SM
       SM=ZERO
       DO 36 M = 0,N-1
   36  SM=SM+V(M)*C(N-M,P)
       WGPLG=SGN(N)*SK+SGN(P)*(SM+V(N+P))
       RETURN
      END IF
 
      IF(X .GT. HALF) THEN
       X1=SGN(0)-X
       H=C1*X1+C2
       ALFA=H+H
       V(0)=SGN(0)
       U(0)=SGN(0)
       V(1)=LOG(DCMPLX(X1,ZERO))
       U(1)=LOG(X)
       DO 23 L = 2,P
   23  V(L)=V(1)*V(L-1)/L
       DO 26 L = 2,N
   26  U(L)=U(1)*U(L-1)/L
       SK=ZERO
       DO 24 K = 0,N-1
       P1=N-K
       R=X1**P1/FCT(P1)
       SM=ZERO
       DO 25 M = 0,P-1
       N1=P-M
       L=INDEX(10*N1+P1-10)
       B1=ZERO
       B2=ZERO
       DO 12 I = NC(L),0,-1
       B0=A(I,L)+ALFA*B1-B2
       B2=B1
   12  B1=B0
       Q=SGN(M)*(B0-H*B2)*R/P1**N1
   25  SM=SM+V(M)*Q
   24  SK=SK+U(K)*(S1(P1,P)-SM)
       WGPLG=SK+SGN(P)*U(N)*V(P)
       RETURN
      END IF
 
      L=INDEX(10*N+P-10)
      H=C1*X+C2
      ALFA=H+H
      B1=ZERO
      B2=ZERO
      DO 11 I = NC(L),0,-1
      B0=A(I,L)+ALFA*B1-B2
      B2=B1
   11 B1=B0
      WGPLG=(B0-H*B2)*X**P/(FCT(P)*P**N)
      RETURN
 1000 FORMAT(/' ***** CERN SUBROUTINE WGPLG ... ILLEGAL VALUES',
     1        '   N = ',I3,'   P = ',I3)
      END
 
 
 
      DOUBLE PRECISION FUNCTION CL2(THETA)
C     ************************************
C CALCULATES CLAUSENS INTEGRAL OF SECOND ORDER FOR REAL ARGUMENTS
C THETA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=30,K=20,DM=1.D-15,PM=1.D-35,VM=1.D35)
      DIMENSION B(N)
      B(1)=.1D1/.6D1
      B(2)=.1D1/.30D2
      B(3)=.1D1/.42D2
      B(4)=.1D1/.30D2
      B(5)=.5D1/.66D2
      B(6)=.691D3/.2730D4
      B(7)=.7D1/.6D1
      B(8)=.3617D4/.510D3
      B(9)=.43867D5/.798D3
      B(10)=.174611D6/.330D3
      B(11)=.854513D6/.138D3
      B(12)=.236364091D9/.2730D4
      B(13)=.8553103D7/.6D1
      B(14)=.23749461029D11/.870D3
      B(15)=.8615841276005D13/.14322D5
      B(16)=.7709321041217D13/.510D3
      B(17)=.2577687858367D13/.6D1
      B(18)=.26315271553053477373D20/.1919190D7
      B(19)=.2929993913841559D16/.6D1
      B(20)=.261082718496449122051D21/.13530D5
      B(21)=.1520097643918070802691D22/.1806D4
      B(22)=.27833269579301024235023D23/.690D3
      B(23)=.596451111593912163277961D24/.282D3
      B(24)=.5609403368997817686249127547D28/.46410D5
      B(25)=.495057205241079648212477525D27/.66D2
      B(26)=.801165718135489957347924991853D30/.1590D4
      B(27)=.29149963634884862421418123812691D32/.798D3
      B(28)=.2479392929313226753685415739663229D34/.870D3
      B(29)=.84483613348880041862046775994036021D35/.354D3
      B(30)=.1215233140483755572040304994079820246041491D43/.56786730D8
*      PI=.4D1*DATAN(.1D1)
      X=DMOD(THETA,.2D1*PI)
      IF (X.GT.PI) X=X-.2D1*PI
      IF (X.LE.-PI) X=X+.2D1*PI
      IF (X.EQ.0.D0) THEN
      CL2=0.D0
      ELSE
      H=.1D1-DLOG(DABS(X))
      XS=X**2
      P=XS
      V=.2D1
      DO 1 I=1,K
      D=B(I)/(DFLOAT(2*I*(2*I+1))*V)*P
      H=H+D
      IF (D.LT.DM.OR.P.LT.PM.OR.V.GT.VM) GOTO 2
      P=P*XS
1     V=V*DFLOAT((2*I+1)*(2*I+2))
2     CL2=X*H
      END IF
      RETURN
C***********************************************************************
      END
 
      DOUBLE PRECISION FUNCTION CL3(THETA)
C     ************************************
C CALCULATES CLAUSENS INTEGRAL OF THIRD ORDER FOR REAL ARGUMENTS THETA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=30,K=20,DM=1.D-35,PM=1.D-35,VM=1.D35)
      DIMENSION B(N)
      B(1)=.1D1/.6D1
      B(2)=.1D1/.30D2
      B(3)=.1D1/.42D2
      B(4)=.1D1/.30D2
      B(5)=.5D1/.66D2
      B(6)=.691D3/.2730D4
      B(7)=.7D1/.6D1
      B(8)=.3617D4/.510D3
      B(9)=.43867D5/.798D3
      B(10)=.174611D6/.330D3
      B(11)=.854513D6/.138D3
      B(12)=.236364091D9/.2730D4
      B(13)=.8553103D7/.6D1
      B(14)=.23749461029D11/.870D3
      B(15)=.8615841276005D13/.14322D5
      B(16)=.7709321041217D13/.510D3
      B(17)=.2577687858367D13/.6D1
      B(18)=.26315271553053477373D20/.1919190D7
      B(19)=.2929993913841559D16/.6D1
      B(20)=.261082718496449122051D21/.13530D5
      B(21)=.1520097643918070802691D22/.1806D4
      B(22)=.27833269579301024235023D23/.690D3
      B(23)=.596451111593912163277961D24/.282D3
      B(24)=.5609403368997817686249127547D28/.46410D5
      B(25)=.495057205241079648212477525D27/.66D2
      B(26)=.801165718135489957347924991853D30/.1590D4
      B(27)=.29149963634884862421418123812691D32/.798D3
      B(28)=.2479392929313226753685415739663229D34/.870D3
      B(29)=.84483613348880041862046775994036021D35/.354D3
      B(30)=.1215233140483755572040304994079820246041491D43/.56786730D8
*      PI=.4D1*DATAN(.1D1)
      ZETA3=DTRILOG(.1D1)
C     ZETA3=.12020569032D1
      X=DMOD(THETA,.2D1*PI)
      IF (X.GT.PI) X=X-.2D1*PI
      IF (X.LE.-PI) X=X+.2D1*PI
      IF (X.EQ.0.D0) THEN
      CL3=ZETA3
      ELSE
      H=.75D0-.5D0*DLOG(DABS(X))
      XS=X**2
      P=XS
      V=.2D1
      DO 1 I=1,K
      D=B(I)/(DFLOAT(2*I*(2*I+1)*(2*I+2))*V)*P
      H=H+D
      IF (D.LT.DM.OR.P.LT.PM.OR.V.GT.VM) GOTO 2
      P=P*XS
1     V=V*DFLOAT((2*I+1)*(2*I+2))
2     CL3=ZETA3-X**2*H
      END IF
      RETURN
C***********************************************************************
      END
 
	 
	double complex function myasin(z)
	implicit none
	double precision x,y
	double complex z
#include "HPL.h"

	x = dreal(z)
	y = dimag(z)
	if (y.eq.0d0) then
	   if (abs(z).lt.1d0) then
	      myasin = log(z+sqrt(1d0-z*z)*dcmplx(0d0,1d0))
	   else
	      print*,'arcsin out of range'
	      stop
	   end if
	else
	myasin = log(z+sqrt(z*z-1d0))

	endif

        myasin = myasin*dcmplx(0d0,1d0)
        myasin = myasin+PI/2d0

	end function myasin

	double complex function myasinh(z)
	implicit none
	double complex z

	myasinh = log(z+cdsqrt(z*z+1))

	end function myasinh


	double complex function myacosh(z)
	implicit none
	double complex z

	myacosh = log(z+sqrt(z*z-1))

	end function myacosh

	double complex function mycosh(z)
	implicit none
	double complex z

	mycosh = 0.5d0*(cdexp(z)+cdexp(-z))

	end function mycosh

	double complex function mysinh(z)
	implicit none
	double complex z

	mysinh = 0.5d0*(cdexp(z)-cdexp(-z))

	end function mysinh

*************************************************************************
	double precision function phaseover2factor(z1,z2)
	implicit none
#include "HPL.h"

	double complex z1,z2
	double precision phi1,phi2,phase

	if (abs(dimag(z1)).gt.EPS) then
	   phi1 = atan(dimag(z1)/dreal(z1))
	else
	   if (dreal(z1).lt.0d0) then
	      phi1 = PI
	   else
	      phi1 = 0d0
	   end if
	end if

	if (abs(dimag(z2)).gt.EPS) then
	   phi2 = atan(dimag(z2)/dreal(z2))
	else
	   if (dreal(z2).lt.0d0) then
	      phi2 = PI
	   else
	      phi2 = 0d0
	   end if
	end if

	phase = phi1+phi2


c$$$	if (phi1+phi2.gt.PI) then
c$$$	   phase = phi1+phi2-2d0*PI
c$$$	else if (((phi1+phi2).lt.PI).and.
c$$$    -           ((phi1+phi2).gt.-PI)) then
c$$$	   phase = phi1+phi2
c$$$	else if (phi1+phi2.lt.-PI) then
c$$$	   phase = phi1+phi2+2d0*PI
c$$$        end if

	phaseover2factor = exp(phase/2d0*dcmplx(0d0,1d0))

	end function

*************************************************************************

	double complex function cquadrilog(x)
	double complex x

	cquadrilog=x+
     &  x**2/16. + x**3/81. + x**4/256. + x**5/625. + x**6/1296. + 
     -  x**7/2401. + x**8/4096. + x**9/6561. + x**10/10000. + 
     -  x**11/14641. + x**12/20736. + x**13/28561. + x**14/38416. + 
     -  x**15/50625. + x**16/65536. + x**17/83521. + x**18/104976. + 
     -  x**19/130321. + x**20/160000. + x**21/194481. + x**22/234256.+ 
     -  x**23/279841. + x**24/331776. + x**25/390625.+ x**26/456976. + 
     -  x**27/531441. + x**28/614656. + x**29/707281.+ x**30/810000. + 
     -  x**31/923521. + x**32/1.048576e6 + x**33/1.185921e6 + 
     -  x**34/1.336336e6 + x**35/1.500625e6 + x**36/1.679616e6 + 
     -  x**37/1.874161e6 + x**38/2.085136e6 + x**39/2.313441e6 + 
     -  x**40/2.56e6


	end function



