      subroutine v0spinor(p,i,f)
c-----subroutine for massless v spinor

C     Weyl representation
C   Gamma0= 
c            [ 0  0   1    0 ]
c            [               ]
c            [ 0  0   0    1 ]
c            [               ]
c            [ 1  0   0   0  ]
c            [               ]
c            [ 0  1   0   0  ]

c  Gamma1=
c            [  0    0   0  -1 ]
c            [                 ]
c            [  0    0   -1  0 ]
c            [                 ]
c            [  0   +1  0  0   ]
c            [                 ]
c            [ +1   0   0  0   ]
C  Gamma2=
c            [  0    0   0   +%i]
c            [                  ]
c            [  0    0   -%i   0]
c            [                  ]
c            [  0    -%i  0    0]
c            [                  ]
c            [ + %i  0   0    0 ]
c Gamma3= 
c            [  0   0  -1  0 ]
c            [               ]
c            [  0   0  0   1 ]
c            [               ]
c            [  1  0  0   0  ]
c            [               ]
c            [  0  -1  0   0 ]
c 
      implicit none
      integer i
      double complex p(4),f(4),fc,czip,im
      parameter(czip=(0d0,0d0),im=(0d0,1d0))
c--- E=p(1), px=p(2), py=p(3), pz=p(4)
      fc=sqrt(p(1)+p(4))
      
      if (abs(fc) .gt. 1d-8) then 

        if (i.eq.-1) then 
          f(1)=fc
          f(2)=(p(2)+im*p(3))/fc
          f(3)=czip
          f(4)=czip
        elseif (i.eq.1) then 
          f(1)=czip
          f(2)=czip
          f(3)=(p(2)-im*p(3))/fc
          f(4)=-fc
        endif 

      else

        if (i.eq.-1) then    ! case for p(1)=-p(4) along -ve z axis 
          f(1)=czip
          f(2)=sqrt(2d0*p(1))
          f(3)=czip
          f(4)=czip
        elseif (i.eq.1) then
          f(1)=czip
          f(2)=czip
          f(3)=sqrt(2d0*p(1))
          f(4)=czip
        endif 
        
      endif 

      return
      end

