      subroutine qqb_wpwp_qqb_str(p,msq,chn,identical)
      use qqqqampl
      use consts_MCFM
      use consts_dp, only: ci
      implicit none
      character(len=3) :: chn 
      logical :: identical 
      integer i,j,k
      double precision msq(-5:5,-5:5,3),p(12,4)
      double precision mqqb(3),mqbq(3),mqqq(3),mqbb(3)
      double precision mtot(3),mtot_bits(3)
      double precision aveqqqq, fac
      double precision sw1, sw2, facprop
      double complex propw1, propw2

c---set msq=0 to initialize

      msq=0d0
      aveqqqq = 1d0/9d0/4d0
      fac = (MCFMgw**8)*(MCFMgsq**2)/4d0
      sw1 =2d0*(p(3,4)*p(4,4)-p(3,1)*p(4,1)-p(3,2)*p(4,2)-p(3,3)*p(4,3))
      sw2 =2d0*(p(5,4)*p(6,4)-p(5,1)*p(6,1)-p(5,2)*p(6,2)-p(5,3)*p(6,3))
      propw1 = sw1/(sw1-MCFMwmass**2+ci*MCFMwwidth*MCFMwmass)
      propw2 = sw2/(sw2-MCFMwmass**2+ci*MCFMwwidth*MCFMwmass)
      facprop = abs(propw1)**2*abs(propw2)**2
      fac = fac*aveqqqq*facprop

      if (chn .eq. 'qqb') then

      call mtotqqqq(p(1:8,:),mtot,mtot_bits,chn)
      mqqb(1) = fac*mtot(1)
      mqqb(2) = fac*mtot(2)
      mqqb(3) = fac*mtot(3)

      msq(2,-3,1) = mqqb(3)  ! u sbar initial state
      msq(4,-1,1) = mqqb(3)  ! c dbar initial state
      
C      msq(2,-1,1) = 2d0*fac*mtot_bits(1)  ! u dbar initial state
C      msq(2,-1,2) = fac*mtot_bits(2)  ! u dbar initial state
C      msq(2,-1,3) = fac*mtot_bits(3)  ! u dbar initial state
      if (identical) then 
         msq(2,-1,1) = fac*mtot_bits(1) ! u dbar initial state
         msq(2,-1,2) = fac*mtot_bits(2) ! u dbar initial state
         msq(2,-1,3) = fac*mtot_bits(3) ! u dbar initial state
      else
         msq(2,-1,1) = fac*mtot_bits(1) ! u dbar initial state
         msq(2,-1,2) = 0d0              ! u dbar initial state
         msq(2,-1,3) = 0d0              ! u dbar initial state
      endif

      msq(4,-3,1) = msq(2,-1,1)  ! c sbar initial state
      msq(4,-3,2) = msq(2,-1,2)  ! c sbar initial state
      msq(4,-3,3) = msq(2,-1,3)  ! c sbar initial state

      elseif (chn .eq. 'qbq') then

      call mtotqqqq(p(1:8,:),mtot,mtot_bits,chn)
      mqbq(1) = fac*mtot(1)
      mqbq(2) = fac*mtot(2)
      mqbq(3) = fac*mtot(3)
      
      msq(-3,2,1) = mqbq(3) ! sbar u initial state
      msq(-1,4,1) = mqbq(3) ! dbar c initial state

C      msq(-1,2,1) =  2d0*fac*mtot_bits(1) ! dbar u initial state
C      msq(-1,2,2) =  fac*mtot_bits(2) ! dbar u initial state
C      msq(-1,2,3) =  fac*mtot_bits(3) ! dbar u initial state
      if (identical) then 
         msq(-1,2,1) =  fac*mtot_bits(1) ! dbar u initial state
         msq(-1,2,2) =  fac*mtot_bits(2) ! dbar u initial state
         msq(-1,2,3) =  fac*mtot_bits(3) ! dbar u initial state
      else
         msq(-1,2,1) =  fac*mtot_bits(1) ! dbar u initial state
         msq(-1,2,2) =  0d0  ! dbar u initial state
         msq(-1,2,3) =  0d0  ! dbar u initial state
      endif

      msq(-3,4,1) =  msq(-1,2,1) ! sbar c intital state
      msq(-3,4,2) =  msq(-1,2,2) ! sbar c intital state
      msq(-3,4,3) =  msq(-1,2,3) ! sbar c intital state

      elseif (chn .eq. 'qqq') then

      call mtotqqqq(p(1:8,:),mtot,mtot_bits,chn)
      mqqq(1) = fac*mtot(1)
      mqqq(2) = fac*mtot(2)
      mqqq(3) = fac*mtot(3)

      msq(2,4,1) = mqqq(2)            ! u c initial state
      msq(4,2,1) = mqqq(2)            ! c u initial state

      msq(2,2,1) = fac*mtot_bits(1)*(1d0/2d0)  ! u u initial state
      msq(2,2,2) = fac*mtot_bits(2)*(1d0/2d0)  ! u u initial state
      msq(2,2,3) = fac*mtot_bits(3)*(1d0/2d0)  ! u u initial state

      msq(4,4,1) = msq(2,2,1)   ! u u initial state
      msq(4,4,2) = msq(2,2,2)   ! u u initial state
      msq(4,4,3) = msq(2,2,3)   ! u u initial state

      elseif (chn .eq. 'qbb') then

      call mtotqqqq(p(1:8,:),mtot,mtot_bits,chn)
      mqbb(1) = fac*mtot(1)
      mqbb(2) = fac*mtot(2)
      mqbb(3) = fac*mtot(3)

      msq(-1,-3,1) = mqbb(2)           ! dbar sbar initial state
      msq(-3,-1,1) = mqbb(2)           ! sbar dbar initial state

      msq(-1,-1,1) = fac*mtot_bits(1)/2d0 ! dbar dbar initial state
      msq(-1,-1,2) = fac*mtot_bits(2)/2d0 ! dbar dbar initial state
      msq(-1,-1,3) = fac*mtot_bits(3)/2d0 ! dbar dbar initial state

      msq(-3,-3,1) = msq(-1,-1,1) ! sbar sbar initial state
      msq(-3,-3,2) = msq(-1,-1,2) ! sbar sbar initial state
      msq(-3,-3,3) = msq(-1,-1,3) ! sbar sbar initial state

      endif
      
      return
      end
