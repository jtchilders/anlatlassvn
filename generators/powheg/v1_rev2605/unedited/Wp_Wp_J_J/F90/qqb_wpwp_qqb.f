      subroutine qqb_wpwp_qqb(p,msq,chn,identical)
      use qqqqampl
      use consts_MCFM
      use consts_dp, only: ci
      implicit none
      character(len=3) :: chn 
      integer i,j,k
      double precision msq(-5:5,-5:5),p(12,4)
      double precision mqqb(3),mqbq(3),mqqq(3),mqbb(3)
      double precision mtot(3),mtot_bits(3)
      double precision aveqqqq, fac
      double precision sw1, sw2, facprop
      double complex propw1, propw2
      logical :: identical 

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

      
!      msq(2,-1) = mqqb(1) + mqqb(2) ! u dbar initial state
!      msq(4,-3) = mqqb(1) + mqqb(2) ! c sbar initial state
      if (identical) then 
         msq(2,-1) = mqqb(1) ! u dbar initial state
         msq(4,-3) = mqqb(1) ! c sbar initial state
      else
         msq(2,-1) = mqqb(2) ! u dbar initial state
         msq(4,-3) = mqqb(2) ! c sbar initial state
      endif

      msq(2,-3) = mqqb(3)           ! u sbar initial state
      msq(4,-1) = mqqb(3)           ! c dbar initial state

      elseif (chn .eq. 'qbq') then
      call mtotqqqq(p(1:8,:),mtot,mtot_bits,chn)
      mqbq(1) = fac*mtot(1)
      mqbq(2) = fac*mtot(2)
      mqbq(3) = fac*mtot(3)
      
!      msq(-1,2) = mqbq(1) + mqbq(2) ! dbar u initial state
!      msq(-3,4) = mqbq(1) + mqbq(2) ! sbar c intital state
      if (identical) then 
         msq(-1,2) = mqbq(1)    ! dbar u initial state
         msq(-3,4) = mqbq(1)    ! dbar u initial state
      else
         msq(-1,2) = mqbq(2)    ! sbar c intital state
         msq(-3,4) = mqbq(2)    ! sbar c intital state
      endif
      msq(-3,2) = mqbq(3)           ! sbar u initial state
      msq(-1,4) = mqbq(3)           ! dbar c initial state

      elseif (chn .eq. 'qqq') then
      call mtotqqqq(p(1:8,:),mtot,mtot_bits,chn)
      mqqq(1) = fac*mtot(1)
      mqqq(2) = fac*mtot(2)
      mqqq(3) = fac*mtot(3)

      msq(2,2) = mqqq(1)*(1d0/2d0)  ! u u initial state
      msq(2,4) = mqqq(2)            ! u c initial state
      msq(4,2) = mqqq(2)            ! c u initial state
      msq(4,4) = mqqq(1)*(1d0/2d0)  ! c c initial state

      elseif (chn .eq. 'qbb') then
      call mtotqqqq(p(1:8,:),mtot,mtot_bits,chn)
      mqbb(1) = fac*mtot(1)
      mqbb(2) = fac*mtot(2)
      mqbb(3) = fac*mtot(3)

      msq(-1,-1) = mqbb(1)*(1d0/2d0) ! dbar dbar initial state
      msq(-1,-3) = mqbb(2)           ! dbar sbar initial state
      msq(-3,-3) = mqbb(1)*(1d0/2d0) ! sbar sbar initial state
      msq(-3,-1) = mqbb(2)           ! sbar dbar initial state

      endif


      
      return
      end

C     returns, in the case of idential quarks, the two matrix
C     elements squared corersponding to the two possible diagrams 
      subroutine qqb_wpwp_qqb_w1w2(p,chn,w1,w2)
      use qqqqampl
      use consts_MCFM
      use consts_dp, only: ci
      implicit none
      character(len=3) :: chn 
      integer i,j,k
      double precision p(12,4),w1,w2
      double precision mtot(3),mtot_bits(3)
      double precision aveqqqq, fac
      double precision sw1, sw2, facprop
      double complex propw1, propw2
      logical :: identical 

      aveqqqq = 1d0/9d0/4d0
      fac = (MCFMgw**8)*(MCFMgsq**2)/4d0
      sw1 =2d0*(p(3,4)*p(4,4)-p(3,1)*p(4,1)-p(3,2)*p(4,2)-p(3,3)*p(4,3))
      sw2 =2d0*(p(5,4)*p(6,4)-p(5,1)*p(6,1)-p(5,2)*p(6,2)-p(5,3)*p(6,3))
      propw1 = sw1/(sw1-MCFMwmass**2+ci*MCFMwwidth*MCFMwmass)
      propw2 = sw2/(sw2-MCFMwmass**2+ci*MCFMwwidth*MCFMwmass)
      facprop = abs(propw1)**2*abs(propw2)**2
      fac = fac*aveqqqq*facprop

      call mtotqqqq(p(1:8,:),mtot,mtot_bits,chn)
      w1 = fac*mtot_bits(1)
      w2 = fac*mtot_bits(2)

      end      
