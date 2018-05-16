      subroutine qqb_wpwp_qqb_g(p,msq,chn,identical)
      use qqqqgampl
      use consts_MCFM
      use consts_dp, only: ci
      implicit none
      character(len=3) :: chn 
      logical :: identical 
      integer i,j,k
      double precision msq(-5:5,-5:5),p(12,4)
      double precision mqqb(3),mqbq(3),mqqq(3),mqbb(3)
      double precision mqgl(3),mglq(3),mqbg(3),mgqb(3)
      double precision mtot(3)
      double precision aveqq, aveqg,fac
      double precision sw1, sw2, facprop
      double complex propw1, propw2

c---set msq=0 to initialize

      msq=0d0

      aveqq = 1d0/9d0/4d0
      aveqg = 1d0/3d0/8d0/4d0
      fac = (MCFMgw**8)*(MCFMgsq**2)/4d0 ! Born factor 
      fac = fac * MCFMgsq * 2d0 ! 2 dues to TA normalization 

      fac = fac/4d0

      

      sw1 =2d0*(p(3,4)*p(4,4)-p(3,1)*p(4,1)-p(3,2)*p(4,2)-p(3,3)*p(4,3))
      sw2 =2d0*(p(5,4)*p(6,4)-p(5,1)*p(6,1)-p(5,2)*p(6,2)-p(5,3)*p(6,3))
      propw1 = sw1/(sw1-MCFMwmass**2+ci*MCFMwwidth*MCFMwmass)
      propw2 = sw2/(sw2-MCFMwmass**2+ci*MCFMwwidth*MCFMwmass)
      facprop = abs(propw1)**2*abs(propw2)**2
      fac = fac*facprop

      call mtotqqqqg(p(1:9,:),mtot,chn)

      if (chn .eq. 'qqb') then
         mqqb(1) = fac*aveqq*mtot(1)
         mqqb(2) = fac*aveqq*mtot(2)
         mqqb(3) = fac*aveqq*mtot(3)
      
C         msq(2,-1) = mqqb(1) + mqqb(2)               ! u dbar initial state
         if (identical) then 
            msq(2,-1) = mqqb(1) ! u dbar initial state
         else
            msq(2,-1) = mqqb(2) ! u dbar initial state
         endif
         msq(2,-3) = mqqb(3)                         ! u sbar initial state
         msq(4,-3) = msq(2,-1)               ! c sbar initial state
         msq(4,-1) = mqqb(3)                         ! c dbar initial state


      elseif (chn .eq. 'qbq') then
         mqbq(1) = fac*aveqq*mtot(1)
         mqbq(2) = fac*aveqq*mtot(2)
         mqbq(3) = fac*aveqq*mtot(3)

C         msq(-1,2) = mqbq(1) + mqbq(2)               ! dbar u initial state
         if (identical) then 
            msq(-1,2) = mqbq(1) ! dbar u initial state
         else
            msq(-1,2) = mqbq(2) ! dbar u initial state
         endif
         msq(-3,2) = mqbq(3)    ! sbar u initial state
         msq(-3,4) = msq(-1,2)  ! sbar c intital state
         msq(-1,4) = mqbq(3)    ! dbar c initial state
         

      elseif (chn .eq. 'qqq') then
         mqqq(1) = fac*aveqq*mtot(1)
         mqqq(2) = fac*aveqq*mtot(2)
         mqqq(3) = fac*aveqq*mtot(3)

         msq(2,2) = mqqq(1)*(1d0/2d0)                ! u u initial state
                                                     ! factor due to idential final states
         msq(2,4) = mqqq(2)                          ! u c initial state
         msq(4,2) = mqqq(2)                          ! c u initial state
         msq(4,4) = mqqq(1)*(1d0/2d0 )               ! c c initial state
                                                     ! factor due to idential final states

      elseif (chn .eq. 'qbb') then
         mqbb(1) = fac*aveqq*mtot(1)
         mqbb(2) = fac*aveqq*mtot(2)
         mqbb(3) = fac*aveqq*mtot(3)
         
         msq(-1,-1) = mqbb(1)*(1d0/2d0)              ! dbar dbar initial state
                                                     ! factor due to idential final states
         msq(-1,-3) = mqbb(2)                        ! dbar sbar initial state
         msq(-3,-1) = mqbb(2)                        ! sbar dbar initial state
         msq(-3,-3) = mqbb(1)*(1d0/2d0)              ! sbar sbar initial state
                                                     ! factor due to idential final states

      elseif (chn .eq. 'qgl') then
         mqgl(1) = fac*aveqg*mtot(1)
         mqgl(2) = fac*aveqg*mtot(2)
         mqgl(3) = fac*aveqg*mtot(3)
         
         msq(1,0) = 0d0                              ! d g initial state
         msq(2,0) = mqgl(1)/2d0 + mqgl(2) ! u g initial state
         if (identical) then 
            msq(2,0) = mqgl(1)/2d0 ! u g initial state
         else
            msq(2,0) = mqgl(2) ! u g initial state
         endif
         msq(3,0) = 0d0                              ! s g initial state
         msq(4,0) = msq(2,0)                         ! c g initial state
         
      elseif (chn .eq. 'glq') then
         mglq(1) = fac*aveqg*mtot(1)
         mglq(2) = fac*aveqg*mtot(2)
         mglq(3) = fac*aveqg*mtot(3)

         msq(0,1) = 0d0                             ! g d initial state
C         msq(0,2) = mglq(1)/2d0 + mglq(2)           ! g u initial state
         if (identical) then 
            msq(0,2) = mglq(1)/2d0 ! g u initial state
         else
            msq(0,2) = mglq(2) ! g u initial state
         endif
         msq(0,3) = msq(0,1)                        ! g s initial state
         msq(0,4) = msq(0,2)                        ! g c initial state
      

      elseif (chn .eq. 'qbg') then
         mqbg(1) = fac*aveqg*mtot(1)
         mqbg(2) = fac*aveqg*mtot(2)
         mqbg(3) = fac*aveqg*mtot(3)

C         msq(-1,0) = mqbg(1)/2d0 + mqbg(2) ! db g initial state
         if (identical) then 
            msq(-1,0) = mqbg(1)/2d0  ! db g initial state
         else
            msq(-1,0) = mqbg(2) ! db g initial state
         endif
         msq(-2,0) = 0d0                            ! ub g initial state
         msq(-3,0) = msq(-1,0)                      ! sb g initial state
         msq(-4,0) = msq(-2,0)                      ! cb g initial state


      elseif (chn .eq. 'gqb') then
         mgqb(1) = fac*aveqg*mtot(1)
         mgqb(2) = fac*aveqg*mtot(2)
         mgqb(3) = fac*aveqg*mtot(3)
          
C     msq(0,-1) = mgqb(1)/2d0 + mgqb(2)           ! g db initial state
         if (identical) then 
            msq(0,-1) = mgqb(1)/2d0 ! g db initial state
         else
            msq(0,-1) = mgqb(2) ! g db initial state
         endif
         msq(0,-2) = 0d0        ! g ub initial state
         msq(0,-3) = msq(0,-1)  ! g sb initial state
         msq(0,-4) = msq(0,-2)  ! g cb initial state
         
      endif
      

      return
      end
