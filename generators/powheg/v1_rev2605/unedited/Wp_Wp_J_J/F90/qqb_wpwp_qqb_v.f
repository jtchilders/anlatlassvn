      subroutine qqb_wpwp_qqb_v(p,msq,chn,polesonly,identical)
      use qqqqampl_v
      use consts_MCFM
      use consts_dp
      use sub_defs_io 
      implicit none
      include 'qcdcouple.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f' 
      character(len=3) :: chn 
      logical, intent(in):: polesonly, identical  
      integer i,j,k
      logical :: justepinv1, justepinv2
      double precision msq(-5:5,-5:5),p(12,4), muu
      double precision mqqb(3,-2:1),mqbq(3,-2:1)
      double precision mqqq(3,-2:1),mqbb(3,-2:1)
      double precision mtot(3,-2:1)
      double precision aveqqqq, fac
      double precision epinv1, epinv1a, epinv2a
      double precision sw1, sw2, facprop
      double complex propw1, propw2

c---set msq=0 to initialize

      msq=0d0

      Mtot = zero
      mqbq = zero
      mqqb = czero

      justepinv1 = log_val_opt('-justepinv1',.false.)
      justepinv2 = log_val_opt('-justepinv2',.false.)
      
      muu = scale
      
      if (justepinv2)then
         epinv1 = 0d0
         epinv1a = epinv 
         epinv2a = epinv2
      elseif (justepinv1) then
         epinv1 = epinv
         epinv1a = 0d0
         epinv2a = 0d0
      else
         epinv1 = epinv
         epinv1a = epinv
         epinv2a = epinv2
      endif

      aveqqqq = 1d0/9d0/4d0
      fac = (MCFMgw**8)*(MCFMgsq**2)/4d0
      sw1 =2d0*(p(3,4)*p(4,4)-p(3,1)*p(4,1)-p(3,2)*p(4,2)-p(3,3)*p(4,3))
      sw2 =2d0*(p(5,4)*p(6,4)-p(5,1)*p(6,1)-p(5,2)*p(6,2)-p(5,3)*p(6,3))
      propw1 = sw1/(sw1-MCFMwmass**2+ci*MCFMwwidth*MCFMwmass)
      propw2 = sw2/(sw2-MCFMwmass**2+ci*MCFMwwidth*MCFMwmass)
      facprop = abs(propw1)**2*abs(propw2)**2
      fac = fac*aveqqqq*facprop
      fac = fac*ason2pi 


      if (chn .eq. 'qqb') then
        
         call mtotqqqq_v(p(1:8,:),muu, chn, polesonly,mtot)
        
         mqqb(1,:) = fac*mtot(1,:)
         mqqb(2,:) = fac*mtot(2,:)
         mqqb(3,:) = fac*mtot(3,:)

      elseif (chn .eq. 'qbq') then
         call mtotqqqq_v(p(1:8,:),muu, chn, polesonly,mtot)
         mqbq(1,:) = fac*mtot(1,:)
         mqbq(2,:) = fac*mtot(2,:)
         mqbq(3,:) = fac*mtot(3,:)

      elseif (chn .eq. 'qqq') then
         call mtotqqqq_v(p(1:8,:),muu, chn, polesonly,mtot)
         mqqq(1,:) = fac*mtot(1,:)
         mqqq(2,:) = fac*mtot(2,:)
         mqqq(3,:) = fac*mtot(3,:)

       
      elseif (chn .eq. 'qbb') then
         call mtotqqqq_v(p(1:8,:),muu, chn, polesonly,mtot)
         mqbb(1,:) = fac*mtot(1,:)
         mqbb(2,:) = fac*mtot(2,:)
         mqbb(3,:) = fac*mtot(3,:)
      endif


      if (justepinv2 .or. justepinv1) then
      
         if (chn == 'qqb' .or. chn == 'all') then

C            msq(2,-1) = mqqb(1,-2)*epinv1a*epinv2a + mqqb(1,-1)*epinv1 +          
C     .           mqqb(2,-2)*epinv1a*epinv2a + mqqb(2,-1)*epinv1 ! u dbar initial state 
            if (identical) then 
               msq(2,-1) = mqqb(1,-2)*epinv1a*epinv2a+mqqb(1,-1)*epinv1 ! u dbar -> ub d   
            else
               msq(2,-1) = mqqb(2,-2)*epinv1a*epinv2a+mqqb(2,-1)*epinv1 ! u dbar -> cb s  
            endif
            msq(2,-3) = mqqb(3,-2)*epinv1a*epinv2a + mqqb(3,-1)*epinv1 ! u sbar initial state
            msq(4,-3) = msq(2,-1) ! c sbar initial state
            msq(4,-1) = msq(2,-3)                                                 ! c dbar initial state
            
         elseif (chn == 'qbq' .or. chn == 'all') then

C            msq(-1,2) = mqbq(1,-2)*epinv1a*epinv2a + mqbq(1,-1)*epinv1 +           
C     .           mqbq(2,-2)*epinv1a*epinv2a + mqbq(2,-1)*epinv1                   ! dbar u initial state

            if (identical) then 
               msq(-1,2) = mqbq(1,-2)*epinv1a*epinv2a+mqbq(1,-1)*epinv1           ! dbar u -> d ub 
            else           
               msq(-1,2) = mqbq(2,-2)*epinv1a*epinv2a+mqbq(2,-1)*epinv1           ! dbar u -> cb s
            endif
            msq(-1,4) = mqbq(3,-2)*epinv1a*epinv2a + mqbq(3,-1)*epinv1            ! dbar c initial state

            msq(-3,4) = msq(-1,2)                                                 ! sbar c  initial state

            msq(-3,2) = msq(-1,4)                                                 ! sbar u initial state

         
         elseif (chn == 'qqq' .or. chn == 'all') then

            msq(2,2) = mqqq(1,-2)*epinv1a*epinv2a + mqqq(1,-1)*epinv1             ! u u initial state 
            msq(2,2) = msq(2,2)/two                                               ! Symmetry factor
            
            msq(2,4) = mqqq(2,-2)*epinv1a*epinv2a + mqqq(2,-1)*epinv1             ! u c initial state

            msq(4,4) = msq(2,2)                                                   ! c c initial state
            
            msq(4,2) = msq(2,4)                                                   ! c u initial state
            
            
         elseif (chn == 'qbb' .or. chn == 'all') then

            msq(-1,-1) = mqbb(1,-2)*epinv1a*epinv2a + mqbb(1,-1)*epinv1           ! dbar dbar initial state 
            msq(-1,-1) = msq(-1,-1)/two                                           ! Symmetry factor
            
            msq(-1,-3) = mqbb(2,-2)*epinv1a*epinv2a + mqbb(2,-1)*epinv1           ! dbar sbar initial state

            msq(-3,-3) = msq(-1,-1)                                               ! sbar sbar initial state
            
            msq(-3,-1) = msq(-1,-3)                                               ! sbar dbar initial state

         endif
         
         

      else            ! if not justepinv1 or justepinv2

            if (chn == 'qqb' .or. chn == 'all') then
               
C               msq(2,-1) = mqqb(1,0)+ mqqb(1,1)+ 
C     .              mqqb(1,-2)*epinv1a*epinv2a + mqqb(1,-1)*epinv1 + 
C     .              mqqb(2,0) + mqqb(2,1) +                                        
C     .              mqqb(2,-2)*epinv1a*epinv2a + mqqb(2,-1)*epinv1                ! u dbar initial state

               if (identical) then 
                  msq(2,-1) = mqqb(1,0)+ mqqb(1,1)+ 
     .                 mqqb(1,-2)*epinv1a*epinv2a + mqqb(1,-1)*epinv1 
               else
                  msq(2,-1) = mqqb(2,0) + mqqb(2,1) +                                        
     .                 mqqb(2,-2)*epinv1a*epinv2a + mqqb(2,-1)*epinv1 ! u dbar initial state
               endif
               msq(2,-3)=  mqqb(3,0) + mqqb(3,1) +                                ! u sabr initial state
     .              mqqb(3,-2)*epinv1a*epinv2a + mqqb(3,-1)*epinv1                 

               msq(4,-3) = msq(2,-1)                                              ! c sbar initial state
               msq(4,-1) = msq(2,-3)                                              ! c dbar initial state
               
            elseif (chn == 'qbq' .or. chn == 'all') then
               
C               msq(-1,2) = mqbq(1,0) + mqbq(1,1) + 
C     .              mqbq(1,-2)*epinv1a*epinv2a + mqbq(1,-1)*epinv1 + 
C     .              mqbq(2,0) + mqbq(2,1) + 
C     .              mqbq(2,-2)*epinv1a*epinv2a + mqbq(2,-1)*epinv1                ! dbar u initial state

               if (identical) then 
                  msq(-1,2) = mqbq(1,0) + mqbq(1,1) + 
     .                 mqbq(1,-2)*epinv1a*epinv2a + mqbq(1,-1)*epinv1
               else
                  msq(-1,2) = mqbq(2,0) + mqbq(2,1) + 
     .                 mqbq(2,-2)*epinv1a*epinv2a + mqbq(2,-1)*epinv1 ! dbar u initial state
               endif
               msq(-1,4)=  mqbq(3,0) + mqbq(3,1) +
     .              mqbq(3,-2)*epinv1a*epinv2a + mqbq(3,-1)*epinv1                ! dbar c initial state
               
               msq(-3,4) = msq(-1,2)                                              ! sbar c initial state

               msq(-3,2) = msq(-1,4)                                              ! sbar u initial state

            elseif (chn == 'qqq' .or. chn == 'all') then
               
               msq(2,2) = mqqq(1,0) + mqqq(1,1) +
     .              mqqq(1,-2)*epinv1a*epinv2a + mqqq(1,-1)*epinv1                ! u u initial state
               msq(2,2) = msq(2,2)/two                                            ! Symmetry factor
               
               msq(2,4)=  mqqq(2,0) + mqqq(2,1) + 
     .              mqqq(2,-2)*epinv1a*epinv2a + mqqq(2,-1)*epinv1                ! u c initial state
               
               msq(4,4) = msq(2,2)                                                ! c c initial state

               msq(4,2) = msq(2,4)                                                ! c u initial state


            elseif (chn == 'qbb' .or. chn == 'all') then
               
               msq(-1,-1) = mqbb(1,0) + mqbb(1,1) + 
     .              mqbb(1,-2)*epinv1a*epinv2a + mqbb(1,-1)*epinv1               ! dbar dbar initial state
               msq(-1,-1) = msq(-1,-1)/two                                       ! Symmetry factor
               
               msq(-1,-3)=  mqbb(2,0) + mqbb(2,1) +
     .              mqbb(2,-2)*epinv1a*epinv2a + mqbb(2,-1)*epinv1               ! dbar sbar initial state

               msq(-3,-3) = msq(-1,-1)                                           ! sbar sbar initial state

               msq(-3,-1) = msq(-1,-3)                                           ! sbar dbar initial state
               
            else
               write(*,*) chn
               stop 'qqb_wpwp_qqb_v: channel not recognised!'
            endif

         endif

      return
      end
