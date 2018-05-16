c this subroutine calls right analysis routines, depending on decay mode
c (fully or semi leptonic)

      subroutine init_hist
      implicit none
      include 'cvecbos.h'
      logical :: first = .true. 
      save first 
      real *8 powheginput 
      external powheginput
      integer fat_jet
      save fat_jet

      if (first) then 
c     decay products of the two vector bosons
         vdecaymodeWp=powheginput('vdecaymodeWp')
         vdecaymodeWm=powheginput('vdecaymodeWm')

         fat_jet = powheginput('#fat_jet')
         
c     identify decay mode (default is fully lept.):
         decmode_lep = .true.
         decmode_slp = .false.
         decmode_slm = .false.      
         if (abs(vdecaymodeWp).lt.11) then ! Wp decays hadronically
            if (abs(vdecaymodeWm).lt.11) then
               stop 'fully hadronic decays are not allowed'
            else 
               decmode_slp = .true.
               decmode_lep = .false.
            endif
         else                   ! Wp decays lept.  
            if (abs(vdecaymodeWm).lt.11) then
               decmode_slm = .true.
               decmode_lep = .false.
               if (fat_jet.eq.1) then
                  print*,'no fat-jet analysis available'
                  print*,'for Wm decaying hadronically'
                  stop
               endif   
            endif 
         endif   
         

         if (decmode_slp) then
            if (fat_jet.eq.1) then 
               write(*,*) 'Calling analysis_slp_fat'
            else
               write(*,*) 'Calling analysis_slp' 
            endif   
         elseif (decmode_slm) then
            write(*,*) 'Calling analysis_slm' 
         else
            write(*,*) 'Calling analysis_lep' 
         endif
         first = .false. 

      endif

      if (decmode_slp) then
            if (fat_jet.eq.1) then 
               call init_hist_slp_fat
            else   
               call init_hist_slp
            endif   
      elseif (decmode_slm) then
         call init_hist_slm
      else
         call init_hist_lep
      endif   
      
      end

      subroutine analysis(dsig0)
      implicit none
      real * 8 dsig0
      include 'cvecbos.h'
      logical :: first = .true. 
      save first 
      real *8 powheginput 
      external powheginput
      integer fat_jet
      save fat_jet

      if (first) then 
c     decay products of the two vector bosons
         vdecaymodeWp=powheginput('vdecaymodeWp')
         vdecaymodeWm=powheginput('vdecaymodeWm')

         fat_jet = powheginput('#fat_jet')
         
c     identify decay mode (default is fully lept.):
         decmode_lep = .true.
         decmode_slp = .false.
         decmode_slm = .false.      
         if (abs(vdecaymodeWp).lt.11) then ! Wp decays hadronically
            if (abs(vdecaymodeWm).lt.11) then
               stop 'fully hadronic decays are not allowed'
            else 
               decmode_slp = .true.
               decmode_lep = .false.
            endif
         else                   ! Wp decays lept.  
            if (abs(vdecaymodeWm).lt.11) then
               decmode_slm = .true.
               decmode_lep = .false.
               if (fat_jet.eq.1) then
                  print*,'no fat-jet analysis available'
                  print*,'for Wm decaying hadronically'
                  stop
               endif   
            endif 
         endif   
         
c we pretend that quarks are not quarks, otherwise POWHEG
c makes them radiate
         if(abs(vdecaymodeWm).le.7)
     .     vdecaymodeWm = sign(1,vdecaymodeWm)*(abs(vdecaymodeWm) + 100)
         if(abs(vdecaymodeWp).le.7)
     .     vdecaymodeWp = sign(1,vdecaymodeWp)*(abs(vdecaymodeWp) + 100)

         if (decmode_slp) then
            if (fat_jet.eq.1) then 
               write(*,*) 'Calling analysis_slp_fat' 
            else
               write(*,*) 'Calling analysis_slp'
            endif   
         elseif (decmode_slm) then
            write(*,*) 'Calling analysis_slm' 
         else
            write(*,*) 'Calling analysis_lep' 
         endif
         first = .false. 
      endif
      
      if (decmode_slp) then
         if (fat_jet.eq.1) then 
            call analysis_slp_fat(dsig0)
         else
            call analysis_slp(dsig0)
         endif   
      elseif (decmode_slm) then
         call analysis_slm(dsig0)
      else
         call analysis_lep(dsig0)
      endif   
      
      end
