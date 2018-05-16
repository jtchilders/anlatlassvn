      function pt2max_regular()
      real * 8 pt2max_regular
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'PhysPars.h'
      integer choice
      parameter (choice=1)
      logical debug
      parameter (debug=.false.) 

      if(choice.eq.1) then
         pt2max_regular=(kn_sreal/4)*(1-kn_y**2)*kn_csi**2
C     Maximum pt as in the ISR case, since we adopted ISR
C     parametrization for phase space of regular contributions.
         if(debug) then
            if((dabs((kn_cmpreal(1,4)**2+kn_cmpreal(2,4)**2)/
     $           pt2max_regular)- 1d0).gt.1d-8) then
               write(*,*) 'ERROR: wrong pt2max_regular'
               call exit(1)
            endif
         endif
      elseif(choice.eq.2) then
         pt2max_regular=ph_topmass
C     hard scale of the process
      else
         write(*,*) 'ERROR: wrong choihce in pt2max_regular'
         call exit(1)
      endif
      if (pt2max_regular.lt.rad_ptsqmin) then
         write(*,*) '****************************************'
         write(*,*) 'WARNING in pt2max_regular'
         write(*,*) 'pt2max_regular < rad_ptsqmin ',
     #        pt2max_regular,' < ',rad_ptsqmin
         write(*,*) (flst_regular(i,rad_realreg),i=1,nlegreal)
         write(*,*) 'To generate this event, use the following seeds'
         call printcurrentrandom
         write(*,*) '****************************************'
      endif
      end
