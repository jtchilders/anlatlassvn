      subroutine init_processes
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
      include "pwhg_flg.h"
      include "pwhg_st.h"
      include "pwhg_par.h"
      integer i
      real * 8 powheginput
      external powheginput

C     uncomment to link with MINLO 
      st_bornorder=2
      if(powheginput("#minlo").eq.1) then
      flg_minlo=.true.
       else
      flg_minlo=.false.
      endif

      if(powheginput("#withnegweights").eq.0) then
         flg_withnegweights = .false.
      else
         flg_withnegweights = .true.
      endif

      if(powheginput("#withdamp").eq.0) then
         flg_withdamp = .false.
         flg_bornzerodamp=.false.
      else
         flg_withdamp=.true.
         flg_bornzerodamp=.true.
      endif

      if (powheginput("#doublefsr").eq.0) then
         flg_doublefsr=.false.
      else
         flg_doublefsr=.true.
      endif
      if(powheginput("#par_diexp").lt.0) par_diexp=2
      if(powheginput("#par_dijexp").lt.0) par_dijexp=2
      if(powheginput("#par_2gsupp").lt.0) par_2gsupp=4


      call init_processes_born_MCFM
      call init_processes_real_MCFM
      call init_couplings

      st_nlight=5
      do i=3,nlegreal
         if (abs(flst_real(i,1)).le.st_nlight) then
            flst_lightpart=i
            exit
         endif
      enddo
 
      return
      end
 
 

 
 
 
