      subroutine minlo_checks
      implicit none
      include 'pwhg_flg.h'
      real * 8 powheginput
      external powheginput
      if(.not.flg_minlo) then
         return
      endif
      if(powheginput("#bornonly").eq.1) then
         flg_minlo_nnll=.false.
      elseif(powheginput("#minlo_nnll").eq.0) then
         flg_minlo_nnll=.false.
      else
         flg_minlo_nnll=.true.
      endif
c This is set to true when computing the real contributions;
c it should start as false
      flg_minlo_real=.false.
c
      if(powheginput("#runningscales").eq.1) then
         write(*,*) 
     1        ' minlo_checks: you cannot use runningscales with minlo'
         write(*,*) 'change it in powheg.input'
         call exit(-1)
      endif
      end
