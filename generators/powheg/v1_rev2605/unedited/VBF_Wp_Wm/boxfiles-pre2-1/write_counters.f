
      subroutine write_counters
      implicit none
      include 'pwhg_rnd.h'
      integer iun
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      call newunit(iun)
      if(rnd_cwhichseed.eq.'none') then
         open(unit=iun,file=pwgprefix(1:lprefix)//'counters.dat'
     1     ,status='unknown')
      else
         open(unit=iun,file=pwgprefix(1:lprefix)//'counters'//
     1        rnd_cwhichseed//'.dat',status='unknown')
      endif
      call printcnt(iun)
      close(iun)
      end

