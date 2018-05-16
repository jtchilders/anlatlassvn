      subroutine compare(j,legs,wgt,wgtx)
      double precision wgt,wgtx
      integer legs(6)
      write(6,*) 'compare',legs(1),legs(2),legs(5),legs(6)
      if (abs(wgtx/wgt-1d0) .gt. 1d-7) then
      write(6,*) 'j,wgt',j,wgt
      write(6,*) 'j,wgt',j,wgtx
      endif
      end
