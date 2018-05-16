      subroutine comparejk(nproc,legs,wgt,wgtx)
      include 'flags.f'
      double precision wgt(6,6),wgtx(6,6)
      integer legs(6),j,k,nproc
      logical writeout
      writeout=.false.
      do j=1,6
      do k=1,6
      if (wgt(j,k) /= 0d0) then
      if (abs(wgtx(j,k)/wgt(j,k)-1d0) .gt. 1d-7) then
      writeout=.true.
      endif
      endif
      enddo
      enddo
      write(6,*) 'comparejk:',nproc,legs(1),legs(2),legs(5),legs(6)
      if (gflag) then
      if (writeout) then
      write(6,*) legs(1),legs(2)
      write(6,*) '1,1',wgt(1,1),wgtx(1,1)
      write(6,*) '2,2',wgt(2,2),wgtx(2,2)
      write(6,*) '5,5',wgt(5,5),wgtx(5,5)
      write(6,*) '6,6',wgt(6,6),wgtx(6,6)
      write(6,*) '1,2',wgt(1,2),wgtx(1,2)
      write(6,*) '1,5',wgt(1,5),wgtx(1,5)
      write(6,*) '1,6',wgt(1,6),wgtx(1,6)
      write(6,*) '2,5',wgt(2,5),wgtx(2,5)
      write(6,*) '2,6',wgt(2,6),wgtx(2,6)
      write(6,*) '5,6',wgt(5,6),wgtx(5,6)
c      write(6,*) '1',wgtx(1,1),wgtx(1,2),wgtx(1,5),wgtx(1,6)
c      write(6,*) '2',wgt(2,1),wgt(2,2),wgt(2,5),wgt(2,6)
c      write(6,*) '2',wgtx(2,1),wgtx(2,2),wgtx(2,5),wgtx(2,6)
c      write(6,*) '5',wgt(5,1),wgt(5,2),wgt(5,5),wgt(5,6)
c      write(6,*) '5',wgtx(5,1),wgtx(5,2),wgtx(5,5),wgtx(5,6)
c      write(6,*) '6',wgt(6,1),wgt(6,2),wgt(6,5),wgt(6,6)
c      write(6,*) '6',wgtx(6,1),wgtx(6,2),wgtx(6,5),wgtx(6,6)
      pause
      endif
      endif
      end
