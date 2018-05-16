      subroutine openstorenlo(filetag)
      implicit none
      character * (*) filetag
      character * 100 fname
      integer iunit
      integer innlounit,outnlounit
      common/cnlounits/innlounit,outnlounit
      call getfilename(filetag,fname)
c      inquire(file=fname,exist=present)
      call newunit(iunit)
      open(unit=iunit,file=fname,access='sequential',
     1     form='unformatted',status='new')
      innlounit=0
      outnlounit=iunit
      end

      subroutine closeloadstorenlo
      implicit none
      integer innlounit,outnlounit
      common/cnlounits/innlounit,outnlounit
      if(innlounit.gt.0) then
         close(innlounit)
         innlounit=0
      endif
      if(outnlounit.gt.0) then
         close(outnlounit)
         outnlounit=0
      endif
      end

      subroutine storenlo(xx)
      implicit none
      include 'hepevt.h'
      real * 8 xx
      integer innlounit,outnlounit
      common/cnlounits/innlounit,outnlounit
      if(outnlounit.gt.0) then
         write(outnlounit) nhep
         write(outnlounit) idhep(1:nhep),isthep(1:nhep),
     1        phep(:,1:nhep),xx
      endif
      end

      subroutine storeendseqnlo
      implicit none
      include 'hepevt.h'
      integer innlounit,outnlounit
      common/cnlounits/innlounit,outnlounit
      if(outnlounit.gt.0) then
c a negative nhep signals the end of a bunch of correlated events
         write(outnlounit) -1
      endif
      end

      subroutine readnlo(xx,iret)
      implicit none
      include 'hepevt.h'
      real * 8 xx
      integer iret
      integer innlounit,outnlounit
      common/cnlounits/innlounit,outnlounit
      if(innlounit.gt.0) then
         read(innlounit) nhep
         if(nhep.eq.-1) then
            iret=1
         else
            read(outnlounit) idhep(1:nhep),isthep(1:nhep),
     1           phep(:,1:nhep),xx
            iret=0
         endif
      endif
      end


      subroutine getfilename(filetag,fname)
      implicit none
      character *(*) filetag, fname
      include 'pwhg_rnd.h'
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer ltag,lname
      ltag=len(filetag)
      lname=len(fname)
 1    if(filetag(ltag:ltag).eq.' ') then
        ltag=ltag-1
        goto 1
      endif
      if(rnd_cwhichseed.eq.'none') then
         fname=pwgprefix(1:lprefix)//filetag(1:ltag)//'.dat'
      else
         fname=pwgprefix(1:lprefix)//filetag(1:ltag)//'-'//
     1        rnd_cwhichseed//'.dat'
      endif
      end
