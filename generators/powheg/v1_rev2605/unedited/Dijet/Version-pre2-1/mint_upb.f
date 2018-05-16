c      implicit none
c      integer ndim
c      character * 20 pwgprefix
c      integer lprefix
c      common/cpwgprefix/pwgprefix,lprefix
c      include 'pwhg_rnd.h'
c      parameter (ndim=19)
c      real * 8 ymax(50,ndim),xint,ymmm,ymin
c      integer k,j
c      lprefix=3
c      pwgprefix='pwg'
c      rnd_cwhichseed='none'
c      xint=8.8d-3
c      call loadmintupb(ndim,'btildeupb',xint,ymax)
c      ymmm=0
c      ymin=1d30
c      do k=1,50
c         do j=1,ndim
c            ymmm=max(ymmm,ymax(k,j))
c            ymin=min(ymin,ymax(k,j))
c         enddo
c         write(*,'(19(d8.2,1x))') (ymax(k,j),j=1,ndim)
c      enddo
c      write(*,*) ymmm, ymin, xint**(1d0/ndim)
c      end


c initialize the storage of values for the determination of the
c upper bounding envelope in MINT
      subroutine startstoremintupb(filetag)
      implicit none
      include 'pwhg_rnd.h'
      character * (*) filetag
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer iunit,l
      logical active
      common/storeubc/iunit,active
      save /storeubc/
      data active/.false./
      active=.true.
      l=len(filetag)
 1    if(filetag(l:l).eq.' ') then
        l=l-1
        goto 1
      endif
      call newunit(iunit)
      if(rnd_cwhichseed.eq.'none') then
         open(unit=iunit,
     1        file=pwgprefix(1:lprefix)//filetag(1:l)//'.dat',
     2        status='unknown')
      else
         open(unit=iunit,
     1        file=pwgprefix(1:lprefix)//filetag(1:l)//'-'//
     1        rnd_cwhichseed//'.dat',status='unknown')
      endif
      end

      subroutine storemintupb(ndim,ncell,imode,f,f0)
      implicit none
      include 'nlegborn.h'
      integer ndim,ncell(ndim),imode
      real * 8 f,f0
      integer iunit
      logical active
      common/storeubc/iunit,active
      character * 30 fmt
      integer k
      logical ini
      data ini/.true./
      save ini,fmt
      if(active) then
         if(ini) then
            fmt='(   i2,2d11.5)'
            write(fmt(2:4),'(i2)') ndim
            ini=.false.
         endif
         if(imode.eq.0) then
            write(iunit,fmt) (ncell(k),k=1,ndim),f,f0
         else
            write(iunit,fmt) (ncell(k),k=1,ndim),f
         endif
      endif
      end

      subroutine stopstoremintupb
      implicit none
      integer iunit
      logical active
      common/storeubc/iunit,active
      close(iunit)
      active=.false.
      end

      subroutine getlinemintupb1(filetag,ndim,cells,f,f0,iret)
      implicit none
      include 'nlegborn.h'
      character *(*) filetag
      integer ndim,cells(ndiminteg),iret
      real * 8 f,f0
      character * 1, dimension(:,:), allocatable, save :: allcells
      real, dimension(:), allocatable, save :: allf
      real, dimension(:), allocatable, save :: allf0
      integer nlines,j
      integer status,index
      data status/0/
      save status,index,nlines
c a call with iret=-10 deallocate all arrays and returns ;
      if(iret.eq.-10) then
         write(*,*) ' deallocating mintupb arrays'
         deallocate(allcells,stat=j)
         deallocate(allf,stat=j)
         deallocate(allf0,stat=j)
         write(*,*) ' end deallocating '
         status=0
         return
      endif
      if(status.eq.0) then
         write(*,*) ' getlinemintupb1: loading file(s)'
c status=0 is the initial call;
         iret=0
         nlines=0
c in this block count the lines in the file (nlines
 1       continue
         call getlinemintupb(filetag,ndim,cells,f,f0,iret)
         if(iret.eq.0) then
            nlines=nlines+1
            goto 1
         endif
c lines counted
c allocates enough stuff for nlines lines
         deallocate(allcells,stat=iret)
         allocate(allcells(ndim,nlines),stat = iret)
         deallocate(allf,stat=iret)
         allocate(allf(nlines),stat=iret)
         deallocate(allf0,stat=iret)
         allocate(allf0(nlines),stat=iret)
c store file content in allocated array
         do j=1,nlines
            call getlinemintupb(filetag,ndim,cells,f,f0,iret)
            call inttochar(cells,allcells(:,j),ndim)
            allf(j)=f
            allf0(j)=f0
         enddo
c this call should return 1 in iret (no more lines)
         call getlinemintupb(filetag,ndim,cells,f,f0,iret)
c Set status to 1 (cells loaded), index to 1 (initiate reading)
         status=1
         index=1
         write(*,*) 'getlinemintupb1: loaded file(s)'
      endif
      if(index.le.nlines) then
         f=allf(index)
         f0=allf0(index)
         call chartoint(allcells(:,index),cells,ndim)
         index=index+1
         iret=0
      else
c end of data
         iret=1
         index=1
      endif
      end
         
      subroutine getlinemintupb(filetag,ndim,cells,f,f0,iret)
c reads a line from the upper bound file or files.
c iret=-1: failure
c iret=0 : success
c iret=1 : end of stream; next call restart reading
c          from the beginning
      implicit none
      character *(*) filetag
      integer ndim,cells(ndim),iret
      real * 8 f,f0
      include 'pwhg_rnd.h'
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer ltag
      character * 100 fname
      integer jfile,k,iunit
      character * 20 fmt
      character * 4 chnum
      logical lpresent
      integer status
      data status/0/
      save status,iunit,jfile,fmt
c initial call
      if(status.eq.0) then
         fmt='(   i2,2d11.5)'
         ltag=len(filetag)
 33      if(filetag(ltag:ltag).eq.' ') then
            ltag=ltag-1
            goto 33
         endif
c The format for reading the bounds file, for ndim dimensions
         write(fmt(2:4),'(i2)') ndim
         call newunit(iunit)
c see if there are files to load
         if(rnd_cwhichseed.eq.'none') then
            fname=pwgprefix(1:lprefix)//filetag(1:ltag)//'.dat'
         else
            do jfile=1,9999
               write(chnum,'(i4)') jfile
               if(chnum(1:1).eq.' ') chnum(1:1)='0'
               if(chnum(2:2).eq.' ') chnum(2:2)='0'
               if(chnum(3:3).eq.' ') chnum(3:3)='0'
               fname=pwgprefix(1:lprefix)//filetag(1:ltag)
     1              //'-'//chnum//'.dat'
               inquire(file=fname,exist=lpresent)
               if(lpresent) goto 10
            enddo
            goto 999
         endif
 10      continue
         open(unit=iunit,file=fname,status='old',err=999)
         write(*,*) ' opened ',trim(fname)
c file opened for reading
         status=1
      endif
 12   continue
      if(filetag.eq.'btildeupb') then
         read(unit=iunit,fmt=fmt,end=11) (cells(k),k=1,ndim),f,f0
      else
         read(unit=iunit,fmt=fmt,end=11) (cells(k),k=1,ndim),f
         f0=-1
      endif
      iret=0
      return
 11   continue
      close(iunit)
      if(rnd_cwhichseed.ne.'none') then
 13      jfile=jfile+1
         if(jfile.lt.9999) then
            write(chnum,'(i4)') jfile
            if(chnum(1:1).eq.' ') chnum(1:1)='0'
            if(chnum(2:2).eq.' ') chnum(2:2)='0'
            if(chnum(3:3).eq.' ') chnum(3:3)='0'
            fname=pwgprefix(1:lprefix)//filetag//'-'//chnum//'.dat'
            inquire(file=fname,exist=lpresent)
            if(lpresent) then
               open(unit=iunit,file=fname,status='old',err=999)
               write(*,*) ' opened ',trim(fname)
               goto 12
            else
               goto 13
            endif
         endif
      endif
      iret=1
      status=0
      return
 999  iret=-1
      end
         





      subroutine loadmintupb(ndim,filetag,ymax,ymaxrat)
      implicit none
      integer ndim
      character *(*) filetag
      real * 8 ymax(50,ndim),xint,xintrat
      real * 8 ymaxrat(50,ndim)
      integer cells(ndim)
      integer kdim,iret,j,iunit
      real * 8 f,f0,prod,prodrat,xless,xmore,xmorerat
      real * 8 fail,tot,ubtot,failrat,totrat,ubtotrat
      integer ipoints
      integer iterations
      logical btildeflg
      iterations=1
      if(filetag.eq.'btildeupb') then
         btildeflg=.true.
      else
         btildeflg=.false.
      endif
c First compute total for f and f/f0 (f/b)
      iret=0
      xint=0
      xintrat=0
      ipoints=0
      do while (iret.eq.0)
         call getlinemintupb1(filetag,ndim,cells,f,f0,iret)
         if(iret.eq.-1) goto 998
         ipoints=ipoints+1
         xint=xint+f
         if(btildeflg) then
            if(f0.gt.0) xintrat=xintrat+f/f0
         endif
      enddo
      ipoints=ipoints-1
      xint=xint/ipoints
      xintrat=xintrat/ipoints/10
      xmore=0.01
      xmorerat=xmore/5
      xless=xmore/2
c The bound is initially set as if the function was uniform
c with the given integral
      do kdim=1,ndim
         do j=1,50
            ymax(j,kdim)=xint**(1d0/ndim)
            if(btildeflg) then
               if(xintrat.gt.0) ymaxrat(j,kdim)=xintrat**(1d0/ndim)
            endif
         enddo
      enddo
c loop for reading u-bound data
 1    continue
      call getlinemintupb1(filetag,ndim,cells,f,f0,iret)
      if(iret.lt.0) then
         write(*,*) ' error while loading bound files'
         call exit(-1)
      endif
      if(iret.eq.0) then
         prod=1
         prodrat=1
         do kdim=1,ndim
            prod=prod*ymax(cells(kdim),kdim)
            if(btildeflg) prodrat=prodrat*ymaxrat(cells(kdim),kdim)
         enddo
         if(f.gt.prod) then
            do kdim=1,ndim
               ymax(cells(kdim),kdim)=
     1              ymax(cells(kdim),kdim)*(f/prod+0.1)**(xmore/ndim)
            enddo
         endif
         if(btildeflg) then
            if(f0.gt.0) then
               if(f/f0.gt.prodrat) then
                  do kdim=1,ndim
                     ymaxrat(cells(kdim),kdim)=
     1                    ymaxrat(cells(kdim),kdim)*
     2                    (f/f0/prodrat+0.1)**(xmorerat/ndim)
                  enddo
               endif
            endif
         endif
         goto 1
      endif
c check if the failure rate is satisfactory
      fail=0
      tot=0
      ubtot=0
      if(btildeflg) then
         failrat=0
         totrat=0
         ubtotrat=0
      endif
 2    continue
      call getlinemintupb1(filetag,ndim,cells,f,f0,iret)
      if(iret.lt.0) goto 998
      if(iret.eq.0) then         
         prod=1
         if(btildeflg) prodrat=1
         do kdim=1,ndim
            prod=prod*ymax(cells(kdim),kdim)
            if(btildeflg) prodrat=prodrat*ymaxrat(cells(kdim),kdim)
         enddo
         if(f.gt.prod) fail=fail+(f-prod)
         tot=tot+f
         ubtot=ubtot+prod
         if(btildeflg) then
            if(f0.ne.0) then
               if(f.gt.f0*prodrat) failrat=failrat+(f-f0*prodrat)
               totrat=totrat+f
               ubtotrat=ubtotrat+f0*prodrat
            endif
         endif
         goto 2
      endif
      if(fail/tot.gt.1d-3.or.(btildeflg.and.failrat/totrat.gt.1d-3))then
c stop updating the rat grid, if satisfactory
         if(btildeflg.and.failrat/totrat.lt.1d-3) then
            btildeflg=.false.
            write(*,*) ' envelope efficiency',totrat/ubtotrat
            write(*,*) 'failure estimate',failrat/totrat
         endif
         if(iterations.lt.4) then
            write(*,*) ' iterating upper bounding envelope formation'
         elseif(iterations.lt.5) then
            write(*,*) ' more iterations needed'
            write(*,*) ' this can take a moment ...'
         endif
         write(*,*) 'failure estimate',fail/tot
         if(btildeflg) then
            write(*,*) 'ratios failure estimate',failrat/totrat
         endif
         iterations=iterations+1
         goto 1
      endif
      if(filetag.eq.'btildeupb') then
         btildeflg=.true.
      endif
      if(btildeflg) then
         write(*,*) 'envelope efficiency: '
         write(*,*)
     1 ' # of generated configurations over number of btilde calls'
         write(*,*) totrat/ubtotrat
         write(*,*) 'failure estimate',failrat/totrat
         write(*,*)
     1 ' # of generated configurations over number of born calls'
         write(*,*) tot/ubtot
         write(*,*) 'failure estimate',fail/tot
      else
         write(*,*) 'envelope efficiency=',tot/ubtot
         write(*,*) 'failure estimate',fail/tot
      endif
      write(*,*) 'processed ',ipoints,' points',iterations,'iterations'
      call newunit(iunit)
      open(unit=iunit,file='testbndrat.top',status='unknown')
      do kdim=1,ndim
         write(iunit,*) 'set limits x 0 50 y 0 5'
         do j=1,50
            write(iunit,*) j, ymaxrat(j,kdim)
         enddo
         write(iunit,*) 'hist'
         write(iunit,*) 'newplot'
      enddo
      close(iunit)
      open(unit=iunit,file='testbnd.top',status='unknown')
      do kdim=1,ndim
         write(iunit,*) 'set limits x 0 50 y 0 5'
         do j=1,50
            write(iunit,*) j, ymax(j,kdim)
         enddo
         write(iunit,*) 'hist'
         write(iunit,*) 'newplot'
      enddo
      close(iunit)
      call getlinemintupb1(filetag,ndim,cells,f,f0,-10)
      return
 998  continue
      write(*,*) ' error while loading bound files'
      call exit(-1)      
      end

      subroutine monitorubound(x,icalls)
      implicit none
      real * 8 x
      integer icalls
      include 'pwhg_rnd.h'
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      character * 80 file
      integer iun
      if(rnd_cwhichseed.eq.'none') then
         file=pwgprefix(1:lprefix)//'boundviolations.dat'
      else
         file=pwgprefix(1:lprefix)//'boundviolations-'//rnd_cwhichseed
     1        //'.dat'
      endif
      call newunit(iun)
      open(unit=iun,file=file,access='append')
      write(iun,*) 'calls=',icalls,'f/ubound',x
      close(iun)
      end


      function strictubound(cells)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flg.h'
      integer cells(ndiminteg)
      real * 8 strictubound
      include 'cgengrids.h'
      real * 8 prodrat
      integer kdim
      prodrat=1
      do kdim=1,ndiminteg
         prodrat=prodrat*ymaxrat(cells(kdim),kdim)
      enddo
      strictubound=prodrat
      end

      subroutine inttochar(int_arr,ch_arr,ndim)
      integer ndim
      integer int_arr(ndim)
      character * 1 ch_arr(ndim)
      integer k
      do k=1,ndim
         ch_arr(k)=char(int_arr(k))
      enddo
      end

      subroutine chartoint(ch_arr,int_arr,ndim)
      integer ndim
      integer int_arr(ndim)
      character * 1 ch_arr(ndim)
      integer k
      do k=1,ndim
         int_arr(k)=ichar(ch_arr(k))
      enddo
      end

