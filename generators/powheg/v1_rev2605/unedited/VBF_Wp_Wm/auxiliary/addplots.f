c combine multiseed powheg output pwgNLO-????.top.
      implicit none
      character * 60 tag, files(500)
      integer npoints,j,k,l,ncalls,npt,nmax,nfiles,firstfile,jj, index 
      parameter (nmax=500)
      real * 8 xxx(nmax),yyy(nmax),eee(nmax)
      real * 8 xx(nmax),yy(nmax),ee(nmax),xmin,xmax
      write(*,*) ' enter number of files'
      read(*,*) nfiles
      write(*,*) ' enter number of calls'
      read(*,*) ncalls
      write(*,*) ' enter number of first file'
      read(*,*) firstfile
      if(nfiles.gt.500) then
         write(*,*) ' number exceeds maximum limit, increase it'
         stop
      endif
      do j=1,nfiles
         files(j)='wwNLO-0001.top'
         jj=j+firstfile-1
         if(jj.lt.10) then
            write(files(j)(10:10),'(i1)') jj
         elseif(jj.lt.100) then
            write(files(j)(9:10),'(i2)') jj
         else
            write(files(j)(8:10),'(i3)') jj
         endif
      enddo
      open(unit=11,file=files(1),status='old')
      open(unit=12,file='addedNLO.top',status='unknown')
      open(unit=112,file='addedNLO-gnu.top',status='unknown')
      npoints=nmax
      index = 0 
 1    do k=1,1000
         call getdata(11,12,tag,xmin,xmax,npoints,xxx,yyy,eee)
         if(npoints.eq.0) goto 999
         call completetop(npoints,nmax,xmin,xmax,xxx,yyy,eee)
         write(*,*) tag
         do j=2,nfiles
            npt=nmax
            call getdatafromtag(files(j),tag,npt,xx,yy,ee)
            call completetop(npt,nmax,xmin,xmax,xx,yy,ee)
            if(npt.ne.npoints) then
               write(*,*) ' error: ',files(j)
               call exit(-1)
            endif
            do l=1,npoints
               if(xx(l).ne.0.or.xxx(l).ne.0) then
                  if(abs(xx(l)-xxx(l))/abs(xx(l)+xxx(l)).gt.1d-13) then
                     write(*,*) ' error: two data sets do not match'
                     call exit(-1)
                  endif
               endif
               eee(l)=sqrt((eee(l)**2*(j-1)**2+ee(l)**2)/j**2
     1              +(j-1)*(yyy(l)-yy(l))**2/(j**3*ncalls))
               yyy(l)=(yyy(l)*(j-1)+yy(l))/j
            enddo
         enddo
         do j=1,npoints
            yyy(j) = nfiles*yyy(j)
            eee(j) = nfiles*eee(j)
            write(12,'(f13.4,2x,2(e14.8,2x))') xxx(j),yyy(j),eee(j)
         enddo
         write(112,'(a9,i2,a)') '# index ', index, tag(4:)  
         do j=1,npoints
            yyy(j) = yyy(j)
            eee(j) = eee(j)
            write(112,'(f13.4,2x,2(e14.8,2x))') xxx(j),yyy(j),eee(j)
         enddo
         write(112,*) 
         write(112,*) 
         index = index +1 
      enddo
 999  end
   

      subroutine writeout(iun,line)
      integer iun
      character *(*) line
      integer n,j
      n=len(line)
      do j=n,1,-1
         if(line(j:j).ne.' ') then
            write(iun,'(a)') line(1:j)
            return
         endif
      enddo
      end



      subroutine getdata(iin,iout,tag,xmin,xmax,npt,x,y,e)
      implicit none
      integer iin,iout,npt
      character *(*) tag
      real * 8 x(npt),y(npt),e(npt),xmin,xmax
      character * 100 line,oline
      real * 8 rarr(10)
      integer narr(10),karr,j,l1,l2,k
      tag=' '
      oline=' '
      do k=1,1000000
         read(unit=iin,fmt='(a)',err=99,end=99) line
         l1=index(line,'"')
         if(l1.gt.0) then
            l2=index(line(l1+1:),'"')
            if(l2.gt.0) tag=line(l1+1:l2-1)
         endif
         call writeout(12,line)
         if(line(1:7).eq.' ( INT=') then
            backspace(iin)
            backspace(iin)
            backspace(iin)
            backspace(iin)
            backspace(iin)
            read(iin,'(a)') oline
            if(oline(1:15).eq.'  SET LIMITS X ') then
               read(oline(16:),*) xmin,xmax
            else
               write(*,*) ' error: not a set limits line'
               call exit(-1)
            endif
            read(iin,*)
            read(iin,*)
            read(iin,'(a)') tag
            read(iin,*)
            do j=1,1000000
               read(unit=iin,fmt='(a)',err=99,end=99) line
               call reads(line,100,rarr,narr,karr)
               if(karr.eq.3) then
                  x(j)=rarr(1)
                  y(j)=rarr(2)
                  e(j)=rarr(3)
               else
                  npt=j-1
                  backspace(iin)
                  return
               endif
            enddo
         endif
         oline=line
      enddo
 99   continue
      npt=0
      end

      subroutine completetop(npt,nmax,xmin,xmax,x,y,e)
      implicit none
      integer npt,nmax
      real * 8 x(nmax),y(nmax),e(nmax),xmin,xmax
      real * 8 binsz
      integer ioffset,nempty,j,jj,k,ntop
      if(npt.eq.0) then
         write(*,*) ' completetop: 0 points!'
         call exit(-1)
      endif
      binsz=1d15
      do j=2,npt
         binsz=min(binsz,abs(x(j)-x(j-1)))
      enddo
      binsz=min(binsz,2*(x(1)-xmin),2*(xmax-x(npt)))
      ioffset=0
      ntop=npt
      do j=1,npt+1
         jj=j+ioffset
         ntop=npt+ioffset
         if(j.eq.1) then
            nempty=nint((x(1)-(xmin-binsz/2))/binsz)-1
         elseif(j.lt.npt+1) then
            nempty=nint((x(jj)-x(jj-1))/binsz)-1
         else
            nempty=nint((xmax+binsz/2-x(ntop))/binsz)-1
         endif
         if(ntop+nempty.gt.nmax) then
            write(*,*) ' completetop: arrays x,y,e too small'
            stop
         endif
         if(nempty.gt.0) then
            do k=ntop,jj,-1
               x(k+nempty)=x(k)
               y(k+nempty)=y(k)
               e(k+nempty)=e(k)
            enddo
            do k=0,nempty-1
               if(jj.eq.1) then
                  x(jj+k)=xmin-binsz/2+(k+1)*binsz
               else
                  x(jj+k)=x(jj-1)+(k+1)*binsz
               endif
               y(jj+k)=0
               e(jj+k)=0
            enddo
         endif
         ioffset=ioffset+nempty
      enddo
      ntop=npt+ioffset
      npt=ntop
      do k=npt+1,nmax
         x(k)=0
         y(k)=0
         e(k)=0
      enddo
      end

      subroutine getdatafromtag(file,tag,npt,x,y,e)
      implicit none
      integer npt
      character *(*) tag, file
      real * 8 x(npt),y(npt),e(npt)
      character * 100 line,oline
      real * 8 rarr(10)
      integer narr(10),karr,j,l1,l2,k
      open(unit=13,file=file,status='old')
      do k=1,1000000
         read(unit=13,fmt='(a)',err=99,end=99) line
         if(line(1:len(tag)).eq.tag) then
            read(13,*)
            do j=1,1000000
               read(unit=13,fmt='(a)',err=99,end=99) line
               call reads(line,100,rarr,narr,karr)
               if(karr.eq.3) then
                  x(j)=rarr(1)
                  y(j)=rarr(2)
                  e(j)=rarr(3)
               else
                  npt=j-1
                  close(13)
                  return
               endif
            enddo
         endif
      enddo
 99   continue
      npt=0
      close(13)
      end

c Program to read numbers from strings
      subroutine reads(string,nstr,rarr,narr,karr)
c string: input string of nstr characters
c it outputs in rarr the real numbers in string,
c karr is the number of real numbers found
      implicit none
      integer nstr,narr,karr
      real * 8 rarr(narr)
      character *(*) string
c
      character * 15 numstr
      integer istr,istart,ios
      karr=0
      istr=1
 1    continue
c no more numbers:
      if(istr.gt.nstr) goto 999
c skip blanks
 10   if(string(istr:istr).eq.' ') then
         if(istr.eq.nstr) goto 999
         istr=istr+1
         goto 10
      endif
c first non-blank
      istart=istr
c find next non-blank
 20   if(string(istr:istr).ne.' ') then
         if(istr.eq.nstr) goto 22
         istr=istr+1
         goto 20
      endif
      istr=istr-1
 22   continue
      if(istr-istart+1.gt.15) then
         write(*,*) ' error: number too long'
         stop
      endif
      numstr=string(istart:istr)
      karr=karr+1
      read(unit=numstr,iostat=ios,fmt='(bn,f15.0)') rarr(karr)
c if not a valid number ignore silently
      if(ios.ne.0) then
         karr=karr-1
         goto 999
      endif
      istr=istr+1
      goto 1
 999  end

c Program to read integers from strings
      subroutine ireads(string,nstr,iarr,narr,karr)
      implicit none
      integer narr
      integer iarr(narr)
      character * (*) string
      integer nstr,karr,j
c local variables
      real * 8 rarr(20)
      call reads(string,nstr,rarr,20,karr)
      if(karr.gt.narr.or.(narr.gt.20.and.karr.eq.20)) then
         write(*,*) ' more than 20 numbers in input?'
         stop
      endif
      do j=1,karr
         iarr(j)=nint(rarr(j))
         if(rarr(j)-iarr(j).ne.0) then
            write(*,*) ' not an integer in input!'
            stop
         endif
      enddo
      end

c program to read numbers from input line
      subroutine iread(iun,iarr,nel,nread)
      dimension iarr(*)
      character * 300 str
      read(iun,'(a)') str
      call ireads(str,300,iarr,nel,nread)
      end

c program to read numbers from input line
      subroutine rread(iun,rarr,nel,nread)
      real * 8  rarr(*)
      character * 300 str
      read(iun,'(a)') str
      call reads(str,300,rarr,nel,nread)
      end
      
c      real * 8 rarr(10)
c      integer iarr(10)
c 1    call rread(5,rarr,10,k)
c      write(*,*) k,' values, ',(rarr(j),j=1,k)
c      goto 1
c      end
