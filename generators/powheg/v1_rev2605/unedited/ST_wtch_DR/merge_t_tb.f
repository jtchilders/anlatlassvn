      program merge_t_tb
      implicit none
      integer input
      read(*,*) input
      if(input.gt.0) then
         call find_totals(input)
      elseif(input.lt.0) then
         call merge(input)
      else
         write(*,*) 'Error: wrong input value in merge_t_tb: ',input
         call exit(1)
      endif
      end

      subroutine find_totals(total_ev)
      implicit none
      integer total_ev
      include 'LesHouches.h'
      integer iun_t,iun_tb,ios,nev_t,nev_tb
      character * 50 name_t,name_tb
      double precision Xsec_t,Xsec_tb,Xsec_tot
      integer id(2)
      double precision eb(2)
      name_t='temp_t-events.lhe'
      name_tb='temp_tb-events.lhe'
      call newunit(iun_t)
      open(unit=iun_t,file=name_t,status='old',iostat=ios)
      if(ios.ne.0) then
         write(*,*)' Problems opening file '//name_t
         call exit(1)
      endif
      call newunit(iun_tb)
      open(unit=iun_tb,file=name_tb,status='old',iostat=ios)
      if(ios.ne.0) then
         write(*,*)' Problems opening file '//name_tb
         call exit(1)
      endif

      Xsec_t  =0d0
      call lhefreadhdr(iun_t)
      Xsec_t = xsecup(1)
      id(1) = idbmup(1)
      id(2) = idbmup(2)
      eb(1) = ebmup(1)
      eb(2) = ebmup(2)
      
      Xsec_tb =0d0
      call lhefreadhdr(iun_tb)
      Xsec_tb = xsecup(1)
c     sanity check
      if (idbmup(1).ne.id(1).and.idbmup(2).ne.id(2).and.
     #     ebmup(1).ne.eb(1).and.ebmup(2).ne.eb(2)) then
         write(*,*) 'The two input files have different beam energy '//
     #        'or beam particles'
         call exit(1)
      endif

c     close opened files
      rewind (iun_t)
      close (iun_t)
      rewind (iun_tb)
      close (iun_tb)

      write(*,*) 'This file is read by the merge_ttb script'
      write(*,*) Xsec_t, '       ! Xsecs t'
      write(*,*) Xsec_tb,'       ! Xsecs tb'
      Xsec_tot=Xsec_t+Xsec_tb
      write(*,*) Xsec_tot,'       ! Xsecs tot'

      nev_t=nint(1.05*(Xsec_t/Xsec_tot)*total_ev)
      nev_tb=nint(1.05*(Xsec_tb/Xsec_tot)*total_ev)
      write(*,*) nev_t,       '       ! 1nev_t'
      write(*,*) nev_tb,      '       ! 2nev_tb'
      write(*,*) nev_t+nev_tb,'       ! nev_t+nev_tb'
      write(*,*) total_ev,    '       ! total_ev'

c     check that total number of events is not lower
c     than the required one
      if((nev_t+nev_tb).lt.total_ev) then
         write(*,*) 'Error in find_totals: nev_t+nev_tb<total_ev '
         call exit(1)
      endif

      end
      

      subroutine merge(total_ev)
      implicit none
      integer total_ev
      include 'LesHouches.h'
      integer iun_t,iun_tb,iun_merge,ios,nev_t,nev_tb,iev,nt,ntb
      character * 50 name_t,name_tb,name_merge
      double precision Xsec_t,Xsec_tb,Xsec_tot,errXsec_t,errXsec_tb
     $     ,errXsec_tot,r,ratio
      integer id(2)
      double precision eb(2)
      double precision random
      external random

      if(total_ev.ge.0) then
         write(*,*) 'Error: input of merge subroutine must be negative'
         call exit(1)
      else
         total_ev=-total_ev
      endif

      name_t='temp_t-events.lhe'
      name_tb='temp_tb-events.lhe'
      name_merge='t_tb_sample-events.lhe'
      call newunit(iun_t)
      open(unit=iun_t,file=name_t,status='old',iostat=ios)
      if(ios.ne.0) then
         write(*,*)' Problems opening file '//name_t
         call exit(1)
      endif
      call newunit(iun_tb)
      open(unit=iun_tb,file=name_tb,status='old',iostat=ios)
      if(ios.ne.0) then
         write(*,*)' Problems opening file '//name_tb
         call exit(1)
      endif
      call newunit(iun_merge)
      open(unit=iun_merge,file=name_merge,status='unknown')   
      write(*,*) '*Merging started'

      Xsec_t  =0d0
      call lhefreadhdr(iun_t)
      Xsec_t = xsecup(1)
      errXsec_t = xerrup(1)
      id(1) = idbmup(1)
      id(2) = idbmup(2)
      eb(1) = ebmup(1)
      eb(2) = ebmup(2)
      
      Xsec_tb =0d0
      call lhefreadhdr(iun_tb)
      Xsec_tb = xsecup(1)
      errXsec_tb = xerrup(1)

      Xsec_tot=Xsec_t+Xsec_tb
      errXsec_tot=sqrt(errXsec_t**2+errXsec_tb**2)
      
      nev_t=nint(1.05*(Xsec_t/Xsec_tot)*total_ev)
      nev_tb=nint(1.05*(Xsec_tb/Xsec_tot)*total_ev)

      write(*,*) ' top sample events: ',nev_t
      write(*,*) ' antitop sample events: ',nev_tb

c     header of merged file
      xsecup(1) = Xsec_tot
      xerrup(1) = errXsec_tot
c     special lprup for merged file
      lprup(1) = 666
      call lhefwritehdr_merged(iun_merge,
     $     19,'temp_t-powheg.input',
     $     20,'temp_tb-powheg.input')
      call flush(iun_merge)
c     creating merged file
      ratio=Xsec_t/Xsec_tot

      nt=0
      ntb=0
      do iev=1,total_ev
         r=random()
         if(r.lt.ratio) then
c     top event
            nt=nt+1
            if(nt.gt.nev_t) then
               write(*,*) '* Merging ended with ',iev-1
     $              ,' events in the sample'
               write(*,*) '  Top events were not enough'
               goto 999
            endif
            call lhefreadev(iun_t)
            if(nup.eq.0) then
               write(*,*) '* Merging ended with ',iev-1
     $              ,' events in the sample'
               write(*,*) '  Top events were not enough'
               goto 999
            endif
            idprup=666
            call lhefwritev(iun_merge)
         else
c     antitop event
            ntb=ntb+1
            if(ntb.gt.nev_tb) then
               write(*,*) '* Merging ended with ',iev-1
     $              ,' events in the sample'
               write(*,*) '  Antitop events were not enough'
               goto 999
            endif
            call lhefreadev(iun_tb)
            if(nup.eq.0) then
               write(*,*) '* Merging ended with ',iev-1
     $              ,' events in the sample'
               write(*,*) '  Antitop events were not enough'
               goto 999
            endif
            idprup=666
            call lhefwritev(iun_merge)
         endif
      enddo

      write(iun_merge,'(a)') '</LesHouchesEvents>'
      write(*,*) '* Merging ended regularly'
      write(*,*) '  Merged sample has ',total_ev,' events'

 999  continue

c     close opened files
      rewind(iun_t)
      close(iun_t)
      rewind(iun_tb)
      close(iun_tb)
      rewind(iun_merge)
      close(iun_merge)

      end

c...lhefheader(nlf)
c...reads initialization information from a les houches events file on unit nlf. 
      subroutine lhefreadhdr(nlf)
      implicit none
      integer nlf
      character * 100 string
      integer ipr
      include 'LesHouches.h'
 1    read(nlf,fmt='(a)',err=998,end=998) string
      if(string(1:5).eq.'<init') then
         read(nlf,*) idbmup(1),idbmup(2),ebmup(1),ebmup(2),
     &        pdfgup(1),pdfgup(2),pdfsup(1),pdfsup(2),idwtup,nprup
         do ipr=1,nprup
            read(nlf,*) xsecup(ipr),xerrup(ipr),xmaxup(ipr),
     &           lprup(ipr)
         enddo
         goto 999
      else
         goto 1
      endif
 998  write(*,*) 'lhefreadhdr: could not find <init> data'
      call exit(1)
 999  end


c...reads event information from a les houches events file on unit nlf. 
      subroutine lhefreadev(nlf)
      implicit none
      integer nlf
      character * 100 string
      include 'LesHouches.h'
      integer i,j
 1    read(nlf,fmt='(a)',err=998,end=998) string
      if(string.eq.'</LesHouchesEvents>') then
         goto 998
      endif
      if(string(1:6).eq.'<event') then
         read(nlf,*) nup,idprup,xwgtup,scalup,aqedup,aqcdup
         do i=1,nup
            read(nlf,*) idup(i),istup(i),mothup(1,i),
     &           mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     &           vtimup(i),spinup(i)
         enddo
         call lhefreadpdfrw(nlf)
         goto 999
      else
         goto 1
      endif
c no event found:
 998  nup=0
 999  end


c...lhefwritehdr(nlf)
c...writes initialization information to a les houches events file on unit nlf. 
      subroutine lhefwritehdr_merged(nlf,lnamet,namet,lnametb,nametb)
      implicit none
      integer nlf
      integer lnamet,lnametb
      character * 50 namet,nametb
      real * 8 version
      common/cversion/version
      data version/1.0/
      integer ipr,iran,n1ran,n2ran
      include 'LesHouches.h'
      write(nlf,'(a)') '<LesHouchesEvents version="1.0">'
      write(nlf,'(a)') '<!--'
      write(nlf,'(a,f3.1)') 'file generated with POWHEG-BOX version ',
     #     version
      write(nlf,'(a)')
     $     '*** Events obtained merging top and antitop samples ***'
      write(nlf,'(a)') '*** Input file for top events was: ***'
      call wrtpowheginput_merged(nlf,lnamet,namet)
      write(nlf,'(a)') '*** Input file for antitop events was: ***'
      call wrtpowheginput_merged(nlf,lnametb,nametb)
      write(nlf,'(a)') '*** End of input files ***'
      call rm48ut(iran,n1ran,n2ran)
      write(nlf,*) 'Random number generator initialized with: ',
     # iran,' ',n1ran,' ',n2ran
      write(nlf,'(a)') '-->'
      write(nlf,'(a)') '<init>'
      write(nlf,110) idbmup(1),idbmup(2),ebmup(1),ebmup(2),
     &pdfgup(1),pdfgup(2),pdfsup(1),pdfsup(2),idwtup,nprup
      do 100 ipr=1,nprup
         write(nlf,120) xsecup(ipr),xerrup(ipr),xmaxup(ipr),
     &        lprup(ipr)
 100  continue
      write(nlf,'(a)') '</init>'
 110  format(1p,2(1x,i8),2(1x,e12.5),6(1x,i6))
 120  format(1p,3(1x,e12.5),1x,i6)
      end

c...lhefeader(nlf)
c...writes event information to a les houches events file on unit nlf. 
      subroutine lhefwritev(nlf)
      implicit none
      integer nlf
      include 'LesHouches.h'
      include 'pwhg_flg.h'
      integer i,j
      write(nlf,'(a)')'<event>'
      write(nlf,210) nup,idprup,xwgtup,scalup,aqedup,aqcdup
      do 200 i=1,nup
         write(nlf,220) idup(i),istup(i),mothup(1,i),
     & mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     & vtimup(i),spinup(i)
 200  continue
      if(flg_pdfreweight) call lhefwritepdfrw(nlf)
      write(nlf,'(a)')'</event>'      
 210  format(1p,2(1x,i6),4(1x,e12.5))
 220  format(1p,i8,5(1x,i5),5(1x,e16.9),1x,e12.5,1x,e10.3)
      end


      subroutine wrtpowheginput_merged(output,lname,filename)
      implicit none
      integer output,lname
      character * 50 filename
      integer iunit,k,ios
      character * 100 line
      call newunit(iunit)
      open(unit=iunit,file=filename(1:lname),status='old',iostat=ios)
      if(ios.ne.0) then
         write(*,*)' Error in wrtpowheginput_merged'
         write(*,*)' while opening file ',filename(1:lname)
         call exit(1)
      endif
 1    continue
      read(unit=iunit,fmt='(a)',iostat=ios,end=999) line
      if(ios.ne.0) then
         write(*,*) ' cannot read ',filename(1:lname)
         call exit(1)
      endif
      k=100
 2    if(k.gt.0) then
         if(line(k:k).eq.' ') then
            k=k-1
            goto 2
         endif
      endif
      if(k.eq.0) k=1
      write(output,'(a)') line(1:k)
      goto 1
      close(iunit)
 999  end



      subroutine lhefreadpdfrw(nlf)
      implicit none
      include 'pwhg_flg.h'
      integer nlf
      integer id1,id2
      real * 8 x1,x2,xf1,xf2,xmufact
      common/cpdfrwinfo/id1,id2,x1,x2,xmufact,xf1,xf2
      character *4 pdftag
      read(nlf,*,err=999,end=999) pdftag,id1,id2,x1,x2,xmufact,xf1,xf2
      if(pdftag.eq.'#pdf') then
         flg_pdfreweight=.true.
      else
         flg_pdfreweight=.false.
      endif
      return
 999  write(*,*) 'Error in lhefreadpdfrw (merging)',pdftag
      call exit(1)
      end


      subroutine lhefwritepdfrw(nlf)
      implicit none
      integer nlf
      integer id1,id2
      real * 8 x1,x2,xf1,xf2,xmufact
      common/cpdfrwinfo/id1,id2,x1,x2,xmufact,xf1,xf2
      write(nlf,111)'#pdf ',id1,id2,x1,x2,xmufact,xf1,xf2
 111  format(a,2(1x,i2),5(1x,d14.8))
      end





