      implicit none
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'pwhg_rnd.h'
      include 'pwhg_flg.h'
      integer j,iun,nev,maxev
      common/cnev/nev
      real * 8 weight,tmp
      real * 8 powheginput
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer ios
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      integer iseed,n1,n2
      logical testplots
      if (powheginput('#testplots').eq.1d0) then
         testplots=.true.
      else
         testplots=.false.
      endif
      nev=powheginput('numevts')
      call newunit(iun)
c The following allows to perform multiple runs with
c different random seeds in the same directory.
c If manyseeds is set to 1, the program asks for an integer j;
c The file 'pwgprefix'seeds.dat at line j is read, and the
c integer at line j is used to initialize the random
c sequence for the generation of the event.
c The event file is called 'pwgprefix'events-'j'.lhe
      if(powheginput("#manyseeds").eq.1) then
         open(unit=iun,status='old',iostat=ios,
     1        file=pwgprefix(1:lprefix)//'seeds.dat')
          if(ios.ne.0) then
             write(*,*) 'option manyseeds required but '
             write(*,*) 'file ',pwgprefix(1:lprefix)/
     $            /'seeds.dat not found'
            call exit(-1)
         endif 
         do j=1,1000000
            read(iun,*,iostat=ios)  rnd_initialseed
            if(ios.ne.0) goto 10
         enddo
 10      continue
         rnd_numseeds=j-1
         write(*,*) 'enter which seed'
         read(*,*) rnd_iwhichseed
         if(rnd_iwhichseed.gt.rnd_numseeds) then
            write(*,*) ' no more than ',rnd_numseeds, ' seeds in ',
     1           pwgprefix(1:lprefix)//'seeds.dat'
            call exit(-1)
         endif
         rewind(iun)
         do j=1,rnd_iwhichseed
c Commented line to be used instead, for testing that manyseed runs
c yield the same results as single seed runs, provided the total number
c of calls is the same.
c     read(iun,*) rnd_initialseed,rnd_i1,rnd_i2
            read(iun,*) rnd_initialseed
            rnd_i1=0
            rnd_i2=0
         enddo
         close(iun)
         write(rnd_cwhichseed,'(i4)') rnd_iwhichseed
         do j=1,4
            if(rnd_cwhichseed(j:j).eq.' ') rnd_cwhichseed(j:j)='0'
         enddo
      else
         rnd_cwhichseed='none'
      endif
      if (testplots) WHCPRG='NLO   '
      call pwhginit
      if(nev.gt.0) then
         if(flg_newweight) then
            if (testplots) then 
               write(*,*) '-------> Warning: testplots has been reset to
     1 false since we are doing reweighting' 
               testplots = .false. 
            endif
            continue
         else
            if(rnd_cwhichseed.ne.'none') then
               write(*,*) pwgprefix(1:lprefix)//'events-'//
     1            rnd_cwhichseed//'.lhe', rnd_iwhichseed,rnd_initialseed
               open(unit=iun,status='new',file=pwgprefix(1:lprefix)
     1              //'events-'//rnd_cwhichseed//'.lhe')
            else
               open(unit=iun,status='new',
     1              file=pwgprefix(1:lprefix)//'events.lhe')
            endif
         endif
      else
         write(*,*) ' No events requested'
         goto 999
      endif
      call lhefwritehdr(iun)
      if (testplots) then
         call init_hist 
c     let the analysis subroutine know that it is run by this program
         WHCPRG='LHE   '
      endif
c if we are using manyseeds, and iseed is given, it means that we want
c to examine that event in particular
      if(rnd_cwhichseed.ne.'none') then
         iseed=powheginput('#iseed')
         n1=powheginput('#rand1')
         n2=powheginput('#rand2')
         if(iseed.ge.0.and.n1.ge.0.and.n2.ge.0)
     1        call setrandom(iseed,n1,n2)
      endif
      call resetcnt
     1       ('upper bound failure in inclusive cross section')
      call resetcnt
     1       ('vetoed calls in inclusive cross section')
      call resetcnt(
     1 'upper bound failures in generation of radiation')
      call resetcnt('vetoed radiation')
      write(*,*)
      write(*,*)' POWHEG: generating events'
      if(flg_newweight) then
         call opencount(maxev)
         call openoutputrw
         if(maxev.ne.nev) then
            write(*,*) ' Warning: powheg.input says ',nev,' events'
            write(*,*) ' the file contains ', maxev, ' events'
            write(*,*) ' Doing ',maxev,' events'
            nev = maxev 
         endif
      endif
      do j=1,nev
         if(flg_newweight) then
            call pwhgnewweight
         else
            call pwhgevent
            if(nup.eq.0) then
               write(*,*) ' nup = 0 skipping event'
               goto 111
            endif
            call lhefwritev(iun)
         endif
         if(idwtup.eq.3) then
            weight=rad_totgen*xwgtup*rad_branching
         elseif(idwtup.eq.-4) then
            weight=xwgtup
         else
            write(*,*) ' only 3 and -4 are allowed for idwtup'
            call exit(-1)
         endif
         if(testplots) then
            call lhtohep
            call analysis(weight)
            call pwhgaccumup
            if (mod(j,20000).eq.0) then
               if(rnd_cwhichseed.eq.'none') then
                  open(unit=99,file=pwgprefix(1:lprefix)//
     1                 'pwhgalone-output.top')
               else
                  open(unit=99,file=pwgprefix(1:lprefix)//
     1                 'pwhgalone-output'//rnd_cwhichseed//'.top')
               endif
               call pwhgsetout
               call pwhgtopout
               close(99)
            endif
         endif
 111     continue
      enddo
      if (testplots) then
         if(rnd_cwhichseed.eq.'none') then
            open(unit=99,file=pwgprefix(1:lprefix)//
     1           'pwhgalone-output.top')
         else
            open(unit=99,file=pwgprefix(1:lprefix)//
     1           'pwhgalone-output'//rnd_cwhichseed//'.top')
         endif
         call pwhgsetout
         call pwhgtopout
         close(99)
      endif
      call lhefwritetrailer(iun)
      close(iun)
 999  continue
      call write_counters
c this causes powheginput to print all unused keywords
c in the powheg.input file; useful to catch mispelled keywords
      tmp=powheginput('print unused tokens')
      end
