      program main_pythia
      implicit none
      include 'LesHouches.h'
      include 'hepevt.h'
      include 'PhysPars.h'
c      integer NMXHEP,NEVHEP,NHEP,ISTHEP,IDHEP,
c     &     JMOHEP,JDAHEP
c      double precision phep,vhep
c      PARAMETER (NMXHEP=4000)
c      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
c     &     JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),
c     &     VHEP(4,NMXHEP)
      real * 8 parp,pari
      integer mstp,msti
      common/pypars/mstp(200),parp(200),msti(200),pari(200)
      integer MSTU,MSTJ
      double precision PARU,PARJ
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      integer maxev
      common/mcmaxev/maxev
      integer MDCY,MDME,KFDP
      double precision brat
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      integer pycomp
      external pycomp
      integer iev,temp,i
      external pydata
      real * 8 powheginput
      external powheginput
      integer hdecaymode
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      
      WHCPRG='PYTHIA'
c     set Higgs boson mass and width
      call init_couplings

c PARAMETERS
c     mstj(41)=3 !No photon radiations off leptons
c     mstp(61)=0 !No IS shower
c     mstp(71)=0 !No FS shower
      mstp(81)=0 !No Multiple interactions (MI increases the execution time)
c     mstp(91)=0  !No Primordial kt
c     mstp(131)=0 !No Pile Up
     mstp(111)=0               !No hadronization
      
c     Make PI0 stable as in herwig default
      mdcy(pycomp(111),1)=0
c     Make tau stable
      mdcy(pycomp(15),1)=0
C---PROCESS; 
c      hdecaymode=powheginput('hdecaymode')
c      if ((hdecaymode.lt.-1).or.(hdecaymode.gt.12)) then
c         write(*,*) "Higgs decay mode not allowed"
c         stop
c      endif   
      hdecaymode=-1

c     choose Higgs decay channel
      if (hdecaymode.eq.-1) then
         mdcy(pycomp(25),1)=0
      else   
         mdcy(pycomp(25),1)=1
         if (hdecaymode.gt.0) then
            do i=210,288
               if (mdme(i,1).ne.-1) mdme(i,1)=0
            enddo
            if (hdecaymode.eq.12) then
               mdme(223,1)=1
            elseif(hdecaymode.eq.11) then
               mdme(225,1)=1
            elseif(hdecaymode.eq.10) then
               mdme(226,1)=1
            elseif(hdecaymode.eq.7) then
               mdme(218,1)=1
            elseif(hdecaymode.eq.8) then
               mdme(219,1)=1
            elseif(hdecaymode.eq.9) then
               mdme(220,1)=1   
            else
               mdme(209+hdecaymode,1)=1
            endif
         endif      
      endif
      
c Set up PYTHIA to accept user processes
      call PYINIT('USER','','',0d0)
      
      call PYABEG
      nevhep=0
      do iev=1,maxev
         call pyevnw
         if(nup.eq.0) goto 111
c     Convert from PYJETS event record to HEPEVT event record
         temp=nevhep
         call pyhepc(1)
         nevhep=temp
C     Print out the event record
         IF (IEV.le.1) THEN 
c     list the event
c     CALL PYLIST(7) ! print the HEPEUP common block
c     CALL PYLIST(5) ! print the HEPEVT common block
            CALL PYLIST(5)      ! print the event
c     call PYLIST(1) ! as PYLIST(2) but with less information
         ENDIF
         
         call PYANAL
         
         IF(MOD(IEV,20000).EQ.0) THEN
            WRITE(*,*)'# of events processed=',NEVHEP
            CALL PYAEND
         ENDIF          
      enddo
 111  continue
      write(*,*) 'At the end NEVHEP is ',nevhep
      call PYAEND
      END


      subroutine UPINIT
      implicit none
      integer maxev
      common/mcmaxev/maxev      
C--- Opens input file and counts number of events, setting MAXEV;
      call opencount(maxev)
      rewind(97)
      call lhefreadhdr(97)
      end


      subroutine countevents
      implicit none
      character * 6 string
      integer maxev
      common/mcmaxev/maxev
 1    continue
      read(unit=97,fmt='(a)',end=2) string
      if(string.eq.'<event') then
         maxev=maxev+1
         goto 1
      endif
      goto 1
 2    continue
      write(*,*) ' found ',maxev,' events in file'
      end


      subroutine UPEVNT
      call lhefreadev(97)
      end


      subroutine upveto
c pythia routine to abort event
      end


      subroutine pyabeg
      call init_hist
      end

c      subroutine pyabeg
c      implicit none
c      include 'hepevt.h'
c      nevhep=0
c      call abegin
c      end

      subroutine pyaend
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      open(unit=99,file=pwgprefix(1:lprefix)//'POWHEG+PYTHIA-output.top'
     #     ,status='unknown')
      call pwhgsetout
      call pwhgtopout
      close(99)
      end


      subroutine pyanal
      implicit none
      include 'LesHouches.h'
      include 'hepevt.h'
c      integer NMXHEP,NEVHEP,NHEP,ISTHEP,IDHEP,
c     &     JMOHEP,JDAHEP
c      double precision phep,vhep
c      PARAMETER (NMXHEP=4000)
c      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
c     &   JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      integer MDCY,MDME,KFDP
      double precision brat
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      integer pycomp
      external pycomp
      real * 8 powheginput,bratio
      external powheginput
      integer hdecaymode

      nevhep=nevhep+1
c      hdecaymode=powheginput('hdecaymode')
      hdecaymode=-1
      if (hdecaymode.eq.0) then
         bratio=1d0
      elseif (hdecaymode.eq.12) then
         bratio=brat(223)
      elseif(hdecaymode.eq.11) then
         bratio=brat(225)
      elseif(hdecaymode.eq.10) then
         bratio=brat(226)
      elseif(hdecaymode.eq.7) then
         bratio=brat(218)
      elseif(hdecaymode.eq.8) then
         bratio=brat(219)
      elseif(hdecaymode.eq.9) then
         bratio=brat(220) 
      elseif(hdecaymode.eq.-1) then
         bratio=1d0  
      else
         bratio=brat(209+hdecaymode)
      endif
      if(idwtup.eq.3) xwgtup=xwgtup*xsecup(1)
      xwgtup=xwgtup*bratio
      call analysis(xwgtup)
      call pwhgaccumup 
      end




