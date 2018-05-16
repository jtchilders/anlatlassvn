      program main_pythia
      implicit none
      include 'LesHouches.h'
      include 'hepevt.h'
 
      integer iev,temp,i
      external pydata
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      integer maxpr
      parameter (maxpr=6)
c     mcmaxev
      integer maxev
      common/mcmaxev/maxev

c     WHCPRG tells the analysis subroutine which program is calling the
c     analysis
      WHCPRG='PYTHIA'

      call getmaxev(maxev)

c Set up tune
      call setup_PYTHIA_tune

c Set up PYTHIA to accept user processes
      call PYINIT('USER','','',0d0)
      
c Set up initial parameter      
      call setup_PYTHIA_parameters
      
      call PYABEG
      nevhep=0
      do iev=1,maxev
         call pyevnt
         if(nup.eq.0) then
            write(*,*) ' no event generated; skipping'
            goto 111
         endif
c     Convert from PYJETS event record to HEPEVT event record
         temp=nevhep
         call pyhepc(1)
         nevhep=temp
C     Print out the event record
         IF (IEV.le.maxpr) THEN 
c     list the event
c            call pystat(2)      ! print cross sections, widths, branchings,...
c            CALL PYLIST(7)      ! print the HEPEUP common block
             CALL PYLIST(5)      ! print the HEPEVT common block
c            CALL PYLIST(2)      ! print the event
c            call PYLIST(1)      ! as PYLIST(2) but with less information
         ENDIF
         
         call PYANAL
          IF(nevhep.gt.0.and.MOD(nevhep,20000).EQ.0) THEN
            WRITE(*,*)'# of events processed=',iev
            WRITE(*,*)'# of events generated=',NEVHEP
            CALL PYAEND
         ENDIF 
      enddo
 111  continue
      write(*,*) 'At the end NEVHEP is ',nevhep
!:      write(*,*) 'At the end: #warnings= ',mstu(27),' #errors= ',mstu(23)
C---USER'S TERMINAL CALCULATIONS
      call PYAEND
      END

