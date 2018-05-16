      program leshouchesanal
      implicit none
      include 'LesHouches.h'
      integer j,nev
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
c     let the analysis subroutine know that it is run by this program
      WHCPRG='LHE   '
      call opencount(nev)
      call upinit
      call init_hist 
      do j=1,nev
         call upevnt
         if(nup.eq.0) then
            write(*,*) ' nup = 0 skipping event'
            goto 111
         endif
         call lhuptohepevt(j)
         if(abs(idwtup).eq.3) xwgtup=xwgtup*xsecup(1)
         call analysis(xwgtup)
         call pwhgaccumup
         if (mod(j,20000).eq.0) then
            write(*,*) "# of events processed =",j
            open(unit=99,file='LHEF_analysis.top')
            call pwhgsetout
            call pwhgtopout
            close(99)
         endif
111     continue
      enddo
      open(unit=99,file='LHEF_analysis.top')
      call pwhgsetout
      call pwhgtopout
      close(99)
      write(*,*) 'EVENTS FOUND : ',nev
      end
      
      subroutine UPINIT
      implicit none
      double precision R_jet,ptmin_jet,powheginput
      external powheginput
      common/cjetdefs/R_jet,ptmin_jet
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   JET DEFINITIONS MANDATORY TO DEFINE THE PROCESS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      R_jet=powheginput('R_jet')
      ptmin_jet=powheginput('ptmin_jet')
      call lhefreadhdr(97)
      end

      subroutine UPEVNT
      call lhefreadev(97)
      end

      subroutine lhuptohepevt(n)
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      integer ihep,mu,n
      
      nhep=nup
      nevhep=n
      do ihep=1,nhep
         isthep(ihep)=istup(ihep)
         idhep(ihep)=idup(ihep)
         do mu=1,2
            jmohep(mu,ihep)=mothup(mu,ihep)
         enddo
         do mu=1,5
            phep(mu,ihep)=pup(mu,ihep)
         enddo
      enddo
      end
