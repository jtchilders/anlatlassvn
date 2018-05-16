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
            call lheanend
         endif
111     continue
      enddo
      call lheanend
      write(*,*) 'EVENTS FOUND : ',nev
      end

      subroutine lheanend
      implicit none
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      include 'pwhg_rnd.h'
      if(rnd_cwhichseed.ne.'none') then
         open(unit=99,file=pwgprefix(1:lprefix)//'LHEF_analysis-'
     1        //rnd_cwhichseed//'.top'
     2     ,status='unknown')
      else
         open(unit=99,file=pwgprefix(1:lprefix)//'LHEF_analysis.top'
     1     ,status='unknown')
      endif
      call pwhgsetout
      call pwhgtopout
      close(99)
      end
      
      subroutine UPINIT
      implicit none
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
