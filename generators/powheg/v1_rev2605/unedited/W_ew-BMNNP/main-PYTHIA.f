      program main_pythia
      implicit none
      include 'LesHouches.h'
      include 'hepevt.h'
 
      integer iev,temp,i
      external pydata
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      integer maxpr
      parameter (maxpr=1)
c     mcmaxev
      integer maxev
      common/mcmaxev/maxev
      integer id_w

      logical ok
      common/ok_photos/ok

      real * 8 powheginput
      external powheginput

      real *8 kt2minqed
      common/showerqed/kt2minqed

      integer iun
c set epsilon
      kt2minqed = powheginput("#kt2minqed")
      if (kt2minqed.le.0d0) kt2minqed  = 0.001d0**2

c     WHCPRG tells the analysis subroutine which program is calling the
c     analysis
      WHCPRG='PYTHIA'

      call getmaxev(maxev)

c Set up tune
      call setup_PYTHIA_tune

c Set up PYTHIA to accept user processes
      call PYINIT('USER','','',0d0)
      call phoini
      
c Set up initial parameter      
      call setup_PYTHIA_parameters
      
      call PYABEG
      nevhep=0
      do iev=1,maxev
         call lhefreadev(97)
         call lhuptophhepevt(iev)
         call seteps
         call find_w(id_w)

         ok = .false.
c PHOTOS UNTIL EVENT OK
         if (.not.ok) then
             call lhuptophhepevt(iev)
             call photos_make(id_w)
             call confronta_pt
         endif
 
         call translate
         call hepevttolhef

c PHYTHIA
         call pyevnt
c         if(nup.eq.0) goto 111
c     Convert from PYJETS event record to HEPEVT event record
         temp=nevhep
         call pyhepc(1)
         nevhep=temp
c     Print out the event record
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
            WRITE(*,*)'# of events processed=',NEVHEP
            CALL PYAEND
         ENDIF 
      enddo
 111  continue
      write(*,*) 'At the end NEVHEP is ',nevhep
!:      write(*,*) 'At the end: #warnings= ',mstu(27),' #errors= ',mstu(23)
C---USER'S TERMINAL CALCULATIONS
      call PYAEND
      END
*
**
*
      subroutine hepevttolhef
      implicit none
      include 'LesHouches.h'
      include 'hepevt.h'
      integer i,mu,nold
      integer icoltmp(2)

      nold=nup
      nup=nhep
      do i=1,nup
         istup(i)=isthep(i)
         idup(i)=idhep(i)
         do mu=1,2
            mothup(mu,i)=jmohep(mu,i)
         enddo
         do mu=1,5
            pup(mu,i)=phep(mu,i)
         enddo
         if (idup(i).eq.22) then
             if (icolup(1,i).ne.0.or.icolup(2,i).ne.0) then
                 icoltmp=icolup(:,i)
                 icolup(:,i)=(/0,0/)
             endif
         endif
      enddo
      do i=1,nup
         if (idup(i).eq.21.or.abs(idup(i)).lt.6) then
             if (icolup(1,i).eq.0.and.icolup(2,i).eq.0) then
                 icolup(:,i)=icoltmp
             endif
         endif
      enddo
      do i=nold,nup
         vtimup(i)=0d0
         spinup(i)=9d0
      enddo
      end

*
**
*
      subroutine lhuptophhepevt(n)
      implicit none
      include 'LesHouches.h'
      integer ihep,mu,n

      integer nmxhep
      parameter (nmxhep=10000)
      integer idhep,isthep,jdahep,jmohep,nevhep,nhep
      real*8 phep,vhep
      common/ph_hepevt/nevhep,nhep,isthep(nmxhep),idhep(nmxhep),
     &jmohep(2,nmxhep),jdahep(2,nmxhep),phep(5,nmxhep),vhep(4,nmxhep)

      logical qedrad
      common/ph_phoqed/qedrad(nmxhep)

      nhep=nup
      nevhep=n
      do ihep=1,nhep
         isthep(ihep)=istup(ihep)
         idhep(ihep)=idup(ihep)
         do mu=1,2
            jmohep(mu,ihep)=mothup(mu,ihep)
            jdahep(mu,ihep)=0
         enddo
         do mu=1,5
            phep(mu,ihep)=pup(mu,ihep)
         enddo
         qedrad(ihep)=.true.
      enddo
      do ihep=1,nhep
         do mu=1,nhep
            if (jmohep(1,mu).eq.ihep) then
                if(jdahep(1,ihep).eq.0)then
                    jdahep(1,ihep)=mu
                else
                    jdahep(2,ihep)=mu
                endif
            elseif (jmohep(2,mu).eq.ihep) then
                if(jdahep(1,ihep).eq.0)then
                    jdahep(1,ihep)=mu
                else
                    jdahep(2,ihep)=mu
                endif
            endif
         enddo
      enddo
      end
*
** convert from ph_hepevt to hepevt
*
      subroutine translate
      implicit none
      include 'hepevt.h'
      integer ihep,mu,n

      integer ph_nmxhep
      parameter (ph_nmxhep=10000)
      integer ph_idhep,ph_isthep,ph_jdahep,ph_jmohep,ph_nevhep,ph_nhep
      real*8 ph_phep,ph_vhep
      common/ph_hepevt/ph_nevhep,ph_nhep,ph_isthep(ph_nmxhep),
     &ph_idhep(ph_nmxhep),ph_jmohep(2,ph_nmxhep),ph_jdahep(2,ph_nmxhep),
     &ph_phep(5,ph_nmxhep),ph_vhep(4,ph_nmxhep)

      nhep=ph_nhep
      do ihep=1,nhep
         isthep(ihep)=ph_isthep(ihep)
         idhep(ihep)=ph_idhep(ihep)
         do mu=1,2
            jmohep(mu,ihep)=ph_jmohep(mu,ihep)
            jdahep(mu,ihep)=ph_jdahep(mu,ihep)
         enddo
         do mu=1,5
            phep(mu,ihep)=ph_phep(mu,ihep)
         enddo
      enddo
      end

*
**
*
      subroutine find_w(id)
      implicit none
      include 'LesHouches.h'
      integer id

      do id=0,nup
          if(abs(idup(id)).eq.24) return
      enddo

      end
*
**
*
      subroutine seteps
      implicit none

      include "LesHouches.h"
      logical interf,isec,itre,iexp,iftop,ifw
      real*8 fint,fsec,expeps
      common /phokey/ fsec,fint,expeps,interf,isec,itre,iexp,iftop,ifw

      integer nmxhep
      parameter (nmxhep=10000)
      integer idhep,isthep,jdahep,jmohep,nevhep,nhep
      real*8 phep,vhep
      common/ph_hepevt/nevhep,nhep,isthep(nmxhep),idhep(nmxhep),
     &jmohep(2,nmxhep),jdahep(2,nmxhep),phep(5,nmxhep),vhep(4,nmxhep)

      real *8 kt2minqed
      common/showerqed/kt2minqed

      real*8 alpha,xphcut
      common/phocop/alpha,xphcut

      integer ihep

      integer vdecaytemp
      logical ini
      data ini/.true./
      save ini

      if (ini) then
          ini=.false.
          vdecaytemp=lprup(1)-10000
      endif
      
c find lepton
      do ihep=1,nhep
          if (idhep(ihep).eq.vdecaytemp.and.isthep(ihep).eq.1) exit
      enddo

c      expeps = sqrt(kt2minqed)/phep(4,ihep)
c      print*,expeps
      xphcut = sqrt(kt2minqed)/phep(4,ihep)

      end
*
**
*
      subroutine confronta_pt
      implicit none
      include "LesHouches.h"
      integer nmxhep
      parameter (nmxhep=10000)
      integer idhep,isthep,jdahep,jmohep,nevhep,nhep
      real*8 phep,vhep
      common/ph_hepevt/nevhep,nhep,isthep(nmxhep),idhep(nmxhep),
     &jmohep(2,nmxhep),jdahep(2,nmxhep),phep(5,nmxhep),vhep(4,nmxhep)

      real *8 kt2minqed
      common/showerqed/kt2minqed

      integer jmo,jphot,ihep

      logical ok
      common/ok_photos/ok
      
      real*8 pt

c find mother of photons
      do ihep=1,nhep
          if (idhep(ihep).eq.22.and.isthep(ihep).eq.1) then
              jmo = jmohep(1,ihep)
              pt = sqrt ( 1d0 - ( (phep(1,ihep)*phep(1,jmo) + 
     &                             phep(2,ihep)*phep(2,jmo) +
     &                             phep(3,ihep)*phep(3,jmo) )/
     &                            (abs(phep(4,ihep)
     &      *sqrt( phep(1,jmo)**2 + phep(2,jmo)**2 + phep(3,jmo)**2  )))
     &                          )**2  ) * phep(4,ihep)

              if (pt**2.lt.kt2minqed) return
              if (pt**2.gt.scalup) return
          endif
      enddo

      ok = .true.

      end
