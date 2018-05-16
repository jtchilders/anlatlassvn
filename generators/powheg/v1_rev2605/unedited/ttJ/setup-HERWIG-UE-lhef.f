      subroutine setup_HERWIG_parameters
      include 'HERWIG65.INC'
      include 'LesHouches.h'
      logical uevent 
      parameter (uevent=.true.)
C---USER CAN RESET PARAMETERS AT
C   THIS POINT, OTHERWISE DEFAULT
C   VALUES IN HWIGIN WILL BE USED.
c      PTMIN=100.
      MAXPR=2

      ptrms=2.5d0
      write(*,*)
      write(*,*) '*******************************************'
      write(*,*) '*******************************************'
      write(*,*) ' INITIAL p_T SPREADING OF ',ptrms,' GEV    '
      write(*,*) '*******************************************'
      write(*,*) '*******************************************'
      write(*,*)
      if(.not.uevent) then
         PRSOF=0d0
         write(*,*)
         write(*,*) '*******************************************'
         write(*,*) '*******************************************'
         write(*,*) ' NO UNDERLYING EVENT  WILL BE GENERATED    '
         write(*,*) '*******************************************'
         write(*,*) '*******************************************'
         write(*,*)   
      else
         write(*,*)
         write(*,*) '*******************************************'
         write(*,*) '*******************************************'
         write(*,*) ' UNDERLYING EVENT  WILL BE GENERATED     '
         write(*,*) '*******************************************'
         write(*,*) '*******************************************'
         write(*,*)   
      endif
c----DO NOT USE SOFT ME CORRECTION     
      SOFTME=.FALSE.
      end

      subroutine setup_HERWIG_decay
c     Needed to force t->e ve b decay.
      CALL HWMODK(6,1D0,100,12,-11,5,0,0)
      CALL HWMODK(-6,1D0,100,-12,11,-5,0,0)
c      call HWIODK(42,3,0)  !to see the HERWIG decay table
      end

      subroutine getmaxev(maxev)
      integer maxev
C---  Opens input file and counts number of events, setting MAXEV
      call opencount(maxev)
      end
      
      subroutine UPINIT
      include 'HERWIG65.INC'
      include 'LesHouches.h'
      double precision R_jet,ptmin_jet,powheginput
      external powheginput
      common/cjetdefs/R_jet,ptmin_jet
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   JET DEFINITIONS MANDATORY TO DEFINE THE PROCESS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      R_jet=powheginput('R_jet')
      ptmin_jet=powheginput('ptmin_jet')

      nevhep=0
c first call lhefreadhdr; this sets lprup;
      call lhefreadhdr(97)
C---CALL HWUSTA TO MAKE ANY PARTICLE STABLE
      CALL HWUSTA('PI0     ')
      CALL HWUSTA('TAU-    ')
      CALL HWUSTA('TAU+    ')
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CALL HWUSTA('B+      ')
      CALL HWUSTA('B-      ')
      CALL HWUSTA('B_D0    ')
      CALL HWUSTA('B_DBAR0 ')
      CALL HWUSTA('B_S0    ')
      CALL HWUSTA('B_SBAR0 ')
      CALL HWUSTA('SIGMA_B+')
      CALL HWUSTA('SIGMA_B-')
      CALL HWUSTA('XI_B0   ')
      CALL HWUSTA('XI_B+   ')
      CALL HWUSTA('XI_B-   ')
      CALL HWUSTA('B_C+    ')
      CALL HWUSTA('B_C-    ')
      CALL HWUSTA('UPSLON1S')
      CALL HWUSTA('SGM_BBR+')
      CALL HWUSTA('SGM_BBR-')
      CALL HWUSTA('LMD_BBR0')
      CALL HWUSTA('OMEGA_B-')
      CALL HWUSTA('XI_BBAR0')
      CALL HWUSTA('OMG_BBR+')
      CALL HWUSTA('LMBDA_B0')



c      call exit(1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end
      
      subroutine colour_conj(icol)
      implicit none
      integer icol(2)
      integer itmp
      itmp=icol(2)
      icol(2)=icol(1)
      icol(1)=itmp
      end

      
c$$$      subroutine displeshouches
c$$$      include 'LesHouches.h'
c$$$      integer mark(200),icol(200),j,k,ilist,icur
c$$$      write(*,*) 'incoming beams:', idbmup(1), idbmup(2)
c$$$      write(*,*) 'number of partons in subprocess',nup
c$$$      write(*,'(a,10(i3,1x))') 'ist:',(istup(j),j=1,nup)
c$$$      write(*,'(a,10(i3,1x))') 'ids:',(idup(j),j=1,nup)
c$$$      do j=1,nup
c$$$         write(*,'(5(d10.4,1x))') (pup(k,j),k=1,5)
c$$$      enddo
c$$$c check that quarks have zero anticolour, and antiquarks
c$$$c have zero colour
c$$$      do j=1,nup
c$$$         if(idup(j).ge.1.and.idup(j).le.6.and.icolup(2,j).ne.0) then
c$$$            write(*,*) ' error: anticolor in quark'
c$$$         elseif(idup(j).ge.-6.and.idup(j).le.-1.and.icolup(1,j).ne.0)
c$$$     #           then
c$$$            write(*,*) ' error: color in quark'
c$$$         endif
c$$$      enddo
c$$$c Print color linked lists
c$$$c first conjugate incoming colors
c$$$      call colour_conj(icolup(1,1))
c$$$      call colour_conj(icolup(1,2))
c$$$      do j=1,nup
c$$$         mark(j)=1
c$$$      enddo
c$$$      imarked=0
c$$$      do j=1,nup
c$$$         if(istup(j).eq.2.or.
c$$$     #        (icolup(1,j).eq.0.and.icolup(2,j).eq.0)) then
c$$$            mark(j)=0
c$$$            imarked=imarked+1
c$$$         endif
c$$$         if(istup(j).ne.2.and.idup(j).gt.0.and.idup(j).le.6) then
c$$$            icur=j
c$$$         endif
c$$$      enddo
c$$$      mark(icur)=0
c$$$      ilist=1
c$$$      imarked=imarked+1
c$$$      icol(ilist)=icur
c$$$ 12   continue
c$$$      do j=1,nup
c$$$         if(istup(j).ne.2.and.
c$$$     #        (icolup(1,j).ne.0.or.icolup(2,j).ne.0)) then
c$$$         if(mark(j).ne.0.and.icolup(1,icur).eq.icolup(2,j)) then
c$$$            ilist=ilist+1
c$$$            imarked=imarked+1
c$$$            icol(ilist)=j
c$$$            mark(j)=0
c$$$            icur=j
c$$$            if(imarked.eq.nup) goto 22
c$$$            if(icolup(1,icur).eq.0) then
c$$$               do k=1,nup
c$$$                  if(mark(k).ne.0.and.icolup(2,k).eq.0) then
c$$$                     ilist=ilist+2
c$$$                     icol(ilist)=k
c$$$                     icol(ilist-1)=0
c$$$                     imarked=imarked+1
c$$$                     icur=k
c$$$                     mark(k)=0
c$$$                     goto 12
c$$$                  endif
c$$$               enddo
c$$$               write(*,*) ' inconsistent colors!'
c$$$               stop
c$$$            endif
c$$$            goto 12
c$$$         endif
c$$$         endif
c$$$      enddo
c$$$      write(*,*) ' inconsistent colors!'
c$$$      stop
c$$$ 22   continue
c$$$      call colour_conj(icolup(1,1))
c$$$      call colour_conj(icolup(1,2))
c$$$      write(*,*) (icol(j),j=1,ilist)
c$$$      end


      subroutine UPEVNT
      implicit none
      include 'LesHouches.h'
      real * 8 scalupfac,powheginput
      external powheginput
      logical ini
      data ini/.true./
      save ini
      call lhefreadev(97)
c      call displeshouches
      scalupfac=powheginput('#scalupfac')
      if (scalupfac.lt.0) scalupfac=1d0
      if ((scalupfac.ne.1d0).and.ini) print *,"######  SCALUPFAC ="
     $     ,scalupfac
      ini=.false.
      scalup=scalup*scalupfac
      end


      subroutine hwabeg
      call init_hist
      end

      subroutine hwaend
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      open(unit=99,file=pwgprefix(1:lprefix)//'POWHEG+HERWIG-UE-output.top'
     #     ,status='unknown')
      call pwhgsetout
      call pwhgtopout
      close(99)
      end
      

      subroutine hwanal
      INCLUDE 'HERWIG65.INC'
      include 'LesHouches.h'
c     check parameters
      logical verbose
      parameter (verbose=.false.)

      if (ierror.ne.0) then
         if(verbose) then
            write(*,*) 'Killed event'
            write(*,*) 'Scalup= ',scalup
            call HWUPUP         !hepeup
            call hwuepr         !all the event
         endif
         return
      endif
      if(abs(idwtup).eq.3) xwgtup=xwgtup*xsecup(1)
      call analysis(xwgtup)
      call pwhgaccumup 
      end

