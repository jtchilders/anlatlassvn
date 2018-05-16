      subroutine setup_PYTHIA_tune
      implicit none
      real * 8 powheginput
      external powheginput
      integer pythiatune
C --------------------------------------------------------------- C
C - N.B. PYTUNE(ITUNE) must be called before the call to PYINIT - C
C --------------------------------------------------------------- C
C -  100         A :  Rick Field's CDF Tune A     (Oct 2002)
C -  103        DW :  Rick Field's CDF Tune DW    (Apr 2006)
C -  320 Perugia 0 :  Perugia update of S0-Pro    (Feb 2009)      
      
      pythiatune=powheginput('#pythiatune')
      if (pythiatune.lt.0) then
         pythiatune=340
      endif         
      call PYTUNE(pythiatune)
      end

      subroutine setup_PYTHIA_parameters
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      double precision parp,pari
      integer mstp,msti
      common/pypars/mstp(200),parp(200),msti(200),pari(200)
      integer MSTU,MSTJ
      double precision PARU,PARJ
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      integer MDCY,MDME,KFDP
      double precision brat
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
c      integer MRPY
c      double precision RRPY
c      COMMON/PYDATR/MRPY(6),RRPY(100)
      integer pycomp
      external pycomp
      integer maxev
      common/mcmaxev/maxev
      real * 8 scalupfact
      common/cscalupfac/scalupfact
      real * 8 powheginput
      external powheginput
c      character * 20 pwgprefix
c      integer lprefix
c      common/cpwgprefix/pwgprefix,lprefix
c      integer iseed,ios

c$$$      if(lprefix.eq.7.and.pwgprefix(1:3).eq.'pwg') then
c$$$         read(pwgprefix(4:7),fmt='(i4)',iostat=ios) iseed
c$$$         if(ios.eq.0) then
c$$$c initialize random seed
c$$$            write(*,*) ' *****************************************'
c$$$            write(*,*) ' initializing PYTHIA random seed to ',iseed
c$$$            write(*,*) ' *****************************************'
c$$$            mrpy(1)=iseed+10
c$$$            mrpy(2)=0
c$$$         endif
c$$$      endif
c$$$

c      mstj(41) =3 ! Photon radiation off leptons. Not recommended to touch
c                  ! this -- causes interference with the pT ordered shower.
c      mstp(61) =0 ! No IS shower
c      mstp(71) =0 ! No FS shower
c      mstp(81) =0 ! No Multiple interactions (MI increases execution time).
c      mstp(91) =0 ! No Primordial kt
c      mstp(131)=0 ! No Pile Up
c      mstp(111)=0 ! No hadronization
      
c     N.B.
c     ====
c     For the case of jet production the following parameter setting
c     limits the transverse momentum of secondary scatterings, due
c     to multiple parton interactions, to be less than that of the
c     primary interaction (see POWHEG Dijet paper arXiv:1012.3380
c     [hep-ph] sec. 4.1 and also the PYTHIA Manual).
      mstp(86)=1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     maximum number of errors before pythia aborts (def=10)
      mstu(22)=1000
c     number of warnings printed on the shell
      mstu(26)=20
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(powheginput("#nohad").gt.0) then
         write(*,*)' ****** switching off hadronization ****'
c     Hadronization off
         mstp(111)=0
c     primordial kt off
         mstp(91)=0
c     No multiple parton interactions
         if(mstp(81).eq.1) then
c     Q2 ordered shower
            mstp(81)=0
         elseif(mstp(81).eq.21) then
c     p_t^2 ordered shower
            mstp(81)=20
         endif
      endif
      end

      subroutine getmaxev(maxev)
      integer maxev
C--- Opens input file and counts number of events, setting MAXEV;
      call opencount(maxev)
      end

      subroutine UPINIT
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      double precision parp,pari
      integer mstp,msti
      common/pypars/mstp(200),parp(200),msti(200),pari(200)
      integer MSTU,MSTJ
      double precision PARU,PARJ
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      integer MDCY,MDME,KFDP
      double precision brat
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      integer pycomp
      external pycomp
      integer hdecaymode,i   
      integer maxev
      common/mcmaxev/maxev
      real * 8 powheginput
      external powheginput
      nevhep=0
c read the header first, so lprup is set
      call lhefreadhdr(97)

c     Make PI0 stable as in herwig default
      mdcy(pycomp(111),1)=0
c     Make tau stable
      mdcy(pycomp(15),1)=0 ! tau stable. Done by Tauola too

C---Decay;
      hdecaymode=powheginput('#hdecaymode')
      if(hdecaymode.lt.-10) hdecaymode= -1
      if ((hdecaymode.lt.-1).or.(hdecaymode.gt.12)) then
         write(*,*) "Higgs decay mode not allowed"
         call pwhg_exit(-1)
      endif   

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
      end

      subroutine UPEVNT
      implicit none
      real * 8 powheginput
      external powheginput      
      logical changescalup,ini
      save changescalup,ini
      data ini/.true./
      if (ini) then
         changescalup = powheginput("#changescalup").eq.1d0 
         if (changescalup) then
            write(*,*) '*********************************'
            write(*,*) ' CHANGE SCALUP activated '
            write(*,*) '*********************************'
         endif
         ini=.false.
      endif
      call lhefreadev(97)
      if (changescalup) then
         call change_scalup
      endif
      end

      subroutine upveto
c pythia routine to abort event
      end

      subroutine pyabeg
      call init_hist
      end

      subroutine pyaend
      implicit none
      include 'pwhg_rnd.h'
       character * 100 filename
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      if(rnd_cwhichseed.eq.'none') then
         filename=pwgprefix(1:lprefix)//
     1        'POWHEG+PYTHIA-output'
      else
         filename=pwgprefix(1:lprefix)//'-'//
     1        rnd_cwhichseed //'-'//
     2        'POWHEG+PYTHIA-output'
      endif
      call pwhgsetout
      call pwhgtopout(filename)
      close(99)
      end


      subroutine pyanal
      implicit none
      include 'LesHouches.h'
      include 'hepevt.h'
      integer MDCY,MDME,KFDP
      double precision brat
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      integer pycomp
      external pycomp
      double precision bratio
      integer hdecaymode
      integer mint
      double precision vint
      COMMON/PYINT1/MINT(400),VINT(400)
c     check parameters
      logical verbose
      parameter (verbose=.false.)
      real * 8 powheginput
      external powheginput
      logical ini
      save ini,hdecaymode
      data ini/.true./
      if (ini) then
         hdecaymode=powheginput('#hdecaymode')
         ini=.false.
      endif

c     non so a che servono le righe seguenti!! C.O.
      if(mint(51).ne.0) then
         if(verbose) then
            write(*,*) 'Killed event'
            write(*,*) 'Scalup= ',scalup
            call pylist(7)      !hepeup
            call pylist(2)      !all the event
         endif
         return
      endif

      nevhep=nevhep+1
      
      if(hdecaymode.lt.-10) hdecaymode=-1
      if ((hdecaymode.eq.0).or.(hdecaymode.eq.-1)) then
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
      else
         bratio=brat(209+hdecaymode)
      endif
      if(abs(idwtup).eq.3) xwgtup=xwgtup*xsecup(1)
      xwgtup=xwgtup*bratio
      call analysis(xwgtup)
      call pwhgaccumup 
      end


      function dotp(p1,p2)
      implicit none
      real * 8 dotp,p1(0:3),p2(0:3)
      dotp = (p1(0)*p2(0) - p1(3)*p2(3)) - p1(1)*p2(1) - p1(2)*p2(2)
      end


      subroutine change_scalup
      implicit none
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      real * 8 ptmin,ptmin2,pcm(0:3,maxnup),beta,vec(3),ptkj2
      integer k,mu,j
      real * 8 scalupfact,powheginput
      external powheginput
      logical ini
      data ini/.true./
      save ini,scalupfact
      real * 8 dotp
      external dotp
      integer npart

      if (ini) then
         write(*,*) '*************************************'
         write(*,*) 'scalup set to the min pt in the event'
         write(*,*) '*************************************'
         scalupfact=powheginput('#scalupfact')         
         if(scalupfact.gt.0) then
            write(*,*) '****************************************'
            write(*,*)'ALL scalup values rescaled by the factor ',
     $           scalupfact
            write(*,*)'*****************************************'
         endif            
         ini=.false.
      endif         

      npart=0
      do k=3,nup
c     only light partons
         if (idup(k).eq.21.or.abs(idup(k)).le.5) then
            npart=npart+1
            do mu=1,3
               pcm(mu,npart)=pup(mu,k)
            enddo
            pcm(0,npart)=pup(4,k)
         endif
      enddo
c     compute min pt of light partons with respect to the incoming beam
      ptmin2=1d30
      do k=1,npart
         ptmin2=min(ptmin2,pcm(1,k)**2+pcm(2,k)**2)
      enddo

c     compute pt's of the final state partons with respect to each other
      beta=-(pup(3,1)+pup(3,2))/(pup(4,1)+pup(4,2))
      vec(1)=0
      vec(2)=0
      vec(3)=1
c     go in the CM frame   
      call mboost(npart,vec,beta,pcm,pcm)
      do k=1,npart-1
         do j=k+1,npart
            ptkj2 = 2*dotp(pcm(0,k),pcm(0,j))*
     $           pcm(0,k)*pcm(0,j)/(pcm(0,k)+pcm(0,j))**2
            ptmin2=min(ptmin2,ptkj2)
         enddo
      enddo
      ptmin=sqrt(ptmin2)
      if(scalup.gt.ptmin) then
         scalup = ptmin
      endif

      if(scalupfact.gt.0) scalup=scalup*scalupfact

      end



