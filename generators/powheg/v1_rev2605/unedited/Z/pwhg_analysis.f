c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include '../pwhg_book.h'
      integer diag
      real * 8 binsize(100)
      common/pwhghistcommon/binsize

      call pwhginihist

c     total cross section sanity check
      diag=1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'total','LOG',binsize(diag),0d0,3d0)

      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'mt(Z)','LOG',binsize(diag),50d0,100d0)

      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'m(Z)','LOG',binsize(diag),60d0,120d0)

      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'pt(l+)','LOG',binsize(diag),0d0
     $     ,100d0)

      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'eta(l+)','LOG',binsize(diag),-3d0,3d0)

      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'pt(Z)','LOG',binsize(diag),0d0,100d0)

      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'y(Z)','LOG',binsize(diag),-3d0,3d0)

      end

      
     
      subroutine analysis(dsig)
      implicit none
      real * 8 dsig
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include  'LesHouches.h'
      real *8 p_lminus(0:3),p_lplus(0:3),pcm(0:3),p_ll(0:3)
      real *8 pt_lplus,pt_lminus,eta_lplus,eta_lminus,
     $delphi,mt_v,mv,ptv,yv
      integer ihep,mu
      logical ini
      data ini/.true./
      save ini
c     binsize
      integer diag
      real * 8 binsize(100)
      common/pwhghistcommon/binsize
c     we need to tell to this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer vdecaytemp
      real *8 Zmass,Zwidth,Zmass2low,Zmass2high
      
      Zmass = 91.1876d0
      Zwidth = 2.4952d0
      Zmass2low = (Zmass-40*Zwidth)**2
      Zmass2high = (Zmass+40*Zwidth)**2


      if (ini) then
         write (*,*)
         write (*,*) '********************************************'
         if(whcprg.eq.'NLO') then
            write(*,*) '       NLO analysis'
         elseif(WHCPRG.eq.'LHE   ') then
            write(*,*) '       LHE analysis'
         elseif(WHCPRG.eq.'HERWIG') then
            write (*,*) '           HERWIG ANALYSIS            '
            write(*,*) 'not implemented analysis'
            write(*,*) 'no plots will be present at the end of the run'
         elseif(WHCPRG.eq.'PYTHIA') then
            write (*,*) '           PYTHIA ANALYSIS            '
            write(*,*) 'not implemented analysis'
            write(*,*) 'no plots will be present at the end of the run'
         endif
         write(*,*) '*****************************'
         vdecaytemp=lprup(1)-10000 ! Z decay product, with positive id
         if(vdecaytemp.eq.11.or.vdecaytemp.eq.13
     $        .or.vdecaytemp.eq.15) then
            continue
         else
            write(*,*) '**************************************'
            write(*,*) ' template analysis works only for e, mu and tau'
            write(*,*) '                 STOP     '
            write(*,*) '**************************************'
            call exit(1)
         endif
         ini=.false.
      endif

      diag=0

      do mu=0,3
         pcm(mu)=0d0
      enddo


      if((WHCPRG.eq.'NLO   ').or.(WHCPRG.eq.'LHE   ')) then
c     find Z decay products
         do ihep=1,nhep
            if(isthep(ihep).eq.1) then
               if(idhep(ihep).eq.vdecaytemp) then
                  p_lminus(0)=phep(4,ihep)
                  do mu=1,3
                     p_lminus(mu)=phep(mu,ihep)
                  enddo
               elseif(idhep(ihep).eq.-vdecaytemp) then
                  p_lplus(0)=phep(4,ihep)
                  do mu=1,3
                     p_lplus(mu)=phep(mu,ihep)
                  enddo
               endif
            endif
         enddo
      else
         return
      endif

c     lminus transverse momentum
      pt_lminus=sqrt(p_lminus(1)**2 + p_lminus(2)**2)
      call get_pseudorap(p_lminus,eta_lminus)
c     lplus transverse momentum
      pt_lplus=sqrt(p_lplus(1)**2 + p_lplus(2)**2)
      call get_pseudorap(p_lplus,eta_lplus)
c     invariant mass of the charged lepton system
      do mu=0,3
         p_ll(mu)=p_lplus(mu)+p_lminus(mu)
      enddo
      call getinvmass(p_ll,mv)
c     azimuthal separation between leptons
      delphi = dabs(atan2(p_lplus(2),p_lplus(1)) - 
     $     atan2(p_lminus(2),p_lminus(1)))
      delphi=min(delphi,2d0*pi-delphi)
c     transverse mass of the charged lepton system
      mt_v=sqrt(2*pt_lplus*pt_lminus*(1d0-dcos(delphi)))
c     rapidity of the charged lepton system
      call getrapidity(p_ll,yv)
      ptv=sqrt((p_lplus(1)+p_lminus(1))**2
     $     + (p_lplus(2)+p_lminus(2))**2)




c     total sigma
      diag=diag+1
      call pwhgfill(diag,1.5d0,dsig/binsize(diag))

c     transverse mass of the charged lepton system
      diag=diag+1
      call pwhgfill(diag,mt_v,dsig/binsize(diag))

c     invariant mass of the charged lepton system
      diag=diag+1
      call pwhgfill(diag,mv,dsig/binsize(diag))

c     pt(l+)
      diag=diag+1
      call pwhgfill(diag,pt_lplus,dsig/binsize(diag))

c     eta(l+)
      diag=diag+1
      call pwhgfill(diag,eta_lplus,dsig/binsize(diag))

c     pt(v)
      diag=diag+1
      call pwhgfill(diag,ptv,dsig/binsize(diag))

c     y(v)
      diag=diag+1
      call pwhgfill(diag,yv,dsig/binsize(diag))
         
      end
      

      subroutine getrapidity(p,y)
      implicit none
      real * 8 p(0:3),y
      y=0.5d0*log((p(0)+p(3))/(p(0)-p(3)))
      end

      subroutine getinvmass(p,m)
      implicit none
      real * 8 p(0:3),m
      m=sqrt(abs(p(0)**2-p(1)**2-p(2)**2-p(3)**2))
      end

      subroutine get_pseudorap(p,eta)
      implicit none
      real*8 p(0:3),eta,pt,th
      real *8 tiny
      parameter (tiny=1.d-5)

      pt=sqrt(p(1)**2+p(2)**2)
      if(pt.lt.tiny.and.abs(p(3)).lt.tiny)then
         eta=sign(1.d0,p(3))*1.d8
      elseif(pt.lt.tiny) then   !: added this elseif
         eta=sign(1.d0,p(3))*1.d8
      else
         th=atan2(pt,p(3))
         eta=-log(tan(th/2.d0))
      endif
      end
