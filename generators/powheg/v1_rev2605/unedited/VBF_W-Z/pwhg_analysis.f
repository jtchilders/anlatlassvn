c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none
      include '../pwhg_book.h'
      include 'PhysPars.h'
      include 'PhysPars_Higgs.h'
!       include 'pwhg_math.h'

      include 'process.inc'
      real * 8 pi,pi2
      parameter(pi = 3.141592653589793D0, pi2 = 9.869604401089358D0)
      real * 8 ptvbcut
      common/cptvbcut/ptvbcut
      character * 10 cut
      integer nsigma,diag
      real * 8 step,invmasslow,invmasshigh,ymax
      real * 8 binsize(100)
      common/pwhghistcommon/binsize
      logical ini
      data ini/.true./
      save ini
      double precision vbfnloinput
      external vbfnloinput
      include 'global.inc'      
      
       procID = vbfnloinput("#PROC_ID")     
      select case (procID)
      case (zjj_l)
      if (ini) then
         write(*,*) '********************************************'
         write(*,*) '********************************************'
         write(*,*) 'ph_Zmass = ',ph_Zmass
         write(*,*) 'ph_Zwidth = ',ph_Zwidth
         write(*,*) '********************************************'
         write(*,*) '********************************************'
         ini=.false.
      endif

      
      nsigma = 30
      invmasslow =0!ph_Zmass-nsigma*ph_Zwidth
      invmasshigh=250!ph_Zmass+nsigma*ph_Zwidth
      step = 2*nsigma*ph_Zwidth/50
      cut = ' WBF cuts '
      ymax = 7.2d0

      call pwhginihist
      diag=1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'abs(yj(1) - yj(2))'//cut,'LIN',
     #     binsize(diag),0d0,10d0)

      diag=diag+1
      binsize(diag) = 0.5d0
      call pwhgbookup(diag,'yj(1) - yj(2)'//cut,'LIN',
     #     binsize(diag),-10d0,10d0)
      diag = diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'sig(all VBF cuts)','LIN',
     &                binsize(diag),0d0,1d0)
      diag=diag+1
      binsize(diag) = 70d0
      call pwhgbookup(diag,'mjj all'//cut,'LIN',binsize(diag),0d0,3500d0)

      diag=diag+1
      binsize(diag) = 2d0
      call pwhgbookup(diag,'ptj(3)'//cut,'LOG',binsize(diag),0d0,100d0)

      diag=diag+1
      binsize(diag) = 5d0
      call pwhgbookup(diag,'ptj(2)'//cut,'LOG',binsize(diag),0d0,400d0)

      diag=diag+1
      binsize(diag) = 10d0
      call pwhgbookup(diag,'ptj(1)'//cut,'LOG',binsize(diag),0d0,500d0)


      diag=diag+1
      binsize(diag) = 10d0
      call pwhgbookup(diag,'ptl(2)'//cut,'LOG',binsize(diag),0d0,400d0)

      diag=diag+1
      binsize(diag) = 10d0
      call pwhgbookup(diag,'ptl(1)'//cut,'LOG',binsize(diag),0d0,400d0)

      diag=diag+1
      binsize(diag) = 0.25d0
      call pwhgbookup(diag,'yj(3)'//cut,'LIN',binsize(diag),-5d0,5d0)
 

      diag=diag+1
      binsize(diag) = 0.25d0
      call pwhgbookup(diag,'yj(2)'//cut,'LIN',binsize(diag),-5d0,5d0)


      diag=diag+1
      binsize(diag) = 0.25d0
      call pwhgbookup(diag,'yj(1)'//cut,'LIN',binsize(diag),-5d0,5d0)

      diag=diag+1
      binsize(diag) = 0.25d0
      call pwhgbookup(diag,'yj(3)-0.5*(yj(1)+yj(2))'//cut,'LIN',
     #     binsize(diag),-5d0,5d0) 

      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'min(abs(yj(1)),abs(yj(2)))'//cut,'LIN',
     #     binsize(diag),0d0,ymax)


      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'max(abs(yj(1)),abs(yj(2)))'//cut,'LIN',
     #     binsize(diag),0d0,ymax)

      diag=diag+1
      binsize(diag) = 3.6d0
      call pwhgbookup(diag,'delphi_jj'//cut,'LIN',
     1  binsize(diag),-180d0,180d0)

      diag=diag+1
      binsize(diag) = 3.6d0
      call pwhgbookup(diag,'delphi_ll'//cut,'LIN',
     1  binsize(diag),-180d0,180d0)

      diag=diag+1
      binsize(diag) = 40d0
      call pwhgbookup(diag,'mZjj ','LIN',binsize(diag),0d0,3600d0)


      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'jet multip'//cut,'LOG',
     #     binsize(diag),1.5d0,8.5d0)

      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'jet mult. btw tag jets'//cut,'LOG',
     #     binsize(diag),0.5d0,8.5d0)

      
      
      case(wpjj,wmjj)

      if (ini) then
         write(*,*) '********************************************'
         write(*,*) '********************************************'
         write(*,*) 'ph_Wmass = ',ph_Wmass
         write(*,*) 'ph_Wwidth = ',ph_Wwidth
         write(*,*) '********************************************'
         write(*,*) '********************************************'
         ini=.false.
      endif      
      
      nsigma = 30
      invmasslow =0
      invmasshigh=250
      step = 5
      cut = ' WBF cuts '
      ymax = 7.2d0


      call pwhginihist
      diag=1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'abs(yj(1) - yj(2))'//cut,'LIN',
     #     binsize(diag),0d0,10d0)

      diag=diag+1
      binsize(diag) = 0.5d0
      call pwhgbookup(diag,'yj(1) - yj(2)'//cut,'LIN',
     #     binsize(diag),-10d0,10d0)

      diag = diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'sig(all VBF cuts)','LIN',
     &                binsize(diag),0d0,1d0)

      diag=diag+1
      binsize(diag) = 70d0
      call pwhgbookup(diag,'mjj all'//cut,'LIN',binsize(diag),0d0,3500d0)

      diag=diag+1
      binsize(diag) = 2d0
      call pwhgbookup(diag,'ptj(3)'//cut,'LOG',binsize(diag),0d0,100d0)

      diag=diag+1
      binsize(diag) = 5d0
      call pwhgbookup(diag,'ptj(2)'//cut,'LOG',binsize(diag),0d0,400d0)

      diag=diag+1
      binsize(diag) = 10d0
      call pwhgbookup(diag,'ptj(1)'//cut,'LOG',binsize(diag),0d0,500d0)


      diag=diag+1
      binsize(diag) = 10d0
      call pwhgbookup(diag,'ptl(2)'//cut,'LOG',binsize(diag),0d0,400d0)

      diag=diag+1
      binsize(diag) = 10d0
      call pwhgbookup(diag,'ptl(1)'//cut,'LOG',binsize(diag),0d0,400d0)

      diag=diag+1
      binsize(diag) = 0.25d0
      call pwhgbookup(diag,'yj(3)'//cut,'LIN',binsize(diag),-5d0,5d0)
 

      diag=diag+1
      binsize(diag) = 0.25d0
      call pwhgbookup(diag,'yj(2)'//cut,'LIN',binsize(diag),-5d0,5d0)


      diag=diag+1
      binsize(diag) = 0.25d0
      call pwhgbookup(diag,'yj(1)'//cut,'LIN',binsize(diag),-5d0,5d0)

      diag=diag+1
      binsize(diag) = 0.25d0
      call pwhgbookup(diag,'yj(3)-0.5*(yj(1)+yj(2))'//cut,'LIN',
     #     binsize(diag),-5d0,5d0) 

      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'min(abs(yj(1)),abs(yj(2)))'//cut,'LIN',
     #     binsize(diag),0d0,ymax)


      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'max(abs(yj(1)),abs(yj(2)))'//cut,'LIN',
     #     binsize(diag),0d0,ymax)

      diag=diag+1
      binsize(diag) = 3.6d0
      call pwhgbookup(diag,'delphi_jj'//cut,'LIN',
     1  binsize(diag),-180d0,180d0)

      diag=diag+1
      binsize(diag) = 3.6d0
      call pwhgbookup(diag,'delphi_lnu'//cut,'LIN',
     1  binsize(diag),-180d0,180d0)


      diag=diag+1
      binsize(diag) = 40d0
      call pwhgbookup(diag,'mWjj ','LIN',binsize(diag),0d0,3600d0)


      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'jet multip'//cut,'LOG',
     #     binsize(diag),1.5d0,8.5d0)

      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'jet mult. btw tag jets'//cut,'LOG',
     #     binsize(diag),0.5d0,8.5d0)


      end select


      end

      
     
      subroutine analysis(dsig)
      implicit none
      include 'nlegborn.h'      
      include 'pwhg_math.h'
      include 'vbfnlo-files/global.inc'
      include 'process.inc'
      include 'pwhg_kn.h'

      real * 8 dsig 
      
!       if(kn_jacborn.eq.0d0) return
      if(dsig.eq.0d0) return
!       print*, procID
      select case(procID)  
      
      case(Zjj_l)
         call analysisZ(dsig)
      case(Wpjj,Wmjj)
         call analysisW(dsig)      
         
      end select
      
      end
      
      subroutine analysisZ(dsig)
      implicit none      
      real * 8 dsig
      include 'hepevt.h'
      include 'pwhg_math.h'  
      include 'nlegborn.h'
      include 'pwhg_kn.h'
c arrays to reconstruct jets
      integer maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real *8 ptrack(4,maxtrack)
      integer maxtagjets
      parameter (maxtagjets=10)
      integer jet,lep,tagjets
      real *8 ptj(maxtagjets),yj(maxtagjets)
      real *8 pjet(4,maxjet) 
      integer mu,jpart,jjet,njj(maxtagjets),found,njets,
     #     ihep,ntracks,ijet,j1
      real * 8 vec(3),pjetout(0:3),beta,ptrel,get_ptrel,
     #     ptrackin(0:3),ptrackout(0:3),partons(0:3,5)
      integer i,diag,njets_passcut
      external get_ptrel
      real * 8 R,ptmin_fastkt, mllmin, Rllmin,rminjj, mll,Rll, rsep_jjmin,rsepjj
      integer jetvec(maxtrack)
      logical ini
      data ini/.true./
      save ini
      integer HZZ,HWW
      integer idH, k, kk, kmax, j, nphotons, npartons
      parameter (nphotons=2, npartons=5)
      real * 8 pH(0:3),ptH,yH,inv_mH,Edec,thl,phil,plepCM(0:3,2),
     $     plep(0:3,2),mod_vecH,pj(0:3,3),pT_lep(2),y_lep(2),pHjj(0:3),
     $     mHjj, isol_phopart, eff_phopart, y_partons(5), phi_lep(2)
      real * 8 random,rsepn,getrapidity0,mjj,azi,mjj2,rllmin2
      external random,rsepn,getrapidity0,mjj,azi,mjj2
      real * 8 Rsep(2,2),Rmin,delphi_jj, Rll_min,dely_jl
      logical pass_cuts,pass_cuts_no_mjjmin,pass_cuts_no_deltayjjmin
      real * 8 ptjetmin,ptalljetmin,yjetmax,deltay_jjmin,ptlepmin,
     #     ylepmax,mjjmin,Rsep_jlmin,ptvetojet,sig2NLO,pjN(0:3)
      logical ylep_between_jets,jet_opphem,exist3rdjet
      logical onlyquarks,Z_exchange
      real * 8 binsize(100), ptl(2)
      common/pwhghistcommon/binsize
c      logical iniptcut
c      save iniptcut
c      data iniptcut/.true./
      logical higgsfound,found_hardest_vetojet    
c     we need to tell to the this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer counter1, counter2, counter3
      data counter1/0/
      save counter1
      data counter2/0/
      save counter2
      data counter3/0/
      save counter3
      double precision RSEPS, powheginput, rsepn4
      external RSEPS, powheginput, rsepn4
      double precision z_three, pt_av(10), r_av(10)  
c     photon isolation stuff
      double precision ptsum,factor
      double precision delta(npartons)
      double precision dlist(npartons,nphotons)
      double precision ptlist(npartons,nphotons)
      double precision flist(npartons,nphotons)
      double precision spt(npartons,nphotons)
      integer counter_partons, countertot,countertot1, passcuts_count
      double precision pt_av2(10), r_av2(10), verhpt ,ptrho, rseprho
      double precision pt_av21(10), pt_av22(10), pt_avt1(10), pt_avt2(10)
      save pt_av2
      save r_av2
      data countertot/0/
      save countertot
      data countertot1/0/
      save countertot1         
!       print*, "HERE"

      
      
c     DISABLE ALL CUTS

      if (ini) then
      ptalljetmin = 0d0            
      ptjetmin = 0d0         
      yjetmax = 100d0           
      mjjmin = 0d0              
      deltay_jjmin = 0d0        

      ptlepmin = 0d0            
      ylepmax = 100d0           
      Rsep_jlmin = 0d0          
      Rsep_jjmin = 0d0   
      dely_jl= 0d0

      ylep_between_jets = .false.
      jet_opphem = .false.
      Rllmin = 0d0
      mllmin=0d0
      
      ptalljetmin=powheginput('#ptalljetmin')
      ptjetmin=powheginput('#ptjetmin')
      mjjmin=powheginput('#mjjmin')  
      deltay_jjmin=powheginput('#deltay_jjmin')
      yjetmax=powheginput('#yjetmax')     
      if (yjetmax.lt.0.1d0) yjetmax=100d0      
      ylepmax=powheginput('#ylepmax')
      if (ylepmax.lt.0.1d0) ylepmax=100d0
      ptlepmin=powheginput('#ptlepmin')
      mllmin=powheginput('#mllmin')
      Rllmin=powheginput('#Rllmin')
      Rsep_jjmin=powheginput('#Rsep_jjmin')
      Rsep_jlmin=powheginput('#Rsep_jlmin')
      dely_jl=powheginput('#dely_jl')
      if(powheginput('#ylep_between_jets').eq.1d0) ylep_between_jets=.true.
      if(powheginput('#jet_opphem').eq.1d0) jet_opphem=.true.      
      
!       ini=.false.
      endif

      

       goto 811
c     ACTIVATE CUTS!!!
      Rsep_jjmin = 0.7d0   
      ptalljetmin = 30d0
      ptjetmin = 30d0
      yjetmax = 5d0
      mjjmin = 600d0
      deltay_jjmin = 4.2d0
      jet_opphem = .true.


      
       goto 811

c     LEPTONIC CUTS
      ptlepmin = 20d0
      ylepmax = 2.5d0
      Rsep_jlmin = 0.6d0
      ylep_between_jets = .true.
      Rllmin = 0.2d0   
      mllmin=20d0

 811  continue

      if (ini) then
         write(*,*) '**************************************************'
         write(*,*) '**************************************************'
         write(*,*) '                ANALYSIS CUTS                     '
         write(*,*) '**************************************************'
         write(*,*) '**************************************************'
         write(*,*) 'Rsep_jjmin = ',Rsep_jjmin
         write(*,*) 'ptalljetmin = ',ptalljetmin
         write(*,*) 'ptjetmin = ',ptjetmin
         write(*,*) 'yjetmax = ',yjetmax
         write(*,*) 'mjjmin = ',mjjmin 
         write(*,*) 'deltay_jjmin = ',deltay_jjmin
         write(*,*) 'ptlepmin = ',ptlepmin 
         write(*,*) 'ylepmax = ',ylepmax 
         write(*,*) 'Rsep_jlmin = ',Rsep_jlmin
         write(*,*) 'ylep_between_jets = ',ylep_between_jets 
         write(*,*) 'jet_opphem = ',jet_opphem 
         write(*,*) '**************************************************'
         write(*,*) '**************************************************'
         ini = .false.
                 
      endif


      goto 143
      onlyquarks = .false.
      if (onlyquarks) then
c     select events with incoming quarks only
         if (.not.(idhep(1).ne.21.and.idhep(2).ne.21)) return
      else
c     select events with incoming gluons
         if (.not.(idhep(1).eq.21.or.idhep(2).eq.21)) return
      endif

 143  continue

      goto 144
      Z_exchange = .false.
c     select Z events
      if (nhep.eq.6) then
         if (idhep(1).eq.21) then
            if (idhep(4) + idhep(6).eq.0) then
               Z_exchange = .true.
            endif
         elseif (idhep(6).eq.21) then
            if (idhep(1).eq.idhep(4)) then
               Z_exchange = .true.
            endif 
         endif
      endif

      if (.not.Z_exchange) then
      else
         return
      endif

 144  continue

!       idH=0
      higgsfound = .false.
c     find Higgs boson
   
      do ihep=1,nhep
         if ((idhep(ihep).eq.11).or.(idhep(ihep).eq.13)) then
            do mu=1,3
               plep(mu,1)=phep(mu,ihep)
            enddo
            plep(0,1)=phep(4,ihep)
         endif
         if ((idhep(ihep).eq.-11).or.(idhep(ihep).eq.-13)) then
            do mu=1,3
               plep(mu,2)=phep(mu,ihep)          
            enddo
            plep(0,2)=phep(4,ihep)
         endif         
            do mu=0,3
               pH(mu)=plep(mu,2)+plep(mu,1)         
            enddo         
      enddo

!       if(dsig.lt.1d-20) return
!  
      do mu=0,3
         pH(mu) = plep(mu,1)+plep(mu,2)
      enddo
      
      ptH = sqrt(pH(1)**2+pH(2)**2)
      yH= getrapidity0(pH)
      call getinvmass(pH,inv_mH)
      

      
c     set up arrays for jet finding
      do jpart=1,maxtrack
         do mu=1,4
            ptrack(mu,jpart)=0d0
         enddo
         jetvec(jpart)=0
      enddo      
      do jjet=1,maxjet
         do mu=1,4
            pjet(mu,jjet)=0d0
         enddo
      enddo
      j1=0
      found=0
      ntracks=0
      njets=0
c     loop over final state particles to find jets 
      do ihep=1,nhep
         if ((isthep(ihep).eq.1).and.
c     exclude leptons, gauge and Higgs bosons
     #        (((abs(idhep(ihep)).le.10).or.(abs(idhep(ihep)).ge.40))
c     but include gluons 
     #        .or.(abs(idhep(ihep)).eq.21))) then
            if(ntracks.eq.maxtrack) then
               write(*,*)
     #              'hwanal: too many particles, increase maxtrack'
               stop
            endif
c     copy momenta to construct jets 
            ntracks=ntracks+1
            do mu=1,4
               ptrack(mu,ntracks)=phep(mu,ihep)
            enddo
         endif
      enddo      

      
      if (ntracks.eq.0) then
         return
      endif
      
************************************************************************
*     siscone algorithm
**********************************************************************
c     R = 0.7  radius parameter
c     f = 0.5  overlapping fraction
c.....run the clustering        
c      call fastjetsiscone(ptrack,ntracks,0.7d0,0.5d0,pjet,njets) 
************************************************************************
*     fastkt algorithm
**********************************************************************
c      R = 0.7  Radius parameter
c.....run the clustering 
      R = Rsep_jjmin          
      ptmin_fastkt = 0d0


      call fastjetppgenkt(ptrack,ntracks,R,-1,ptmin_fastkt,pjet,njets,
     $                        jetvec)     

      


     
c     now we have the jets
      if (njets.gt.0) then
c     find the first THREE hardest jets, if any
         call find_hardest_jets(njets,pjet,3,tagjets,njj)    
       
c     at least TWO tagging jets to continue
         if (tagjets.le.1) then
            return
         endif
         
         do ijet=1,tagjets
            do mu=1,3
               pj(mu,ijet)=pjet(mu,njj(ijet))
               partons(mu,2+ijet)=pjet(mu,njj(ijet))
            enddo
            pj(0,ijet)=pjet(4,njj(ijet))
            partons(0,2+ijet)=pjet(4,njj(ijet))
         enddo

c     get pt's and rapidities of the jets
         do ijet=1,tagjets
            ptj(ijet) = sqrt(pj(1,ijet)**2 + pj(2,ijet)**2)
            yj(ijet) = getrapidity0(pj(0,ijet))
         enddo
         do i=1,tagjets+2
             y_partons=getrapidity0(partons(0,i))
         enddo
         
         do i=1,2
            call legoy(plep(0,i), pT_lep(i),y_lep(i),phi_lep(i))
         enddo

c     compute min R_jj separations

         Rmin = 1d50
         do jet=1,tagjets-1
            do lep=jet+1,tagjets
               Rsepjj=rsepn(pj(0,jet),pj(0,lep))
               Rminjj=min(Rmin,Rsepjj)
            enddo
         enddo

c     compute min R_jlep separations
         Rmin = 1d50
         do jet=1,2
            do lep=1,2
               Rsep(jet,lep)=rsepn(pj(0,jet),plep(0,lep))
               Rmin=min(Rmin,Rsep(jet,lep))
            enddo
         enddo
         Rllmin2 = 1d50
         Rll=rsepn(plep(0,1),plep(0,2))
         rll=min(rllmin2,rll)
         
         mll=mjj(plep(0,1),plep(0,2))

         if (yj(1).gt.yj(2)) then
            delphi_jj = azi(pj(0,1))-azi(pj(0,2))
         else
            delphi_jj = azi(pj(0,2))-azi(pj(0,1))
         endif
 
c         delphi_jj =abs(azi(pj(0,1))-azi(pj(0,2)))
c      IF (DELPHI.GT.PI) THEN
c         DELPHI = 2*PI-DELPHI
c      ENDIF         
         if (delphi_jj.gt.Pi) then
            delphi_jj = delphi_jj - 2d0*pi
         else if (delphi_jj.lt.-Pi) then
            delphi_jj = delphi_jj + 2d0*pi
         endif
         delphi_jj = delphi_jj * 180d0 / Pi
 
c     compute invariant mass of the Hjj system
         do mu=0,3
            pHjj(mu) = pH(mu)+pj(mu,1)+pj(mu,2)
         enddo
         mHjj = sqrt(abs(phjj(0)**2-phjj(1)**2-phjj(2)**2-phjj(3)**2))




       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCC               APPLY CUTS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         pass_cuts = 
     $        (min(pT_lep(1),pT_lep(2)).gt.ptlepmin) .and.
     $        (max(abs(y_lep(1)),abs(y_lep(2))).lt.ylepmax) .and.
     $        (min(pTj(1),pTj(2)).gt.ptjetmin) .and.
     $        (max(abs(yj(1)),abs(yj(2))).lt.yjetmax) .and.
     $        (Rmin.gt.Rsep_jlmin).and.(rminjj.gt.rsep_jjmin)
!          print*, '1', pass_cuts
         if (ylep_between_jets) then
            pass_cuts = pass_cuts .and.
     $           ((min(yj(1),yj(2))+dely_jl).lt.y_lep(1)) .and.
     $           ((max(yj(1),yj(2))-dely_jl).gt.y_lep(1)) .and.
     $           ((min(yj(1),yj(2))+dely_jl).lt.y_lep(2)) .and.
     $           ((max(yj(1),yj(2))-dely_jl).gt.y_lep(2))
         endif
    
         pass_cuts = pass_cuts .and. (mll.gt.mllmin).and. (Rll.gt.Rllmin)
         
!          print*, '2', pass_cuts
         if (jet_opphem) then
            pass_cuts = pass_cuts .and. 
     $           (yj(1)*yj(2).lt.0)
         endif
!          print*, '3', pass_cuts
         if (tagjets.lt.3) then
            exist3rdjet = .false.
         else
            exist3rdjet = (pTj(3).gt.ptalljetmin) .and.
     $           (abs(yj(3)).lt.yjetmax)
         endif

         pass_cuts_no_mjjmin = pass_cuts .and. 
     $        (abs(yj(1)-yj(2)).gt.deltay_jjmin)
         
         pass_cuts_no_deltayjjmin = pass_cuts .and.
     $        (mjj(pj(0,1),pj(0,2)).gt.mjjmin)
         
         pass_cuts = pass_cuts .and.
     $        (mjj(pj(0,1),pj(0,2)).gt.mjjmin) .and.
     $        (abs(yj(1)-yj(2)).gt.deltay_jjmin)
!          print*, '4', pass_cuts
         
         if(exist3rdjet) then
           do lep=1,2
           Rsepjj=rsepn(pj(0,3),plep(0,lep))
           pass_cuts=pass_cuts.and.(Rsepjj.gt.Rsep_jlmin)
           enddo
         endif            

         if (pass_cuts) then    
c            call increasecnt('pass WBF cuts')

            diag=1
            call pwhgfill(diag,abs(yj(1) - yj(2)),dsig/binsize(diag))

            diag=diag+1
            call pwhgfill(diag,yj(1) - yj(2),dsig/binsize(diag))

            diag =diag+ 1   
            call pwhgfill(diag,0.5d0,dsig)

            diag=diag+1
            call pwhgfill(diag,mjj(pj(0,1),pj(0,2)),dsig/binsize(diag))

            diag=diag+1            
            if (exist3rdjet) then
               call pwhgfill(diag,ptj(3),dsig/binsize(diag))
            endif
            
            
            diag=diag+1
            call pwhgfill(diag,ptj(2),dsig/binsize(diag))

            diag=diag+1
            call pwhgfill(diag,ptj(1),dsig/binsize(diag))

            diag=diag+1
            call pwhgfill(diag,pt_lep(2),dsig/binsize(diag))

            diag=diag+1
            call pwhgfill(diag,pt_lep(1),dsig/binsize(diag))
          
c            if (iniptcut) then
c               iniptcut = .false.
c               write(*,*) '****************************'
c               write(*,*) '****************************'
c               write(*,*) 'pt cut on yj(3) in place!!  '
c               write(*,*) '****************************'
c               write(*,*) '****************************'
c            endif

            diag=diag+1
            if (exist3rdjet)  then
               call pwhgfill(diag,yj(3),dsig/binsize(diag))
            endif

            diag=diag+1
            call pwhgfill(diag,yj(2),dsig/binsize(diag))

            diag=diag+1
            call pwhgfill(diag,yj(1),dsig/binsize(diag))

            diag=diag+1
            if (exist3rdjet) then 
               call pwhgfill(diag,yj(3)-0.5*(yj(1)+yj(2)),
     $              dsig/binsize(diag))
            endif

            diag=diag+1
            call pwhgfill(diag,min(abs(yj(1)),abs(yj(2))),
     $           dsig/binsize(diag))

            diag=diag+1
            call pwhgfill(diag,max(abs(yj(1)),abs(yj(2))),
     $           dsig/binsize(diag))

            diag=diag+1
            call pwhgfill(diag,delphi_jj,dsig/binsize(diag))
!             if (toobig)  print*, 'delphi_jj  ', delphi_jj
 

         if (y_lep(1).gt.y_lep(2)) then
            delphi_jj = phi_lep(1)-phi_lep(2)
         else
            delphi_jj = phi_lep(2)-phi_lep(1)
         endif
 
c         delphi_jj =abs(azi(pj(0,1))-azi(pj(0,2)))
c      IF (DELPHI.GT.PI) THEN
c         DELPHI = 2*PI-DELPHI
c      ENDIF         
         if (delphi_jj.gt.Pi) then
            delphi_jj = delphi_jj - 2d0*pi
         else if (delphi_jj.lt.-Pi) then
            delphi_jj = delphi_jj + 2d0*pi
         endif
         delphi_jj = delphi_jj * 180d0 / Pi
c         if (abs(delphi_jj).lt.5d0) then
c             print*, plep(:,1), pt_lep(1), y_lep(1), phi_lep(1)
c             print*, plep(:,2), pt_lep(2), y_lep(2), phi_lep(2)
c         endif

            diag=diag+1
            call pwhgfill(diag,delphi_jj,dsig/binsize(diag))
!             if (toobig)  print*, 'delphi_gaga  ', delphi_jj            

            diag=diag+1
            call pwhgfill(diag,mHjj,dsig/binsize(diag))  


c     now count how many jets have pt > ptalljetmin
c     find the first 10 hardest jets, if any
            call find_hardest_jets(njets,pjet,10,tagjets,njj)
c            write(*,*) 'tagjets ',tagjets            
c     two jets surely already exists at this step
            njets_passcut = 2
            do ijet=3,tagjets
               ptj(ijet)=
     #              sqrt(pjet(1,njj(ijet))**2 + pjet(2,njj(ijet))**2)
               call getrapidity(pjet(1,njj(ijet)),yj(ijet))    
           
               if (ptj(ijet).gt.ptalljetmin .and.
     #             abs(yj(ijet)).lt.yjetmax) then
                   do mu=1,3
                      pjN(mu)=pjet(mu,njj(ijet))
                   enddo
                   pjN(0)=pjet(4,njj(ijet))
                   if(min(rsepn(pjN(0),plep(0,1)),rsepn(pjN(0),plep(0,2))).gt.Rsep_jlmin) 
     &                   njets_passcut = njets_passcut + 1
               endif
            enddo
c            write(*,*) 'njets_passcut ',njets_passcut
            diag=diag+1
            call pwhgfill(diag,njets_passcut*1d0,dsig/binsize(diag))

c     jet multeplicity between the two tagging jets
            njets_passcut = 0
            do ijet=3,tagjets
               if (min(yj(1),yj(2)).lt.yj(ijet) .and.
     #                 yj(ijet).lt.max(yj(1),yj(2)) .and.
     #                 ptj(ijet).gt.ptalljetmin) then
                   do mu=1,3
                      pjN(mu)=pjet(mu,njj(ijet))
                   enddo
                   pjN(0)=pjet(4,njj(ijet))
                   if(min(rsepn(pjN(0),plep(0,1)),rsepn(pjN(0),plep(0,2))).gt.Rsep_jlmin) 
     &                   njets_passcut = njets_passcut + 1
               endif
            enddo            
            diag=diag+1
            call pwhgfill(diag,njets_passcut*1d0,dsig/binsize(diag))
     
         
    
      
         
      endif
      endif
      end




      subroutine getrapidity(p,y)
      implicit none
      real * 8 p(4),y
      y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      end

      function getrapidity0(p)
      implicit none
      real * 8 p(0:3),getrapidity0
      getrapidity0=0.5d0*log((p(0)+p(3))/(p(0)-p(3)))
      end

      subroutine getinvmass(p,m)
      implicit none
      real * 8 p(0:3),m
      m=sqrt(abs(p(0)**2-p(1)**2-p(2)**2-p(3)**2))
      end






      function get_ptrel(pin,pjet)
      implicit none
      real * 8 get_ptrel,pin(0:3),pjet(0:3)
      real * 8 pin2,pjet2,cth2,scalprod
      pin2  = pin(1)**2 + pin(2)**2 + pin(3)**2
      pjet2 = pjet(1)**2 + pjet(2)**2 + pjet(3)**2
      scalprod = pin(1)*pjet(1) + pin(2)*pjet(2) + pin(3)*pjet(3)
      cth2 = scalprod**2/pin2/pjet2
      get_ptrel = sqrt(pin2*abs(1d0 - cth2))
      end
      

      function azi(p)
      implicit none
      include 'pwhg_math.h'  
      real * 8 azi,p(0:3)
      azi = atan(p(2)/p(1))
      if (p(1).lt.0d0) then
         if (azi.gt.0d0) then               
            azi = azi - pi
         else
            azi = azi + pi
         endif
      endif    
      end

c     calculate the separation in the lego plot between the two momenta
c     p1 and p2
      function rsepn(p1,p2)
      implicit none
      include 'pwhg_math.h'  
      real * 8 rsepn,p1(0:3),p2(0:3)
      real * 8 y1,phi1,y2,phi2
      real * 8 delphi
      real * 8 getrapidity0,azi
      external getrapidity0,azi

      phi1 = azi(p1)   
      phi2 = azi(p2)
      y1 = getrapidity0(p1)
      y2 = getrapidity0(p2)

      delphi = abs(phi1-phi2)
      if (delphi.gt.pi) then
         delphi = 2*pi-delphi
      endif
      if (delphi.lt.0 .or. delphi.gt.pi) then
         print*,' problem in rsepn. delphi = ',delphi
      endif
      rsepn = sqrt( (y1-y2)**2 + delphi**2 )
      end

c     calculate the separation in the lego plot between the two momenta
c     p1 and p2 for pj(1:4)
      function rsepn4(p1,p2)
      implicit none
      include 'pwhg_math.h'  
      real * 8 rsepn4,p1(4),p2(4)
      real * 8 y1,phi1,y2,phi2
      real * 8 delphi
      real * 8 azi4
      external azi4

      phi1 = azi4(p1)   
      phi2 = azi4(p2)
      call getrapidity(p1,y1)
      call getrapidity(p2,y2)

      delphi = abs(phi1-phi2)
      if (delphi.gt.pi) then
         delphi = 2*pi-delphi
      endif
      if (delphi.lt.0 .or. delphi.gt.pi) then
         print*,' problem in rsepn. delphi = ',delphi
      endif
      rsepn4 = sqrt( (y1-y2)**2 + delphi**2 )
      end
      
      
      
      function azi4(p)
      implicit none
      include 'pwhg_math.h'  
      real * 8 azi4,p(4)
      azi4 = atan(p(2)/p(1))
      if (p(1).lt.0d0) then
         if (azi4.gt.0d0) then               
            azi4 = azi4 - pi
         else
            azi4 = azi4 + pi
         endif
      endif    
      end


c mjj^2 = (p1+p2)^2 = p1^2 + p2^2 + 2*dotp(p1,p2)
      function mjj(p1,p2)
      implicit none
      real * 8 mjj,p1(0:3),p2(0:3)
      real * 8 p(0:3)
      integer mu
      do mu=0,3
         p(mu)=p1(mu)+p2(mu)
      enddo
      mjj = sqrt(abs(p(0)**2-p(1)**2-p(2)**2-p(3)**2))
      end



c     find the first "nhardjets" hardest jets in pjet (that contains njets)
c     and return their position.
c     foundhardjets is the number of found hard jets (.le.nhardjets)
      subroutine find_hardest_jets(njets,pjet,nhardjets,
     #     foundhardjets,jj)
      implicit none
      integer njets
      real *8 pjet(4,njets) 
      integer nhardjets,jj(nhardjets)
      real * 8 ptj(nhardjets),pt
      integer ijet,hjet,foundhardjets,i
      logical is_i_in_array
      external is_i_in_array

      if (njets.eq.0) then
         write(*,*) 'WARNING!!!!!!!!!!!  EMPTY  PJET ARRAY'
         nhardjets=0
         return
      endif

      do hjet=1,nhardjets
         jj(hjet)=0d0
         ptj(hjet)=0d0
      enddo
      foundhardjets=1
      do ijet=1,njets   
         pt=sqrt(pjet(1,ijet)**2 + pjet(2,ijet)**2)
         do hjet=1,min(foundhardjets,nhardjets)
            if (pt.gt.ptj(hjet).and.
     $           .not.is_i_in_array(nhardjets,ijet,jj)) then
               foundhardjets = foundhardjets + 1
               do i=nhardjets,hjet+1,-1
                  ptj(i)=ptj(i-1)
                  jj(i)=jj(i-1)
               enddo
               ptj(hjet)=pt
               jj(hjet)=ijet
            endif
         enddo
      enddo
c     set number of jets found
      foundhardjets = min(foundhardjets-1,nhardjets)
      end

      function is_i_in_array(nhardjets,i,jj)
      implicit none
      logical is_i_in_array
      integer nhardjets,i,jj(nhardjets)
      integer j
      is_i_in_array = .false.
      do j=1,nhardjets
         if (i.eq.jj(j)) then
            is_i_in_array = .true.
            return
         endif
      enddo
      end



! !       subroutine particle_identif(HWW,HZZ)
! !       implicit none
! !       integer pdg_Higgs,pdg_Z,pdg_W,HZZ,HWW
! !       pdg_Higgs = 25
! !       pdg_Z = 23
! !       pdg_W = 24      
! ! c     build an identifier for Higgs production in WW and ZZ fusion 
! !       HWW = 10000*pdg_W + pdg_Higgs
! !       HZZ = 10000*pdg_Z + pdg_Higgs
! !       end


      subroutine analysisW(dsig)
      implicit none
      real * 8 dsig
      include 'hepevt.h'
      include 'pwhg_math.h'  
      include 'nlegborn.h'
      include 'pwhg_kn.h'
c arrays to reconstruct jets
      integer maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real *8 ptrack(4,maxtrack)
      integer maxtagjets
      parameter (maxtagjets=10)
      integer jet,lep,tagjets
      real *8 ptj(maxtagjets),yj(maxtagjets)
      real *8 pjet(4,maxjet) 
      integer mu,jpart,jjet,njj(maxtagjets),found,njets,
     #     ihep,ntracks,ijet,j1
      real * 8 vec(3),pjetout(0:3),beta,ptrel,get_ptrel,
     #     ptrackin(0:3),ptrackout(0:3),partons(0:3,5)
      integer i,diag,njets_passcut
      external get_ptrel
      real * 8 R,ptmin_fastkt, mllmin, Rllmin,rminjj, mll,Rll, rsep_jjmin,rsepjj
      integer jetvec(maxtrack)
      logical ini
      data ini/.true./
      save ini
      integer HZZ,HWW
      integer idH, k, kk, kmax, j, nphotons, npartons
      parameter (nphotons=2, npartons=5)
      real * 8 pH(0:3),ptH,yH,inv_mH,Edec,thl,phil,plepCM(0:3,2),pjN(0:3),
     $     plep(0:3,2),mod_vecH,pj(0:3,3),pT_lep(2),y_lep(2),pHjj(0:3),
     $     mHjj, isol_phopart, eff_phopart, y_partons(5), phi_lep(2)
      real * 8 random,rsepn,getrapidity0,mjj,azi,mjj2,rsepn4, verhpt
      external random,rsepn,getrapidity0,mjj,azi,mjj2,rsepn4
      real * 8 Rsep(8,2),Rmin,delphi_jj, Rll_min,dely_jl,actcuts
      logical pass_cuts,pass_cuts_no_mjjmin,pass_cuts_no_deltayjjmin
      real * 8 ptjetmin,ptalljetmin,yjetmax,deltay_jjmin,ptlepmin,
     #     ylepmax,mjjmin,Rsep_jlmin,ptvetojet,sig2NLO
      logical ylep_between_jets,jet_opphem,exist3rdjet
      logical onlyquarks,Z_exchange,activate_cuts
      real * 8 binsize(100), ptl(2), pj4(0:3)
      common/pwhghistcommon/binsize
c      logical iniptcut
c      save iniptcut
c      data iniptcut/.true./
      logical higgsfound,found_hardest_vetojet, passcuts4j
c     we need to tell to the this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer counter1, counter2, counter3
      data counter1/0/
      save counter1
      data counter2/0/
      save counter2
      data counter3/0/
      save counter3
      double precision RSEPS, powheginput,z_three, pt_av(10), r_av(10)
      external RSEPS, powheginput
c     photon isolation stuff
      double precision ptsum,factor,rseprho, ptrho
      double precision delta(npartons)
      double precision dlist(npartons,nphotons)
      double precision ptlist(npartons,nphotons)
      double precision flist(npartons,nphotons)
      double precision spt(npartons,nphotons)
      integer counter_partons, countertot,countertot1, passcuts_count
      double precision pt_av2(10), r_av2(10)
      double precision pt_av21(10), pt_av22(10), pt_avt1(10), pt_avt2(10)
      save pt_av2
      save r_av2
      data countertot/0/
      save countertot
      data countertot1/0/
      save countertot1      

      data passcuts_count/0/
      save passcuts_count


c     DISABLE ALL CUTS

      if (ini) then
      ptalljetmin = 0d0            
      ptjetmin = 0d0         
      yjetmax = 100d0           
      mjjmin = 0d0              
      deltay_jjmin = 0d0        

      ptlepmin = 0d0            
      ylepmax = 100d0           
      Rsep_jlmin = 0d0          
      Rsep_jjmin = 0.1d0 !For jet algorithm   
      dely_jl= 0d0

      ylep_between_jets = .false.
      jet_opphem = .false.
      passcuts4j = .false.
      Rllmin = 0d0
      mllmin=0d0
      
      actcuts=powheginput('#activate_cuts')
      if (actcuts.eq.1d0) then
      ptalljetmin=powheginput('#ptalljetmin')
      ptjetmin=powheginput('#ptjetmin')
      mjjmin=powheginput('#mjjmin')  
      deltay_jjmin=powheginput('#deltay_jjmin')
      yjetmax=powheginput('#yjetmax')     
      if (yjetmax.lt.0.1d0) yjetmax=100d0      
      ylepmax=powheginput('#ylepmax')
      if (ylepmax.lt.0.1d0) ylepmax=100d0
      ptlepmin=powheginput('#ptlepmin')
      mllmin=powheginput('#mllmin')
      Rllmin=powheginput('#Rllmin')
      Rsep_jjmin=powheginput('#Rsep_jjmin')
      Rsep_jlmin=powheginput('#Rsep_jlmin')
      dely_jl=powheginput('#dely_jl')
      if(powheginput('#ylep_between_jets').eq.1d0) ylep_between_jets=.true.
      if(powheginput('#jet_opphem').eq.1d0) jet_opphem=.true.   
      if(powheginput('#passcuts4j').eq.1d0) passcuts4j=.true.   
      endif
      
!       ini=.false.
      endif
!       Rllmin = 0.2d0 
!       Rsep_jlmin = 0.6d0
!       dely_jl=0.6
!       ptjetmin = 15d0
!       yjetmax = 4.5d0
!       mjjmin = 100d0
! c      deltay_jjmin = 4.2d0
! c      jet_opphem = .true.
!       ptlepmin = 10d0
!       ylepmax = 4.5d0
!       Rsep_jlmin = 0.2d0
      

       goto 811
c     ACTIVATE CUTS!!!
      Rsep_jjmin = 0.7d0   
      ptalljetmin = 30d0
      ptjetmin = 30d0
      yjetmax = 5d0
      mjjmin = 600d0
      deltay_jjmin = 4.2d0
      jet_opphem = .true.


      
       goto 811

c     LEPTONIC CUTS
      ptlepmin = 20d0
      ylepmax = 2.5d0
      Rsep_jlmin = 0.6d0
      ylep_between_jets = .true.
      Rllmin = 0.2d0   
      mllmin=20d0

 811  continue

      if (ini) then
         write(*,*) '**************************************************'
         write(*,*) '**************************************************'
         write(*,*) '                ANALYSIS CUTS                     '
         write(*,*) '**************************************************'
         write(*,*) '**************************************************'
         write(*,*) 'Rsep_jjmin = ',Rsep_jjmin
         write(*,*) 'ptalljetmin = ',ptalljetmin
         write(*,*) 'ptjetmin = ',ptjetmin
         write(*,*) 'yjetmax = ',yjetmax
         write(*,*) 'mjjmin = ',mjjmin 
         write(*,*) 'deltay_jjmin = ',deltay_jjmin
         write(*,*) 'ptlepmin = ',ptlepmin 
         write(*,*) 'ylepmax = ',ylepmax 
         write(*,*) 'Rsep_jlmin = ',Rsep_jlmin
         write(*,*) 'ylep_between_jets = ',ylep_between_jets 
         write(*,*) 'jet_opphem = ',jet_opphem 
         write(*,*) '**************************************************'
         write(*,*) '**************************************************'
         ini = .false.

      endif



      do ihep=1,nhep
         if ((abs(idhep(ihep)).eq.11).or.(abs(idhep(ihep)).eq.13)) then
            do mu=1,3
               plep(mu,1)=phep(mu,ihep)
            enddo
            plep(0,1)=phep(4,ihep)
         endif
         if ((abs(idhep(ihep)).eq.12).or.(abs(idhep(ihep)).eq.14)) then
            do mu=1,3
               plep(mu,2)=phep(mu,ihep)          
            enddo
            plep(0,2)=phep(4,ihep)
         endif         
            do mu=0,3
               pH(mu)=plep(mu,2)+plep(mu,1)          
            enddo
      enddo

      do mu=0,3
         pH(mu) = plep(mu,1)+plep(mu,2)
      enddo
      
      ptH = sqrt(pH(1)**2+pH(2)**2)
      yH= getrapidity0(pH)
      call getinvmass(pH,inv_mH)

      
c     set up arrays for jet finding
      do jpart=1,maxtrack
         do mu=1,4
            ptrack(mu,jpart)=0d0
         enddo
         jetvec(jpart)=0
      enddo      
      do jjet=1,maxjet
         do mu=1,4
            pjet(mu,jjet)=0d0
         enddo
      enddo
      j1=0
      found=0
      ntracks=0
      njets=0
c     loop over final state particles to find jets 
      do ihep=1,nhep
         if ((isthep(ihep).eq.1).and.
c     exclude leptons, gauge and Higgs bosons
     #        (((abs(idhep(ihep)).le.10).or.(abs(idhep(ihep)).ge.40))
c     but include gluons 
     #        .or.(abs(idhep(ihep)).eq.21))) then
            if(ntracks.eq.maxtrack) then
               write(*,*)
     #              'hwanal: too many particles, increase maxtrack'
               stop
            endif
c     copy momenta to construct jets 
            ntracks=ntracks+1
            do mu=1,4
               ptrack(mu,ntracks)=phep(mu,ihep)
            enddo
         endif
      enddo



      

      
      if (ntracks.eq.0) then
         return
      endif
      
************************************************************************
*     siscone algorithm
**********************************************************************
c     R = 0.7  radius parameter
c     f = 0.5  overlapping fraction
c.....run the clustering        
c      call fastjetsiscone(ptrack,ntracks,0.7d0,0.5d0,pjet,njets) 
************************************************************************
*     fastkt algorithm
**********************************************************************
c      R = 0.7  Radius parameter
c.....run the clustering 
      R = Rsep_jjmin          
      ptmin_fastkt = 0d0

     
      call fastjetppgenkt(ptrack,ntracks,R,-1,ptmin_fastkt,pjet,njets,
     $                        jetvec)

      


     
c     now we have the jets
      if (njets.gt.0) then
c     find the first THREE hardest jets, if any
         call find_hardest_jets(njets,pjet,3,tagjets,njj)

       
c     at least TWO tagging jets to continue
         if (tagjets.le.1) then
            return
         endif
         
         do ijet=1,tagjets
            do mu=1,3
               pj(mu,ijet)=pjet(mu,njj(ijet))
            enddo
            pj(0,ijet)=pjet(4,njj(ijet))
         enddo

         do ijet=1,tagjets
            ptj(ijet) = sqrt(pj(1,ijet)**2 + pj(2,ijet)**2)
            yj(ijet) = getrapidity0(pj(0,ijet))
         enddo

         
         do i=1,2
            call legoy(plep(0,i), pT_lep(i),y_lep(i),phi_lep(i))
         enddo

c     compute min R_jlep separations

         Rmin = 1d50
cccc already done by kt-algorithm
         do jet=1,2   !tagjets-1
            do lep=jet+1,tagjets
               Rsepjj=rsepn(pj(0,jet),pj(0,lep))
               Rminjj=min(Rmin,Rsepjj)
            enddo
         enddo


         Rmin = 1d50
         do jet=1,2!!njets
            do lep=1,1
               Rsep(jet,lep)=rsepn(pj(0,jet),plep(0,lep))
               Rmin=min(Rmin,Rsep(jet,lep))
            enddo
         enddo

         Rll=rsepn(plep(0,1),plep(0,2))
         
         mll=mjj(plep(0,1),plep(0,2))

         if (yj(1).gt.yj(2)) then
            delphi_jj = azi(pj(0,1))-azi(pj(0,2))
         else
            delphi_jj = azi(pj(0,2))-azi(pj(0,1))
         endif
 
      
         if (delphi_jj.gt.Pi) then
            delphi_jj = delphi_jj - 2d0*pi
         else if (delphi_jj.lt.-Pi) then
            delphi_jj = delphi_jj + 2d0*pi
         endif
         delphi_jj = delphi_jj * 180d0 / Pi
 
c     compute invariant mass of the Hjj system
         do mu=0,3
            pHjj(mu) = pH(mu)+pj(mu,1)+pj(mu,2)
         enddo
         mHjj = sqrt(abs(phjj(0)**2-phjj(1)**2-phjj(2)**2-phjj(3)**2))

       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCC               APPLY CUTS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         pass_cuts = ((pT_lep(1)).gt.ptlepmin) .and.
     $        (abs(y_lep(1)).lt.ylepmax) .and.
     $        (min(pTj(1),pTj(2)).gt.ptjetmin) .and.
     $        (max(abs(yj(1)),abs(yj(2))).lt.yjetmax) .and.
     $        (Rmin.gt.Rsep_jlmin).and.(rminjj.gt.rsep_jjmin)

         if (ylep_between_jets) then
            pass_cuts = pass_cuts .and.
     $           ((min(yj(1),yj(2))+dely_jl).lt.y_lep(1)) .and.
     $           ((max(yj(1),yj(2))-dely_jl).gt.y_lep(1)) !.and.
         endif
         

         if (jet_opphem) then
            pass_cuts = pass_cuts .and. 
     $           (yj(1)*yj(2).lt.0)
         endif

         if (tagjets.lt.3) then
            exist3rdjet = .false.
         else
            exist3rdjet = (pTj(3).gt.ptalljetmin) .and.
     $           (abs(yj(3)).lt.yjetmax)
         endif

         

         pass_cuts_no_mjjmin = pass_cuts .and. 
     $        (abs(yj(1)-yj(2)).gt.deltay_jjmin)
         
         pass_cuts_no_deltayjjmin = pass_cuts .and.
     $        (mjj(pj(0,1),pj(0,2)).gt.mjjmin)
         
         pass_cuts = pass_cuts .and.
     $        (mjj(pj(0,1),pj(0,2)).gt.mjjmin) .and.
     $        (abs(yj(1)-yj(2)).gt.deltay_jjmin)

         


            
         if(njets.gt.3) then            
            ptj(4) = sqrt(pjet(1,4)**2 + pjet(2,4)**2)
            call getrapidity(pjet(1,4),yj(4) )         
            if (passcuts4j) then
               pass_cuts = pass_cuts .and.(pTj(4).gt.ptalljetmin) .and.
     $           (abs(yj(4)).lt.yjetmax)
            endif
         endif
         
         if(exist3rdjet) then
           Rsepjj=rsepn(pj(0,3),plep(0,1))
           pass_cuts=pass_cuts.and.(Rsepjj.gt.Rsep_jlmin)
         endif

         if (pass_cuts) then    

            diag=1
            call pwhgfill(diag,abs(yj(1) - yj(2)),dsig/binsize(diag))

            diag=diag+1
            call pwhgfill(diag,yj(1) - yj(2),dsig/binsize(diag))
            
            diag =diag+ 1   
            call pwhgfill(diag,0.5d0,dsig)

            diag=diag+1
            call pwhgfill(diag,mjj(pj(0,1),pj(0,2)),dsig/binsize(diag))

            diag=diag+1            
            if (exist3rdjet) then
               call pwhgfill(diag,ptj(3),dsig/binsize(diag))
            endif
            
            
            diag=diag+1
            call pwhgfill(diag,ptj(2),dsig/binsize(diag))

            diag=diag+1
            call pwhgfill(diag,ptj(1),dsig/binsize(diag))

            diag=diag+1
            call pwhgfill(diag,pt_lep(2),dsig/binsize(diag))

            diag=diag+1
            call pwhgfill(diag,pt_lep(1),dsig/binsize(diag))
          

            diag=diag+1
            if (exist3rdjet)  then
               call pwhgfill(diag,yj(3),dsig/binsize(diag))
            endif

            diag=diag+1
            call pwhgfill(diag,yj(2),dsig/binsize(diag))

            diag=diag+1
            call pwhgfill(diag,yj(1),dsig/binsize(diag))

            diag=diag+1
            if (exist3rdjet) then 
               call pwhgfill(diag,yj(3)-0.5*(yj(1)+yj(2)),
     $              dsig/binsize(diag))
            endif

            diag=diag+1
            call pwhgfill(diag,min(abs(yj(1)),abs(yj(2))),
     $           dsig/binsize(diag))

            diag=diag+1
            call pwhgfill(diag,max(abs(yj(1)),abs(yj(2))),
     $           dsig/binsize(diag))

            diag=diag+1
            call pwhgfill(diag,delphi_jj,dsig/binsize(diag))
 

         if (y_lep(1).gt.y_lep(2)) then
            delphi_jj = phi_lep(1)-phi_lep(2)
         else
            delphi_jj = phi_lep(2)-phi_lep(1)
         endif
 
  
         if (delphi_jj.gt.Pi) then
            delphi_jj = delphi_jj - 2d0*pi
         else if (delphi_jj.lt.-Pi) then
            delphi_jj = delphi_jj + 2d0*pi
         endif
         delphi_jj = delphi_jj * 180d0 / Pi

         diag=diag+1
         call pwhgfill(diag,delphi_jj,dsig/binsize(diag))       




         diag=diag+1
         call pwhgfill(diag,mHjj,dsig/binsize(diag))   


c     now count how many jets have pt > ptalljetmin
c     find the first 10 hardest jets, if any
            call find_hardest_jets(njets,pjet,10,tagjets,njj)   
c     two jets surely already exists at this step
            njets_passcut = 2
            do ijet=3,tagjets
               ptj(ijet)=
     #              sqrt(pjet(1,njj(ijet))**2 + pjet(2,njj(ijet))**2)
               call getrapidity(pjet(1,njj(ijet)),yj(ijet))    
           
               if (ptj(ijet).gt.ptalljetmin .and.
     #             abs(yj(ijet)).lt.yjetmax) then
                   do mu=1,3
                      pjN(mu)=pjet(mu,njj(ijet))
                   enddo
                   pjN(0)=pjet(4,njj(ijet))
                   Rsepjj=rsepn(pjN(0),plep(0,1))
                   if(Rsepjj.gt.Rsep_jlmin) njets_passcut = njets_passcut + 1
               endif
            enddo
c            write(*,*) 'njets_passcut ',njets_passcut
            diag=diag+1
            call pwhgfill(diag,njets_passcut*1d0,dsig/binsize(diag))

c     jet multeplicity between the two tagging jets
            njets_passcut = 0
            do ijet=3,tagjets
               if (min(yj(1),yj(2)).lt.yj(ijet) .and.
     #                 yj(ijet).lt.max(yj(1),yj(2)) .and.
     #                 ptj(ijet).gt.ptalljetmin) then
                   do mu=1,3
                      pjN(mu)=pjet(mu,njj(ijet))
                   enddo
                   pjN(0)=pjet(4,njj(ijet))
                   Rsepjj=rsepn(pjN(0),plep(0,1))
                   if(Rsepjj.gt.Rsep_jlmin) njets_passcut = njets_passcut + 1
               endif
            enddo            
            diag=diag+1
            call pwhgfill(diag,njets_passcut*1d0,dsig/binsize(diag))

      endif
      endif
      

      
      
      
      end


