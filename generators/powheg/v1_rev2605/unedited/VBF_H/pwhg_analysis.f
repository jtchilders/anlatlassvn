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

      if (ini) then
         write(*,*) '********************************************'
         write(*,*) '********************************************'
         write(*,*) 'inv Higgs boson mass plot done assuming the '
         write(*,*) 'following values'
         write(*,*) 'ph_Hmass = ',ph_Hmass
         write(*,*) 'ph_Hwidth = ',ph_Hwidth
         write(*,*) '********************************************'
         write(*,*) '********************************************'
         ini=.false.
      endif

      nsigma = 3
      invmasslow =ph_Hmass-nsigma*ph_Hwidth
      invmasshigh=ph_Hmass+nsigma*ph_Hwidth
      step = 2*nsigma*ph_Hwidth/50
      cut = ' WBF cuts '
      ymax = 7.2d0

      call pwhginihist
      diag=1
      binsize(diag) = 10d0
      call pwhgbookup(diag,'pt H ','LIN',binsize(diag),0d0,400d0)

      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'y H ','LIN',binsize(diag),-ymax,ymax)

      diag=diag+1
      binsize(diag) = step
      call pwhgbookup(diag,'inv mass ','LOG',binsize(diag),invmasslow,
     #     invmasshigh)

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'yj(1) - yj(2) ','LIN',binsize(diag),
     #     -ymax,ymax)

      diag=diag+1
      binsize(diag) = 40d0
      call pwhgbookup(diag,'mjj ','LIN',binsize(diag),0d0,3600d0)

      diag=diag+1
      binsize(diag) = 40d0
      call pwhgbookup(diag,'mjj'//cut,'LIN',binsize(diag),0d0,3600d0)

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'abs(yj(1) - yj(2))'//cut,'LIN',
     #     binsize(diag),0d0,2*ymax)

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'yj(1) - yj(2)'//cut,'LIN',
     #     binsize(diag),-ymax,ymax)

      diag=diag+1
      binsize(diag) = 40d0
      call pwhgbookup(diag,'mjj'//cut,'LIN',binsize(diag),0d0,3600d0)

      diag=diag+1
      binsize(diag) = 10d0
      call pwhgbookup(diag,'ptj(3)'//cut,'LOG',binsize(diag),0d0,150d0)

      diag=diag+1
      binsize(diag) = 10d0
      call pwhgbookup(diag,'ptj(2)'//cut,'LOG',binsize(diag),0d0,250d0)

      diag=diag+1
      binsize(diag) = 10d0
      call pwhgbookup(diag,'ptj(1)'//cut,'LOG',binsize(diag),0d0,400d0)

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'yj(3)'//cut,'LIN',binsize(diag),-ymax,ymax)
 

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'yj(2)'//cut,'LIN',binsize(diag),-ymax,ymax)


      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'yj(1)'//cut,'LIN',binsize(diag),-ymax,ymax)

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'yj(3)-0.5*(yj(1)+yj(2))'//cut,'LIN',
     #     binsize(diag),-ymax,ymax) 

      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'min(abs(yj(1)),abs(yj(2)))'//cut,'LIN',
     #     binsize(diag),0d0,ymax)


      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'max(abs(yj(1)),abs(yj(2)))'//cut,'LIN',
     #     binsize(diag),0d0,ymax)

      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'delphi_jj'//cut,'LIN',
     1  binsize(diag),0d0,3.2d0)


      diag=diag+1
      binsize(diag) = 40d0
      call pwhgbookup(diag,'mHjj ','LIN',binsize(diag),0d0,3600d0)

      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'ptrel j1'//cut,'LOG',binsize(diag),0d0,30d0)

      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'ptrel j2'//cut,'LOG',binsize(diag),0d0,30d0)

      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'jet multip'//cut,'LOG',
     #     binsize(diag),1.5d0,8.5d0)

      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'jet mult. btw tag jets'//cut,'LOG',
     #     binsize(diag),0.5d0,8.5d0)

      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'Pveto'//cut,'LIN',binsize(diag),10d0,60d0)


      end

      
     
      subroutine analysis(dsig)
      implicit none
      real * 8 dsig
      include 'hepevt.h'
      include 'pwhg_math.h'      
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
     #     ptrackin(0:3),ptrackout(0:3)
      integer i,diag,njets_passcut
      external get_ptrel
      real * 8 R,ptmin_fastkt
      integer jetvec(maxtrack)
      logical ini
      data ini/.true./
      save ini
      integer HZZ,HWW
      integer idH
      real * 8 pH(0:3),ptH,yH,inv_mH,Edec,thl,phil,plepCM(0:3,2),
     $     plep(0:3,2),mod_vecH,pj(0:3,3),pT_lep(2),y_lep(2),pHjj(0:3),
     $     mHjj
      real * 8 random,rsepn,getrapidity0,mjj,azi
      external random,rsepn,getrapidity0,mjj,azi
      real * 8 Rsep(2,2),Rmin,delphi_jj
      logical pass_cuts,pass_cuts_no_mjjmin,pass_cuts_no_deltayjjmin
      real * 8 ptjetmin,ptalljetmin,yjetmax,deltay_jjmin,ptlepmin,
     #     ylepmax,mjjmin,Rsep_jlmin,ptvetojet,sig2NLO
      logical ylep_between_jets,jet_opphem,exist3rdjet
      logical onlyquarks,Z_exchange
      real * 8 binsize(100)
      common/pwhghistcommon/binsize
c      logical iniptcut
c      save iniptcut
c      data iniptcut/.true./
      logical higgsfound,found_hardest_vetojet    
c     we need to tell to the this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      
      
c     DISABLE ALL CUTS
      ptalljetmin = 0d0            
      ptjetmin = 0d0         
      yjetmax = 100d0           
      mjjmin = 0d0              
      deltay_jjmin = 0d0        

      ptlepmin = 0d0            
      ylepmax = 100d0           
      Rsep_jlmin = 0d0          

      ylep_between_jets = .false.
      jet_opphem = .false.

c      goto 811
c     ACTIVATE CUTS!!!
      ptalljetmin = 20d0
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

 811  continue

      if (ini) then
         write(*,*) '**************************************************'
         write(*,*) '**************************************************'
         write(*,*) '                ANALYSIS CUTS                     '
         write(*,*) '**************************************************'
         write(*,*) '**************************************************'
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

      idH=0
      higgsfound = .false.
c     find Higgs boson
      call particle_identif(HWW,HZZ)      
      do ihep=1,nhep
c     set up different searching stategies according to the shower Monte Carlo
c     program used      
         if ((WHCPRG.eq.'HERWIG').and.(idhep(ihep).eq.25).and.
     $        (isthep(ihep).eq.155)) then 
            higgsfound = .true.      
         elseif ((WHCPRG.eq.'PYTHIA').and.(idhep(ihep).eq.25).and.
     #           (ISTHEP(ihep).eq.1)) then 
            higgsfound = .true.      
         elseif (( (WHCPRG.eq.'POWHEG').or.(WHCPRG.eq.'NLO   ').or.
     $           (WHCPRG.eq.'LHE   ')).and.
     #           ((idhep(ihep).eq.HWW).or.(idhep(ihep).eq.HZZ).or.
     #           (idhep(ihep).eq.25) )) then
            higgsfound = .true.            
         endif
         if (higgsfound) then
            idH=ihep
            idhep(ihep)=25
            goto 333
         endif
      enddo
 333  continue
      if (idH.eq.0) then
         write(*,*) 'HIGGS NOT FOUND'
         call exit(2)
      endif

      do mu=1,3
         pH(mu) = phep(mu,idH)
      enddo
      pH(0) = phep(4,idH)
      
      ptH = sqrt(pH(1)**2+pH(2)**2)
      call getrapidity(phep(1,idH),yH)
      call getinvmass(phep(1,idH),inv_mH)
      
      diag=1
      call pwhgfill(diag,ptH,dsig/binsize(diag))
      diag=diag+1
      call pwhgfill(diag,yH,dsig/binsize(diag))
      diag=diag+1
      call pwhgfill(diag,inv_mH,dsig/binsize(diag))

      
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
      R = 0.7d0          
      ptmin_fastkt = 0d0
      call fastjetktwhich(ptrack,ntracks,ptmin_fastkt,R,
     #     pjet,njets,jetvec) 
      
c     if (njets.eq.3) then
c     write(*,*) 'njets ',njets
c     endif     
c     njets=5
c     do i=1,njets
c     do mu=0,3
c     pjet(mu,i) = random()
c     enddo
c     enddo

     
c     now we have the jets
      if (njets.gt.0) then
c     find the first THREE hardest jets, if any
         call find_hardest_jets(njets,pjet,3,tagjets,njj)
c     write(*,*) 'tagjets ',tagjets
c     write(*,*) 'jet found new',njj
         
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
         
c     get pt's and rapidities of the jets
         do ijet=1,tagjets
            ptj(ijet) = sqrt(pj(1,ijet)**2 + pj(2,ijet)**2)
            yj(ijet) = getrapidity0(pj(0,ijet))
         enddo
         
c     decay the Higgs boson in the Higgs boson CM in order to impose
c     cuts on the decay products, reffered to as massless leptons here
         Edec = inv_mH/2
         thl = pi*random()
         phil = 2*pi*random()
         plepCM(1,1)=Edec*sin(thl)*cos(phil)
         plepCM(2,1)=Edec*sin(thl)*sin(phil)
         plepCM(3,1)=Edec*cos(thl)
         plepCM(0,1)=Edec
         plepCM(0,2)=Edec
         do mu=1,3
            plepCM(mu,2)=-plepCM(mu,1)
         enddo
         mod_vecH=sqrt(pH(1)**2+pH(2)**2+pH(3)**2)
         do mu=1,3
            vec(mu) = pH(mu)/mod_vecH
         enddo
         beta = mod_vecH/pH(0)
         call mboost(2,vec,beta,plepCM(0,1),plep(0,1))
c     debug
c     do mu=1,3
c     vec(mu) = -vec(mu)
c     enddo
c     call mboost(1,vec,beta,pH(0),ptmp(0))
c     write(*,*) 'only time component!!'
c     write(*,*) ptmp
c     do mu=0,3
c     write(*,*) pH(mu)-plep(mu,1)-plep(mu,2)
c     enddo
         
         do i=1,2
            pT_lep(i) = sqrt(plep(1,i)**2+plep(2,i)**2)               
            y_lep(i) = getrapidity0(plep(0,i))
         enddo
c     compute min R_jlep separations
         Rmin = 1d50
         do jet=1,2
            do lep=1,2
               Rsep(jet,lep)=rsepn(pj(0,jet),plep(0,lep))
               Rmin=min(Rmin,Rsep(jet,lep))
            enddo
         enddo
 
         delphi_jj = abs(azi(pj(0,1))-azi(pj(0,2)))
         if (delphi_jj.gt.pi) then
            delphi_jj = 2*pi-delphi_jj
         endif
 
c     compute invariant mass of the Hjj system
         do mu=0,3
            pHjj(mu) = pH(mu)+pj(mu,1)+pj(mu,2)
         enddo
         mHjj = sqrt(abs(phjj(0)**2-phjj(1)**2-phjj(2)**2-phjj(3)**2))

         diag=diag+1
         call pwhgfill(diag,yj(1) - yj(2),dsig/binsize(diag))
         diag=diag+1
         call pwhgfill(diag,mjj(pj(0,1),pj(0,2)),dsig/binsize(diag))


       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCC               APPLY CUTS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         pass_cuts = 
     $        (min(pT_lep(1),pT_lep(2)).gt.ptlepmin) .and.
     $        (max(abs(y_lep(1)),abs(y_lep(2))).lt.ylepmax) .and.
     $        (min(pTj(1),pTj(2)).gt.ptjetmin) .and.
     $        (max(abs(yj(1)),abs(yj(2))).lt.yjetmax) .and.
     $        (Rmin.gt.Rsep_jlmin)
         
         if (ylep_between_jets) then
            pass_cuts = pass_cuts .and.
     $           (min(yj(1),yj(2)).lt.y_lep(1)) .and.
     $           (max(yj(1),yj(2)).gt.y_lep(1)) .and.
     $           (min(yj(1),yj(2)).lt.y_lep(2)) .and.
     $           (max(yj(1),yj(2)).gt.y_lep(2))
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

         
         diag=diag+1
         if (pass_cuts_no_mjjmin) then
            call pwhgfill(diag,mjj(pj(0,1),pj(0,2)),dsig/binsize(diag))
         endif

         diag=diag+1
         if (pass_cuts_no_deltayjjmin) then 
            call pwhgfill(diag,abs(yj(1) - yj(2)),dsig/binsize(diag))
         endif
            

         if (pass_cuts) then    
c            call increasecnt('pass WBF cuts')
            diag=diag+1
            call pwhgfill(diag,yj(1) - yj(2),dsig/binsize(diag))

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

            diag=diag+1
            call pwhgfill(diag,mHjj,dsig/binsize(diag))

c     loop on the hardest and next-to-hardest jet
            do ijet=1,2
               vec(1)=0d0
               vec(2)=0d0
               vec(3)=1d0
               beta = -pj(3,ijet)/pj(0,ijet)
               call mboost(1,vec,beta,pj(0,ijet),pjetout)  
c               write(*,*) pjetout
               ptrel = 0
               do i=1,ntracks
                  if (jetvec(i).eq.njj(ijet)) then
                     do mu=1,3
                        ptrackin(mu) = ptrack(mu,i)
                     enddo
                     ptrackin(0) = ptrack(4,i)
                     call mboost(1,vec,beta,ptrackin,ptrackout) 
                     ptrel = ptrel + get_ptrel(ptrackout,pjetout)
c                     write(*,*) ijet,ptrel
                  endif
               enddo
               diag=diag+1
               if (ijet.eq.1) then 
                  call pwhgfill(diag,ptrel,dsig/binsize(diag))
               endif
c               diag=diag+1
               if (ijet.eq.2) then 
                  call pwhgfill(diag,ptrel,dsig/binsize(diag))
               endif
            enddo         


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
                  njets_passcut = njets_passcut + 1
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
                  njets_passcut = njets_passcut + 1
               endif
            enddo            
            diag=diag+1
            call pwhgfill(diag,njets_passcut*1d0,dsig/binsize(diag))

c     Pveto as in eq. (3.14) of arXiv:0710.5621v2
            found_hardest_vetojet = .false.
            do ijet=3,tagjets               
               if (.not.found_hardest_vetojet .and. 
     #                 (min(yj(1),yj(2)).lt.yj(ijet) .and.
     #                 yj(ijet).lt.max(yj(1),yj(2)))) then
                  found_hardest_vetojet = .true.
                  ptvetojet = ptj(ijet)
               endif
            enddo            
            diag=diag+1
            sig2NLO = 0.72331d0
 111        continue
c     NB: NO division by binsize(diag) since we want the integral
            call pwhgfill(diag,ptvetojet,dsig/sig2NLO)
            ptvetojet = ptvetojet - binsize(diag)
            if (ptvetojet.gt.10d0) goto 111
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
      real * 8 p(4),m
      m=sqrt(abs(p(4)**2-p(1)**2-p(2)**2-p(3)**2))
      end



c      subroutine pwhgfillup(n,x,y)
c      implicit none
c      real * 8 x,y
c      integer n
c      call pwhgfill(n,x,y)
c      call pwhgfill(n+100,x,y*y)
c      end





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



      subroutine particle_identif(HWW,HZZ)
      implicit none
      integer pdg_Higgs,pdg_Z,pdg_W,HZZ,HWW
      pdg_Higgs = 25
      pdg_Z = 23
      pdg_W = 24      
c     build an identifier for Higgs production in WW and ZZ fusion 
      HWW = 10000*pdg_W + pdg_Higgs
      HZZ = 10000*pdg_Z + pdg_Higgs
      end




c      subroutine topout
c      implicit none
c      include 'hepevt.h'
c      character * 50 title
c      integer i
c      integer maxnumplot
c      common/cmaxnumplot/maxnumplot
c      logical lin(100)
c      common/lin_scale/lin
c      character * 3 scale
cc     
cc     If histogram I contains accumulated weights and
cc     histogram I+100 contains its squared values,
cc     then a temporary copy of both is made in order
cc     to safely have intermediate results
c
c      do i=1,maxnumplot
c	call pwhgfinal(i)
c        call pwhgcopy(i,i+200)
c        call pwhgcopy(i+100,i+300)
c        call pwhgopera(i+200,'F',i+200,i+200,1d0/dble(nevhep),0d0)
c        call pwhgerror(i+200,i+300,dble(nevhep))
c        call pwhgfinal(i+200)
c        call pwhgfinal(i+300)
c      enddo
c      do i=1,maxnumplot
c         call pwhggettitle(i+200,title)
c         if (lin(i)) then
c            scale = 'LIN'
c         else
c            scale = 'LOG'
c         endif
c         call pwhgmultitop(i+200,i+300,2,3,title,' ',scale)
c      enddo
c      end            





