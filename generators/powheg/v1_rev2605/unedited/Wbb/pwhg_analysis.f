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
c     binsize
      real * 8 bsz(100)
      common/pwhghistcommon/bsz
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      integer nptWcutmax,nptWcut
      parameter(nptWcutmax = 10)
      real * 8 ptWcuts(0:nptWcutmax-1)
      common/cptvbcut/ptWcuts,nptWcut,numplots
      integer ncut,numplots
      
      character * 10 cut

      ptWcuts(0) = 0d0
      ptWcuts(1) = 50d0
c     number of pt W cuts implemented
      nptWcut=2

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      numplots = 34  ! <========== DO NOT FORGET TO SET THIS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      

      call pwhginihist

      do ncut=0,nptWcut-1
      write(unit=cut,fmt="(f5.2)") ptWcuts(ncut)

      diag = 1
      bsz(diag) = 1d0
      call pwhgbookup(diag+numplots*ncut,'total Xsec ptW>'//cut,'LIN',
     $     bsz(diag),0d0,1d0)
   
      diag=2
      bsz(diag) = 10d0
      call pwhgbookup(diag+numplots*ncut,'ptb leading ptW>'//cut,'LOG',
     $     bsz(diag),0d0,400d0)

      diag=3
      bsz(diag) = 10d0
      call pwhgbookup(diag+numplots*ncut,'ptb sublead ptW>'//cut,'LOG',
     $     bsz(diag),0d0,400d0)
      
      diag=4
      bsz(diag) = 0.4d0
      call pwhgbookup(diag+numplots*ncut,'etab leading ptW>'//cut,'LOG',
     $     bsz(diag),-5d0,5d0)

      diag=5
      bsz(diag) = 0.4d0
      call pwhgbookup(diag+numplots*ncut,'etab sublead ptW>'//cut,'LOG',
     $     bsz(diag),-5d0,5d0)

      diag=6
      bsz(diag) = 10d0
      call pwhgbookup(diag+numplots*ncut,'ptW ptW>'//cut,'LOG',
     $     bsz(diag),0d0,400d0)

      diag=7
      bsz(diag) = 0.4d0
      call pwhgbookup(diag+numplots*ncut,'etaW ptW>'//cut,'LOG',
     $     bsz(diag),-5d0,5d0)

      diag=8
      bsz(diag) = 3d0
      call pwhgbookup(diag+numplots*ncut,'mbb ptW>'//cut,'LOG',
     $     bsz(diag),0d0,300d0)

      diag=9
      bsz(diag) = 0.2d0
      call pwhgbookup(diag+numplots*ncut,'Rbb ptW>'//cut,'LOG',
     $     bsz(diag),0d0,6d0)

      diag=10
      bsz(diag) = 0.1d0
      call pwhgbookup(diag+numplots*ncut,'dphi(bb,W) ptW>'//cut,'LOG',
     $     bsz(diag),0d0,3.2d0)

      diag=11
      bsz(diag) = 0.4d0
      call pwhgbookup(diag+numplots*ncut,'deta(bb,W) ptW>'//cut,'LOG',
     $     bsz(diag),-5d0,5d0)

      diag=12
      bsz(diag) = 10d0
      call pwhgbookup(diag+numplots*ncut,'ptl ptW>'//cut,'LOG',
     $     bsz(diag),0d0,400d0)

      diag=13
      bsz(diag) = 0.4d0
      call pwhgbookup(diag+numplots*ncut,'etal ptW>'//cut,'LOG',
     $     bsz(diag),-5d0,5d0)

      diag=14
      bsz(diag) = 10d0
      call pwhgbookup(diag+numplots*ncut,'ptj ptW>'//cut,'LOG',
     $     bsz(diag),0d0,400d0)

      diag=15
      bsz(diag) = 10d0
      call pwhgbookup(diag+numplots*ncut,'ptj yj< yjcut ptW>'//cut,
     $     'LOG',bsz(diag),0d0,400d0)

      diag=16
      bsz(diag) = 0.4d0
      call pwhgbookup(diag+numplots*ncut,'yj ptW>'//cut,'LOG',
     $     bsz(diag),-5d0,5d0)

      diag=17
      bsz(diag) = 0.4d0
      call pwhgbookup(diag+numplots*ncut,'yj ptj>ptjcut ptW>'//cut,
     $     'LOG',bsz(diag),-5d0,5d0)

      diag=18
      bsz(diag) = 0.4d0
      call pwhgbookup(diag+numplots*ncut,'ybbW-yj ptW>'//cut,
     $     'LOG',bsz(diag),-5d0,5d0)
   
      diag=19
      bsz(diag) = 1d0
      call pwhgbookup(diag+numplots*ncut,'inv mass W ptW>'//cut,
     $     'LOG',bsz(diag),50d0,110d0)
  
      diag=20
      bsz(diag) = 0.4d0
      call pwhgbookup(diag+numplots*ncut,'deta(b,b) ptW>'//cut,
     $     'LOG',bsz(diag),-5d0,5d0)
 
      diag=21
      bsz(diag) = 0.1d0
      call pwhgbookup(diag+numplots*ncut,'dphi(b,b) ptW>'//cut,
     $     'LOG',bsz(diag),0d0,3.2d0)
  
      diag=22
      bsz(diag) = 1d0
      call pwhgbookup(diag+numplots*ncut,'ptj zoom ptW>'//cut,'LOG',
     $     bsz(diag),0d0,50d0)

      diag=23
      bsz(diag) = 10d0
      call pwhgbookup(diag+numplots*ncut,'ptj qg ptW>'//cut,'LOG',
     $     bsz(diag),0d0,400d0)

      diag=24
      bsz(diag) = 10d0
      call pwhgbookup(diag+numplots*ncut,'ptj qq ptW>'//cut,'LOG',
     $     bsz(diag),0d0,400d0)

      diag=25
      bsz(diag) = 10d0
      call pwhgbookup(diag+numplots*ncut,'sqrt(shat) ptW>'//cut,'LOG',
     $     bsz(diag),0d0,2000d0)

      diag=26
      bsz(diag) = 10d0
      call pwhgbookup(diag+numplots*ncut,'m_bbW ptj>0 ptW>'//cut,'LOG',
     $     bsz(diag),0d0,2000d0)

      diag=27
      bsz(diag) = 10d0
      call pwhgbookup(diag+numplots*ncut,'m_bbW ptj>100 ptW>'//cut,
     $     'LOG',bsz(diag),0d0,2000d0)

      diag=28
      bsz(diag) = 10d0
      call pwhgbookup(diag+numplots*ncut,'m_bbW ptj>200 ptW>'//cut,
     $     'LOG',bsz(diag),0d0,2000d0)

      diag=29
      bsz(diag) = 10d0
      call pwhgbookup(diag+numplots*ncut,'m_bbW ptj>300 ptW>'//cut,
     $     'LOG',bsz(diag),0d0,2000d0)

      diag=30
      bsz(diag) = 10d0
      call pwhgbookup(diag+numplots*ncut,'m_bbW ptj>400 ptW>'//cut,
     $     'LOG',bsz(diag),0d0,2000d0)

      diag=31
      bsz(diag) = 10d0
      call pwhgbookup(diag+numplots*ncut,'m_bbW no light jet ptW>'//cut,
     $     'LOG',bsz(diag),0d0,2000d0)

      diag=32
      bsz(diag) = 10d0
      call pwhgbookup(diag+numplots*ncut,'sqrt(shat) btl ptW>'//cut,
     $     'LOG',bsz(diag),0d0,2000d0)

      diag=33
      bsz(diag) = 10d0
      call pwhgbookup(diag+numplots*ncut,'sqrt(shat) rmn ptW>'//cut,
     $     'LOG',bsz(diag),0d0,2000d0)

      diag=34
      bsz(diag) = 10d0
      call pwhgbookup(diag+numplots*ncut,'m_bbW part ptW>'//cut,
     $     'LOG',bsz(diag),0d0,2000d0)


      enddo
      end 

      
      subroutine analysis(dsig)
      implicit none
      real * 8 dsig
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
c arrays to reconstruct jets
      integer maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real *8 pjet(4,maxjet)
      integer ihep,mu,jetvec(maxtrack)
      logical ini
      data ini/.true./
      save ini
c     bsz
      integer diag
      real * 8 bsz(100)
      common/pwhghistcommon/bsz
c     we need to tell to this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer maxnumlep
      parameter (maxnumlep=100)
      integer lepvec(maxnumlep),nuvec(maxnumlep)
      integer iheparr_w(maxnumlep),iheparr_l(maxnumlep),
     $     iheparr_v(maxnumlep),nbjet_array(maxnumlep),
     $     nbbarjet_array(maxnumlep)
      integer ilep,inu,lep,nu,nlep,nnu,jlep,jnu,i,j
      real *8 mV2ref,mV2
      real *8 Wmass,Wwidth,Wmass2low,Wmass2high,mW,mW2
      logical foundlep
      real * 8 etabmax,etabmin,mbb,ptw,etaw,pb(4),pbb(4),pw(4),pbbb(4),
     $     ptbmin,ptbmax,ptb,ptbb,pj(4),ptj
      character * 20 process
      integer ist_w,nw,nl,nv,nb,nbbar,ist,id,njets,ihepin1,ihepin2
      logical is_B_hadron,is_BBAR_hadron
      real * 8 p_w(4,maxnumlep),p_l(4,maxnumlep),p_b(4,maxnumlep),
     $     p_bbar(4,maxnumlep),p_v(4,maxnumlep)
      external is_B_hadron,is_BBAR_hadron
      real * 8 powheginput
      external powheginput
      integer idvecbos,vdecaymode,chargelept,neutralept,sign
      save idvecbos,vdecaymode,chargelept,neutralept
      integer invalidisthep
      parameter (invalidisthep = 23751)                
      integer ihepb, ihepbbar
      real * 8 pl(1:4),pv(1:4),ptl,ptv,tmp
      integer nbjet,nbbarjet,ihepl,ihepv,ihepw,jetinfo(maxjet)
      integer ntracks,ihep_of_track(maxtrack)
      real * 8 ptrack(4,maxtrack)
      common/cptrack/ihep_of_track,ptrack,ntracks
      real * 8 pt
      logical found_hardjet,found_bjet1,found_bjet2
      real * 8 phardjet(4),pbjet1(4),pbjet2(4),pin1(4),pin2(4)
      integer typeb1,typeb2
      integer kill,kept
      save kill,kept
      data kill/0/
      data kept/0/
      integer nptWcutmax,nptWcut
      parameter(nptWcutmax = 10)
      real * 8 ptWcuts(0:nptWcutmax-1)
      common/cptvbcut/ptWcuts,nptWcut,numplots
      integer ncut,numplots
      real * 8 ptb1,ptb2,etab1,etab2,pbbjet(4),Rbb,dphi_bbW,
     #     etabb,deta_bbW,etal,yj,ybbW,dy_bbWj,pbbjetW(4),invmW,
     $     deta_bb,dphi_bb,ptot(4),mbbW,sqrshat
      character * 3 machine
      real * 8 ptbcut,etabcut,etaWcut,ptlcut,etalcut,ptvcut,yjcut,
     $     ptjcut 
      logical Wdecay
      real * 8 ebeam1
      save  ptbcut,etabcut,etaWcut,ptlcut,etalcut,ptvcut,yjcut,ptjcut
      integer wrong_bb_sequence
      save  wrong_bb_sequence
      data  wrong_bb_sequence/0/
      logical apply_cuts
      parameter (apply_cuts=.false.)
      real * 8 ptmin_fastkt
      common/cptmin_fastkt/ptmin_fastkt
      data rad_type/1/

      if (ini) then
         ebeam1=powheginput('ebeam1')
         if (ebeam1.eq.980d0) then
            machine='TEV'         
         elseif (ebeam1.eq.7000d0.or.ebeam1.eq.3500d0) then
            machine='LHC'
         else
            write(*,*) 'UNKNOWN accelerator energy. Program stops!'
            call exit(1)
         endif
         if (apply_cuts) then
            if (machine.eq.'LHC') then
               write(*,*) '********************************'            
               write(*,*) '********************************'            
               write(*,*) '************   LHC   ***********'            
               write(*,*) '********************************'            
               write(*,*) '********************************'            
               ptbcut  = 30d0
               etabcut =  2.5d0
               etaWcut =  2.5d0
               ptlcut  = 25d0
               etalcut =  1.5d0
               ptvcut  = 25d0
               yjcut   =  2.8d0
               ptjcut  = 30d0
               ptmin_fastkt = 5d0
            elseif (machine.eq.'TEV') then
               write(*,*) '********************************'            
               write(*,*) '********************************'            
               write(*,*) '************   TEV   ***********'            
               write(*,*) '********************************'            
               write(*,*) '********************************'            
               ptbcut  = 20d0
               etabcut =  2d0
               etaWcut =  2d0
               ptlcut  = 20d0
               etalcut =  1.1d0
               ptvcut  = 25d0
               yjcut   =  2.1d0
               ptjcut  = 20d0
               ptmin_fastkt = 5d0
            else
               write(*,*) 'UNKNOWN machine'
               call exit(1)
            endif
         else
            write(*,*) '********************************'            
            write(*,*) '********************************'            
            write(*,*) '**********  NO CUTS   **********'            
            write(*,*) '********************************'            
            write(*,*) '********************************'            
            ptbcut  =  0d0
            etabcut =  1d10
            etaWcut =  1d10
            ptlcut  =  0d0
            etalcut =  1d10
            ptvcut  =  0d0
            yjcut   =  1d10
            ptjcut  =  0d0
            ptmin_fastkt = 5d0         
         endif
         write(*,*) '********************************'            
         write(*,*) '********************************'            
         write(*,*) '**********    CUTS    **********'                     
         write(*,*) 'ptbcut =', ptbcut  
         write(*,*) 'etabcut =', etabcut
         write(*,*) 'etaWcut =', etaWcut
         write(*,*) 'ptlcut =', ptlcut  
         write(*,*) 'etalcut =',etalcut 
         write(*,*) 'ptvcut =', ptvcut   
         write(*,*) 'yjcut =', yjcut    
         write(*,*) 'ptjcut =', ptjcut    
         write(*,*) 'ptmin_fastkt =', ptmin_fastkt
         write(*,*) '********************************'            
         write(*,*) '********************************'            

         write (*,*)
         write (*,*) '********************************************'
         if(WHCPRG.eq.'NLO   ') then
            write (*,*) '           NLO ANALYSIS CALLED        '
         elseif(WHCPRG.eq.'LHE   ') then
            write (*,*) '           LHE ANALYSIS CALLED        '
         elseif(WHCPRG.eq.'HERWIG') then
            write (*,*) '           HERWIG ANALYSIS CALLED     '
         elseif(WHCPRG.eq.'PYTHIA') then
            write (*,*) '           PYTHIA ANALYSIS CALLED     '
         endif
         write (*,*) '********************************************'
         write (*,*)
         idvecbos=powheginput('idvecbos')
         vdecaymode=powheginput('vdecaymode')
         if (idvecbos.eq.24) then
c     W+ production
            sign=+1
         else
c     W- production            
            sign=-1
         endif
         if (vdecaymode.eq.0) then
c     decay into electron
            chargelept = -11*sign
            neutralept =  12*sign 
         elseif (vdecaymode.eq.1) then
c     decay into electron
            chargelept = -11*sign
            neutralept =  12*sign 
         elseif (vdecaymode.eq.2) then
c     decay into muon
            chargelept = -13*sign
            neutralept =  14*sign 
         elseif (vdecaymode.eq.3) then
c     decay into tau
            chargelept = -15*sign
            neutralept =  16*sign 
         endif
         ini=.false.
      endif

      nw=0
      nl=0
      nv=0
      nb=0
      nbbar=0
      ihepb=0
      ihepbbar=0
c     assume as default that the W decay products are present
      Wdecay=.true.


      if(WHCPRG.eq.'HERWIG') then
c         ist_w=195              !: for undecayed events
         ist_w=155              !: for decayed events
c the 4rd and 5th particle in HERWIG are the two incoming ones         
         ihepin1=4
         ihepin2=5
      elseif(WHCPRG.eq.'PYTHIA') then
         ist_w=3
c the 3rd and 4th particle in PYTHIA are the two incoming ones         
         ihepin1=3
         ihepin2=4
      elseif(WHCPRG.eq.'LHE   ') then
c     lhef_analysis: W present in the list of particles with ist=2,
c     since it is decayed
         ist_w=2
c the 1st and 2nd particle in LHE are the two incoming ones         
         ihepin1=1
         ihepin2=2
c     here we have two option: 
c     1) if we want to compare with NLO distributions, then Wdecay
c        should be set to false
c     2) if we want to compare with the full showered result, then 
c        Wdecay should be set to true
c         Wdecay=.false.
         Wdecay=.true.
      elseif(WHCPRG.eq.'NLO   ') then
c     NLO analysis: W NOT decayed
         Wdecay=.false.
         ist_w=1
c the 1st and 2nd particle in NLO are the two incoming ones         
         ihepin1=1
         ihepin2=2
      endif
      
c     assign incoming particles momenta
      do mu=1,4
         pin1(mu)=phep(mu,ihepin1)
         pin2(mu)=phep(mu,ihepin2)
      enddo

      do ihep=1,nhep
         ist=isthep(ihep)
         id=idhep(ihep)     
c     W
         if(ist.eq.ist_w.and.id.eq.idvecbos) then
            nw=nw+1
            iheparr_w(nw)=ihep
            do mu=1,4
               p_w(mu,nw)=phep(mu,ihep)
            enddo
c     charged lepton               
         elseif(ist.eq.1.and.id.eq.chargelept) then
            nl=nl+1
            iheparr_l(nl)=ihep
            do mu=1,4
               p_l(mu,nl)=phep(mu,ihep)
            enddo
c     neutral lepton
         elseif(ist.eq.1.and.id.eq.neutralept) then
            nv=nv+1
            iheparr_v(nv)=ihep
            do mu=1,4
               p_v(mu,nv)=phep(mu,ihep)
            enddo
c     b hadrons
         elseif((ist.eq.1).and.is_B_hadron(id)) then
c     write(*,*) 'ihep B hadron =',ihep
            nb=nb+1
            do mu=1,4
               p_b(mu,nb)=phep(mu,ihep)
            enddo
c            ihepb=ihep
c     bbar hadrons
         elseif((ist.eq.1).and.is_BBAR_hadron(id)) then
c     write(*,*) 'ihep Bbar hadron =',ihep
            nbbar=nbbar+1
            do mu=1,4
               p_bbar(mu,nbbar)=phep(mu,ihep)
            enddo
c            ihepbbar=ihep
         endif
c     Upsilon meson b bbar
         if (id.eq.553) then
            write(*,*) 'WARNING: ****** FOUND UPSILON ******'
c     write(*,*) '************ FIX CODE ************'     
c     write(*,*) 'ihepb,ihepbbar ',ihepb,ihepbbar  
c     CALL PYLIST(4) 
c     stop
            return
         endif   
      enddo         

c     write(*,*) '======>',nb,nbbar
      
      if (nb*nbbar.eq.0) then
         if (((nb.eq.0).and.(nbbar.ne.0)).or.
     $        ((nb.ne.0).and.(nbbar.eq.0))) then
            write(*,*) 'SEVERE WARNING: ***** One b is missing ******' 
            write(*,*) 'ihepb,ihepbbar ',ihepb,ihepbbar            
c     CALL PYLIST(4) 
            return
         endif
      endif

      
c      if (nw*nl*nv*nb*nbbar.ge.2) then
c         write(*,*) 'summary'      
c         write(*,*) nw,nl,nv,nb,nbbar
c      endif
      

      if (Wdecay.and.nl*nv.eq.0) then
c     not found at least one charged lepton and one neutral lepton
         write(*,*) 'WARNING: nl=',nl,' nv=',nv
         return
      endif

c     find highest-pt charged lepton
      ptl=0d0
      ihepl=0
      do i=1,nl
         tmp=sqrt(p_l(1,i)**2+sqrt(p_l(2,i)**2))
         if (tmp.gt.ptl) then
            ptl=tmp
            ihepl=iheparr_l(i)
         endif         
      enddo
c     now ihepl contains the position of the highest-pt charged lepton
      if (ihepl.ne.0) then
c     change status of this lepton and save momentum in pl
         isthep(ihepl)=invalidisthep
         do mu=1,4
            pl(mu)=phep(mu,ihepl)
         enddo  
      endif
      

c     find highest-pt neutral lepton
      ptv=0d0
      ihepv=0
      do i=1,nv
         tmp=sqrt(p_v(1,i)**2+sqrt(p_v(2,i)**2))
         if (tmp.gt.ptv) then
            ptv=tmp
            ihepv=iheparr_v(i)
         endif         
      enddo
c     now ihepv contains the position of the highest-pt neutral lepton
      if (ihepv.ne.0) then
c     change status of this lepton and save momentum in pv
         isthep(ihepv)=invalidisthep
         do mu=1,4
            pv(mu)=phep(mu,ihepv)
         enddo  
      endif
      
           

c     find highest-pt W
      ptw=0d0
      ihepw=0
      do i=1,nw
         tmp=sqrt(p_w(1,i)**2+sqrt(p_w(2,i)**2))
         if (tmp.gt.ptw) then
            ptw=tmp
            ihepw=iheparr_w(i)
         endif         
      enddo
c     now ihepw contains the position of the highest-pt W
      if (ihepw.ne.0) then
c     change status of this lepton and save momentum in pw
         isthep(ihepw)=invalidisthep
         do mu=1,4
            pw(mu)=phep(mu,ihepw)
         enddo  
      endif

c      do mu=1,4
c         if (phep(mu,ihepw)-phep(mu,ihepl)-phep(mu,ihepv).gt.1d-6) then
c            write(*,*) '================='
c         endif
c      enddo

c      write(*,*) (phep(mu,ihepw)-phep(mu,ihepl)-phep(mu,ihepv),mu=1,4)
     

      process = "antikt"
      call buildjets(process,njets,pjet,jetvec)

c      write(*,*) 'total # jets',njets
c      write(*,*) (sqrt(pjet(1,i)**2+pjet(2,i)**2),i=1,njets)


c     find in which ptrack the B hadrons ended up
      nbjet=0
      nbbarjet=0
c     loop over tracks
      do i=1,ntracks
         id=idhep(ihep_of_track(i))
         if (is_B_hadron(id)) then
            nbjet=nbjet+1
            nbjet_array(nbjet)=jetvec(i)            
         elseif (is_BBAR_hadron(id)) then   
            nbbarjet=nbbarjet+1
            nbbarjet_array(nbbarjet)=jetvec(i)                        
         endif
      enddo

      if (nbjet*nbbarjet.eq.0) return

c     Now there are nbjet b-jets and nbbarjet bbar-jets. The ith b jet
c     is the nbjet_array(i)th jet. The ith bbar jet is the
c     nbbarjet_array(i)th jet
      
c      do i=1, nbjet
c         write(*,*) 'b jets'
c         write(*,*) nbjet_array(i)
c      enddo
c       do i=1, nbbarjet
c         write(*,*) 'bbar jets'
c         write(*,*) nbbarjet_array(i)
c      enddo
 
c     find the hardest and next-to-hardest b-flavored jet (both b or
c     bbar) and hardest non-b jet

c     jets are ordered in decreasing pt. Set up array of info on jets
c     if jetinfo=0 then non-b jet
c     if jetinfo=5 then b jet
c     if jetinfo=-5 then bbar jet
      do i=1,njets
          jetinfo(i)=0
       enddo
      do i=1,njets
         do j=1,nbjet
            if (i.eq.nbjet_array(j)) then
               jetinfo(i)=5
            endif
         enddo
         do j=1,nbbarjet
            if (i.eq.nbbarjet_array(j)) then
               jetinfo(i)=-5
            endif
         enddo
      enddo
c      write(*,*) 'nbjet nbbarjet', nbjet,nbbarjet
c      write(*,*) 'positions of the first twos',
c     $     nbjet_array(1),nbbarjet_array(1)
c      write(*,*) (jetinfo(i),i=1,njets)

      found_hardjet=.false.
      found_bjet1=.false.
      found_bjet2=.false.
      typeb1=0
      typeb2=0
      do i=1,njets
         if (jetinfo(i).eq.0.and..not.found_hardjet) then
            found_hardjet=.true.
            do mu=1,4
               phardjet(mu)=pjet(mu,i)
            enddo            
         elseif (abs(jetinfo(i)).eq.5) then
            if (.not.found_bjet1) then
               found_bjet1=.true.
               do mu=1,4
                  pbjet1(mu)=pjet(mu,i)
               enddo 
c     keep track of b flavor of the 1st jet (quark or antiquark)
               typeb1=jetinfo(i)
               goto 111
            endif
            if (.not.found_bjet2.and.found_bjet1) then
               found_bjet2=.true.
               do mu=1,4
                  pbjet2(mu)=pjet(mu,i)
               enddo 
c     keep track of b flavor of the 2nd jet (quark or antiquark)
               typeb2=jetinfo(i)
            endif
            if (found_bjet1.and.found_bjet2) then
c     they must come from a b-bbar couple. Otherwise, return
               if (typeb1*typeb2.gt.0) then
                  wrong_bb_sequence = wrong_bb_sequence + 1                  
                  if (mod(wrong_bb_sequence,200).eq.0) then
                     write(*,*) 'WARNING: 2 b or 2 bbar in sequence ', 
     $                    wrong_bb_sequence
                  endif
c                  write(*,*) 'nbjet nbbarjet', nbjet,nbbarjet
c                  write(*,*) (jetinfo(i),i=1,njets)
                  return
               endif
            endif
         endif
 111     continue
      enddo


c     if there is only one b jet, then return
      if (.not.(found_bjet1.and.found_bjet2)) then
         return
      endif
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     The analysis starts here!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      

      call pwhg_getpt(pW,ptW)
      call pwhg_getpseudorapidity(pW,etaW)
      if (abs(etaW).gt.etaWcut) return
      call pwhg_getinvmass(pW,invmW)
      call pwhg_getpt(pbjet1,ptb1)
      if (ptb1.lt.ptbcut) return
      call pwhg_getpt(pbjet2,ptb2)
      if (ptb2.lt.ptbcut) return
      call pwhg_getpseudorapidity(pbjet1,etab1)
      if (abs(etab1).gt.etabcut) return
      call pwhg_getpseudorapidity(pbjet2,etab2)
      if (abs(etab2).gt.etabcut) return
      do mu=1,4
         pbbjet(mu)= pbjet1(mu)+pbjet2(mu)
      enddo
      call pwhg_getinvmass(pbbjet,mbb)
      call pwhg_getR_phiy(pbjet1,pbjet2,Rbb)

      call pwhg_getdelta_azi(pbbjet,pw,dphi_bbW)
      call pwhg_getpseudorapidity(pbbjet,etabb)
      deta_bbW=etabb-etaW
      
      deta_bb=etab1-etab2
      call pwhg_getdelta_azi(pbjet1,pbjet2,dphi_bb)

      if (Wdecay) then
         call pwhg_getpt(pv,ptv)
         if (ptv.lt.ptvcut) return
         call pwhg_getpt(pl,ptl)
         if (ptl.lt.ptlcut) return
         call pwhg_getpseudorapidity(pl,etal)
         if (abs(etal).gt.etalcut) return
      endif
      
c      ptj=0d0
c      yj=1d8
c      ybbW=1d8
      do mu=1,4
         pbbjetW(mu)= pbbjet(mu) + pW(mu)
      enddo         
      if (found_hardjet) then
         call pwhg_getpt(phardjet,ptj)
         call pwhg_getrapidity(phardjet,yj)
         call pwhg_getrapidity(pbbjetW,ybbW)
         dy_bbWj=ybbW-yj
      endif
      call pwhg_getinvmass(pbbjetW,mbbW)
      
      do mu=1,4
         ptot(mu) = pin1(mu) + pin2(mu)
      enddo
      call pwhg_getinvmass(ptot,sqrshat)


      do ncut=0,nptWcut-1     

      if (ptW.gt.ptWcuts(ncut)) then
         
      diag=1
      call pwhgfill(diag+numplots*ncut,0.5d0,dsig/bsz(diag))
      diag=2
      call pwhgfill(diag+numplots*ncut,ptb1,dsig/bsz(diag))
      diag=3
      call pwhgfill(diag+numplots*ncut,ptb2,dsig/bsz(diag))
      diag=4
      call pwhgfill(diag+numplots*ncut,etab1,dsig/bsz(diag))
      diag=5
      call pwhgfill(diag+numplots*ncut,etab2,dsig/bsz(diag))
      diag=6
      call pwhgfill(diag+numplots*ncut,ptW,dsig/bsz(diag))
      diag=7
      call pwhgfill(diag+numplots*ncut,etaW,dsig/bsz(diag))
      diag=8
      call pwhgfill(diag+numplots*ncut,mbb,dsig/bsz(diag))
      diag=9
      call pwhgfill(diag+numplots*ncut,Rbb,dsig/bsz(diag))
      diag=10
      call pwhgfill(diag+numplots*ncut,dphi_bbW,dsig/bsz(diag))
      diag=11
      call pwhgfill(diag+numplots*ncut,deta_bbW,dsig/bsz(diag))

      if (Wdecay) then
         diag=12
         call pwhgfill(diag+numplots*ncut,ptl,dsig/bsz(diag))
         diag=13
         call pwhgfill(diag+numplots*ncut,etal,dsig/bsz(diag))
      endif
      
      if (.not.found_hardjet) then
         diag=31
         call pwhgfill(diag+numplots*ncut,mbbW,dsig/bsz(diag))
      endif


      if (found_hardjet) then
         diag=14
         call pwhgfill(diag+numplots*ncut,ptj,dsig/bsz(diag))
         
         if (idhep(1).eq.21.or.idhep(2).eq.21) then
            diag=23
            call pwhgfill(diag+numplots*ncut,ptj,dsig/bsz(diag))
         else
            diag=24
            call pwhgfill(diag+numplots*ncut,ptj,dsig/bsz(diag))
         endif

         if (abs(yj).lt.yjcut) then
            diag=15
            call pwhgfill(diag+numplots*ncut,ptj,dsig/bsz(diag))
            diag=22
            call pwhgfill(diag+numplots*ncut,ptj,dsig/bsz(diag))
         endif
         diag=16
         call pwhgfill(diag+numplots*ncut,yj,dsig/bsz(diag))
         if (ptj.gt.ptjcut) then
            diag=17
            call pwhgfill(diag+numplots*ncut,yj,dsig/bsz(diag))
         endif
         diag=18
         call pwhgfill(diag+numplots*ncut,dy_bbWj,dsig/bsz(diag))
         
         diag=26
         if (ptj.gt.0d0) then
            call pwhgfill(diag+numplots*ncut,mbbW,dsig/bsz(diag))
         endif
         diag=27
         if (ptj.gt.100d0) then
            call pwhgfill(diag+numplots*ncut,mbbW,dsig/bsz(diag))
         endif
         diag=28
         if (ptj.gt.200d0) then
            call pwhgfill(diag+numplots*ncut,mbbW,dsig/bsz(diag))
         endif
         diag=29
         if (ptj.gt.300d0) then
            call pwhgfill(diag+numplots*ncut,mbbW,dsig/bsz(diag))
         endif
         diag=30
         if (ptj.gt.400d0) then
            call pwhgfill(diag+numplots*ncut,mbbW,dsig/bsz(diag))
         endif
      endif
               
      diag=19
      call pwhgfill(diag+numplots*ncut,invmW,dsig/bsz(diag))
      
      diag=20
      call pwhgfill(diag+numplots*ncut,deta_bb,dsig/bsz(diag))
      
      diag=21
      call pwhgfill(diag+numplots*ncut,dphi_bb,dsig/bsz(diag))
      
      diag=25
      call pwhgfill(diag+numplots*ncut,sqrshat,dsig/bsz(diag))

      if (rad_type.eq.1) then
         diag=32
         call pwhgfill(diag+numplots*ncut,sqrshat,dsig/bsz(diag))
      elseif (rad_type.eq.2) then
         diag=33
         call pwhgfill(diag+numplots*ncut,sqrshat,dsig/bsz(diag))
      else
         write(*,*) 'Unexpected rad_type value: ',rad_type
      endif
      

      if ((WHCPRG.eq.'NLO   ').or.(WHCPRG.eq.'LHE   ')) then
         diag=34
         do mu=1,4
            ptot(mu) = phep(mu,3)+phep(mu,4)+phep(mu,5)
         enddo
         call pwhg_getinvmass(ptot,mbbW)
         call pwhgfill(diag+numplots*ncut,mbbW,dsig/bsz(diag))
      endif



      endif
      enddo

      end



      subroutine buildjets(process,njets,pjet,jetvec)
c     arrays to reconstruct jets
      implicit none
      include 'hepevt.h'
      integer njets
      real * 8 pjet(4,*)
      character * 20 process
      integer maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real *8 ptrack(4,maxtrack)
      real *8 r,ptmin_fastkt
      integer jetvec(maxtrack),ihep_of_track(maxtrack)
      integer ihep,j1,ntracks,jpart,jjet,mu
      real * 8 found
      real * 8 f,caf,sf
      integer seed
      data seed/1/
      save seed
      common/cptrack/ihep_of_track,ptrack,ntracks
      common/cptmin_fastkt/ptmin_fastkt
      logical ini
      save ini
      data ini/.true./

      if (ini) then
         write(*,*) '*******************************'
         write(*,*) '   ptmin_fastkt = ',ptmin_fastkt
         write(*,*) '*******************************'
         ini=.false.
      endif
c     get valid tracks
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
c     "stable" particles (particles that we detect)
         if (isthep(ihep).eq.1) then
            if(ntracks.eq.maxtrack) then
               write(*,*)
     #              'analyze: too many particles, increase maxtrack'
               stop
            endif
c     copy momenta to construct jets 
            ntracks=ntracks+1
            ihep_of_track(ntracks)=ihep
            do mu=1,4
               ptrack(mu,ntracks)=phep(mu,ihep)
            enddo
         endif
      enddo
      if (ntracks.eq.0) then
         njets=0
         return
      endif
c     siscone algorithm
c*********************************************************************
c      R = 0.7  radius parameter
c      f = 0.5  overlapping fraction
c.....run the clustering
c      call fastjetsiscone(ptrack,ntracks,0.7d0,0.5d0,pjet,njets) 
c*********************************************************************
c     fastkt algorithm
c*********************************************************************
c      R = 0.7  Radius parameter
c.....run the clustering 
c      R = 0.5d0          
c      ptmin_fastkt = 0d0
c      call fastjetktwhich(ptrack,ntracks,ptmin_fastkt,R,
c     #     pjet,njets,jetvec)
c     now we have the jets

c      call fastjetd0runiicone(ptrack,ntracks,0.7d0,6d0,0.5d0,pjet,njets)
c      return

c     kt algo
c      R=0.4
c      f=0.75
c      sf=1
c      caf=1
c      call fastjetcdfmidpoint(ptrack,ntracks,r,f,sf,caf,pjet,njets)

c      if     (process.eq."antikt R04") then
c         R=0.4d0 
c      elseif (process.eq."antikt R06") then
c         R=0.6d0
c      elseif (process.eq."antikt R02") then
c         R=0.2d0
c      elseif (process.eq."antikt R08") then
c         R=0.8d0
c      else
c         write(*,*) 'JET ANALYSIS TO USE UNKNOWN:',process
c         call exit(1)
c      endif
      R=0.4d0
      if (process.eq."antikt") then
         call fastjetantikt(ptrack,ntracks,ptmin_fastkt,R,
     $        pjet,njets,jetvec)
      elseif (process.eq."kt") then
         call fastjetktwhich(ptrack,ntracks,ptmin_fastkt,R,
     $        pjet,njets,jetvec)
      else
         write(*,*) 'JET ANALYSIS TO USE UNKNOWN:',process
         call exit(1)
      endif
      end





      function is_B_hadron(id)
      implicit none
      logical is_B_hadron
      integer id
      is_B_hadron=((id.gt.-600).and.(id.lt.-500)).or.
     $     ((id.gt.5000).and.(id.lt.6000)).or.(id.eq.5)
      end

      function is_BBAR_hadron(id)
      implicit none
      logical is_BBAR_hadron
      integer id
      is_BBAR_hadron=((id.gt.500).and.(id.lt.600)).or.
     $     ((id.gt.-6000).and.(id.lt.-5000)).or.(id.eq.-5)
      end


      
