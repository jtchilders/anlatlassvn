c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none
      include '../pwhg_book.h'
      real * 8 pi,pi2
      parameter(pi = 3.141592653589793D0, pi2 = 9.869604401089358D0)
      integer nptWcutmax,nptWcut
      parameter(nptWcutmax = 10)
      real * 8 ptWcuts(0:nptWcutmax-1)
      common/cptvbcut/ptWcuts,nptWcut,numplots
      integer ncut,numplots
      
      character * 10 cut
      integer i
c     binsize
      real * 8 bsz(100)
      common/pwhghistcommon/bsz

      
      ptWcuts(0) = 0d0
      ptWcuts(1) = 5d0
      ptWcuts(2) = 10d0
      ptWcuts(3) = 20d0
c      ptWcuts(4) = 30d0
c     number of pt W cuts implemented
      nptWcut=4

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      numplots = 69  ! <========== DO NOT FORGET TO SET THIS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call pwhginihist

      

c     loop on ptW cut
      do ncut=0,nptWcut-1
      write(unit=cut,fmt="(f5.2)") ptWcuts(ncut)
      bsz(1) = 20d0      
      call pwhgbookup(1+numplots*ncut,'pt W ptW>'//cut,'LOG',
     #     bsz(1),0d0,800d0)
      bsz(2) = 20d0      
      call pwhgbookup(2+numplots*ncut,'pt J1 ptW>'//cut,'LOG',
     #     bsz(2),0d0,800d0)
      bsz(3) = 20d0
      call pwhgbookup(3+numplots*ncut,'pt J2 ptW>'//cut,'LOG',
     #     bsz(3),0d0,800d0)
      bsz(4) = 1d0
      call pwhgbookup(4+numplots*ncut,'inv mass W ptW>'//cut,'LOG',
     #     bsz(4),60d0,120d0)

      bsz(5) = 0.4d0
      call pwhgbookup(5+numplots*ncut,'y W, ptW>10 ptW>'//cut,'LOG',
     #     bsz(5),-5d0,5d0)
      bsz(6) = 0.4d0
      call pwhgbookup(6+numplots*ncut,'y W, ptW>20 ptW>'//cut,'LOG',
     #     bsz(6),-5d0,5d0)
      bsz(7) = 0.4d0
      call pwhgbookup(7+numplots*ncut,'y W, ptW>40 ptW>'//cut,'LOG',
     #     bsz(7),-5d0,5d0)
      bsz(8) = 0.4d0
      call pwhgbookup(8+numplots*ncut,'y W, ptW>60 ptW>'//cut,'LOG',
     #     bsz(8),-5d0,5d0)
      bsz(9) = 0.4d0
      call pwhgbookup(9+numplots*ncut,'y W, ptW>80 ptW>'//cut,'LOG',
     #     bsz(9),-5d0,5d0)
      bsz(10) = 0.4d0
      call pwhgbookup(10+numplots*ncut,'y W, ptW>100 ptW>'//cut,'LOG',
     #     bsz(10),-5d0,5d0)
      bsz(11) = 0.4d0
      call pwhgbookup(11+numplots*ncut,'y W, ptW>150 ptW>'//cut,'LOG',
     #     bsz(11),-5d0,5d0)
      bsz(12) = 0.4d0
      call pwhgbookup(12+numplots*ncut,'y W, ptW>200 ptW>'//cut,'LOG',
     #     bsz(12),-5d0,5d0)
      bsz(13) = 0.4d0
      call pwhgbookup(13+numplots*ncut,'y W, ptW>300 ptW>'//cut,'LOG',
     #     bsz(13),-5d0,5d0)

      bsz(14) = 0.4d0 
      call pwhgbookup(14+numplots*ncut,'yJ1, ptJ1>10 ptW>'//cut,'LOG',
     #     bsz(14),-5d0,5d0)
      bsz(15) = 0.4d0
      call pwhgbookup(15+numplots*ncut,'yJ1, ptJ1>20 ptW>'//cut,'LOG',
     #     bsz(15),-5d0,5d0)
      bsz(16) = 0.4d0
      call pwhgbookup(16+numplots*ncut,'yJ1, ptJ1>40 ptW>'//cut,'LOG',
     #     bsz(16),-5d0,5d0)
      bsz(17) = 0.4d0
      call pwhgbookup(17+numplots*ncut,'yJ1, ptJ1>60 ptW>'//cut,'LOG',
     #     bsz(17),-5d0,5d0)
      bsz(18) = 0.4d0
      call pwhgbookup(18+numplots*ncut,'yJ1, ptJ1>80 ptW>'//cut,'LOG',
     #     bsz(18),-5d0,5d0)
      bsz(19) = 0.4d0
      call pwhgbookup(19+numplots*ncut,'yJ1, ptJ1>100 ptW>'//cut,'LOG',
     #     bsz(19),-5d0,5d0)
      bsz(20) = 0.4d0
      call pwhgbookup(20+numplots*ncut,'yJ1, ptJ1>150 ptW>'//cut,'LOG',
     #     bsz(20),-5d0,5d0)
      bsz(21) = 0.4d0
      call pwhgbookup(21+numplots*ncut,'yJ1, ptJ1>200 ptW>'//cut,'LOG',
     #     bsz(21),-5d0,5d0)
      bsz(22) = 0.4d0
      call pwhgbookup(22+numplots*ncut,'yJ1, ptJ1>300 ptW>'//cut,'LOG',
     #     bsz(22),-5d0,5d0)

      bsz(23) = 0.4d0
      call pwhgbookup(23+numplots*ncut,'dy WJ1-j2,ptJ2>10 ptW>'//cut,
     #     'LOG',bsz(23),-10d0,10d0)
      bsz(24) = 0.4d0
      call pwhgbookup(24+numplots*ncut,'dy WJ1-j2, ptJ2>20 ptW>'//cut,
     #     'LOG',bsz(24),-10d0,10d0)
      bsz(25) = 0.4d0
      call pwhgbookup(25+numplots*ncut,'dy WJ1-j2, ptJ2>40 ptW>'//cut,
     #     'LOG',bsz(25),-10d0,10d0)
      bsz(26) = 0.4d0
      call pwhgbookup(26+numplots*ncut,'dy WJ1-j2, ptJ2>60 ptW>'//cut,
     #     'LOG',bsz(26),-10d0,10d0)
      bsz(27) = 0.4d0
      call pwhgbookup(27+numplots*ncut,'dy WJ1-j2, ptJ2>80 ptW>'//cut,
     #     'LOG',bsz(27),-10d0,10d0)
      bsz(28) = 0.4d0
      call pwhgbookup(28+numplots*ncut,'dy WJ1-j2, ptJ2>100 ptW>'//cut,
     #     'LOG',bsz(28),-10d0,10d0)
      bsz(29) = 0.4d0
      call pwhgbookup(29+numplots*ncut,'dy WJ1-j2, ptJ2>150 ptW>'//cut,
     #     'LOG',bsz(29),-10d0,10d0)
      bsz(30) = 0.4d0
      call pwhgbookup(30+numplots*ncut,'dy WJ1-j2, ptJ2>200 ptW>'//cut,
     #     'LOG',bsz(30),-10d0,10d0)
      bsz(31) = 0.4d0
      call pwhgbookup(31+numplots*ncut,'dy WJ1-j2, ptJ2>300 ptW>'//cut,
     #     'LOG',bsz(31),-10d0,10d0)

      bsz(32) = 20d0
      call pwhgbookup(32+numplots*ncut,'pt lep ptW>'//cut,'LOG',
     #     bsz(32),0d0,800d0)
      bsz(33) = 20d0
      call pwhgbookup(33+numplots*ncut,'pt vl ptW>'//cut,'LOG',
     #     bsz(33),0d0,800d0)
      bsz(34) = 0.4d0
      call pwhgbookup(34+numplots*ncut,'y lep, pt lep >10 ptW>'//cut,
     #     'LOG',bsz(34),-5d0,5d0)
      bsz(35) = 0.4d0
      call pwhgbookup(35+numplots*ncut,'y lep, pt lep >20 ptW>'//cut,
     #     'LOG', bsz(35),-5d0,5d0)
      bsz(36) = 0.4d0
      call pwhgbookup(36+numplots*ncut,'y lep, pt lep >40 ptW>'//cut,
     #     'LOG',bsz(36),-5d0,5d0)
      bsz(37) = 0.4d0
      call pwhgbookup(37+numplots*ncut,'y lep, pt lep >60 ptW>'//cut,
     #     'LOG',bsz(37),-5d0,5d0)
      bsz(38) = 0.4d0
      call pwhgbookup(38+numplots*ncut,'y lep, pt lep >80 ptW>'//cut,
     #     'LOG',bsz(38),-5d0,5d0)
      bsz(39) = 0.4d0
      call pwhgbookup(39+numplots*ncut,'y lep, pt lep >100 ptW>'//cut,
     #     'LOG',bsz(39),-5d0,5d0)

      bsz(40) = 0.4d0
      call pwhgbookup(40+numplots*ncut,'y vl, pt vl >10 ptW>'//cut,
     #     'LOG',bsz(40),-5d0,5d0)
      bsz(41) = 0.4d0
      call pwhgbookup(41+numplots*ncut,'y vl, pt vl >20 ptW>'//cut,
     #     'LOG',bsz(41),-5d0,5d0)
      bsz(42) = 0.4d0
      call pwhgbookup(42+numplots*ncut,'y vl, pt vl >40 ptW>'//cut,
     #     'LOG',bsz(42),-5d0,5d0)
      bsz(43) = 0.4d0
      call pwhgbookup(43+numplots*ncut,'y vl, pt vl >60 ptW>'//cut,
     #     'LOG',bsz(43),-5d0,5d0)
      bsz(44) = 0.4d0
      call pwhgbookup(44+numplots*ncut,'y vl, pt vl >80 ptW>'//cut,
     #     'LOG',bsz(44),-5d0,5d0)
      bsz(45) = 0.4d0
      call pwhgbookup(45+numplots*ncut,'y vl, pt vl >100 ptW>'//cut,
     #     'LOG',bsz(45),-5d0,5d0)     

      bsz(46) = 0.4d0
      call pwhgbookup(46+numplots*ncut,'y lep, ptW>'//cut,'LOG',
     #     bsz(46),-5d0,5d0)
      bsz(47) = 0.4d0
      call pwhgbookup(47+numplots*ncut,'y vl, ptW>'//cut,'LOG',
     #     bsz(47),-5d0,5d0)

      bsz(48) = 4d0
      call pwhgbookup(48+numplots*ncut,'pt lep, zoom ptW>'//cut,'LOG',
     #     bsz(48),0d0,100d0)
      bsz(49) = 4d0
      call pwhgbookup(49+numplots*ncut,'pt vl, zoom ptW>'//cut,'LOG',
     #     bsz(49),0d0,100d0)
      bsz(50) = 2d0
      call pwhgbookup(50+numplots*ncut,'pt W, zoom ptW>'//cut,'LOG',
     #     bsz(50),0d0,100d0)
      bsz(51) = 0.5d0
      call pwhgbookup(51+numplots*ncut,'pt W, zoom2 ptW>'//cut,'LOG',
     #     bsz(51),0d0,30d0)

      bsz(52) = 0.5d0
      call pwhgbookup(52+numplots*ncut,'pt J1, zoom ptW>'//cut,'LOG',
     #     bsz(52),0d0,30d0)
      bsz(53) = 0.5d0
      call pwhgbookup(53+numplots*ncut,'pt J2, zoom ptW>'//cut,'LOG',
     #     bsz(53),0d0,30d0)

      bsz(54) = 0.5d0
      call pwhgbookup(54+numplots*ncut,'pt_rel J1 ptW>'//cut,'LOG',
     #     bsz(54),0d0,15d0)
      bsz(55) = 0.5d0
      call pwhgbookup(55+numplots*ncut,'pt_rel J2 ptW>'//cut,'LOG',
     #     bsz(55),0d0,15d0)

      bsz(56) = 0.2d0
      call pwhgbookup(56+numplots*ncut,'azimuth W ptW>'//cut,'LIN',
     #     bsz(56),-3.2d0,3.2d0)    
      bsz(57) = 1d0
      call pwhgbookup(57+numplots*ncut,'total Xsec ptW>'//cut,'LOG',
     #     bsz(57),0d0,1d0)

      bsz(58) = 0.2d0
      call pwhgbookup(58+numplots*ncut,'deltaR_12 ptW>'//cut,'LIN',
     #     bsz(58),0d0,10d0)


      bsz(59) = 3d0
      call pwhgbookup(59+numplots*ncut,'pt J1, Rjj>1, ptW>'
     #     //cut,'LOG',bsz(59),0d0,300d0)
      bsz(60) = 3d0
      call pwhgbookup(60+numplots*ncut,'pt J1, Rjj>2, ptW>'
     #     //cut,'LOG',bsz(60),0d0,300d0)
      bsz(61) = 3d0
      call pwhgbookup(61+numplots*ncut,'pt J1, Rjj>3, ptW>'
     #     //cut,'LOG',bsz(61),0d0,300d0)
      bsz(62) = 3d0
      call pwhgbookup(62+numplots*ncut,'pt J1, Rjj>4, ptW>'
     #     //cut,'LOG',bsz(62),0d0,300d0)
      bsz(63) = 3d0
      call pwhgbookup(63+numplots*ncut,'pt J1, Rjj>5, ptW>'
     #     //cut,'LOG',bsz(63),0d0,300d0)

      bsz(64) = 3d0
      call pwhgbookup(64+numplots*ncut,'pt J2, Rjj>1, ptW>'
     #     //cut,'LOG',bsz(64),0d0,300d0)
      bsz(65) = 3d0
      call pwhgbookup(65+numplots*ncut,'pt J2, Rjj>2, ptW>'
     #     //cut,'LOG',bsz(65),0d0,300d0)
      bsz(66) = 3d0
      call pwhgbookup(66+numplots*ncut,'pt J2, Rjj>3, ptW>'
     #     //cut,'LOG',bsz(66),0d0,300d0)
      bsz(67) = 3d0
      call pwhgbookup(67+numplots*ncut,'pt J2, Rjj>4, ptW>'
     #     //cut,'LOG',bsz(67),0d0,300d0)
      bsz(68) = 3d0
      call pwhgbookup(68+numplots*ncut,'pt J2, Rjj>5, ptW>'
     #     //cut,'LOG',bsz(68),0d0,300d0)

      bsz(69) = 2d0
      call pwhgbookup(69+numplots*ncut,'pt j1, zoom2 ptW>'//cut,'LOG',
     #     bsz(69),0d0,100d0)
      enddo
      end

      
     
      subroutine analysis(dsig)
      implicit none
      real * 8 dsig
      include '../include/hepevt.h'
      include '../include/LesHouches.h'
c arrays to reconstruct jets
      integer maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real *8 ptrack(4,maxtrack)
      real *8 ptj1,ptj2,yj1,yj2,ptvb,yvb,yvbj1
      real *8 ptvl,ptlep,yvl,ylep
      real *8 pjet(4,maxjet) 
      real * 8 mvb,pvb(4),pvbj1(4),tmp

      integer nptWcutmax,nptWcut,numplots
      parameter(nptWcutmax = 10)
      real * 8 ptWcuts(0:nptWcutmax-1)
      common/cptvbcut/ptWcuts,nptWcut,numplots
      integer ncut

      integer mu,jpart,jjet,j1,j2,found,njets,
     1     nvl,nlep,ihep,ntracks,ijet
      logical buildjets
      parameter (buildjets=.true.)

      real * 8 vec(3),pjetin(0:3),pjetout(0:3),beta,ptrel,get_ptrel,
     #     ptrackin(0:3),ptrackout(0:3)
      integer i
      external get_ptrel
      real * 8 R,ptmin_fastkt
      integer jetvec(maxtrack),jj
      logical ini
      data ini/.true./
      save ini
      integer maxnumlep
      parameter (maxnumlep=10)
      integer lepvec(maxnumlep),vlvec(maxnumlep),ivl,ilep,vl,lep
      integer nleps,nvls
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data   WHCPRG/'NLO   '/
      logical is_W
      integer nfoundjets
      integer maxjets
      parameter (maxjets=10)
      integer njj(maxjets)      
      real *8 ptj(maxjets),yj(maxjets),pj(0:3,maxjets)
c     binsize
      real * 8 bsz(100)
      common/pwhghistcommon/bsz
      real * 8 getrapidity0
      external getrapidity0
      real * 8 rsep,rsepn_p
      external rsepn_p
      logical pass_lept_cuts
      external pass_lept_cuts
      real * 8 dist
      real * 8 yl1,mV2
      real * 8 pvl(4,maxnumlep),plep(4,maxnumlep)
      common/clepton_cuts/yl1
      integer vdecaytemp,vdecay2temp,idvecbos
      save vdecaytemp,vdecay2temp,idvecbos

      if (ini) then
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
         
!     id of the charged decay product of the W
         vdecaytemp=lprup(1)-10000
         if(vdecaytemp.lt.0) then 
            vdecay2temp=-vdecaytemp+1 
            idvecbos=24
         elseif(vdecaytemp.gt.0) then
            vdecay2temp=-(vdecaytemp+1)
            idvecbos=-24
         else
            write(*,*) 'Error in decay mode in pwhg_analysis'
            call exit(1)
         endif
         if(abs(vdecaytemp).eq.11.or.abs(vdecaytemp).eq.13
     $        .or.abs(vdecaytemp).eq.15) then
            continue
         else
            write(*,*) '**************************************'
            write(*,*) ' Analysis works only for e, mu or tau decays'
            write(*,*) '                 STOP     '
            write(*,*) '**************************************'
            call exit(1)
         endif
          write (*,*) '********************************************'
         write (*,*)
         ini=.false.
      endif

      if ((WHCPRG.eq.'NLO   ').or.(WHCPRG.eq.'LHE   ')) then 
         nleps=0
         nvls=0
         do i=1,maxnumlep
            lepvec(i) = 0
            vlvec(i) = 0
         enddo
         do ihep=1,nhep
            if (isthep(ihep).eq.1) then
               if(idhep(ihep).eq.vdecaytemp) then
                  nleps=nleps+1
                  lepvec(nleps)=ihep
               elseif(idhep(ihep).eq.vdecay2temp) then
                  nvls=nvls+1
                  vlvec(nvls)=ihep
               endif
            endif         
         enddo
         
         if (nleps.ne.1.or.nvls.ne.1) then
            write(*,*) "Too many leptons found. PROGRAM ABORT"
            call exit(1)
         else 
            ilep=lepvec(1)
            ivl=vlvec(1)
         endif
      else
         return
      endif

      
      do i=1,maxnumlep
         lepvec(i) = 0
         vlvec(i) = 0
      enddo

c     W momentum
      do mu=1,4
         pvb(mu)=phep(mu,ilep)+phep(mu,ivl)         
      enddo
      mV2 = pvb(4)**2-pvb(1)**2-pvb(2)**2-pvb(3)**2
c      write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>',sqrt(mV2)

      ptvb=sqrt(pvb(1)**2+pvb(2)**2)

      call getinvmass(pvb,mvb)
      call getrapidity(pvb,yvb)
      ptvl=sqrt(phep(1,ivl)**2+phep(2,ivl)**2)
      call getrapidity(phep(1,ivl),yvl)
      ptlep=sqrt(phep(1,ilep)**2+phep(2,ilep)**2)
      call getrapidity(phep(1,ilep),ylep)

c$$$      call build_jets(njets,pjet)

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

      found=0
      ntracks=0
      njets=0
c     Loop over final state particles to find jets 
      do ihep=1,nhep
         if ((isthep(ihep).eq.1).and.
c     exclude leptons, gauge and higgs bosons 
     1        (((abs(idhep(ihep)).le.10).or.(abs(idhep(ihep)).ge.40))
c     but include gluons
     2        .or.(abs(idhep(ihep)).eq.21))) then
           if (ntracks.eq.maxtrack) then
              write(*,*) 'Too many particles. Increase maxtrack.'//
     #             ' PROGRAM ABORTS'
              call exit(1)
           endif
c     copy momenta to construct jets 
           ntracks=ntracks+1
           do mu=1,4
              ptrack(mu,ntracks)=phep(mu,ihep)
           enddo
        endif
      enddo
      
       if(buildjets.and.ntracks.gt.0) then
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
         ptmin_fastkt = 20d0
         call fastjetktwhich(ptrack,ntracks,ptmin_fastkt,R,
     #        pjet,njets,jetvec) 
c     
c     ... now we have the jets
         if (njets.gt.0) then
c     find the first 2 hardest jets, if any
            call find_hardest_jets(njets,pjet,2,nfoundjets,njj)
         else
            dsig=0
            return
         endif
      endif


c     loop on ptW cut
      do ncut=0,nptWcut-1


      if (ptvb.lt.ptWcuts(ncut)) goto 987

      call pwhgfill( 1+numplots*ncut,ptvb,dsig/bsz(1))
      call pwhgfill(50+numplots*ncut,ptvb,dsig/bsz(50))
      call pwhgfill(51+numplots*ncut,ptvb,dsig/bsz(51))
      
      call pwhgfill(4+numplots*ncut,mvb,dsig/bsz(4))

      if(buildjets.and.ntracks.gt.0) then
         if (njets.gt.0) then
            do ijet=1,nfoundjets
               do mu=1,3
                  pj(mu,ijet)=pjet(mu,njj(ijet))
               enddo
               pj(0,ijet)=pjet(4,njj(ijet))
            enddo
            
c     get pt's and rapidities of the jets
            do ijet=1,nfoundjets
               ptj(ijet) = sqrt(pj(1,ijet)**2 + pj(2,ijet)**2)
               yj(ijet) = getrapidity0(pj(0,ijet))
            enddo

            j1 = 0
            j2 = 0            
            if (nfoundjets.ge.1) then
               j1 = njj(1)
               ptj1 = ptj(1)
               yj1 = yj(1)
               do mu=1,4
                  pvbj1(mu)=pvb(mu)+pjet(mu,j1)
               enddo
               call getrapidity(pvbj1,yvbj1)
            endif
            if (nfoundjets.ge.2) then
               j2 = njj(2)
               ptj2 = ptj(2)
               yj2 = yj(2)
            endif

            if(j1.gt.0) then
               do mu=1,4
                  pvbj1(mu)=pvb(mu)+pjet(mu,j1)
               enddo
               call getrapidity(pvbj1,yvbj1)
            endif
            
            if(j1.gt.0) then
               call pwhgfill( 2+numplots*ncut,ptj1,dsig/bsz(2))
               call pwhgfill(52+numplots*ncut,ptj1,dsig/bsz(52))
               call pwhgfill(69+numplots*ncut,ptj1,dsig/bsz(69))
            endif
            rsep = 0d0
            if(j2.gt.0) then
               call pwhgfill( 3+numplots*ncut,ptj2,dsig/bsz(3))
               call pwhgfill(53+numplots*ncut,ptj2,dsig/bsz(53))
c     compute the separation in the pseudorapidity-phi plane
               rsep = rsepn_p(pj(0,1),pj(0,2))
               call pwhgfill(58+numplots*ncut,rsep,dsig/bsz(58))
               if (rsep.gt.1d0) then
                  call pwhgfill(59+numplots*ncut,ptj1,dsig/bsz(59))
                  call pwhgfill(64+numplots*ncut,ptj2,dsig/bsz(64))
               endif
               if (rsep.gt.2d0) then
                  call pwhgfill(60+numplots*ncut,ptj1,dsig/bsz(60))
                  call pwhgfill(65+numplots*ncut,ptj2,dsig/bsz(65))
               endif
               if (rsep.gt.3d0) then
                  call pwhgfill(61+numplots*ncut,ptj1,dsig/bsz(61))
                  call pwhgfill(66+numplots*ncut,ptj2,dsig/bsz(66))
               endif
               if (rsep.gt.4d0) then
                  call pwhgfill(62+numplots*ncut,ptj1,dsig/bsz(62))
                  call pwhgfill(67+numplots*ncut,ptj2,dsig/bsz(67))
               endif
               if (rsep.gt.5d0) then
                  call pwhgfill(63+numplots*ncut,ptj1,dsig/bsz(63))
                  call pwhgfill(68+numplots*ncut,ptj2,dsig/bsz(68))
               endif
            endif
         
            if(ptvb.gt. 10) call pwhgfill( 5+numplots*ncut,yvb,
     #           dsig/bsz(5))
            if(ptvb.gt. 20) call pwhgfill( 6+numplots*ncut,yvb,
     #           dsig/bsz(6))
            if(ptvb.gt. 40) call pwhgfill( 7+numplots*ncut,yvb,
     #           dsig/bsz(7))
            if(ptvb.gt. 60) call pwhgfill( 8+numplots*ncut,yvb,
     #           dsig/bsz(8))
            if(ptvb.gt. 80) call pwhgfill( 9+numplots*ncut,yvb,
     #           dsig/bsz(9))
            if(ptvb.gt.100) call pwhgfill(10+numplots*ncut,yvb,
     #           dsig/bsz(10))
            if(ptvb.gt.150) call pwhgfill(11+numplots*ncut,yvb,
     #           dsig/bsz(11))
            if(ptvb.gt.200) call pwhgfill(12+numplots*ncut,yvb,
     #           dsig/bsz(12))
            if(ptvb.gt.300) call pwhgfill(13+numplots*ncut,yvb,
     #           dsig/bsz(13))
            
            if(j1.gt.0) then
               if(ptj1.gt. 10) call pwhgfill(14+numplots*ncut,yj1,
     #              dsig/bsz(14))
               if(ptj1.gt. 20) call pwhgfill(15+numplots*ncut,yj1,
     #              dsig/bsz(15))
               if(ptj1.gt. 40) call pwhgfill(16+numplots*ncut,yj1,
     #              dsig/bsz(16))
               if(ptj1.gt. 60) call pwhgfill(17+numplots*ncut,yj1,
     #              dsig/bsz(17))
               if(ptj1.gt. 80) call pwhgfill(18+numplots*ncut,yj1,
     #              dsig/bsz(18))
               if(ptj1.gt.100) call pwhgfill(19+numplots*ncut,yj1,
     #              dsig/bsz(19))
               if(ptj1.gt.150) call pwhgfill(20+numplots*ncut,yj1,
     #     dsig/bsz(20))
               if(ptj1.gt.200) call pwhgfill(21+numplots*ncut,yj1,
     #              dsig/bsz(21))
               if(ptj1.gt.300) call pwhgfill(22+numplots*ncut,yj1,
     #              dsig/bsz(22))
            endif
            if(j2.gt.0) then
               if(ptj2.gt. 10) call pwhgfill(23+numplots*ncut,yvbj1-yj2,
     #              dsig/bsz(23))
               if(ptj2.gt. 20) call pwhgfill(24+numplots*ncut,yvbj1-yj2,
     #              dsig/bsz(24))
               if(ptj2.gt. 40) call pwhgfill(25+numplots*ncut,yvbj1-yj2,
     #              dsig/bsz(25))
               if(ptj2.gt. 60) call pwhgfill(26+numplots*ncut,yvbj1-yj2,
     #              dsig/bsz(26))
               if(ptj2.gt. 80) call pwhgfill(27+numplots*ncut,yvbj1-yj2,
     #              dsig/bsz(27))
               if(ptj2.gt.100) call pwhgfill(28+numplots*ncut,yvbj1-yj2,
     #              dsig/bsz(28))
               if(ptj2.gt.150) call pwhgfill(29+numplots*ncut,yvbj1-yj2,
     #              dsig/bsz(29))
               if(ptj2.gt.200) call pwhgfill(30+numplots*ncut,yvbj1-yj2,
     #              dsig/bsz(30))
               if(ptj2.gt.300) call pwhgfill(31+numplots*ncut,yvbj1-yj2,
     #              dsig/bsz(31))
            endif
            
c     loop on the hardest and next-to-hardest jet
            do ijet=1,min(njets,2)
               if (ijet.eq.1) then
                  jj=j1 
               else 
                  jj=j2 
               endif            
               do mu=1,3
                  pjetin(mu) = pjet(mu,jj)
               enddo
               pjetin(0) = pjet(4,jj)         
               vec(1)=0d0
               vec(2)=0d0
               vec(3)=1d0
               beta = -pjet(3,jj)/pjet(4,jj)
               call mboost(1,vec,beta,pjetin,pjetout)         
c     write(*,*) pjetout
               ptrel = 0
               do i=1,ntracks
                  if (jetvec(i).eq.jj) then
                     do mu=1,3
                        ptrackin(mu) = ptrack(mu,i)
                     enddo
                     ptrackin(0) = ptrack(4,i)
                     call mboost(1,vec,beta,ptrackin,ptrackout) 
                     ptrel = ptrel + get_ptrel(ptrackout,pjetout)
                  endif
               enddo
               if (ijet.eq.1) then 
                  call pwhgfill(54+numplots*ncut,ptrel,dsig/bsz(54))
               else
                  call pwhgfill(55+numplots*ncut,ptrel,dsig/bsz(55))
               endif
            enddo         
         endif            
c     endif buildjets
      endif
      
      call pwhgfill(32+numplots*ncut,ptlep,dsig/bsz(32))
      call pwhgfill(33+numplots*ncut,ptvl,dsig/bsz(33))
      
      if(ptlep.gt. 10) call pwhgfill(34+numplots*ncut,ylep,dsig/bsz(34))
      if(ptlep.gt. 20) call pwhgfill(35+numplots*ncut,ylep,dsig/bsz(35))
      if(ptlep.gt. 40) call pwhgfill(36+numplots*ncut,ylep,dsig/bsz(36))
      if(ptlep.gt. 60) call pwhgfill(37+numplots*ncut,ylep,dsig/bsz(37))
      if(ptlep.gt. 80) call pwhgfill(38+numplots*ncut,ylep,dsig/bsz(38))
      if(ptlep.gt.100) call pwhgfill(39+numplots*ncut,ylep,dsig/bsz(39))
      
      if(ptvl.gt. 10) call pwhgfill(40+numplots*ncut,yvl,dsig/bsz(40))
      if(ptvl.gt. 20) call pwhgfill(41+numplots*ncut,yvl,dsig/bsz(41))
      if(ptvl.gt. 40) call pwhgfill(42+numplots*ncut,yvl,dsig/bsz(42))
      if(ptvl.gt. 60) call pwhgfill(43+numplots*ncut,yvl,dsig/bsz(43))
      if(ptvl.gt. 80) call pwhgfill(44+numplots*ncut,yvl,dsig/bsz(44))
      if(ptvl.gt.100) call pwhgfill(45+numplots*ncut,yvl,dsig/bsz(45))
      
      call pwhgfill(46+numplots*ncut,ylep,dsig/bsz(46))
      call pwhgfill(47+numplots*ncut,yvl,dsig/bsz(47))
      
      call pwhgfill(48+numplots*ncut,ptlep,dsig/bsz(48))
      call pwhgfill(49+numplots*ncut,ptvl,dsig/bsz(49))
      
      call pwhgfill(56+numplots*ncut,atan2(pvb(2),pvb(1)),dsig/bsz(56))

      call pwhgfill(57+numplots*ncut,0.5d0,dsig/bsz(57))

 987  continue
c     end of loop on ptW cuts
      enddo


      end



      subroutine build_jets(njets,pjet)
c     arrays to reconstruct jets
      implicit none
      include '../include/hepevt.h'
      integer njets
      real * 8 pjet(4,*)
      integer maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real *8 ptrack(4,maxtrack)
      real *8 pp
      integer jetvec(maxtrack)
      integer ihep,j,j1,ntracks,jpart,jjet,mu
      real * 8 found
      real * 8 random
      integer seed
      data seed/1/
      save seed
      real *8 R,sf,f,caf,Etmin
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
         if (isthep(ihep).eq.1) then
            if(ntracks.eq.maxtrack) then
               write(*,*)
     #              'analyze: too many particles, increase maxtrack'
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
c      R = 0.7d0          
c      ptmin_fastkt = 20d0
c      call fastjetktwhich(ptrack,ntracks,ptmin_fastkt,R,
c     #     pjet,njets,jetvec)
c     now we have the jets
c MidPoint CDF
c$$$      R=0.7
c$$$      f=0.75
c$$$      sf=1
c$$$      caf=1
c$$$      call fastjetCDFMidPoint(ptrack,ntracks,R,f,sf,caf,pjet,njets)
c D0RunII Cone
c$$$      R=0.7
c$$$      f=0.5
c$$$      Etmin=15d0
c$$$      call fastjetD0RunIICone(ptrack,ntracks,R,Etmin,f,pjet,njets)

      end

      subroutine getktetaphi(njets,pjet,ktjet,etajet,phijet)
      implicit none
      integer njets
      real * 8 pjet(4,njets),ktjet(njets),etajet(njets),phijet(njets),pp
      integer j
      do j=1,njets
         ktjet(j)=sqrt(pjet(1,j)**2+pjet(2,j)**2)
         pp = sqrt(ktjet(j)**2+pjet(3,j)**2)
         etajet(j)=0.5d0*log((pp+pjet(3,j))/(pp-pjet(3,j)))
         phijet(j)=atan2(pjet(2,j),pjet(1,j))
      enddo
      end
      


      subroutine getrapidity(p,y)
      implicit none
      real * 8 p(4),y
      y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      end

      subroutine getinvmass(p,m)
      implicit none
      real * 8 p(4),m
      m=sqrt(abs(p(4)**2-p(1)**2-p(2)**2-p(3)**2))
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

      function getrapidity0(p)
      implicit none
      real * 8 p(0:3),getrapidity0
      getrapidity0=0.5d0*log((p(0)+p(3))/(p(0)-p(3)))
      end



c     calculate the separation in the lego plot between the two momenta
c     p1 and p2 in azi and pseudorapidity
      function rsepn_p(p1,p2)
      implicit none
      real * 8 pi,pi2
      parameter(pi = 3.141592653589793D0, pi2 = 9.869604401089358D0)
      real * 8 rsepn_p,p1(0:3),p2(0:3)
      real * 8 eta1,phi1,eta2,phi2
      real * 8 delphi
      real * 8 pseudorapidity,azi
      external pseudorapidity,azi

      phi1 = azi(p1)   
      phi2 = azi(p2)
      eta1 = pseudorapidity(p1)
      eta2 = pseudorapidity(p2)

      delphi = abs(phi1-phi2)
      if (delphi.gt.pi) then
         delphi = 2*pi-delphi
      endif
      if (delphi.lt.0 .or. delphi.gt.pi) then
         print*,' problem in rsepn. delphi = ',delphi
      endif
      rsepn_p = sqrt( (eta1-eta2)**2 + delphi**2 )
      end

      function azi(p)
      implicit none
      real * 8 pi,pi2
      parameter(pi = 3.141592653589793D0, pi2 = 9.869604401089358D0)
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

      function pseudorapidity(p)
      implicit none
      real * 8 p(0:3),pseudorapidity
      real * 8 mod, costh
      mod = sqrt(p(1)**2+p(2)**2+p(3)**2)
      costh = p(3)/mod
      pseudorapidity=0.5*log((1+costh)/(1-costh))
      end




      function pass_lept_cuts(p)
      implicit none
      logical pass_lept_cuts
      real * 8 p(1:4)
      real * 8 yl1
      common/clepton_cuts/yl1
      pass_lept_cuts = .true.
      end
