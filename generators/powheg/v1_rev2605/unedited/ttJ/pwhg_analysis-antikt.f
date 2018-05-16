c     The next subroutines, open some histograms and prepare them 
c     to receive data. The relevant routines are 
c     init_hist   :  opens the histograms
c     pwhgfill    :  fills the histograms with data
c     pwhgaccumup :  accumulates results and clear temporary histograms
c     pwhgsetout  :  performs statistical analysis and close histograms
c     pwhgtopout  :  output the results on a .top file
c     You can substitute any of these with your favourite ones

      subroutine init_hist
      implicit none
      include 'pwhg_book.h'
      double precision pi,pi2
      parameter(pi = 3.141592653589793D0, pi2 = 9.869604401089358D0)
      integer icut,idiag,diag
      character * 3 jetcut
      character * 4 attcut
c     binsize
      double precision bsz(nmh)
      common/pwhghistcommon/bsz
c
      integer nptttbcutmax,nptttbcut
      parameter(nptttbcutmax = 10)
      integer ptttbcuts(0:nptttbcutmax-1)
      common/cptvbcut/ptttbcuts,nptttbcut,numplots
      integer ncut,numplots
      
      character * 3 ttbcut
      integer i

      ptttbcuts(0) = 0
      ptttbcuts(1) = 10
      ptttbcuts(2) = 20
      ptttbcuts(3) = 50
c     number of pt cuts implemented one more than the last previous one
      nptttbcut=4

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      numplots = 109  ! <========== DO NOT FORGET TO SET THIS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      write(*,*) '**************************************************'
      write(*,*) '**************************************************'
      write(*,*) '             INITIALIZING HISTOGRAMS          '
      write(*,*) '**************************************************'
      

      call pwhginihist

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C LOOP ON TTBAR PT CUTS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do ncut=0,nptttbcut-1 
         write(unit=ttbcut,fmt="(i3)") ptttbcuts(ncut)

c total cross section sanity check
         bsz(1) = 1d0
         call pwhgbookup(1+numplots*ncut,'total pt(ttb)>'//ttbcut,'LOG'
     $        ,bsz(1),0d0,1d0)
c yt asymmetry
         bsz(2) = 1d0
         call pwhgbookup(2+numplots*ncut,'yt asym. pt(ttb)>'//ttbcut,
     $        'LOG',bsz(2),0d0,1d0)
c yt asymmetry |yt|<1
         bsz(3) = 1d0
         call pwhgbookup(3+numplots*ncut,'yt asym. |Dy|<1 pt(ttb)>'
     $        //ttbcut,'LOG',bsz(3),0d0,1d0)
c total |yt|< 1
         bsz(4) = 1d0
         call pwhgbookup(4+numplots*ncut,'total |Dy|<1 pt(ttb)>'
     $        //ttbcut,'LOG',bsz(4),0d0,1d0)
c yt asymmetry |yt|> 1
         bsz(5) = 1d0
         call pwhgbookup(5+numplots*ncut,'yt asym. |Dy|>1 pt(ttb)>'
     $        //ttbcut,'LOG',bsz(5),0d0,1d0)
c total |yt|> 1
         bsz(6) = 1d0
         call pwhgbookup(6+numplots*ncut,'total |Dy|>1 pt(ttb)>'
     $        //ttbcut,'LOG',bsz(6),0d0,1d0)
c yt asymmetry Mtt<450
         bsz(7) = 1d0
         call pwhgbookup(7+numplots*ncut,'yt asym. Mtt<450 pt(ttb)>'
     $        //ttbcut,'LOG',bsz(7),0d0,1d0)
c total Mtt < 450
         bsz(8) = 1d0
         call pwhgbookup(8+numplots*ncut,'total Mtt<450 pt(ttb)>'
     $        //ttbcut,'LOG',bsz(8),0d0,1d0)
c yt asymmetry Mtt>450
         bsz(9) = 1d0
         call pwhgbookup(9+numplots*ncut,'yt asym. Mtt>450 pt(ttb)>'
     $        //ttbcut,'LOG',bsz(9),0d0,1d0)
c total Mtt > 450
         bsz(10) = 1d0
         call pwhgbookup(10+numplots*ncut,'total Mtt>450 pt(ttb)>'
     $        //ttbcut,'LOG',bsz(10),0d0,1d0)
c top antitop charge asymmetry
         bsz(11) = 1d0
         call pwhgbookup(11+numplots*ncut,'t-tbar charge asym. pt(ttb)>'
     $        //ttbcut,'LOG',bsz(11),0d0,1d0)
c mass of the pair
         bsz(12)=20d0
         call pwhgbookup(12+numplots*ncut,'ttbar mass pt(ttb)>'//ttbcut,
     $        'LOG',bsz(12),0d0,1000d0)
c rapidity of the pair
         bsz(13)=0.2d0
         call pwhgbookup(13+numplots*ncut,'ttbar y pt(ttb)>'//ttbcut,
     $        'LOG',bsz(13),-4d0,4d0)

c pt of the pair
         bsz(14)=4d0
         call pwhgbookup(14+numplots*ncut,'ttbar pt pt(ttb)>'//ttbcut
     $        ,'LOG',bsz(14),0d0,800d0)
c pt of the pair (zoom)
         bsz(15)=0.5d0
         call pwhgbookup(15+numplots*ncut,'ttbar pt zoom pt(ttb)>'/
     $        /ttbcut,'LOG',bsz(15),0d0,100d0)
c pt of top
         bsz(16)=5d0
         call pwhgbookup(16+numplots*ncut,'t pt pt(ttb)>'//ttbcut,'LOG'
     $        ,bsz(16),0d0,500d0)
c pt of antitop
         bsz(17)=5d0
         call pwhgbookup(17+numplots*ncut,'tbar pt pt(ttb)>'//ttbcut
     $        ,'LOG',bsz(17),0d0,500d0)
c top rapidity
         bsz(18)=0.2d0
         call pwhgbookup(18+numplots*ncut,'t y pt(ttb)>'//ttbcut,'LOG'
     $        ,bsz(18),-4d0,4d0)
c antitop rapidity
         bsz(19)=0.2d0
         call pwhgbookup(19+numplots*ncut,'tbar y pt(ttb)>'//ttbcut
     $        ,'LOG',bsz(19),-4d0,4d0)
c top-antitop rapidity asymmetry
         bsz(20)=0.2
         call pwhgbookup(20+numplots*ncut,'t-tbar,y asym.  pt(ttb)>'/
     $        /ttbcut,'LOG',bsz(20),-4d0,4d0)
c top-antitop pseudo-rapidity asymmetry
         bsz(21) = 0.2
         call pwhgbookup(21+numplots*ncut,'t-tbar,eta asym.  pt(ttb)>'/
     $        /ttbcut,'LOG',bsz(21), -4d0,4d0)

c top-antitop charge asymmetry (using y)
         bsz(22)=0.2
         call pwhgbookup(22+numplots*ncut
     $        ,'t-tbar, ch asym.(y)  pt(ttb)>'//ttbcut,'LOG',bsz(22),
     $        -4d0,4d0)

c top-antitop charge asymmetry (using eta)
         bsz(23)=0.2
         call pwhgbookup(23+numplots*ncut
     $        ,'t-tbar, ch asym.(eta)  pt(ttb)>'//ttbcut,'LOG',bsz(23),
     $        -4d0,4d0)

c top-antitop azimuthal distance
         bsz(24) = 0.02d0
         call pwhgbookup(24+numplots*ncut,'t-tbar,Delta_phi pt(ttb)>'/
     $        /ttbcut,'LOG',bsz(24),0d0,1d0)
c top-antitop Delta R distance
         bsz(25) = 0.25d0
         call pwhgbookup(25+numplots*ncut,'t-tbar,Delta_R   pt(ttb)>'/
     $        /ttbcut,'LOG',bsz(25),0 .d0,10.d0)
         

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c  PLOTS INVOLVING JETS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         
c N jets
         bsz(29) = 1d0
         call pwhgbookup(29+numplots*ncut,'Njets   pt(ttb)>'/
     $        /ttbcut,'LOG',bsz(29),0d0,6d0)
         
         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c  PLOTS INVOLVING HARDEST JET
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c (Leading non-b) jet pt
         bsz(30)=5
         call pwhgbookup(30+numplots*ncut,'pt1stjet pt(ttb)>'//ttbcut ,
     $        'LOG',bsz(30),0d0,500d0)

c (Leading non-b) jet pt (zoom)
         bsz(31)=1
         call pwhgbookup(31+numplots*ncut,'pt1stjet zoom pt(ttb)>'/
     $        /ttbcut,'LOG',bsz(31),0d0,100d0)

c (Leading non-b) jet pt_rel
         bsz(32)=1
         call pwhgbookup(32+numplots*ncut,'pt_rel 1stjet pt(ttb)>'/
     $        /ttbcut,'LOG',bsz(32),0d0,60d0)

c (Leading non-b) jet y
         idiag=0
         do icut=20,100,20
            write(jetcut,'(i3)') icut
            bsz(33+idiag)=0.2
            call pwhgbookup(33+idiag+numplots*ncut,'y1stjet, pt(jet)>'/
     $           /jetcut//' pt(ttb)>'//ttbcut,'LOG',bsz(33+idiag),-4d0
     $           ,4d0)
            idiag=idiag+1
         enddo
c (Leading non-b) jet y-yttbar
         idiag=0
         do icut=20,100,20
            write(jetcut,'(i3)') icut
            bsz(38+idiag)=0.2
            call pwhgbookup(38+idiag+numplots*ncut
     $           ,'y1stjet-y(ttbar),pt(jet)>'//jetcut//' pt(ttb)>'/
     $           /ttbcut,'LOG',bsz(38+idiag),-4d0,4d0)
            idiag=idiag+1
         enddo
c (Leading non-b) jet phi-phittbar
         idiag=0
         do icut=20,100,20
            write(jetcut,'(i3)') icut
            bsz(43+idiag)=0.02d0
            call pwhgbookup(43+idiag+numplots*ncut
     $           ,'phi1stjet-phi(ttbar),pt(jet)>'//jetcut//' pt(ttb)>'/
     $           /ttbcut,'LOG',bsz(43+idiag),0d0,1d0)
            idiag=idiag+1
         enddo
c (Leading non-b) jet delta R jet-ttbar
         idiag=0
         do icut=20,100,20
            write(jetcut,'(i3)') icut
            bsz(48+idiag)=0.25d0
            call pwhgbookup(48+idiag+numplots*ncut
     $           ,'deltaR 1stjet-(ttbar),pt(jet)>' //jetcut//' pt(ttb)>'
     $           //ttbcut,'LOG',bsz(48+idiag),0d0,10d0)
            idiag=idiag+1
         enddo      

c (Leading non-b) jet ttbar system pt(log)
         bsz(53)=0.06
         call pwhgbookup(53+numplots*ncut
     $        ,'log10(pt) 1stjet-ttbar pt(ttb)>'//ttbcut ,'LOG',bsz(53)
     $        ,0d0,3d0)
         
c (Leading non-b) jet ttbar system pt
         bsz(54)=3
         call pwhgbookup(54+numplots*ncut ,'pt 1stjet-ttbar pt(ttb)>'/
     $        /ttbcut,'LOG',bsz(54),0d0 ,300d0)
    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c  PLOTS INVOLVING NEXT-TO-HARDEST JET
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c (Next-to-Leading non-b) jet pt
         bsz(60)=10
         call pwhgbookup(60+numplots*ncut,'pt2ndjet pt(ttb)
     $        >'//ttbcut ,'LOG',bsz(60),0d0,600d0)

c (Next-to-Leading non-b) jet pt  (zoom)    
         bsz(61)=1
         call pwhgbookup(61+numplots*ncut,'pt2ndjet (zoom) pt(ttb)>'/
     $        /ttbcut,'LOG',bsz(61),0d0,100d0)

c (Next-to-Leading non-b) jet pt_rel 
         bsz(62)=1
         call pwhgbookup(62+numplots*ncut,'pt_rel 2ndjet pt(ttb)
     $        >'//ttbcut ,'LOG',bsz(62),0d0,60d0)

c (Next-to-Leading non-b) jet y
         idiag=0
         do icut=20,100,20
            write(jetcut,'(i3)') icut
            bsz(63+idiag)=0.2
            call pwhgbookup(63+idiag+numplots*ncut,'y2ndjet,pt(jet)>'/
     $           /jetcut//' pt(ttb)>'//ttbcut, 'LOG',bsz(63+idiag),-4d0
     $           ,4d0)
            idiag=idiag+1
         enddo
c (Next-to-Leading non-b) jet y-yttbar
         idiag=0
         do icut=20,100,20
            write(jetcut,'(i3)') icut
            bsz(68+idiag)=0.2
            call pwhgbookup(68+idiag+numplots*ncut
     $           ,'y2ndjet-y(ttbar),pt(jet)>'//jetcut//' pt(ttb)>'/
     $           /ttbcut,'LOG',bsz(68+idiag),-4d0,4d0)
            idiag=idiag+1
         enddo

c (Next-to-Leading non-b) jet delta R jet-ttbar
         idiag=0
         do icut=20,100,20
            write(jetcut,'(i3)') icut
            bsz(73+idiag)=0.25d0
            call pwhgbookup(73+idiag+numplots*ncut
     $           ,'deltaR 2ndjet-(ttbar),pt(jet)>'//jetcut//' pt(ttb)>'/
     $           /ttbcut,'LOG',bsz(73+idiag),0d0,10d0)
            idiag=idiag+1
         enddo

c  delta R (Leading non-b) jet - (Next-to-Leading non-b) jet
         idiag=0
         do icut=20,100,20
            write(jetcut,'(i3)') icut
            bsz(78+idiag)=0.25d0
            call pwhgbookup(78+idiag+numplots*ncut
     $           ,'deltaR 1stjet-2ndjet,pt(jet)>'//jetcut//' pt(ttb)>'/
     $           /ttbcut,'LOG',bsz(78+idiag),0d0,10d0)
            idiag=idiag+1
         enddo


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c  PLOTS INVOLVING LEPTONS FROM TOP DECAY PRODUCTS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c Lepton+
c Lepton energy in top rest frame 
         bsz(90)=2
         call pwhgbookup(90+numplots*ncut,'eem pt(ttb)>'//ttbcut,'LOG'
     $        ,bsz(90),0d0,100d0)
c Lepton cos theta in top rest frame 
         bsz(91)=0.1
         call pwhgbookup(91+numplots*ncut ,'costh1 costh2 pt(ttb)>'/
     $        /ttbcut,'LOG',bsz(91),-1d0,1d0)
c Lepton azimuth in top rest frame. The origin is the top
c transverse direction
         bsz(92)=0.02d0
         call pwhgbookup(92+numplots*ncut ,'phi dec pt(ttb)
     $        >'//ttbcut ,'LOG',bsz(92),-1d0,1d0)
         
c  Hardest lepton+ pt   
         bsz(93)=3
         call pwhgbookup(93+numplots*ncut,'lep+ pt(ttb)>'//ttbcut,'LOG'
     $        ,bsz(93),0d0,300d0)

c  Hardest lepton- pt   
         bsz(94)=3
         call pwhgbookup(94+numplots*ncut,'lep- pt(ttb)>'//ttbcut,'LOG'
     $        ,bsz(94),0d0,300d0)

c  Hardest lepton+ rapidity   
         bsz(95)=0.2
         call pwhgbookup(95+numplots*ncut ,'lep+ y pt(ttb)>'//ttbcut,
     $        'LOG',bsz(95),-4d0,4d0)

c  Hardest lepton- rapidity   
         bsz(96)=0.2
         call pwhgbookup(96+numplots*ncut ,'lep- y pt(ttb)>'//ttbcut,
     $        'LOG',bsz(96),-4d0,4d0)

c  dilepton invariant mass   
         bsz(97)=5
         call pwhgbookup(97+numplots*ncut ,'dilep inv. mass  pt(ttb)
     $        >'//ttbcut ,'LOG',bsz(97),0d0,500d0)

c dilepton azimuthal distance
         bsz(98)=0.02
         call pwhgbookup(98+numplots*ncut ,'dilep Delta phi  pt(ttb)
     $        >'//ttbcut ,'LOG',bsz(98),0d0,1d0)

c dilepton azimuthal distance with mtt cut
         bsz(99)=0.02
         call pwhgbookup(99+numplots*ncut
     $        ,'dilep Delta phi  mtt < 400 pt(ttb)>'//ttbcut ,'LOG'
     $        ,bsz(99),0d0,1d0)

         

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C PLOTS FOR TRACKING NEGATIVE WEIGHTS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    

      bsz(100) = 2d0      
      call pwhgbookup(100+numplots*ncut,'pt(ttb) pos-|neg| pt(ttb)>'/
     $     /ttbcut,'LOG',bsz(100),0d0,200d0)
      bsz(101) = 2d0      
      call pwhgbookup(101+numplots*ncut,'pt(ttb) pos pt(ttb)
     $     >'//ttbcut ,'LOG',bsz(101),0d0,200d0)
      bsz(102) = 2d0      
      call pwhgbookup(102+numplots*ncut,'pt(ttb) neg pt(ttb)
     $     >'//ttbcut ,'LOG',bsz(102),0d0,200d0)
      bsz(103) = 2d0      
      call pwhgbookup(103+numplots*ncut,'pt J1 pos-|neg| pt(ttb)>'/
     $     /ttbcut,' LOG', bsz(103),0d0,200d0)
      bsz(104) = 2d0      
      call pwhgbookup(104+numplots*ncut,'pt J1 pos pt(ttb)>'//ttbcut
     $     ,'LOG',bsz(104),0d0,200d0)
      bsz(105) = 2d0      
      call pwhgbookup(105+numplots*ncut,'pt J1 neg pt(ttb)>'//ttbcut
     $     ,'LOG',bsz(105),0d0,200d0)
      bsz(106) = 2d0
      call pwhgbookup(106+numplots*ncut,'pt J2 pos-|neg| pt(ttb) >'/
     $     /ttbcut,'LOG', bsz(106),0d0,200d0)
      bsz(107) = 2d0
      call pwhgbookup(107+numplots*ncut,'pt J2 pos pt(ttb)>'//ttbcut
     $     ,'LOG',bsz(107),0d0,200d0)
      bsz(108) = 2d0
      call pwhgbookup(108+numplots*ncut,'pt J2 neg pt(ttb)>'//ttbcut
     $     ,'LOG',bsz(108),0d0,200d0)


      enddo                     ! end of pt (ttb) cuts loop

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C PLOTS IN ARXIV 0810.0452
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      diag=numplots*nptttbcut

      diag=diag+1
      bsz(diag)=25d0
      call pwhgbookup(diag,'DUW jet pt bin4TeV','LOG',bsz(diag)
     $     ,0d0,325d0)
      
      diag=diag+1
      bsz(diag)=50d0
      call pwhgbookup(diag,'DUW jet pt bin4LHC','LOG',bsz(diag)
     $     ,0d0,700d0)

      diag=diag+1
      bsz(diag)=25d0
      call pwhgbookup(diag,'DUW ttbar pt bin4TeV','LOG'
     $     ,bsz(diag),0d0,325d0)

      diag=diag+1
      bsz(diag)=50d0
      call pwhgbookup(diag,'DUW ttbar pt bin4LHC','LOG'
     $     ,bsz(diag),0d0,700d0)

      diag=diag+1
      bsz(diag)=25d0
      call pwhgbookup(diag,'DUW t pt bin4TeV','LOG',bsz(diag),0d0
     $     ,400d0)

      diag=diag+1
      bsz(diag)=50d0
      call pwhgbookup(diag,'DUW t pt bin4LHC','LOG',bsz(diag),0d0
     $     ,700d0)

      diag=diag+1
      bsz(diag)=0.4d0
      call pwhgbookup(diag,'DUW t eta bin4TeV','LOG',bsz(diag),
     $     -4d0,4d0)

      diag=diag+1
      bsz(diag)=0.6d0
      call pwhgbookup(diag,'DUW t eta bin4LHC','LOG',bsz(diag),
     $     -6d0,6d0)

      diag=diag+1
      bsz(diag)=0.4d0
      call pwhgbookup(diag,'DUW t y bin4TeV','LOG',bsz(diag),-2d0
     $     ,2d0)

      diag=diag+1
      bsz(diag)=0.6d0
      call pwhgbookup(diag,'DUW t y bin4LHC','LOG',bsz(diag),
     $     -4.2d0,4.2d0)

      diag=diag+1
      bsz(diag)=0.4d0
      call pwhgbookup(diag,'DUW jet y bin4TeV','LOG',bsz(diag),
     $     -4d0,4d0)

      diag=diag+1
      bsz(diag)=0.6d0
      call pwhgbookup(diag,'DUW jet y bin4LHC','LOG',bsz(diag),-5.4d0,5
     $     .4d0)
     


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C PLOTS FOR ASYMMETRY
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
      diag=diag+1
      
c (Leading non-b) jet pt with inv mass cuts
      idiag=0
      do icut=0,1000,200
         write(attcut,'(i4)') icut
         bsz(diag+idiag)=50
         call pwhgbookup(diag+idiag,'pt 1st jet , mtt>'//attcut,'LOG'
     $        ,bsz(diag+idiag),0d0,1000d0)
         idiag=idiag+1
      enddo

c inv mass with jet cuts
      diag=diag+idiag
      idiag=0
      do icut=0,1000,200
         write(attcut,'(i4)') icut
         bsz(diag+idiag)=50
         call pwhgbookup(diag+idiag,'mtt,  pt 1stjet>'//attcut,'LOG'
     $        ,bsz(diag+idiag),0d0,1000d0)
         idiag=idiag+1
      enddo
      diag=diag+idiag
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C PLOTS FOR MASS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
      diag=diag+1
      
      bsz(diag)=0.05d0
      call pwhgbookup(diag,'m*','LIN',bsz(diag),0d0
     $     ,1d0)

      end


     
      subroutine analysis(dsig)
      implicit none
      double precision dsig
      include 'hepevt.h'
      double precision pi
      parameter(pi = 3.141592653589793D0)
c arrays to reconstruct jets
      integer maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=50)
      double precision pjet(4,maxjet),pj(4,maxjet)
      integer njets,j
      double precision ktjet(maxjet),yjet(maxjet),phijet(maxjet)
     $     ,ptreljet(maxjet)
      logical ini
      data ini/.true./
      save ini
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
c     binsize
      include 'pwhg_book.h'
      double precision bsz(nmh)
      common/pwhghistcommon/bsz
      integer diag,idiag
c
      integer nptttbcutmax,nptttbcut
      parameter(nptttbcutmax = 10)
      integer ptttbcuts(0:nptttbcutmax-1)
      common/cptvbcut/ptttbcuts,nptttbcut,numplots
      integer ncut,numplots
c
      double precision deltaphi,rsep_azi_y,get_ptrel_jet
      external deltaphi,rsep_azi_y,get_ptrel_jet
      character * 20 jetalgo
      integer icut
      integer it,itbar,iem,iep,iit,iitbar,iheptop,iheptbar,iiem,iiep
      double precision phit,phitbar,deltaRttbar,deltaphittbar
     $     ,deltaRj1ttbar,deltaphij1ttbar,deltaRj1j2,deltaRj2ttbar
      double precision ppairttbar(4),ptt,yt,etat,pttbar,ytbar,etatbar,
     $     ptttbar,mttbar,yttbar,phittbar,etattbar,eep,cthep,phiep,eem
     $     ,cthem,phiem,dphi,ptep,ptem,mepem,yep,yem,etaep,etaem,
     $     ppairepem(4),deltaphiepem,diffabseta,psysttbarj(4)
     $     ,pt2systtbarj
      integer mu,ihep
      double precision maxptep,maxptem
      double precision m0,mstar
      integer ibtag(maxjet),itags
      logical sonofidhep,is_lepton,is_antilepton
      external sonofidhep,is_lepton,is_antilepton
      if (ini) then
         write(*,*) '*****************************'
         if(WHCPRG.eq.'NLO   ') then
            write(*,*) '       NLO ANALYSIS'
         elseif(WHCPRG.eq.'LHE   ') then
            write (*,*) '           LHE ANALYSIS            '
         elseif(WHCPRG.eq.'HERWIG') then
            write (*,*) '           HERWIG ANALYSIS            '
         elseif(WHCPRG.eq.'PYTHIA') then
            write (*,*) '           PYTHIA ANALYSIS            '
         endif
         write(*,*) '*****************************'
         ini=.false.
      endif
c RESET ALL COUNTERS
      it=0
      iit=0
      iheptop=0
      itbar=0
      iitbar=0
      iheptbar=0
      iem=0
      iiem=0
      maxptem=0d0
      iep=0
      iiep=0
      maxptep=0d0
c Parton level analysis
      if((whcprg.eq.'NLO   ').or.(whcprg.eq.'LHE   ')) then
         do ihep=3,nhep
            if(idhep(ihep).eq.6) then
               it=ihep
               iit=iit+1
            elseif(idhep(ihep).eq.-6) then
               itbar=ihep
               iitbar=iitbar+1
            elseif(is_lepton(idhep(ihep))) then
               iem=ihep
               iiem=iiem+1
            elseif(is_antilepton(idhep(ihep)))  then
               iep=ihep
               iiep=iiep+1
            endif
         enddo
c check that exactly 1 top and 1 anti-top are selected
         if (iit.ne.1) then
            write(*,*) "Error in pwhg_analisys: ",iit," tops"
            call printleshouches
            call exit(1)
         endif
         if (iitbar.ne.1) then
            write(*,*) "Error in pwhg_analisys: ",iitbar," anti-tops"
            call printleshouches
            call exit(1)
         endif
c Hadron level analysis
      else
         do ihep=1,nhep
            if(isthep(ihep).eq.1) then
               if(is_antilepton(idhep(ihep))) then
                  if (sonofidhep(6,ihep,j)) then
                     it=j
c     check that the positrons that come from a top come from the same one 
                     if (iit.eq.0) then
                        iheptop=it
                        iit=iit+1
                        iiep=iiep+1
                     endif
                     if ((iit.gt.0).and.(it.ne.iheptop)) then
                        iit=iit+1
                        iiep=iiep+1
                     endif
c     select the hardest positron coming from the same top
                     if(phep(1,ihep)**2+phep(2,ihep)**2.gt.maxptep) then
                        iep=ihep
                        maxptep=phep(1,iep)**2+phep(2,iep)**2
                     endif
                  endif
               elseif(is_lepton(idhep(ihep))) then
                  if(sonofidhep(-6,ihep,j)) then
                     itbar=j
c     check that the electrons that come from a tbar come from the same one 
                     if (iitbar.eq.0) then
                        iheptbar=itbar
                        iitbar=iitbar+1
                        iiem=iiem+1
                     endif
                     if ((iitbar.gt.0).and.(itbar.ne.iheptbar)) then
                        iitbar =iitbar+1
                        iiem=iiem+1
                     endif
c     select the hardest electron coming from the same tbar
                     if(phep(1,ihep)**2+phep(2,ihep)**2.gt.maxptem) then
                        iem=ihep
                        maxptem=phep(1,iem)**2+phep(2,iem)**2
                     endif
                  endif
               endif
            endif
         enddo
c check that at least one top and one anti-top are selected
c during showering there may be more than one entry with PDG=6,
c with different IST. The sonof procedure select the first back
c in shower history.
         if (iit.lt.1) then
            write(*,*) "Error in pwhg_analisys: ",iit," tops"
            call printleshouches
            call exit(1)
         endif
         if (iitbar.lt.1) then
            write(*,*) "Error in pwhg_analisys: ",iitbar," anti-tops"
            call printleshouches
            call exit(1)
         endif
c check that at least one antilepton and one lepton are selected
         if (iiep.lt.1) then
            write(*,*) "Error in pwhg_analisys: ",iep," leptons"
            call printleshouches
            call exit(1)
         endif
         if (iiem.lt.1) then
            write(*,*) "Error in pwhg_analisys: ",iem," anti-leptons"
            call printleshouches
            call exit(1)
         endif
      endif
c     evaluate quantities related to top-antitop 
      call ptyetaphi(phep(1,it),ptt,yt,etat,phit)
      call ptyetaphi(phep(1,itbar),pttbar,ytbar,etatbar,phitbar)
      diffabseta=abs(etat)-abs(etatbar)
      do mu=1,4
         ppairttbar(mu)=phep(mu,it)+phep(mu,itbar)
      enddo
      call getinvmass(ppairttbar,mttbar)
      call ptyetaphi(ppairttbar,ptttbar,yttbar,etattbar,phittbar)

      deltaphittbar=deltaphi(phit,phitbar)
      deltaRttbar=rsep_azi_y(phep(1,it),phep(1,itbar))

      
c     build jets with the inclusive kt algo (TOPS and b-jet NOT INCLUDED)
      jetalgo="antikt"
c     njets is the maximum number of jets to be searched for (buildjets
c     then sets it to the actual value found)
      if((whcprg.eq.'NLO   ').or.(whcprg.eq.'LHE   ')) then
         njets=2
      else
         njets=10
      endif

      call build_non_b_jets(njets,pjet,jetalgo,ptreljet)
      call getktyphi(njets,pjet,ktjet,yjet,phijet)


c     if there are no jets skip the analysis
      if (njets.lt.1) then
         dsig=0
         return
      endif

c     evaluate quantities related to leptons (if present) 
      if(iiep.gt.0) then
c     Compute lepton variables in top rest frame 
         call decvariables(phep(1,iep),phep(1,it),eep,cthep,phiep)
c     Compute lepton variables in laboratory frame  
         call ptyetaphi(phep(1,iep),ptep,yep,etaep,phiep)
      endif
      if(iiem.gt.0) then
c     Compute lepton variables in top rest frame 
         call decvariables(phep(1,iem),phep(1,itbar),eem,cthem,phiem)
c     Compute lepton variables in laboratory frame        
         call ptyetaphi(phep(1,iem),ptem,yem,etaem,phiem)
      endif
      if((iiep.gt.0).and.(iiem.gt.0)) then   
         dphi=(phiep-phiem)-nint((phiep-phiem)/(2*pi))*2*pi
         do mu=1,4
            ppairepem(mu)=phep(mu,iep)+phep(mu,iem)
         enddo
         call getinvmass(ppairepem,mepem)
         deltaphiepem=deltaphi(phiep,phiem)
      endif


c     fill histograms
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C LOOP ON TTBAR PT CUTS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCC
      do ncut=0,nptttbcut-1
         
         if (ptttbar.lt.ptttbcuts(ncut)) goto 999
c total sigma
         call pwhgfill(1+numplots*ncut,0.5d0,dsig/bsz(1))
c yt asymmetry
         if(yt.gt.0d0) call pwhgfill(2+numplots*ncut,0.5d0,dsig/bsz(2))
         if(yt.lt.0d0) call pwhgfill(2+numplots*ncut,0.5d0,-dsig/bsz(2))
c yt asymmetry |Dy|<1
         if((yt.gt.0d0).and.(abs(yt).lt.1d0)) call pwhgfill(3+numplots
     $        *ncut,0.5d0,dsig/bsz(3))
         if((yt.lt.0d0).and.(abs(yt).lt.1d0)) call pwhgfill(3+numplots
     $        *ncut,0.5d0,-dsig/bsz(3))
c total |Dy|<1
         if(abs(yt).lt.1d0) call pwhgfill(4+numplots*ncut,0.5d0,dsig
     $        /bsz(4))
    
c yt asymmetry |Dy|>1
         if((yt.gt.0d0).and.(abs(yt).ge.1d0)) call pwhgfill(5+numplots
     $        *ncut,0.5d0,dsig/bsz(5))
         if((yt.lt.0d0).and.(abs(yt).ge.1d0)) call pwhgfill(5+numplots
     $        *ncut,0.5d0,-dsig/bsz(5))

c total |Dy|>1
         if(abs(yt).ge.1d0) call pwhgfill(6+numplots*ncut,0.5d0,dsig
     $        /bsz(6))

c yt asymmetry Mtt<450
         if((yt.gt.0d0).and.(mttbar.lt.450d0)) call pwhgfill(7+numplots
     $        *ncut,0.5d0,dsig/bsz(7))
         if((yt.lt.0d0).and.(mttbar.lt.450d0)) call pwhgfill(7+numplots
     $        *ncut,0.5d0,-dsig/bsz(7))

c total Mtt<450
         if(mttbar.lt.450d0) call pwhgfill(8+numplots*ncut,0.5d0,dsig
     $        /bsz(8))


c yt asymmetry Mtt>450
         if((yt.gt.0d0).and.(mttbar.ge.450d0)) call pwhgfill(9+numplots
     $        *ncut,0.5d0,dsig/bsz(9))
         if((yt.lt.0d0).and.(mttbar.ge.450d0)) call pwhgfill(9+numplots
     $        *ncut,0.5d0,-dsig/bsz(9))

c total Mtt<450
         if(mttbar.ge.450d0) call pwhgfill(10+numplots*ncut,0.5d0,dsig
     $        /bsz(10))

c top antitop charge asymmetry
         if(diffabseta.gt.0d0) call pwhgfill(11+numplots*ncut,0.5d0,dsig
     $        /bsz(11))
         if(diffabseta.lt.0d0) call pwhgfill(11+numplots*ncut,0.5d0,
     $        -dsig/bsz(11))

c mass of the pair
         call pwhgfill(12+numplots*ncut,mttbar,dsig/bsz(12))
c rapidity of the pair
         call pwhgfill(13+numplots*ncut,yttbar,dsig/bsz(13))
c pt of the pair
         call pwhgfill(14+numplots*ncut,ptttbar,dsig/bsz(14))
c pt of the pair zoom
         call pwhgfill(15+numplots*ncut,ptttbar,dsig/bsz(15))
c pt of top
         call pwhgfill(16+numplots*ncut,ptt,dsig/bsz(16))
c pt of antitop
         call pwhgfill(17+numplots*ncut,pttbar,dsig/bsz(17))
c top rapidity
         call pwhgfill(18+numplots*ncut,yt,dsig/bsz(18))
c antitop rapidity
         call pwhgfill(19+numplots*ncut,ytbar,dsig/bsz(19))
c top - antitop rapidity asymmetry
         call pwhgfill(20+numplots*ncut,yt,dsig/bsz(20))
         call pwhgfill(20+numplots*ncut,ytbar,-dsig/bsz(20))
c top - antitop pseudo-rapidity asymmetry
         call pwhgfill(21+numplots*ncut,etat,dsig/bsz(21))
         call pwhgfill(21+numplots*ncut,etatbar,-dsig/bsz(21))
c top - antitop charge asymmetry (using y)
         call pwhgfill(22+numplots*ncut,abs(yt)-abs(ytbar),dsig/bsz(22))
c top - antitop charge asymmetry (using eta)
         call pwhgfill(23+numplots*ncut,diffabseta,dsig/bsz(23))

c top - antitop azimuthal distance
         call pwhgfill(24+numplots*ncut,deltaphittbar/pi,dsig/bsz(24))
c top - antitop Delta R   
         call pwhgfill(25+numplots*ncut,deltaRttbar,dsig/bsz(25))

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          PLOTS INVOLVING JETS      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c N jets
         call pwhgfill(29+numplots*ncut,njets+0.5d0,dsig/bsz(29))

         deltaphij1ttbar=deltaphi(phijet(1),phittbar)
   
C--- consistency check, only works if tops decay semileptonically
c         if((whcprg.eq.'NLO   ').or.(whcprg.eq.'LHE   ')) then
c            if (abs(deltaphij1ttbar)/pi.lt.0.5d0) then
c               print *,"*************************************"
c               print *,"Event : ",nevhep
c               print *,"deltaphij1ttbar < pi/2 with only 1 jets"
c               print *," numjets ",njets
c               print *,deltaphij1ttbar,phijet(1),phittbar
c               call printleshouches
c               print *,"*************************************"
c            endif
c         endif

         deltaRj1ttbar=rsep_azi_y(pjet(1,1),ppairttbar)

         do mu=1,4
            psysttbarj(mu)=ppairttbar(mu)+pjet(mu,1)
         enddo
         pt2systtbarj=psysttbarj(1)**2+psysttbarj(2)**2
         if(pt2systtbarj.le.0d0)  pt2systtbarj=1d-10
         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c  PLOTS INVOLVING HARDEST JET
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c (Leading non-b) jet pt
         call pwhgfill(30+numplots*ncut,ktjet(1),dsig/bsz(30))
         call pwhgfill(31+numplots*ncut,ktjet(1),dsig/bsz(31))
         
c (Leading non-b) jet pt_rel
         call pwhgfill(32+numplots*ncut,ptreljet(1),dsig/bsz(32))
c (Leading non-b) jet y
         idiag=0
         do icut=20,100,20
            if(ktjet(1).gt.icut) then
               call pwhgfill(33+idiag+numplots*ncut,yjet(1),dsig/bsz(33
     $              +idiag))
            endif
            idiag=idiag+1
         enddo
c (Leading non-b) jet y-yttbar
         idiag=0
         do icut=20,100,20
            if(ktjet(1).gt.icut) then
               call pwhgfill(38+idiag+numplots*ncut,yjet(1)-yttbar,dsig
     $              /bsz(38+idiag))
            endif
            idiag=idiag+1
         enddo
c (Leading non-b) jet phi-phittbar
         idiag=0
         do icut=20,100,20
            if(ktjet(1).gt.icut) then
               call pwhgfill(43+idiag+numplots*ncut,deltaphij1ttbar/pi
     $              ,dsig/bsz(43+idiag))
            endif
            idiag=idiag+1
         enddo
c (Leading non-b) jet delta R jet-ttbar
         idiag=0
         do icut=20,100,20
            if(ktjet(1).gt.icut) then
               call pwhgfill(48+idiag+numplots*ncut,deltaRj1ttbar,dsig
     $              /bsz(48+idiag))
            endif
            idiag=idiag+1
         enddo

c (Leading non-b) jet +ttbar system pt
         call pwhgfill(53+numplots*ncut,0.5d0*log10(pt2systtbarj),dsig
     $        /bsz(53))
         call pwhgfill(54+numplots*ncut,sqrt(pt2systtbarj),dsig/bsz(54))

         if (njets.lt.2) goto 777
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c  PLOTS INVOLVING NEXT-TO-HARDEST JET
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         deltaRj2ttbar=rsep_azi_y(pjet(1,2),ppairttbar)
         deltaRj1j2=rsep_azi_y(pjet(1,1),pjet(1,2))
        
c (Next-to-Leading non-b) jet pt
         call pwhgfill(60+numplots*ncut,ktjet(2),dsig/bsz(60))
         call pwhgfill(61+numplots*ncut,ktjet(2),dsig/bsz(61))

c (Next-to-Leading non-b) jet pt_rel
         call pwhgfill(62+numplots*ncut,ptreljet(2),dsig/bsz(62))

c (Next-to-Leading non-b) jet y
         idiag=0
         do icut=20,100,20
            if(ktjet(2).gt.icut) then
               call pwhgfill(63+idiag+numplots*ncut,yjet(2),dsig/bsz(63
     $              +idiag))
            endif
            idiag=idiag+1
         enddo
c (Next-to-Leading non-b) jet y-yttbar
         idiag=0
         do icut=20,100,20
            if(ktjet(2).gt.icut) then
               call pwhgfill(68+idiag+numplots*ncut,yjet(2)-yttbar,dsig
     $              /bsz(68+idiag))
            endif
            idiag=idiag+1
         enddo
c (Next-to-Leading non-b) jet delta R jet-ttbar
         idiag=0
         do icut=20,100,20
            if(ktjet(2).gt.icut) then
               call pwhgfill(73+idiag+numplots*ncut,deltaRj2ttbar,dsig
     $              /bsz(73+idiag))
            endif
            idiag=idiag+1
         enddo

c  delta R (Leading non-b) jet - (Next-to-Leading non-b) jet
         idiag=0
         do icut=20,100,20 
            if(ktjet(2).gt.icut) then
               call pwhgfill(78+idiag+numplots*ncut,deltaRj1j2,dsig
     $              /bsz(78+idiag))
            endif
            idiag=idiag+1
         enddo

         
 777     continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C PLOTS FOR TRACKING NEGATIVE WEIGHTS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  ttbar pt         
         call pwhgfill(100+numplots*ncut,ptttbar,dsig/bsz(100))
         if (dsig.gt.0d0) then
            call pwhgfill(101+numplots*ncut,ptttbar,dsig/bsz(101))
         else
            call pwhgfill(102+numplots*ncut,ptttbar,abs(dsig)/bsz(102))
         endif
c  hardest jet pt         
         call pwhgfill(103+numplots*ncut,ktjet(1),dsig/bsz(103))
         if (dsig.gt.0d0) then
            call pwhgfill(104+numplots*ncut,ktjet(1),dsig/bsz(104))
         else
            call pwhgfill(105+numplots*ncut,ktjet(1),abs(dsig)/bsz(105))
         endif

         if (njets.lt.2) goto 888
c next-to-hardest jet pt
         call pwhgfill(106+numplots*ncut,ktjet(2),dsig/bsz(106))
         if (dsig.gt.0d0) then
            call pwhgfill(107+numplots*ncut,ktjet(2),dsig/bsz(107))
         else
            call pwhgfill(108+numplots*ncut,ktjet(2),abs(dsig)/bsz(108))
         endif
         
 888     continue
c leptons may not be there, in case we are called by the NLO
c section
         if(iiep.eq.0.and.iiem.eq.0) goto 999        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c  PLOTS INVOLVING LEPTONS FROM TOP DECAY PRODUCTS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c Lepton energy in top rest frame 
         if (iiep.gt.0) call pwhgfill(90+numplots*ncut,eep,dsig/bsz(90))
c Lepton cos theta in top rest frame 
         if (iiep.gt.0) call pwhgfill(91+numplots*ncut,cthep*cthem,dsig
     $        /bsz(91))
c Lepton azimuth in top rest frame. The origin is the top
c transverse direction
         if (iiep*iiem.gt.0) call pwhgfill(92+numplots*ncut,dphi/pi,dsig
     $        /bsz(92))

c Hardest lepton+ pt in lab. frame
         if (iiep.gt.0) call pwhgfill(93+numplots*ncut,ptep,dsig/bsz(93))
c Hardest lepton- pt in lab. frame
         if (iiem.gt.0) call pwhgfill(94+numplots*ncut,ptem,dsig/bsz(94))
c Hardest lepton+ y in lab. frame
         if (iiep.gt.0) call pwhgfill(95+numplots*ncut,yep,dsig/bsz(95))
c Hardest lepton- y in lab. frame
         if (iiem.gt.0) call pwhgfill(96+numplots*ncut,yem,dsig/bsz(96))
c Dilepton inv. mass
         if (iiep*iiem.gt.0)call pwhgfill(97+numplots*ncut,mepem,dsig
     $        /bsz(97))
c Dilepton azimuthal sep.
         if (iiep*iiem.gt.0)call pwhgfill(98+numplots*ncut,deltaphiepem
     $        /pi,dsig/bsz(98))
c Dilepton azimuthal sep. with mttcut
         if((iiep*iiem.gt.0).and.(mttbar.le.400)) then
            call pwhgfill(99+numplots*ncut,deltaphiepem/pi,dsig/bsz(99))
         endif
 999     continue
      enddo ! end of ptttbcut do loop

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PLOTS IN 0810.0452
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      diag=numplots*nptttbcut

      diag=diag+1
      call pwhgfill(diag,ktjet(1),dsig/bsz(diag))
      diag=diag+1
      call pwhgfill(diag,ktjet(1),dsig/bsz(diag))
      diag=diag+1
      call pwhgfill(diag,ptttbar,dsig/bsz(diag))
      diag=diag+1
      call pwhgfill(diag,ptttbar,dsig/bsz(diag))
      diag=diag+1
      call pwhgfill(diag,ptt,dsig/bsz(diag))
      diag=diag+1
      call pwhgfill(diag,ptt,dsig/bsz(diag))
      diag=diag+1
      call pwhgfill(diag,etat,dsig/bsz(diag))
      diag=diag+1
      call pwhgfill(diag,etat,dsig/bsz(diag))
      diag=diag+1
      call pwhgfill(diag,yt,dsig/bsz(diag))
      diag=diag+1
      call pwhgfill(diag,yt,dsig/bsz(diag))
      diag=diag+1
      call pwhgfill(diag,yjet(1),dsig/bsz(diag))
      diag=diag+1
      call pwhgfill(diag,yjet(1),dsig/bsz(diag))

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PLOTS FOR ASYMMETRY
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      diag=diag+1
      
      idiag=0
      do icut=0,1000,200
         if(mttbar.gt.icut) then
            call pwhgfill(diag+idiag,ktjet(1),dsig/bsz(diag+idiag))
         endif
         idiag=idiag+1
      enddo

      diag=diag+idiag

      idiag=0
      do icut=0,1000,200
         if(ktjet(1).gt.icut) then
            call pwhgfill(diag+idiag,mttbar,dsig/bsz(diag+idiag))
         endif
         idiag=idiag+1
      enddo
      
      diag=diag+idiag
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C PLOTS FOR MASS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      diag=diag+1
      m0= 170d0
      mstar= 2d0* m0/sqrt(ppairttbar(1)*pjet(1,1)+ppairttbar(2)*pjet(2
     $     ,1)+ppairttbar(3)*pjet(3,1)+ppairttbar(4)*pjet(4,1))
      call pwhgfill(diag,mstar,dsig/bsz(diag))

      end

CCCCCCCCCC  Ancillary routines CCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine build_non_b_jets(njets,pjet,process,ptreljet)
c     arrays to reconstruct jets
      implicit none
      include 'hepevt.h'
      integer njets,requestedjets,nonbjets
      double precision pjet(4,*)
      character * 20 process
      integer maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=50)
      double precision ptrack(4,maxtrack),pj(4,maxjet),ptrj(maxjet)
      double precision R_jet,ptmin_jet
      common/cjetdefs/R_jet,ptmin_jet
      integer jetvec(maxtrack),itrackhep(maxtrack)
      integer ihep,j1,ntracks,jpart,jjet,mu,k
      integer ibtag(njets)
      double precision ptreljet(njets)
      double precision tmp,get_ptrel_jet
      external get_ptrel_jet,bson
      logical ini,bson,condition
      data ini/.true./
      save ini
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
c     set up maximum number of requested jets
      requestedjets=njets
c
c     get valid tracks
c     set up arrays for jet finding
      do jpart=1,maxtrack
         do mu=1,4
            ptrack(mu,jpart)=0d0
         enddo
         jetvec(jpart)=0
      enddo      
      do jjet=1,requestedjets
         do mu=1,4
            pjet(mu,jjet)=0d0
            pj(mu,jjet)=0d0
         enddo
         ptreljet(jjet)=0d0
         ptrj(jjet)=0d0
      enddo
      j1=0
      ntracks=0
      njets=0
c     loop over final state particles to find jets 
      do ihep=1,nhep
         if((isthep(ihep).eq.1).and.
c     exclude leptons, gauge and higgs bosons and TOP quarks 
     1        (((abs(idhep(ihep)).lt.6).or.(abs(idhep(ihep)).ge.40))
c     but include gluons
     2        .or.(abs(idhep(ihep)).eq.21))) then
         
              condition=.true.
              if(WHCPRG.eq.'LHE   ') then
                 if((abs(idhep(ihep)).eq.5)
     4              .and.(abs(idhep(jmohep(1,ihep))).eq.6)) then
c     exclude final state b quarks whose mother is a top 
c     (needed in LHEF analysis of decayed event)          
                    condition=.false.
                 endif
              endif
              if (condition) then
                 if (ntracks.eq.maxtrack) then
                   write(*,*) 'Too many particles. Increase maxtrack.'//
     $                   ' PROGRAM ABORTS'
                    call exit(1)
                 endif
c     copy momenta to construct jets 
                 ntracks=ntracks+1
                 do mu=1,4
                    ptrack(mu,ntracks)=phep(mu,ihep)
                 enddo
                 itrackhep(ntracks)=ihep
              endif
         endif
      enddo
      if (ntracks.eq.0) then
         njets=0
         return
      endif
      if (process.eq."antikt") then
         if (ini) then
      write(*,*) '**************************************************'
      write(*,*) '**************************************************'
      write(*,*) '                JET PARAMETERS                    '
      write(*,*) '**************************************************'
      write(*,*) '**************************************************'
      write(*,*) '   inclusive anti-kt (FASTJET implementation): '
      write(*,*) '      jet radius ',  R_jet
      write(*,*) '      jet ptmin ',   ptmin_jet
      write(*,*) '**************************************************'
      write(*,*) '**************************************************'
      if((R_jet.le.0d0).or.(ptmin_jet.le.0d0)) then
         write(*,*) 
         write(*,*) '**************************************************'
         write(*,*) '**************************************************'
         write(*,*) '   ERROR: JET ALGORITHM NOT CORRECTLY DEFINED     '
         write(*,*) '    BOTH R_JET AND PTMIN_JET MUST BE POSITIVE     '
         write(*,*) '**************************************************'
         write(*,*) '**************************************************'
         call exit(-1)
      endif
      ini=.false.
      endif
c     FastJet wrapper is contained in the the file fastjet_wrap.cpp in the same directory of libvirtual
         call fastjet_kt(ptrack,ntracks,ptmin_jet,R_jet,
     $        -1d0,3,0,0,njets,pj,jetvec)

      else
         write(*,*) 'JET ANALYSIS TO USE UNKNOWN:',process
         call exit(1)
      endif

c  Adjust the jet number (if larger than requested)
      njets=min(njets,requestedjets)
c  otherwise skip analysis if no jets found

      if (njets.eq.0)  return
      
c    Evaluate pt_rel          
      do jjet=1,njets
         ptrj(jjet)=get_ptrel_jet(jjet,ntracks,ptrack,pj,jetvec)
      enddo

      do jjet=1,njets
         ibtag(jjet)=0
      enddo
      if(whcprg.ne.'NLO   ') then
c Find b decay products among tracks and fill non-b jet arrays
         do k=1,ntracks
            if(jetvec(k).gt.0) then
               if(bson(itrackhep(k))) then
                  ibtag(jetvec(k))=ibtag(jetvec(k))+1
               endif
            endif
         enddo
      endif
      nonbjets=0
      do jjet=1,njets
         if(ibtag(jjet).eq.0) then
            nonbjets=nonbjets+1
            do mu=1,4
                pjet(mu,nonbjets)=pj(mu,jjet)
            enddo
            ptreljet(nonbjets)=ptrj(jjet)
         endif
      enddo
      
      
c     set output value of njets
      njets=nonbjets

      

c$$$c check consistency ( CAVEAT it only works for of E-scheme recombination )
c$$$      do k=1,ntracks
c$$$         if(jetvec(k).gt.0) then
c$$$            do mu=1,4
c$$$               pj(mu,jetvec(k))=pj(mu,jetvec(k))+ptrack(mu,k)
c$$$            enddo
c$$$         endif
c$$$      enddo
c$$$      tmp=0
c$$$      do jjet=1,njets
c$$$         do mu=1,4
c$$$            tmp=tmp+abs(pj(mu,jjet)-pjet(mu,jjet))
c$$$         enddo
c$$$      enddo
c$$$      if(tmp.gt.1d-4) then
c$$$         write(*,*) ' bug!',tmp,njets
c$$$         do jjet=1,njets
c$$$            write(*,*)
c$$$            write(*,'(5(d10.4,1x))') (pjet(k,jjet),k=1,4)
c$$$            write(*,*) "--------------------------------"
c$$$            write(*,'(5(d10.4,1x))') (pj(k,jjet),k=1,4)
c$$$            write(*,*)
c$$$         enddo
c$$$         stop
c$$$      endif

      end
      
 
      function idigit(k,l)
      implicit none
      integer idigit,k,l
      idigit=abs(mod(l,10**k)/10**(k-1))
      end

      subroutine getrapidity(p,y)
      implicit none
      double precision p(4),y
      y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      end

      function getrapidity0(p)
      implicit none
      double precision p(0:3),getrapidity0
      getrapidity0=0.5d0*log((p(0)+p(3))/(p(0)-p(3)))
      end


      
      subroutine getinvmass(p,m)
      implicit none
      double precision p(4),m
      double precision m2
      m2 = p(4)**2-p(1)**2-p(2)**2-p(3)**2
      if (m2.ge.0d0) then
         m = sqrt(abs(m2))
      else
         m = -sqrt(abs(m2))
      endif
      end


      function get_ptrel_jet(ijet,ntracks,ptrack,pjet,jetvec)
      implicit none
      include 'hepevt.h'
      double precision pjet(4,*),pjetin(0:3),pjetout(0:3),beta,vec(3)
      integer maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=50)
      double precision ptrack(4,maxtrack),ptrackin(0:3),ptrackout(0:3)
      double precision ptrel,get_ptrel,get_ptrel_jet
      external get_ptrel
      integer jetvec(maxtrack),ntracks
      integer ijet,mu,i
      do mu=1,3
         pjetin(mu) = pjet(mu,ijet)
      enddo
      pjetin(0) = pjet(4,ijet)         
      vec(1)=0d0
      vec(2)=0d0
      vec(3)=1d0
      beta = -pjet(3,ijet)/pjet(4,ijet)
      call mboost(1,vec,beta,pjetin,pjetout)         
c     write(*,*) pjetout
      ptrel = 0
      do i=1,ntracks
         if (jetvec(i).eq.ijet) then
            do mu=1,3
               ptrackin(mu) = ptrack(mu,i)
            enddo
            ptrackin(0) = ptrack(4,i)
            call mboost(1,vec,beta,ptrackin,ptrackout) 
            ptrel = ptrel + get_ptrel(ptrackout,pjetout)
         endif
      enddo
      get_ptrel_jet=ptrel
      end


      function get_ptrel(pin,pjet)
      implicit none
      double precision get_ptrel,pin(0:3),pjet(0:3)
      double precision pin2,pjet2,cth2,scalprod
      pin2  = pin(1)**2 + pin(2)**2 + pin(3)**2
      pjet2 = pjet(1)**2 + pjet(2)**2 + pjet(3)**2
      scalprod = pin(1)*pjet(1) + pin(2)*pjet(2) + pin(3)*pjet(3)
      cth2 = scalprod**2/pin2/pjet2
      get_ptrel = sqrt(pin2*abs(1d0 - cth2))
      end


      subroutine ptyetaphi(p,pt,y,eta,phi)
      implicit none
      double precision p(4),pt,y,eta,phi
      double precision pp,tiny
      parameter (tiny=1d-12)
      pt=sqrt(p(1)**2+p(2)**2)
      y=log((p(4)+p(3))/(p(4)-p(3)))/2
      pp=sqrt(pt**2+p(3)**2)*(1+tiny)
      eta=log((pp+p(3))/(pp-p(3)))/2
      phi=atan2(p(2),p(1))
      end

      subroutine getktyphi(njets,pjet,ktjet,yjet,phijet)
      implicit none
      integer njets,j
      double precision pjet(4,njets),ktjet(njets),yjet(njets)
     $     ,phijet(njets)
      do j=1,njets
         ktjet(j)=sqrt(pjet(1,j)**2+pjet(2,j)**2)        
         yjet(j)=0.5d0*log((pjet(4,j)+pjet(3,j))/(pjet(4,j)-pjet(3,j)))
         phijet(j)=atan2(pjet(2,j),pjet(1,j))
      enddo
      end

      subroutine decvariables(pdec0,ppart0,edec,cthdec,phidec)
      implicit none
      include 'pwhg_math.h'
      double precision pdec0(4),ppart0(4),edec,cthdec,phidec
      double precision pdec(0:3),ppart(0:3)
      double precision vec(3),beta,pt
      integer mu
      do mu=1,3
         pdec(mu)=pdec0(mu)
         ppart(mu)=ppart0(mu)
      enddo
      pdec(0)=pdec0(4)
      ppart(0)=ppart0(4)
      vec(1)=0
      vec(2)=0
      vec(3)=-1
      beta=ppart(3)/ppart(0)
      call mboost(1,vec,beta,pdec,pdec)
      call mboost(1,vec,beta,ppart,ppart)
      pt=sqrt(ppart(1)**2+ppart(2)**2)
      vec(1)=ppart(1)/pt
      vec(2)=ppart(2)/pt
      vec(3)=0
      beta=-pt/ppart(0)
      call mboost(1,vec,beta,pdec,pdec)
      call mboost(1,vec,beta,ppart,ppart)
      edec=pdec(0)
      cthdec=pdec(3)/sqrt(pdec(1)**2+pdec(2)**2+pdec(3)**2)
      phidec=atan2(pdec(2),pdec(1))-atan2(vec(2),vec(1))
c bring it back between -pi and pi
      phidec=phidec-nint(phidec/(2*pi))*2*pi
      end


      logical function is_lepton(j)
      implicit none
      integer j
      is_lepton=(j.eq.11).or.(j.eq.13).or.(j.eq.15)
      end

      logical function is_antilepton(j)
      implicit none
      integer j
      is_antilepton=(j.eq.-11).or.(j.eq.-13).or.(j.eq.-15)
      end
      

      logical function bson(j)
      implicit none
      integer J
      include 'hepevt.h'
      integer jcurr,oldjcurr
      logical bhadr
      jcurr=j
      if(jcurr.eq.0) then
         bson=.false.
         return
      endif
c     check that the mother of the b quark is a top quark to avoid
c     mistaking the b quarks coming from the hard scattering. this
c     may happen in LHEF level analysis with decayed top
      if((abs(idhep(jcurr)).eq.5).and.
     $     (abs(idhep(jmohep(1,jcurr))).eq.6)) then
         bson=.true.
c     check it should never enter here after a shower
         if (isthep(jcurr).ne.1) stop "PROBLEM IN bson"
         return
      endif
 1    continue
      if(jcurr.eq.0) then
         bson=.false.
         return
      endif
      bson=.false.
      if(bhadr(idhep(jcurr))) then
         bson=.true.
         return
      endif
c      print *,jcurr,idhep(jcurr),jmohep(1,jcurr)
      oldjcurr=jcurr
      jcurr=jmohep(1,jcurr)  
      if (jmohep(1,jcurr).eq.oldjcurr) then
c     avoid entering an infinite loop 
c     (it happens with HERWIG )
         bson=.false.
         return
      endif
      if (jmohep(1,jcurr).eq.0) then
          bson=.false.
         return
      endif
      goto 1
      end

      logical function bhadr(idhep)
      implicit none
      integer idhep
      integer i1,i2,idigit
      i1=mod(idigit(1,idhep),2)
c if the rightmost digit is an odd number it is a bottomed meson (1,3,5,7)
c otherwise it is a bottomed baryon
      if(i1.eq.1) then
c if is a bottomed meson the 3rd rightmost digit should be 5
         i2=idigit(3,idhep)
      elseif(i1.eq.0) then
c if is a bottomed barion the 4th rightmost digit should be 5
         i2=idigit(4,idhep)
      endif
      if(i2.eq.5) then
         bhadr=.true.
      else
         bhadr=.false.
      endif
      end

      
      logical function sonofidhep(id,j,posj)
      implicit none
      integer id,j,posj
      include 'hepevt.h'
      integer jcurr,oldjcurr
      logical bhadr
      jcurr=j
      if(jcurr.eq.0) then
         sonofidhep=.false.
         posj=0
         return
      endif
 1    continue
      if(jcurr.eq.0) then
         sonofidhep=.false.
         posj=0
         return
      endif
      sonofidhep=.false.
      if(idhep(jcurr).eq.id) then
         sonofidhep=.true.
         posj=jcurr
         return
      endif
      oldjcurr=jcurr
      jcurr=jmohep(1,jcurr)  
      if (jmohep(1,jcurr).eq.oldjcurr) then
c     avoid entering an infinite loop
c     (it happens with HERWIG)
         sonofidhep=.false.
         return
      endif
      goto 1
      end


c     calculate the separation in the lego plot between the two momenta p1 and p2 
c     in azi and rapidity
      function rsep_azi_y(p1,p2)
      implicit none
      double precision rsep_azi_y,p1(*),p2(*)
      double precision y1,phi1,y2,phi2,kt1,kt2
      double precision deltaphi,eta1,eta2
      external deltaphi
      call ptyetaphi(p1,kt1,y1,eta1,phi1)
      call ptyetaphi(p2,kt2,y2,eta2,phi2)
      rsep_azi_y = sqrt( (y1-y2)**2 + (deltaphi(phi1,phi2))**2 )
      end


      function deltaphi(phi1,phi2)
      implicit none
      double precision deltaphi,phi1,phi2
      double precision pi,pi2
      parameter(pi = 3.141592653589793D0, pi2 = 9.869604401089358D0)
      deltaphi = abs(phi1-phi2)
c make sure it is below 2 pi
c      delphi=delphi-2*pi*int(delphi/(2*pi))
      if (deltaphi.gt.pi) then
         deltaphi = 2*pi-deltaphi
      endif
      if (deltaphi.lt.0 .or. deltaphi.gt.pi) then
         print*,' problem in deltaphi = ',deltaphi
      endif
      end

      function azi(p)
      implicit none
      double precision pi,pi2
      parameter(pi = 3.141592653589793D0, pi2 = 9.869604401089358D0)
      double precision azi,p(0:3)
      azi = atan(p(2)/p(1))
      if (p(1).lt.0d0) then
         if (azi.gt.0d0) then               
            azi = azi - pi
         else
            azi = azi + pi
         endif
      endif    
      end

      

