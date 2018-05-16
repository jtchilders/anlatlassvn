c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include 'pwhg_math.h'
      integer j,k,i
      real * 8 dy,dylep,dpt,dr
      character * 1 cnum(9)
      data cnum/'1','2','3','4','5','6','7','8','9'/
      integer maxjet
      parameter (maxjet=3)
      integer nptmin
      parameter (nptmin=3)
      character * 4 cptmin(nptmin)
      real * 8 ptminarr(nptmin)
      data cptmin/  '-001',  '-010',  '-020'/
      data ptminarr/   1d0,     10d0,    20d0/
      common/infohist/ptminarr,cnum,cptmin
      save /infohist/
      real * 8 Hmass,Hwidth,powheginput
      external powheginput

      call inihists

      dy=0.5d0
      dylep=0.4d0
      dpt=10d0
      dr=0.2d0

      Hmass = powheginput('hmass')
      Hwidth = powheginput('hwidth')

      
      do i=1,nptmin
c     total cross section sanity check

      call bookupeqbins('sigtot'//cptmin(i),1d0,0.5d0,1.5d0)
      
      call bookupeqbins('Nevents'//cptmin(i),2.5d0,0,2d2)

      call bookupeqbins('Njet'//cptmin(i),1d0,-0.5d0,5.5d0)

      call bookupeqbins('H-y'//cptmin(i),dy,-5d0,5d0)
      call bookupeqbins('H-eta'//cptmin(i),dy,-5d0,5d0)
      call bookupeqbins('H-pt'//cptmin(i),dpt,0d0,400d0)
      call bookupeqbins('H-m'//cptmin(i),Hwidth,Hmass-20*Hwidth,
     $     Hmass+20*Hwidth)
c      call bookupeqbins('H-m'//cptmin(i),0.2d-2,124.98d0,125.020d0)

      call bookupeqbins('W-y'//cptmin(i),dy,-5d0,5d0)
      call bookupeqbins('W-eta'//cptmin(i),dy,-5d0,5d0)
      call bookupeqbins('W-pt'//cptmin(i),dpt,0d0,400d0)
      call bookupeqbins('W-m'//cptmin(i),dpt,0d0,200d0)

      call bookupeqbins('lept-eta'//cptmin(i),dylep,-4d0,4d0)
      call bookupeqbins('lept-pt'//cptmin(i),dpt,0d0,500d0)
      call bookupeqbins('miss-pt'//cptmin(i),dpt,0d0,500d0)

      call bookupeqbins('HW-y'//cptmin(i),dy,-5d0,5d0)
      call bookupeqbins('HW-eta'//cptmin(i),dy,-5d0,5d0)
      call bookupeqbins('HW-pt'//cptmin(i),dpt,0d0,400d0)
      call bookupeqbins('HW-ptzoom'//cptmin(i),2d0,1d0,151d0)
      call bookupeqbins('HW-ptzoom2'//cptmin(i),0.5d0,0d0,20d0)
      call bookupeqbins('HW-ptzoom3'//cptmin(i),0.001d0,0d0,0.2d0)
      call bookupeqbins('HW-m'//cptmin(i),dpt,0d0,400d0)


      do j=1,maxjet
         call bookupeqbins('j'//cnum(j)//'-y'//cptmin(i),dy,-5d0,5d0)
         call bookupeqbins('j'//cnum(j)//'-eta'//cptmin(i),dy,-5d0,5d0)
         call bookupeqbins('j'//cnum(j)//'-pt'//cptmin(i),dpt,0d0,400d0)
         call bookupeqbins('j'//cnum(j)//'-ptzoom'//cptmin(i),
     $        2d0,1d0,151d0)
         call bookupeqbins('j'//cnum(j)//'-m'//cptmin(i),dpt,0d0,400d0) 
         call bookupeqbins('j'//cnum(j)//'-ptzoom2'//cptmin(i),
     $        0.5d0,0d0,20d0)
      enddo


      do j=1,maxjet-1
         do k=j+1,maxjet
            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1           '-y'//cptmin(i),dy,-5d0,5d0)  
            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1           '-eta'//cptmin(i),dy,-5d0,5d0)
            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1           '-pt'//cptmin(i),dpt,0d0,400d0)
            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1           '-m'//cptmin(i),dpt,0d0,400d0)  
         enddo
      enddo
  
 
      do j=1,maxjet-1
         do k=j+1,maxjet
            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1           '-dy'//cptmin(i),dy,-5d0,5d0)  
            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1           '-deta'//cptmin(i),dy,-5d0,5d0)
            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1           '-delphi'//cptmin(i),pi/20,0d0,pi)
            call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1           '-dr'//cptmin(i),dr,0d0,10d0)  
         enddo
      enddo
  
      do j=1,maxjet-1
            call bookupeqbins('j'//cnum(j)//'lept'//
     1           '-dy'//cptmin(i),dy,-5d0,5d0)  
            call bookupeqbins('j'//cnum(j)//'lept'//
     1           '-deta'//cptmin(i),dy,-5d0,5d0)
            call bookupeqbins('j'//cnum(j)//'lept'//
     1           '-delphi'//cptmin(i),pi/20,0d0,pi)
            call bookupeqbins('j'//cnum(j)//'lept'//
     1           '-dr'//cptmin(i),dr,0d0,10d0)  
      enddo

    
c      do j=1,maxjet-1
c      do k=j+1,maxjet
c         call bookupeqbins('Hj'//cnum(j)//'-j'//cnum(k)//
c     1        '-dy'//cptmin(i),dy,-5d0,5d0)  
c         call bookupeqbins('Hj'//cnum(j)//'-j'//cnum(k)//
c     1        '-deta'//cptmin(i),dy,-5d0,5d0)
c         call bookupeqbins('Hj'//cnum(j)//'-j'//cnum(k)//
c     1        '-delphi'//cptmin(i),pi/20,0d0,pi)
c         call bookupeqbins('Hj'//cnum(j)//'-j'//cnum(k)//
c     1        '-dr'//cptmin(i),dr,0d0,20d0)  
c      enddo
c      enddo

c      if(maxjet.ge.3) then
c         call bookupeqbins('Hj1j2-j3-dy'//cptmin(i),dy,-5d0,5d0)  
c         call bookupeqbins('Hj1j2-j3-deta'//cptmin(i),dy,-5d0,5d0)
c         call bookupeqbins('Hj1j2-j3-delphi'//cptmin(i),pi/20,0d0,pi)
c         call bookupeqbins('Hj1j2-j3-dr'//cptmin(i),dr,0d0,20d0)
c      endif

c      do j=1,maxjet
c         call bookupeqbins('ptrel'//cnum(j)//cptmin(i),0.5d0,0d0,20d0)
c      enddo      
c$$$
c$$$      do j=1,maxjet
c$$$         call bookupeqbins('ptrel'//cnum(j)//'qqqq'//cptmin(i),
c$$$     $        0.5d0,0d0,20d0)
c$$$         call bookupeqbins('ptrel'//cnum(j)//'qqgg'//cptmin(i),
c$$$     $        0.5d0,0d0,20d0)
c$$$         call bookupeqbins('ptrel'//cnum(j)//'ggqq'//cptmin(i),
c$$$     $        0.5d0,0d0,20d0)
c$$$         call bookupeqbins('ptrel'//cnum(j)//'gggg'//cptmin(i),
c$$$     $        0.5d0,0d0,20d0)
c$$$         call bookupeqbins('ptrel'//cnum(j)//'qgqg'//cptmin(i),
c$$$     $        0.5d0,0d0,20d0)
c$$$      enddo      
c$$$

      enddo
      end
     
      subroutine analysis(dsig0)
      implicit none
      real * 8 dsig0
      include 'hepevt.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h' 
      include 'pwhg_rad.h' 
      include 'pwhg_weights.h'
c      include 'pwhg_flg.h'
c      include 'LesHouches.h'
      integer isthep_loc(NMXHEP)  ! local copy of isthep
      logical ini
      data ini/.true./
      save ini
      integer   maxjet,mjets,njets,numjets,ntracks
      parameter (maxjet=2048)
      real * 8  ktj(maxjet),etaj(maxjet),rapj(maxjet),
     1    phij(maxjet),pj(4,maxjet),rr,ptrel(4)
      integer maxtrack
      parameter (maxtrack=2048)
      real * 8  ptrack(4,maxtrack)
      integer   jetvec(maxtrack),itrackhep(maxtrack)
      character * 1 cnum(9)
      integer nptmin
      parameter (nptmin=3)
      character * 4 cptmin(nptmin)
      real * 8 ptminarr(nptmin)      
      common/infohist/ptminarr,cnum,cptmin
      save /infohist/
      integer j,k,i,jj
c     we need to tell to this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer ih,il,inu
      real * 8 ph(4),pl(4),pnu(4),pw(4)
      real * 8 httot,y,eta,pt,m
      real * 8 dy,deta,delphi,dr
      integer ihep
      real * 8 powheginput,dotp
      external powheginput,dotp
      real * 8 ptmin
      integer idvecbos,Vdecmod,idl,idnu
      save idvecbos,Vdecmod,idl,idnu
      integer maxnumlep
      parameter (maxnumlep=10)
      real * 8 pvl(4,maxnumlep),plep(4,maxnumlep)
      integer mu,ilep,ivl,nlep,nvl
      logical is_W
      real * 8 mV2,ptvb,mvb,ptlep,ptminfastjet,ptvl,R,ylep,yvb,yvl
      real * 8 Wmass,Wwidth,Wmasslow,Wmasshigh
      integer jpart, jjet
      real * 8 palg
      integer ii
      integer  minlo
      save minlo
      data minlo/0/
      character * 20 processid
c      real * 8 rescfac1,rescfac2
c      common /crescfac/rescfac1,rescfac2
      real * 8 dsig(7)
      integer nweights
      logical inimulti
      data inimulti/.true./
      save inimulti


c      call reweightifneeded(dsig0,dsig)

      if(inimulti) then
         if(weights_num.eq.0) then
            call setupmulti(1)
         else
            call setupmulti(weights_num)
         endif
         inimulti=.false.
      endif

      dsig=0
      if(weights_num.eq.0) then
         dsig(1)=dsig0
         nweights=1
      else
         dsig(1:weights_num)=weights_val(1:weights_num)
          nweights=weights_num
      endif

      if(sum(abs(dsig)).eq.0) return

c      if(dsig.eq.0) return


      if (ini) then
         idvecbos=powheginput('idvecbos')
         Vdecmod=powheginput('vdecaymode')

         if (WHCPRG.ne.'NLO   ') then
            if (Vdecmod.eq.1) then
               idl=-11
               idnu=12
            elseif (Vdecmod.eq.2) then
               idl=-13
               idnu=14
            elseif (Vdecmod.eq.3) then
               idl=-15
               idnu=16
            endif
         else
            idl=-11
            idnu=12           
         endif
c     if idvecbos=24 idl and idnu are ok
         if (idvecbos.eq.-24) then
            idl = -idl
            idnu= -idnu
         endif

         minlo=powheginput('#minlo')
         if (minlo.eq.1) then
            processid='HW'
         else
            include 'pwhg_processid.h'
         endif
         ini=.false.
      endif


      ilep=0
      ih=0
      ivl=0

      do ihep=1,nhep  
         isthep_loc(ihep) = isthep(ihep)
      enddo
      
      if ((WHCPRG.eq.'NLO   ').or.(WHCPRG.eq.'LHE   ')) then 
         do ihep=1,nhep            
            if(idhep(ihep).eq.idl) then
               ilep=ihep
               do mu=1,4
                  plep(mu,1)=phep(mu,ihep)
               enddo
            elseif(idhep(ihep).eq.idnu) then
               ivl=ihep
               do mu=1,4
                  pvl(mu,1)=phep(mu,ihep)
               enddo              
            elseif(idhep(ihep).eq.25) then
               ih=ihep
               do mu=1,4
                  ph(mu)=phep(mu,ihep)
               enddo              
            endif
         enddo
      endif


c     Analysis after MC shower
      if((WHCPRG.eq.'HERWIG').or.(WHCPRG.eq.'PYTHIA')) then
c     Loop again over final state particles to find products of W decay, by
c     looking into the shower branchings.
         nlep=0
         nvl=0
         do ihep=1,nhep
c     works for POWHEG+HERWIG, POWHEG+PYHIA, HERWIG, PYTHIA and real in
c     MC@NLO
            if(idhep(ihep).eq.25.and.isthep_loc(ihep).eq.1) then
               ph=phep(1:4,ihep)
               ih=ihep
            endif
            if (isthep_loc(ihep).eq.1.and.(idhep(ihep).eq.idl .or.
     $           idhep(ihep).eq.idnu))
     1           then
c               is_W = idhep(jmohep(1,jmohep(1,ihep))).eq.idvecbos .or.
c     $         idhep(jmohep(1,jmohep(1,jmohep(1, ihep)))).eq.idvecbos
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     if the hadronization is switched off, then the electrons and neutrinos
c     that are found come only from the W decay
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               is_W = .true.
               if (is_W) then
c     find first decay product
                  if(idhep(ihep).eq.idl) then
                     ilep=ihep
                     nlep=nlep+1
c     find second decay product
                  elseif(idhep(ihep).eq.idnu) then
                     ivl=ihep
                     nvl=nvl+1
                  endif
               endif
            endif
         enddo
         if(nvl.ne.1.or.nlep.ne.1) then
            write(*,*) 'Problems with leptons from W decay'
c            write(*,*) 'PROGRAM ABORT'
            write(*,*) 'nvl= ',nvl, 'nlep= ',nlep
c            call exit(1)
            return
         endif
            
         do mu=1,4
            plep(mu,nlep)=phep(mu,ilep)
            pvl(mu,nvl)=phep(mu,ivl)
         enddo
      endif

      if (ilep*ih*ivl.eq.0) then
         write(*,*) 
     $        'ERROR... have NOT found the electron/neutrino/Higgs'
         write(*,*) 'idhep'
         write(*,*) (idhep(ihep),ihep=1,nhep)
         return
      endif
         
c     change status of l vu and Higgs
      isthep_loc(ilep)=10000
      isthep_loc(ivl)=10000
      isthep_loc(ih)=10000

c     W momentum
      do mu=1,4
         pw(mu)=plep(mu,1) + pvl(mu,1)
      enddo

c     set up arrays for jet finding
c      do jpart=1,maxtrack
c         do mu=1,4
c            ptrack(mu,jpart)=0d0
c         enddo
c         jetvec(jpart)=0
c      enddo      
c      do jjet=1,maxjet
c         do mu=1,4
c            pj(mu,jjet)=0d0
c         enddo
c      enddo

      ntracks=0
      mjets=0
c     Loop over final state particles to find jets 
      do ihep=1,nhep
         if (isthep_loc(ihep).eq.1) then
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
      if (ntracks.eq.0) then
         numjets=0
      else
c     palg=1 is standard kt, -1 is antikt
         palg = -1d0
         R = 0.5d0              ! radius parameter
         ptminfastjet = 1d0
         call fastjetppgenkt(ptrack,ntracks,R,palg,ptminfastjet,
     $        pj,numjets,jetvec)
c         call fastjetktwhich(ptrack,ntracks,ptminfastjet,R,
c     $        pj,mjets,jetvec) 
      endif

      do i=1,nptmin        
         njets=0
         
         do j=1,min(3,numjets)
            ktj(j) = sqrt(pj(1,j)**2 + pj(2,j)**2 )
            if (ktj(j).gt.ptminarr(i)) then
               njets=njets+1
            endif
         enddo
         
c     since ptminarr(1) is the smallest value, the following return is correct
         if (processid.eq.'HWJ') then
            if(njets.eq.0) return
         endif
         
         call filld('sigtot'//cptmin(i),1d0,dsig)
         
         call filld('Nevents'//cptmin(i),abs(dsig),1d0)
         
         if(njets.eq.0) then
            call filld('Njet'//cptmin(i),0d0,dsig)
         elseif(njets.eq.1) then
            call filld('Njet'//cptmin(i),1d0,dsig)
         elseif(njets.eq.2) then
            call filld('Njet'//cptmin(i),2d0,dsig)
         elseif(njets.eq.3) then
            call filld('Njet'//cptmin(i),3d0,dsig)
         elseif(njets.eq.4) then
            call filld('Njet'//cptmin(i),4d0,dsig)
         elseif(njets.eq.5) then
            call filld('Njet'//cptmin(i),5d0,dsig)
         else
c     write(*,*) ' Njet?',mjets
         endif
c     Higgs
         call getyetaptmass(ph,y,eta,pt,m)
         call filld('H-y'//cptmin(i),    y, dsig)
         call filld('H-eta'//cptmin(i),eta, dsig)
         call filld('H-pt'//cptmin(i),  pt, dsig)
         call filld('H-m'//cptmin(i), m, dsig)
c     W
         call getyetaptmass(pw,y,eta,pt,m)
         call filld('W-y'//cptmin(i),    y, dsig)
         call filld('W-eta'//cptmin(i),eta, dsig)
         call filld('W-pt'//cptmin(i),  pt, dsig)
         call filld('W-m'//cptmin(i), m, dsig)
c     lepton
         call getyetaptmass(plep(:,1),y,eta,pt,m)
         call filld('lept-eta'//cptmin(i),eta, dsig)
         call filld('lept-pt'//cptmin(i),  pt, dsig)
c     neutrino
         call getyetaptmass(pvl(:,1),y,eta,pt,m)
         call filld('miss-pt'//cptmin(i),  pt, dsig)
c     HW
         call getyetaptmass(ph+pw,y,eta,pt,m)

c         if (pt.gt.0.09d0 .and. pt.lt.0.11d0) then
c            write(*,*) 'ci siamo'
c            write(*,*) rescfac1,rescfac2
c         endif


         call filld('HW-y'//cptmin(i),    y, dsig)
         call filld('HW-eta'//cptmin(i),eta, dsig)
         call filld('HW-pt'//cptmin(i),  pt, dsig)
         call filld('HW-ptzoom'//cptmin(i),  pt, dsig)
         call filld('HW-ptzoom2'//cptmin(i),  pt, dsig)
         call filld('HW-ptzoom3'//cptmin(i),  pt, dsig)
         call filld('HW-m'//cptmin(i), m, dsig)
         
c     jets
         mjets=min(njets,2)
         
         do j=1,mjets
            call getyetaptmass(pj(:,j),y,eta,pt,m)
            call filld('j'//cnum(j)//'-y'//cptmin(i),     y, dsig)
            call filld('j'//cnum(j)//'-eta'//cptmin(i), eta, dsig)
            call filld('j'//cnum(j)//'-pt'//cptmin(i),   pt, dsig)
            call filld('j'//cnum(j)//'-ptzoom'//cptmin(i),   pt, dsig)
            call filld('j'//cnum(j)//'-ptzoom2'//cptmin(i),   pt, dsig)
            call filld('j'//cnum(j)//'-m'//cptmin(i),     m, dsig)
c     call filld('ptrel'//cnum(j)//cptmin(i),ptrel(j), dsig)         
         enddo
         
         do j=1,mjets
            do k=j+1,mjets
               call getyetaptmass(pj(:,j)+pj(:,k),y,eta,pt,m)
               call filld('j'//cnum(j)//'j'//cnum(k)//'-y'//cptmin(i),
     $              y, dsig)
               call filld('j'//cnum(j)//'j'//cnum(k)//'-eta'//cptmin(i),
     $              eta, dsig)
               call filld('j'//cnum(j)//'j'//cnum(k)//'-pt'//cptmin(i),
     $              pt, dsig)
               call filld('j'//cnum(j)//'j'//cnum(k)//'-m'//cptmin(i), 
     $              m, dsig)
            enddo
         enddo
         
         do j=1,mjets
            do k=j+1,mjets
               call deltaplot(pj(:,j),pj(:,k),dsig,
     1              'j'//cnum(j)//'j'//cnum(k),cptmin(i))
            enddo
         enddo
         
         do j=1,mjets
            call deltaplot(pj(:,j),plep(:,1),dsig,
     1           'j'//cnum(j)//'lept',cptmin(i))
         enddo
         
      enddo
      end




c      subroutine yetaptmassplot(p,dsig,prefix)
c      implicit none
c      real * 8 p(4),dsig
c      character *(*) prefix
c      real * 8 y,eta,pt,m
c      call getyetaptmass(p,y,eta,pt,m)
c      call filld(prefix//'-y',y,dsig)
c      call filld(prefix//'-eta',eta,dsig)
c      call filld(prefix//'-pt',pt,dsig)
c      call filld(prefix//'-m',m,dsig)
c      end

      subroutine deltaplot(p1,p2,dsig,prefix,postfix)
      implicit none
      real * 8 p1(4),p2(4),dsig(7)
      character *(*) prefix,postfix
      real * 8 dy,deta,delphi,dr
      call getdydetadphidr(p1,p2,dy,deta,delphi,dr)
      call filld(prefix//'-dy'//postfix,dy,dsig)
      call filld(prefix//'-deta'//postfix,deta,dsig)
      call filld(prefix//'-delphi'//postfix,delphi,dsig)
      call filld(prefix//'-dr'//postfix,dr,dsig)
      end

      subroutine getyetaptmass(p,y,eta,pt,mass)
      implicit none
      real * 8 p(4),y,eta,pt,mass
      call pwhg_getrapidity(p,y)      
      pt=sqrt(p(1)**2+p(2)**2)
      call pwhg_getpseudorapidity(p,eta)
      call pwhg_getinvmass(p,mass)
      end

      subroutine getdydetadphidr(p1,p2,dy,deta,dphi,dr)
      implicit none
      include 'pwhg_math.h' 
      real * 8 p1(*),p2(*),dy,deta,dphi,dr
      real * 8 y1,eta1,pt1,mass1,phi1
      real * 8 y2,eta2,pt2,mass2,phi2
      call getyetaptmass(p1,y1,eta1,pt1,mass1)
      call getyetaptmass(p2,y2,eta2,pt2,mass2)
      dy=y1-y2
      deta=eta1-eta2
      phi1=atan2(p1(1),p1(2))
      phi2=atan2(p2(1),p2(2))
      dphi=abs(phi1-phi2)
      dphi=min(dphi,2d0*pi-dphi)
      dr=sqrt(deta**2+dphi**2)
      end


      subroutine sortbypt(n,iarr)
      implicit none
      integer n,iarr(n)
      include 'hepevt.h'
      integer j,k
      real * 8 tmp,pt(nmxhep)
      logical touched
      do j=1,n
         pt(j)=sqrt(phep(1,iarr(j))**2+phep(2,iarr(j))**2)
      enddo
c bubble sort
      touched=.true.
      do while(touched)
         touched=.false.
         do j=1,n-1
            if(pt(j).lt.pt(j+1)) then
               k=iarr(j)
               iarr(j)=iarr(j+1)
               iarr(j+1)=k
               tmp=pt(j)
               pt(j)=pt(j+1)
               pt(j+1)=tmp
               touched=.true.
            endif
         enddo
      enddo
      end

      function islept(j)
      implicit none
      logical islept
      integer j
      if(abs(j).ge.11.and.abs(j).le.15) then
         islept = .true.
      else
         islept = .false.
      endif
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

      subroutine pwhgfinalopshist
      end
