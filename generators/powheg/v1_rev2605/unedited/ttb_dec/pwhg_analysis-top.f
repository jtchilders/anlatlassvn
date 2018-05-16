C ----------------------------------------------------------------- C
C - This is a parton level only analysis for ttbar production     - C
C - tops bs Ws etc are constructed using MC truth only - no jets  - C
C - The MC truth reconstruction has been tested (see sanity check - C
C - code in the analysis below).                                  - C
C - Since it is parton level you need to comment out from         - C
C - CALL HWDHOB down to CALL HWDHOB inclusive in main-HERWIG.f .  - C
C - Also at some point, the showering went into what looked like  - C
C - an infinite loop after 137K events - gdb said it was in       - C
C - HWHGUP. The same glitch did not occur with *** herwig6520.f *** C
C - I also eliminated the analysis as a possible cause (it occurs - C
C - with HWANAL removed). I did not see anything fishy with the   - C
C - Tevatron, semileptonic event that got caught.                 - C
C ----------------------------------------------------------------- C

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
      integer j,l
      character * 20 prefix
      integer nbins
      parameter (nbins=11)
      real * 8 pT_tt_bins(nbins+1)
      data pT_tt_bins/  0d0, 10d0, 25d0, 50d0,100d0,
     1     150d0,200d0,250d0,300d0,400d0,600d0, 900d0/          
      real * 8 m_tt_bins(nbins+1)
      data m_tt_bins/ 320d0,360d0,400d0,450d0,500d0,
     1     550d0,600d0,650d0,700d0,800d0,900d0,1000d0/          
      character * 2 digit(20)
      data digit/'01','02','03','04','05','06','07','08','09','10',
     1           '11','12','13','14','15','16','17','18','19','20'/
      integer lenocc
      external lenocc

      call inihists

      call bookupeqbins('Njets-pt10'  ,1d0,-0.5d0,10.5d0)
      call bookupeqbins('Njets-pt20'  ,1d0,-0.5d0,10.5d0)
      call bookupeqbins('Njets-pt30'  ,1d0,-0.5d0,10.5d0)
      call bookupeqbins('Njets-pt40'  ,1d0,-0.5d0,10.5d0)

      do j=1,9
         if(j.eq.1) then
            prefix='t'
         elseif(j.eq.2) then
            prefix='tb'
         elseif(j.eq.3) then
            prefix='btop'
         elseif(j.eq.4) then
            prefix='bbtop'
         elseif(j.eq.5) then
            prefix='lwp'
         elseif(j.eq.6) then
            prefix='lwm'
         elseif(j.eq.7) then
            prefix='ttb'
         elseif(j.eq.8) then
            prefix='ttb-radPY'
         elseif(j.eq.9) then
            prefix='ttb-radPW'
         endif
         l=lenocc(prefix)
         call bookupeqbins(prefix(1:l)//'_y'  ,0.2d0,-4d0,4d0)
         call bookupeqbins(prefix(1:l)//'_eta',0.2d0,-4d0,4d0)
         call bookupeqbins(prefix(1:l)//'_pt' ,2d0,0d0,400d0)
         call bookupeqbins(prefix(1:l)//'_m'  ,2.5d0,0d0,500d0)
      enddo
      
      call bookupeqbins('m_lp_lm',2d0,0d0,400d0)
      call bookupeqbins('mT_lp_MET',2d0,0d0,400d0)
      call bookupeqbins('mT_lm_MET',2d0,0d0,400d0)
      call bookupeqbins('m_lp_jb',2d0,0d0,400d0)
      call bookupeqbins('m_lm_jbbar',2d0,0d0,400d0)
      call bookupeqbins('m_wp_b',2d0,0d0,400d0)
      call bookupeqbins('m_wm_bb',2d0,0d0,400d0)
      call bookupeqbins('m_wp_bj',2d0,0d0,400d0)
      call bookupeqbins('m_wm_bbj',2d0,0d0,400d0)
      call bookupeqbins('bfrag',0.01d0,0d0,1d0)
      call bookupeqbins('bptdec',1d0,0d0,90d0)
      call bookupeqbins('bbptdec',1d0,0d0,90d0)
      call bookupeqbins('bbfrag',0.01d0,0d0,1d0)
      call bookupeqbins('wpmom',0.01d0,0d0,1d0)
      call bookupeqbins('wmmom',0.01d0,0d0,1d0)
      call bookupeqbins('lwp-lwm-dy',0.2d0,-4d0,4d0)
      call bookupeqbins('lwp-lwm-deta',0.2d0,-4d0,4d0)
      call bookupeqbins('lwp-lwm-delphi',pi/10,0d0,pi)
      call bookupeqbins('lwp-lwm-dr',0.2d0,0d0,8d0)
      do j=1,10
         call bookupeqbins('cth1cth2-'//digit(j),0.2d0,-1d0,1d0)
      enddo
      end

      subroutine analysis(dsig0)
      implicit none
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include 'LesHouches.h'
      character * 6 whcprg      
      common/cwhcprg/whcprg
      integer jpref
      character * 20 prefix(18)
      common/ccccprefix/jpref,prefix
      data whcprg/'NLO   '/
      real * 8  dsig0,dsig
      logical   ini
      data      ini/.true./
      save      ini
      integer   ihep                ! HEPEVT index.
      real * 8 p_top(4),p_tb(4),p_wp(4),p_wm(4),p_lwp(4),p_lwm(4),
     1         p_nuwp(4),p_nuwm(4),p_b(4),p_bb(4),y,eta,pt,mass,
     2         ptzmf(4),plzmf(4)
      integer   maxtracks,maxjets
      parameter (maxtracks=nmxhep,maxjets=20)
      integer mjets,jetvec(maxtracks)
      logical   isForClustering(maxtracks)
      real * 8 j_kt(maxjets),j_eta(maxjets),j_rap(maxjets),
     1     j_phi(maxjets),j_p(4,maxjets)
      integer j,id,i_top,i_atop,i_bfromtop,i_abfromatop,
     1     i_wp,i_wm,i_lwp,i_lwm,i_nuwp,i_nuwm,i_bjet,i_abjet,jhep,
     1     i_part,njets20,njets30,njets40
      real * 8 mtop,mtb,mwp,mwm,mb,mbb,p_bmax,e_bmax,xb,
     1     p_bbmax,e_bbmax,xbb,ewp,pwp,ewm,pwm,xw,
     2     dy,deta,dphi,dr,cth1,cth2
      integer jcth1
      real * 8 w(4),pb(4),ptb
      real * 8 prodvec2,powheginput
      logical sonofid
      external sonofid
      integer in_jet
      external in_jet
      integer ngenerations,inotfound,iprodrad
      common/cngenerations/ngenerations
      character * 2 digit(20)
      data digit/'01','02','03','04','05','06','07','08','09','10',
     1           '11','12','13','14','15','16','17','18','19','20'/
      integer id1,id2
      ngenerations = powheginput("#ngenerations")
      if(ngenerations.lt.0) ngenerations = 4

      dsig  = dsig0

      i_top = 0
      i_atop = 0
      i_wp = 0
      i_wm = 0
      i_lwp = 0
      i_lwm = 0
      i_bfromtop = 0
      i_abfromatop = 0

      if(whcprg.eq.'NLO') then
         i_top = 3
         i_atop = 4
         i_wp = 5
         i_wm = 6

         i_lwp = 7
         i_lwm = 9
         i_nuwp = 8
         i_nuwm = 10
         i_bfromtop = 11
         i_abfromatop = 12

         IsForClustering = .false.
         IsForClustering(13) = .true.
         IsForClustering(i_bfromtop) = .true.
         IsForClustering(i_abfromatop) = .true.
C --------------------------------------------- C
C - LHE PARTICLE TOP RECONSTRUCTION: MC TRUTH - C
C --------------------------------------------- C
      else
c Build top MC; find the last top (tbar)
c in the event record, i.e. before decays
         do jhep=1,nhep
            id=idhep(jhep)
            if(idhep(jhep).eq.6) i_top = jhep
            if(idhep(jhep).eq.-6) i_atop = jhep
            id=abs(id)
            if(id.eq.5.or.id.eq.24) then
            if(sonofid(6,jhep)) then
               if(idhep(jhep).eq.5) i_bfromtop = jhep
               if(idhep(jhep).eq.-5) i_abfromatop = jhep
               if(idhep(jhep).eq.24) i_wp = jhep
               if(idhep(jhep).eq.-24) i_wm = jhep
            endif
            endif
            if(id.ge.11.and.id.le.14) then
            if(sonofid(24,jhep)) then
               if(idhep(jhep).eq.-11.or.idhep(jhep).eq.-13) i_lwp = jhep
               if(idhep(jhep).eq.11.or.idhep(jhep).eq.13) i_lwm = jhep
               if(idhep(jhep).eq.-12.or.idhep(jhep).eq.-14)i_nuwm = jhep
               if(idhep(jhep).eq.12.or.idhep(jhep).eq.14) i_nuwp = jhep
            endif
            endif
c for jets, using only final state particles excluding leptons
            if(isthep(jhep).eq.1.and.
     1           (idhep(jhep).lt.11.or.idhep(jhep).gt.16)) then
               IsForClustering(jhep) = .true.
            else
               IsForClustering(jhep) = .false.
            endif
         enddo
      endif

      inotfound = 0
      if(i_top.eq.0) then
         write(*,*) 'top not found'
         inotfound = inotfound + 1
      endif
      if(i_atop.eq.0) then
         write(*,*) 'antitop not found'
         inotfound = inotfound + 1
      endif
      if(i_wp.eq.0) then
         write(*,*) 'wp not found'
         inotfound = inotfound + 1
      endif
      if(i_wm.eq.0) then
         write(*,*) 'wm not found'
         inotfound = inotfound + 1
      endif
      if(i_lwp.eq.0) then
         write(*,*) 'lwp not found'
         inotfound = inotfound + 1
      endif
      if(i_lwm.eq.0) then
         write(*,*) 'lwm not found'
         inotfound = inotfound + 1
      endif
      if(i_nuwp.eq.0) then
         write(*,*) 'nuwp not found'
         inotfound = inotfound + 1
      endif
      if(i_nuwm.eq.0) then
         write(*,*) 'nuwm not found'
         inotfound = inotfound + 1
      endif
      if(i_bfromtop.eq.0) then
         write(*,*) 'b from top not found'
         inotfound = inotfound + 1
      endif
      if(i_abfromatop.eq.0) then
         write(*,*) 'bbar from tbar not found'
         inotfound = inotfound + 1
      endif

      if(inotfound.gt.0) return

c 
      if(whcprg.ne.'NLO'.and.whcprg.ne.'LHE') then
c Setup a flag:
c 1 for events with no production radiation in the LHE
c 2 for events with production radiation
c 0 otherwise
         if((pup(1,3)+pup(1,4))**2+(pup(2,3)+pup(2,4))**2.lt.1d-2) then
            iprodrad=1
c            if(scalup.gt.46) then
c               write(*,*) ' warning: scalup = ',scalup,
c     1              ' in radiation in decay'
c            endif
         else
            iprodrad=2
         endif
      else
         iprodrad = 0
      endif

      id1=idup(1)
      id2=idup(2)
      if(id1.eq.21) id1=0
      if(id2.eq.21) id2=0

      if(whcprg.eq.'LHE') then
         if(powheginput('#subprocess').eq.1) then
            if(id1.eq.0.and.id2.eq.0) then
               continue
            else
               return
            endif
         elseif(powheginput('#subprocess').eq.2) then
            if(id1.gt.0.and.id2.lt.0) then
               continue
            else
               return
            endif
         elseif(powheginput('#subprocess').eq.3) then
            if(id1.lt.0.and.id2.gt.0) then
               continue
            else
               return
            endif
         elseif(powheginput('#subprocess').eq.4) then
            if(id1.gt.0.and.id2.eq.0) then
               continue
            else
               return
            endif
         elseif(powheginput('#subprocess').eq.5) then
            if(id1.lt.0.and.id2.eq.0) then
               continue
            else
               return
            endif
         elseif(powheginput('#subprocess').eq.6) then
            if(id1.eq.0.and.id2.gt.0) then
               continue
            else
               return
            endif
         elseif(powheginput('#subprocess').eq.7) then
            if(id1.eq.0.and.id2.lt.0) then
               continue
            else
               return
            endif
         endif
      endif
               

      p_top=phep(1:4,i_top)
      p_tb=phep(1:4,i_atop)
      p_wp=phep(1:4,i_wp)
      p_wm=phep(1:4,i_wm)
      p_lwp=phep(1:4,i_lwp)
      p_lwm=phep(1:4,i_lwm)
      p_nuwp=phep(1:4,i_nuwp)
      p_nuwm=phep(1:4,i_nuwm)
      p_b=phep(1:4,i_bfromtop)
      p_bb=phep(1:4,i_abfromatop)

      mjets = maxjets
      call buildjets(mjets,j_kt,j_eta,j_rap,j_phi,j_p,jetvec,
     1     isForClustering)

      i_bjet = in_jet(i_bfromtop,jetvec)
      i_abjet = in_jet(i_abfromatop,jetvec)

      njets20 = 0
      njets30 = 0
      njets40 = 0
      do j=1,mjets
         if(j_kt(j).gt.20) then
            njets20 = njets20 + 1
         endif
         if(j_kt(j).gt.30) then
            njets30 = njets30 + 1
         endif
         if(j_kt(j).gt.40) then
            njets40 = njets40 + 1
         endif
      enddo

      call filld('Njets-pt10',dble(mjets),dsig)
      call filld('Njets-pt20',dble(njets20),dsig)
      call filld('Njets-pt30',dble(njets30),dsig)
      call filld('Njets-pt40',dble(njets40),dsig)

      call yetaptmassplot(p_top,dsig,'t')

      call yetaptmassplot(p_tb,dsig,'tb')

      call yetaptmassplot(p_b,dsig,'btop')

      call yetaptmassplot(p_bb,dsig,'bbtop')

      call yetaptmassplot(p_lwp,dsig,'lwp')

      call yetaptmassplot(p_lwm,dsig,'lwm')


      call yetaptmassplot(p_top+p_tb,dsig,'ttb')


      if(iprodrad.eq.1) then
         call yetaptmassplot(p_top+p_tb,dsig,'ttb-radPY')
      elseif(iprodrad.eq.2) then
         call yetaptmassplot(p_top+p_tb,dsig,'ttb-radPW')
      endif

      call getyetaptmass(p_lwp+p_lwm,y,eta,pt,mass)
      call filld('m_lp_lm',mass,dsig)

      call getyetaptmass(p_lwp+p_nuwp,y,eta,pt,mass)
      call filld('mT_lp_MET',mass,dsig)

      call getyetaptmass(p_lwm+p_nuwm,y,eta,pt,mass)
      call filld('mT_lm_MET',mass,dsig)

      if(i_bjet.ne.0) then
         call getyetaptmass(p_lwp+j_p(:,i_bjet),y,eta,pt,mass)
         call filld('m_lp_jb',mass,dsig)
      endif

      if(i_abjet.ne.0) then
         call getyetaptmass(p_lwm+j_p(:,i_abjet),y,eta,pt,mass)
         call filld('m_lm_jbbar',mass,dsig)
      endif

c b W mass
      call getyetaptmass(p_wp+p_b,y,eta,pt,mass)
      call filld('m_wp_b',mass,dsig)

c bb W- mass
      call getyetaptmass(p_wm+p_bb,y,eta,pt,mass)
      call filld('m_wm_bb',mass,dsig)

c b-jet W mass
      if(i_bjet.ne.0) then
         call getyetaptmass(p_wp+j_p(:,i_bjet),y,eta,pt,mass)
         call filld('m_wp_bj',mass,dsig)
      endif

c bb-jet W- mass
      if(i_abjet.ne.0) then
         call getyetaptmass(p_wm+j_p(:,i_abjet),y,eta,pt,mass)
         call filld('m_wm_bbj',mass,dsig)
      endif
      
c b frag: p_top.p_b/(p_top.p_b_max)
      mtop=sqrt(p_top(4)**2-p_top(1)**2-p_top(2)**2-p_top(3)**2)
      mwp=sqrt(p_wp(4)**2-p_wp(1)**2-p_wp(2)**2-p_wp(3)**2)
      mb=sqrt(abs(p_b(4)**2-p_b(1)**2-p_b(2)**2-p_b(3)**2))
      p_bmax=sqrt((mtop**2-(mb+mwp)**2)*(mtop**2-(mb-mwp)**2))/(2*mtop)
      e_bmax=sqrt(p_bmax**2+mb**2)
      xb=(p_top(4)*p_b(4)-p_top(1)*p_b(1)
     1     -p_top(2)*p_b(2)-p_top(3)*p_b(3))/(mtop*e_bmax)
      call filld('bfrag',xb,dsig)
      call boost2reson4(p_top,1,p_wp,w)
      call boost2reson4(p_top,1,p_b,pb)
      ptb=sqrt( abs( (pb(1)**2+pb(2)**2+pb(3)**2)
     1     -(pb(1)*w(1)+pb(2)*w(2)+pb(3)*w(3))**2/
     1     (w(1)**2+w(2)**2+w(3)**2)))
      call filld('bptdec',ptb,dsig)

c bbar frag: p_top.p_bb/(p_top.p_bb_max)
      mtb=sqrt(p_tb(4)**2-p_tb(1)**2-p_tb(2)**2-p_tb(3)**2)
      mwm=sqrt(p_wm(4)**2-p_wm(1)**2-p_wm(2)**2-p_wm(3)**2)
      mbb=sqrt(abs(p_bb(4)**2-p_bb(1)**2-p_bb(2)**2-p_bb(3)**2))
      p_bbmax=sqrt((mtb**2-(mbb+mwm)**2)*(mtb**2-(mbb-mwm)**2))/(2*mtb)
      e_bbmax=sqrt(p_bbmax**2+mbb**2)
      xbb=(p_tb(4)*p_bb(4)-p_tb(1)*p_bb(1)
     1     -p_tb(2)*p_bb(2)-p_tb(3)*p_bb(3))/(mtb*e_bbmax)
      call filld('bbfrag',xbb,dsig)
      call boost2reson4(p_tb,1,p_wm,w)
      call boost2reson4(p_tb,1,p_bb,pb)
      ptb=sqrt( abs( (pb(1)**2+pb(2)**2+pb(3)**2)
     1     -(pb(1)*w(1)+pb(2)*w(2)+pb(3)*w(3))**2/
     1     (w(1)**2+w(2)**2+w(3)**2)) )
      call filld('bbptdec',ptb,dsig)

c W momentum relative to its maximum in top rest frame
      ewp=(p_top(4)*p_wp(4)-p_top(1)*p_wp(1)
     1     -p_top(2)*p_wp(2)-p_top(3)*p_wp(3))/mtop
      pwp=sqrt(ewp**2-mwp**2)
      xw=pwp/p_bmax
      call filld('wpmom',xw,dsig)

c W momentum relative to its maximum in top rest frame
      ewm=(p_tb(4)*p_wm(4)-p_tb(1)*p_wm(1)
     1     -p_tb(2)*p_wm(2)-p_tb(3)*p_wm(3))/mtb
      pwm=sqrt(ewm**2-mwm**2)
      xw=pwm/p_bbmax
      call filld('wmmom',xw,dsig)

      call deltaplot(p_lwp,p_lwm,dsig,'lwp-lwm','')

c theta1 and theta2 from Bernreuther
      call boost2reson4(p_top+p_tb,1,p_top,ptzmf)
      call boost2reson4(p_top+p_tb,1,p_lwp,plzmf)
      call boost2reson4(ptzmf,1,plzmf,w)
      cth1 = (w(1)*ptzmf(1)+w(2)*ptzmf(2)+w(3)*ptzmf(3))/
     1 sqrt((w(1)**2+w(2)**2+w(3)**2)
     2     *(ptzmf(1)**2+ptzmf(2)**2+ptzmf(3)**2))


      call boost2reson4(p_top+p_tb,1,p_tb,ptzmf)
      call boost2reson4(p_top+p_tb,1,p_lwm,plzmf)
      call boost2reson4(ptzmf,1,plzmf,w)
      cth2 = (w(1)*ptzmf(1)+w(2)*ptzmf(2)+w(3)*ptzmf(3))/
     1 sqrt((w(1)**2+w(2)**2+w(3)**2)
     2     *(ptzmf(1)**2+ptzmf(2)**2+ptzmf(3)**2))


      jcth1 = int((cth1+1) * 5) + 1
      if(jcth1.eq.0) jcth1=1
      if(jcth1.eq.11) jcth1=10
      call filld("cth1cth2-"//digit(jcth1),cth2,dsig)

      end


      function in_jet(i_part,jetvec)
      implicit none
      include 'hepevt.h'
      integer   maxtracks,maxjets
      parameter (maxtracks=nmxhep,maxjets=20)
      integer in_jet,jetvec(maxtracks),i_part
      integer j
      logical sonofhep
      external sonofhep
      do j=1,nhep
         if(jetvec(j).ne.0) then
            if(sonofhep(i_part,j)) then
               in_jet = jetvec(j)
               return
            endif
         endif
      enddo
      in_jet = 0
      end

      function prodvec2(vec1,vec2)
      implicit none
      real * 8 prodvec2,vec1(4),vec2(4)
      prodvec2=vec1(4)*vec2(4)-vec1(1)*vec2(1)
     1 -vec1(2)*vec2(2)-vec1(3)*vec2(3)
      end


      subroutine getyetaptmass(p,y,eta,pt,mass)
      implicit none
      real * 8 p(4),y,eta,pt,mass,pv
      real *8 tiny
      parameter (tiny=1.d-5)
      y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      pt=sqrt(p(1)**2+p(2)**2)
      pv=sqrt(pt**2+p(3)**2)
      if(pt.lt.tiny)then
         eta=sign(1.d0,p(3))*1.d8
      else
         eta=0.5d0*log((pv+p(3))/(pv-p(3)))
      endif
      mass=sqrt(abs(p(4)**2-pv**2))
      end


      subroutine yetaptmassplot(p,dsig,prefix)
      implicit none
      real * 8 p(4),dsig
      character *(*) prefix
      real * 8 y,eta,pt,m
      call getyetaptmass(p,y,eta,pt,m)
      call filld(prefix//'_y',y,dsig)
      call filld(prefix//'_eta',eta,dsig)
      call filld(prefix//'_pt',pt,dsig)
      call filld(prefix//'_m',m,dsig)
      end


      subroutine pwhgfinalopshist
      implicit none
      include 'pwhg_bookhist-new.h'
      integer  f_idx,b_idx,x_idx,ixx,jxx,j,k,l
      integer  indexhist
      real * 8 num,numerr,den,denerr

      f_idx=indexhist('dF-dpT')
      b_idx=indexhist('dB-dpT')
      x_idx=indexhist('dAFB-dpT')
      call get_afb_hist(f_idx,b_idx,x_idx)

      f_idx=indexhist('dF-dpT')
      b_idx=indexhist('dB-dpT')
      x_idx=indexhist('F-gt-pT')
      call integratehist(f_idx,x_idx,1)
      x_idx=indexhist('B-gt-pT')
      call integratehist(b_idx,x_idx,1)
      x_idx=indexhist('F-lt-pT')
      call integratehist(f_idx,x_idx,-1)
      x_idx=indexhist('B-lt-pT')
      call integratehist(b_idx,x_idx,-1)

      f_idx=indexhist('F-gt-pT')
      b_idx=indexhist('B-gt-pT')
      x_idx=indexhist('AFB-gt-pT')
      call get_afb_hist(f_idx,b_idx,x_idx)

      f_idx=indexhist('F-lt-pT')
      b_idx=indexhist('B-lt-pT')
      x_idx=indexhist('AFB-lt-pT')
      call get_afb_hist(f_idx,b_idx,x_idx)

      f_idx=indexhist('dF-dm_tt')
      b_idx=indexhist('dB-dm_tt')
      x_idx=indexhist('dAFB-dm_tt')
      call get_afb_hist(f_idx,b_idx,x_idx)

      f_idx=indexhist('dF-dm_tt')
      b_idx=indexhist('dB-dm_tt')
      x_idx=indexhist('F-gt-m_tt')
      call integratehist(f_idx,x_idx,1)
      x_idx=indexhist('B-gt-m_tt')
      call integratehist(b_idx,x_idx,1)
      x_idx=indexhist('F-lt-m_tt')
      call integratehist(f_idx,x_idx,-1)
      x_idx=indexhist('B-lt-m_tt')
      call integratehist(b_idx,x_idx,-1)

      f_idx=indexhist('F-gt-m_tt')
      b_idx=indexhist('B-gt-m_tt')
      x_idx=indexhist('AFB-gt-m_tt')
      call get_afb_hist(f_idx,b_idx,x_idx)

      f_idx=indexhist('F-lt-m_tt')
      b_idx=indexhist('B-lt-m_tt')
      x_idx=indexhist('AFB-lt-m_tt')
      call get_afb_hist(f_idx,b_idx,x_idx)

      end


      subroutine integratehist(integrand,integral,direction)
      implicit none
      include 'pwhg_bookhist-new.h'
      integer  integrand,integral,direction
      integer  indexhist,ixx,jxx
      external indexhist
      if(direction.eq.1) then
        do ixx=1,nbins(integrand)
          do jxx=ixx,nbins(integrand)
            yhistarr2(ixx,integral)=yhistarr2(ixx,integral)
     1         +yhistarr2(jxx,integrand)
     1         *(xhistarr(jxx+1,integrand)-xhistarr(jxx,integrand))
            errhistarr2(ixx,integral)=errhistarr2(ixx,integral)
     1         +(errhistarr2(jxx,integrand)
     1         *(xhistarr(jxx+1,integrand)-xhistarr(jxx,integrand)))**2
          enddo
          errhistarr2(ixx,integral)=sqrt(errhistarr2(ixx,integral))
        enddo
      else if(direction.eq.-1) then
        do ixx=1,nbins(integrand)
          do jxx=ixx,1,-1
            yhistarr2(ixx,integral)=yhistarr2(ixx,integral)
     1         +yhistarr2(jxx,integrand)
     1         *(xhistarr(jxx+1,integrand)-xhistarr(jxx,integrand))
            errhistarr2(ixx,integral)=errhistarr2(ixx,integral)
     1         +(errhistarr2(jxx,integrand)
     1         *(xhistarr(jxx+1,integrand)-xhistarr(jxx,integrand)))**2
          enddo
          errhistarr2(ixx,integral)=sqrt(errhistarr2(ixx,integral))
        enddo
      else
         write(*,*) 'subroutine integratehist: error!'
         write(*,*) '--------------------------------'
         write(*,*) 'Invalid direction specified for histogram given.'
         write(*,*) 'direction=1/-1 only - integral from lowest bin to'
         write(*,*) 'highest bin or highest bin to lowest bin.'
         call exit(-1)
      endif

      end


      subroutine get_afb_hist(f_idx,b_idx,afb_idx)
      implicit none
      include 'pwhg_bookhist-new.h'
      real*8   f,ef,b,eb
      integer  ixx,f_idx,b_idx,afb_idx
      
      do ixx=1,nbins(f_idx)
         f=yhistarr2(ixx,f_idx)
         ef=errhistarr2(ixx,f_idx)
         b=yhistarr2(ixx,b_idx)
         eb=errhistarr2(ixx,b_idx)
         if((f+b).gt.0d0) then         ! Guard against division by zero.
            yhistarr2(ixx,afb_idx)=(f-b)/(f+b)
            errhistarr2(ixx,afb_idx)=2*sqrt((f*ef)**2+(b*eb)**2)
     1                              /(f+b)**2
         else
            yhistarr2(ixx,afb_idx)=0d0
            errhistarr2(ixx,afb_idx)=0d0
         endif
      enddo

      end


      function islept(j)
      implicit none
      logical islept
      integer j
      if(abs(j).ge.11.and.abs(j).le.16) then
         islept = .true.
      else
         islept = .false.
      endif
      end

      function phepDot(p_A,p_B)
      implicit none
      real * 8  phepDot
      real * 8  p_A(4),p_B(4)
      phepDot=p_A(4)*p_B(4)-p_A(1)*p_B(1)
     1       -p_A(2)*p_B(2)-p_A(3)*p_B(3)
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


      subroutine buildjets(mjets,kt,eta,rap,phi,pjet,jetvechep,
     1                                               isForClustering)
c     arrays to reconstruct jets
      implicit  none
      include  'hepevt.h'
      integer   maxtracks,maxjets
      parameter (maxtracks=nmxhep,maxjets=20)
      integer   mjets,jetvechep(maxtracks)
      real * 8  kt(maxjets),eta(maxjets),rap(maxjets),
     1     phi(maxjets),pjet(4,maxjets)
      logical   isForClustering(maxtracks)
      real * 8  ptrack(4,maxtracks),pj(4,maxjets)
      integer   jetvec(maxtracks),itrackhep(maxtracks)
      integer   ntracks,njets
      integer   j,k,mu
      real * 8  r,palg,ptmin,pp,tmp
      logical sonofid
      external sonofid
C - Initialize arrays and counters for output jets
      ptrack = 0
      jetvec = 0
      ntracks=0
      pjet = 0
      pj = 0
C - Extract final state particles to feed to jet finder
      do j=1,nhep
         if(.not.isForClustering(j)) cycle
         if(ntracks.eq.maxtracks) then
            write(*,*) 'analyze: need to increase maxtracks!'
            write(*,*) 'ntracks: ',ntracks
            call exit(-1)
         endif
         ntracks=ntracks+1
         ptrack(:,ntracks) = phep(1:4,j)
         itrackhep(ntracks)=j
      enddo
      if (ntracks.eq.0) then
         mjets=0
         return
      endif
C --------------- C
C - Run FastJet - C
C --------------- C
C - R = 0.7   radius parameter
C - f = 0.75  overlapping fraction
      palg  = -1
      r     = 0.5d0
      ptmin = 10d0
      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin,pjet,njets,jetvec)
      mjets=min(mjets,njets)
      if(njets.eq.0) return
c check consistency
      do k=1,ntracks
         if(jetvec(k).gt.0) then
            do mu=1,4
               pj(mu,jetvec(k))=pj(mu,jetvec(k))+ptrack(mu,k)
            enddo
         endif
      enddo
      tmp=0
      do j=1,mjets
         do mu=1,4
            tmp=tmp+abs(pj(mu,j)-pjet(mu,j))
         enddo
      enddo
      if(tmp.gt.1d-4) then
         write(*,*) ' bug!'
      endif
C --------------------------------------------------------------------- C
C - Computing arrays of useful kinematics quantities for hardest jets - C
C --------------------------------------------------------------------- C
      do j=1,mjets
         kt(j)=sqrt(pjet(1,j)**2+pjet(2,j)**2)
         pp = sqrt(kt(j)**2+pjet(3,j)**2)
         eta(j)=0.5d0*log((pp+pjet(3,j))/(pp-pjet(3,j)))
         rap(j)=0.5d0*log((pjet(4,j)+pjet(3,j))/(pjet(4,j)-pjet(3,j)))
         phi(j)=atan2(pjet(2,j),pjet(1,j))
      enddo
      jetvechep = 0
      do j=1,ntracks
         jetvechep(itrackhep(j))=jetvec(j)
      enddo
      end

      function sonofid(m,k)
      implicit none
      logical sonofid
      integer m,k
      include  'hepevt.h'
      integer j,kcurr
      integer ngenerations
      common/cngenerations/ngenerations
      kcurr=k
      do j=1,ngenerations
         if(abs(idhep(kcurr)).eq.m) then
            sonofid = .true.
            return
         endif
         kcurr = jmohep(1,kcurr)
         if(kcurr.eq.0) then
            sonofid = .false.
            return
         endif
      enddo
      sonofid=.false.
      end


      function sonofhep(m,k)
      implicit none
      logical sonofhep
      integer m,k
      include  'hepevt.h'
      integer j,kcurr
      integer ngenerations
      common/cngenerations/ngenerations
      kcurr=k
      do j=1,ngenerations
         if(kcurr.eq.m) then
            sonofhep = .true.
            return
         endif
         kcurr = jmohep(1,kcurr)
         if(kcurr.eq.0) then
            sonofhep = .false.
            return
         endif
      enddo
      sonofhep = .false.
      end

      subroutine boost2reson4(pres,nm,pin,pout)
      implicit none
      integer nm
      real * 8 pres(4),pin(4,nm),pout(4,nm)
      real * 8 vec(3),beta
      beta=sqrt(pres(1)**2+pres(2)**2+pres(3)**2)/pres(4)
      vec(1)=pres(1)/(beta*pres(4))
      vec(2)=pres(2)/(beta*pres(4))
      vec(3)=pres(3)/(beta*pres(4))
      call mboost4(nm,vec,-beta,pin,pout)
      end


      
      subroutine mboost4(m,vec,beta,vin,vout)
c     boosts the m vectors vin(4,m) into the vectors vout(4,m) (that can
c     be the same) in the direction of vec(3) (|vec|=1) with velocity
c     beta.  Lorents convention: (t,x,y,z).
      implicit none
      integer m
      real * 8 vec(3),beta,vin(4,m),vout(4,m)
      real * 8 betav,gamma
      real * 8 vdotb
      integer ipart,idim
      gamma=1/sqrt(1-beta**2)
      do ipart=1,m
         vdotb=vin(1,ipart)*vec(1)
     #         +vin(2,ipart)*vec(2)+vin(3,ipart)*vec(3)
         do idim=1,3
            vout(idim,ipart)=vin(idim,ipart)
     #           +vec(idim)*((gamma-1)*vdotb
     #           +gamma*beta*vin(4,ipart))
         enddo
         vout(4,ipart)=gamma*(vin(4,ipart)+vdotb*beta)
      enddo
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

      subroutine deltaplot(p1,p2,dsig,prefix,postfix)
      implicit none
      real * 8 p1(4),p2(4),dsig
      character *(*) prefix,postfix
      real * 8 dy,deta,delphi,dr
      call getdydetadphidr(p1,p2,dy,deta,delphi,dr)
      call filld(prefix//'-dy'//postfix,dy,dsig)
      call filld(prefix//'-deta'//postfix,deta,dsig)
      call filld(prefix//'-delphi'//postfix,delphi,dsig)
      call filld(prefix//'-dr'//postfix,dr,dsig)
      end

