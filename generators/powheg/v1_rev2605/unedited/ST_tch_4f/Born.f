      subroutine setborn(p,bflav,born,bornjk,bmunu)
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'
      double precision CKM(1:6,1:6)
      common/cckm/CKM
      double precision p(0:3,nlegborn),bornjk(nlegborn,nlegborn)
      integer bflav(nlegborn)
      double precision kb_mad(0:3,nlegborn),amp2mad
      double precision bmunu(0:3,0:3,nlegborn),born,colcf
      integer ileg,j,k,mu,nu
      double precision born_bmunu
      double precision ckm_light,ckm_tb
      logical ini
      data ini/.true./
      save ini
      include 'coupl.inc'

cccccccccccccccccccccccccccccccc
c     check with FC
c$$$      double precision k3k5,k1k4,k2k4,k1k2,mt,mb,mw,PropWW,PropTT,PropBB
c$$$     $     ,PropBT,bornFC,k1k3,k1k5,k2k3,k2k5,k3k4,k4k5
c$$$      double precision dotp
c$$$      external dotp
cccccccccccccccccccccccccccccccc



      if(ini) then 
c     setting madgraph parameters (needed for madgraph subroutines)
         call set_madgraph_parameters
         ini=.false.
         write(*,*) 'Relevant couplings used by Madgraph'
         write(*,*) 'e/(sinw*sqrt(2))     = ',gwf(1)
         write(*,*) 'MCFM gwsq (from Mad) = ',(gwf(1)*sqrt(2.))**2
c         stop
      endif
c     set madgraph parameters that can change on an event-by-event basis
      call mad_setparam
      do ileg=1,nlegborn
         do mu=0,3
            kb_mad(mu,ileg)=p(mu,ileg)
         enddo
      enddo
c     to avoid bugs in HELAS, restore exact masslessness of  incoming partons 
      kb_mad(0,1)=dabs(kb_mad(3,1))
      kb_mad(0,2)=dabs(kb_mad(3,2))

      if(ttype.eq.1) then
         call compborn(kb_mad,bflav,amp2mad)
      elseif(ttype.eq.-1) then
         call compborn_tb(kb_mad,bflav,amp2mad)
      else
         write(*,*) 'wrong ttype in Born.f'
         call exit(-1)
      endif
      born=amp2mad




c     Colour factors for colour-correlated Born amplitudes;
      do j=1,nlegborn
         if(abs(bflav(j)).le.6) then
            do k=j+1,nlegborn
               if(abs(bflav(k)).le.6) then
c     light current
                  if(abs(bflav(j)).le.4.and.abs(bflav(k)).le.4.and.
     $                 bflav(j)*bflav(k).ne.0) then
                     colcf=   cf
c     heavy antenna, t-b
                  elseif((abs(bflav(j)).eq.5.and.abs(bflav(k)).eq.6).or.
     $                    (abs(bflav(j)).eq.6.and.abs(bflav(k)).eq.5)) then
                     colcf=   - (ca-2.*cf)/2.
c     heavy antenna, t-g or b-g
                  elseif((abs(bflav(j)).eq.0.and.abs(bflav(k)).ge.5).or.
     $                    (abs(bflav(k)).eq.0.and.abs(bflav(j)).ge.5)) then
                     colcf=   ca/2.
                  else
c     non correlated
                     colcf=0.
                  endif
                  bornjk(j,k)=born*colcf
                  bornjk(k,j)=bornjk(j,k)
               endif
            enddo
         endif
      enddo
      
      ckm_tb=CKM(abs(bflav(3)),abs(bflav(4)))**2
      
c     !ER: with ttype as follows, it works.
      if(bflav(2).eq.0) then
c     since Bmunu is computed from scratch, I have to supply the
c     proper CKM factor.
         ckm_light=CKM(abs(bflav(1)),abs(bflav(5)))**2

         if(ttype*bflav(1).gt.0) then
            do mu=0,3
               do nu=0,3
                  bmunu(mu,nu,2)= born_bmunu(mu,nu,
     $                 kb_mad(0,1),kb_mad(0,2),
     $                 kb_mad(0,3),kb_mad(0,4),kb_mad(0,5))
     $                 /3./8./4.
     $                 *(st_alpha*4.*pi)
     $                 *(nc*tf*ca)
     $                 *(ph_alphaem*4.*pi)**2 /64./ph_sthw2**2
     $                 *8./3.
     $                 *ckm_light *ckm_tb
               enddo
            enddo
         elseif(ttype*bflav(1).lt.0) then
c     flip 1-5, since light current is flipped
            do mu=0,3
               do nu=0,3
                  bmunu(mu,nu,2)= born_bmunu(mu,nu,
     $                 kb_mad(0,5),kb_mad(0,2),
     $                 kb_mad(0,3),kb_mad(0,4),kb_mad(0,1))
     $                 /3./8./4.
     $                 *(st_alpha*4.*pi)
     $                 *(nc*tf*ca)
     $                 *(ph_alphaem*4.*pi)**2 /64./ph_sthw2**2
     $                 *8./3.
     $                 *ckm_light *ckm_tb
               enddo
            enddo
         else
            write(*,*) 'Problem in Born.f'
            call exit(1)
         endif
c$$$      write(*,*) '--------------------------'
c$$$      do mu=0,3
c$$$         write(*,*) (bmunu(mu,nu,2), nu=0,3)
c$$$      enddo
      elseif(bflav(1).eq.0) then
         ckm_light=CKM(abs(bflav(2)),abs(bflav(5)))**2

c     EXCHANGE 1 and 2, since bmunu assumes GLUON IN
c     POSITION 2
         if(ttype*bflav(2).gt.0) then
            do mu=0,3
               do nu=0,3
                  bmunu(mu,nu,1)= born_bmunu(mu,nu,
     $                 kb_mad(0,2),kb_mad(0,1),
     $                 kb_mad(0,3),kb_mad(0,4),kb_mad(0,5))
     $                 /3./8./4.
     $                 *(st_alpha*4.*pi)
     $                 *(nc*tf*ca)
     $                 *(ph_alphaem*4.*pi)**2 /64./ph_sthw2**2
     $                 *8./3.
     $                 *ckm_light *ckm_tb
               enddo
            enddo
         elseif(ttype*bflav(2).lt.0) then
c     flip 2-5, since light current is flipped
            do mu=0,3
               do nu=0,3
                  bmunu(mu,nu,1)= born_bmunu(mu,nu,
     $                 kb_mad(0,5),kb_mad(0,1),
     $                 kb_mad(0,3),kb_mad(0,4),kb_mad(0,2))
     $                 /3./8./4.
     $                 *(st_alpha*4.*pi)
     $                 *(nc*tf*ca)
     $                 *(ph_alphaem*4.*pi)**2 /64./ph_sthw2**2
     $                 *8./3.
     $                 *ckm_light *ckm_tb
               enddo
            enddo
         else
            write(*,*) 'Problem in Born.f'
            call exit(1)
         endif
      endif

      

cccccccccccccccccccccccccccccccccc
c$$$c     Born with Feyncalc, contracting Bmunu with gmunu
c$$$      if(bflav(2).eq.0.and.bflav(1).gt.0) then
c$$$         k3k5=dotp(p(0,3),p(0,5))
c$$$         k1k4=dotp(p(0,1),p(0,4))
c$$$         k2k4=dotp(p(0,2),p(0,4))
c$$$         k1k2=dotp(p(0,1),p(0,2))
c$$$
c$$$         k1k3=dotp(p(0,1),p(0,3))
c$$$         k1k5=dotp(p(0,1),p(0,5))
c$$$         k2k3=dotp(p(0,2),p(0,3))
c$$$
c$$$         k2k5=dotp(p(0,2),p(0,5))
c$$$         k3k4=dotp(p(0,3),p(0,4))
c$$$         k4k5=dotp(p(0,4),p(0,5))
c$$$
c$$$         mt=topmass_pow
c$$$         mb=bmass_pow
c$$$         mw=ph_Wmass
c$$$         PropWW=1./(-2.*k1k5-mw**2)**2
c$$$         PropTT=1./(-2.*k2k3)**2
c$$$         PropBB=1./(-2.*k2k4)**2
c$$$         PropBT=1./(4.*k2k3*k2k4)
c$$$
c$$$         bornFC=(-1024*k3k5*(-(k1k4*mb**2) + k1k2*(k2k4 + mb**2))*PropBB + 
c$$$     -    1024*(-(k1k3*k2k4*k3k5) + k1k2*k3k4*k3k5 + 
c$$$     -       k1k4*(k2k5*k3k4 + k2k3*k3k5 + k2k4*k3k5 - 2*k3k4*k3k5 - 
c$$$     -          k2k3*k4k5))*PropBT - 
c$$$     -    1024*k1k4*(k2k3*k2k5 + (k2k5 - k3k5)*mt**2)*PropTT)*PropWW
c$$$
c$$$         bornFC=-bornFC
c$$$     $        /3./8./4.
c$$$     $        *(st_alpha*4.*pi)
c$$$     $        *(nc*tf*ca)
c$$$     $        *(ph_alphaem*4.*pi)**2 /64./ph_sthw2**2
c$$$     $        *8./3.
c$$$
c$$$c$$$         print*, mt,mb,mw
c$$$c$$$         stop
c$$$
c$$$         print*, bflav,"--> ",born/bornFC
c$$$
c$$$      endif
cccccccccccccccccccccccccccccccccc



      end





      subroutine borncolour_lh
c Sets up the colour for the given flavour configuration
c already filled in the Les Houches interface.
c In case there are several colour structure, one
c should pick one with a probability proportional to
c the value of the corresponding cross section, for the
c kinematics defined in the Les Houches interface
      implicit none
      include 'LesHouches.h'
      include 'PhysPars.h'
      integer ileg,tmp
      integer tcol,bcol,lcol
      data tcol/501/
      data bcol/502/
      data lcol/510/


cccccccccccccccccccccccccccccc
c     t-channel subprocesses (4f)
c     3rd is the TOP, 4th the BOTTOM
cccccccccccccccccccccccccccccc

      if(abs(idup(3)).ne.6) then
         write(*,*)'error in borncolour_lh, top leg'
      endif
      icolup(1,3)=tcol
      icolup(2,3)=0


      if(abs(idup(4)).ne.5) then
         write(*,*)'error in borncolour_lh, b leg'
      endif
      icolup(1,4)=0
      icolup(2,4)=bcol

c     gu / gd
      if((idup(1).eq.21).and.
     $     (ttype*idup(2).gt.0).and.
     $     (ttype*idup(2).le.4)) then
c               ===>=== 3,t
c              /
c     1,g ggggg
c              \---<--- 4,bbar
c
c     2,q ------>------ 5,qp
         icolup(1,1)=tcol
         icolup(2,1)=bcol
         icolup(1,2)=lcol
         icolup(2,2)=0
         icolup(1,5)=icolup(1,2)
         icolup(2,5)=icolup(2,2)
c     gubar / gdbar
      elseif((idup(1).eq.21).and.
     $     (ttype*idup(2).lt.0).and.
     $     (ttype*idup(2).ge.-4)) then
c               ===>=== 3,t
c              /
c     1,g ggggg
c              \---<--- 4,bbar
c
c     2,qbar  ------<------ 5,qbarp
         icolup(1,1)=tcol
         icolup(2,1)=bcol
         icolup(1,2)=0
         icolup(2,2)=lcol
         icolup(1,5)=icolup(1,2)
         icolup(2,5)=icolup(2,2)
c     ug /dg
      elseif((idup(2).eq.21).and.
     $     (ttype*idup(1).gt.0).and.
     $     (ttype*idup(1).le.4)) then
c               ===>=== 3,t
c              /
c     2,g ggggg
c              \---<--- 4,bbar
c
c     1,q ------>------ 5,qp
         icolup(1,1)=lcol
         icolup(2,1)=0
         icolup(1,2)=tcol
         icolup(2,2)=bcol
         icolup(1,5)=icolup(1,1)
         icolup(2,5)=icolup(2,1)
c     ubarg / dbarg
      elseif((idup(2).eq.21).and.
     $     (ttype*idup(1).lt.0).and.
     $     (ttype*idup(1).ge.-4)) then
c               ===>=== 3,t
c              /
c     2,g ggggg
c              \---<--- 4,bbar
c
c     1,qbar  ------<------ 5,qbarp
         icolup(1,1)=0
         icolup(2,1)=lcol
         icolup(1,2)=tcol
         icolup(2,2)=bcol
         icolup(1,5)=icolup(1,1)
         icolup(2,5)=icolup(2,1)
      else
         write(*,*)'error in borncolor_lh, invalid flavs ',ttype
         write(*,*)'program stops'
         call exit(-1)
      endif



c     if charge-conjugated process, exchange color
c     with anticolors
      if(ttype.lt.0) then
         do ileg=1,5
            tmp=icolup(1,ileg)
            icolup(1,ileg)=icolup(2,ileg)
            icolup(2,ileg)=tmp
         enddo
      endif

      end

      subroutine finalize_lh
      implicit none
c     Set up the resonance whose mass must be preserved
c     on the Les Houches interface.
c     resonance Z -> e-(3) e+(4)
c$$$      call add_resonance(23,3,4)
c     give masses to final-state products
c$$$      call lhefinitemasses
      end
