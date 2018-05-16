      subroutine setborn(p,bflav,bres,bornjk,bmunu)
      implicit none
      include 'nlegborn.h'
c      include 'pwhg_math.h'
      include 'pwhg_st.h'
c      include 'pwhg_kn.h'
      include 'constants.f'
      include 'plabel.f'
      include 'process.f'
      include 'scale.f'
      include 'qcdcouple.f'
      include 'msq_cs.f'
      include 'mcfmtopwhg.f'
      integer nlegs
      parameter (nlegs=nlegborn)
c For t tbar with decay, nlegborn=12
      real * 8 p(0:3,nlegs),bres,bornjk(nlegs,nlegs),
     &     bmunu(0:3,0:3,nlegs)
      integer bflav(nlegs),ro
c mxpart is in constants.f
      real * 8 mcfmp(mxpart,4)
      double precision msq(-nf:nf,-nf:nf)
      integer iem,iep,inu,inub,ib,ibb,fl,j,k,myin(nlegs)
      logical isodd,isaquark
      external qqb_QQbdk_gvec
      real * 8 brcorr, brcorrect
      external brcorrect
c set scale in mcfm blocks
      scale=sqrt(st_muren2)
      musq=st_muren2
      as=st_alpha
      gsq=4d0*pi*as
      ason2pi=as/(2d0*pi)
c identify semileptonic decays
      do j=3,nlegs
c W decay products are non-b, non-t fermions, odd are down type (e or d,s)
c even are up type (nu or u,c)
         fl=bflav(j)
c first non-b products
         if(abs(fl).ne.5) then
            if(abs(fl).ne.6.and.abs(fl).lt.20) then
c now only non b t dec. products are allowed
               if(isodd(fl)) then
c e mu tau, d,s
                  if(fl.gt.0) then
                     iem=j
                  else
                     iep=j
                  endif
               else
c nu, u,c
                  if(fl.gt.0) then
                     inu=j
                  else
                     inub=j
                  endif
               endif
            endif
         else
c b
            if(fl.gt.0) then
               ib=j
            else
               ibb=j
            endif
         endif
      enddo
C     setup plabels for MCFM
      plabel(:)='ig'
      if (isaquark(bflav(inu))) then
      plabel(3)='pp'
      case='tt_bbh'
      endif
      if (isaquark(bflav(iem))) then
      plabel(7)='pp'
      case='tt_bbh'
      endif
c Pedantic e nu example
c if(isodd(bflav(j)).and.bflav(j).gt.0.and.bflav(j).lt.20) iem=j
c if(iseven(bflav(j)).and.bflav(j).gt.0.and.bflav(j).lt.20.and.bflav(j).ne.6) inu=j
c         if(bflav(j).eq.11) iem=j
c         if(bflav(j).eq.-12) inub=j
c         if(bflav(j).eq.-11) iep=j
c         if(bflav(j).eq.12) inu=j
c         if(bflav(j).eq.5) ib=j
c         if(bflav(j).eq.-5) ibb=j
c
c iem, inub, ibb come from tb      
c iep, inu, ib come from t
c incoming partons

C     vector for redirection of powheg vector onto mcfm
      myin(:)=-1
      myin(1)=1
      myin(2)=2
      myin(inu)=3
      myin(iep)=4
      myin(ib)=5
      myin(ibb)=6
      myin(iem)=7
      myin(inub)=8


      do j=1,nlegs
      do ro=1,4
      if (myin(j) .gt. 0) then
      if (j .le. 2) then
c incoming partons
      mcfmp(myin(j),ro)=-p(pwhg(ro),j)
      else
      mcfmp(myin(j),ro)=+p(pwhg(ro),j)
      endif
      endif
      enddo
      enddo


C-----Setup Born matrix element
      call qqb_QQbdk(mcfmp,msq)
      bres=msq(bflav(1),bflav(2))


      bornjk(:,:)=zip
      bmunu(:,:,:)=zip

      if ((bflav(1).eq.0).and.(bflav(2).eq.0)) then
C     gg-> QQb case
C-----Remember total cross section is 
C-----msq_cs(1,j,k)+msq_cs(2,j,k)+msq_cs(0,j,k)
C-----where
C-----msq_cs(1,j,k) propto A1^2
C-----msq_cs(2,j,k) propto A2^2
C-----msq_cs(0,j,k) propto -1/xn^2*|(A1+A2)|^2

         bornjk(1,2)=xn/2d0*(msq_cs(1,0,0)+msq_cs(2,0,0))
         bornjk(1,4)=xn/2d0*(msq_cs(1,0,0)+msq_cs(0,0,0))
         bornjk(1,3)=xn/2d0*(msq_cs(2,0,0)+msq_cs(0,0,0))
         bornjk(2,3)=bornjk(1,4)
         bornjk(2,4)=bornjk(1,3)
         bornjk(3,4)=-0.5d0/xn*(msq_cs(1,0,0)+msq_cs(2,0,0)
     &        +msq_cs(0,0,0)*(1d0+xn**2))
         
         
         call setupBmunu(mcfmp,myin,qqb_QQbdk_gvec,nlegs,bflav,Bmunu)
         
      elseif (bflav(2).eq.-bflav(1)) then
         bornjk(1,2)=-0.5d0/xn*bres
         if (bflav(1) .gt. 0d0) then
C--   qqb case
            bornjk(1,3)=(0.5d0*xn-1d0/xn)*bres
            bornjk(1,4)=+1d0/xn*bres
         else
C--   qbq case
            bornjk(1,4)=(0.5d0*xn-1d0/xn)*bres
            bornjk(1,3)=+1d0/xn*bres
         endif
         bornjk(2,3)=bornjk(1,4)
         bornjk(2,4)=bornjk(1,3)
         bornjk(3,4)=bornjk(1,2)
      endif

      bornjk(3,ib)=cf*bres
      bornjk(4,ibb)=cf*bres

c For hadronic W decays must add colour ordered amplitudes

      if(abs(bflav(iep)).lt.11) then
         bornjk(iep,inu)=cf*bres
         bornjk(inu,iep)=cf*bres
      endif

      if(abs(bflav(iem)).lt.11) then
         bornjk(iem,inub)=cf*bres
         bornjk(inub,iem)=cf*bres
      endif

C fill other non-zero values
      do j=1,nlegs-1
         do k=j+1,nlegs
            bornjk(k,j)=bornjk(j,k)
         enddo
      enddo
c Supply strong correction to branching ratio, if needed
      brcorr = brcorrect(p)
      bres   = bres   * brcorr
      bornjk = bornjk * brcorr
      bmunu  = bmunu  * brcorr

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
c      include 'nlegborn.h'
c      include 'pwhg_flst.h'
c      include 'pwhg_kn.h'
      integer iclabel
      common/ciclabel/iclabel
      real * 8 random
      external random
      iclabel=500
c--- Set up color labels for basic processes, gg -> tt~ and qq~ -> tt~      
      if(idup(1).eq.21) then
c gg
         if(random().gt.0.5d0) then
            call clinkqgga
     1           (icolup(1,3),icolup(1,1),icolup(1,2),icolup(1,4))
         else
            call clinkqgga
     1           (icolup(1,3),icolup(1,2),icolup(1,1),icolup(1,4))
         endif
      elseif(idup(1).gt.0) then
         call clinkqa(icolup(1,3),icolup(1,1))
         call clinkqa(icolup(1,2),icolup(1,4))
      else
         call clinkqa(icolup(1,3),icolup(1,2))
         call clinkqa(icolup(1,1),icolup(1,4))
      endif
c---  Copy color labels from t and t~ to b and b~
      icolup(:,11)=icolup(:,3)
      icolup(:,12)=icolup(:,4)
c--- Color-less W+ and W-
      icolup(:,5)=0
      icolup(:,6)=0
c--- W+ decay products
      if (abs(idup(7)) .le. 6) then
c---    jets
        if (idup(7) .gt. 0) then
           call clinkqa(icolup(1,7),icolup(1,8))
	else
           call clinkqa(icolup(1,8),icolup(1,7))
	endif
      else
c---    leptons
         icolup(:,7)=0
         icolup(:,8)=0
      endif            
c--- W- decay products
      if (abs(idup(9)) .le. 6) then
c---    jets
        if (idup(9) .gt. 0) then
           call clinkqa(icolup(1,9),icolup(1,10))
	else
           call clinkqa(icolup(1,10),icolup(1,9))
	endif
      else
c---    leptons
         icolup(:,9)=0
         icolup(:,10)=0
      endif            
           
c     1 and 2 are incoming! conjugate color
      call conjcolor(icolup(1,1))
      call conjcolor(icolup(1,2))
      end


      subroutine conjcolor(cl)
      integer cl(2),i
      i=cl(1)
      cl(1)=cl(2)
      cl(2)=i
      end

c Subroutine to link colours for
c quark - gluon -gluon  -gluon - aquark
c in planar order
      subroutine clinkqggga(ic1,ic2,ic3,ic4,ic5)
      integer iclabel
      common/ciclabel/iclabel
      integer ic1(2),ic2(2),ic3(2),ic4(2),ic5(2)
c ic1 is a quark: has colour, zero anticolor
      ic1(1)=iclabel+2
      ic1(2)=0
c ic2 is a gluon: link to quark
      ic2(1)=iclabel+3
      ic2(2)=iclabel+2
c ic3 is a gluon
      ic3(1)=iclabel+4
      ic3(2)=iclabel+3
c ic4 is an gluon
      ic4(1)=iclabel+5
      ic4(2)=iclabel+4
c ic5 is an antiquark
      ic5(1)=0
      ic5(2)=iclabel+5
      iclabel=iclabel+10
      end

c Subroutine to link colours for
c quark - gluon -gluon - aquark
c in planar order
      subroutine clinkqgga(ic1,ic2,ic3,ic4)
      integer iclabel
      common/ciclabel/iclabel
      integer ic1(2),ic2(2),ic3(2),ic4(2)
c ic1 is a quark: has colour, zero anticolor
      ic1(1)=iclabel+2
      ic1(2)=0
c ic2 is a gluon: link to quark
      ic2(1)=iclabel+3
      ic2(2)=iclabel+2
c ic3 is a gluon
      ic3(1)=iclabel+4
      ic3(2)=iclabel+3
c ic4 is an anti-quark
      ic4(1)=0
      ic4(2)=iclabel+4
      iclabel=iclabel+10
      end

c Subroutine to link colours for
c quark - gluon - aquark
c in planar order
      subroutine clinkqga(ic1,ic2,ic3)
      integer iclabel
      common/ciclabel/iclabel
      integer ic1(2),ic2(2),ic3(2)
c ic1 is a quark: has colour, zero anticolor
      ic1(1)=iclabel+2
      ic1(2)=0
c ic2 is a gluon: link to quark
      ic2(1)=iclabel+3
      ic2(2)=iclabel+2
c ic3 is an antiquark
      ic3(1)=0
      ic3(2)=iclabel+3
      iclabel=iclabel+10
      end

c Subroutine to link colours for
c quark - aquark
c in planar order
      subroutine clinkqa(ic1,ic2)
      integer iclabel
      common/ciclabel/iclabel
      integer ic1(2),ic2(2)
c ic1 is a quark: has colour, zero anticolor
      ic1(1)=iclabel+2
      ic1(2)=0
c ic2 is an antiquark
      ic2(1)=0
      ic2(2)=iclabel+2
      iclabel=iclabel+10
      end



      subroutine finalize_lh
c Thing to do: in case of hadronic decay of a W, set
c the right proportion of u dbar, u sbar, c sbar and c dbar
      implicit none
      include 'LesHouches.h'
      include 'PhysPars.h'
      integer i
      real * 8 random
      external random
      real * 8 sin2cb
      logical ini,nospincorr
      data ini/.true./
      save ini,nospincorr
      real * 8 powheginput
      external powheginput
      sin2cb = ph_CKM(1,2)**2
      do i=7,9,2
         if(abs(idup(i)).lt.3.and.abs(idup(i)).ge.1) then
            if(abs(idup(i)).eq.1) then
               if(random().gt.0.5d0) then
                  idup(i)=sign(3,idup(i))
                  if(random().gt.sin2cb) then
                     idup(i+1)=sign(4,idup(i+1))
                  else
                     idup(i+1)=sign(2,idup(i+1))
                  endif
               else
                  if(random().gt.sin2cb) then
                     idup(i+1)=sign(2,idup(i+1))
                  else
                     idup(i+1)=sign(4,idup(i+1))
                  endif
               endif
            elseif(abs(idup(i)).eq.2) then
               if(random().gt.0.5d0) then
                  idup(i)=sign(4,idup(i))
                  if(random().gt.sin2cb) then
                     idup(i+1)=sign(3,idup(i+1))
                  else
                     idup(i+1)=sign(1,idup(i+1))
                  endif
               else
                  if(random().gt.sin2cb) then
                     idup(i+1)=sign(1,idup(i+1))
                  else
                     idup(i+1)=sign(3,idup(i+1))
                  endif
               endif
            else
               write(*,*) ' finalize_lh: something wrong'
               call exit(-1)
            endif
         endif
      enddo
      if(ini) then
         nospincorr = powheginput("#nospincorr").eq.1
      endif
      if(nospincorr) then
c perform a random rotation of the t (tbar) decay products in the
c t (tbar) rest frame
         if(ini) then
            write(*,*) ' rotating randomly the t and tbar systems'
            ini = .false.
         endif
         call randomrotate(3)
         call randomrotate(4)
      endif
      end

      subroutine randomrotate(ind)
      implicit none
      integer ind
      include 'LesHouches.h'
      real * 8 pres(5),vec(3),beta,r(3,3)
      logical sonof
      integer j
      beta=sqrt(pup(1,ind)**2+pup(2,ind)**2+pup(3,ind)**2)/pup(4,ind)
      vec(1)=pup(1,ind)/(beta*pup(4,ind))
      vec(2)=pup(2,ind)/(beta*pup(4,ind))
      vec(3)=pup(3,ind)/(beta*pup(4,ind))
      call uniformrot(r)
      do j=3,nup
         if(sonof(ind,j)) then
            call mboost5(1,vec,-beta,pup(:,j),pup(:,j))
            call matrixmultvec(r,pup(1:3,j))
            call mboost5(1,vec,beta,pup(:,j),pup(:,j))
         endif
      enddo
      end

      function sonof(m,k)
      implicit none
      logical sonof
      integer m,k
      include  'LesHouches.h'
      integer j,kcurr
      integer ngenerations
      parameter (ngenerations=4)
      kcurr=mothup(1,k)
      do j=1,ngenerations
         if(kcurr.eq.m) then
            sonof = .true.
            return
         endif
         kcurr = mothup(1,kcurr)
         if(kcurr.eq.0) then
            sonof = .false.
            return
         endif
      enddo
      sonof=.false.
      end


      subroutine mboost5(m,vec,beta,vin,vout)
c     boosts the m vectors vin(4,m) into the vectors vout(4,m) (that can
c     be the same) in the direction of vec(3) (|vec|=1) with velocity
c     beta.  Lorents convention: (t,x,y,z).
      implicit none
      integer m
      real * 8 vec(3),beta,vin(5,m),vout(5,m)
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
         vout(5,ipart)=vin(5,ipart)
      enddo
      end


      subroutine matrixmultvec(r,v)
      implicit none
      real * 8 r(3,3),v(3),res(3)
      integer j
      do j=1,3
         res(j)=r(j,1)*v(1)+r(j,2)*v(2)+r(j,3)*v(3)
      enddo
      v = res
      end


      subroutine uniformrot(R)
      implicit none
c     Generate a uniformly distributed rotation
      real * 8 r(3,3)
      real * 8 pi
      parameter (pi=3.141592653589793d0)
      real * 8 costh,sinth,phi,gamma,sing,cosg,norm
      real * 8 random
      external random
      costh=2*random()-1
      sinth=sqrt(abs(1-costh**2))
      phi=2*pi*random()
c First axis in random direction
      r(1,1)=costh
      r(2,1)=sinth*sin(phi)
      r(3,1)=sinth*cos(phi)
c now pick a vector orthogonal to the first axis
      if(costh.gt.0.5d0) then
         norm=sqrt(r(1,1)**2+r(2,1)**2)
         r(1,2)=r(2,1)/norm
         r(2,2)=-r(1,1)/norm
         r(3,2)=0
      else
         norm=sqrt(r(2,1)**2+r(3,1)**2)
         r(1,2)=0
         r(2,2)=r(3,1)/norm
         r(3,2)=-r(2,1)/norm
      endif
c Now totate r(:,2) around r(:,1) of an arbitrary angle
      gamma = 2*pi * random()
      sing = sin(gamma)
      cosg = cos(gamma)
      call mrotate(r(:,1),sing,cosg,r(:,2))
c Last axis is cross product of 1 and 2
      r(1,3)=r(2,1)*r(3,2)-r(3,1)*r(2,2)
      r(2,3)=r(3,1)*r(1,2)-r(1,1)*r(3,2)
      r(3,3)=r(1,1)*r(2,2)-r(2,1)*r(1,2)
      end

