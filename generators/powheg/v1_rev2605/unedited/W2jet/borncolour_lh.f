      subroutine borncolour_lh
c Sets up the colour for the given flavour configuration
c already filled in the Les Houches interface.
c In case there are several colour structure, one
c should pick one with a probability proportional to
c the value of the corresponding cross section, for the
c kinematics defined in the Les Houches interface.
c Here we assume all particles to be outgoing, and
c assign colour according to the corresponding colour amplitudes.
c At the end, the colour of incoming partons are conjugated.
      implicit none
      include 'nlegborn.h'
      include 'LesHouches.h'
      include 'pwhg_kn.h'
      include 'constants.f'
      include 'qq_cs.f'
      integer iq,ia,iq1,iq2,ia1,ia2,ig1,ig2,j,cpoint(4)
      integer bflav(nlegborn)
      double precision mcfmp(mxpart,4),bres
      data cpoint/1,2,5,6/
      
c     transform to MCFM momenta
      
      mcfmp(1,1:3)=-kn_cmpborn(1:3,1)
      mcfmp(1,4)=-kn_cmpborn(0,1)
      
      mcfmp(2,1:3)=-kn_cmpborn(1:3,2)
      mcfmp(2,4)=-kn_cmpborn(0,2)
      
      mcfmp(3,1:3)=kn_cmpborn(1:3,3)
      mcfmp(3,4)=kn_cmpborn(0,3)
      
      mcfmp(4,1:3)=kn_cmpborn(1:3,4)
      mcfmp(4,4)=kn_cmpborn(0,4)
      
      mcfmp(5,1:3)=kn_cmpborn(1:3,5)
      mcfmp(5,4)=kn_cmpborn(0,5)
      
      mcfmp(6,1:3)=kn_cmpborn(1:3,6)
      mcfmp(6,4)=kn_cmpborn(0,6)
      
      bflav(1)=idup(1)
      bflav(2)=idup(2)
      bflav(5)=idup(5)
      bflav(6)=idup(6)
      bflav(3)=12               ! not used 
      bflav(4)=-11              ! not used 
      do j=1,4
         if (bflav(cpoint(j)) .eq. 21) bflav(cpoint(j))=0
      enddo
!      write(*,*) 'bflav', bflav 
      call qqb_w2jet_pwhg(mcfmp,bflav,bres)
      
      
c     q qb g g or permutations-crossing
      if(idup(1).eq.21.or.idup(2).eq.21
     1     .or.idup(5).eq.21.or.idup(6).eq.21) then
c     find the quarks and gluons
         ig1=-1
         do j=1,4
            if(idup(cpoint(j)).eq.21) then
               if(ig1.lt.0) then
                  ig1=cpoint(j)
               else
                  ig2=cpoint(j)
               endif
            elseif(idup(cpoint(j))*istup(cpoint(j)).gt.0) then
               iq=cpoint(j)
            elseif(idup(cpoint(j))*istup(cpoint(j)).lt.0) then
               ia=cpoint(j)
            else
               write(*,*) 'borncolour_lh: should not be here!'
               call exit(1)
            endif
         enddo
         
         call borncolourqagg(icolup(1,iq),icolup(1,ia),
     1        icolup(1,ig1),icolup(1,ig2))
         
      else
c q q qb qb, or q Q qb Qb, plus permutations-crossing
         iq1=-1
         iq2=-1
         ia1=-1
         ia2=-1
         do j=1,4
            if(idup(cpoint(j))*istup(cpoint(j)).gt.0) then
               if(iq1.lt.0) then
                  iq1=cpoint(j)
               else
                  iq2=cpoint(j)
               endif
            else
               if(ia1.lt.0) then
                  ia1=cpoint(j)
               else
                  ia2=cpoint(j)
               endif
            endif
         enddo
         call borncolour4q(icolup(1,iq1),icolup(1,iq2),
     1        icolup(1,ia1),icolup(1,ia2))
      endif
c     Conjugate incoming colours
      call colour_conj(icolup(1,1))
      call colour_conj(icolup(1,2))   
      end
      
      subroutine borncolour4q(icol1,icol2,icol3,icol4)
c                             q     q     qbar  qbar
      implicit none
      include 'qq_cs.f'
      integer icol1(2),icol2(2),icol3(2),icol4(2)
      real * 8 r
      real * 8 random
c q Q->q Q
      r=random()*(qq_cs(1)+qq_cs(2))
      if(r.lt.qq_cs(1)) then
c t channel gluon (23 channel; thus colour is exchanged
c 1->4 and 2->3
         call colourjoin4q(icol1,icol2,icol3,icol4)
      else
c u channel gluon (23 channel; thus colour is exchanged
c 1->3 and 2->4
         call colourjoin4q(icol1,icol2,icol4,icol3)
      endif
      end

      subroutine colourjoin4q(icol1,icol2,icol3,icol4)
c                             q     q     qbar  qbar
c perform a planar colour connection on the planar sequence
c q q qbar qbar
      implicit none
      integer icol1(2),icol2(2),icol3(2),icol4(2)
      integer newcolor
      icol1(2)=0
      icol2(2)=0
      icol3(1)=0
      icol4(1)=0
      call getnewcolor(newcolor)
      icol1(1)=newcolor
      icol4(2)=newcolor
      call getnewcolor(newcolor)
      icol2(1)=newcolor
      icol3(2)=newcolor
      end

