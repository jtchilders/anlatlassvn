c      program testfindregions
c      implicit none
c      integer n
c      parameter (n=6)
c      integer rflav(n),nregions,iregions(2,n*(n-1)/2),i
c      character  * 30 process
c
c     process = 'u~ u  -> e- e+ u  u~'
c     call from_madgraph_to_number(process,rflav) 
c     call find_regions(n,rflav,nregions,iregions)
c      
c     write(*,*) nregions
c     do i=1,nregions
c        write(*,*) iregions(1,i),iregions(2,i)
c     enddo
c      call genflavreglist
c      end



      subroutine find_regions(a,ares,atags,indexreal,nregions,iregions)
c Finds all singular regions for the real graph indexed by indexreal.
c It returns:
c integer nregions
c integer iregion(2,nregions): the indices of particles forming singular
c                              regions i,j (i<j).
c                              For initial state singularities, if the
c                              emitter can be both of the initial state 
c                              particles, 
c                              and if the radiated particle is a gluon,
c                              only one region is generated with first
c                              index equal to zero.
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer a(nlegreal,*),ares(nlegreal,*),atags(nlegreal,*)
      integer indexreal,nregions,iregions(2,maxregions)
      logical ireg(2)
      logical validBorn
      external validBorn
      integer itag,iret,i,j,k,ibornfl,iborn,flrad,
     1     bflav(nlegborn),res(nlegborn),tags(nlegborn)
      nregions=0
c final state regions
      do i=flst_lightpart,nlegreal
         do j=i+1,nlegreal
c find if they can arise from the same splitting
            call same_splitting(a,ares,atags,
     1           indexreal,i,j,ibornfl,itag,iret)
c cannot come from the same splitting
            if(iret.lt.0) then
               goto 10
            endif
c build the underlying born flavour structure in bflav
            iborn=0
            do k=1,nlegreal
               if(k.eq.i) then
                  iborn=iborn+1
                  bflav(iborn)=ibornfl
                  res(iborn)=ares(k,indexreal)
                  tags(iborn)=itag
               elseif(k.ne.j) then
                  iborn=iborn+1
                  bflav(iborn)=a(k,indexreal)
                  res(iborn)=ares(k,indexreal)
                  tags(iborn)=atags(k,indexreal)
               endif
            enddo
            if(validBorn(bflav,tags,res)) then
               nregions=nregions+1
               iregions(1,nregions)=i
               iregions(2,nregions)=j
            endif
 10         continue
         enddo
      enddo
c initial state region
      do j=flst_lightpart,nlegreal
         do i=1,2
            ireg(i)=.false.
            call same_splitting(a,ares,atags,
     1           indexreal,i,j,ibornfl,itag,iret)
            if(iret.lt.0) then
               goto 11
            endif
            iborn=0
            do k=1,nlegreal
               if(k.eq.i) then
                  iborn=iborn+1
                  bflav(iborn)=ibornfl
                  res(iborn)=0
                  tags(iborn)=itag
               elseif(k.ne.j) then
                  iborn=iborn+1
                  bflav(iborn)=a(k,indexreal)
                  res(iborn)=ares(k,indexreal)
                  tags(iborn)=atags(k,indexreal)
               endif
            enddo
            if(validBorn(bflav,tags,res)) then
               ireg(i)=.true.
            endif
 11         continue
         enddo
         flrad=a(nlegreal,indexreal)
         if(ireg(1).and.ireg(2).and.(flrad.eq.0.or.flrad.eq.22))
     1        then
c if both regions are singular and the radiated parton is a gluon
c or a photon emit a single region with emitter 0
            nregions=nregions+1
            iregions(1,nregions)=0
            iregions(2,nregions)=j
         else
            if(ireg(1)) then
               nregions=nregions+1
               iregions(1,nregions)=1
               iregions(2,nregions)=j
            endif
            if(ireg(2)) then
               nregions=nregions+1
               iregions(1,nregions)=2
               iregions(2,nregions)=j
            endif
         endif
      enddo
      end

      subroutine ubornflav(alr)
      implicit none
      integer alr
      include 'nlegborn.h'
      include 'pwhg_flst.h'
c finds the underlying Born flavour
c integer n: number of legs in real graph
c integer rflav(n): flavours of legs in real graph
c                   (1 and 2 incoming)
      integer itag,ibornfl,iret,l,n,j
c j:  singularity in region j,n
c     j=0 (1 and 2), j=1, j=2: initial state sing.
c     j>2 final state sing.
      j=flst_emitter(alr)
c if it is both 1 and 2, pretend it is 1
      if(j.eq.0) j=1
      iret=0
      n=nlegreal
c this is only to find itag and ibornfl. We already now that j and n
c come from the same splitting.
      if(j.eq.1.or.j.eq.2) then
         call same_splitting0('isr',flst_alr(j,alr),flst_alr(n,alr),
     1       flst_alrtags(j,alr),flst_alrtags(n,alr),itag,ibornfl,iret)
      elseif(j.gt.2) then
         call same_splitting0('fsr',flst_alr(j,alr),flst_alr(n,alr),
     1       flst_alrtags(j,alr),flst_alrtags(n,alr),itag,ibornfl,iret)
      endif
      if(iret.lt.0) then
         write(*,*) ' ubornflav: error'
         write(*,*) ' j: ',j
         call print_lists(nlegreal,flst_alr(l,alr),
     1              flst_alrres(l,alr),flst_alrtags(l,alr))
         write(*,*) ' emitter:',flst_emitter(alr)
         call exit(-1)
      endif
      do l=1,nlegborn
         if(l.eq.j) then
            flst_uborn(l,alr)=ibornfl
            flst_uborntags(l,alr)=itag
            flst_ubornres(l,alr)=flst_alrres(l,alr)
         else
            flst_uborn(l,alr)=flst_alr(l,alr)
            flst_uborntags(l,alr)=flst_alrtags(l,alr)
            flst_ubornres(l,alr)=flst_alrres(l,alr)
         endif
      enddo
      return
 998  continue
      end


      function flavequiv(m,n,ja,jb,arr1,arr2,arr3)
c returns true if the flavour structures aflav and bflav are
c equivalent up to a permutation of the final state lines,
c false otherwise.
      implicit none
      logical flavequiv
      integer m,n,ja,jb,arr1(m,*),arr2(m,*),arr3(m,*)
c we need the parameter nlegreal
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer j,k,itmp,i1(nlegreal),i2(nlegreal),i3(nlegreal)
      call intassign(n,arr1(1,jb),i1)
      call intassign(n,arr2(1,jb),i2)
      call intassign(n,arr3(1,jb),i3)
      do j=1,n
         if(
     1         arr1(j,ja).ne.i1(j)
     2     .or.arr2(j,ja).ne.i2(j)
     2     .or.arr3(j,ja).ne.i3(j) ) then
            if(j.le.2) then
               flavequiv=.false.
               return
            endif
            do k=j+1,n
               if(
     1          arr1(j,ja).eq.i1(k)
     2     .and.arr2(j,ja).eq.i2(k)
     2     .and.arr3(j,ja).eq.i3(k) ) then
                  call exchange_ind(j,k,i1,i2,i3)
                  goto 10
               endif
            enddo
            flavequiv=.false.
            return
         endif
 10      continue
      enddo
      flavequiv=.true.
      end

      function validBorn(bflav,tags,res)
c Find if the flavour structure bflav is equivalent to an element
c in the list of Born processes. Equivalence means that it can be
c made identical with a permutation of final state particles.
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer bflav(nlegborn),tags(nlegborn),res(nlegborn)
      logical validBorn
      integer j,k,kb,itmp,ib(nlegborn),it(nlegborn),ir(nlegborn)
      call intassign(nlegborn,bflav,ib)
      call intassign(nlegborn,tags,it)
      call intassign(nlegborn,res,ir)
      do kb=1,flst_nborn
         do j=1,nlegborn
c recursive exit
            if(j.eq.nlegborn
     1           .and.flst_born(j,kb).eq.ib(j)
     2           .and.flst_borntags(j,kb).eq.it(j)
     3           .and.flst_bornres(j,kb).eq.ir(j) ) then
               validBorn=.true.
               return
            endif
            if(
     1                flst_born(j,kb).ne.ib(j)
     2           .or.flst_borntags(j,kb).ne.it(j)
     3           .or.flst_bornres(j,kb).ne.ir(j) ) then
               if(j.le.2) then
                  goto 999
               endif
               do k=j+1,nlegborn
                  if(
     1                flst_born(k,kb).eq.ib(j)
     2           .and.flst_borntags(k,kb).eq.it(j)
     3           .and.flst_bornres(k,kb).eq.ir(j) ) then
                     call exchange_ind(j,k,ib,ir,it)
                     goto 10
                  endif
               enddo
               goto 999
            endif
 10         continue
         enddo         
 999     continue
      enddo
c      write(*,*) ' validBorn false,',bflav
      validBorn=.false.
      end

c      -6  -5  -4  -3  -2  -1  0  1  2  3  4  5  6                    
c      t~  b~  c~  s~  u~  d~  g  d  u  s  c  b  t                    
      subroutine from_madgraph_to_number(stringa,ferm_flav)
      implicit none
      integer nmax
      parameter (nmax=30)
      character stringa(nmax)
      integer ferm_flav(*)

      integer i, parton
      character *2 flav(-5:5)
      real * 8 charge(-5:5)
      common/flav_ordering/charge,flav

      parton = 0
      do i=1,nmax
         if (stringa(i).eq.'g') then
            parton = parton + 1
            ferm_flav(parton) = 0
         elseif (stringa(i).eq.'H') then
            parton = parton + 1
            ferm_flav(parton) = 503
         elseif (stringa(i).eq.'d') then
            parton = parton + 1
            ferm_flav(parton) = +1
         elseif (stringa(i).eq.'u') then
            parton = parton + 1
            ferm_flav(parton) = +2
         elseif (stringa(i).eq.'s') then
            parton = parton + 1
            ferm_flav(parton) = +3
         elseif (stringa(i).eq.'c') then
            parton = parton + 1
            ferm_flav(parton) = +4
         elseif (stringa(i).eq.'b') then
            parton = parton + 1
            ferm_flav(parton) = +5
         elseif (stringa(i).eq.'t') then
            parton = parton + 1
            ferm_flav(parton) = +6
         elseif (stringa(i).eq.'~') then
            ferm_flav(parton) = -ferm_flav(parton)
         elseif (stringa(i).eq.' ') then
         elseif (stringa(i).eq.'Z') then
            parton = parton + 1
            ferm_flav(parton) = +10
            parton = parton + 1
            ferm_flav(parton) = -10
          elseif (stringa(i).eq.'e') then
            parton = parton + 1
            ferm_flav(parton) = +10
          elseif (stringa(i).eq.'+') then
            ferm_flav(parton) = -ferm_flav(parton)           
          elseif (stringa(i).eq.'/') then
             return
        endif            
      enddo

      end
      function isalightparton(ipart)
      implicit none
      include 'pwhg_st.h'
      logical isalightparton
      integer ipart
      if(abs(ipart).le.st_nlight) then
         isalightparton=.true.
         return
      endif
      if(ipart.eq.22) then
         isalightparton=.true.
         return
      endif
      if(abs(ipart).ge.11.and.abs(ipart).le.16) then
         isalightparton=.true.
         return
      endif
      isalightparton=.false.
      end

      subroutine genflavreglist
      implicit none
      include 'pwhg_flg.h'
      include 'pwhg_st.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer nregions,iregions(2,maxregions)
      integer iflregl,k,l,ipart,j,itmp,nreg,iret,tmpfl,fl1,fl2
      logical equalintlists
      external equalintlists
      logical verbose
      parameter (verbose=.true.)
      logical flavequiv,isalightparton,equiv_entry_alr_real,equal_lists
      external flavequiv,isalightparton,equiv_entry_alr_real,equal_lists
c check that there are no coloured light partons before flst_lightpart
      do j=1,flst_nreal
         do ipart=3,flst_lightpart -1
            if(isalightparton(flst_real(ipart,j))) then
               write(*,*) 
     1      ' genflavreglist: light parton before flst_lightpart'
c               stop
            endif
         enddo
         do ipart=flst_lightpart,nlegreal
            if(.not.isalightparton(flst_real(ipart,j))) then
               write(*,*) 
     1      ' genflavreglist: not a light parton after flst_lightpart'
c               stop
            endif
         enddo
      enddo
      do j=1,flst_nborn
         do ipart=3,flst_lightpart-1
            if(isalightparton(flst_born(ipart,j))) then
               write(*,*) 
     1      ' genflavreglist: light parton before flst_lightpart'
c               stop
            endif
         enddo
         do ipart=flst_lightpart,nlegborn
            if(.not.isalightparton(flst_born(ipart,j))) then
               write(*,*) 
     1      ' genflavreglist: not a light parton after flst_lightpart'
c               stop
            endif
         enddo
      enddo
c sanity check on real flavour configurations;
c they should all be inequivalent
      do j=1,flst_nreal
         do k=j+1,flst_nreal
            if(flavequiv(nlegreal,nlegreal,j,k,flst_real,
     1           flst_realres,flst_realtags)) then               
               write(*,*)'found two equivalent real flavour processes:'
               write(*,*)'process',j
               call print_lists(nlegreal,flst_real(l,j)
     1              ,flst_realres(l,j),flst_realtags(l,j))
               stop
            endif
         enddo
      enddo
c sanity check on Born flavour configurations;
c they should all be inequivalent
      do j=1,flst_nborn
         do k=j+1,flst_nborn
            if(flavequiv(nlegborn,nlegborn,j,k,flst_born,flst_bornres,
     1           flst_borntags)) then
               write(*,*)'found two equivalent Born flavour processes:'
                write(*,*)'process',j,', flavours '
               call print_lists(nlegborn,flst_born(l,j),
     1              flst_bornres(l,j),flst_borntags(l,j))
               stop
            endif
         enddo
      enddo
c Start search for regions (i.e. alr)
c current number of alr found
      iflregl=0
      flst_nregular=0
      if(flst_nreal.gt.maxprocreal) then
         write(*,*)' genflavreglist: increase maxprocreal'
         stop
      endif
      flg_withreg=.false.
      do k=1,flst_nreal
         call find_regions(flst_real,flst_realres,flst_realtags,
     1        k,nregions,iregions)
         if(nregions.eq.0) then
            flst_nregular=flst_nregular+1
c There are remnants! set up the appropriate flag:
            flg_withreg=.true.
            call intassign
     #(nlegreal,flst_real(1,k),flst_regular(1,flst_nregular))
            call intassign
     #(nlegreal,flst_realres(1,k),flst_regularres(1,flst_nregular))
            call intassign
     #(nlegreal,flst_realtags(1,k),flst_regulartags(1,flst_nregular))
         endif
         do l=1,nregions
            iflregl=iflregl+1
            if(iregions(1,l).le.2) then
               flst_emitter(iflregl)=iregions(1,l)
            else
               flst_emitter(iflregl)=nlegreal-1
            endif
            if(iflregl.ge.maxalr) then
               write(*,*)' genflavreglist: increase maxalr'
               stop
            endif
            ipart=0
c final state singularity
            if(iregions(1,l).gt.2) then
               do j=1,nlegreal
                  if(j.ne.iregions(1,l)
     #                  .and.j.ne.iregions(2,l)) then
                     ipart=ipart+1
                     flst_alr(ipart,iflregl)=flst_real(j,k)
                     flst_alrres(ipart,iflregl)=flst_realres(j,k)
                     flst_alrtags(ipart,iflregl)=flst_realtags(j,k)
                  endif
               enddo
               ipart=ipart+1
               flst_alr(ipart,iflregl)=flst_real(iregions(1,l),k)
               flst_alrres(ipart,iflregl)=flst_realres(iregions(1,l),k)
               flst_alrtags(ipart,iflregl)=
     1              flst_realtags(iregions(1,l),k)
               ipart=ipart+1
               flst_alr(ipart,iflregl)=flst_real(iregions(2,l),k)
               flst_alrres(ipart,iflregl)=flst_realres(iregions(2,l),k)
               flst_alrtags(ipart,iflregl)=
     1              flst_realtags(iregions(2,l),k)
c put always in the order q g and q q~, i.e. fl(i)>fl(j)
               fl1=flst_alr(nlegreal-1,iflregl)
               fl2=flst_alr(nlegreal,iflregl)
               if(  (fl2.ne.22 .and. fl2.ne.0)  .and.
     1           (  (fl1.eq.0 .or. fl1.eq.22) .or.
     1              (fl1.lt.fl2)    )  ) then
                  call exchange_ind(nlegreal,nlegreal-1,
     1                 flst_alr(1,iflregl),flst_alrres(1,iflregl),
     2                 flst_alrtags(1,iflregl))
               endif
            else
c initial state singularity
               do j=1,nlegreal
                  if(j.ne.iregions(2,l)) then
                     ipart=ipart+1
                     flst_alr(ipart,iflregl)=flst_real(j,k)
                     flst_alrres(ipart,iflregl)=flst_realres(j,k)
                     flst_alrtags(ipart,iflregl)=flst_realtags(j,k)
                  endif
               enddo
               ipart=ipart+1
               flst_alr(ipart,iflregl)=flst_real(iregions(2,l),k)
               flst_alrres(ipart,iflregl)=flst_realres(iregions(2,l),k)
               flst_alrtags(ipart,iflregl)=
     1              flst_realtags(iregions(2,l),k)
            endif
c            write(*,*) (flst_alr(ipart,iflregl),ipart=1,nlegreal),
c     #     '   em:',flst_emitter(iflregl)
         enddo
      enddo
      nreg=iflregl
      flst_nalr=nreg
      call pretty_print_flst
c bunch together identical elements, increasing their multiplicities
      do j=1,nreg
         flst_mult(j)=1
      enddo
      do j=1,nreg
         if(flst_mult(j).gt.0) then
            do k=j+1,nreg
c Previously was:
c               if(flst_emitter(j).eq.flst_emitter(k).and.
c     #  equalintlists(nlegreal,flst_alr(1,j),flst_alr(1,k))) then
c now accounts for equivalence by permutation of final state lines.
c Notice: identity of emitter and radiated parton must be valid
c  without permutations
               if(flst_mult(k).ne.0) then
                  if(flst_emitter(j).eq.flst_emitter(k).and.
c     ISR: is ISR, has same radiated parton, is equivalent
c          (excluding the radiated parton)
     1           ( (flst_emitter(j).le.2 .and.
     2              equiv_entry_alr_real(nlegreal,j,k).and.
     3              flavequiv(nlegreal,nlegreal-1,j,k,
     4              flst_alr,flst_alrres,flst_alrtags))
     5                   .or.
c     FSR: has the same radiated and emitter parton, is equivalent
c          (excluding emitter and emitted parton)
     6             (flst_emitter(j).gt.2 .and.
     7              equiv_entry_alr_real(nlegreal,j,k).and.
     8              equiv_entry_alr_real(nlegreal-1,j,k).and.
     9              flavequiv(nlegreal,nlegreal-2,j,k,
     1              flst_alr,flst_alrres,flst_alrtags))
     2           )) then
                     call print_lists(nlegreal,flst_alr(1,j),
     1                 flst_alrres(1,j),flst_alrtags(1,j))
                     call print_lists(nlegreal,flst_alr(1,k),
     1                 flst_alrres(1,k),flst_alrtags(1,k))
                     flst_mult(j)=flst_mult(j)+flst_mult(k)
                     flst_mult(k)=0
                  endif
               endif
            enddo
         endif
      enddo
c browse the list, put together identical elements, compute
c associated underlying Born
      flst_nalr=nreg
      call pretty_print_flst
      iflregl=0
      do j=1,nreg
         if(flst_mult(j).gt.0) then
            iflregl=iflregl+1
            if(j.gt.iflregl) then
               flst_emitter(iflregl)=flst_emitter(j)
               call alr_move(j,iflregl)
               flst_mult(iflregl) = flst_mult(j)
            endif
            call ubornflav(iflregl)
         endif
      enddo
      flst_nalr=iflregl
      call pretty_print_flst
c
c Build unique list of underlying Born; reorder flavours in alpha_r, uborn, emitter
c so that the underlying Born matches exactly a Born flavour structure in the flst_born array
c flavour structures arising as underlying Born
      do j=1,flst_nalr
         do k=1,flst_nborn
c are they the same permutation?
            call reorder_regions(j,k,iret)
c            if(iret.eq.1) write(*,*) ' reordering took place'
            if(iret.ne.-1) goto 11
         enddo
c they are inequivalent
         write(*,*) ' error: underlying born not present in born list'
         call print_lists(nlegborn,flst_uborn(1,j),flst_ubornres(1,j),
     1        flst_uborntags(1,j))
         stop
 11      continue
      enddo
      call pretty_print_flst

c Build pointers from alpha_r -> born
      do j=1,flst_nalr
         do k=1,flst_nborn
            if(equal_lists(nlegborn,j,k,
     1           flst_uborn,flst_ubornres,flst_uborntags(1,j),
     2           flst_born,flst_bornres,flst_borntags)) then
               flst_alr2born(j)=k
            endif
         enddo
      enddo
c Build pointers from born -> alpha_r
      do j=1,flst_nborn
         flst_born2alr(0,j)=0
         do k=1,flst_nalr
            if(equal_lists(nlegborn,k,j,
     1           flst_uborn,flst_ubornres,flst_uborntags(1,j),
     2           flst_born,flst_bornres,flst_borntags)) then
               flst_born2alr(0,j)=flst_born2alr(0,j)+1
               flst_born2alr(flst_born2alr(0,j),j)=k
            endif
         enddo
c Sanity check: each Born should be the underlying Born of some alr
         if(flst_born2alr(0,j).eq.0) then
            write(*,*) ' Born graph ',j,' is never the underlying Born'
     #           //' of some alr'
            call print_lists(nlegborn,flst_born(1,j),
     1           flst_bornres(1,j),flst_borntags(1,j))
c            stop
         endif
      enddo
c Find regions for each alpha_r
      do j=1,flst_nalr
         call find_regions(flst_alr,flst_alrres,flst_alrtags,
     1        j,nregions,iregions)
         do k=1,nregions
            flst_allreg(1,k,j)=iregions(1,k)
            flst_allreg(2,k,j)=iregions(2,k)
         enddo
         flst_allreg(1,0,j)=nregions
      enddo
c For each region, compute the underlying Born multiplicity
      do j=1,flst_nalr
         if(flst_emitter(j).gt.2) then
            flst_ubmult(j)=0
c find flavour of emitter IN THE UNDERLYING BORN
            do k=3,nlegborn
               if(flst_uborn(k,j).eq.flst_uborn(flst_emitter(j),j)
     1  .and.  flst_ubornres(k,j).eq.flst_ubornres(flst_emitter(j),j)
     2  .and.  flst_uborntags(k,j).eq.flst_uborntags(flst_emitter(j),j))
     3              then
                  flst_ubmult(j)=flst_ubmult(j)+1
               endif
            enddo
         else
            flst_ubmult(j)=1
         endif
      enddo
c     debug information
      if (verbose) then
         call pretty_print_flst
      endif
      end

      function equal_lists(n,j,k,a,ares,atags,b,bres,btags)
      logical equal_lists
      integer n,j,k,a(n,*),ares(n,*),atags(n,*),
     1     b(n,*),bres(n,*),btags(n,*)
      integer l
      do l=1,n
         if(a(l,j).ne.b(l,k) .or.
     1      ares(l,j).ne.bres(l,k) .or.
     1      atags(l,j).ne.btags(l,k)) then
            equal_lists=.false.
            return
         endif
      enddo
      equal_lists=.true.
      end

      subroutine from_number_to_madgraph(n,flav,emitter,string)
      implicit none
      integer n,flav(n),emitter
      include 'nlegborn.h'
      character * (*) string
      integer min_partnames,max_partnames
      parameter (min_partnames=-25)
      parameter (max_partnames=25)
      character * 3 partnames(min_partnames:max_partnames)
      data partnames/ 
     &     '   ','W- ','   ','   ','   ','   ','   ','   ','   ',
     $     'vt~','ta+','vm~','mu+','ve~','e+',' ','  ',' ',
     $     '  ','t~','b~','c~','s~','u~','d~','g ','d ','u ','s ' ,'c ',
     $     'b ','t ','  ','  ','',' ','e-','ve','mu-','vmu','ta-','vta',
     $     '  ','  ','  ','  ','  ','gam','   ','W+ ','H  '/
      integer lastnb,j,next
      external lastnb
      string=' '
      if(emitter.eq.0) then
         string='('//partnames(flav(1))//' '//partnames(flav(2))//').'
      elseif(emitter.eq.1) then
         string='('//partnames(flav(1))//')'//partnames(flav(2))//' .'
      elseif(emitter.eq.2) then
         string=' '//partnames(flav(1))//'('//partnames(flav(2))//').'
      else
         string=' '//partnames(flav(1))//' '//partnames(flav(2))//' .'
      endif
      next=lastnb(string)
      string(next:)=' ==> .'
      next=lastnb(string)
      do j=3,n
         if(emitter.eq.j) then
            string(next:)='('//partnames(flav(j))//').'
         else
            string(next:)=' '//partnames(flav(j))//' .'
         endif
         next=lastnb(string) 
      enddo
      string(next:)='    |'
      end

      function lastnb(string)
      implicit none
      integer lastnb
      character *(*) string
      integer ll,l
      ll=len(string)
      do l=ll,1,-1
         if(string(l:l).ne.' ') then
            lastnb=l
            return
         endif
      enddo
      end

      subroutine pretty_print_flst
      implicit none
      include 'nlegborn.h'
      character * 200 string,stringb
      include 'pwhg_flst.h'
      integer j,k,l,iun,lstring,lstringb,lastnb
      external lastnb
c     logical ini
c     data ini/.true./
c     save ini,iun
c     if(ini) then
         call newunit(iun)
         open(unit=iun,file='FlavRegList',status='unknown')
c        ini=.false.
c     endif
c     write(unit=iun,fmt=*) ' index= ',index
      do j=1,flst_nalr
         call from_number_to_madgraph
     #         (nlegreal,flst_alr(1,j),flst_emitter(j),string)
         call from_number_to_madgraph
     #         (nlegborn,flst_uborn(1,j),-1,stringb)
         lstring=lastnb(string)
         lstringb=lastnb(stringb)
         if(flst_alrres(flst_emitter(j),j).ne.0) then
            write(string(lstring:),'(a,i2)')
     1           'res. ',flst_alrres(flst_emitter(j),j)
         endif
         lstring=lastnb(string)
         write(iun,'(a,i3)') string(1:lstring)//' mult=', flst_mult(j)
         write(iun,'(a,i3)') stringb(1:lstringb)//' uborn, mult=',
     1        flst_ubmult(j)
         write(iun,'(20(1x,2(1x,i2)))')
     #   ((flst_allreg(l,k,j),l=1,2),k=1,flst_allreg(1,0,j))
      enddo
      close(iun)
      end


      subroutine intassign(n,iarr1,iarr2)
      implicit none
      integer n,iarr1(n),iarr2(n)
      integer j
      do j=1,n
         iarr2(j)=iarr1(j)
      enddo
      end

      function equalintlists(n,iarr1,iarr2)
      implicit none
      integer n,iarr1(n),iarr2(n)
      logical equalintlists
      integer j
      do j=1,n
         if(iarr2(j).ne.iarr1(j)) then
            equalintlists=.false.
            return
         endif
      enddo
      equalintlists=.true.
      end

      subroutine reorder_regions(alr,iborn,iret)
c      subroutine reorder_regions(n,uborn0,uborn,rflav,emit,iret)
c reorders the amplitude in uborn according to the amplitude in uborn0
c It also reorders rflav, and sets the emitter to its appropriate value
c If the amplitudes have inequivalent flavour structures, it returns -1
c without any other action.
c If the flavour structures are identical, it returns 0, with no other action.
c If the flavour structures have been reordered, it returns 1.
      implicit none
      integer alr,iborn,iret
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer emit
      integer j,k,itmp,ib(nlegborn),ir(nlegborn),it(nlegborn)
      iret=0
      emit=flst_emitter(alr)
      call intassign(nlegborn,flst_uborn(1,alr),ib)
      call intassign(nlegborn,flst_ubornres(1,alr),ir)
      call intassign(nlegborn,flst_uborntags(1,alr),it)
      do j=1,nlegborn
         if(
     1        flst_born(j,iborn).ne.ib(j)    .or.
     2        flst_bornres(j,iborn).ne.ir(j) .or.
     3        flst_borntags(j,iborn).ne.it(j) ) then
            if(j.le.2) then
               iret=-1
               goto 999
            endif
            iret=1
            do k=j+1,nlegborn
               if(
     1              flst_born(j,iborn).eq.ib(k)    .and.
     2              flst_bornres(j,iborn).eq.ir(k) .and.
     3              flst_borntags(j,iborn).eq.it(k) ) then

                  call exchange_ind(j,k,ib,ir,it)

                  goto 10
               endif
            enddo
c they differ in flavour content
            iret=-1
            goto 999
 10         continue
         endif
      enddo
c they are identical; no reordering needed
      if(iret.eq.0) return
c reorder
      do j=3,nlegborn
         
         if(
     1        flst_born(j,iborn).ne.flst_uborn(j,alr)   .or.
     2        flst_bornres(j,iborn).ne.flst_ubornres(j,alr) .or.
     3        flst_borntags(j,iborn).ne.flst_uborntags(j,alr) ) then
            iret=1
            do k=j+1,nlegborn
               if(
     1              flst_born(j,iborn).eq.flst_uborn(k,alr)   .and.
     2              flst_bornres(j,iborn).eq.flst_ubornres(k,alr) .and.
     3              flst_borntags(j,iborn).eq.flst_uborntags(k,alr))then

                  call exchange_ind(j,k,flst_uborn(1,alr),
     1                 flst_ubornres(1,alr),flst_uborntags(1,alr))
                  call exchange_ind(j,k,flst_alr(1,alr),
     1                 flst_alrres(1,alr),flst_alrtags(1,alr))

                  if(flst_emitter(alr).eq.j) then
                     flst_emitter(alr)=k
                  elseif(flst_emitter(alr).eq.k) then
                     flst_emitter(alr)=j
                  endif
                  goto 11
               endif
            enddo
            write(*,*) ' should never get here'
            call exit(-1)
 11         continue
         endif
      enddo
 999  continue
      end

      subroutine exchange_ind(j,k,a,ares,atags)
      implicit none
      integer j,k,a(*),ares(*),atags(*)
      integer itmp
      itmp=a(j)
      a(j)=a(k)
      a(k)=itmp
      itmp=ares(j)
      ares(j)=ares(k)
      ares(k)=itmp
      itmp=atags(j)
      atags(j)=atags(k)
      atags(k)=itmp
      end

      function valid_emitter(j)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      logical valid_emitter
      integer j,alr
      integer validarr(0:nlegborn)
      data validarr/nlegreal*0/
      save validarr
      if(validarr(j).eq.1) then
         valid_emitter=.true.
         return
      elseif(validarr(j).eq.-1) then
         valid_emitter=.false.
         return
      else
         do alr=1,flst_nalr
            if(j.eq.flst_emitter(alr)) then
               valid_emitter=.true.
               validarr(j)=1
               return
            endif
         enddo
         valid_emitter=.false.
         validarr(j)=-1
         return
      endif
      end


      subroutine same_splitting(a,ares,atags,
     1     indexreal,i1,i2,ibornfl,itag,iret)
c returns iret=1 if partons i1,i2 in real indexreal come from the
c same splitting
      implicit none
      character * 3 iof
      integer indexreal,i1,i2,ibornfl,itag,iret
      integer fl1,fl2,tag1,tag2
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer a(nlegreal,*),ares(nlegreal,*),atags(nlegreal,*)
      if(ares(i1,indexreal).ne.
     1     ares(i2,indexreal)) then
c do not come from the same resonance
         iret=-1
         return
      endif
      fl1=a(i1,indexreal)
      fl2=a(i2,indexreal)
      tag1=atags(i1,indexreal)
      tag2=atags(i2,indexreal)
      if(i1.le.2) then
         call same_splitting0('isr',fl1,fl2,tag1,tag2,itag,ibornfl,iret)
      elseif(i2.le.2) then
         call same_splitting0('isr',fl2,fl1,tag2,tag1,itag,ibornfl,iret)
      else
         call same_splitting0('fsr',fl1,fl2,tag1,tag2,itag,ibornfl,iret)
      endif
      end

      subroutine same_splitting0
     1     (iof,fl1,fl2,tag1,tag2,itag,ibornfl,iret)
      implicit none
      character * 3 iof
      integer fl1,fl2,tag1,tag2,itag,ibornfl,iret
      logical is_charged, is_coloured
      external is_charged, is_coloured
      if(iof.eq.'isr'.and.fl2.ne.22) then
c in the isr case, 2 is the outgoing parton
         fl2=-fl2
      endif
      iret=1
      if(fl1+fl2.eq.0.and.is_coloured(fl1).and.tag1.eq.tag2) then
         ibornfl=0
         itag=0
      elseif(fl2.eq.0.and.is_coloured(fl1)) then
         ibornfl=fl1
         itag=tag1
      elseif(fl1.eq.0.and.is_coloured(fl2)) then
         ibornfl=fl2
         itag=tag2
      elseif(fl1.eq.22.and.is_charged(fl2)) then
         ibornfl=fl2
         itag=tag2
      elseif(fl2.eq.22.and.is_charged(fl1)) then
         ibornfl=fl1
         itag=tag1
      else
c cannot come from the same splitting
         iret=-1
      endif
      if(iof.eq.'isr'.and.fl2.ne.22) then
         fl2=-fl2
      endif
      end

      function is_charged(fl)
      implicit none
      logical is_charged
      integer fl
      if(fl.eq.0) then
         is_charged=.false.
      elseif(abs(fl).le.6) then
         is_charged=.true.
      elseif(abs(fl).ge.11.and.abs(fl).le.15.and.2*(fl/2).ne.fl) then
         is_charged=.true.
      else
         is_charged=.false.
      endif
      end

      function is_coloured(fl)
      implicit none
      logical is_coloured
      integer fl
      if(abs(fl).le.6) then
         is_coloured=.true.
      else
         is_coloured=.false.
      endif
      end


      function equiv_entry_alr_real(j,alr1,alr2)
      implicit none 
      logical equiv_entry_alr_real
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer j,alr1,alr2
      equiv_entry_alr_real=
     1     flst_alr(j,alr1).eq.flst_alr(j,alr2)   .and.
     2     flst_alrres(j,alr1).eq.flst_alrres(j,alr2) .and.
     3     flst_alrtags(j,alr1).eq.flst_alrtags(j,alr2)
      end


      subroutine alr_move(j,k)
      implicit none 
      integer j,k
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer i
      do i=1,nlegreal
         flst_alr(i,k) = flst_alr(i,j)
         flst_alrres(i,k) = flst_alrres(i,j)
         flst_alrtags(i,k) = flst_alrtags(i,j)
      enddo
      end

      subroutine print_lists(n,a,ares,atags)
      implicit none
      integer n,a(n),ares(n),atags(n)
      integer j
      character * 50 format
      format = '(         (i3,1x))'
      write(format(2:10),'(i9)') n
      write(*,*)' Flavours:'
      write(*,format) (a(j),j=1,n)
      write(*,*)' Resonance mapping:'
      write(*,format) (ares(j),j=1,n)
      write(*,*)' Tags:'
      write(*,format) (atags(j),j=1,n)
      end
