      subroutine reorderreal(nwz,i1,i2,i5in,i6in,i7in,ii5,ii6,ii7)
      implicit none
C-----reordering rules.
C     if final state has qqb, q always comes first, except for gg initial state
C     if final state has (qqb or q) and g  then g always comes last
C     if final state has (q,q) or (qqb,qqb) ordering is matched 
C-----with identical q or qb in initial state
      integer i1,i2,i5,i6,i7,ii5,ii6,ii7,temp,imax,imin,imid
      integer i5in,i6in,i7in,nwz
      i5=i5in
      i6=i6in
      i7=i7in
         if (i5 .eq. 0) then
         temp=i5         
         i5=i7         
         i7=temp
         endif
         if (i6 .eq. 0) then
         temp=i6         
         i6=i7         
         i7=temp
         endif
         imax=max(i5,i6,i7)
         imin=min(i5,i6,i7)
         imid=i5+i6+i7-imax-imin
         if ((imid .gt. 0) .and. (imax .gt.0) .and. (i7 .ne. 0)) then
         i7=imin
         i6=imid
         i5=imax
         endif
         if ((imid .lt. 0) .and. (imin .lt.0) .and. (i7 .ne. 0)) then
         i7=imax
         i6=imid
         i5=imin
         endif

      if ((i1 .eq. 0) .and. (i2 .eq. 0)) then
          ii5=min(i5,i6)
          ii6=max(i5,i6)
          ii7=i7
c          write(6,*) '1,ii5,ii6,ii7',ii5,ii6,ii7
      elseif (((i6 .gt. 0) .and. (i5 .lt. 0)) .or. (i5 .eq. 0)) then
          ii5=i6
          ii6=i5
          ii7=i7
c          write(6,*) '2,ii5,ii6,ii7',ii5,ii6,ii7
      elseif (((i5 .gt. 0) .and. (i6 .lt. 0)) .or. (i6 .eq. 0)) then
          ii5=i5 
          ii6=i6 
          ii7=i7
c          write(6,*) '3,ii5,ii6,ii7',ii5,ii6,ii7
      elseif ((i5 .gt. 0) .and. (i6 .gt. 0)) then
          if (nwz .eq. 1) then 
              if ((i5 .eq. i2) .or. (i6 .eq. i1)) then
                ii5=i6
                ii6=i5
                ii7=i7
              elseif ((i5 .eq. i1) .or. (i6.eq.i2)) then
                ii5=i5
                ii6=i6
                ii7=i7
              else
                ii5=max(i5,i6)
                ii6=min(i5,i6)
                ii7=i7
              endif
          elseif (nwz .eq. -1) then 
               if ((i5 .eq. i1) .or. (i6.eq.i2)) then
                 ii5=i5
                 ii6=i6
                 ii7=i7
               elseif ((i5 .eq. i2) .or. (i6 .eq. i1)) then
                 ii5=i6
                 ii6=i5
                 ii7=i7
               else
                 ii5=max(i5,i6)
                 ii6=min(i5,i6)
                 ii7=i7
               endif
          endif
      elseif ((i5 .lt. 0) .and. (i6 .lt. 0))  then
          if (nwz .eq. 1) then 
              if ((i5 .eq. i1) .or. (i6 .eq. i2)) then
                ii5=i5
                ii6=i6
                ii7=i7
              elseif ((i5 .eq. i2) .or. (i6.eq.i1)) then
                ii5=i6
                ii6=i5
                ii7=i7
              else
                ii5=max(i5,i6)
                ii6=min(i5,i6)
                ii7=i7
              endif
          elseif (nwz .eq. -1) then 
               if ((i5 .eq. i2) .or. (i6.eq.i1)) then
                 ii5=i6
                 ii6=i5
                 ii7=i7
               elseif ((i5 .eq. i1) .or. (i6 .eq. i2)) then
                 ii5=i5
                 ii6=i6
                 ii7=i7
               else
                 ii5=max(i5,i6)
                 ii6=min(i5,i6)
                 ii7=i7
               endif
          endif
      else
          ii5=max(i5,i6)
          ii6=min(i5,i6)
          ii7=i7
c          write(6,*) '5,ii5,ii6,ii7',ii5,ii6,ii7
      endif

C------adjust special cases
      if ((i1.eq.-3).and.(i2.eq.0)
     & .and.(ii5.eq.-2).and.(ii6.eq.-4).and.(ii7.eq.2)) then
      ii5=-4
      ii6=-2
      elseif ((i1.eq.-3).and.(i2.eq.0)
     & .and.(ii5.eq.-1).and.(ii6.eq.-4).and.(ii7.eq.1)) then
      ii5=-4
      ii6=-1
      elseif ((i1.eq.0).and.(i2.eq.-3)
     & .and.(ii5.eq.-4).and.(ii6.eq.-5).and.(ii7.eq.5)) then
      ii5=-5
      ii6=-4
      elseif ((i1.eq.0).and.(i2.eq.-1)
     & .and.(ii5.eq.-2).and.(ii6.eq.-5).and.(ii7.eq.5)) then
      ii5=-5
      ii6=-2
      elseif ((i1.eq.0).and.(i2.eq.-1)
     & .and.(ii5.eq.-2).and.(ii6.eq.-4).and.(ii7.eq.4)) then
      ii5=-4
      ii6=-2
      elseif ((i1.eq.0).and.(i2.eq.-1)
     & .and.(ii5.eq.-2).and.(ii6.eq.-3).and.(ii7.eq.3)) then
      ii5=-3
      ii6=-2
      elseif ((i1.eq.0).and.(i2.eq.4)
     & .and.(ii5.eq.3).and.(ii6.eq.1).and.(ii7.eq.-1)) then
      ii5=1
      ii6=3
      elseif ((i1.eq.0).and.(i2.eq.4)
     & .and.(ii5.eq.3).and.(ii6.eq.2).and.(ii7.eq.-2)) then
      ii5=2
      ii6=3
      elseif ((i1.eq.2).and.(i2.eq.0)
     & .and.(ii5.eq.5).and.(ii6.eq.1).and.(ii7.eq.-5)) then
      ii5=1
      ii6=5
      elseif ((i1.eq.2).and.(i2.eq.0)
     & .and.(ii5.eq.4).and.(ii6.eq.1).and.(ii7.eq.-4)) then
      ii5=1
      ii6=4
      elseif ((i1.eq.2).and.(i2.eq.0)
     & .and.(ii5.eq.3).and.(ii6.eq.1).and.(ii7.eq.-3)) then
      ii5=1
      ii6=3
      elseif ((i1.eq.4).and.(i2.eq.0)
     & .and.(ii5.eq.+5).and.(ii6.eq.3).and.(ii7.eq.-5)) then
      ii5=3
      ii6=5
      endif

      if ((i1.eq.-4).and.(i2.eq.0)
     & .and.(ii5.eq.-2).and.(ii6.eq.-3).and.(ii7.eq.2)) then
      ii5=-3
      ii6=-2
      elseif ((i1.eq.-4).and.(i2.eq.0)
     & .and.(ii5.eq.-1).and.(ii6.eq.-3).and.(ii7.eq.1)) then
      ii5=-3
      ii6=-1
      elseif ((i1.eq.0).and.(i2.eq.-4)
     & .and.(ii5.eq.-3).and.(ii6.eq.-5).and.(ii7.eq.5)) then
      ii5=-5
      ii6=-3
      elseif ((i1.eq.0).and.(i2.eq.-2)
     & .and.(ii5.eq.-1).and.(ii6.eq.-5).and.(ii7.eq.5)) then
      ii5=-5
      ii6=-1
      elseif ((i1.eq.0).and.(i2.eq.-2)
     & .and.(ii5.eq.-1).and.(ii6.eq.-4).and.(ii7.eq.4)) then
      ii5=-4
      ii6=-1
      elseif ((i1.eq.0).and.(i2.eq.-2)
     & .and.(ii5.eq.-1).and.(ii6.eq.-3).and.(ii7.eq.3)) then
      ii5=-3
      ii6=-1
      elseif ((i1.eq.0).and.(i2.eq.3)
     & .and.(ii5.eq.4).and.(ii6.eq.2).and.(ii7.eq.-2)) then
      ii5=2
      ii6=4
      elseif ((i1.eq.0).and.(i2.eq.3)
     & .and.(ii5.eq.4).and.(ii6.eq.1).and.(ii7.eq.-1)) then
      ii5=1
      ii6=4
      elseif ((i1.eq.1).and.(i2.eq.0)
     & .and.(ii5.eq.5).and.(ii6.eq.2).and.(ii7.eq.-5)) then
      ii5=2
      ii6=5
      elseif ((i1.eq.1).and.(i2.eq.0)
     & .and.(ii5.eq.4).and.(ii6.eq.2).and.(ii7.eq.-4)) then
      ii5=2
      ii6=4
      elseif ((i1.eq.1).and.(i2.eq.0)
     & .and.(ii5.eq.3).and.(ii6.eq.2).and.(ii7.eq.-3)) then
      ii5=2
      ii6=3
      elseif ((i1.eq.3).and.(i2.eq.0)
     & .and.(ii5.eq.5).and.(ii6.eq.4).and.(ii7.eq.-5)) then
      ii5=4
      ii6=5
      endif
      return
      end

