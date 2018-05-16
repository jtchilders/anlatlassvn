      subroutine reorderlord(nwz,i1,i2,i5in,i6in,ii5,ii6)
      implicit none
C-----reordering rules.
C     if final state has qqb, q always comes first, except for gg initial state
C     if final state has (qqb or q) and g  then g always comes last
C     if final state has (q,q) or (qqb,qqb) ordering is matched 
C-----with identical q or qb in initial state
      integer i1,i2,i5,i6,i5in,i6in,ii5,ii6,nwz
      i5=i5in
      i6=i6in
      if ((i1 .eq. 0) .and. (i2 .eq. 0)) then
          ii5=min(i5,i6)
          ii6=max(i5,i6)
      elseif (((i6 .gt. 0) .and. (i5 .lt. 0)) .or. (i5 .eq. 0)) then
          ii5=i6
          ii6=i5
      elseif (((i5 .gt. 0) .and. (i6 .lt. 0)) .or. (i6 .eq. 0)) then
          ii5=i5 
          ii6=i6 
      elseif (((i5 .gt. 0) .and. (i6 .gt. 0)) 
     &   .or. ((i5 .lt. 0) .and. (i6 .lt. 0)))  then
          if (nwz .eq. 1) then 
              if ((i5 .eq. i1) .or. (i6 .eq. i2)) then
                ii5=i5
                ii6=i6
              elseif ((i5 .eq. i2) .or. (i6.eq.i1)) then
                ii5=i6
                ii6=i5
              endif
          elseif (nwz .eq. -1) then 
               if ((i5 .eq. i2) .or. (i6.eq.i1)) then
                 ii5=i6
                 ii6=i5
               elseif ((i5 .eq. i1) .or. (i6 .eq. i2)) then
                 ii5=i5
                 ii6=i6
               endif
          endif
      else
          ii5=max(i5,i6)
          ii6=min(i5,i6)
      endif

      return
      end
