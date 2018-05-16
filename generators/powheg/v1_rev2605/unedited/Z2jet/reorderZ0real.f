      subroutine reorderZ0real(i1,i2,i5in,i6in,i7in,ii5,ii6,ii7)
      implicit none
C-----reordering rules.
C     if final state has qqb, q always comes first, except for gg initial state
C     if final state has (qqb or q) and g  then g always comes last
C     if final state has (q,q) or (qqb,qqb) ordering is matched 
C-----with identical q or qb in initial state
      integer i1,i2,i5,i6,i7,ii5,ii6,ii7,temp,imax,imin,imid,j
      integer i5in,i6in,i7in,nwz
      logical qqa,qaa
      data j/0/
      j=j+1
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
      
C-at this point gluon is last
      imax=max(i5,i6,i7)
      imin=min(i5,i6,i7)
      imid=i5+i6+i7-imax-imin
      qqa=(imax > 0) .and. (imid > 0) .and. (imin < 0)
      qaa=(imax > 0) .and. (imin < 0) .and. (imid < 0)
      if (i7 .ne. 0) then
            if (qqa) then
               if (imax == -imin) then
                  ii5=imax
                  ii6=imin
                  ii7=imid
               elseif (imid == -imin) then
                  ii5=imid
                  ii6=imin
                  ii7=imax
               endif
            elseif (qaa) then
              if (imax == -imin) then
                  ii5=imax
                  ii6=imin
                  ii7=imid
              elseif (imax == -imid) then
                  ii5=imax
                  ii6=imid
                  ii7=imin
              endif
            endif
      elseif (i7 .eq. 0) then
          ii7=i7
          if ((i1 == 0) .and. (i2 == 0)) then
          ii5=min(i5,i6)
          ii6=max(i5,i6)
          else
          if ((i5 .gt. 0) .and. (i6 .gt. 0)) then
              if ((i5 .eq. i2) .or. (i6 .eq. i1)) then
                ii5=i6
                ii6=i5
              elseif ((i5 .eq. i1) .or. (i6.eq.i2)) then
                ii5=i5
                ii6=i6
              else
                ii5=max(i5,i6)
                ii6=min(i5,i6)
              endif
          elseif ((i5 .lt. 0) .and. (i6 .lt. 0))  then
              if ((i5 .eq. i2) .or. (i6 .eq. i1)) then
                ii5=i6
                ii6=i5
              elseif ((i5 .eq. i1) .or. (i6.eq.i2)) then
                ii5=i5
                ii6=i6
              else
                ii5=max(i5,i6)
                ii6=min(i5,i6)
              endif
          elseif ((i5 .gt. 0) .and. (i6 .lt. 0))  then
                ii5=max(i5,i6)
                ii6=min(i5,i6)
          elseif ((i5 .lt. 0) .and. (i6 .gt. 0))  then
                ii5=max(i5,i6)
                ii6=min(i5,i6)
          else
              ii5=i5
              ii6=i6
          endif
          endif
      endif
      return
      end

