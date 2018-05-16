      implicit none
      character * 60 file
      character * 100 string
      integer j,k,lt
      integer nobl
      character * 3 tag
      write(*,*) ' enter 1 for LO, 2 for NLO'
      read(*,*) lt
      if(lt.eq.1) then
         tag='LO'
         lt=2
      else
         tag='NLO'
         lt=3
      endif
      file='pwg'//tag(1:lt)//'.top'
      open(unit=10,file=file,status='old')
      open(unit=20,file='cmp.top',status='unknown')
      do j=1,1000000
         read(unit=10,fmt='(a)',end=999,err=999) string
         write(20,'(a)') string(1:nobl(string))
         if(string.eq.' (  PT jet 1') then
            file='pt7_'//tag(1:lt)//'_150.top'
         elseif(string.eq.' (  PT jet 2') then
            file='pt8_'//tag(1:lt)//'_150.top'
         elseif(string.eq.' (  Eta jet 1') then
            file='eta7_'//tag(1:lt)//'_150.top'
         elseif(string.eq.' (  Eta jet 2') then
            file='eta8_'//tag(1:lt)//'_150.top'
         elseif(string.eq.' (  Etaj1-Etaj2') then
            file='eta78_'//tag(1:lt)//'_150.top'
         elseif(string.eq.' (  HT,TOT') then
            file='HTTOT_'//tag(1:lt)//'_150.top'
         elseif(string.eq.' (  PT LEPT') then
            file='pt4_'//tag(1:lt)//'_150.top'
         elseif(string.eq.' (  PT MISS') then
            file='etmiss_'//tag(1:lt)//'_150.top'
         elseif(string.eq.' (  Eta_lept') then
            file='eta4_'//tag(1:lt)//'_150.top'
         elseif(string.eq.' (  Eta_l1-eta_l2') then
            file='eta45_'//tag(1:lt)//'_150.top'
         elseif(string.eq.' (  M(l1l2)') then
            file='mll_'//tag(1:lt)//'_150.top'
         elseif(string.eq.' (  MT_WW') then
            file='mtww_'//tag(1:lt)//'_150.top'
         elseif(string.eq.' (  DelR(j1lep)') then
            file='Rlj1_'//tag(1:lt)//'_150.top'
         elseif(string.eq.' (  DelR(j2lept)') then
            file='Rlj2_'//tag(1:lt)//'_150.top'
         elseif(string.eq.' (  Phi(l1l2)') then
            file='phill_'//tag(1:lt)//'_150.top'
         else
            goto 30
         endif
         open(unit=11,file=file,status='old')
         write(20,*) ' set order x y'
         do k=1,10000
            read(unit=11,fmt='(a)',end=22,err=22) string
            write(20,'(a)') string(1:nobl(string))
         enddo
 22      write(20,*) ' plot red symbol=5O'
 30      continue
      enddo
 999  end

               
            
      function nobl(string)
      implicit none
      integer nobl
      character *(*) string
      integer l
      l=len(string)
      do nobl=l,1,-1
        if(string(nobl:nobl).ne.' ') return
      enddo
      end
