c     various shared helpers
c     2012-07 AvM


c     decode two PDG ids for MSSM fermions from a single integer
      subroutine decode_sfermion_pair(combination,ida,idb)
      implicit none
      integer combination,ida,idb,idaa,idbb
      idaa = combination/1000
      idbb = mod(combination,1000)
      ida =  (idaa/100)*1000000 + mod(idaa,100)
      idb = -(idbb/100)*1000000 - mod(idbb,100)
      end



c     encode two PDG ids for MSSM fermions into a single integer
c     (encodes 2000013,-200013 to 213213 etc, opposite signs for A,B)
      integer function encode_sfermion_pair(ida,idb)
      implicit none
      integer ida,idb,idaa,idbb,combination,idar,idbr
      idaa = (ida / 1000000)*100 + mod(ida,100)
      idbb = (idb / 1000000)*100 + mod(idb,100)
      combination = idaa - 1000*idbb
      ! check if reconstruction works, then we are fine in any case
      call decode_sfermion_pair(combination,idar,idbr)
      if ((ida.ne.idar).or.(idb.ne.idbr))  stop "invalid particle ID"
      encode_sfermion_pair = combination
      end
