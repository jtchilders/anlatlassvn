      subroutine init_hist
      call pwhginihist
      end

      subroutine analysis(dsig)
      real * 8 dsig
      end

      subroutine particle_identif(HWW,HZZ)
      implicit none
      integer pdg_Higgs,pdg_Z,pdg_W,HZZ,HWW
      pdg_Higgs = 25
      pdg_Z = 23
      pdg_W = 24      
c     build an identifier for Higgs production in WW and ZZ fusion 
      HWW = 10000*pdg_W + pdg_Higgs
      HZZ = 10000*pdg_Z + pdg_Higgs
      end


