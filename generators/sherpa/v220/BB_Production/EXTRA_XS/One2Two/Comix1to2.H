#include "METOOLS/Main/Spin_Structure.H"

namespace METOOLS {
  class Current;
  class Vertex;
}

namespace PHASIC {
  class Color_Integrator;
}

namespace EXTRAXS {

  class Comix1to2 : public METOOLS::Spin_Amplitudes {
    std::vector<METOOLS::Current*> m_cur;
    METOOLS::Current *m_fcur;
    std::vector<METOOLS::Vertex*> m_v1;

    std::vector<METOOLS::Current*> m_anticur;
    METOOLS::Current *m_antifcur;
    std::vector<METOOLS::Vertex*> m_antiv1;

    std::vector<size_t> m_nhel;

    PHASIC::Color_Integrator* p_ci;

    size_t NHel(const ATOOLS::Flavour& fl);
  public:
    Comix1to2(const std::vector<ATOOLS::Flavour>& flavs);
    ~Comix1to2();
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
  }; 

}

