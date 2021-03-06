#ifndef ATOOLS__Org__My_MPI_H
#define ATOOLS__Org__My_MPI_H

#include "ATOOLS/Org/CXXFLAGS.H"
#ifdef USING__MPI
#include "mpi.h"
#endif

#include <vector>

namespace ATOOLS {

  class Data_Reader;

  class My_MPI {
  private:

#ifdef USING__MPI
    MPI::Intracomm *p_comm;
#endif

    void SetMPIRecv(std::vector<int> r);

  public:

    My_MPI();

    ~My_MPI();

    void SetUpSendRecv(Data_Reader *const read);

#ifdef USING__MPI
    inline MPI::Intracomm *MPIComm() { return p_comm; }
#endif

  };// end of class My_MPI

  extern My_MPI *mpi;

}// end of namespace ATOOLS

#endif
