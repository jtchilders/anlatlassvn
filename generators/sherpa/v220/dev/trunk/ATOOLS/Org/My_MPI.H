#ifndef ATOOLS__Org__My_MPI_H
#define ATOOLS__Org__My_MPI_H

#include "ATOOLS/Org/CXXFLAGS.H"
#ifdef USING__MPI
#include "mpi.h"
#endif

#include <vector>
#include <string>

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

    static int GetRank(){
      int rank = 0;
#ifdef USING__MPI
      rank = MPI::COMM_WORLD.Get_rank();
#endif
      return rank;
    }

    static int GetWorldSize(){
      int size = 1;
#ifdef USING__MPI
      size = MPI::COMM_WORLD.Get_size();
#endif
      return size;
    }

#ifdef USING__MPI
    static void broadcast_string(MPI::Intracomm comm, int root, std::string& strg) {
      int numprocs = 1, localproc = 0;
      numprocs = comm.Get_size();
      localproc = comm.Get_rank();
      if (numprocs==1)
        return;
      
      if (localproc == root) {
        int length = strg.size();
        MPI_Bcast(&length, 1, MPI::INT, root, comm);
        char* cstr = const_cast<char*>(strg.c_str());
        MPI_Bcast(cstr, strg.size(), MPI::CHAR, root, comm);
      }
      else {
        int length = 0;
        MPI_Bcast(&length, 1, MPI::INT, root, comm);
        char* charstr = new char[length+1];
        MPI_Bcast(charstr, length, MPI::CHAR, root, comm);
        charstr[length] = '\0';
        strg = charstr;
        delete [] charstr;
      }
      return;
    }

    static void broadcast_string_vector(MPI::Intracomm comm, int root, std::vector<std::string>& strg) {
      int numprocs = 1, localproc = 0;
      numprocs = comm.Get_size();
      localproc = comm.Get_rank();
      if (numprocs==1)
        return;

      std::string long_string;
      std::string seperator = ">:>:>:<:<?<";
      
      // add strings together into one string
      if (localproc == root) {
        // catch empty vectors
        if(strg.size() != 0){        
          long_string = strg[0];
          // if the vector is size == 1 then this isn't called
          for(int i=1;i<strg.size();++i)
            long_string = long_string + seperator + strg[i];
        }
      }
      // broadcast the long string
      broadcast_string(comm,root,long_string);
      
      if(long_string.size() == 0) return;

      if( localproc == root ){
        int vector_length = 0;
        vector_length = strg.size();
        MPI_Bcast(&vector_length,1,MPI_INT,root,comm);
      }
      else{
        strg.clear();
        int vector_length = 0;
        MPI_Bcast(&vector_length,1,MPI_INT,root,comm);
        strg.resize(vector_length);

        int substr_start_index = 0;
        int substr_end_index = 0;
        int string_index = 0;
        while((substr_end_index = long_string.find(seperator,substr_start_index)) != std::string::npos){
          strg[string_index++] = long_string.substr(substr_start_index,substr_end_index-substr_start_index);
          substr_start_index = substr_end_index + seperator.size();
        }
        // if the vector was only size == 1 this should be sufficient
        // if the vector is longer, then the last entry is this one.
        strg[string_index] = long_string.substr(substr_start_index);
      }
     return;
    }
#endif

#ifdef USING__MPI
    inline MPI::Intracomm *MPIComm() { return p_comm; }
#endif

  };// end of class My_MPI

  extern My_MPI *mpi;

}// end of namespace ATOOLS

#endif
