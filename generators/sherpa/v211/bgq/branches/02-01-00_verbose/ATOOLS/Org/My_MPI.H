#ifndef ATOOLS__Org__My_MPI_H
#define ATOOLS__Org__My_MPI_H

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/MyTiming.H"
#include "ATOOLS/Org/Message.H"
#ifdef USING__MPI
#include "mpi.h"
#endif

#include <vector>
#include <typeinfo>
#include <cxxabi.h>
#include <sstream>

namespace ATOOLS {

  class Data_Reader;

  class My_MPI {
  private:

    bool m_hassend, m_hasrecv;

#ifdef USING__MPI
    MPI::Intracomm m_send, m_recv;
#endif

    void SetMPIRecv(std::vector<int> r);

  public:

    My_MPI();

    ~My_MPI();

    void SetUpSendRecv(Data_Reader *const read);

    bool HasMPISend() const;
    bool HasMPIRecv() const;

#ifdef USING__MPI
    inline MPI::Intracomm &MPISend() { return m_send; }
    inline MPI::Intracomm &MPIRecv() { return m_recv; }

    static int num_bcast_string;
    static int num_bcast;
    static int num_bcast_vect;
    static int num_bcast_matrix;

    template <class Read_Type>
    static void Bcast(Read_Type& val){
       num_bcast++;
       const std::type_info  &ti = typeid(val);
       int status = 0;
       char *realname = abi::__cxa_demangle(ti.name(), 0, 0, &status);
       std::stringstream ss;
       ss << "type=" << realname << " calls=" << num_bcast << " starting val=" << val;
       My_MPI::PrintMessage(__PRETTY_FUNCTION__,ss.str());
       MPI::COMM_WORLD.Bcast(&val,1,GetDataType<Read_Type>(),0);
       ss.str("");
       ss << "type=" << realname << " calls=" << num_bcast << " ending val=" << val;
       My_MPI::PrintMessage(__PRETTY_FUNCTION__,ss.str());
    }

    template <class Vector_Type>
    static void Bcast(std::vector<Vector_Type>& vect){
      num_bcast_vect++;
      const std::type_info  &ti = typeid(vect);
      int status = 0;
      char *realname = abi::__cxa_demangle(ti.name(), 0, 0, &status);
      int vector_size = vect.size();
      
      std::stringstream ss;
      ss << "vector_type=" << realname << " calls=" << num_bcast_vect << " starting size=" << vector_size;
      My_MPI::PrintMessage(__PRETTY_FUNCTION__,ss.str());
      Bcast<int>(vector_size);
      ss.str("");
      ss << "vector_type=" << realname << " calls=" << num_bcast_vect << " middle size=" << vector_size;
      My_MPI::PrintMessage(__PRETTY_FUNCTION__,ss.str());
      
      for(unsigned int i=0;i<vector_size;++i){
         if(MPI::COMM_WORLD.Get_rank()!=0) vect.push_back(Vector_Type());
         Bcast<Vector_Type>(vect[i]);
      }
      ss.str("");
      ss << "vector_type=" << realname << " calls=" << num_bcast_vect << " ending size=" << vect.size();
      My_MPI::PrintMessage(__PRETTY_FUNCTION__,ss.str());
      
    }

    template <class Matrix_Type>
    static void Bcast(std::vector<std::vector<Matrix_Type> >& matrix){
      num_bcast_matrix++;
      const std::type_info  &ti = typeid(matrix);
      int status = 0;
      char *realname = abi::__cxa_demangle(ti.name(), 0, 0, &status);
      int matrix_size = matrix.size();

      std::stringstream ss;
      ss << "matrix_type=" << realname << " calls=" << num_bcast_matrix << " starting size=" << matrix_size;
      My_MPI::PrintMessage(__PRETTY_FUNCTION__,ss.str());
      Bcast<int>(matrix_size);
      ss.str("");
      ss << "matrix_type=" << realname << " calls=" << num_bcast_matrix << " middle size=" << matrix_size;
      My_MPI::PrintMessage(__PRETTY_FUNCTION__,ss.str());
      
      for(unsigned int i=0;i<matrix_size;++i){
         if(MPI::COMM_WORLD.Get_rank()!=0) matrix.push_back(std::vector<Matrix_Type>());
         Bcast<Matrix_Type>(matrix[i]);
      }
      ss.str("");
      ss << "matrix_type=" << realname << " calls=" << num_bcast_matrix << " ending size=" << matrix.size();
      My_MPI::PrintMessage(__PRETTY_FUNCTION__,ss.str());
      
    }
    

    template <class Read_Type>
    static MPI::Datatype GetDataType(){
      MPI::Datatype datatype = MPI::INT;
      if (typeid(int) == typeid(Read_Type))
         datatype = MPI::INT;
      else if(typeid(short) == typeid(Read_Type))
         datatype = MPI::SHORT;
      else if(typeid(long) == typeid(Read_Type))
         datatype = MPI::LONG;
      else if(typeid(unsigned) == typeid(Read_Type))
         datatype = MPI::UNSIGNED;
      else if(typeid(unsigned short) == typeid(Read_Type))
         datatype = MPI::UNSIGNED_SHORT;
      else if(typeid(unsigned long) == typeid(Read_Type))
         datatype = MPI::UNSIGNED_LONG;
      else if(typeid(float) == typeid(Read_Type))
         datatype = MPI::FLOAT;
      else if(typeid(double) == typeid(Read_Type))
         datatype = MPI::DOUBLE;
      else if(typeid(char) == typeid(Read_Type))
         datatype = MPI::CHAR;
      else if(typeid(unsigned char) == typeid(Read_Type))
         datatype = MPI::UNSIGNED_CHAR;
      return datatype;
    };

    static void PrintMessage(std::string class_method,std::string msg);

    static int GetRank();

#endif

  };// end of class My_MPI

  template <>
  void My_MPI::Bcast<std::string>(std::string& val);

  extern My_MPI *mpi;

}// end of namespace ATOOLS

#endif
