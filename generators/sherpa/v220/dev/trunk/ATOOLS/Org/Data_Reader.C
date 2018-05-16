#include "ATOOLS/Org/Data_Reader.H"

#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "My_MPI.H"
#include <typeinfo>
#include <ctype.h>

using namespace ATOOLS;

Data_Reader::Data_Reader(): 
  Read_Write_Base(2,0), m_allowunits(false)
{
  SetInFileMode(fom::permanent);
}

Data_Reader::~Data_Reader() 
{
  CloseInFile(0);
  CloseInFile(1);
}

Data_Reader::Data_Reader(const std::string &wordsep,const std::string &linesep,
			 const std::string &comment,const std::string &ignore):
  Read_Write_Base(2,0,wordsep,linesep,comment,ignore), m_allowunits(false)
{
  SetInFileMode(fom::permanent);
}

void Data_Reader::SetString(const std::string string)
{ 
  m_string=string; 
  FileContent(1).clear();
  AddFileContent(m_string,1);
#ifdef DEBUG__Data_Reader
  msg_IODebugging()<<METHOD<<"(): Read string content '"<<m_string<<"' {\n";
  for (size_t j(0);j<FileContent(1).size();++j) {
    msg_IODebugging()<<"  ";
    for (size_t k(0);k<FileContent(1)[j].size();++k)
      msg_IODebugging()<<"'"<<FileContent(1)[j][k]<<"' ";
    msg_IODebugging()<<"\n";
  }
  msg_IODebugging()<<"}\n";
#endif
}

template <class Read_Type>
Read_Type Data_Reader::Convert(std::string cur) const
{
  if (cur==nullstring) return Default<Read_Type>();
  cur=ReplaceTags(cur);
  Read_Type value;
  if (typeid(value)==typeid(int) || typeid(value)==typeid(unsigned int) ||
      typeid(value)==typeid(long) ||
      typeid(value)==typeid(float) ||	typeid(value)==typeid(double)) {
    if (!AllowNans()) 
      if (cur=="nan" || cur=="inf" || cur=="NAN" || cur=="INF") cur="1";
    if (AllowUnits()) cur=ReplaceUnits(cur);
    if (Interprete()) cur=Interpreter()->Interprete(StripEscapes(cur));
  }
  return ATOOLS::ToType<Read_Type>(cur);
}

template <class Read_Type>
Read_Type Data_Reader::ReadValue(const std::string &parameter,
				 const size_t &file)
{
  std::cerr << METHOD << ":" << My_MPI::GetRank() << ": Start parameter = " << parameter << "\n";
  std::string cur;
#ifdef USING__MPI
  int buffer_size = 0;
  char *buffer = 0;
  if (My_MPI::GetRank() == 0){
    std::cerr << METHOD << ":" << My_MPI::GetRank() << ": Rank 0 moving on to Read.\n";
#endif
  if (file==0) OpenInFile(file);
  for (size_t i(0);i<FileContent(file).size();++i)
    for (size_t j(0);j<FileContent(file)[i].size();++j) {
      size_t pos(0), length(0);
      std::string par(ReplaceTags(FileContent(file)[i][j]));
      if (parameter==nullstring ||
	  (pos=Find(par,parameter,length))!=std::string::npos) {
	cur=par;
	if ((cur=cur.substr(pos+length)).length()==0) {
	  if (j<FileContent(file)[i].size()-1) cur=FileContent(file)[i][j+1];
	  else cur=nullstring;
	}
	if (cur!=nullstring) break;
      }
    }
#ifdef USING__MPI
  } //end if rank 0
  std::cerr << METHOD << ":" << My_MPI::GetRank() << ": cur = " << cur << "\n";
  My_MPI::broadcast_string(MPI::COMM_WORLD,0,cur);
  std::cerr << METHOD << ":" << My_MPI::GetRank() << ": cur = " << cur << "\n";
#endif
  return Convert<Read_Type>(cur);
}

template <class Read_Type> std::vector<Read_Type> 
Data_Reader::ReadVector(const std::string &parameter,const size_t &file)
{
  std::cerr << METHOD << ":" << My_MPI::GetRank() << ": Start parameter = " << parameter << "\n";
  std::vector<Read_Type> values;
  
#ifdef USING__MPI
  std::vector<std::string> str_values;
  if (My_MPI::GetRank() == 0){
    std::cerr << METHOD << ":" << My_MPI::GetRank() << ": Rank 0 is reading file...\n";
#endif

  if (file==0) OpenInFile(file);
  Read_Type value;
  size_t last(0);
  bool found(false);
  for (size_t i(0);i<FileContent(file).size();++i) {
#ifdef USING__MPI
    std::string found_cur;
#endif
    for (size_t j(0);j<FileContent(file)[i].size();++j) {
      size_t pos(0), length(0);
      std::string cur(ReplaceTags(FileContent(file)[i][j]));
      if (parameter==nullstring ||
	  (pos=Find(cur,parameter,length))!=std::string::npos) {
      	if ((cur=cur.substr(pos+length)).length()==0) {
      	  ++j;
      	  if (j<FileContent(file)[i].size()) cur=FileContent(file)[i][j];
      	  else cur="";
      	}
      	if (VectorType()==vtc::vertical) {
#ifdef USING__MPI
          found_cur = cur;
#else
       	  value=Convert<Read_Type>(cur);
#endif
      	  found=true;
      	}
      	else {
      	  if (i>last) values.clear();
#ifdef USING__MPI
          str_values.push_back(cur);
#else
      	  values.push_back(Convert<Read_Type>(cur));
#endif
      	  for (++j;j<FileContent(file)[i].size();++j) 
#ifdef USING__MPI
            str_values.push_back(FileContent(file)[i][j]);
#else
      	    values.push_back(Convert<Read_Type>(FileContent(file)[i][j]));
#endif
      	  last=i;
      	}
      } // end if(par==nullstring || pos!=npos)
    } // end for FileContent(file)[i].size()
    if (VectorType()==vtc::vertical) {
#ifdef USING__MPI
      if (found) str_values.push_back(found_cur);
#else
      if (found) values.push_back(value);
#endif
      found=false;
    }
  } // end for FileContent(file).size()

#ifdef USING__MPI
    std::cerr << METHOD << ":" << My_MPI::GetRank() << ": Rank 0 is done reading the file.\n";
  } // end if(rank == 0)
  std::cerr << METHOD << ":" << My_MPI::GetRank() << ": Ranks will now collect the str_values\n";
  My_MPI::broadcast_string_vector(MPI::COMM_WORLD,0,str_values);
  std::cerr << METHOD << ":" << My_MPI::GetRank() << ": str_values.size() = " << str_values.size() << "\n";
  values.resize(str_values.size());
  for(int i=0;i<str_values.size();++i){
    values[i] = Convert<Read_Type>(str_values[i]);
    std::cerr << METHOD << ":" << My_MPI::GetRank() << ": values[" << i << "] = " << values[i] << "\n";
  }
  std::cerr << METHOD << ":" << My_MPI::GetRank() << ": Done.\n";
#endif

  return values;
}

template <class Read_Type> std::vector< std::vector<Read_Type> > 
Data_Reader::ReadMatrix(const std::string &parameter,const size_t &file)
{
  std::cerr << METHOD << ":" << My_MPI::GetRank() << ": Start.\n";
  if (file==0) OpenInFile(file);
  std::vector< std::vector<Read_Type> > values;
  for (size_t i(0);i<FileContent(file).size();++i) {
    for (size_t j(0);j<FileContent(file)[i].size();++j) {
      size_t pos(0), length(0);
      std::string cur(ReplaceTags(FileContent(file)[i][j]));
      if (parameter==nullstring ||
	  (pos=Find(cur,parameter,length))!=std::string::npos) {
	if ((cur=cur.substr(pos+length)).length()==0) {
	  ++j;
	  if (j<FileContent(file)[i].size()) cur=FileContent(file)[i][j];
	  else cur="";
	}
	values.push_back(std::vector<Read_Type>(1,Convert<Read_Type>(cur)));
	for (++j;j<FileContent(file)[i].size();++j) 
	  values.back().push_back(Convert<Read_Type>
				  (FileContent(file)[i][j]));
      }
    }
  }
  if (MatrixType()==mtc::transposed ||
      values.empty()) return values;
  std::vector< std::vector<Read_Type> > nvalues;
  size_t max(std::numeric_limits<int>::max());
  for (size_t j(0);j<values.size();++j) 
    if (max>values[j].size()) max=values[j].size();
  nvalues.resize(max,std::vector<Read_Type>
		 (values.size(),Default<Read_Type>()));
  for (size_t i(0);i<max;++i)  
    for (size_t j(0);j<values.size();++j) 
      nvalues[i][j]=values[j][i];
  return nvalues;
}

template <class Read_Type > bool Data_Reader::
ReadFromFile(Read_Type &result,std::string parameter) 
{ 
  if ((result=ReadValue<Read_Type>(parameter,0))!=
      Default<Read_Type>()) return true; 
  return false; 
}

template <class Read_Type > bool Data_Reader::
ReadFromString(Read_Type &result,std::string parameter) 
{ 
  if ((result=ReadValue<Read_Type>(parameter,1))!=
      Default<Read_Type>()) return true; 
  return false; 
}

template <class Read_Type > bool Data_Reader::
VectorFromFile(std::vector<Read_Type> &result,std::string parameter) 
{ 
  if ((result=ReadVector<Read_Type>(parameter,0)).size()!=0) return true; 
  return false; 
}

template <class Read_Type> bool Data_Reader::
VectorFromString(std::vector<Read_Type> &result,std::string parameter) 
{ 
  if ((result=ReadVector<Read_Type>(parameter,1)).size()!=0) return true; 
  return false; 
}

template <class Read_Type> bool Data_Reader::
MatrixFromFile(std::vector<std::vector<Read_Type> > &result,
	       std::string parameter) 
{ 
  if ((result=ReadMatrix<Read_Type>(parameter,0)).size()!=0) return true; 
  else return false; 
}

template <class Read_Type> bool Data_Reader::
MatrixFromString(std::vector<std::vector<Read_Type> > &result,
		 std::string parameter) 
{ 
  if ((result=ReadMatrix<Read_Type>(parameter,1)).size()!=0) return true; 
  else return false; 
}

namespace ATOOLS {

  template bool Data_Reader::ReadFromFile<int>
  (int &,std::string);
  template bool Data_Reader::ReadFromFile<unsigned int>
  (unsigned int &,std::string);
  template bool Data_Reader::ReadFromFile<long int>
  (long int &,std::string);
  template bool Data_Reader::ReadFromFile<unsigned long int>
  (unsigned long int &,std::string);
  template bool Data_Reader::ReadFromFile<float>
  (float &,std::string);
  template bool Data_Reader::ReadFromFile<double>
  (double &,std::string);
  template bool Data_Reader::ReadFromFile<std::string>
  (std::string &,std::string);

  template bool Data_Reader::ReadFromString<int>
  (int &,std::string);
  template bool Data_Reader::ReadFromString<unsigned int>
  (unsigned int &,std::string);
  template bool Data_Reader::ReadFromString<long int>
  (long int &,std::string);
  template bool Data_Reader::ReadFromString<unsigned long int>
  (unsigned long int &,std::string);
  template bool Data_Reader::ReadFromString<float>
  (float &,std::string);
  template bool Data_Reader::ReadFromString<double>
  (double &,std::string);
  template bool Data_Reader::ReadFromString<std::string>
  (std::string &,std::string);

  template bool Data_Reader::VectorFromFile<int>
  (std::vector<int> &,std::string);
  template bool Data_Reader::VectorFromFile<unsigned int>
  (std::vector<unsigned int> &,std::string);
  template bool Data_Reader::VectorFromFile<long int>
  (std::vector<long int> &,std::string);
  template bool Data_Reader::VectorFromFile<unsigned long int>
  (std::vector<unsigned long int> &,std::string);
  template bool Data_Reader::VectorFromFile<float>
  (std::vector<float> &,std::string);
  template bool Data_Reader::VectorFromFile<double>
  (std::vector<double> &,std::string);
  template bool Data_Reader::VectorFromFile<std::string>
  (std::vector<std::string> &,std::string);

  template bool Data_Reader::VectorFromString<int>
  (std::vector<int> &,std::string);
  template bool Data_Reader::VectorFromString<unsigned int>
  (std::vector<unsigned int> &,std::string);
  template bool Data_Reader::VectorFromString<long int>
  (std::vector<long int> &,std::string);
  template bool Data_Reader::VectorFromString<unsigned long int>
  (std::vector<unsigned long int> &,std::string);
  template bool Data_Reader::VectorFromString<float>
  (std::vector<float> &,std::string);
  template bool Data_Reader::VectorFromString<double>
  (std::vector<double> &,std::string);
  template bool Data_Reader::VectorFromString<std::string>
  (std::vector<std::string> &,std::string);

  template bool Data_Reader::MatrixFromFile<int>
  (std::vector<std::vector<int> > &,std::string); 
  template bool Data_Reader::MatrixFromFile<unsigned int>
  (std::vector<std::vector<unsigned int> > &,std::string);
  template bool Data_Reader::MatrixFromFile<long int>
  (std::vector<std::vector<long int> > &,std::string);
  template bool Data_Reader::MatrixFromFile<unsigned long int>
  (std::vector<std::vector<unsigned long int> > &,std::string);
  template bool Data_Reader::MatrixFromFile<float>
  (std::vector<std::vector<float> > &,std::string);
  template bool Data_Reader::MatrixFromFile<double>
  (std::vector<std::vector<double> > &,std::string);
  template bool Data_Reader::MatrixFromFile<std::string>
  (std::vector<std::vector<std::string> > &,std::string);

  template bool Data_Reader::MatrixFromString<int>
  (std::vector<std::vector<int> > &,std::string);
  template bool Data_Reader::MatrixFromString<unsigned int>
  (std::vector<std::vector<unsigned int> > &,std::string);
  template bool Data_Reader::MatrixFromString<long int>
  (std::vector<std::vector<long int> > &,std::string);
  template bool Data_Reader::MatrixFromString<unsigned long int>
  (std::vector<std::vector<unsigned long int> > &,std::string);
  template bool Data_Reader::MatrixFromString<float>
  (std::vector<std::vector<float> > &,std::string);
  template bool Data_Reader::MatrixFromString<double>
  (std::vector<std::vector<double> > &,std::string);
  template bool Data_Reader::MatrixFromString<std::string>
  (std::vector<std::vector<std::string> > &,std::string);

}
