#include "ATOOLS/Org/My_MPI.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"

#include <stddef.h>
#include <cstring>
#include <unistd.h>

using namespace ATOOLS;

My_MPI *ATOOLS::mpi(NULL);

My_MPI::My_MPI():
  m_hassend(false), m_hasrecv(false)
{
}

My_MPI::~My_MPI()
{
#ifdef USING__MPI
  if (m_hassend) m_send.Free();
  if (m_hasrecv) m_recv.Free();
#endif  
}

void My_MPI::SetUpSendRecv(Data_Reader *const read)
{
#ifdef USING__MPI
  int size=MPI::COMM_WORLD.Get_size();
  if (size>1) {
    std::string rname=read->GetValue<std::string>
      ("MPI_NODE_NAME",std::string(".*"));
    double starttime=rpa->gen.Timer().RealTime();
    msg_Info()<<METHOD<<"(): Analyzing MPI environment {\n";
    int rank=MPI::COMM_WORLD.Get_rank(), pid(getpid()), hlen;
    char host[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(host,&hlen);
    std::vector<std::string> hres(RegExMatch(host,rname));
    if (hres.size()) strcpy(host,hres.front().c_str());
    // Rank 0 keeps track of all the other ranks
    if (rank==0) {
      msg_Info()<<"  Rank "<<rank<<", pid "<<pid
		<<" running on "<<host<<".\n";
      char mhost[MPI_MAX_PROCESSOR_NAME];
      strcpy(mhost,host);
      //std::vector<std::string> hres(RegExMatch(mhost,rname));
      //if (hres.size()) strcpy(mhost,hres.front().c_str());
      std::vector<std::string> hosts(size,host); // the hostname indexed by rank
      std::map<std::string,size_t> procs; // the number of ranks for each host, indexed by host 
      std::vector<std::vector<int> > recv(size); // vector of ranks, indexed by rank
      std::map<std::string,std::vector<int> > recvmap; // vector of rank numbers for each host, indexed by host
      std::vector<int> mrecv(1,-1); // vector of rank numbers, not sure what this is 
      procs[mhost]=1;
      std::map<std::string,size_t>::const_iterator mhost_procs_itr = procs.find(mhost);
      for (int tag=1;tag<size;++tag) { // loop over ranks
	MPI::COMM_WORLD.Recv(&pid,1,MPI::INT,MPI::ANY_SOURCE,tag);
	MPI::COMM_WORLD.Recv(host,MPI_MAX_PROCESSOR_NAME,
			     MPI::CHAR,MPI::ANY_SOURCE,tag);
	msg_Info()<<"  Rank "<<tag<<", pid "<<pid
		  <<" running on "<<host<<"."<<std::endl;
	hosts[tag]=host;
	recvmap[host].push_back(tag);
	if (procs.find(host)!=procs.end()) { // if not a new host received
	  if (strcmp(host,mhost)==0 && mhost_procs_itr->second==1) // if rank is on mhost and first additional rank on this host ( host == mhost and proc[mhost] == 1)
	    mrecv.push_back(tag);
	  ++procs[host];
	}
	else { // if new host received
	  mrecv.push_back(tag);
	  procs[host]=1;
	}
      } // end for(tag)
      // if there is only one host machine
      if (procs.size()==1) {
           std::map<std::string,std::vector<int> >::iterator recvmap_itr = recvmap.find(mhost);
           if(recvmap_itr != recvmap.end()){
   // add main host rank numbers at the end of mrecv
	mrecv.insert(mrecv.end(),
		     recvmap_itr->second.begin(),
		     recvmap_itr->second.end());
   // resize host rank numbers to one, erase old ones.
	recvmap_itr->second.resize(1);
           }
      }
      std::set<std::string> send;
      for (int tag=1;tag<size;++tag) { // loop over ranks
	std::vector<int> recv(1,recvmap[hosts[tag]][0]); // 0th rank for the host on which this rank (tag) is located. (overrided previous 'recv' variable)
	if (send.find(hosts[tag])==send.end()) { // if host is new
	  send.insert(hosts[tag]); // add hostname only once
	  recv=recvmap[hosts[tag]]; // get rank's host's vector of ranks
	  recv[0]=0; // set first rank to 0
	}
	if (recv.size()>1)
	  msg_Info()<<"  Rank "<<tag<<" send/recv "<<recv<<".\n";
	int nrecv(recv.size());
	MPI::COMM_WORLD.Send(&nrecv,1,MPI::INT,tag,size+tag);
	MPI::COMM_WORLD.Send(&recv.front(),nrecv,MPI::INT,tag,size+tag);
   for(int aa=0;aa<recv.size();aa++)
      msg_MPIDebugging() << " rank " << tag << " recv[ " << aa << "] = " << recv[aa] << std::endl;
      } // end for(tag)
      for(int aa=0;aa<mrecv.size();aa++)
         msg_MPIDebugging() << " mrecv[ " << aa << "] = " << mrecv[aa] << std::endl;
      SetMPIRecv(mrecv);
      double diff=rpa->gen.Timer().RealTime()-starttime;
      msg_Info()<<"} -> "<<FormatTime(size_t(diff))<<" elapsed"<<std::endl;
    }
    else {
      MPI::COMM_WORLD.Send(&pid,1,MPI::INT,0,rank);
      MPI::COMM_WORLD.Send(host,MPI_MAX_PROCESSOR_NAME,MPI::CHAR,0,rank);
      int nrecv;
      MPI::COMM_WORLD.Recv(&nrecv,1,MPI::INT,0,size+rank);
      std::vector<int> recv(nrecv);
      MPI::COMM_WORLD.Recv(&recv.front(),nrecv,MPI::INT,0,size+rank);
      SetMPIRecv(recv);
      msg_MPIDebugging() << " rank = " << rank << "\n";
      for(int aa=0;aa<recv.size();++aa)
         msg_MPIDebugging() << "    recv[" << aa << "] = " << recv[aa] << "\n";
    }
  }
#endif
}

void My_MPI::SetMPIRecv(std::vector<int> r)
{
#ifdef USING__MPI
  int rank=MPI::COMM_WORLD.Get_rank();
  if (rank==0) { // rank 0 controls everything
    m_hasrecv=true;
    m_recv=MPI::COMM_WORLD.Split(rank,rank);
    m_send=MPI::COMM_WORLD.Split(MPI_UNDEFINED,rank);
  }
  else {
    if (r[0]==0) { // rank 0 of each node (or rank 1 in the case of the host that shares rank 0)
      m_hassend=m_hasrecv=true;
      m_send=MPI::COMM_WORLD.Split(r[0],rank); // create a subcomm of all ranks on one host
      m_recv=MPI::COMM_WORLD.Split(rank,rank); // create a subcomm of all just this rank
    }
    else { // all other ranks
      m_hassend=true; 
      m_recv=MPI::COMM_WORLD.Split(MPI_UNDEFINED,rank); // set m_recv to be a NULL COMM
      m_send=MPI::COMM_WORLD.Split(r[0],rank); // create a subcomm of all ranks on one host (same as above)
    }
  }
#endif
}

bool My_MPI::HasMPISend() const
{
  return m_hassend;
}

bool My_MPI::HasMPIRecv() const
{
#ifdef USING__MPI
  if (m_hasrecv) return m_recv.Get_size()>1;
#endif
  return false;
}


#ifdef USING__MPI
  int My_MPI::num_bcast_string = 0;
  int My_MPI::num_bcast = 0;
  int My_MPI::num_bcast_vect = 0;
  int My_MPI::num_bcast_matrix = 0;
  template <>
    void My_MPI::Bcast<std::string>(std::string& val){
       num_bcast_string++;
       int string_size =  val.size();
       std::stringstream ss;
       ss << "type=std::string calls=" << num_bcast_string << " starting size=" << string_size << "string=" << val;
       My_MPI::PrintMessage(__PRETTY_FUNCTION__,ss.str());
       MPI::COMM_WORLD.Bcast(&string_size,1,MPI::INT,0);
       ss.str("");
       ss << "type=std::string calls=" << num_bcast_string << " middle size=" << string_size << "string=" << val;
       My_MPI::PrintMessage(__PRETTY_FUNCTION__,ss.str());
       char buffer[string_size+1];
       strcpy(buffer,val.c_str());
       MPI::COMM_WORLD.Bcast(buffer,string_size,MPI::CHAR,0);
       buffer[string_size] = '\0';
       val = std::string(buffer);
       ss.str("");
       ss << "type=std::string calls=" << num_bcast_string << " ending size=" << val.size() << "string=" << val;
       My_MPI::PrintMessage(__PRETTY_FUNCTION__,ss.str());
     }
  void My_MPI::PrintMessage(std::string class_method,std::string msg){
    msg_MPIDebugging() << "MPI=-=" << GetRank() << "=-=" << class_method << "=-=" << MyTiming::GetCurrentTimeString() << "=-=" << msg << "\n";
  }
  int  My_MPI::GetRank(){
    return MPI::COMM_WORLD.Get_rank();
  } 
#endif

