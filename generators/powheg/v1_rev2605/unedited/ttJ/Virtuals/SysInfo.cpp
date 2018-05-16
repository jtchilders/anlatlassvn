#include "SysInfo.h"
#include <cstdlib>
#include <iostream>
#include <sys/resource.h>

#define LINE "# -------------------------------------------------------------"

using namespace std;

SysInfo::SysInfo(){

  time(&start);
  nicelevel = getpriority(PRIO_PROCESS, 0);
  pid = getpid();
  
  cout << "#" << endl;
  cout << LINE << endl;
  cout << "#    Started: " << ctime(&start);
  cout << "#        Pid: " << pid << endl;
  cout << "# Nice level: " << nicelevel << endl; 
  
  jobid = getenv("JOB_ID");
  if ( jobid != 0 )
    cout << "#\n# qsub job id: " << string(jobid) << "\n#" << endl;
  

  int retrycounter = 5;
  int err;
  do {
    err = ::gethostname(hostname, HOSTNAME_LENGTH);
    if ( 0 == err ) {
      cout << "# Running on: " << hostname << std::endl;
    } else {
      cout << "# Cannot retrieve hostname, error = " 
	   << err << endl;
      cout << "# retry..."  << endl;
      sleep(1);
      retrycounter--;
    }
  } while ( ( retrycounter > 0 ) && ( err != 0 ) ); 
  cout << LINE << endl;
  cout << "#" << endl;
}

SysInfo::~SysInfo(){
  time_t end;
  time(&end);
  time_t elapsed = end - start;
  
  cout << LINE << endl;
  cout << LINE << endl;
  cout << "# Finished: " << ctime(&end);
  cout << "# Total time elapsed " << elapsed
       << " seconds" << endl << "#" << endl; 
  
  if ( 0 == getrusage(0,&rusage)  ) {
    
    double utime = rusage.ru_utime.tv_sec
      + rusage.ru_utime.tv_usec/1000000.0;
    double stime = rusage.ru_stime.tv_sec
      +rusage.ru_stime.tv_usec/1000000.0 ;
    
    cout << "# Page faults (soft): " << rusage.ru_minflt << endl;
    cout << "# Page faults (hard): " << rusage.ru_majflt << endl;
    cout << "# Context switches  : " << rusage.ru_nivcsw << endl
	 << "#" << endl;
    
    cout << "# User time used since program start: " 
	 << utime << " seconds " << endl;
    cout << "# System time used since program start: " 
	 << stime << " seconds "<< endl;  
    
    if ( elapsed > 0 ) 
      cout << "# CPU fraction = " << utime/elapsed*100. << " %" << endl;
    
  } else {
    cout << "# Could not get information from getrusage." << endl;
  }
  cout << "# ...done." << endl;
}

void SysInfo::renice(){
  if ( 0 == setpriority(PRIO_PROCESS, 0, 19) ) {
    cout << "# Changed nice level to 19 \n";
  } else {
    cout << "# Failed to change nice level \n";
  }
  nicelevel = getpriority(PRIO_PROCESS, 0);
  cout << "# Current nice level: " << nicelevel << endl << "#\n"; 
}

string SysInfo::gethostname(){
  return(string(hostname));
}

string SysInfo::getjobid(){
  if ( jobid != 0 ) 
    return(string(jobid));
  return(string("not available"));
}

#ifdef MAIN

int main(){
  
  SysInfo sysinfo;
}

#endif
