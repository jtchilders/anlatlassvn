/* $Modified: Fri Oct  5 13:52:05 2007 by uwer $ */
#ifndef SYSINFO_H_
#define SYSINFO_H_
#include <time.h>  
#include <string> 
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>
#include <iostream>
  

class SysInfo {
  enum { HOSTNAME_LENGTH = 150};
  int pid;
  char *jobid;
  int nicelevel;
  time_t start;
  char hostname[HOSTNAME_LENGTH];
  struct rusage rusage;
 public:
  SysInfo();
  ~SysInfo();
  void renice();
  std::string gethostname();
  std::string getjobid();
};
#endif
