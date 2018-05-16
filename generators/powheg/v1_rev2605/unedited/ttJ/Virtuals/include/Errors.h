/* $Modified: Wed Aug  5 13:58:55 2009 by uwer $ */
#ifndef ERRORS_H_
#define ERRORS_H_
#include <cstring>
//#define VERBOSE

class ErrorCounter {
private: 
  int count;
  int mod;
  std::string message;
public:
  ErrorCounter(std::string s) : message(s), count(0), mod(10000) {}
    ErrorCounter(std::string s, int hlp) : message(s), count(0), mod(hlp) {}
  ~ErrorCounter(){
#ifdef VERBOSE
    std::cout << "# WARNING: ErrorCounter: The message: " 
	      << message << std::endl
	      << "# was produced " << count << " times." << std::endl;
#endif
  }
  void addError(){
    count++;
#ifdef VERBOSE
    if ( (count % mod) == 0 ){
      std::cout << "# WARNING: " << message 
		<< ", ErrorCount = " << count << std::endl;
    }
#endif
  }
  inline int getErrorCount(){
    return(count);
  }
};

#endif // ERRORS_H_
