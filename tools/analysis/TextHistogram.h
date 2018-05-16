#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

class Bin1D{
public:
   Bin1D(double low_edge, double content = 0.0){
      m_low_edge = low_edge;
      m_content = content;
   }
   double content(){return m_content;}
   void increment(){m_content++;}
   void reset(){m_content=0.0}
   double low_edge(){return m_low_edge;}
private:
   double m_low_edge;
   double m_content;
};

class TextHistogram1D{
public:
   TextHistogram1D(double bin_width,double bin_modulus = 0.0)
   {
      m_bin_width       = bin_width;
      m_bin_modulus     = bin_modulus;
      m_lowest_bin_edge = 0.0;
      m_bins            = new std::vector<Bin1D>;
   }
   ~TextHistogram1D(){
      delete m_bins; m_bins = 0;
   }

   int Fill(double value){
      // loop over the vector to find the proper bin
      for(std::vector<double>::iterator it = m_bin_content->begin();
          it != m_bin_content->end();++it){

      }
   }


private:
   double m_bin_width;
   double m_bin_modulus;
   double m_lowest_bin_edge;
   std::vector<Bin1D>* m_bins;

   int AddBin()
};