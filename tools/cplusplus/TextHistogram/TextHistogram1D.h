#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include <exception>
#include <iomanip>

class Bin1D{
public:
   Bin1D(double bin_low_edge, double bin_width,double bin_content = 0.){
      m_bin_low_edge = bin_low_edge;
      m_bin_width = bin_width;
      m_bin_content = bin_content;
   }


   static Bin1D BuildBin(const double value,const double bin_intercept,const double bin_width){
      double bin_low_edge = int((value - bin_intercept)/bin_width) * bin_width + bin_intercept;
      // std::cout << "BuildBin: " << value << " " << bin_intercept << " " << bin_width << " " << bin_low_edge << "\n";
      return Bin1D(bin_low_edge,bin_width);
   }

   void Fill(double weight = 1.0){
      m_bin_content += weight;
   }
   double bin_content(){return m_bin_content;}
   double bin_low_edge(){return m_bin_low_edge;}
   double bin_high_edge(){return m_bin_low_edge + m_bin_width;}
   double bin_width(){return m_bin_width;}
   bool contains(double value){
      if(m_bin_low_edge <= value && value < (m_bin_low_edge + m_bin_width))
         return true;
      return false;
   }
private:
   double m_bin_content;
   double m_bin_low_edge;
   double m_bin_width;
};

class TextHistogram1D{
public:
   TextHistogram1D(double bin_width,double bin_intercept):m_bins(*(new std::vector<Bin1D>)){
      m_bin_width = bin_width;
      m_bin_intercept = bin_intercept;
   }
   ~TextHistogram1D(){delete &m_bins;}

   void Fill(double value,double weight = 1.0){
      // Find bin to fill
      unsigned int bin_number = 0;
      try{
         // std::cout << " Calling FindBin: " << value << "\n";
         bin_number = FindBin(value);
      } // no bin found so let's add one
      catch(std::string message){
         try{
            // std::cout << " Calling AddBin: " << value << "\n";
            bin_number = AddBin(value);
         } // failed to add bin
         catch(std::string message){
            std::cerr << "ERROR Adding bin for value: " << value << "\n";
            return;
         }
      }
      // std::cout << " After Find/AddBin, Filling: " << bin_number << " " << value << "\n";
      m_bins[bin_number].Fill(weight);
   }

   unsigned int FindBin(double value){
      // std::cout << "In FindBin: " << value << "\n";
      for(std::vector<Bin1D>::iterator it = m_bins.begin();
          it != m_bins.end();++it){
         
         // if the value falls within the boundaries of this bin return the bin number
         if( it->contains(value) ){
            unsigned int bin_number = (unsigned int)(it - m_bins.begin());
            // std::cout << " Found Bin: " << value << " " << it->bin_low_edge() << " " << it->bin_width() << " " << bin_number << "\n";
            return bin_number;
         }
      }

      std::stringstream ss;
      ss << "Did not find a bin for value: " << value;
      throw ss.str();
      return 0;
   }

   unsigned int AddBin(double value){
      // std::cout << " In AddBin: " << value << "\n";
      // find the place in the vector to insert a bin
      unsigned int new_bin_number = 0;
      // if there are no bins, add the first one
      if(m_bins.size() == 0){
         // std::cout << " This is the first bin created.\n";
         m_bins.push_back(Bin1D::BuildBin(value,m_bin_intercept,m_bin_width));
         return 0;
      }
      // this should catch all other cases
      else{
         //std::cout << " Looping over bins\n";
         for(std::vector<Bin1D>::iterator it = m_bins.begin();
             it != m_bins.end();++it){
            // if the value falls between the previous bin and the current bin we need to add a new bin in between
            if( it != m_bins.begin() ){
               if( (it-1)->bin_high_edge() <= value && value < it->bin_low_edge()){
                  //std::cout << " Found the location at which to add the bin.\n";
                  // add bin
                  m_bins.insert(it,Bin1D::BuildBin(value,m_bin_intercept,m_bin_width));
                  return it - m_bins.begin();
               }
            }
            // if this is the first bin, check if value is below it.
            else{
               //std::cout << " checking if bin is below the first bin.\n";
               if(value < it->bin_low_edge()){
                  m_bins.insert(it,Bin1D::BuildBin(value,m_bin_intercept,m_bin_width));
                  return it - m_bins.begin();
               }
            }
         }
         // if the value falls after the last bin, check for that
         if(m_bins.back().bin_high_edge() < value){
            //std::cout << " This value is above the last bin.\n";
            m_bins.push_back(Bin1D::BuildBin(value,m_bin_intercept,m_bin_width));
            return m_bins.size() - 1;
         }
      }
      std::stringstream ss;
      ss << "Did not find a bin for value: " << value;
      throw ss.str();
      return 0;
   }

   std::string str(){

      std::stringstream ss;

      ss << "|              Histogram Contents             |\n";
      ss << " -------------- --------------- --------------\n";
      ss << "| Bin Low Edge | Bin High Edge |  bin content |\n";
      for(std::vector<Bin1D>::iterator it = m_bins.begin();it != m_bins.end();++it)
         ss << "| " << std::setw(12) << it->bin_low_edge() << " | " << std::setw(13) << it->bin_high_edge() << " | " << std::setw(12) << it->bin_content() << " |\n";
      ss << " -------------- --------------- --------------\n";
      return ss.str();
   }

private:
   double m_bin_width;
   double m_bin_intercept;
   std::vector<Bin1D>& m_bins;

};