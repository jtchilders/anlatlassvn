#include "TextHistogram1D.h"


int main(){

   TextHistogram1D h(1.,0.);

   std::cout << "FILL: 1\n";
   h.Fill(1.);
   std::cout << h.str();
   std::cout << "FILL 1.5\n";
   h.Fill(1.5);
   std::cout << h.str();
   std::cout << "FILL 3\n";
   h.Fill(3.);
   std::cout << h.str();
   std::cout << "FILL 9\n";
   h.Fill(9.);
   std::cout << h.str();
   std::cout << "FILL 15.0001\n";
   h.Fill(15.0001);
   std::cout << h.str();

   std::cout << h.str();


   return 0;
}