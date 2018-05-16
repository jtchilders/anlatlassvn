// To compile: 
//  g++ -o alproot alproot.cc \
// -I/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/5.34.14-x86_64-slc6-gcc4.7/include

#define MAXJETS 6

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

#include "TFile.h"
#include "TNtuple.h"
#include "TVector3.h"

void rootify(float x[], TNtuple *h1){

   bool unswapped;
   int i,j;
   float y[3*MAXJETS+12];
   float t0,t1,t2;

   // Processes a flattened version of the unweighted files

   // Process the leptons and jets
   for(i=0;i<MAXJETS+2; i++) {
      TVector3 *v = new TVector3(x[3*i],x[3*i+1],x[3*i+2]);
      y[3*i]=v->Perp();
      y[3*i+1]=v->Phi();
      y[3*i+2]=v->Eta();
   }

   // Sort the jets by pT using an inefficient algorithm

   unswapped=true;
   while(unswapped) {
      unswapped = false;
      for(i=2;i<MAXJETS+1; i++) {
         for(j=i+1;j<MAXJETS+2;j++) { 
            if(y[3*j] > y[3*i] ) {

               t0 = y[3*i];
               t1 = y[3*i+1];
               t2 = y[3*i+2];

               y[3*i]   = y[3*j];
               y[3*i+1] = y[3*j+1];
               y[3*i+2] = y[3*j+2];

               y[3*j]   = t0;
               y[3*j+1] = t1;
               y[3*j+2] = t2;

            } // if
         } // for
      }  // for
   } // while

   // Process the other ntuple values
   for(i=3*MAXJETS+6;i<3*MAXJETS+12;i++) {
      y[i]=x[i];
   }

   h1->Fill(y);

}

void alproot() {

   // Flattens the unweighted files

   int i;

   char *inputstring;
   int bytes_read;
   size_t max_line_length = 132;

   float outdata[3*MAXJETS+12];

   unsigned int eventnum;
   unsigned int processId;
   unsigned int nparticles;
   float weight;
   float scaleInGeV;

   int particle_code;
   int color1;
   int color2;
   float px,py,pz,mass;

   double e,p;


   // Initalization

   inputstring = (char *) malloc (max_line_length + 1);
   for(i=0; i< 3*MAXJETS+12; i++) outdata[i]=0.001;

   std::ifstream inputFile("alpout.unw");
   if(inputFile.fail()) {
      fprintf(stderr,"Input file does not exist.\n");
      exit(1);
   }

   TFile f("alpout.root","recreate");

   // Booking (this still requires adjustments when changing MAXJETS)
   TNtuple *h1 = new TNtuple("Alpgen","Alpgen",
   "l1_pt:l1_phi:l1_eta:l2_pt:l2_phi:l2_eta:j1_pt:j1_phi:j1_eta:j2_pt:j2_phi:j2_eta:j3_pt:j3_phi:j3_eta:j4_pt:j4_phi:j4_eta:j5_pt:j5_phi:j5_eta:j6_pt:j6_phi:j6_eta:ll_mass:processId:weight:scaleInGeV:color1:color2");

   // Process the file

   while(inputFile.getline(inputstring, (std::streamsize) 132)) {

      // Jets

      sscanf(inputstring,"%d %d %d %f %f",
         &eventnum, &processId, &nparticles, &weight, &scaleInGeV);

      inputFile.getline(inputstring, (std::streamsize) 132);
      inputFile.getline(inputstring, (std::streamsize) 132);
      for(i=0; i < nparticles-4; i++) {  // Jets

         inputFile.getline(inputstring, (std::streamsize) 132);
         sscanf(inputstring,"%d %d %d %f %f %f %f",
            &particle_code,&color1,&color2,
            &px,&py,&pz,&mass);
         if(mass!=0.) {
            fprintf(stderr,"Jet-lepton confusion.  Exiting.\n");
            exit(1);
         }
         outdata[6+3*i]=px;
         outdata[7+3*i]=py;
         outdata[8+3*i]=pz;
      }


      // Leptons

      inputFile.getline(inputstring, (std::streamsize) 132);   
      sscanf(inputstring,"%d %d %d %f %f %f %f",
         &particle_code,&color1,&color2,
         &px,&py,&pz,&mass);
      if(mass==0.) {
         fprintf(stderr,"Lepton-jet confusion.  Exiting.\n");
         exit(1);
      }
      outdata[0]=px;
      outdata[1]=py;
      outdata[2]=pz;

      inputFile.getline(inputstring, (std::streamsize) 132);             
      sscanf(inputstring,"%d %d %d %f %f %f %f",
         &particle_code,&color1,&color2,
         &px,&py,&pz,&mass);
      if(mass==0.) {
         fprintf(stderr,"Lepton-jet confusion.  Exiting.\n");
         exit(1);
      }
      outdata[3]=px;
      outdata[4]=py;
      outdata[5]=pz;
 
      // Calculated Quantities

      e = sqrt(pow(mass,2) + pow(outdata[0],2) + pow(outdata[1],2) + pow(outdata[2],2)) +
          sqrt(pow(mass,2) + pow(outdata[3],2) + pow(outdata[4],2) + pow(outdata[5],2));
      p = sqrt( pow(outdata[0] + outdata[3],2) +
                pow(outdata[1] + outdata[4],2) +
                pow(outdata[2] + outdata[5],2) );
      outdata[3*MAXJETS+6]=sqrt(e*e-p*p);
      outdata[3*MAXJETS+7]=processId;
      outdata[3*MAXJETS+8]=weight;
      outdata[3*MAXJETS+9]=scaleInGeV;
      outdata[3*MAXJETS+10]=color1;
      outdata[3*MAXJETS+11]=color2;

      rootify(outdata, h1);

   } // while

   // Finalize and clean up

   inputFile.close();
   h1->Write();
   f.Close();

}


// Magic for running inside ROOT
 
# ifndef __CINT__
int main()
{
  alproot();
  return 0;
}
# endif
