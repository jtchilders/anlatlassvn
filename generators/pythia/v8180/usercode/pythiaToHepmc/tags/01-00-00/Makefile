# ROOT Stuff
ROOTC = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)

COPTS_OPTIONAL=
LOPTS_OPTIONAL=-static

HEPMCPATH=/users/hpcusers/svn/tools/HepMC-2.06.08
HEPMC_COPTS=-I $(HEPMCPATH)
HEPMC_LOPTS=$(HEPMCPATH)/src/.libs/*.o
#-lHepMC -L$(HEPMCPATH)/lib
#$(HEPMCPATH)/src/.libs/*.o

PYTHIAPATH=/users/hpcusers/svn/generators/pythia/v8180/bgq/trunk/
PYTHIA_COPTS=-I $(PYTHIAPATH)/include -I $(PYTHIAPATH)/examples 
PYTHIA_LOPTS= $(PYTHIAPATH)/tmp/*.o
#-L$(PYTHIAPATH)/lib -lpythia8 -llhapdfdummy -lpythia8tohepmc

TCLAPPATH=/users/hpcusers/svn/tools/tclap-1.2.1
TCLAP_COPTS=-I $(TCLAPPATH)/include

OPENMPI_INCPATH=/usr/include/openmpi-x86_64/
OPENMPI_COPTS=-I $(OPENMPI_INCPATH)
OPENMPI_LIBPATH=/usr/lib64/openmpi/lib/
OPENMPI_LOPTS=-L $(OPENMPI_LIBPATH) -lmpi

COPTS=-c $(COPTS_OPTIONAL) $(PYTHIA_COPTS) $(HEPMC_COPTS) $(TCLAP_COPTS)
LOPTS=$(PYTHIA_LOPTS) $(HEPMC_LOPTS) $(LOPTS_OPTIONAL) $(OPENMPI_LOPTS)

GCC=mpic++

BINARY=pythia

target: $(BINARY)

# Link binary executables
$(BINARY): $(BINARY).o
	@echo " "
	@echo "--------->>>>>> Linking" $<
	$(GCC) $(LOPTS) $< -o $@
	@echo "----->>> Done."


# Build binary Source
$(BINARY).o: $(BINARY).cpp
	@echo " "
	@echo "--------->>>>>> Compiling" $*.o
	$(GCC) -MD $(COPTS) $< -o $@
	@echo "----->>> Done."

clean:
	@echo " "
	@echo "--------->>>>>> Cleaning"
	rm -f $(BINARY).o $(BINARY)
	@echo "----->>> Done."

