# ROOT Stuff
ROOTC = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)

COPTS_OPTIONAL= -I/home/projects/ATLASQ/tools/zlib-1.2.8
LOPTS_OPTIONAL=-lz -L/home/projects/ATLASQ/tools/zlib-1.2.8 

HEPMCPATH=/users/hpcusers/svn/tools/HepMC-2.06.08
HEPMC_COPTS=-I $(HEPMCPATH)
HEPMC_LOPTS=-L$(HEPMCPATH)/local/lib -lHepMC -lHepMCfio

PYTHIAPATH=/users/hpcusers/svn/generators/pythia/v8180/bgq/trunk
PYTHIA_COPTS=-I $(PYTHIAPATH)/include -I $(PYTHIAPATH)/examples
PYTHIA_LOPTS=-L $(PYTHIAPATH)/lib/archive -lpythia8 -lpythia8tohepmc

LHAPDFPATH=/users/hpcusers/svn/tools/lhapdf-5.9.1/local
LHAPDF_LOPTS=-L $(LHAPDFPATH)/lib -lLHAPDF -lgfortran

TCLAPPATH=/users/hpcusers/svn/tools/tclap-1.2.1
TCLAP_COPTS=-I $(TCLAPPATH)/include

COPTS=-c $(COPTS_OPTIONAL) $(PYTHIA_COPTS) $(HEPMC_COPTS) $(TCLAP_COPTS)
LOPTS=$(PYTHIA_LOPTS) $(HEPMC_LOPTS) $(LHAPDF_LOPTS) $(LOPTS_OPTIONAL)

GCC=mpic++

SOURCES=Gzip_Stream.cpp
OBJECTS=$(SOURCES:.cpp=.o)
BINARY=pythia

target: $(BINARY) $(SOURCES)

# Link binary executables
$(BINARY): $(BINARY).o $(OBJECTS)
	@echo " "
	@echo "--------->>>>>> Linking" $<
	$(GCC) $(OBJECTS) $(BINARY).o $(LOPTS) $(LOPTS) -o $@
	@echo "----->>> Done."


# Build binary Source
$(BINARY).o: $(BINARY).cxx
	@echo " "
	@echo "--------->>>>>> Compiling" $*.o
	$(GCC)  $(COPTS) $< -o $@
	@echo "----->>> Done."

# Build Generic Object Files
.cpp.o:
	@echo " "
	@echo "--------->>>>>> Compiling" $*.o
	$(GCC) $(COPTS) $< -o $@
	@echo "----->>> Done."


link:
	ln -s $(PYTHIAPATH)/xmldoc

clean:
	@echo " "
	@echo "--------->>>>>> Cleaning"
	rm -f *.o $(BINARY)
	@echo "----->>> Done."

