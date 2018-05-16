
COMPILER=gfortran
FC=$(COMPILER)
ifeq ("$(FC)","gfortran")
F77 = gfortran
FFLAGS = -fno-automatic -ffixed-line-length-none -O2
endif
ifeq ("$(FC)","g77")
F77= g77 -fno-automatic -ffixed-line-length-none -O2
endif
ifeq ("$(FC)","ifort")
F77= ifort -save -132 -O3
endif


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
#
# object files for qq->qq(g)W+W+ to be implemented in POWHEG;
OBJECTS_PH = ./loops/ggf_loops/loops_aux.o \
          ./vvjj_amp_aux.o \
          ./anomHiggs.o \
          ./anomV.o \
          ./anomal.o\
          ./brakets.o \
	  ./qqZWjqqi.o \
	  ./qqbqq.o \
	  ./init_tensors.o \
          ./help.o \
          ./kinematics_aux.o \
          ./tbv.o \
          ./jvirt5.o \
          ./jvirt5-isig.o \
          ./HtoWW.o \
          ./dotproducts.o \
          ./parameters.o \
          ./helas/ixxxxx.o \
          ./helas/oxxxxx.o \
          ./helas/fvixxx.o \
          ./helas/fvoxxx.o \
          ./helas/jvvxxx.o \
          ./helas/vvvxxx.o \
          ./helas/jioxxx.o \
          ./helas/iovxxx.o \
          ./helas/hvvxxx.o \
          ./helas/jvsxxx.o \
          ./helas/w3w3xx.o \
          ./helas/wwwwxx.o \
          ./helas/vvsxxx.o \
          ./helas/vxxxxx.o \
          ./helas/jwwwxx.o \
	  ./loops/aux_loops.o \
              ./loops/Integrals/A0.o \
              ./loops/Integrals/B0C0D01MASG.o \
              ./loops/Integrals/C0functions.o  \
              ./loops/Integrals/D0functionsDiv.o \
              ./loops/Integrals/dilog.o        \
              ./loops/Integrals/F0.o  \
               ./loops/Integrals/E0functions.o  \
              ./loops/Integrals/scalarints.o\
              ./loops/Integrals/Aux_Functions.o  \
              ./loops/Integrals/B0functions.o  \
              ./loops/Integrals/cdilog.o       \
              ./loops/Integrals/D0functions.o \
              ./loops/Utilities/bra_ket_subroutines.o  ./loops/Utilities/GaugeTest.o  ./loops/Utilities/MatrixElements.o  \
              ./loops/Utilities/pent45_noAbe.o  ./loops/Utilities/Pent_to_boxes.o  ./loops/Utilities/Wardtest.o
# 	          ./helas/anomal6.o\
#           ./helas/anomal8.o\  
TARGETS_PH = libvbfnlo.a

# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# rules for vbfnlo_lib:
#

 libvbfnlo.a:  $(OBJECTS_PH)
	$(AR) r libvbfnlo.a $(OBJECTS_PH) 

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

# for g77 or gfortran: 
.o.o:
	$(F77)  -c -o $*.o  $(FFLAGS) $*.o

# ----------------------------------------------------------------------------

clean:
	rm -f *.o */*.o */*/*.o */*/*/*.o core *~ *.x


## 	  ./tensor_redux/B0C0D0fin.o \
