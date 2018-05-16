FC = gfortran
FFLAGS = -g -fno-automatic -ffixed-line-length-none -O2

# ----------------------------------------------------------------------------

# object files for qq->qqZ with leptonic decays

OBJECTS =  ./li2.o \
	  ./help.o \
	  ./braket_vbf.o \
	  ./init_tensors_z.o \
          ./qqzqq.o \
          ./qqzqqj.o \
          ./dhelas/ixxxxx.o \
          ./dhelas/oxxxxx.o \
          ./dhelas/fvixxx.o \
          ./dhelas/fvoxxx.o \
          ./dhelas/jvvxxx.o \
          ./dhelas/vvvxxx.o \
          ./dhelas/jioxxx.o \
          ./dhelas/iovxxx.o \
          ./dhelas/hvvxxx.o \
          ./dhelas/jvsxxx.o \
          ./dhelas/w3w3xx.o \
          ./dhelas/wwwwxx.o \
          ./dhelas/vvsxxx.o \
          ./dhelas/vxxxxx.o \
          ./dhelas/jwwwxx.o \
          ./dhelas/jw3wxx.o \
          ./toll-hel.o \
	  ./boxline-coeffs-zjj-full.o \
	  ./boxline-coeffs-cross-zjj-full.o \
	  ./tensors/c0-int-mu.o \
	  ./tensors/d0-int-mu.o \
	  ./tensors/scalar-dint-order.o \
	  ./tensors/tensor-c-full.o \
	  ./tensors/tensor-d-quarks.o   

TARGETS = libnlo.a


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# rules for nlo_lib:
#
 libnlo.a:  $(OBJECTS)
	$(AR) r libnlo.a $(OBJECTS) 

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

# for g77 or gfortran: 
.f.o:
	$(FC) -c -o $*.o  $(FFLAGS) $*.f

# ----------------------------------------------------------------------------

clean:
	rm -f *.o core *~ *.x libnlo.a



