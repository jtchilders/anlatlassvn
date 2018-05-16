FC = gfortran
INCLUDE1=$(PWD)
INCLUDE2=$(PWD)/..
INCLUDE3=$(PWD)/../boxfiles-pre2-1
FFLAGS = -g -fno-automatic -ffixed-line-length-none -O2 -I$(INCLUDE1) -I$(INCLUDE2) 

# ----------------------------------------------------------------------------

# object files for qq->qqW+W- with leptonic decays

OBJECTS =  ./amp_aux.o \
	  ./help.o \
	  ./braket_vbf.o \
	  ./init_tensors_wpm.o \
	  ./provide_tensors_wpm.o \
          ./qqwwqq.o \
          ./qqwwqqj.o \
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
          ./dhelas/hioxxx.o \
          ./dhelas/fsixxx.o \
          ./dhelas/fsoxxx.o \
          ./toww.o \
          ./toww_slp.o \
          ./toww_slm.o \
          ./jvirt5.o \
          ./tensor_redux/tri_box_penta_tens_red.o \
	  ./tensor_redux/bra_ket_subroutines_vbf.o\
	  ./tensor_redux/B0C0D0fin.o \
          ./tensor_redux/tensor_redux_Cij.o \
          ./tensor_redux/tensor_redux_Dij.o \
          ./tensor_redux/tensor_redux_Eij.o \
	  ./tensor_redux/dilog.o \
          ./tensor_redux/virtual1.o \
	  ./tensor_redux/virtual2.o  

TARGETS = libvbfnlo.a


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# rules for vbfnlo_lib:
#
 libvbfnlo.a:  $(OBJECTS)
	$(AR) r libvbfnlo.a $(OBJECTS) 

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

# for g77 or gfortran: 
.f.o:
	$(FC) -c -o $*.o  $(FFLAGS) $*.f

# ----------------------------------------------------------------------------

clean:
	rm -f *.o core *~ *.x libvbfnlo.a



