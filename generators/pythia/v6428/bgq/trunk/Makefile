

SRC=pythia-6.4.28.f
OBJ=$(SRC:.f=.o)
GFF=gfortran

target: ${OBJ}

${OBJ}: ${SRC}
	@echo ' >>>> compiling ' ${SRC}
	${GFF} -cpp -c ${SRC}

