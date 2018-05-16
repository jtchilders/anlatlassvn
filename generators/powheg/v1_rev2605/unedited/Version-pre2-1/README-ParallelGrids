In several case, the formation of the importance sampling grid
has proven to be a bottelneck in POWHEG event generation. This
is because, as it stands, it cannot be parallelized. The files
in this directory provide a way to parallelize also the grid
generation. In order to use them, edit the Makefile of the
process you need, and modify the line

VPATH= ./:../:obj/

into

VPATH= ../Parallel-grids/:./:../:obj/

so that the file in this directory override the default ones.

In order to drive the parallel grid calculation, refer to the
example files in this directory:

run-with-parallel-grids
powheg.input-save

that are ment for running on a multi-core computer.
