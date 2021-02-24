BASE_DIR=$(pwd)
CC = ../mpich/bin/mpicxx
HEADER = evp.h evpFlags.h V3math.h pf_xrx.h PFgrid.h
CFLAGS = -Wall -O3
SEEKFLAG = -DMPICH_IGNORE_CXX_SEEK
MPI_CFLAGS = -I../mpich/include
MPI_LIBS = -L../mpich/lib
FFTW_CFLAGS = -I../fftw/include
FFTW_LIBS = -L../fftw/lib
LIBS = -lfftw_mpi -lfftw -lgsl -lgslcblas -lm -lhdf5
SRCS = evp.c io.c init.c kinematics.c evolution.c constitutive.c integrated.c pf_xrx.c  
OBJS = $(SRCS:.c=.o)
MAIN = forge
.PHONY: depend clean
#all: $(MAIN)

$(MAIN): $(OBJS) Makefile
	$(CC) $(SEEKFLAG) $(FFTW_CFLAGS) $(MPI_CFLAGS)  $(CFLAGS)  -o $(MAIN) $(OBJS) $(FFTW_LIBS) $(MPI_LIBS) $(LIBS)

%.o: %.c $(HEADER) Makefile
	$(CC) $(SEEKFLAG) $(FFTW_CFLAGS) $(MPI_CFLAGS) $(CFLAGS)  -c $< -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $^
