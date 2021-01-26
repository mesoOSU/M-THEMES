BASE_DIR=$(pwd)
cd $BASE_DIR/mpich/src
rm -r mpich-3.2.1
tar xf mpich-3.2.1.tar.gz
cd $BASE_DIR/mpich/src/mpich-3.2.1
./configure --prefix=$BASE_DIR/mpich  --disable-fortran
make
make install
cd $BASE_DIR/fftw/src
rm -r fftw-2.1.5
tar xf fftw-2.1.5.tar.gz
cd $BASE_DIR/fftw/src/fftw-2.1.5
./configure --prefix=$BASE_DIR/fftw --enable-mpi MPICC=$BASE_DIR/mpich/bin/mpicc
make
make install
cd $BASE_DIR/src
rm *.o
make
cd $BASE_DIR