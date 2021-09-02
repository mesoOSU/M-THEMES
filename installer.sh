BASE_DIR=$(pwd)
# install mpich
cd $BASE_DIR/mpich/src
rm -r mpich-3.2.1
tar xf mpich-3.2.1.tar.gz
cd $BASE_DIR/mpich/src/mpich-3.2.1
./configure --prefix=$BASE_DIR/mpich
make
make install

# install fftw
cd $BASE_DIR/fftw/src
rm -r fftw-2.1.5
tar xf fftw-2.1.5.tar.gz
cd $BASE_DIR/fftw/src/fftw-2.1.5
./configure --prefix=$BASE_DIR/fftw --enable-mpi MPICC=$BASE_DIR/mpich/bin/mpicc
make
make install

# install HDF
cd $BASE_DIR/hdf/src
rm -r "hdf5-1.10.5"
tar xf "hdf5-1.10.5.tar.gz"
cd "$BASE_DIR/hdf/src/hdf5-1.10.5"
./configure --enable-parallel --disable-fortran --prefix=$BASE_DIR/hdf
make
make install

# build M-THEMES
cd $BASE_DIR/src
rm *.o
make
cd $BASE_DIR
