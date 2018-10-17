##Build Macaulay2 1.12.0.1

module load compilers/gcc/5.3
module load libs/gcc/5.2/boost/1.59
module load libs/gcc/lapack/3.3.0
module load apps/xzutils/5.2.4/gcc-5.3
module load apps/yasm/1.3.0/gcc-5.3
module load apps/make/3.82/gcc-5.3
module load apps/gcc/5.2/git/2.5

workers=4

git clone https://github.com/Macaulay2/M2
cd M2
git checkout release 1.12.0.1

cd M2

make -f Makefile get-tools
export CPPFLAGS="$CPPFLAGS $CFLAGS"
export CFLAGS= 
make
./configure --enable-download --enable-ntl-wizard --prefix=/usr/local/packages6/apps/macaulay2/1.12.0.1/gcc-5.3
make -j${workers}
make check
make install

##Notes
##cohomcalg fails to build due to -g3 o2 appearing in CFLAGS, which is not recognised
##edit libraries/cohomcalg/makefile.in
##set CFLAGS='-std=gnu11'
##build completes