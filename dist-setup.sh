## Create workspace
## set environment
export ROOT=~/build/fst-gl/static-build

mkdir -p $ROOT
cd $ROOT


## Build libdeflate
cd $ROOT
git clone https://github.com/ebiggers/libdeflate.git
cd libdeflate 
git checkout v1.25 
cmake -B build -DCMAKE_INSTALL_PREFIX=$ROOT
cmake --build build
cmake --install build


## build hts-lib
cd $ROOT
git clone https://github.com/samtools/htslib.git
cd htslib
git checkout 1.23.1   # stable version recommended
git submodule update --init --recursive
autoreconf -i
./configure --prefix=$ROOT \
            --disable-libcurl \
            --with-libdeflate \
            CPPFLAGS="-I$ROOT/include" \
            LDFLAGS="-L$ROOT/lib64"
make && make install


cd $ROOT
export NIM_PREFIX=/opt/nim
git clone https://github.com/nim-lang/Nim.git
cd Nim
git checkout v2.2.4   # or whatever version you want

sh build_all.sh
bin/nim c koch
./koch boot -d:release
./koch install /opt

mkdir -p $NIM_PREFIX
rsync -av bin/ $NIM_PREFIX/bin/
rsync -av lib/ $NIM_PREFIX/lib/

export PATH=$NIM_PREFIX/bin:$PATH

export NIM_PATH=$NIM_PREFIX/lib
export NIM_LIB_PREFIX=$NIM_PREFIX/lib
export NIMBLE_PATH=$NIM_PREFIX

cd $ROOT

export LD_LIBRARY_PATH=$ROOT/lib:$ROOT/lib64:$LD_LIBRARY_PATH