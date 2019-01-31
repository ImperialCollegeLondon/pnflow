
all: hypre zlib  tiff 

psInstall ?= $(abspath ..)
 $(info psInstall set to $(psInstall))

ifndef PS3rdParty
 # set PS3rdParty to this file directory
 PS3rdParty := $(abspath $(lastword $(MAKEFILE_LIST))/..)
 $(info PS3rdParty set to $(PS3rdParty))
endif
CmakeOpts=cmake  -DBUILD_SHARED_LIBS=False  -DCMAKE_INSTALL_PREFIX=$(psInstall) \
            -DHYPRE_SEQUENTIAL=true -DCMAKE_TOOLCHAIN_FILE=$(PS3rdParty)/toolchain$(OPT).txt

#hypre:
#	ls hypre-2.11.2/src/hypre/lib/libHYPRE.a || (cd hypre-2.11.2/src && ./configure --without-MPI   --without-fei  CC=$(preCC)gcc  F77=$(preCC)gfortran CPP=$(preCC)g++  CXX=$(preCC)g++  AR=$(preCC)ar\ -rcu   RANLIB=$(preCC)ranlib --host=$(psHost) ) || echo 
#	ls hypre-2.11.2/src/hypre/lib/libHYPRE.a || (cd hypre-2.11.2/src && make -j)

#hypreDbg:
#	(cp -r hypre-2.11.2 hypre-2.11.2Dbg)
#	(cd hypre-2.11.2Dbg/src && ./configure --enable-debug --without-MPI   --without-fei  CC=gcc-6  F77=gfortran CPP=g++-6  CXX=g++  AR=ar\ -rcu   RANLIB=ranlib ) 
#	(cd hypre-2.11.2Dbg/src && make)

#hypreWin:	 
#	(cp -r hypre-2.11.2 hypre-2.11.2MinGW)
#	(cd hypre-2.11.2MinGW/src && ./configure --without-MPI   --without-fei  CC=x86_64-w64-mingw32-gcc  F77=x86_64-w64-mingw32-gfortran CPP=x86_64-w64-mingw32-g++  CXX=x86_64-w64-mingw32-g++  AR=x86_64-w64-mingw32-ar\ -rcu   RANLIB=x86_64-w64-mingw32-ranlib --host x86_64-w64-mingw32 )
#	(cd hypre-2.11.2MinGW/src && make) 

psBuild=$(psInstall)/build


biltT=$(shell [ -d $(PS3rdParty)/libtiff ] && [ ! -f $(psInstall)/lib/libtiffxx.a ] && echo Tiff)
tiff:
	[ "Tiff" != "$(biltT)" ] || (mkdir -p $(psBuild)/libtif )
	[ "Tiff" != "$(biltT)" ] || (cd       $(psBuild)/libtif && $(CmakeOpts)  $(PS3rdParty)/libtiff)
	[ "Tiff" != "$(biltT)" ] || (cd       $(psBuild)/libtif && make -j install)


biltH=$(shell [ -d $(PS3rdParty)/hypre-2.11.2/src ] && [ ! -f $(psInstall)/lib/libHYPRE.a ] && echo Hypre)
hypre:
	[ "Hypre" != "$(biltH)" ] || (mkdir -p $(psBuild)/hypre )
	[ "Hypre" != "$(biltH)" ] || (cd       $(psBuild)/hypre &&  $(CmakeOpts)  $(PS3rdParty)/hypre-2.11.2/src )
	[ "Hypre" != "$(biltH)" ] || (cd       $(psBuild)/hypre && make -j install)

biltZ=$(shell [ -d $(PS3rdParty)/zlib-1.2.11 ] && [ ! -f $(psInstall)/lib/libzxx.a ] && echo Zlib)
zlib:
	[ "Zlib"  != "$(biltZ)" ] || (mkdir -p $(psBuild)/zlib )
	[ "Zlib"  != "$(biltZ)" ] || (cd       $(psBuild)/zlib &&  $(CmakeOpts)  $(PS3rdParty)/zlib-1.2.11)
	[ "Zlib"  != "$(biltZ)" ] || (cd       $(psBuild)/zlib && make -j install)




distclean: 


clean: 
#	(rm -rf libtiff/build)
#	(rm -rf zlib-1.2.11/build)
#	(rm -rf hypre-2.11.2/src/build)
