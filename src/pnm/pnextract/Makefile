


all: directories Linux
#all: Win
#all: winVS



ifeq "$(OPT)" "dbg"
 COPTS     = -g -fno-inline  -D_debugCompile_  $(CFLAGSExtra)
 BuildDir = ../../build/pnextractdbg
else
 COPTS     = -O3  $(CFLAGSExtra)
endif

ifdef USE_CORRELATIONS



endif






USE_ZLIB=1
USE_TIFF=1
USE_VXLIMG=1
psDir ?= $(abspath ..)
include  $(psDir)/Makefile.in


psBinDir ?= $(psDir)/../bin
BuildDir ?= $(psDir)/../build/pnextract

directories:
	@mkdir -p $(psBinDir)
	@mkdir -p $(BuildDir)




psCXXWin 	?= x86_64-w64-mingw32-g++
psCXX 		?= g++


psCFLAGS := $(psCFLAGS) $(COPTS)   -Wall -pedantic -MMD -MP    -I.
psLFLAGS := $(psLFLAGS) -static-libgcc -static-libstdc++ 


srcs = \
	blockNet.cpp  medialSurf.cpp  nextract.cpp \
	writers_vtk.cpp  writers_vxl.cpp $(gfsrcs)


objs = $(addprefix $(BuildDir)/,$(subst /,__,$(srcs:.cpp=.o)))
deps = $(addprefix $(BuildDir)/,$(subst /,__,$(srcs:.cpp=.d)))

# delete the following if nmake (from visual studio) complains, uncomment for linux compillation
Linux: $(objs)
	@echo Building $(psBinDir)/pnextract$(OPT):
	$(psCXX)  $(objs) -o $(psBinDir)/pnextract$(OPT) $(psLFLAGS)

-include $(deps)


$(BuildDir)/%.o: $(subst $(BuildDir)/,,$(subst __,/,$(subst .o,.cpp,$@)))
	$(psCXX)  $(psCFLAGS)  $(EXTRACFLAGS) -c $(subst $(BuildDir)/,,$(subst __,/,$(subst .o,.cpp,$@))) -o $@



.PHONY: clean $(deps)

clean:
	$(RM) -f $(objs) $(deps) nelmTo*Par.o



Win:
	$(psCXXWin)  $(psCFLAGS)  $(srcs) -o $(psBinDir)/pnextract$(OPT).exe   $(psLFLAGS)


winVS:
	@cl /nologo /MT /EHsc /Ox  $(srcs)  /link -LARGEADDRESSAWARE /out:..\bin\pnextract$(OPT).exe
	@del *.obj





