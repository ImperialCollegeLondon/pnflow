


all: directories Linux
#all: Win32
#all: Win64
#all: winVS
# "c:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\bin\amd64_x86\vcvarsamd64_x86.bat"
# nmake

USE_HYPRE=1
psDir ?= $(abspath ..)
include  $(psDir)/Makefile.in

#COPTS     = -g  -Wall  -fno-inline 
COPTS     = -O3 -Wall -std=c++11


psBinDir ?= ../../bin
BuildDir ?= ../../build/pnflow

directories:
	@mkdir -p $(psBinDir)
	@mkdir -p $(BuildDir)

#psCXXWin 		= i686-w64-mingw32-g++
psCXXWin 	?= x86_64-w64-mingw32-g++-6
psCXX 		?= g++

#EXTRACFLAGS= -Wall -pedantic
#MMD is for creating dependancy
psCFLAGS +=  -IporeModel -I. -IpostProcess 


srcs =  \
	poreModel/layerApex.cpp  poreModel/cornerApex.cpp  \
	poreModel/polygon.cpp  poreModel/polygon_calcR.cpp  \
	poreModel/polygonDrain.cpp  poreModel/polygonImb.cpp \
	ElementConstSolver.cpp  ElementPore.cpp    elem_Model.cpp  \
	postProcess/netStatistics.cpp Element.cpp  \
	ElementThroat.cpp    \
	inputData.cpp   netsimSP.cpp  netsim.cpp  netsim_drain.cpp   \
	netsim_imbibe.cpp   utility.cpp   hypre.cpp pnflow.cpp  postProcess/vtuWriter.cpp


objs = $(addprefix $(BuildDir)/,$(subst /,__,$(srcs:.cpp=.o)))
deps = $(addprefix $(BuildDir)/,$(subst /,__,$(srcs:.cpp=.d)))

# delete the following if nmake (from visual studio) complains, uncomment for linux compillation
Linux: Hypre $(objs)
	$(psCXX)  $(objs) -o $(psBinDir)/pnflow$(OPT) $(psLFLAGS)


-include $(deps)


$(BuildDir)/%.o: $(subst $(BuildDir)/,,$(subst __,/,$(subst .o,.cpp,$@)))
	$(psCXX)  $(psCFLAGS)  $(EXTRACFLAGS) -c $(subst $(BuildDir)/,,$(subst __,/,$(subst .o,.cpp,$@))) -o $@



.PHONY: clean $(deps)

clean:
	$(RM) $(objs) $(deps)




Hypre:
	echo Hypre
	(cd $(psDir)/../thirdparty && make hypre)



winVS:
	@cl -I..\libf2c /nologo /MT /EHsc /Ox  $(srcs) ..\..\bin\vcf2c.lib  /link -LARGEADDRESSAWARE /out:..\bin\pnflow.exe
	@del *.obj


