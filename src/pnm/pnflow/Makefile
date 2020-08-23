
all:  msApp

#need to compile hypre with fPIC for this to work
#all: msLib
#msCFLAGS += -DMAIN -fPIC


msSrc ?= $(abspath ../..)

EXTRACFLAGS= -Wall -pedantic
#MMD is for creating dependancy
msCFLAGS +=  -I. -IpsPrc -Ipolygon -MMD -MP


srcs = \
	polygon/layerApex.cpp  polygon/cornerApex.cpp  \
	polygon/polygon.cpp  polygon/polygon_calcR.cpp  \
	polygon/polygonDrain.cpp  polygon/polygonImb.cpp \
	ElementConstSolver.cpp  ElementPore.cpp    elem_Model.cpp  \
	psPrc/netStatistics.cpp Element.cpp  \
	ElementThroat.cpp    cnflow.cpp \
	inputData.cpp   FlowData.cpp  FlowDomain.cpp  FlowDomain_drain.cpp   \
	FlowDomain_imbibe.cpp   NetworkTransform.cpp   hypreSolver.cpp psPrc/vtuWriter.cpp

tsts = tests/readCAs.py





USE_msTEST=1
USE_HYPRESEQ=1
USE_msMAKE=1
include  ${msSrc}/script/Makefile.in

