QT -= gui

TARGET = pnflow_tom
CONFIG += c++11 console
CONFIG -= app_bundle
CONFIG -= app_bundle
TEMPLATE = app

QMAKE_CC  = gcc
QMAKE_CXX = g++

CONFIG += debug_and_release

DESTDIR = ./../../bin
#LIBS += -L/usr/local/lib -lmath
INCLUDEPATH = ./poreModel/ ./postProcess/ ./  ./../../include/
LIBS += -L./../../lib -lHYPRE

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0
 
SOURCES += \
    ./poreModel/cornerApex.cpp \
    ./poreModel/layerApex.cpp \
    ./poreModel/polygon.cpp \
    ./poreModel/polygon_calcR.cpp \
    ./poreModel/polygonDrain.cpp \
    ./poreModel/polygonImb.cpp \
    ./postProcess/netStatistics.cpp \
    ./postProcess/vtuWriter.cpp \
    ./elem_Model.cpp \
    ./Element.cpp \
    ./ElementConstSolver.cpp \
    ./ElementPore.cpp \
    ./ElementThroat.cpp \
    ./hypre.cpp \
    ./inputData.cpp \
    ./pnflow.cpp \
    ./netsim.cpp \
    ./netsim_drain.cpp \
    ./netsim_imbibe.cpp \
    ./netsimSP.cpp \
    ./utility.cpp




DISTFILES += \
    ./pnflow_input.dat \
    ./Makefile

HEADERS += \
    ./poreModel/apex.h \
    ./poreModel/cornerApex.h \
    ./poreModel/layerApex.h \
    ./poreModel/polygon.h \
    ./postProcess/netStatistics.h \
    ./postProcess/typses.h \
    ./postProcess/vtuWriter.h \
    ./compareFuncs.h \
    ./elem_Model.h \
    ./Element.h \
    ./fluid.h \
    ./hypre.h \
    ./inputData.h \
    ./inputFile.h \
    ./netsim.h \
    ./netsim_data.h \
    ./NetworkTransform.h \
    ./solver.h \
    ./sortedEvents.h \
    ./threeSome.h \
    ./utility.h



