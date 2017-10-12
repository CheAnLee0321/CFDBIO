TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    ddmodel.cpp \
    cfd.cpp \
    montecarlo.cpp \
    fvmesh.cpp \
    cdmodel.cpp

HEADERS += \
    ddmodel.h \
    cfd.h \
    montecarlo.h \
    fvmesh.h \
    Parameter.h \
    cdmodel.h

QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp
