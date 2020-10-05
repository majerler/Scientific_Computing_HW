TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        ../Lab1/grid2d.cpp \
        cf_2.cpp \
        eno_advection.cpp \
        main.cpp

QMAKE_CXXFLAGS += -Xpreprocessor -fopenmp -lomp -I/usr/local/include
QMAKE_LFLAGS += -lomp
LIBS += -L /usr/local/lib /usr/local/lib/libomp.dylib

HEADERS += \
    ../Lab1/grid2d.h \
    cf_2.h \
    eno_advection.h
