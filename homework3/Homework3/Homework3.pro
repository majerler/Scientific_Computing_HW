TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -Xpreprocessor -fopenmp -lomp -I/usr/local/include
QMAKE_LFLAGS += -lomp
LIBS += -L /usr/local/lib /usr/local/lib/libomp.dylib

SOURCES += \
        ../classes/cf_2.cpp \
        ../classes/godunov.cpp \
        ../classes/grid2d.cpp \
        ../classes/math_tools.cpp \
        ../classes/sl_method.cpp \
        main.cpp

HEADERS += \
    ../classes/cf_2.h \
    ../classes/godunov.h \
    ../classes/grid2d.h \
    ../classes/math_tools.h \
    ../classes/sl_method.h
