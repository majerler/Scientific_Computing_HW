TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        autoencoder.cpp \
        main.cpp \
        net_tools.cpp

HEADERS += \
    autoencoder.h \
    net_tools.h
