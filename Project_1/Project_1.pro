TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    project1_problems.cpp


HEADERS += \
    project1_problems.h

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib
LIBS += -larmadillo -llapack -lblas
