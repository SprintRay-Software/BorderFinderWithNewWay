#-------------------------------------------------
#
# Project created by QtCreator 2020-10-27T09:24:33
#
#-------------------------------------------------

QT       += gui

TARGET = BorderFinder
TEMPLATE = lib

DEFINES += _USE_MATH_DEFINES

DEFINES += OM_STATIC_BUILD

DEFINES += BORDERFINDER_LIBRARY

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
        borderfinder.cpp \
    clipper.cpp \
    grahamscan.cpp \
    node.cpp

HEADERS += \
        borderfinder.h \
        borderfinder_global.h \ 
    clipper.h \
    borderfinderconstants.h \
    grahamscan.h \
    node.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}

macx {
    INCLUDEPATH += $$PWD/include_Mac

    LIBS += -L$$PWD/lib_Mac -lOpenMeshCore
    LIBS += -L$$PWD/lib_Mac -lOpenMeshTools

}

#CONFIG += debug_and_release
win32 {
#    CONFIG(debug,debug|release)
#    {
#        INCLUDEPATH += $$PWD/include
#        LIBS += -L$$PWD/libd -lOpenMeshCored
#        LIBS += -L$$PWD/libd -lOpenMeshToolsd
#    }
    CONFIG(release,debug|release)
    {
        INCLUDEPATH += $$PWD/include
        LIBS += -L$$PWD/lib -lOpenMeshCore
        LIBS += -L$$PWD/lib -lOpenMeshTools
    }
}

#win32 {
#    INCLUDEPATH += $$PWD/include

#    LIBS += -L$$PWD/lib -lOpenMeshCore
#    LIBS += -L$$PWD/lib -lOpenMeshCored
#    LIBS += -L$$PWD/lib -lOpenMeshTools
#    LIBS += -L$$PWD/lib -lOpenMeshToolsd
#}
