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
    clipper.cpp

HEADERS += \
        borderfinder.h \
        borderfinder_global.h \ 
    clipper.h \
    borderfinderconstants.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}

macx {
    INCLUDEPATH += $$PWD/include_Mac

    LIBS += -L$$PWD/lib_Mac -lOpenMeshCore
    LIBS += -L$$PWD/lib_Mac -lOpenMeshTools

}
win32 {
    CONFIG(release,debug|release){
message(release)
        INCLUDEPATH += $$PWD/include
        LIBS += -L$$PWD/lib -lOpenMeshCore
        LIBS += -L$$PWD/lib -lOpenMeshTools
        target.path = D:\QTProjects\MR\Moonray\PlanU\Moonray-0319\Moonray-0319\AutoPlace\Borderfinder\lib\
        target.files += D:\QTProjects\MR\Nest_DLL\build-BorderFinder-Desktop_Qt_5_14_2_MSVC2017_64bit-Release\release\BorderFinder.dll
    }else{
message(debug)
        INCLUDEPATH += $$PWD/include
        LIBS += -L$$PWD/lib -lOpenMeshCored
        LIBS += -L$$PWD/lib -lOpenMeshToolsd
        target.path = D:\QTProjects\MR\Moonray\PlanU\Moonray-0319\Moonray-0319\AutoPlace\Borderfinder\lib\
        target.files += D:\QTProjects\MR\Nest_DLL\build-BorderFinder-Desktop_Qt_5_14_2_MSVC2017_64bit-Debug\debug\BorderFinder.dll
    }

}
INSTALLS += target
