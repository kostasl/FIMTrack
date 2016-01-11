#-------------------------------------------------
#
# Project created by QtCreator 2013-02-13T15:27:57
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

#Added for linux build
QT += printsupport

TARGET = FIMTrack
TEMPLATE = app

include(Algorithm/algorithm.pri)
include(Data/data.pri)
include(Control/control.pri)
include(Utility/utility.pri)
include(Configuration/configuration.pri)
include(GUI/GUI.pri)
include(Main/main.pri)

CONFIG   += app_bundle

##macx {
# #   QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.9
##    QMAKE_CXXFLAGS += -stdlib=libc++
##    INCLUDEPATH += /usr/local/include
#    LIBS += -L/usr/local/lib \
#            -lopencv_core \
#            -lopencv_highgui \
#            -lopencv_imgproc

##    QMAKE_CXXFLAGS_WARN_ON = -Wno-unused-variable -Wno-reorder
#}

#Works with OPENCV 2.4 - 3.0 has a problem with persistence.hpp
unix {

    QMAKE_CXXFLAGS += -std=c++11

#    INCLUDEPATH += /usr/include
    INCLUDEPATH +=  /home/klagogia/opencv3git-install/include

    LIBS += -L/home/klagogia/opencv3git-install/lib \
            -lopencv_core \
            -lopencv_highgui \
            -lopencv_imgproc \
            -lopencv_imgcodecs


#    LIBS += -L/usr/lib/x86_64-linux-gnu \
#            -lopencv_core \
#            -lopencv_highgui \
#            -lopencv_imgproc


    QMAKE_CXXFLAGS_WARN_ON = -Wno-unused-variable -Wno-reorder
}

win32{

    INCLUDEPATH += D:\\Arbeit\\Bibliotheken\\OpenCV\\OpenCV-2.4.3\\include

    CONFIG(debug,debug|release){
        LIBS += D:\\Arbeit\\Bibliotheken\\OpenCV\\OpenCV-2.4.3\\build\\x86\\vc10\\lib\\opencv_core243d.lib
        LIBS += D:\\Arbeit\\Bibliotheken\\OpenCV\\OpenCV-2.4.3\\build\\x86\\vc10\\lib\\opencv_highgui243d.lib
        LIBS += D:\\Arbeit\\Bibliotheken\\OpenCV\\OpenCV-2.4.3\\build\\x86\\vc10\\lib\\opencv_imgproc243d.lib

    }

    CONFIG(release,debug|release){
        LIBS += D:\\Arbeit\\Bibliotheken\\OpenCV\\OpenCV-2.4.3\\build\\x86\\vc10\\lib\\opencv_core243.lib
        LIBS += D:\\Arbeit\\Bibliotheken\\OpenCV\\OpenCV-2.4.3\\build\\x86\\vc10\\lib\\opencv_highgui243.lib
        LIBS += D:\\Arbeit\\Bibliotheken\\OpenCV\\OpenCV-2.4.3\\build\\x86\\vc10\\lib\\opencv_imgproc243.lib
    }
}
