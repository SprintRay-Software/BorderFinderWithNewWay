#pragma once
#ifndef BORDERFINDERCONSTANTS_H
#define BORDERFINDERCONSTANTS_H

#include <QString>
#include <QColor>
#include <qdir.h>
#include <QCoreApplication>

namespace BORDERFINDERCONSTANTS {
    // Location of the slicer and output image (UPDATE THIS FOR YOUR SYSTEM!)
#ifdef __APPLE__
    const QString SQUASHED_IMG_PATH = QCoreApplication::applicationDirPath() + "/Borderfinder";
    const QString SLICER_PATH = QCoreApplication::applicationDirPath() + "/Borderfinder/DLPSlicerSingleMac";
#else
    const QString SQUASHED_IMG_PATH = QCoreApplication::applicationDirPath() + "/Borderfinder";
    const QString SLICER_PATH = QCoreApplication::applicationDirPath() + "/Borderfinder/Command Line/DLPSlicerSingleWin.exe";
#endif
    const QString CLONE_STL_PATH = QCoreApplication::applicationDirPath() + "/Borderfinder";

    // Color of the background. All other colrs are considered to be the object.
    const QColor EMPTY_PIXEL = QColor(0, 0, 0);
    // The color to draw when adding lines to merge islands
    const QColor FILLED_PIXEL = QColor(255, 255, 255);
    // The color to draw border pixels as when outputting image
    const QColor BORDER_PIXEL = QColor(255, 255, 0);

    // identity matrix
    const QString DEFAULT_TRANSFORMATION_MATRIX = "1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1";

    const QString FILE_NOT_FOUND = "BorderFinder: The file %s does not exist. Please check that the filepath is correct when calling BorderFinder's constructor";
    const QString FILE_TYPE_INCORRECT = "BorderFinder: The input %s is not an stl, obj, or valid image for QImage. (Note that spr does not work on this branch.) Please check that " + SQUASHED_IMG_PATH + " is correct in borderfinderConstants.h and that the filepath is correct when calling BorderFinder's constructor";
    const QString NEXT_EDGE_NOT_FOUND = "outlineShape: error in finding next edge. This is a bug, please send the input file to the maintainer of this algorithm. Current point is %i, %i.";
};

#endif // BORDERFINDERCONSTANTS_H


