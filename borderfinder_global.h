#ifndef BORDERFINDER_GLOBAL_H
#define BORDERFINDER_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(BORDERFINDER_LIBRARY)
#  define BORDERFINDERSHARED_EXPORT Q_DECL_EXPORT
#else
#  define BORDERFINDERSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // BORDERFINDER_GLOBAL_H
