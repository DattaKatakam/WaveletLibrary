#ifndef WAVELETLIB_GLOBAL_H
#define WAVELETLIB_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(WAVELETLIB_LIBRARY)
#  define WAVELETLIB_EXPORT Q_DECL_EXPORT
#else
#  define WAVELETLIB_EXPORT Q_DECL_IMPORT
#endif

#endif // WAVELETLIB_GLOBAL_H
