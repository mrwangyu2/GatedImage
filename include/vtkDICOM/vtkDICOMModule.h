
#ifndef VTKDICOM_EXPORT_H
#define VTKDICOM_EXPORT_H

#ifdef VTKDICOM_STATIC_DEFINE
#  define VTKDICOM_EXPORT
#  define VTKDICOM_NO_EXPORT
#else
#  ifndef VTKDICOM_EXPORT
#    ifdef vtkDICOM_EXPORTS
        /* We are building this library */
#      define VTKDICOM_EXPORT __declspec(dllexport)
#    else
        /* We are using this library */
#      define VTKDICOM_EXPORT __declspec(dllimport)
#    endif
#  endif

#  ifndef VTKDICOM_NO_EXPORT
#    define VTKDICOM_NO_EXPORT 
#  endif
#endif

#ifndef VTKDICOM_DEPRECATED
#  define VTKDICOM_DEPRECATED __declspec(deprecated)
#  define VTKDICOM_DEPRECATED_EXPORT VTKDICOM_EXPORT __declspec(deprecated)
#  define VTKDICOM_DEPRECATED_NO_EXPORT VTKDICOM_NO_EXPORT __declspec(deprecated)
#endif

#define DEFINE_NO_DEPRECATED 0
#if DEFINE_NO_DEPRECATED
# define VTKDICOM_NO_DEPRECATED
#endif

/* AutoInit dependencies.  */
#include "vtkIOImageModule.h"
#include "vtkIOSQLModule.h"

#endif
