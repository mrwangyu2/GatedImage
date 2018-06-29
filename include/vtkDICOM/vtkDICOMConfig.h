/*=========================================================================

  Program: DICOM for VTK

  Copyright (c) 2012-2013 David Gobbi
  All rights reserved.
  See Copyright.txt or http://www.cognitive-antics.net/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkDICOMConfig_h
#define vtkDICOMConfig_h

/* Configuration information. */
#define DICOM_BUILD_SHARED_LIBS
#define DICOM_BUILD_TESTING
/* #undef DICOM_USE_GDCM */
/* #undef DICOM_USE_DCMTK */
#define DICOM_USE_VTKZLIB

/* Version number. */
#define DICOM_MAJOR_VERSION 0
#define DICOM_MINOR_VERSION 7
#define DICOM_PATCH_VERSION 10
#define DICOM_SHORT_VERSION "0.7"
#define DICOM_VERSION "0.7.10"

/* Legacy (for backwards compatibility) */
#define DICOM_BUILD_VERSION DICOM_PATCH_VERSION

/* For compatibility with VTK 7.1 */
#ifndef VTK_DELETE_FUNCTION
#if (__cplusplus >= 201103L) || (defined(_MSC_VER) && _MSC_VER >= 1800)
#define VTK_DELETE_FUNCTION =delete
#else
#define VTK_DELETE_FUNCTION
#endif
#endif

#endif
