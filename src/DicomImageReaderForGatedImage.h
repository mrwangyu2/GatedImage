#ifndef DicomImageReaderForGatedImage_h
#define DicomImageReaderForGatedImage_h

#include <windows.h>
#include <vector>
#include <algorithm>
#include <filesystem>

#include <vtkSmartPointer.h>
#include <vtkDICOMImageReader.h>
#include <vtkStringArray.h>

#include "NaturalSort.h"


class DicomImageReaderForGatedImage
{
public:
	DicomImageReaderForGatedImage();
	~DicomImageReaderForGatedImage();
	void setDirectory(const char* dn);

	int getNumberOfTimeSlots(const char* directory);

	void updateAllSeries();

	std::vector<vtkSmartPointer<vtkDICOMImageReader>> _readerGroup;

private:
	std::string _dataDirectory;
	std::vector<std::string> _dataDirectoryGroup;
	std::vector<vtkSmartPointer<vtkStringArray>> _fileNameGroup;
	vtkSmartPointer<vtkDICOMImageReader> _reader;



public:
	std::vector<vtkSmartPointer<vtkDICOMImageReader>>& getReadGroup();
};


#endif

