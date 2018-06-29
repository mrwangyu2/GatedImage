#include "DicomImageReaderForGatedImage.h"

DicomImageReaderForGatedImage::DicomImageReaderForGatedImage()
{
}

DicomImageReaderForGatedImage::~DicomImageReaderForGatedImage()
{
}

void DicomImageReaderForGatedImage::setDirectory(const char* directory)
{
	int numberOfTimeSlots = getNumberOfTimeSlots(directory);
	int fileNamesSize = this->_reader->GetNumberOfDICOMFileNames();
	int sizeInFileNameArray = fileNamesSize / numberOfTimeSlots;

	for (int i = 0; i < numberOfTimeSlots; i++)
	{
		vtkSmartPointer<vtkDICOMImageReader> reader = vtkSmartPointer<vtkDICOMImageReader>::New();
		vtkSmartPointer<vtkStringArray> fileNameArray = vtkSmartPointer<vtkStringArray>::New();
		fileNameArray->Resize(sizeInFileNameArray);

		this->_readerGroup.push_back(reader);
		this->_fileNameGroup.push_back(fileNameArray);
	}

	int group = 0;
	for (int i = 0; i < fileNamesSize; i++)
	{
		vtkStdString fileName = this->_reader->GetDICOMFileName(i);

		if ((i / sizeInFileNameArray) > group)
		{
			group++;
		}

		this->_fileNameGroup[group]->InsertNextValue(fileName.c_str());
	}

}

int DicomImageReaderForGatedImage::getNumberOfTimeSlots(const char* directory)
{
	this->_reader = vtkSmartPointer<vtkDICOMImageReader>::New();
	this->_reader->SetDirectoryName(directory);
	this->_reader->Update();

    // Sort natural file name in the directory
	vtkDICOMImageReaderVector* dicomFileNames = this->_reader->DICOMFileNames;
	std::sort(dicomFileNames->begin(), dicomFileNames->end(), stlnatstrlt);

	return this->_reader->GetNumberOfTimeSlots();
}

void DicomImageReaderForGatedImage::updateAllSeries()
{
	char currentDir[MAX_PATH];
	GetCurrentDirectory(MAX_PATH, currentDir);

	for (int i = 0; i < this->_fileNameGroup.size(); i++)
	{
		auto fileNames = this->_fileNameGroup.at(i);

		int size = fileNames->GetSize();
		std::string directory(currentDir);
		this->_dataDirectory = directory + "\\temp\\";

		RemoveDirectory(this->_dataDirectory.c_str());
		CreateDirectory(this->_dataDirectory.c_str(), NULL);

		std::string dataDirectory = this->_dataDirectory + std::to_string(i) + "\\";
		this->_dataDirectoryGroup.push_back(dataDirectory);

		CreateDirectory(dataDirectory.c_str(), NULL);

		for (int j = 0; j < size; j++)
		{
			std::string filePath = fileNames->GetValue(j);
			filePath.replace(filePath.find("/"), 1, "\\");

			int pos = filePath.find_last_of('\\');
			std::string fileName(filePath.substr(pos + 1));

			std::string destinationFile = dataDirectory + fileName;
			CopyFile(filePath.c_str(), destinationFile.c_str(), TRUE);
			cout << GetLastError();
		}
	}

	for (int i = 0; i < this->_readerGroup.size(); i++)
	{
		this->_readerGroup.at(i)->SetDirectoryName(this->_dataDirectoryGroup.at(i).c_str());
		this->_readerGroup.at(i)->Update();
	}
}


std::vector<vtkSmartPointer<vtkDICOMImageReader>>& DicomImageReaderForGatedImage::getReadGroup()
{
	return this->_readerGroup;
}
