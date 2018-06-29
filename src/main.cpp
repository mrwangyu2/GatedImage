#include "GatedImageViewer.h"

int main(int argc, char* argv[])
{
	// Verify input arguments
	if (argc != 2)
	{
		std::cout << "Usage: " << argv[0]
			<< " FolderName" << std::endl;
		return EXIT_FAILURE;
	}

	std::string folder = argv[1];
	//std::string folder = "C:\\VTK\\vtkdata-5.8.0\\Data\\DicomTestImages";

	GatedImageViewer gatedImagerViewer;
	gatedImagerViewer.parserData(folder);
	gatedImagerViewer.buildPipelineForFreedomBreathe();
	return EXIT_SUCCESS;
}

