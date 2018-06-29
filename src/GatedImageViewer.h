#ifndef GatedImageViewer_h
#define GatedImageViewer_h

// some standard vtk headers
#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkActor.h>

// headers needed for this example
#include <vtkImageActor.h>
#include <vtkCamera.h>
#include <vtkDICOMImageReader.h>
#include <vtkInteractorStyleImage.h>
#include <vtkActor2D.h>
#include <vtkTextProperty.h>
#include <vtkTextMapper.h>
// needed to easily convert int to std::string
#include <sstream>

#include "ImageViewerForGatedImage.h"
#include "DicomImageReaderForGatedImage.h"
#include "MyVtkInteractorStyleImage.h"


class vtkImageActor;
class vtkCamera;
class ImageViewerForGatedImage;
class myVtkInteractorStyleImage;

class  GatedImageViewer
{
public:
	GatedImageViewer();
	~GatedImageViewer();

	void parserData(const std::string& folder);
	void buildPipelineForFreedomBreathe();

	void createImageViewerGroup();

	void setInputConnectionToImageView();

private:
	DicomImageReaderForGatedImage _dicomImageReader;

	vtkImageMapToWindowLevelColors  *_windowLevel;
	vtkRenderWindow                 *_renderWindow;
	vtkRenderer                     *_renderer;
	vtkImageActor                   *_imageActor;
	vtkRenderWindowInteractor       *_interactor;
	myVtkInteractorStyleImage       *_interactorStyle;
	int _currentActorNumber; 
	bool _sequence;

	std::vector<vtkSmartPointer<ImageViewerForGatedImage>> _imageViewerGroup;

	void installPipeline();
	void createVTKPipelineObject();
	void setImageViewer(vtkSmartPointer<ImageViewerForGatedImage>& imageViewer, vtkSmartPointer<vtkDICOMImageReader>& reader);
	void setVTKPipelineOject();
public:
	void displayNextActor();

	void displayActor();

	void calculatedActorNumber();

};

#endif
