#include "GatedImageViewer.h"
#include "vtkAutoInit.h"
#include "ImageViewerForGatedImage.h"
VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK was built with vtkRenderingOpenGL2
VTK_MODULE_INIT(vtkInteractionStyle);


GatedImageViewer::GatedImageViewer()
	: _windowLevel(nullptr)
	, _renderWindow(nullptr)
	, _renderer(nullptr)
	, _imageActor(nullptr)
	, _interactor(nullptr)
	, _interactorStyle(nullptr)
	, _currentActorNumber(0)
	, _sequence(true)
{
}


GatedImageViewer::~GatedImageViewer()
{
}


void GatedImageViewer::parserData(const std::string& folder)
{
	DicomImageReaderForGatedImage dicomImageReader;
	dicomImageReader.setDirectory(folder.c_str());
	dicomImageReader.updateAllSeries();
	this->_dicomImageReader = dicomImageReader;
}


void GatedImageViewer::buildPipelineForFreedomBreathe()
{
	this->createVTKPipelineObject();
	this->setVTKPipelineOject();

	this->createImageViewerGroup();

	this->installPipeline();

	this->_interactorStyle->StartTimer();
	this->_interactor->Start();
}


void GatedImageViewer::createVTKPipelineObject()
{
	this->_renderWindow = vtkRenderWindow::New();
	this->_renderer = vtkRenderer::New();
	this->_interactor = vtkRenderWindowInteractor::New();
	this->_interactorStyle = myVtkInteractorStyleImage::New();
}


void GatedImageViewer::setVTKPipelineOject()
{
	this->_interactorStyle->setGatedImageViewer(this);
	this->_interactorStyle->UseTimersOn();
	this->_interactorStyle->SetTimerDuration(100);
}


void GatedImageViewer::createImageViewerGroup()
{
	std::vector<vtkSmartPointer<vtkDICOMImageReader>> readerGroup = this->_dicomImageReader.getReadGroup();

	for (auto reader : readerGroup)
	{
		vtkSmartPointer<ImageViewerForGatedImage> imageViewer = vtkSmartPointer<ImageViewerForGatedImage>::New();
		this->setImageViewer(imageViewer, reader);
		this->_imageViewerGroup.push_back(imageViewer);
	}
}


void GatedImageViewer::setImageViewer(vtkSmartPointer<ImageViewerForGatedImage>& imageViewer, vtkSmartPointer<vtkDICOMImageReader>& reader)
{
	imageViewer->SetRenderer(this->_renderer);
	imageViewer->SetRenderWindow(this->_renderWindow);

	imageViewer->SetInputConnection(reader->GetOutputPort());
	imageViewer->SetupInteractor(this->_interactor);

	imageViewer->SetSliceOrientation(ImageViewerForGatedImage::SLICE_ORIENTATION_XY);
	//imageViewer->SetSliceOrientation(ImageViewerForGatedImage::SLICE_ORIENTATION_YZ);
	//imageViewer->SetSliceOrientation(ImageViewerForGatedImage::SLICE_ORIENTATION_XZ);

	imageViewer->Rotate180();
	imageViewer->SetPetLookupTable();

	imageViewer->SetColorWindow(30000);
	imageViewer->SetColorLevel(15000);

	imageViewer->GetRenderWindow()->SetSize(600, 600);
	imageViewer->GetRenderer()->SetBackground(1, 1, 1);
	imageViewer->Render();
	imageViewer->GetRenderer()->ResetCamera();
	imageViewer->Render();

}


void GatedImageViewer::setInputConnectionToImageView()
{

}


void GatedImageViewer::installPipeline()
{
	if (this->_renderWindow && this->_renderer)
	{
		this->_renderWindow->AddRenderer(this->_renderer);
	}

	if (this->_interactor)
	{
		//if (!this->_interactorStyle)
		//{
		//	this->_interactorStyle = vtkInteractorStyleImage::New();
		//	ImageViewerForGatedImageCallback *cbk = ImageViewerForGatedImageCallback::New();
		//	cbk->IV = this;
		//	this->_interactorStyle->AddObserver(
		//		vtkCommand::WindowLevelEvent, cbk);
		//	this->_interactorStyle->AddObserver(
		//		vtkCommand::StartWindowLevelEvent, cbk);
		//	this->_interactorStyle->AddObserver(
		//		vtkCommand::ResetWindowLevelEvent, cbk);
		//	cbk->Delete();
		//}

		this->_interactor->SetInteractorStyle(this->_interactorStyle);
		this->_interactor->SetRenderWindow(this->_renderWindow);
	}

	if (this->_renderer && this->_imageViewerGroup.size() > 0)
	{
		for (auto imageViewer : this->_imageViewerGroup)
		{
			imageViewer->addActorToRender();
			imageViewer->addOutputPortToMapper();


		}
	}



}


void GatedImageViewer::displayNextActor()
{
	calculatedActorNumber();
	displayActor();
}


void GatedImageViewer::calculatedActorNumber()
{
	if (this->_sequence == true)
	{
		this->_currentActorNumber++;
		if (this->_currentActorNumber >= this->_imageViewerGroup.size())
		{
			this->_currentActorNumber = this->_currentActorNumber - 2;
			this->_sequence = false;
		}
	}
	else
	{
		this->_currentActorNumber--;

		if (this->_currentActorNumber < 0)
		{
			this->_currentActorNumber++;
			this->_sequence = true;
		}
	}
}


void GatedImageViewer::displayActor()
{
	vtkSmartPointer<ImageViewerForGatedImage> imageViewer;
	for (int i = 0; i < this->_imageViewerGroup.size(); i++)
	{
		imageViewer = this->_imageViewerGroup.at(i);
		if (i == this->_currentActorNumber)
		{
			imageViewer->GetImageActor()->VisibilityOn();
		}
		else
		{
			imageViewer->GetImageActor()->VisibilityOff();
		}
	}
	imageViewer->Render();
}


