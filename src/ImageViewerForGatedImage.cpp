#include "ImageViewerForGatedImage.h"

#include "vtkCamera.h"
#include "vtkCommand.h"
#include "vtkImageActor.h"
#include "vtkImageData.h"
#include "vtkImageMapper3D.h"
#include "vtkImageMapToWindowLevelColors.h"
#include "vtkInformation.h"
#include "vtkInteractorStyleImage.h"
#include "vtkObjectFactory.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkStreamingDemandDrivenPipeline.h"


vtkStandardNewMacro(ImageViewerForGatedImage);


ImageViewerForGatedImage::ImageViewerForGatedImage()
{
	this->RenderWindow = NULL;
	this->Renderer = NULL;
	this->ImageActor = vtkImageActor::New();
	this->WindowLevel = vtkImageMapToWindowLevelColors::New();
	this->Interactor = NULL;
	this->InteractorStyle = NULL;

	this->Slice = 0;
	this->FirstRender = 1;
	this->SliceOrientation = ImageViewerForGatedImage::SLICE_ORIENTATION_XY;

	// Setup the pipeline

	vtkRenderWindow *renwin = vtkRenderWindow::New();
	this->SetRenderWindow(renwin);
	renwin->Delete();

	vtkRenderer *ren = vtkRenderer::New();
	this->SetRenderer(ren);
	ren->Delete();

	this->InstallPipeline();
}


ImageViewerForGatedImage::~ImageViewerForGatedImage()
{
	if (this->WindowLevel)
	{
		this->WindowLevel->Delete();
		this->WindowLevel = NULL;
	}

	if (this->ImageActor)
	{
		this->ImageActor->Delete();
		this->ImageActor = NULL;
	}

	if (this->Renderer)
	{
		this->Renderer->Delete();
		this->Renderer = NULL;
	}

	if (this->RenderWindow)
	{
		this->RenderWindow->Delete();
		this->RenderWindow = NULL;
	}

	if (this->Interactor)
	{
		this->Interactor->Delete();
		this->Interactor = NULL;
	}

	if (this->InteractorStyle)
	{
		this->InteractorStyle->Delete();
		this->InteractorStyle = NULL;
	}
}


void ImageViewerForGatedImage::SetupInteractor(vtkRenderWindowInteractor *arg)
{
	if (this->Interactor == arg)
	{
		return;
	}

	this->UnInstallPipeline();

	if (this->Interactor)
	{
		this->Interactor->UnRegister(this);
	}

	this->Interactor = arg;

	if (this->Interactor)
	{
		this->Interactor->Register(this);
	}

	this->InstallPipeline();

	if (this->Renderer)
	{
		this->Renderer->GetActiveCamera()->ParallelProjectionOn();
	}
}


void ImageViewerForGatedImage::SetRenderWindow(vtkRenderWindow *arg)
{
	if (this->RenderWindow == arg)
	{
		return;
	}

	this->UnInstallPipeline();

	if (this->RenderWindow)
	{
		this->RenderWindow->UnRegister(this);
	}

	this->RenderWindow = arg;

	if (this->RenderWindow)
	{
		this->RenderWindow->Register(this);
	}

	this->InstallPipeline();
}


void ImageViewerForGatedImage::SetRenderer(vtkRenderer *arg)
{
	if (this->Renderer == arg)
	{
		return;
	}

	this->UnInstallPipeline();

	if (this->Renderer)
	{
		this->Renderer->UnRegister(this);
	}

	this->Renderer = arg;

	if (this->Renderer)
	{
		this->Renderer->Register(this);
	}

	this->InstallPipeline();
	this->UpdateOrientation();
}


void ImageViewerForGatedImage::SetSize(int a, int b)
{
	this->RenderWindow->SetSize(a, b);
}


int* ImageViewerForGatedImage::GetSize()
{
	return this->RenderWindow->GetSize();
}


void ImageViewerForGatedImage::GetSliceRange(int &min, int &max)
{
	vtkAlgorithm *input = this->GetInputAlgorithm();
	if (input)
	{
		input->UpdateInformation();
		int *w_ext = input->GetOutputInformation(0)->Get(
			vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT());
		min = w_ext[this->SliceOrientation * 2];
		max = w_ext[this->SliceOrientation * 2 + 1];
	}
}


int* ImageViewerForGatedImage::GetSliceRange()
{
	vtkAlgorithm *input = this->GetInputAlgorithm();
	if (input)
	{
		input->UpdateInformation();
		return input->GetOutputInformation(0)->Get(
			vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()) +
			this->SliceOrientation * 2;
	}
	return NULL;
}


int ImageViewerForGatedImage::GetSliceMin()
{
	int *range = this->GetSliceRange();
	if (range)
	{
		return range[0];
	}
	return 0;
}


int ImageViewerForGatedImage::GetSliceMax()
{
	int *range = this->GetSliceRange();
	if (range)
	{
		return range[1];
	}
	return 0;
}


void ImageViewerForGatedImage::SetSlice(int slice)
{
	int *range = this->GetSliceRange();
	if (range)
	{
		if (slice < range[0])
		{
			slice = range[0];
		}
		else if (slice > range[1])
		{
			slice = range[1];
		}
	}

	if (this->Slice == slice)
	{
		return;
	}

	this->Slice = slice;
	this->Modified();

	this->UpdateDisplayExtent();
	this->Render();
}


void ImageViewerForGatedImage::SetSliceOrientation(int orientation)
{
	if (orientation < ImageViewerForGatedImage::SLICE_ORIENTATION_YZ ||
		orientation > ImageViewerForGatedImage::SLICE_ORIENTATION_XY)
	{
		vtkErrorMacro("Error - invalid slice orientation " << orientation);
		return;
	}

	if (this->SliceOrientation == orientation)
	{
		return;
	}

	this->SliceOrientation = orientation;

	// Update the viewer

	int *range = this->GetSliceRange();
	if (range)
	{
		this->Slice = static_cast<int>((range[0] + range[1]) * 0.5);
	}

	this->UpdateOrientation();
	this->UpdateDisplayExtent();

	if (this->Renderer && this->GetInput())
	{
		double scale = this->Renderer->GetActiveCamera()->GetParallelScale();
		this->Renderer->ResetCamera();
		this->Renderer->GetActiveCamera()->SetParallelScale(scale);
	}

	this->Render();
}


void ImageViewerForGatedImage::UpdateOrientation()
{
	// Set the camera position

	vtkCamera *cam = this->Renderer ? this->Renderer->GetActiveCamera() : NULL;
	if (cam)
	{
		switch (this->SliceOrientation)
		{
		case ImageViewerForGatedImage::SLICE_ORIENTATION_XY:
			cam->SetFocalPoint(0, 0, 0);
			cam->SetPosition(0, 0, 1); // -1 if medical ?
			cam->SetViewUp(0, 1, 0);
			break;

		case ImageViewerForGatedImage::SLICE_ORIENTATION_XZ:
			cam->SetFocalPoint(0, 0, 0);
			cam->SetPosition(0, -1, 0); // 1 if medical ?
			cam->SetViewUp(0, 0, 1);
			break;

		case ImageViewerForGatedImage::SLICE_ORIENTATION_YZ:
			cam->SetFocalPoint(0, 0, 0);
			cam->SetPosition(1, 0, 0); // -1 if medical ?
			cam->SetViewUp(0, 0, 1);
			break;
		}
	}
}


void ImageViewerForGatedImage::UpdateDisplayExtent()
{
	vtkAlgorithm *input = this->GetInputAlgorithm();
	if (!input || !this->ImageActor)
	{
		return;
	}

	input->UpdateInformation();
	vtkInformation* outInfo = input->GetOutputInformation(0);
	int *w_ext = outInfo->Get(
		vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT());

	// Is the slice in range ? If not, fix it

	int slice_min = w_ext[this->SliceOrientation * 2];
	int slice_max = w_ext[this->SliceOrientation * 2 + 1];
	if (this->Slice < slice_min || this->Slice > slice_max)
	{
		this->Slice = static_cast<int>((slice_min + slice_max) * 0.5);
	}

	// Set the image actor

	switch (this->SliceOrientation)
	{
	case ImageViewerForGatedImage::SLICE_ORIENTATION_XY:
		this->ImageActor->SetDisplayExtent(
			w_ext[0], w_ext[1], w_ext[2], w_ext[3], this->Slice, this->Slice);
		break;

	case ImageViewerForGatedImage::SLICE_ORIENTATION_XZ:
		this->ImageActor->SetDisplayExtent(
			w_ext[0], w_ext[1], this->Slice, this->Slice, w_ext[4], w_ext[5]);
		break;

	case ImageViewerForGatedImage::SLICE_ORIENTATION_YZ:
		this->ImageActor->SetDisplayExtent(
			this->Slice, this->Slice, w_ext[2], w_ext[3], w_ext[4], w_ext[5]);
		break;
	}

	// Figure out the correct clipping range

	if (this->Renderer)
	{
		if (this->InteractorStyle &&
			this->InteractorStyle->GetAutoAdjustCameraClippingRange())
		{
			this->Renderer->ResetCameraClippingRange();
		}
		else
		{
			vtkCamera *cam = this->Renderer->GetActiveCamera();
			if (cam)
			{
				double bounds[6];
				this->ImageActor->GetBounds(bounds);
				double spos = bounds[this->SliceOrientation * 2];
				double cpos = cam->GetPosition()[this->SliceOrientation];
				double range = fabs(spos - cpos);
				double *spacing = outInfo->Get(vtkDataObject::SPACING());
				double avg_spacing =
					(spacing[0] + spacing[1] + spacing[2]) / 3.0;
				cam->SetClippingRange(
					range - avg_spacing * 3.0, range + avg_spacing * 3.0);
			}
		}
	}
}


void ImageViewerForGatedImage::SetPosition(int a, int b)
{
	this->RenderWindow->SetPosition(a, b);
}


int* ImageViewerForGatedImage::GetPosition()
{
	return this->RenderWindow->GetPosition();
}


void ImageViewerForGatedImage::SetDisplayId(void *a)
{
	this->RenderWindow->SetDisplayId(a);
}


void ImageViewerForGatedImage::SetWindowId(void *a)
{
	this->RenderWindow->SetWindowId(a);
}


void ImageViewerForGatedImage::SetParentId(void *a)
{
	this->RenderWindow->SetParentId(a);
}


double ImageViewerForGatedImage::GetColorWindow()
{
	return this->WindowLevel->GetWindow();
}


double ImageViewerForGatedImage::GetColorLevel()
{
	return this->WindowLevel->GetLevel();
}


void ImageViewerForGatedImage::SetColorWindow(double s)
{
	this->WindowLevel->SetWindow(s);
}


void ImageViewerForGatedImage::SetColorLevel(double s)
{
	this->WindowLevel->SetLevel(s);
}


class ImageViewerForGatedImageCallback : public vtkCommand
{
public:
	static ImageViewerForGatedImageCallback *New() { return new ImageViewerForGatedImageCallback; }

	void Execute(vtkObject *caller,
		unsigned long event,
		void *vtkNotUsed(callData)) VTK_OVERRIDE
	{
		if (this->IV->GetInput() == NULL)
		{
			return;
		}

		// Reset

		if (event == vtkCommand::ResetWindowLevelEvent)
		{
			this->IV->GetInputAlgorithm()->UpdateWholeExtent();
			double *range = this->IV->GetInput()->GetScalarRange();
			this->IV->SetColorWindow(range[1] - range[0]);
			this->IV->SetColorLevel(0.5 * (range[1] + range[0]));
			this->IV->Render();
			return;
		}

		// Start

		if (event == vtkCommand::StartWindowLevelEvent)
		{
			this->InitialWindow = this->IV->GetColorWindow();
			this->InitialLevel = this->IV->GetColorLevel();
			return;
		}

		// Adjust the window level here

		vtkInteractorStyleImage *isi =
			static_cast<vtkInteractorStyleImage *>(caller);

		int *size = this->IV->GetRenderWindow()->GetSize();
		double window = this->InitialWindow;
		double level = this->InitialLevel;

		// Compute normalized delta

		double dx = 4.0 *
			(isi->GetWindowLevelCurrentPosition()[0] -
			isi->GetWindowLevelStartPosition()[0]) / size[0];
		double dy = 4.0 *
			(isi->GetWindowLevelStartPosition()[1] -
			isi->GetWindowLevelCurrentPosition()[1]) / size[1];

		// Scale by current values

		if (fabs(window) > 0.01)
		{
			dx = dx * window;
		}
		else
		{
			dx = dx * (window < 0 ? -0.01 : 0.01);
		}
		if (fabs(level) > 0.01)
		{
			dy = dy * level;
		}
		else
		{
			dy = dy * (level < 0 ? -0.01 : 0.01);
		}

		// Abs so that direction does not flip

		if (window < 0.0)
		{
			dx = -1 * dx;
		}
		if (level < 0.0)
		{
			dy = -1 * dy;
		}

		// Compute new window level

		double newWindow = dx + window;
		double newLevel;
		newLevel = level - dy;

		// Stay away from zero and really

		if (fabs(newWindow) < 0.01)
		{
			newWindow = 0.01*(newWindow < 0 ? -1 : 1);
		}
		if (fabs(newLevel) < 0.01)
		{
			newLevel = 0.01*(newLevel < 0 ? -1 : 1);
		}

		this->IV->SetColorWindow(newWindow);
		this->IV->SetColorLevel(newLevel);
		this->IV->Render();
	}

	ImageViewerForGatedImage *IV;
	double InitialWindow;
	double InitialLevel;
};


void ImageViewerForGatedImage::InstallPipeline()
{
	if (this->RenderWindow && this->Renderer)
	{
		this->RenderWindow->AddRenderer(this->Renderer);
	}

	if (this->Interactor)
	{
		if (!this->InteractorStyle)
		{
			this->InteractorStyle = vtkInteractorStyleImage::New();
			ImageViewerForGatedImageCallback *cbk = ImageViewerForGatedImageCallback::New();
			cbk->IV = this;
			this->InteractorStyle->AddObserver(
				vtkCommand::WindowLevelEvent, cbk);
			this->InteractorStyle->AddObserver(
				vtkCommand::StartWindowLevelEvent, cbk);
			this->InteractorStyle->AddObserver(
				vtkCommand::ResetWindowLevelEvent, cbk);
			cbk->Delete();
		}

		this->Interactor->SetInteractorStyle(this->InteractorStyle);
		this->Interactor->SetRenderWindow(this->RenderWindow);
	}

	if (this->Renderer && this->ImageActor)
	{
		this->Renderer->AddViewProp(this->ImageActor);
	}

	if (this->ImageActor && this->WindowLevel)
	{
		this->ImageActor->GetMapper()->SetInputConnection(
			this->WindowLevel->GetOutputPort());
	}
}


void ImageViewerForGatedImage::UnInstallPipeline()
{
	if (this->ImageActor)
	{
		this->ImageActor->GetMapper()->SetInputConnection(NULL);
		//this->ImageActor->VisibilityOff();
	}

	if (this->Renderer && this->ImageActor)
	{
		this->Renderer->RemoveViewProp(this->ImageActor);
	}

	if (this->RenderWindow && this->Renderer)
	{
		this->RenderWindow->RemoveRenderer(this->Renderer);
	}

	if (this->Interactor)
	{
		this->Interactor->SetInteractorStyle(NULL);
		this->Interactor->SetRenderWindow(NULL);
	}
}


void ImageViewerForGatedImage::Render()
{
	if (this->FirstRender)
	{
		// Initialize the size if not set yet

		vtkAlgorithm *input = this->GetInputAlgorithm();
		if (input)
		{
			input->UpdateInformation();
			int *w_ext = this->GetInputInformation()->Get(
				vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT());
			int xs = 0, ys = 0;

			switch (this->SliceOrientation)
			{
			case ImageViewerForGatedImage::SLICE_ORIENTATION_XY:
			default:
				xs = w_ext[1] - w_ext[0] + 1;
				ys = w_ext[3] - w_ext[2] + 1;
				break;

			case ImageViewerForGatedImage::SLICE_ORIENTATION_XZ:
				xs = w_ext[1] - w_ext[0] + 1;
				ys = w_ext[5] - w_ext[4] + 1;
				break;

			case ImageViewerForGatedImage::SLICE_ORIENTATION_YZ:
				xs = w_ext[3] - w_ext[2] + 1;
				ys = w_ext[5] - w_ext[4] + 1;
				break;
			}

			// if it would be smaller than 150 by 100 then limit to 150 by 100
			if (this->RenderWindow->GetSize()[0] == 0)
			{
				this->RenderWindow->SetSize(
					xs < 150 ? 150 : xs, ys < 100 ? 100 : ys);
			}

			if (this->Renderer)
			{
				this->Renderer->ResetCamera();
				this->Renderer->GetActiveCamera()->SetParallelScale(
					xs < 150 ? 75 : (xs - 1) / 2.0);
			}
			this->FirstRender = 0;
		}
	}
	if (this->GetInput())
	{
		this->RenderWindow->Render();
	}
}


const char* ImageViewerForGatedImage::GetWindowName()
{
	return this->RenderWindow->GetWindowName();
}


void ImageViewerForGatedImage::SetOffScreenRendering(int i)
{
	this->RenderWindow->SetOffScreenRendering(i);
}


int ImageViewerForGatedImage::GetOffScreenRendering()
{
	return this->RenderWindow->GetOffScreenRendering();
}


void ImageViewerForGatedImage::SetInputData(vtkImageData *in)
{
	this->WindowLevel->SetInputData(in);
	this->UpdateDisplayExtent();
}


vtkImageData* ImageViewerForGatedImage::GetInput()
{
	return vtkImageData::SafeDownCast(this->WindowLevel->GetInput());
}


vtkInformation* ImageViewerForGatedImage::GetInputInformation()
{
	return this->WindowLevel->GetInputInformation();
}


vtkAlgorithm* ImageViewerForGatedImage::GetInputAlgorithm()
{
	return this->WindowLevel->GetInputAlgorithm();
}


void ImageViewerForGatedImage::SetInputConnection(vtkAlgorithmOutput* input)
{
	this->WindowLevel->SetInputConnection(input);
	this->UpdateDisplayExtent();
}


void ImageViewerForGatedImage::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);

	os << indent << "RenderWindow:\n";
	this->RenderWindow->PrintSelf(os, indent.GetNextIndent());
	os << indent << "Renderer:\n";
	this->Renderer->PrintSelf(os, indent.GetNextIndent());
	os << indent << "ImageActor:\n";
	this->ImageActor->PrintSelf(os, indent.GetNextIndent());
	os << indent << "WindowLevel:\n" << endl;
	this->WindowLevel->PrintSelf(os, indent.GetNextIndent());
	os << indent << "Slice: " << this->Slice << endl;
	os << indent << "SliceOrientation: " << this->SliceOrientation << endl;
	os << indent << "InteractorStyle: " << endl;
	if (this->InteractorStyle)
	{
		os << "\n";
		this->InteractorStyle->PrintSelf(os, indent.GetNextIndent());
	}
	else
	{
		os << "None";
	}
}


void ImageViewerForGatedImage::addActorToRender()
{
	this->GetRenderer()->AddActor(this->ImageActor);
		//this->_renderer->AddViewProp(this->ImageActor);
}


void ImageViewerForGatedImage::addOutputPortToMapper()
{
	if (this->ImageActor && this->WindowLevel)
	{
		this->ImageActor->GetMapper()->SetInputConnection(
			this->WindowLevel->GetOutputPort());
	}
}


unsigned char GetR(unsigned int v)
{
	return v >> 16;
}


unsigned char GetG(unsigned int  v)
{
	v = v << 16;
	v = v >> 24;
	return v;
}


unsigned char GetB(unsigned int v)
{
	v = v << 24;
	v = v >> 24;
	return v;
}


unsigned int RGBA(unsigned char r, unsigned char g, unsigned char b)
{
	unsigned int cc = r;
	cc = cc << 8;
	cc = cc + g;
	cc = cc << 8;
	cc = cc + b;

	return cc;
}


void ImageViewerForGatedImage::SetPetLookupTable()
{
	//this->WindowLevel->GetMTime();
	//vtkScalarsToColors* scalarsToColors = this->WindowLevel->GetLookupTable();
	//CLookupTableWrapper* lookupTableWrapper = dynamic_cast<CLookupTableWrapper*>(this->WindowLevel->GetLookupTable());

	//this->ImageActor->GetProperty()->GetLookupTable();

	//CLookupTableWrapper* lookupTableWrapper = dynamic_cast<CLookupTableWrapper*>(this->ImageActor->GetProperty()->GetLookupTable());
	CLookupTableWrapper* lookupTableWrapper = CLookupTableWrapper::New();

	if (nullptr != lookupTableWrapper)
	{
		int sizeColorRGBA = 256;
		unsigned int colorRGBA[256];

		for (int i = 0; i < 256; i++)
		{
			colorRGBA[i] = RGBA(255 - i, 255 - i, 255 - i);
		}
		
		CHybirdData hybird;
		hybird.Push(colorRGBA);
		hybird.Push(sizeColorRGBA);
		
		hybird.Reset();

		unsigned int* color;
		int number;
		color = (unsigned int*)hybird.Pop(HY_COLORREF_POINTER);
		number = *((int*)hybird.Pop(HY_INT));

		for (int i = 0; i < 256; i++)
		{
			lookupTableWrapper->SetTableValue(i, GetR(color[i]) / 255.0, GetG(color[i]) / 255.0, GetB(color[i]) / 255.0, 1);
		}
		lookupTableWrapper->Build();

		//this->WindowLevel->SetLookupTable(lookupTableWrapper);
		this->ImageActor->GetProperty()->SetLookupTable(lookupTableWrapper);
	}
}


void ImageViewerForGatedImage::Rotate180()
{
	double* focal = this->Renderer->GetActiveCamera()->GetFocalPoint();
	double distance = this->Renderer->GetActiveCamera()->GetDistance();
	double* project = this->Renderer->GetActiveCamera()->GetDirectionOfProjection();
	double* viewUp = this->Renderer->GetActiveCamera()->GetViewUp();

	double up[3];
	for (int i = 0; i < 3; i++)
	{
		up[i] = 0 - viewUp[i];
	}

	double pos[3];
	for (int i = 0; i < 3; i++)
	{
		pos[i] = focal[i] + project[i] * (float)distance;
	}

	this->Renderer->GetActiveCamera()->SetPosition(pos);
	this->Renderer->GetActiveCamera()->SetViewUp(up);
}
