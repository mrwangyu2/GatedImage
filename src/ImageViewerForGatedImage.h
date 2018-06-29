#ifndef ImageViewerForGatedImage_h
#define ImageViewerForGatedImage_h

#include <vtkInteractionImageModule.h> // For export macro
#include <vtkImageProperty.h>
#include "vtkObject.h"
#include "LookupTableWrapper.h"
#include "HybirdData.h"

class vtkAlgorithm;
class vtkAlgorithmOutput;
class vtkImageActor;
class vtkImageData;
class vtkImageMapToWindowLevelColors;
class vtkInformation;
class vtkInteractorStyleImage;
class vtkRenderWindow;
class vtkRenderer;
class vtkRenderWindowInteractor;
class vtkImageProperty;

class  ImageViewerForGatedImage : public vtkObject
{
public:
	static ImageViewerForGatedImage *New();
	vtkTypeMacro(ImageViewerForGatedImage, vtkObject);
	void PrintSelf(ostream& os, vtkIndent indent);

	/**
	* Get the name of rendering window.
	*/
	virtual const char *GetWindowName();

	/**
	* Render the resulting image.
	*/
	virtual void Render(void);

	//@{
	/**
	* Set/Get the input image to the viewer.
	*/
	virtual void SetInputData(vtkImageData *in);
	virtual vtkImageData *GetInput();
	virtual void SetInputConnection(vtkAlgorithmOutput* input);
	//@}

	/**
	* Set/get the slice orientation
	*/

	enum
	{
		SLICE_ORIENTATION_YZ = 0,
		SLICE_ORIENTATION_XZ = 1,
		SLICE_ORIENTATION_XY = 2
	};

	vtkGetMacro(SliceOrientation, int);
	virtual void SetSliceOrientation(int orientation);
	virtual void SetSliceOrientationToXY()
	{
		this->SetSliceOrientation(ImageViewerForGatedImage::SLICE_ORIENTATION_XY);
	};
	virtual void SetSliceOrientationToYZ()
	{
		this->SetSliceOrientation(ImageViewerForGatedImage::SLICE_ORIENTATION_YZ);
	};
	virtual void SetSliceOrientationToXZ()
	{
		this->SetSliceOrientation(ImageViewerForGatedImage::SLICE_ORIENTATION_XZ);
	};


	void ImageViewerForGatedImage::Rotate180();

	//@{
	/**
	* Set/Get the current slice to display (depending on the orientation
	* this can be in X, Y or Z).
	*/
	vtkGetMacro(Slice, int);
	virtual void SetSlice(int s);
	//@}

	/**
	* Update the display extent manually so that the proper slice for the
	* given orientation is displayed. It will also try to set a
	* reasonable camera clipping range.
	* This method is called automatically when the Input is changed, but
	* most of the time the input of this class is likely to remain the same,
	* i.e. connected to the output of a filter, or an image reader. When the
	* input of this filter or reader itself is changed, an error message might
	* be displayed since the current display extent is probably outside
	* the new whole extent. Calling this method will ensure that the display
	* extent is reset properly.
	*/
	virtual void UpdateDisplayExtent();

	//@{
	/**
	* Return the minimum and maximum slice values (depending on the orientation
	* this can be in X, Y or Z).
	*/
	virtual int GetSliceMin();
	virtual int GetSliceMax();
	virtual void GetSliceRange(int range[2])
	{
		this->GetSliceRange(range[0], range[1]);
	}
	virtual void GetSliceRange(int &min, int &max);
	virtual int* GetSliceRange();
	//@}

	//@{
	/**
	* Set window and level for mapping pixels to colors.
	*/
	virtual double GetColorWindow();
	virtual double GetColorLevel();
	virtual void SetColorWindow(double s);
	virtual void SetColorLevel(double s);
	//@}

	//@{
	/**
	* These are here when using a Tk window.
	*/
	virtual void SetDisplayId(void *a);
	virtual void SetWindowId(void *a);
	virtual void SetParentId(void *a);
	//@}

	//@{
	/**
	* Set/Get the position in screen coordinates of the rendering window.
	*/
	virtual int* GetPosition();
	virtual void SetPosition(int a, int b);
	virtual void SetPosition(int a[2]) { this->SetPosition(a[0], a[1]); }
	//@}

	//@{
	/**
	* Set/Get the size of the window in screen coordinates in pixels.
	*/
	virtual int* GetSize();
	virtual void SetSize(int a, int b);
	virtual void SetSize(int a[2]) { this->SetSize(a[0], a[1]); }
	//@}

	//@{
	/**
	* Get the internal render window, renderer, image actor, and
	* image map instances.
	*/
	vtkGetObjectMacro(RenderWindow, vtkRenderWindow);
	vtkGetObjectMacro(Renderer, vtkRenderer);
	vtkGetObjectMacro(ImageActor, vtkImageActor);
	//vtkImageActor* GetImageActor(){ return this->ImageActor; };
	vtkGetObjectMacro(WindowLevel, vtkImageMapToWindowLevelColors);
	vtkGetObjectMacro(InteractorStyle, vtkInteractorStyleImage);
	//@}

	//@{
	/**
	* Set your own renderwindow and renderer
	*/
	virtual void SetRenderWindow(vtkRenderWindow *arg);
	virtual void SetRenderer(vtkRenderer *arg);
	//@}

	/**
	* Attach an interactor for the internal render window.
	*/
	virtual void SetupInteractor(vtkRenderWindowInteractor*);

	//@{
	/**
	* Create a window in memory instead of on the screen. This may not
	* be supported for every type of window and on some windows you may
	* need to invoke this prior to the first render.
	*/
	virtual void SetOffScreenRendering(int);
	virtual int GetOffScreenRendering();
	vtkBooleanMacro(OffScreenRendering, int);
	//@}

protected:
	ImageViewerForGatedImage();
	~ImageViewerForGatedImage();

	virtual void InstallPipeline();
	virtual void UnInstallPipeline();

	vtkImageMapToWindowLevelColors  *WindowLevel;
	vtkRenderWindow                 *RenderWindow;
	vtkRenderer                     *Renderer;
	vtkImageActor                   *ImageActor;
	vtkRenderWindowInteractor       *Interactor;
	vtkInteractorStyleImage         *InteractorStyle;

	int SliceOrientation;
	int FirstRender;
	int Slice;

	virtual void UpdateOrientation();

	vtkAlgorithm* GetInputAlgorithm();
	vtkInformation* GetInputInformation();

	friend class ImageViewerForGatedImageCallback;

private:
	ImageViewerForGatedImage(const ImageViewerForGatedImage&) VTK_DELETE_FUNCTION;
	void operator=(const ImageViewerForGatedImage&)VTK_DELETE_FUNCTION;
public:
	void addActorToRender();
	void addOutputPortToMapper();
	void setImageViewer();
	void SetPetLookupTable();
};

#endif
