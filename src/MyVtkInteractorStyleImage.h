
#ifndef MyVtkInteractorStyleImage_h
#define MyVtkInteractorStyleImage_h

#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkImageActor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkTextMapper.h>

#include "ImageViewerForGatedImage.h"
#include "GatedImageViewer.h"

class GatedImageViewer;

class myVtkInteractorStyleImage : public vtkInteractorStyleImage
{
public:
	myVtkInteractorStyleImage();
	~myVtkInteractorStyleImage();

	static myVtkInteractorStyleImage* New();
	vtkTypeMacro(myVtkInteractorStyleImage, vtkInteractorStyleImage);

protected:
	ImageViewerForGatedImage* _ImageViewer;
	vtkTextMapper* _StatusMapper;
	GatedImageViewer* _gatedImageViewer;

	int _Slice;
	int _MinSlice;
	int _MaxSlice;

public:
	void setGatedImageViewer(GatedImageViewer* viewer);
	void SetImageViewer(ImageViewerForGatedImage* imageViewer);
	void SetStatusMapper(vtkTextMapper* statusMapper);


protected:
	void MoveSliceForward();
	void MoveSliceBackward();
	virtual void OnKeyDown();
	virtual void OnMouseWheelForward();
	virtual void OnMouseWheelBackward();

	virtual void OnTimer();
};

#endif
