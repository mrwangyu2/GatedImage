#include "MyVtkInteractorStyleImage.h"
vtkStandardNewMacro(myVtkInteractorStyleImage);

// helper class to format slice status message
class StatusMessage {
public:
	static std::string Format(int slice, int maxSlice) {
		std::stringstream tmp;
		tmp << "Slice Number  " << slice + 1 << "/" << maxSlice + 1;
		return tmp.str();
	}
};

myVtkInteractorStyleImage::myVtkInteractorStyleImage()
{
}


myVtkInteractorStyleImage::~myVtkInteractorStyleImage()
{
}


void myVtkInteractorStyleImage::SetImageViewer(ImageViewerForGatedImage* imageViewer) {
	_ImageViewer = imageViewer;
	_MinSlice = imageViewer->GetSliceMin();
	_MaxSlice = imageViewer->GetSliceMax();
	_Slice = _MinSlice;
	cout << "Slicer: Min = " << _MinSlice << ", Max = " << _MaxSlice << std::endl;
}


void myVtkInteractorStyleImage::SetStatusMapper(vtkTextMapper* statusMapper) {
	_StatusMapper = statusMapper;
}


void myVtkInteractorStyleImage::MoveSliceForward() {
	if (_Slice < _MaxSlice) {
		_Slice += 1;
		cout << "MoveSliceForward::Slice = " << _Slice << std::endl;
		//_ImageViewer->SetSlice(_Slice);
		std::string msg = StatusMessage::Format(_Slice, _MaxSlice);
		//_StatusMapper->SetInput(msg.c_str());
		_ImageViewer->Render();
	}
	else
	{
		_Slice = 0;
	}
}


void myVtkInteractorStyleImage::MoveSliceBackward() {
	if (_Slice > _MinSlice) {
		_Slice -= 1;
		cout << "MoveSliceBackward::Slice = " << _Slice << std::endl;
		//_ImageViewer->SetSlice(_Slice);
		std::string msg = StatusMessage::Format(_Slice, _MaxSlice);
		//_StatusMapper->SetInput(msg.c_str());
		_ImageViewer->Render();
	}
}


void myVtkInteractorStyleImage::OnKeyDown() {
	std::string key = this->GetInteractor()->GetKeySym();
	if (key.compare("Up") == 0) {
		//cout << "Up arrow key was pressed." << endl;
		MoveSliceForward();
	}
	else if (key.compare("Down") == 0) {
		//cout << "Down arrow key was pressed." << endl;
		MoveSliceBackward();
	}
	// forward event
	vtkInteractorStyleImage::OnKeyDown();
}


void myVtkInteractorStyleImage::OnMouseWheelForward() {
	//std::cout << "Scrolled mouse wheel forward." << std::endl;
	//MoveSliceForward();
	// don't forward events, otherwise the image will be zoomed 
	// in case another interactorstyle is used (e.g. trackballstyle, ...)
	// vtkInteractorStyleImage::OnMouseWheelForward();

	std::cout << "display next actor" << std::endl;
	this->_gatedImageViewer->displayNextActor();
	//this->_ImageViewer->displayNextActor();
}


void myVtkInteractorStyleImage::OnMouseWheelBackward() {
	//std::cout << "Scrolled mouse wheel backward." << std::endl;
	if (_Slice > _MinSlice) {
		//MoveSliceBackward();
	}
	// don't forward events, otherwise the image will be zoomed 
	// in case another interactorstyle is used (e.g. trackballstyle, ...)
	// vtkInteractorStyleImage::OnMouseWheelBackward();
}

void myVtkInteractorStyleImage::setGatedImageViewer(GatedImageViewer* viewer)
{
	this->_gatedImageViewer = viewer;
}

void myVtkInteractorStyleImage::OnTimer()
{
	//_ImageViewer->displayNextActor();
	//MoveSliceForward();
	this->_gatedImageViewer->displayNextActor();
}


