
#include "LookupTableWrapper.h"


#include "vtkBitArray.h"
#include "vtkObjectFactory.h"
#include "vtkMath.h"
#include "vtkMathConfigure.h"
#include <assert.h>

vtkStandardNewMacro(CLookupTableWrapper);

// Construct with range=(0,1); and hsv ranges set up for rainbow color table 
// (from red to blue).
CLookupTableWrapper::CLookupTableWrapper(int sze, int ext)
{
	this->NumberOfColors = sze;
	m_maxIndex=NumberOfColors-1;
	this->Table = vtkUnsignedCharArray::New();
	this->Table->Register(this);
	this->Table->Delete();
	this->Table->SetNumberOfComponents(4);
	this->Table->Allocate(4*sze,4*ext);
	m_tableDatPointer=this->Table->GetPointer(0);

	this->HueRange[0] = 0.0;
	this->HueRange[1] = 0.66667;

	this->SaturationRange[0] = 1.0;
	this->SaturationRange[1] = 1.0;

	this->ValueRange[0] = 1.0;
	this->ValueRange[1] = 1.0;

	this->AlphaRange[0] = 1.0;
	this->AlphaRange[1] = 1.0;
	this->Alpha = 1.0;

	this->NanColor[0] = 0.5;
	this->NanColor[1] = 0.0;
	this->NanColor[2] = 0.0;
	this->NanColor[3] = 1.0;

	

	this->Ramp = VTK_RAMP_SCURVE;
	this->Scale = VTK_SCALE_LINEAR;

	SetTableRange(0,1);

	this->OpaqueFlag=1;

	bLogInverse=false;
}

//----------------------------------------------------------------------------
CLookupTableWrapper::~CLookupTableWrapper()
{
	this->Table->UnRegister(this);
	this->Table = NULL;
	m_tableDatPointer=NULL;
}

//----------------------------------------------------------------------------
// Description:
// Return true if all of the values defining the mapping have an opacity
// equal to 1. Default implementation return true.
int CLookupTableWrapper::IsOpaque()
{
	if(this->OpaqueFlagBuildTime<this->GetMTime())
	{
		int opaque=1;
		if (this->NanColor[3] < 1.0) { opaque = 0; }
		int size=this->Table->GetNumberOfTuples();
		int i=0;
		unsigned char *ptr=this->Table->GetPointer(0);
		while(opaque && i<size)
		{
			opaque=ptr[3]==255;
			ptr+=4;
			++i;
		}
		this->OpaqueFlag=opaque;
		this->OpaqueFlagBuildTime.Modified();
	}
	return this->OpaqueFlag;
}

//----------------------------------------------------------------------------
// Scalar values greater than maximum range value are clamped to maximum
// range value.
void CLookupTableWrapper::SetTableRange(double r[2])
{
	this->SetTableRange(r[0],r[1]);
}

//----------------------------------------------------------------------------
// Set the minimum/maximum scalar values for scalar mapping. Scalar values
// less than minimum range value are clamped to minimum range value.
// Scalar values greater than maximum range value are clamped to maximum
// range value.
void CLookupTableWrapper::SetTableRange(double rmin, double rmax)
{
	if (this->Scale == VTK_SCALE_LOG10 && 
		((rmin > 0 && rmax < 0) || (rmin < 0 && rmax > 0)))
	{
		//vtkErrorMacro("Bad table range for log scale: ["<<rmin<<", "<<rmax<<"]");
		//return;
	}
	if (rmax < rmin)
	{
		vtkErrorMacro("Bad table range: ["<<rmin<<", "<<rmax<<"]");
		return;
	}

	if (this->TableRange[0] == rmin && this->TableRange[1] == rmax)
	{
		return;
	}

	this->TableRange[0] = rmin;
	this->TableRange[1] = rmax;

	if (this->Scale == VTK_SCALE_LOG10)
	{   // handle logarithmic scale
		vtkWrapperLookupTableLogRange(this->TableRange, m_logRange);
		m_shift = -m_logRange[0];
		if (m_logRange[1] <= m_logRange[0])
		{
			m_scale = VTK_DOUBLE_MAX;
		}
		else
		{
			m_scale = (m_maxIndex + 1)/(m_logRange[1] - m_logRange[0]);
		}
	}
	else
	{   // plain old linear
		m_shift = -this->TableRange[0];
		if (this->TableRange[1] <= this->TableRange[0])
		{
			m_scale = VTK_DOUBLE_MAX;
		}
		else
		{
			m_scale = (m_maxIndex + 1)/(this->TableRange[1] - this->TableRange[0]);
		}
	}

	this->Modified();
}

//----------------------------------------------------------------------------
// Have to be careful about the range if scale is logarithmic
void CLookupTableWrapper::SetScale(int scale)
{
	if (this->Scale == scale)
	{
		return;
	}
	this->Scale = scale;
	this->Modified();

	double rmin = this->TableRange[0];
	double rmax = this->TableRange[1];

	if (this->Scale == VTK_SCALE_LOG10 && 
		((rmin > 0 && rmax < 0) || (rmin < 0 && rmax > 0)))
	{
		SetTableRange(1,10);
		vtkErrorMacro("Bad table range for log scale: ["<<rmin<<", "<<rmax<<"], "
			"adjusting to [1, 10]");
		return;
	}
	SetTableRange(TableRange[0],TableRange[1]);
}

//----------------------------------------------------------------------------
// Allocate a color table of specified size.
int CLookupTableWrapper::Allocate(int sz, int ext)
{
	this->NumberOfColors = sz;
	int a = this->Table->Allocate(4*this->NumberOfColors,4*ext);
	m_maxIndex=NumberOfColors-1;
	m_tableDatPointer=this->Table->GetPointer(0);
	this->Modified();
	return a;
}

//----------------------------------------------------------------------------
// Force the lookup table to rebuild
void CLookupTableWrapper::ForceBuild()
{
	int i;
	double hue, sat, val, hinc, sinc, vinc, ainc;
	double rgba[4], alpha;
	unsigned char *c_rgba;

	int maxIndex = this->NumberOfColors - 1;

	if( maxIndex )
	{
		hinc = (this->HueRange[1] - this->HueRange[0])/maxIndex;
		sinc = (this->SaturationRange[1] - this->SaturationRange[0])/maxIndex;
		vinc = (this->ValueRange[1] - this->ValueRange[0])/maxIndex;
		ainc = (this->AlphaRange[1] - this->AlphaRange[0])/maxIndex;
	}
	else
	{
		hinc = sinc = vinc = ainc = 0.0; 
	}

	for (i = 0; i <= maxIndex; i++) 
	{
		hue = this->HueRange[0] + i*hinc;
		sat = this->SaturationRange[0] + i*sinc;
		val = this->ValueRange[0] + i*vinc;
		alpha = this->AlphaRange[0] + i*ainc;

		vtkMath::HSVToRGB(hue, sat, val, &rgba[0], &rgba[1], &rgba[2]);
		rgba[3] = alpha;

		c_rgba = this->Table->WritePointer(4*i,4);

		switch(this->Ramp)
		{
		case VTK_RAMP_SCURVE:
			{
				c_rgba[0] = static_cast<unsigned char> 
					(127.5*(1.0+cos((1.0-static_cast<double>(rgba[0]))*3.141593)));
				c_rgba[1] = static_cast<unsigned char> 
					(127.5*(1.0+cos((1.0-static_cast<double>(rgba[1]))*3.141593)));
				c_rgba[2] = static_cast<unsigned char> 
					(127.5*(1.0+cos((1.0-static_cast<double>(rgba[2]))*3.141593)));
				c_rgba[3] = static_cast<unsigned char> (alpha*255.0);
				/* same code, but with rounding for correctness
				c_rgba[0] = static_cast<unsigned char>
				(127.5*(1.0 + cos((1.0 - rgba[0])*3.141593)) + 0.5);
				c_rgba[1] = static_cast<unsigned char>
				(127.5*(1.0 + cos((1.0 - rgba[1])*3.141593)) + 0.5);
				c_rgba[2] = static_cast<unsigned char>
				(127.5*(1.0 + cos((1.0 - rgba[2])*3.141593)) + 0.5);
				c_rgba[3] = static_cast<unsigned char>(alpha*255.0 + 0.5);
				*/
			}
			break;
		case VTK_RAMP_LINEAR:
			{
				c_rgba[0] = static_cast<unsigned char>(rgba[0]*255.0 + 0.5);
				c_rgba[1] = static_cast<unsigned char>(rgba[1]*255.0 + 0.5);
				c_rgba[2] = static_cast<unsigned char>(rgba[2]*255.0 + 0.5);
				c_rgba[3] = static_cast<unsigned char>(rgba[3]*255.0 + 0.5);
			}
			break;
		case VTK_RAMP_SQRT:
			{
				c_rgba[0] = static_cast<unsigned char>(sqrt(rgba[0])*255.0 + 0.5);
				c_rgba[1] = static_cast<unsigned char>(sqrt(rgba[1])*255.0 + 0.5);
				c_rgba[2] = static_cast<unsigned char>(sqrt(rgba[2])*255.0 + 0.5);
				c_rgba[3] = static_cast<unsigned char>(sqrt(rgba[3])*255.0 + 0.5);
			}
			break;
		default:
			assert("check: impossible case." && 0); // reaching this line is a bug.
			break;
		}
	}
	this->BuildTime.Modified();
}

//----------------------------------------------------------------------------
// Generate lookup table from hue, saturation, value, alpha min/max values. 
// Table is built from linear ramp of each value.
void CLookupTableWrapper::Build()
{
	if (this->Table->GetNumberOfTuples() < 1 ||
		(this->GetMTime() > this->BuildTime && 
		this->InsertTime <= this->BuildTime))
	{
		this->ForceBuild();
	}
}

//----------------------------------------------------------------------------
// get the color for a scalar value
void CLookupTableWrapper::GetColor(double v, double rgb[3])
{
	unsigned char *rgb8 = this->MapValue(v);

	rgb[0] = rgb8[0]/255.0;
	rgb[1] = rgb8[1]/255.0;
	rgb[2] = rgb8[2]/255.0;
}

//----------------------------------------------------------------------------
// get the opacity (alpha) for a scalar value
double CLookupTableWrapper::GetOpacity(double v)
{
	unsigned char *rgb8 = this->MapValue(v);

	return rgb8[3]/255.0;
}

//----------------------------------------------------------------------------
// There is a little more to this than simply taking the log10 of the
// two range values: we do conversion of negative ranges to positive
// ranges, and conversion of zero to a 'very small number'
void CLookupTableWrapper::vtkWrapperLookupTableLogRange(const double range[2], double logRange[2])
{
	logRange[0] = 0;
	logRange[1] = 1;
}

//----------------------------------------------------------------------------
// Apply log to value, with appropriate constraints.
double CLookupTableWrapper::vtkWrapperApplyLogScale(double v, const double range[2], 
							   const double logRange[2])
{
	if(v<range[0]) v=range[0];
	if(v>range[1]) v=range[1];
	if(bLogInverse)
	{
		return 1.0-log10(1+(range[1]-v)/(range[1]-range[0])*9.0);
	}else
	{
		return log10(1+(v-range[0])/(range[1]-range[0])*9.0);
	}
}

//----------------------------------------------------------------------------
// Apply shift/scale to the scalar value v and do table lookup.
unsigned char *CLookupTableWrapper::vtkWrapperLinearLookupMain(double v,
										  unsigned char *table,
										  double maxIndex,
										  double shift, double scale)
{
	double findx = (v + shift)*scale;

	// do not change this code: it compiles into min/max opcodes
	findx = (findx > 0 ? findx : 0);
	findx = (findx < maxIndex ? findx : maxIndex);

	return &table[4*static_cast<unsigned int>(findx)];
}

template<class T>
unsigned char *CLookupTableWrapper::vtkWrapperLinearLookup(
	T v, unsigned char *table, double maxIndex, double shift, double scale,
	unsigned char *vtkNotUsed(nanColor))
{
	return vtkWrapperLinearLookupMain(v, table, maxIndex, shift, scale);
}

//----------------------------------------------------------------------------
// Check for not-a-number when mapping double or float
unsigned char *CLookupTableWrapper::vtkWrapperLinearLookup(
	double v, unsigned char *table, double maxIndex, double shift, double scale,
	unsigned char *nanColor)
{
	// calling isnan() instead of vtkMath::IsNan() improves performance
#ifdef VTK_HAS_ISNAN
	if (isnan(v))
#else
	if (vtkMath::IsNan(v))
#endif
	{
		return nanColor;
	}

	return vtkWrapperLinearLookupMain(v, table, maxIndex, shift, scale);
}

unsigned char *CLookupTableWrapper::vtkWrapperLinearLookup(
	float v, unsigned char *table, double maxIndex, double shift, double scale,
	unsigned char *nanColor)
{
	return vtkWrapperLinearLookup(static_cast<double>(v), table, maxIndex, shift, scale,
		nanColor);
}

//----------------------------------------------------------------------------
void CLookupTableWrapper::GetLogRange(const double range[2], double log_range[2])
{
	vtkWrapperLookupTableLogRange(range, log_range);
}

//----------------------------------------------------------------------------
double CLookupTableWrapper::ApplyLogScale(double v, const double range[2],
									 const double log_range[2])
{
	return vtkWrapperApplyLogScale(v, range, log_range);
}

//----------------------------------------------------------------------------
// Given a scalar value v, return an index into the lookup table
vtkIdType CLookupTableWrapper::GetIndex(double v)
{
	if (this->Scale == VTK_SCALE_LOG10)
	{   // handle logarithmic scale
		
		v = vtkWrapperApplyLogScale(v, this->TableRange, m_logRange);
	}

	// map to an index
	int findx =(int)( (v + m_shift)*m_scale);
	if (findx < 0)
	{
		findx = 0;
	}
	if (findx > m_maxIndex)
	{
		findx = m_maxIndex;
	}
	return findx;
}

//----------------------------------------------------------------------------
// Given a table, set the internal table and set the number of colors.
void CLookupTableWrapper::SetTable(vtkUnsignedCharArray *table)
{
	if (table != this->Table && table != NULL)
	{
		// Check for incorrect arrays.
		if (table->GetNumberOfComponents() != this->Table->GetNumberOfComponents())
		{
			vtkErrorMacro(<<"Number of components in given table (" 
				<< table->GetNumberOfComponents()
				<< ") is incorrect, it should have "
				<< this->Table->GetNumberOfComponents() 
				<< "." );
			return;
		}
		this->Table->UnRegister(this);
		this->Table = table;
		this->Table->Register(this);
		this->NumberOfColors = this->Table->GetNumberOfTuples();
		// If InsertTime is not modified the array will be rebuilt.  So we
		// use the same approach that the SetTableValue function does.
		this->InsertTime.Modified();
		m_maxIndex=NumberOfColors-1;
		m_tableDatPointer=this->Table->GetPointer(0);
		this->Modified();
	}
}

//----------------------------------------------------------------------------
// Given a scalar value v, return an rgba color value from lookup table.
unsigned char *CLookupTableWrapper::MapValue(double v)
{
	if (this->Scale == VTK_SCALE_LOG10)
	{   // handle logarithmic scale
		
		v = vtkWrapperApplyLogScale(v, this->TableRange, m_logRange);
	}

	// map to an index
	int findx =(int)( (v + m_shift)*m_scale);
	if (findx < 0)
	{
		findx = 0;
	}
	if (findx > m_maxIndex)
	{
		findx = m_maxIndex;
	}
	
	return (m_tableDatPointer + 4*findx);
}

//----------------------------------------------------------------------------
template<class T>
void CLookupTableWrapperMapData(CLookupTableWrapper *self, T *input, 
						   unsigned char *output, int length, 
						   int inIncr, int outFormat)
{
	int i = length;
	double *range = self->GetTableRange();
	double maxIndex = self->GetNumberOfColors() - 1;
	double shift, scale;
	unsigned char *table = self->GetPointer(0);
	unsigned char *cptr;
	double alpha;

	unsigned char nanColor[4];
	const double *nanColord = self->GetNanColor();
	for (int c = 0; c < 4; c++)
	{
		double v = nanColord[c];
		if (v < 0.0) { v = 0.0; }
		else if (v > 1.0) { v = 1.0; }
		nanColor[c] = static_cast<unsigned char>(v*255.0 + 0.5);
	}

	if ( (alpha=self->GetAlpha()) >= 1.0 ) //no blending required 
	{
		if (self->GetScale() == VTK_SCALE_LOG10)
		{
			double val;
			double logRange[2];
			self->vtkWrapperLookupTableLogRange(range, logRange);
			shift = -logRange[0];
			if (logRange[1] <= logRange[0])
			{
				scale = VTK_DOUBLE_MAX;
			}
			else
			{
				scale = (maxIndex + 1)/(logRange[1] - logRange[0]);
			}
			if (outFormat == VTK_RGBA)
			{
				while (--i >= 0) 
				{
					val = self->vtkWrapperApplyLogScale(*input, range, logRange);
					cptr = self->vtkWrapperLinearLookup(val, table, maxIndex, shift, scale, nanColor);
					output[0] = cptr[0];
					output[1] = cptr[1];
					output[2] = cptr[2];
					output[3] = cptr[3];
					input += inIncr;
					output += 4;
				}
			}
			else if (outFormat == VTK_RGB)
			{
				while (--i >= 0) 
				{
					val = self->vtkWrapperApplyLogScale(*input, range, logRange);
					cptr = self->vtkWrapperLinearLookup(val, table, maxIndex, shift, scale, nanColor);
					output[0] = cptr[0];
					output[1] = cptr[1];
					output[2] = cptr[2];
					input += inIncr;
					output += 3;
				}
			}
			else if (outFormat == VTK_LUMINANCE_ALPHA)
			{
				while (--i >= 0) 
				{
					val = self->vtkWrapperApplyLogScale(*input, range, logRange);
					cptr = self->vtkWrapperLinearLookup(val, table, maxIndex, shift, scale, nanColor);
					output[0] = static_cast<unsigned char>(cptr[0]*0.30 + cptr[1]*0.59 +
						cptr[2]*0.11 + 0.5);
					output[1] = cptr[3];
					input += inIncr;
					output += 2;
				}
			}
			else // outFormat == VTK_LUMINANCE
			{
				while (--i >= 0) 
				{
					val = self->vtkWrapperApplyLogScale(*input, range, logRange);
					cptr = self->vtkWrapperLinearLookup(val, table, maxIndex, shift, scale, nanColor);
					*output++ = static_cast<unsigned char>(cptr[0]*0.30 + cptr[1]*0.59 + 
						cptr[2]*0.11 + 0.5);
					input += inIncr;
				}
			}
		}//if log scale

		else //not log scale
		{
			shift = -range[0];
			if (range[1] <= range[0])
			{
				scale = VTK_DOUBLE_MAX;
			}
			else
			{
				scale = (maxIndex + 1)/(range[1] - range[0]);
			}

			if (outFormat == VTK_RGBA)
			{
				while (--i >= 0) 
				{
					cptr = self->vtkWrapperLinearLookup(*input, table, maxIndex, shift, scale,
						nanColor);
					output[0] = cptr[0];
					output[1] = cptr[1];
					output[2] = cptr[2];
					output[3] = cptr[3];
					input += inIncr;
					output += 4;
				}
			}
			else if (outFormat == VTK_RGB)
			{
				while (--i >= 0) 
				{
					cptr = self->vtkWrapperLinearLookup(*input, table, maxIndex, shift, scale,
						nanColor);
					output[0] = cptr[0];
					output[1] = cptr[1];
					output[2] = cptr[2];
					input += inIncr;
					output += 3;
				}
			}
			else if (outFormat == VTK_LUMINANCE_ALPHA)
			{
				while (--i >= 0) 
				{
					cptr = self->vtkWrapperLinearLookup(*input, table, maxIndex, shift, scale,
						nanColor);
					output[0] = static_cast<unsigned char>(cptr[0]*0.30 + cptr[1]*0.59 +
						cptr[2]*0.11 + 0.5);
					output[1] = cptr[3];
					input += inIncr;
					output += 2;
				}
			}
			else // outFormat == VTK_LUMINANCE
			{
				while (--i >= 0) 
				{
					cptr = self->vtkWrapperLinearLookup(*input, table, maxIndex, shift, scale,
						nanColor);
					*output++ = static_cast<unsigned char>(cptr[0]*0.30 + cptr[1]*0.59 + 
						cptr[2]*0.11 + 0.5);
					input += inIncr;
				}
			}
		}//if not log lookup
	}//if blending not needed

	else //blend with the specified alpha
	{
		if (self->GetScale() == VTK_SCALE_LOG10)
		{
			double val;
			double logRange[2];
			self->vtkWrapperLookupTableLogRange(range, logRange);
			shift = -logRange[0];
			if (logRange[1] <= logRange[0])
			{
				scale = VTK_DOUBLE_MAX;
			}
			else
			{
				scale = (maxIndex + 1)/(logRange[1] - logRange[0]);
			}
			if (outFormat == VTK_RGBA)
			{
				while (--i >= 0) 
				{
					val = self->vtkWrapperApplyLogScale(*input, range, logRange);
					cptr = self->vtkWrapperLinearLookup(val, table, maxIndex, shift, scale, nanColor);
					output[0] = cptr[0];
					output[1] = cptr[1];
					output[2] = cptr[2];
					output[3] = static_cast<unsigned char>(cptr[3]*alpha + 0.5);
					input += inIncr;
					output += 4;
				}
			}
			else if (outFormat == VTK_RGB)
			{
				while (--i >= 0) 
				{
					val = self->vtkWrapperApplyLogScale(*input, range, logRange);
					cptr = self->vtkWrapperLinearLookup(val, table, maxIndex, shift, scale, nanColor);
					output[0] = cptr[0];
					output[1] = cptr[1];
					output[2] = cptr[2];
					input += inIncr;
					output += 3;
				}
			}
			else if (outFormat == VTK_LUMINANCE_ALPHA)
			{
				while (--i >= 0) 
				{
					val = self->vtkWrapperApplyLogScale(*input, range, logRange);
					cptr = self->vtkWrapperLinearLookup(val, table, maxIndex, shift, scale, nanColor);
					output[0] = static_cast<unsigned char>(cptr[0]*0.30 + cptr[1]*0.59 +
						cptr[2]*0.11 + 0.5);
					output[1] = static_cast<unsigned char>(alpha*cptr[3] + 0.5);
					input += inIncr;
					output += 2;
				}
			}
			else // outFormat == VTK_LUMINANCE
			{
				while (--i >= 0) 
				{
					val = self->vtkWrapperApplyLogScale(*input, range, logRange);
					cptr = self->vtkWrapperLinearLookup(val, table, maxIndex, shift, scale, nanColor);
					*output++ = static_cast<unsigned char>(cptr[0]*0.30 + cptr[1]*0.59 + 
						cptr[2]*0.11 + 0.5);
					input += inIncr;
				}
			}
		}//log scale with blending

		else //no log scale with blending
		{
			shift = -range[0];
			if (range[1] <= range[0])
			{
				scale = VTK_DOUBLE_MAX;
			}
			else
			{
				scale = (maxIndex + 1)/(range[1] - range[0]);
			}

			if (outFormat == VTK_RGBA)
			{
				while (--i >= 0) 
				{
					cptr = self->vtkWrapperLinearLookup(*input, table, maxIndex, shift, scale,
						nanColor);
					output[0] = cptr[0];
					output[1] = cptr[1];
					output[2] = cptr[2];
					output[3] = static_cast<unsigned char>(cptr[3]*alpha + 0.5);
					input += inIncr;
					output += 4;
				}
			}
			else if (outFormat == VTK_RGB)
			{
				while (--i >= 0) 
				{
					cptr = self->vtkWrapperLinearLookup(*input, table, maxIndex, shift, scale,
						nanColor);
					output[0] = cptr[0];
					output[1] = cptr[1];
					output[2] = cptr[2];
					input += inIncr;
					output += 3;
				}
			}
			else if (outFormat == VTK_LUMINANCE_ALPHA)
			{
				while (--i >= 0) 
				{
					cptr = self->vtkWrapperLinearLookup(*input, table, maxIndex, shift, scale,
						nanColor);
					output[0] = static_cast<unsigned char>(cptr[0]*0.30 + cptr[1]*0.59 +
						cptr[2]*0.11 + 0.5);
					output[1] = static_cast<unsigned char>(cptr[3]*alpha + 0.5);
					input += inIncr;
					output += 2;
				}
			}
			else // outFormat == VTK_LUMINANCE
			{
				while (--i >= 0) 
				{
					cptr = self->vtkWrapperLinearLookup(*input, table, maxIndex, shift, scale,
						nanColor);
					*output++ = static_cast<unsigned char>(cptr[0]*0.30 + cptr[1]*0.59 + 
						cptr[2]*0.11 + 0.5);
					input += inIncr;
				}
			}
		}//no log scale
	}//alpha blending
}

//----------------------------------------------------------------------------
void CLookupTableWrapper::MapScalarsThroughTable2(void *input, 
											 unsigned char *output,
											 int inputDataType, 
											 int numberOfValues,
											 int inputIncrement,
											 int outputFormat)
{
	switch (inputDataType)
	{
	case VTK_BIT:
		{
			vtkIdType i, id;
			vtkBitArray *bitArray = vtkBitArray::New();
			bitArray->SetVoidArray(input,numberOfValues,1);
			vtkUnsignedCharArray *newInput = vtkUnsignedCharArray::New();
			newInput->SetNumberOfValues(numberOfValues);
			for (id=i=0; i<numberOfValues; i++, id+=inputIncrement)
			{
				newInput->SetValue(i, bitArray->GetValue(id));
			}
			CLookupTableWrapperMapData(this,
				static_cast<unsigned char*>(newInput->GetPointer(0)),
				output,numberOfValues,
				inputIncrement,outputFormat);
			newInput->Delete();
			bitArray->Delete();
		}
		break;

		vtkTemplateMacro(
			CLookupTableWrapperMapData(this,static_cast<VTK_TT*>(input),output,
			numberOfValues,inputIncrement,outputFormat)
			);
	default:
		vtkErrorMacro(<< "MapImageThroughTable: Unknown input ScalarType");
		return;
	}
}  

//----------------------------------------------------------------------------
// Specify the number of values (i.e., colors) in the lookup
// table. This method simply allocates memory and prepares the table
// for use with SetTableValue(). It differs from Build() method in
// that the allocated memory is not initialized according to HSVA ramps.
void CLookupTableWrapper::SetNumberOfTableValues(vtkIdType number)
{
	if (this->NumberOfColors == number)
	{
		return;
	}
	this->Modified();
	this->NumberOfColors = number;
	m_maxIndex=NumberOfColors-1;
	this->Table->SetNumberOfTuples(number);
	m_tableDatPointer=this->Table->GetPointer(0);
}

//----------------------------------------------------------------------------
// Directly load color into lookup table. Use [0,1] double values for color
// component specification. Make sure that you've either used the
// Build() method or used SetNumberOfTableValues() prior to using this method.
void CLookupTableWrapper::SetTableValue(vtkIdType indx, double rgba[4])
{
	// Check the index to make sure it is valid
	if (indx < 0)
	{
		vtkErrorMacro("Can't set the table value for negative index " << indx);
		return;
	}
	if (indx >= this->NumberOfColors)
	{
		vtkErrorMacro("Index " << indx << 
			" is greater than the number of colors " << 
			this->NumberOfColors);
		return;
	}

	unsigned char *_rgba = this->Table->WritePointer(4*indx,4);

	_rgba[0] = static_cast<unsigned char>(rgba[0]*255.0 + 0.5);
	_rgba[1] = static_cast<unsigned char>(rgba[1]*255.0 + 0.5);
	_rgba[2] = static_cast<unsigned char>(rgba[2]*255.0 + 0.5);
	_rgba[3] = static_cast<unsigned char>(rgba[3]*255.0 + 0.5);

	this->InsertTime.Modified();
	this->Modified();
}

//----------------------------------------------------------------------------
// Directly load color into lookup table. Use [0,1] double values for color 
// component specification.
void CLookupTableWrapper::SetTableValue(vtkIdType indx, double r, double g, double b, 
								   double a)
{
	double rgba[4];
	rgba[0] = r; rgba[1] = g; rgba[2] = b; rgba[3] = a;
	this->SetTableValue(indx,rgba);
}

//----------------------------------------------------------------------------
// Return a rgba color value for the given index into the lookup Table. Color
// components are expressed as [0,1] double values.
void CLookupTableWrapper::GetTableValue(vtkIdType indx, double rgba[4])
{
	unsigned char *_rgba;

	indx = (indx < 0 ? 0 : (indx >= this->NumberOfColors ? 
		this->NumberOfColors-1 : indx));

	_rgba = this->Table->GetPointer(indx*4);

	rgba[0] = _rgba[0]/255.0;
	rgba[1] = _rgba[1]/255.0;
	rgba[2] = _rgba[2]/255.0;
	rgba[3] = _rgba[3]/255.0;
}

// Return a rgba color value for the given index into the lookup table. Color
// components are expressed as [0,1] double values.
double *CLookupTableWrapper::GetTableValue(vtkIdType indx)
{
	this->GetTableValue(indx, this->RGBA);
	return this->RGBA;
}

//----------------------------------------------------------------------------
void CLookupTableWrapper::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);

	os << indent << "TableRange: (" << this->TableRange[0] << ", "
		<< this->TableRange[1] << ")\n";
	os << indent << "Scale: "
		<< (this->Scale == VTK_SCALE_LOG10 ? "Log10\n" : "Linear\n");
	os << indent << "HueRange: (" << this->HueRange[0] << ", "
		<< this->HueRange[1] << ")\n";
	os << indent << "SaturationRange: (" << this->SaturationRange[0] << ", "
		<< this->SaturationRange[1] << ")\n";
	os << indent << "ValueRange: (" << this->ValueRange[0] << ", "
		<< this->ValueRange[1] << ")\n";
	os << indent << "AlphaRange: (" << this->AlphaRange[0] << ", "
		<< this->AlphaRange[1] << ")\n";
	os << indent << "NanColor: (" << this->NanColor[0] << ", "
		<< this->NanColor[1] << ", " << this->NanColor[2] << ", "
		<< this->NanColor[3] << ")\n";
	os << indent << "NumberOfTableValues: "
		<< this->GetNumberOfTableValues() << "\n";
	os << indent << "NumberOfColors: " << this->NumberOfColors << "\n";
	os << indent << "Ramp: "
		<< (this->Ramp == VTK_RAMP_SCURVE ? "SCurve\n" : "Linear\n");
	os << indent << "InsertTime: " <<this->InsertTime.GetMTime() << "\n";
	os << indent << "BuildTime: " <<this->BuildTime.GetMTime() << "\n";
	os << indent << "Table: ";
	if( this->Table )
	{
		this->Table->PrintSelf(os << "\n", indent.GetNextIndent());
	}
	else
	{
		// Should not happen
		os << "(none)\n";
	}
}

//----------------------------------------------------------------------------
void CLookupTableWrapper::DeepCopy(vtkScalarsToColors *obj)
{
	if (!obj)
	{
		return;
	}

	CLookupTableWrapper *lut = CLookupTableWrapper::SafeDownCast(obj);

	if (!lut)
	{
		vtkErrorMacro("Cannot DeepCopy a " << obj->GetClassName()
			<< " into a CLookupTableWrapper.");
		return;
	}

	this->Scale               = lut->Scale;
	this->TableRange[0]       = lut->TableRange[0];
	this->TableRange[1]       = lut->TableRange[1];
	this->HueRange[0]         = lut->HueRange[0];
	this->HueRange[1]         = lut->HueRange[1];
	this->SaturationRange[0]  = lut->SaturationRange[0];
	this->SaturationRange[1]  = lut->SaturationRange[1];
	this->ValueRange[0]       = lut->ValueRange[0];
	this->ValueRange[1]       = lut->ValueRange[1];
	this->AlphaRange[0]       = lut->AlphaRange[0];
	this->AlphaRange[1]       = lut->AlphaRange[1];
	this->NumberOfColors      = lut->NumberOfColors;
	this->Ramp                = lut->Ramp;
	this->InsertTime          = lut->InsertTime;
	this->BuildTime           = lut->BuildTime;
	this->Table->DeepCopy(lut->Table);

	m_maxIndex=NumberOfColors-1;
	m_tableDatPointer=this->Table->GetPointer(0);
	SetTableRange(TableRange[0],TableRange[1]);
	this->Superclass::DeepCopy(obj);
}

//----------------------------------------------------------------------------
vtkIdType CLookupTableWrapper::GetNumberOfAvailableColors()
{
	return this->Table->GetNumberOfTuples();
}
