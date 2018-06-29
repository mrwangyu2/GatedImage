
#include "HybirdData.h"


CHybirdData::CHybirdData(void)
{
	m_pData=new char[256];
	m_pCurData=nullptr;
	Reset();
}


CHybirdData::~CHybirdData(void)
{
	delete []m_pData;
}

void CHybirdData::Reset()
{
	m_pCurData=m_pData;
}
void CHybirdData::Push(unsigned int* color)
{
	unsigned int** ptr=(unsigned int**)m_pCurData;
	unsigned int** tmp=ptr;
	*tmp=color;
	m_pCurData+=sizeof(unsigned int*);
}
void CHybirdData::Push(int value)
{
	int* tmp=(int*)m_pCurData;
	*tmp=value;
	m_pCurData+=sizeof(int);
}
void CHybirdData::Push(unsigned long long value)
{
	unsigned long long* tmp=(unsigned long long*)m_pCurData;
	*tmp=value;
	m_pCurData+=sizeof(unsigned long long);
}
void CHybirdData::Push(double value)
{
	double* tmp=(double*)m_pCurData;
	*tmp=value;
	m_pCurData+=sizeof(double);
}

void CHybirdData::Push(char* value)
{
	char** ptr=(char**)m_pCurData;
	char** tmp=ptr;
	*tmp=value;
	m_pCurData+=sizeof(char*);
}
void CHybirdData::Push(double* value)
{
	double** ptr=(double**)m_pCurData;
	double** tmp=ptr;
	*tmp=value;
	m_pCurData+=sizeof(double*);
}
char* CHybirdData::Pop(HYBIRDDATA_TYPE type)
{
	char* tmp=m_pCurData;

	unsigned int** colorPtr;
	double** doublePtr;
	char** charPtr;
	switch (type)
	{
	case HY_DOUBLE_POINTER:
		doublePtr=(double**)m_pCurData;
		tmp=(char*)(*doublePtr);
		m_pCurData+=sizeof(double*);
		break;
	case HY_CHAR_POINTER:
		charPtr=(char**)m_pCurData;
		tmp=(char*)(*charPtr);
		m_pCurData+=sizeof(char*);
		break;
	case HY_COLORREF_POINTER:
		colorPtr=(unsigned int**)m_pCurData;
		tmp=(char*)(*colorPtr);
		m_pCurData+=sizeof(unsigned int*);
		break;
	case HY_INT:
		m_pCurData+=sizeof(int);
		break;
	case HY_UINT64:
		m_pCurData+=sizeof(unsigned long long);
		break;
	case HY_DOUBLE:
		m_pCurData+=sizeof(double);
		break;
	default:
		break;
	}
	return (char*)tmp;
}

