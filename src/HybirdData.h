#ifndef __NS_CHybirdData_H__
#define __NS_CHybirdData_H__

enum HYBIRDDATA_TYPE
{
	HY_COLORREF_POINTER,
	HY_INT,
	HY_UINT64,
	HY_DOUBLE,
	HY_CHAR_POINTER,
	HY_DOUBLE_POINTER
};

class  CHybirdData
{
public:
	CHybirdData(void);
	~CHybirdData(void);

	void Push(char* value);
	void Push(double* value);
	void Push(unsigned int* color);
	void Push(int value);
	void Push(unsigned long long value);
	void Push(double value);

	char* Pop(HYBIRDDATA_TYPE type);

	void Reset();
private:
	
	char* m_pData;
	char* m_pCurData;
};

#endif