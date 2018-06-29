// MyMath.h: wrapper class for basic math calculation 
//
//////////////////////////////////////////////////////////////////////


#if !defined(AFX_MYMATH_H__4AC3D190_94EC_11D2_B21D_00C04FBEDB8F__INCLUDED_)
#define AFX_MYMATH_H__4AC3D190_94EC_11D2_B21D_00C04FBEDB8F__INCLUDED_

#define TEST_DEMO
#define MAX_TRINUM    500000
#define EPSLON        1.0e-6


struct PT_RPY {
	double x;
	double y;
	double z;
	double roll;
	double pitch;
	double yaw;
} ;

class CStdMath 
{
public:
	void AxisRotation(double xyz[], double rot_ang, double nx, double ny, double nz, double xyz_out[]);
	void AxisRotation(double xyz[3], double rot_ang, double axis[6], double xyz_out[3]);
    void LineFitting2D(int     numPnt, 
                       double *coefficients, 
                       double *x, 
                       double *y, 
                       double *stdev);
	void MinEigenvector(double **A, int nrow, double *eigenVector);
	void PlaneFitting(int num_point, double  *a[3], double *nx, double *ny, double *nz, double *d, double *stdDev, double *maxDev);
	void SampleCovarianceAndMean(int numSample, int dim, double *sampleData[], double *mean, double **cov);
	void LineFitting3D(int num_point, double *a[], double *kx, double *ky, double *kz, double *std, double *maxErr);
	double DistanceFromPointToLine(double x, double y, double z, double *lineDirection, double x0, double y0, double z0);

	void MaxEigenvector(double **A, int nrow, double *eigenVector);
	void LineFitting3D(int num_point, double *a[], 
                        double *kx, double *ky, double *kz, 
                        double *x0, double *y0, double *z0, 
                        double *std, double *maxErr);
	void EllipseFitting3D(int numPoint, 
	                      double *ellipse[3], 
						  int fittedPoint, 
						  double *newEllipse[3]);
	void EllipseFitting3D(int     numPoint, 
							   double *ellipse[3], 
							   int     fittedPoint,
							   double *newEllipse[3],
							   double *cx,
							   double *cy,
							   double *cz,
							   double *xAxi, 
							   double *yAxi, 
							   double *Theta);
	void Transfer2Dto3D(int numPoint, double *Point2D[2], double* Point3D[3]);
	void EllipseFitting_2D(int numPoint, double *ellipse[2], int fittedPoint, double *newEllipse[2]);
	void EllipseFitting_2D(int     numPoint, 
	                       double *ellipse[2], 
						   int     fittedPoint,
						   double *newEllipse[2],
						   double  *xc, 
						   double  *yc,  
						   double  *xAxi, 
						   double  *yAxi, 
						   double  *Theta);
 
	void CopyMatrix(double m1[], int startrow1,int startcol1, double m2[],
					int startrow2, int startcol2, int rows, int cols);

	void MultiplyMatrixFloat(double m1[], int row1, int col1, double m2[], 
		int col2, double m3[]);

	void MultiplyMatrix(double m1[], int row1, int col1, double m2[], 
						int col2, double m3[]);

	void LinearLeastSquare(int numRow, int numCol, double **A, double *b, double *x, double *stdDev = 0, double *maxError = 0);
 
	void NonLinearLeastSquare(double** x, double y[], double sig[], int ndata, 
			double a[], int ma, int lista[], int mfit, double** covar, 
			double** alpha, double* chisq, double* alamda, int index);
	void NonLinearLeastSquare(double** x, double y[], double sig[], int ndata, 
			double a[], int ma, int lista[], int mfit, int index);

	void NonLinearLeastSquare_calib(double** x, double y[], double sig[], int ndata, 
			double a[], int ma, int lista[], int mfit, int index, double *pXYZw);

	void NonLinearLeastSquare_calib_I(double** x, double y[], double sig[], int ndata, 
			double a[], int ma, int lista[], int mfit, int index, double *pXYZw);

//	void NonLinearLeastSquare_modelfit(double** x, double y[], double sig[], int ndata, 
//		double a[], double *fiterror, int ma, int lista[], int mfit, int num_line, double **line_pars);

	void NonLinearLeastSquare_modelfit(double** x, double y[], double sig[], int ndata, 
		double a[], double fiterror[2], int ma, int lista[], int mfit, int num_line, double **line_pars);
		
	PT_RPY PTIJK_to_PTRPY( double P1[], double d[], double n[]);

	void MatrixInverse(int numCol, double **A);
	double Atan360(double num, double den);

	void MatrixMultiply33(double m1[3][3],double m2[3][3], double m3[3][3]);
	void MatrixMultiply44(double m1[4][4],double m2[4][4], double m3[4][4]);
	void MatrixMultiplyVector44(double m1[4][4],double v1[4], double v2[4]);
	void MatrixMultiplyVector33(double m1[3][3],double v1[3], double v2[3]);

    // The following functions are used to support 3d circle fitting 
    void CircleFitting(int num_point, double *cir[3], double *xc, double *yc, double *zc, double *R, double *Nx, double *Ny, double *Nz);
    void CircleFitting(int num_point, double *cir[3], double *xc, double *yc, double *zc, double *R, double *Nx, double *Ny, double *Nz, double *stdDev, double *maxDev);
	void GetPlaneNormal(int num_point, double  *a[3], double *nx, double *ny, double *nz, double *d);
	void GetPlaneNormal(int num_point, double  *a[3], double *nx, double *ny, double *nz, double *d, double *stdDev, double *maxDev);
	void GetSphere(int num_point, double  *p[3], double init_par[4], double *x0, double *y0, double *z0, double *R, double *stdDev, double *maxDev);
	void GetSphere_fixRadius(int num_point, double  *p[3], double init_par[4], double *x0, double *y0, double *z0, double *R, double *stdDev, double *maxDev);
	void Get2DCircleFit_NL(int num_point, double  *p[2], double init_par[3], double *x0, double *y0, double *R, double *stdDev, double *maxDev, bool bFixedRadius);
	void GetCylinder(int num_point, double  *p[3], double *nx, double *ny, double *nz, double *x0, double *y0, double *z0, double *R, double *stdDev, double *maxDev);
    void GetXShiftParams(int num_point, double  *p[5], double init_par[3], double *k1, double *k2, double *k3, double *stdDev, double *maxDev);
	void GetXShiftParams(int num_point, double  *p[5], double init_par[6], double *k1, double *k2, double *k3, double *x0, double *y0, double *z0, double *stdDev, double *maxDev);
    void GetXShiftParams2(int num_point, double  *p[5], double init_par[6], double *k1, double *k2, double *k3, double *x0, double *y0, double *z0, double *stdDev, double *maxDev);
	void GetXShiftAndRotateParams(int num_point, double  *p[6], double init_par[8], double *k1, double *k2, double *k3, double *x0, double *y0, double *z0, double *rot_x, double *rot_z, double *stdDev, double *maxDev);
	void GetXShiftAndRotateParams2(int num_point, double  *p[6], double init_par[8], double *k1, double *k2, double *k3, double *x0, double *y0, double *z0, double *rot_x, double *rot_z, double *stdDev, double *maxDev);
    void GetXShiftAndRotateParams3(int num_point, double  *p[6], double init_par[8], double *k1, double *k2, double *k3, double *x0, double *y0, double *z0, double *rot_x, double *rot_z, double *stdDev, double *maxDev);
	void GetXShiftAndRotateParams4(int num_point, double  *p[6], double init_par[9], double *k1, double *k2, double *k3, double *x0, double *y0, double *z0, double *rot_x, double *rot_z, double *del_theta, double *stdDev, double *maxDev);
	void GetXShiftAndRotateParams5(int num_point, double  *p[6], double init_par[8], double *k1, double *k2, double *k3, double *x0, double *y0, double *z0, double *rot_x, double *rot_z, double *stdDev, double *maxDev);
	void CircleFitting_2D (int num_point, double *cir[2], double *xc, double *yc, double *R);
	void CircleFitting_2D (int num_point, double *cir[2], double *xc, double *yc, double *R, double *stdDev, double *maxDev);
	void CircleFitting_3D_NL(int num_point, double *cir[3], double *xc, double *yc, double *zc, double *R, double *Nx, double *Ny, double *Nz, double *stdDev, double *maxDev, bool bFixedRadius);
		
	
    int CalculateRotAxis(double m[3][3], double n[3]);
    
	int  CalculateCoordinatefrom2Vectors(double vector_y[3], double orig[3], double vector_x[3], double t[4][4]);
	int CalculateRotCenter(double p1[3], double p2[3], double n[3], double c[3], double beta);
	double CalculateRotAxisBy2Mat(double m1[3][3], double m2[3][3], double n[3]);
  
	// Basic transform
    void RotateYAxis(double xyz_in[3], double angle, double xyz_out[3]);

    // basic geometry (projection)
	void PointProjection(double x1, double y1, double z1, double nx, double ny, double nz, double d, double *x2, double *y2, double *z2);
    void PointLineProjection(double x0, double y0, double z0, double x1, double y1, double z1, double nx, double ny, double nz,double *xp, double *yp, double *zp);
	void PointLineDistance(double x0, double y0, double z0, double x1, double y1, double z1, double nx, double ny, double nz,double *d);
	void PointLineDistance_2D(double x0, double y0, double line_mode, double a, double b, double *d);
	void LinePlaneIntersection(double xl, double yl, double zl, double nlx, double nly, double nlz, double npx, double npy, double npz, double d, double *xi, double *yi, double *zi);
	//  line1_par = [x1, y1, z1, nx1, ny1, nz1]
	//  line2_par = [x2, y2, z2, nx2, ny2, nz2]
	void Line2LineIntersect(double line1_par[6], double line2_par[6], double intersect[3], double *err);
		
	bool TransposeMatrix(double  m1[], int row, int col, double m2[]);
	bool CalculateTransform(int iNum, double T1[], double T2[], double M[16]);
	bool CalculateTransform2(int iNum, double T1[], double T2[], double M[16]);
	bool CalculateTransform3(double T1[], double T2[], double M[16]);

    void GetRotationAxisByPlanes(int num_point, double  *p[8], double init_par[6], double cal_par[6], double *stdDev, double *maxDev);
	

	// Misc.
	bool Normalize(double &nx,double &ny, double &nz);
	double ran0(int *idum);
	double gaussdev(int *idum);
	void ConvertAngle2Matrix(double alphe, double beta, double gamma, double tx, double ty, double tz, double t[4][4]);		
	
		
	void CreateCirclePoints(double xc, double yc, double zc, double R, double nx, double ny, double nz, int N, int rpy_flag, double cir[][6]);
    void GetCylinder(int num_point, double  *p[3], double int_ori[3], double *nx, double *ny, double *nz, double *x0, double *y0, double *z0, double *R, double *stdDev, double *maxDev);


	// Rotate the input points Point3D[i][3]
	// into the plane that is in parellel with xy plane (Z value are the same)
	// input: points in 3D space (they may lay on a plane): Point3D[]
	// output: points after transform: *Point3D[] in which all Z are the same
	//         its rotation angles are *rpy[3]
	void Rotate3DinZ(int numPoint, double *Point3D[], double rpy[3]);

	// Rotate 3D points Point3D[i][3]
	// with angles rpy[3]
	// input/out: Point3D[i][3]
	void Rotate3D(int numPoint, double *Point3D[], double rpy[3]);
	
	//Given a circle parameters (xc,yc,zc,R,nx,ny,nz)
	// create Arc points arc[N][6]
	// arc portion is defined by three sample point samplePt[9][3];
	void CreateArcPoints(double xc, double yc, double zc, double R, double nx, double ny, double nz, int N, double samplePt[9], int rpy_flag, double arc[][6]);
		

	void RegisterationPoints2Points(double **PtsFrom, double **PtsTo, long numPts, double **matrix, double *ms);
	

	// The following functions are used to calculate the distance between
	// a point to a triangle

	bool IsLeagleTriangle(double* pdTriangle, double dMinDistance);

	/*Function: Calculate square distance of a point to a triangle by
	input three point coordinate.
	return : double as square distance.*/
	double gCalcSqrDistPtToTriangle(double* pdPt, double* pdTriangle,
									double *pdS, double* pdT);

	/*Function: Calculate square distance of a point to a triangle by
	input base point coordinate, edge0 and edge1 coordinate.
	return : double as square distance.*/
	double gCalcSqrDistPtToTriangle1(double* pdPt, double* pdTriangle,
									 double *pdS, double* pdT);

	/*Function: Calculate distance of a point to a triangle  by
	input three point coordinate.*/
	double gCalcDistPtToTriangle(double* pdPt, double* pdTriangle);

	/*Function: Calculate distance of a point to a triangle and output
			   the matched point's coordinate in the triangle by input 
			   three triangle's point coordinate.*/
	double gCalcDistPtToTriangleWithMp(double* pdPt, double* pdTriangle,
									   double *pdMp);

	/*Function: Calculate max-edge length of all input triangle.*/
	double gMaxLengthOfAllTriangle(int iTriNum, double* pdTriangle,double *pdAveLen);

	/*Function: Calculate a point is located at which side of a triangle.*/
	double gPtWithTriangle(double *pdPt,double* pdTriangle);

	/*Function: Calculate the errors of all input point to the surface 
	which is composed of all input triangles.*/
	bool gCalcError(int iPtNum, int iTriNum, double *pdPtTab,double* pdTriTab,
					double dMaxErr,double dHitValue);

	void GetConic(int num_point, double  *p[3], double *x0, double *y0, double *z0, double *nx,double *ny,double *nz, double *stdDev, double *maxDev);
private:
	void Transfer3Dto2D(int numPoint, double* Point3D[3], double* Point2D[2]);

	double nx,ny,nz,d;
	double pitch,yaw,roll;
};

namespace stdMath {

    unsigned char GetR(unsigned int);
    unsigned char GetG(unsigned int);
    unsigned char GetB(unsigned int);
#ifdef WIN32
#else
    unsigned char GetRValue(unsigned int);
    unsigned char GetGValue(unsigned int);
    unsigned char GetBValue(unsigned int);
#endif
    unsigned int RGBA(unsigned char r,unsigned char g,unsigned char b);


}


using namespace stdMath;


#endif // !defined(AFX_MYMATH_H__4AC3D190_94EC_11D2_B21D_00C04FBEDB8F__INCLUDED_)
