 
// Pgmath.cpp: implementation of the CStdMath class.
//
//////////////////////////////////////////////////////////////////////
#include "stdMath.h"
#include <math.h>
#include <stdlib.h>
#include "Mrqmin.h"
#include "svdcmp.h"
#include "polyfit.h"
#include "mv.h"

using namespace mv_math;

#ifndef PI
#define PI		3.14159265358979
#endif 

#ifdef WIN32
#pragma warning(disable : 4101) // Identifier was truncated to '255' characters in the debug information

#endif

#ifndef DELTA
#define	DELTA 0.0005
#endif 

#define SPHERE_FIT 6
#define CIRCLE_FIT 13
  


//for double multiplication
void CStdMath::MultiplyMatrixFloat(double m1[], int row1, int col1, double m2[], int col2, double m3[])
{
	for(int i=0; i<row1; i++)
		for( int j=0; j<col2; j++)
		{
			m3[(col2*i)+j]=0;
			for(int t=0; t<col1; t++)
				m3[(col2*i)+j] += m1[(col1*i)+t]*m2[(col2*t)+j];
		}
		return;
}


/*************************************************	

Function Name: MultipleMatrix

Functionality:	multiply matrix

Input: double[][] m1,	-- matrix 1
		int row1		-- row no. of matrix 1
		int col1		-- col no. of matrix 1
		double[][] m2	-- matrix 2
		int col2		-- col no. of matrix 2
		double m3		-- output matrix
		Note: col 1= row 2
Return: XForm

***************************************************/

void CStdMath::MultiplyMatrix(double m1[], int row1, int col1, double m2[], int col2, double m3[])
{
	for(int i=0; i<row1; i++)
		for( int j=0; j<col2; j++)
		{
			m3[(col2*i)+j]=0;
			for(int t=0; t<col1; t++)
				m3[(col2*i)+j] += m1[(col1*i)+t]*m2[(col2*t)+j];
		}
	return;
}



   
void CStdMath::LinearLeastSquare(int numRow, int numCol, double **A, double *b, double *x, double *stdDev, double *maxError)
{  
 
	int i, j, k;

 	double**	tmpA;
	double**	tmpB;

  	tmpA = new double*[numCol];
 	tmpB = new double*[numCol];
	for(int i = 0; i < numCol; i++) {
		tmpA[i] = new double[numCol];
		tmpB[i] = new double[numRow];
	}
 
	// linear square		
  	int		rl, rh, cl, ch;
 	double** u;
	double*	d;
	double** v;
	
 	rl = 1;
	rh = numRow;
	cl = 1;
	ch = numCol;

	u = matrix(rl, rh, cl, ch);
	d = vector(cl, ch);
	v = matrix(rl, rh, cl, ch);

	for(k = rl; k <= rh; k++) {
		for(j = cl; j <= ch; j++) {
			u[k][j] = (double)A[k-1][j-1];    
		}
	}
	
	// decomposition the Jacobian matrix
	// a = u*s*v
	// a is replaced by u
	svdcmp(u, rh, ch, d, v);

//	tmp.Format(_T("svd --d[]:\n d1:%.4lf d2:%.4lf ...d%d:%.4lf d%d:%.4lf"), 
//				d[1], d[2], numCol-1, d[numCol-1], numCol, d[numCol]);
//	AfxMessageBox(tmp);


	// [A].[x] = [b]
	// [A] = [U].[diag(w)][Vt]
	// [x] = [V].[diag(1/w)].[Ut].[b]
 	for(k = 0; k < numCol; k++) {
		for(int i = 0; i < numCol; i++) {
			tmpA[k][i] =(double) ( v[k+1][i+1] * (1.0/d[i+1]));
		}
 	}

		
	for(k = 0; k < numRow; k++) {
		for(int i = 0; i < numCol; i++) {
				tmpB[i][k] = 0.0;
		}
	} 
  
	
	for(k = 0; k < numRow; k++) {
		for(int i = 0; i < numCol; i++) {
			for (int m=0; m<numCol; m++) {
				tmpB[i][k] += tmpA[i][m]*(double)u[k+1][m+1];
            }
		}

	}

	for(int i = 0; i < numCol; i++) {
		x[i] = 0.0;
	} 

	for(k = 0; k < numRow; k++) {
		for(int i = 0; i < numCol; i++) {
			x[i] += tmpB[i][k] * b[k];
		}
	}


    // calculate the fitting erro 
	// || AX - B||
    // 
	double * error = new double[numRow];
	double maxErr = -1000, sum = 0;
	for(int i = 0; i < numRow; i++) {
		double temp = 0;
		for(j = 0; j < numCol; j++) {
			temp += A[i][j]*x[j];
		}
		error[i] = fabs(temp - b[i]);


		if (maxErr < error[i]) 
			maxErr = error[i];
		sum += error[i]*error[i];

	}

	if (stdDev)
	    *stdDev = sqrt(sum/numRow);
	if (maxError)
	*maxError = maxErr;

	delete error;

 

 
//	tmp.Format(_T("SVD X[]:\n x[0]:%.4lf x[1]:%.4lf.. x[%d]:%.4lf x[%d]:%.4lf"),
//		     x[0], x[1], numCol-1, x[numCol-1],numCol,x[numCol]);
//	AfxMessageBox(tmp);

	// clean up
	free_matrix(u, rl, rh, cl, ch);
	free_matrix(v, rl, rh, cl, ch);
	free_vector(d, cl, ch);
	for(int i = 0; i < numCol; i++) {
		delete []tmpA[i];
		delete []tmpB[i];
	}
	delete []tmpA;
	delete []tmpB;




}


void CStdMath::NonLinearLeastSquare(double** x, double y[], double sig[], int ndata, 
			double a[], int ma, int lista[], int mfit, double** covar, 
			double** alpha, double* chisq, double* alamda, int index)
{
	//CString tmp;
	//CMainFrame* pMainWnd = (CMainFrame*)AfxGetMainWnd();

	// on the first call provide an initial guess for the parameters a, and set alamda < 0
	mrqmin(x, y, sig, ndata, a, ma, lista, mfit, covar, alpha, chisq, alamda, index);

	int count = 0;	// use count to control max loop
	int maxLoop = 100;//10000

	while(fabs(*chisq) > 0.01 && count < maxLoop) {
//	while(fabs(*chisq) > 0.000001 && count < maxLoop) {
		mrqmin(x, y, sig, ndata, a, ma, lista, mfit, covar, alpha, chisq, alamda, index);

		count++;
	}

	//tmmpp.Format(_T("chisq: %lf"), *chisq);
	//pMainWnd->m_wndOutput.m_BuildList.AddString(tmp);

	// make one final call with alamda = 0, 
	// so that covar returns the covariance matrix, and alpha the curvature matrix
	*alamda = 0.0f;
 	mrqmin(x, y, sig, ndata, a, ma, lista, mfit, covar, alpha, chisq, alamda, index);
 

}



void CStdMath::NonLinearLeastSquare(double** x, double y[], double sig[], int ndata, 
			double a[], int ma, int lista[], int mfit, int index)
{

 	double** covar;				// covar[1..ma][1..ma]
	double** alpha;				// alpha[1..ma][1..ma]
	double* chisq;	
	double* alamda;
	int maxLoop;
	double  chiTol;

  	covar = new double*[ma+1];
	alpha = new double*[ma+1];

 	for(int i = 1; i <= ma; i++) {
		covar[i] = new double[ma+1];                                                   
		alpha[i] = new double[ma+1];                                                   
	}

	chisq = new double();
	alamda = new double();

	// initialization
	*alamda = -100.0f;			// must be initialized to less than zero
	*chisq = 0.0f;

	// on the first call provide an initial guess for the parameters a, and set alamda < 0
	mrqmin(x, y, sig, ndata, a, ma, lista, mfit, covar, alpha, chisq, alamda, index);

	int count = 0;	// use count to control max loop
	if (index == 3){ // for TCP calibration -- 
	  maxLoop = 20000;
	  chiTol = 0.0000000001;
	} else if (index ==2){  // for cylinder fitting --  06/01/04
	  maxLoop = 40;
	  chiTol = 0.1;
	} else if (index ==6){  // for sphere fitting --  06/01/04
	  maxLoop = 20;
	  chiTol = 0.1;
	}else {
      maxLoop = 10000;
	  chiTol = 0.00000001;
	}


    // orginal chisq < 0.000001
	while(fabs(*chisq) > chiTol && count < maxLoop) {
		mrqmin(x, y, sig, ndata, a, ma, lista, mfit, covar, alpha, chisq, alamda, index);

		count++;
	}


	// make one final call with alamda = 0, 
	// so that covar returns the covariance matrix, and alpha the curvature matrix
	*alamda = 0.0f;
 	mrqmin(x, y, sig, ndata, a, ma, lista, mfit, covar, alpha, chisq, alamda, index);


	// clean up memory  -- QT 6/5/01
 	for(int i = 1; i <= ma; i++) {
		delete []covar[i];                                                   
		delete []alpha[i];                                                   
	}

	delete []covar;
	delete []alpha;

	delete chisq;
	delete alamda;

}

void CStdMath::NonLinearLeastSquare_calib(double** x, double y[], double sig[], int ndata, 
			double a[], int ma, int lista[], int mfit, int index, double *pXYZw)
{

 	double** covar;				// covar[1..ma][1..ma]
	double** theta;				// theta[1..ma][1..ma]
	double* chisq;	
	double* alamda;
	int maxLoop;
	double  chiTol;

  	covar = new double*[ma+1];
	theta = new double*[ma+1];

 	for(int i = 1; i <= ma; i++) {
		covar[i] = new double[ma+1];                                                   
		theta[i] = new double[ma+1];                                                   
	}

	chisq = new double();
	alamda = new double();

	// initialization
	*alamda = -100.0f;			// must be initialized to less than zero
	*chisq = 0.0f;

	// on the first call provide an initial guess for the parameters a, and set alamda < 0
	mrqmin2(x, y, sig, ndata, a, ma, lista, mfit, covar, theta, chisq, alamda, index, pXYZw);

	int count = 0;	// use count to control max loop
	if (index == 3){ // for TCP calibration -- 
	  maxLoop = 20000;
	  chiTol = 0.0000000001;
	} else {
      maxLoop = 10000;
	  chiTol = 0.00000001;
	}


    // orginal chisq < 0.000001
	while(fabs(*chisq) > chiTol && count < maxLoop) {
		mrqmin2(x, y, sig, ndata, a, ma, lista, mfit, covar, theta, chisq, alamda, index, pXYZw);

		count++;
	}


	// make one final call with alamda = 0, 
	// so that covar returns the covariance matrix, and alpha the curvature matrix
	*alamda = 0.0f;
 	mrqmin2(x, y, sig, ndata, a, ma, lista, mfit, covar, theta, chisq, alamda, index, pXYZw);


	// clean up memory  -- QT 6/5/01
 	for(int i = 1; i <= ma; i++) {
		delete []covar[i];                                                   
		delete []theta[i];                                                   
	}

	delete []covar;
	delete []theta;

	delete chisq;
	delete alamda;

}




void CStdMath::NonLinearLeastSquare_calib_I(double** x, double y[], double sig[], int ndata, 
			double a[], int ma, int lista[], int mfit, int index, double *pPlanVector)
{

 	double** covar;				// covar[1..ma][1..ma]
	double** theta;				// theta[1..ma][1..ma]
	double* chisq;	
	double* alamda;
	int maxLoop;
	double  chiTol;

  	covar = new double*[ma+1];
	theta = new double*[ma+1];

 	for(int i = 1; i <= ma; i++) {
		covar[i] = new double[ma+1];                                                   
		theta[i] = new double[ma+1];                                                   
	}

	chisq = new double();
	alamda = new double();

	// initialization
	*alamda = -100.0f;			// must be initialized to less than zero
	*chisq = 0.0f;

	// on the first call provide an initial guess for the parameters a, and set alamda < 0
	mrqmin3(x, y, sig, ndata, a, ma, lista, mfit, covar, theta, chisq, alamda, index, pPlanVector);

	int count = 0;	// use count to control max loop
	if (index == 3){ // for TCP calibration -- 
	  maxLoop = 20000;
	  chiTol = 0.0000000001;
	} else {
      maxLoop = 10000;
	  chiTol = 0.00000001;
	}


    // orginal chisq < 0.000001
	while(fabs(*chisq) > chiTol && count < maxLoop) {
		mrqmin3(x, y, sig, ndata, a, ma, lista, mfit, covar, theta, chisq, alamda, index, pPlanVector);

		count++;
	}


	// make one final call with alamda = 0, 
	// so that covar returns the covariance matrix, and alpha the curvature matrix
	*alamda = 0.0f;
 	mrqmin3(x, y, sig, ndata, a, ma, lista, mfit, covar, theta, chisq, alamda, index, pPlanVector);


	// clean up memory  -- QT 6/5/01
 	for(int i = 1; i <= ma; i++) {
		delete []covar[i];                                                   
		delete []theta[i];                                                   
	}

	delete []covar;
	delete []theta;

	delete chisq;
	delete alamda;

}




/* *****************************************************************************
 * MatrixInverse is used to calculate the inverse of a square matrix
 * Input: 
 *      A -- pointer to a  square matrix 
 * Output: 
 *      A  -- pointer to a square matrix after inverse 
 * Return Void
 *
 * Note:
 * Calculation is based on SDV algorithm
 *  A = [U].[diag(w)][Vt] by SDV algorithm
 * inv(A)  = [V].[diag(1/w)].[Ut]
 *      
 * The previous function InverseMatrix() is actually a transpose of a orthogonal matrix.
 *
 *  Author: 
 *         11.11.00
 * **********************************************************************************/

void CStdMath::MatrixInverse(int numCol, double **A)
{  
 
	int i, j, k;

 	double**	tmpA;

  	tmpA = new double*[numCol];
	for(int i = 0; i < numCol; i++) {
		tmpA[i] = new double[numCol];
	}
 
	// linear square		
  	int		rl, rh, cl, ch;
 	double** u;
	double*	d;
	double** v;
	
 	rl = 1;
	rh = numCol;
	cl = 1;
	ch = numCol;

	u = matrix(rl, rh, cl, ch);
	d = vector(cl, ch);
	v = matrix(rl, rh, cl, ch);

	for(k = rl; k <= rh; k++) {
		for(j = cl; j <= ch; j++) {
			u[k][j] = (double)A[k-1][j-1];    
		}
	}
	
	// decomposition the Jacobian matrix
	svdcmp(u, rh, ch, d, v);

	// A = [U].[diag(w)][Vt]
	// inv(A)  = [V].[diag(1/w)].[Ut]
 	for(k = 0; k < numCol; k++) {
		for(int i = 0; i < numCol; i++) {
			tmpA[k][i] = v[k+1][i+1] * (1/d[i+1]);
        }
	}

  
    for(k = 0; k < numCol; k++) {
		for(int i = 0; i < numCol; i++) {
			A[i][k] = 0; // initialize A
			for (int m=0; m<numCol; m++) {
				A[i][k] += tmpA[i][m]*(double)u[k+1][m+1];
            }
		}

	}



	for(int i = 0; i < numCol; i++) {
		delete []tmpA[i];
	}
	delete []tmpA;

	// clean up
	free_matrix(u, rl, rh, cl, ch);
	free_matrix(v, rl, rh, cl, ch);
	free_vector(d, cl, ch);

}


// m3[3][3] = m1[3][3]*m2[3][3]
void CStdMath::MatrixMultiply33(double m1[3][3],double m2[3][3], double m3[3][3])
{
   
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
             m3[i][j] = 0.0;

			 for (int k=0; k<3; k++){
			    m3[i][j] = m3[i][j] + m1[i][k]*m2[k][j];
             }

		}
	}
		
		
}

// m3[4][4] = m1[4][4]*m2[4][4]
void CStdMath::MatrixMultiply44(double m1[4][4],double m2[4][4], double m3[4][4])
{
   
	for (int i=0; i<4; i++) {
		for (int j=0; j<4; j++) {
             m3[i][j] = 0.0;

			 for (int k=0; k<4; k++){
			    m3[i][j] = m3[i][j] + m1[i][k]*m2[k][j];
             }

		}
	}
		
		
}


// v2[4] = m1[4][4]*v1[4]
void CStdMath::MatrixMultiplyVector44(double m1[4][4],double v1[4], double v2[4])
{
   
	for (int i=0; i<4; i++) {
             v2[i] = 0.0;
			 for (int j=0; j<4; j++){
			    v2[i] = v2[i] + m1[i][j]*v1[j];
             }
	}
}


// v2[3] = m1[3][3]*v1[3]
void CStdMath::MatrixMultiplyVector33(double m1[3][3],double v1[3], double v2[3])
{
   
	for (int i=0; i<3; i++) {
             v2[i] = 0.0;
			 for (int j=0; j<3; j++){
			    v2[i] = v2[i] + m1[i][j]*v1[j];
             }
	}
}




/*	Purpose : Calculate mean and variance of date set
	Input:	a  ---  pointer to data array
			n  ---	Dimension of Data Array
	Output: med --- mean of Data Array
			sgama---variance of Data Array*/
template<class T> void meanvar(T *a, int n, double &med, double &sgama)
{
	double sum = 0, sumsq = 0;
	if (n < 2)
	{
		med = a[0];
		return;
	}

	for (int i = 1; i < n; i++)
		sum += a[i];
	med = sum / (n - 1);
	for (int i = 1; i < n; i++)
		sumsq += (a[i] - med) * (a[i] - med);
	sgama = sqrt(sumsq / (n-1));
}





/* *****************************************************************************
 * Atan360() is used to calculate arctan() value with the range 0 - 360
 * It is the modification for standard function atan() that results in the value from -PI/2 to PI/2
 *
 * Input: 
 *      numerator (num)
 *      denominator (den) to form triangle
 * Output: 
 *      arc tangent value with value from 0 to 360 degree
 * Return 
 *       Output
 *
 * Author: 
 *         12.04.00
 * **********************************************************************************/

double CStdMath::Atan360(double num, double den)
{
double A;

// double A0
//A0 = atan(num/den);
//if (num>=0 && den >0)           // first Quarter
//   A = A0;
//else if (num>=0 && den < 0)     // second quarter
//   A = PI + A0;                 //  A0 is negative ang
//else if (num<=0 && den < 0)     // Third quarter
//   A = PI + A0;                 // A0 is positive
//else                            // num<0 & den >0 fourth quarter
//   A = 2*PI + A0;               // A0 is negative*/

A = atan2(num, den);
if (A < 0)
	A = 2 * PI + A;

 return A;
}



/* *****************************************************************************/
// Rotate a point xyz_in[3] along its Y axis
// input: a 3d point xyz_in[3]
//        rotataion angle [angle]
// output: xyz_out[3];
// , 9/19/01
/*****************************************************************************/
void CStdMath::RotateYAxis(double xyz_in[3], double angle, double xyz_out[3])
{
    double rot[3][3];

// Rotate along Y axis the rotation matrix is 
    //rot = [cos(angle) 0 -sin(angle)
  // 		0          1       0
  //		sin(angle) 0       cos(angle)];
  
 
// Cam_Pos = rot*Base_Cam_Pos;
   angle = -angle*PI/180;
 
   rot[0][0] = cos(angle);
   rot[0][1] = 0;
   rot[0][2] = -sin(angle);
   rot[1][0] = 0;
   rot[1][1] = 1;
   rot[1][2] = 0;
   rot[2][0] = sin(angle);
   rot[2][1] = 0;
   rot[2][2] = cos(angle);

   xyz_out[0] = cos(angle)*xyz_in[0] - sin(angle)*xyz_in[2];
   xyz_out[1] = xyz_in[1];
   xyz_out[2] = sin(angle)*xyz_in[0] + cos(angle)*xyz_in[2];
 
}




/*******************************************************************************************
 * CircleFitting() calculates the circle parameters based on input points on the circle 
 * Inputs: 
 *        cir[][3]; containing a array of [x,y,z];
 *        num_point; number of pints of the input; 
 * 
 * Outputs:  
 *           poimnters of Circle parameters including 
 *           center position (xc, yc,zc), 
 *           radius R
 *           orientation (Nx, Ny, Nz) in space
 *          
 * return:
 *         void
 * 
 * Notes for approach
 * 1. first step is to calculate the normal vector for the plane that input points lie on.
 * 2. Rotate the input points (or its plane) such that they lie on the plane that is in parellel with xy plane. 
 *    that is,  z is constant.
 * 3. Do two dimensional circle fitting for rotated points. (considering x & y only) to find center and R of the circle
 * 4. Rotate the circle back and determine the real center (xc,yc,zc); 
 * 

 * Author: 
 *          6.5.00 
 */
void  CStdMath::CircleFitting(int num_point, double *cir[3], double *xc, double *yc, double *zc, double *R, double *Nx, double *Ny, double *Nz)
{
  
   	double c_xc, c_yc, c_radius;
    double **cir_2d;

	cir_2d = new double*[num_point];

	for (int k=0; k<num_point; k++)
	{
		 cir_2d[k] = new double[2];
    }

	Transfer3Dto2D(num_point, cir, cir_2d);

	// do circle fitting in two dimensional to find center and radius
	 CircleFitting_2D (num_point,cir_2d,&c_xc, &c_yc, &c_radius);
       	 
	 cir_2d[0][0] = c_xc;
	 cir_2d[0][1] = c_yc;

	 Transfer2Dto3D(1, cir_2d, cir);
	 
	  *xc = cir[0][0];
	  *yc = cir[0][1];
	  *zc = cir[0][2];

	 
	  *R = c_radius;


       for (int k=0; k<num_point; k++)
	   {
		delete  []cir_2d[k];
	   }

       delete []cir_2d;
}



// The same as CircleFitting() except for
// returning stdDev and maxDev

void  CStdMath::CircleFitting(int num_point, double *cir[3], double *xc, double *yc, double *zc, double *R, double *Nx, double *Ny, double *Nz, double *stdDev, double *maxDev)
{
  
    double temp_x,temp_y,temp_z,x1,y1,z1,x2,y2,z2,x,y,z, pl_stdDev, pl_maxDev;
 	double nx, ny, nz, dis, c_xc, c_yc, c_radius;
    double real_nx, real_ny, real_nz, roll, pitch, yaw;
    double **cir_2d;
	double p[3], d[3], n[3];
    PT_RPY	Coord;

	cir_2d = new double*[num_point];

	for (int k=0; k<num_point; k++){
		 cir_2d[k] = new double[2];

    }


    
	GetPlaneNormal(num_point, cir, &nx, &ny, &nz, &dis, &pl_stdDev, &pl_maxDev);
    // normalize the orientation vector in case they are not normalized
	real_nx = nx/sqrt(nx*nx + ny*ny+ nz*nz);
	real_ny = ny/sqrt(nx*nx + ny*ny+ nz*nz);
	real_nz = nz/sqrt(nx*nx + ny*ny+ nz*nz);



	/* 
	 * transform the input points into the plane that is in parellel with xy plane
	 */
	// assume that this vector pass origin
	// since we only consider the orientation
	p[0] = 0.0;
	p[1] = 0.0;
	p[2] = 0.0;

	n[0] = real_nx;
	n[1] = real_ny;
	n[2] = real_nz;

	/* 
	 * Need determine (give) direction vector to create RPY or coordinate
	 * since ijy info is not enogh to determine coordinate
	 */
	 
	// Assume direction vector d[] is normal to n[] or (i,j,k)
	if (real_nx != 0.0){
	  d[0] = -real_ny*0.5/real_nx;
	  d[1] = 0.5;
	  d[2] = 0.0;
	} else if (real_ny != 0.0){
        d[0] = 0.5;
		d[1] = -real_nx*0.5/real_ny;
	    d[2] = 0.0;
    } else { // direction is x 
		d[0] = 1.0;
	    d[1] = 0.0;
		d[2] = 0.0;
	}


	// Define a coordinate for circle based on above normal and direction vector
	Coord = PTIJK_to_PTRPY(p, d, n);


    /*
	 * roll, pitch & yaw indicate how the designed circle coordinate transfers 
	 *into world coordinate. 
	 * From world to circle coordinate reverse the direction and angle
	 */

	roll = Coord.roll*PI/180.0;
	pitch = Coord.pitch*PI/180;
	yaw = Coord.yaw*PI/180.0;
	
	for (int i=0; i< num_point; i++){
	 
	  /* 
	   * Notes:
	   * Rotate it to align with [nx, ny, nz]
	   * rotate along x axis 
	   * y' = cos()*y + sin()*z;
	   * z' = -sin()*y + cos()*z;
	   * rotate along y axis
	   * x' = cos()*x - sin()*z;
	   * z' = sin()*x + cos()*z;
	   * rotate along z axis
	   * x' = cos()*x + sin()*y;
	   * y' = -sin()*x + cos()*y;
	   */
	   
	  /*
	   * Note: 
	   * from circle to  world (forward transform)
	   * a. rotate along z with roll degree
	   * b. rotate along y with pitxh degree
	   * c. rotate along x with yaw degree
	   * From world to circle (inverse trnasform)
	   * a. rotate along x with -yaw degree
	   * b. rotate along y with -pitch degree
	   * c. rotate along z with -roll degree
	   * Here roll pitch, and yaw have been revesed already
	   */


      temp_x = cir[i][0];
	  temp_y = cir[i][1];
	  temp_z = cir[i][2];


	  // first rotate along z axis with roll degree
	  x1 = cos(roll)*temp_x + sin(roll)*temp_y;
	  y1 = -sin(roll)*temp_x + cos(roll)*temp_y;
      z1 = temp_z;
      
	  // rotate along y axis with pitch degree
	  x2 = cos(pitch)*x1 - sin(pitch)*z1;
	  z2 = sin(pitch)*x1 + cos(pitch)*z1;
	  y2 = y1;

      
	  // last rotate along x axis w yaw degree
	  y = cos(yaw)*y2+ sin(yaw)*z2;
	  z = -sin(yaw)*y2 + cos(yaw)*z2;
	  x = x2;
	

      // rotated circle points will be
	  // z should be close to constatnt
	  cir[i][0] = x;
	  cir[i][1] = y;
	  cir[i][2] = z;
	
	  cir_2d[i][0] = cir[i][0];
	  cir_2d[i][1] = cir[i][1];
   
	}

	// do circle fitting in two dimensional to find center and radius
	 CircleFitting_2D(num_point,cir_2d,&c_xc, &c_yc, &c_radius);
       	 

	 /************** 
	  * rotate the circle center back to reflect the resiult in the world coordinate
	  ***************/

      temp_x = c_xc;
	  temp_y = c_yc;
	  temp_z = cir[0][2];   // all z should be very close   -- better to take average 

	  
	  // revese the rotation direction 
	   roll = -roll;
	   pitch = -pitch;
	   yaw = -yaw;
     

	  // first rotate along x axis (yaw)
	  y1 = cos(yaw)*temp_y + sin(yaw)*temp_z;
	  z1 = -sin(yaw)*temp_y + cos(yaw)*temp_z;
	  x1 = temp_x;


	  // rotate along y axis  (pitch)
	  x2 = cos(pitch)*x1- sin(pitch)*z1;
	  z2 = sin(pitch)*x1 + cos(pitch)*z1;
	  y2 = y1;

	  // rotate along z axis (roll)
      x = cos(roll)*x2 + sin(roll)*y2;
	  y = -sin(roll)*x2 + cos(roll)*y2;
      z = z2;

	  	  	  
      // final result after rotaion back
	  *xc = x;
	  *yc = y;
	  *zc = z;
	  *Nx = nx;
	  *Ny = ny;
	  *Nz = nz;
	  *R = c_radius;


	  // calculate starndard deveation and max devation
	  
		  double maxError, sum, delta;
		  maxError = -100000;
		  sum = 0;

		  for (int	i=0; i< num_point; i++){
			   dis = sqrt((cir_2d[i][0] - c_xc)*(cir_2d[i][0]-c_xc) + (cir_2d[i][1] - c_yc)*(cir_2d[i][1]-c_yc));
			   delta = fabs(dis-c_radius);
			   if (maxError < delta) maxError = delta;
			   sum = sum + delta;
		  }
		  *stdDev = sum/num_point;
		  *maxDev = maxError;
	
   	  // cleanup memory
   	for (int k=0; k<num_point; k++){
		delete  []cir_2d[k];
    }
	delete []cir_2d;
 

}


// Circle fitting with nonlinear LS approch
void  CStdMath::CircleFitting_3D_NL(int num_point, double *cir[3], double *xc, double *yc, double *zc, double *R, double *Nx, double *Ny, double *Nz, double *stdDev, double *maxDev, bool bFixedRadius)
{
    double temp_x,temp_y,temp_z,x1,y1,z1,x2,y2,z2,x,y,z, pl_stdDev, pl_maxDev;
 	double nx, ny, nz, dis, c_xc, c_yc, c_radius;
    double real_nx, real_ny, real_nz, roll, pitch, yaw;
    double **cir_2d;
	double p[3], d[3], n[3];
    PT_RPY	Coord;

	cir_2d = new double*[num_point];

	for (int k=0; k<num_point; k++){
		 cir_2d[k] = new double[2];

    }


    
	GetPlaneNormal(num_point, cir, &nx, &ny, &nz, &dis, &pl_stdDev, &pl_maxDev);
    // normalize the orientation vector in case they are not normalized
	real_nx = nx/sqrt(nx*nx + ny*ny+ nz*nz);
	real_ny = ny/sqrt(nx*nx + ny*ny+ nz*nz);
	real_nz = nz/sqrt(nx*nx + ny*ny+ nz*nz);



	/* 
	 * transform the input points into the plane that is in parellel with xy plane
	 */
	// assume that this vector pass origin
	// since we only consider the orientation
	p[0] = 0.0;
	p[1] = 0.0;
	p[2] = 0.0;

	n[0] = real_nx;
	n[1] = real_ny;
	n[2] = real_nz;

	/* 
	 * Need determine (give) direction vector to create RPY or coordinate
	 * since ijy info is not enogh to determine coordinate
	 */
	 
	// Assume direction vector d[] is normal to n[] or (i,j,k)
	if (real_nx != 0.0){
	  d[0] = -real_ny*0.5/real_nx;
	  d[1] = 0.5;
	  d[2] = 0.0;
	} else if (real_ny != 0.0){
        d[0] = 0.5;
		d[1] = -real_nx*0.5/real_ny;
	    d[2] = 0.0;
    } else { // direction is x 
		d[0] = 1.0;
	    d[1] = 0.0;
		d[2] = 0.0;
	}


	// Define a coordinate for circle based on above normal and direction vector
	Coord = PTIJK_to_PTRPY(p, d, n);


    /*
	 * roll, pitch & yaw indicate how the designed circle coordinate transfers 
	 *into world coordinate. 
	 * From world to circle coordinate reverse the direction and angle
	 */

	roll = Coord.roll*PI/180.0;
	pitch = Coord.pitch*PI/180;
	yaw = Coord.yaw*PI/180.0;
	
	for (int i=0; i< num_point; i++){

      temp_x = cir[i][0];
	  temp_y = cir[i][1];
	  temp_z = cir[i][2];


	  // first rotate along z axis with roll degree
	  x1 = cos(roll)*temp_x + sin(roll)*temp_y;
	  y1 = -sin(roll)*temp_x + cos(roll)*temp_y;
      z1 = temp_z;
      
	  // rotate along y axis with pitch degree
	  x2 = cos(pitch)*x1 - sin(pitch)*z1;
	  z2 = sin(pitch)*x1 + cos(pitch)*z1;
	  y2 = y1;

      
	  // last rotate along x axis w yaw degree
	  y = cos(yaw)*y2+ sin(yaw)*z2;
	  z = -sin(yaw)*y2 + cos(yaw)*z2;
	  x = x2;
	

      // rotated circle points will be
	  // z should be close to constatnt
	  cir[i][0] = x;
	  cir[i][1] = y;
	  cir[i][2] = z;
	
	  cir_2d[i][0] = cir[i][0];
	  cir_2d[i][1] = cir[i][1];
   
	}

	// do circle fitting in two dimensional to find center and radius
	 CircleFitting_2D(num_point,cir_2d,&c_xc, &c_yc, &c_radius);
    
   	 // do circle fitting wirh nonlinear method
	 double init_par[3];
	 init_par[0] = c_xc;
	 init_par[1] = c_yc;
	 init_par[2] = c_radius;
	 double stdDev1, maxDev1;
	 Get2DCircleFit_NL(num_point, cir_2d, init_par, &c_xc, &c_yc, &c_radius, &stdDev1, &maxDev1, bFixedRadius);
		 
	 /************** 
	  * rotate the circle center back to reflect the resiult in the world coordinate
	  ***************/

      temp_x = c_xc;
	  temp_y = c_yc;
	  temp_z = cir[0][2];   // all z should be very close   -- better to take average 

	  
	  // revese the rotation direction 
	   roll = -roll;
	   pitch = -pitch;
	   yaw = -yaw;
     

	  // first rotate along x axis (yaw)
	  y1 = cos(yaw)*temp_y + sin(yaw)*temp_z;
	  z1 = -sin(yaw)*temp_y + cos(yaw)*temp_z;
	  x1 = temp_x;


	  // rotate along y axis  (pitch)
	  x2 = cos(pitch)*x1- sin(pitch)*z1;
	  z2 = sin(pitch)*x1 + cos(pitch)*z1;
	  y2 = y1;

	  // rotate along z axis (roll)
      x = cos(roll)*x2 + sin(roll)*y2;
	  y = -sin(roll)*x2 + cos(roll)*y2;
      z = z2;

	  	  	  
      // final result after rotaion back
	  *xc = x;
	  *yc = y;
	  *zc = z;
	  *Nx = nx;
	  *Ny = ny;
	  *Nz = nz;
	  *R = c_radius;


	  // calculate starndard deveation and max devation
	  
		  double maxError, sum, delta;
		  maxError = -100000;
		  sum = 0;

		  for (int		i=0; i< num_point; i++){
			   dis = sqrt((cir_2d[i][0] - c_xc)*(cir_2d[i][0]-c_xc) + (cir_2d[i][1] - c_yc)*(cir_2d[i][1]-c_yc));
			   delta = fabs(dis-c_radius);
			   if (maxError < delta) maxError = delta;
			   sum = sum + delta;
		  }
		  *stdDev = sum/num_point;
		  *maxDev = maxError;
	
   	  // cleanup memory
   	for (int k=0; k<num_point; k++){
		delete  []cir_2d[k];
    }
	delete []cir_2d;
 

}



/*******************************************************************************************
 * GetPlaneNormal() calculates the palne normal vector based on input points 
 * by using least sqaure method
 * Inputs: 
 *        a[][3]; containing a array of [x,y,z];
 *        num_point; number of pints of the input; 
 * Outputs:  
 *           plane normal vector (nx, ny nz) and d (distance from origin) 
 *                     
 * return:
 *         void
 * 
 * Notes:
 * The plane equation is given by 
 * x1*x + x2*y + x3*z = 1
 * The parameter X = [x1 x2 x3]' can be found by using LS eqation
 * A*X = 1 where A is the input point a[][] containing (x,y,z) values;
 * Assuming  k = sqrt(x1*x1 + x2*x2 + x3*x3);
 * The we can get standard plane equation
 * (x1/k)*x + (x2/k)*y + (x3/k)*z = 1/k;
 * or nx*x + ny*y + nz*z = d;
 * where nx = x1/k, ny = x2/k, nz = x3/k, d = 1/k;
 *
 * Author: 
 *          6.9.00 
 */

void CStdMath::GetPlaneNormal(int num_point, double  *a[3], double *nx, double *ny, double *nz, double *d)
{
  // solve LS equation with form AX = B
/*
	double *B;
    double X[3];
	double k;
    int num_col = 3;

	B = new double[num_point];

	for (int i=0; i<num_point; i++){
		B[i] = 1;
	}

	LinearLeastSquare(num_point, num_col, a, B, X); 
    
	k = sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
	
	*nx = X[0]/k;
    *ny = X[1]/k;
	*nz = X[2]/k;
    *d = 1/k;


	delete []B;
*/

	// PlaneFitting can handle special case; 03/31/04 
    double stdDev, maxDev;

	PlaneFitting(num_point, a, nx, ny, nz, d, &stdDev, &maxDev);


}


/*******************************************************************************************
 * GetPlaneNormal() calculates the palne normal vector based on input points 
 * by using least sqaure method
 * Inputs: 
 *        a[][3]; containing a array of [x,y,z];
 *        num_point; number of pints of the input; 
 * Outputs:  
 *           plane normal vector (nx, ny nz) and d (dustance from origin) 
 *           standard deviation (stdDev) and max deviation (maxDev)  of the fitting 
 *                     
 * return:
 *         void
 * 
 * Notes:
 * The plane equation is given by 
 * x1*x + x2*y + x3*z = 1
 * The parameter X = [x1 x2 x3]' can be found by using LS eqation
 * A*X = 1 where A is the input point a[][] containing (x,y,z) values;
 * Assuming  k = sqrt(x1*x1 + x2*x2 + x3*x3);
 * The we can get standard plane equation
 * (x1/k)*x + (x2/k)*y + (x3/k)*z = 1/k;
 * or nx*x + ny*y + nz*z = d;
 * where nx = x1/k, ny = x2/k, nz = x3/k, d = 1/k;
 *
 * Author: 
 *          6.9.00 
 */

void CStdMath::GetPlaneNormal(int num_point, double  *a[3], double *nx, double *ny, double *nz, double *d, double *stdDev, double *maxDev)
{
  // solve LS equation with form AX = B
/*
	double *B;
    double X[3];
	double k;
    int num_col = 3;

	B = new double[num_point];

	for (int i=0; i<num_point; i++){
		B[i] = 1;
	}

	LinearLeastSquare(num_point, num_col, a, B, X); 
	    
	k = sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
	
    // find the fitting error
    // distance from measured point to plane
	double dis;
	double max_d = 0.0;
	double dev = 0.0;
	for (i=0; i<num_point; i++){
		  // distance from measured points to the fitted plane
		  dis = fabs((X[0]*(a[i][0]) + X[1]*(a[i][1]) + X[2]*(a[i][2]) -1.0)/k);
		  if (dis > max_d) max_d = dis;
		  dev = dev + dis * dis;
    }
    dev = sqrt(dev/double(num_point));
	*nx = X[0]/k;
    *ny = X[1]/k;
	*nz = X[2]/k;
    *d = 1/k;
    *stdDev = dev;
	*maxDev = max_d;

	delete []B; 

*/

	// PlaneFitting can handle special case;03/31/04 

	PlaneFitting(num_point, a, nx, ny, nz, d, stdDev, maxDev);
	
	

}

/*******************************************************************************************
 * GetConic() calculates the conic center line (x0, y0, z0,nx,ny,nz) its angle (theta) based on input points 
 * by using nonlinear least sqaure method
 * Inputs: 
 *        p[][3]; containing a measured point array of [x,y,z] on the surface of the conic
 *        num_point; number of pints of the input; 
 * Outputs:  
 *         (x0,y0,z0,nx,ny,nz):center line.
		   theta: angle of fitted conic.
 *         standard devaiation (stdDev) and maxmium deviation (maxDev) of fitting 
 *                    
 * return:
 *         void
 * 
 * Notes:
 * 
 * 
 *         NingSuijun 2011.12.14
 */
void CStdMath::GetConic(int num_point, double  *p[3], double *x0, double *y0, double *z0, double *nx,double *ny,double *nz, double *stdDev, double *maxDev)
{
  
	// declaration & initialization
	int ndata = num_point;		// number of input points x[i][k]; i = 1,... ndata
	int nd = 3;					// nd: number of given x[i][k],  k = 1,2 ,3;
	int ma = 6;					// number of parameters to be determined after fitting
	int mfit = 6;				// = ma

	int* lista;					// lista[1..ma]
	double ** x;					// x[1..ndata][1..nd]	
	double * y;					// y[1..ndata]
	double * sig;					// sig[1..ndata]
	double * a;					// a[1..ma]

//	CString tmp;
//	CMainFrame* pMainWnd = (CMainFrame*)AfxGetMainWnd();
	
	// memory allocation
	lista = new int[ma+1];
	 
	x = new double*[ndata+1];
	for(int i = 1; i <= ndata; i++) {
		x[i] = new double[nd+1];                                                   		 
	}

	y = new double[ndata+1];
	sig = new double[ndata+1];
	a = new double[ma+1];
  
	for(int k = 1; k <= ma; k++) {
		lista[k] = k;			// numbers the parameters a such that 
 						// the first mfit elements correspond to values actually being adjusted
	}
  
	for(int k = 1; k <= ndata; k++) {
		sig[k] = 1;				// standard deviations, this is not experiment data, so set to 1
	    y[k] = 0.0;             // y = 0
	}

	double dsum[3];
	dsum[0]=dsum[1]=dsum[2]=0;
    // assign x value
    for(int i = 1; i <= ndata; i++) {
		for (int j = 1; j <= nd; j++){
			x[i][j] = (double) p[i-1][j-1];
			dsum[j-1]+=x[i][j];
		}
    	
	}



   // initial value of parameters
	double init_par[4];
	a[1] = dsum[0]/ndata;
	a[2] = dsum[1]/ndata;
	a[3] = dsum[2]/ndata;
	a[4] = *nx;
	a[5] = *ny;
	a[6] = *nz;
	//a[7] = *theta;
 
	NonLinearLeastSquare(x, y, sig, ndata, a, ma, lista, mfit, 14);

    
	// calculate the fitting error 
	double max_delta = -10000.0;
    double dd, delta, delta_sq,delta_sum;
	
	delta_sum = 0.0;
	delta_sum = sqrt(a[4]*a[4]+a[5]*a[5]+a[6]*a[6]);
	a[4] /=delta_sum;
	a[5] /=delta_sum;
	a[6] /=delta_sum;
	
	delta_sum = 0.0;
   
    for (int i = 0; i< num_point; i++){
	    
            // fitting error to sphare center
			dd = sqrt((p[i][0]-a[1])*(p[i][0]-a[1])+(p[i][1]-a[2])*(p[i][1]-a[2]) + (p[i][2]-a[3])*(p[i][2]-a[3])); 
			delta = dd*sin(PI/4);
			PointLineDistance(p[i][0],p[i][1],p[i][2],a[1],a[2],a[3],a[4],a[5],a[6],&delta_sq);
			delta_sq = fabs(delta - delta_sq);
            delta_sum = delta_sum + delta_sq;
			delta=delta_sq;
			if (max_delta < delta) max_delta = delta;


	}
	  *stdDev = sqrt(delta_sum/(double)num_point);
      *maxDev = max_delta;
   

	  // normalize the vector
       *x0 = a[1];
	   *y0 = a[2];
	   *z0 = a[3];
	   *nx = a[4];
	   *ny = a[5];
	   *nz = a[6];
	   //*theta = a[7]; 

    delete lista;
	for(int i = 1; i <= ndata; i++) {
		delete x[i];                                                   		 
	}
    delete x;
    delete y;
    delete sig;
    delete a;

 
//	tmp.Format(_T("Fitted results for Sphere"));
//	pMainWnd->m_wndOutput.m_BuildList.AddString(tmp);

//		tmp.Format(_T("X0:%lf  Y0:%lf  Z0:%lf R:%1f, StdDev: %.3f maxDev:%.2f"),
//			*x0,*y0, *z0, *R, *stdDev, *maxDev);
//	pMainWnd->m_wndOutput.m_BuildList.AddString(tmp);
}

/*******************************************************************************************
 * GetSphere() calculates the sphere center (x0, y0, z0) its radius (R) based on input points 
 * by using nonlinear least sqaure method
 * Inputs: 
 *        p[][3]; containing a measured point array of [x,y,z] on the surface of the sphere
 *        num_point; number of pints of the input; 
 *        int_par[4]; intial paramters (init_x0, init_y0, init_z0, init_R)
 * Outputs:  
 *         (x0,y0,z0) , R center and radius of fitted sphere.
 *         standard devaiation (stdDev) and maxmium deviation (maxDev) of fitting 
 *                    
 * return:
 *         void
 * 
 * Notes:
 * Since data fitting is obtained by nonlinear least method that requires initial values
 * The initial paramters can be the average of measurements, for exampe.
 * 
 * 
 *          07.29.01 
 */

void CStdMath::GetSphere(int num_point, double  *p[3], double init_par[4], double *x0, double *y0, double *z0, double *R, double *stdDev, double *maxDev)
{
  
	// declaration & initialization
	int ndata = num_point;		// number of input points x[i][k]; i = 1,... ndata
	int nd = 3;					// nd: number of given x[i][k],  k = 1,2 ,3;
	int ma = 4;					// number of parameters to be determined after fitting
	int mfit = 4;				// = ma

	int* lista;					// lista[1..ma]
	double ** x;					// x[1..ndata][1..nd]	
	double * y;					// y[1..ndata]
	double * sig;					// sig[1..ndata]
	double * a;					// a[1..ma]

//	CString tmp;
//	CMainFrame* pMainWnd = (CMainFrame*)AfxGetMainWnd();
	
	// memory allocation
	lista = new int[ma+1];
	 
	x = new double*[ndata+1];
	for(int i = 1; i <= ndata; i++) {
		x[i] = new double[nd+1];                                                   		 
	}

	y = new double[ndata+1];
	sig = new double[ndata+1];
	a = new double[ma+1];
  
	for(int k = 1; k <= ma; k++) {
		lista[k] = k;			// numbers the parameters a such that 
 						// the first mfit elements correspond to values actually being adjusted
	}
  
	for(int k = 1; k <= ndata; k++) {
		sig[k] = 1;				// standard deviations, this is not experiment data, so set to 1
	    y[k] = 0.0;             // y = 0
	}

    // assign x value
    for(int i = 1; i <= ndata; i++) {
		for (int j = 1; j <= nd; j++){
			x[i][j] = (double) p[i-1][j-1];
		}
    	
	}



   // initial value of parameters
	a[1] = init_par[0];
	a[2] = init_par[1];
	a[3] = init_par[2];
	a[4] = init_par[3];
 
     // call Nonlinear function
// 	CMyMath* m_math = new CMyMath();

	NonLinearLeastSquare(x, y, sig, ndata, a, ma, lista, mfit, SPHERE_FIT);
    // HOW TO USE NLLS:
	// (x,y): input point pair y= f(ax); normally y = 0;
	// sig: weighting factors. normally sig = 1;
	// ndata: number of data point.
	// a: coefficients to be determined -- outcome
	// ma: number of a
	// mfit: number of a to be fitted
	// lista: array of a [1 ...ma];

    
	// calculate the fitting error 
	double max_delta = -10000.0;
    double dd, delta, delta_sq, delta_sum;
	
	delta_sum = 0.0;

   
    for (int i = 0; i< num_point; i++){
	    
            // fitting error to sphare center
			dd = sqrt((p[i][0]-a[1])*(p[i][0]-a[1])+(p[i][1]-a[2])*(p[i][1]-a[2]) + (p[i][2]-a[3])*(p[i][2]-a[3])); 
			delta = fabs(dd - a[4]);
			delta_sq = (dd - a[4])*(dd - a[4]);
            delta_sum = delta_sum + delta_sq;
			if (max_delta < delta) max_delta = delta;


	}
	  *stdDev = sqrt(delta_sum/(double)num_point);
      *maxDev = max_delta;
   

	  // normalize the vector
       *x0 = a[1];
	   *y0 = a[2];
	   *z0 = a[3];
	   *R = a[4];

    delete lista;
	for(int i = 1; i <= ndata; i++) {
		delete x[i];                                                   		 
	}
    delete x;
    delete y;
    delete sig;
    delete a;

 
//	tmp.Format(_T("Fitted results for Sphere"));
//	pMainWnd->m_wndOutput.m_BuildList.AddString(tmp);

//		tmp.Format(_T("X0:%lf  Y0:%lf  Z0:%lf R:%1f, StdDev: %.3f maxDev:%.2f"),
//			*x0,*y0, *z0, *R, *stdDev, *maxDev);
//	pMainWnd->m_wndOutput.m_BuildList.AddString(tmp);
}


// Same as GetSphere() except having fixed radius
void CStdMath::GetSphere_fixRadius(int num_point, double  *p[3], double init_par[4], double *x0, double *y0, double *z0, double *R, double *stdDev, double *maxDev)
{
  
	// declaration & initialization
	int ndata = num_point;		// number of input points x[i][k]; i = 1,... ndata
	int nd = 3;					// nd: number of given x[i][k],  k = 1,2 ,3;
	int ma = 4;					// number of parameters to be determined after fitting
	int mfit = 3;				// = ma

	int* lista;					// lista[1..ma]
	double ** x;					// x[1..ndata][1..nd]	
	double * y;					// y[1..ndata]
	double * sig;					// sig[1..ndata]
	double * a;					// a[1..ma]

//	CString tmp;
//	CMainFrame* pMainWnd = (CMainFrame*)AfxGetMainWnd();
	
	// memory allocation
	lista = new int[ma+1];
	 
	x = new double*[ndata+1];
	for(int i = 1; i <= ndata; i++) {
		x[i] = new double[nd+1];                                                   		 
	}

	y = new double[ndata+1];
	sig = new double[ndata+1];
	a = new double[ma+1];
  
	for(int k = 1; k <= ma; k++) {
		lista[k] = k;			// numbers the parameters a such that 
 						// the first mfit elements correspond to values actually being adjusted
	}
  
	for(int k = 1; k <= ndata; k++) {
		sig[k] = 1;				// standard deviations, this is not experiment data, so set to 1
	    y[k] = 0.0;             // y = 0
	}

    // assign x value
    for(int i = 1; i <= ndata; i++) {
		for (int j = 1; j <= nd; j++){
			x[i][j] = (double) p[i-1][j-1];
		}
    	
	}


   // initial value of parameters
	a[1] = init_par[0];
	a[2] = init_par[1];
	a[3] = init_par[2];
	a[4] = init_par[3];
 
     // call Nonlinear function
// 	CMyMath* m_math = new CMyMath();

	NonLinearLeastSquare(x, y, sig, ndata, a, ma, lista, mfit, SPHERE_FIT);
    // HOW TO USE NLLS:
	// (x,y): input point pair y= f(ax); normally y = 0;
	// sig: weighting factors. normally sig = 1;
	// ndata: number of data point.
	// a: coefficients to be determined -- outcome
	// ma: number of a
	// mfit: number of a to be fitted
	// lista: array of a [1 ...ma];

    
	// calculate the fitting error 
	double max_delta = -10000.0;
    double dd, delta, delta_sq, delta_sum;
	
	delta_sum = 0.0;

   
    for (int i = 0; i< num_point; i++){
	    
            // fitting error to sphare center
			dd = sqrt((p[i][0]-a[1])*(p[i][0]-a[1])+(p[i][1]-a[2])*(p[i][1]-a[2]) + (p[i][2]-a[3])*(p[i][2]-a[3])); 
			delta = fabs(dd - a[4]);
			delta_sq = (dd - a[4])*(dd - a[4]);
            delta_sum = delta_sum + delta_sq;
			if (max_delta < delta) max_delta = delta;


	}
	  *stdDev = sqrt(delta_sum/(double)num_point);
      *maxDev = max_delta;
   

	  // normalize the vector
       *x0 = a[1];
	   *y0 = a[2];
	   *z0 = a[3];
	   *R = a[4];

    delete lista;
	for(int i = 1; i <= ndata; i++) {
		delete x[i];                                                   		 
	}
    delete x;
    delete y;
    delete sig;
    delete a;

 
//	tmp.Format(_T("Fitted results for Sphere"));
//	pMainWnd->m_wndOutput.m_BuildList.AddString(tmp);

//		tmp.Format(_T("X0:%lf  Y0:%lf  Z0:%lf R:%1f, StdDev: %.3f maxDev:%.2f"),
//			*x0,*y0, *z0, *R, *stdDev, *maxDev);
//	pMainWnd->m_wndOutput.m_BuildList.AddString(tmp);
}



// 2D circle fitting with fixed circle radius
// by using nonlinear least square method
void CStdMath::Get2DCircleFit_NL(int num_point, double  *p[2], double init_par[3], double *x0, double *y0, double *R, double *stdDev, double *maxDev, bool bFixedRadius)
{
  
	// declaration & initialization
	int ndata = num_point;		// number of input points x[i][k]; i = 1,... ndata
	int nd = 2;					// nd: number of given x[i][k],  k = 1,2;
	int ma = 3;					// number of parameters to be determined after fitting
	int mfit;				   // = ma

	int* lista;					// lista[1..ma]
	double ** x;					// x[1..ndata][1..nd]	
	double * y;					// y[1..ndata]
	double * sig;					// sig[1..ndata]
	double * a;					// a[1..ma]


	if (bFixedRadius == true)
		mfit = 2;
	else 
		mfit = 3;
	
	// memory allocation
	lista = new int[ma+1];
	 
	x = new double*[ndata+1];
	for(int i = 1; i <= ndata; i++) {
		x[i] = new double[nd+1];                                                   		 
	}

	y = new double[ndata+1];
	sig = new double[ndata+1];
	a = new double[ma+1];
  
	for(int k = 1; k <= ma; k++) {
		lista[k] = k;			// numbers the parameters a such that 
 						// the first mfit elements correspond to values actually being adjusted
	}
  
	for(int k = 1; k <= ndata; k++) {
		sig[k] = 1;				// standard deviations, this is not experiment data, so set to 1
	    y[k] = 0.0;             // y = 0
	}

    // assign x value
    for(int i = 1; i <= ndata; i++) {
		for (int j = 1; j <= nd; j++){
			x[i][j] = (double) p[i-1][j-1];
		}
    	
	}


   // initial value of parameters
	a[1] = init_par[0];
	a[2] = init_par[1];
	a[3] = init_par[2];

	NonLinearLeastSquare(x, y, sig, ndata, a, ma, lista, mfit, CIRCLE_FIT);
    
	// calculate the fitting error 
	double max_delta = -10000.0;
    double dd, delta, delta_sq, delta_sum;
	
	delta_sum = 0.0;

   
    for (int i = 0; i< num_point; i++){
	    
            // fitting error to sphare center
			dd = sqrt((p[i][0]-a[1])*(p[i][0]-a[1])+(p[i][1]-a[2])*(p[i][1]-a[2])); 
			delta = fabs(dd - a[3]);
			delta_sq = (dd - a[3])*(dd - a[3]);
            delta_sum = delta_sum + delta_sq;
			if (max_delta < delta) max_delta = delta;


	}
	  *stdDev = sqrt(delta_sum/(double)num_point);
      *maxDev = max_delta;
   

	  // normalize the vector
       *x0 = a[1];
	   *y0 = a[2];
	   *R = a[3];

    delete lista;
	for(int i = 1; i <= ndata; i++) {
		delete x[i];                                                   		 
	}
    delete x;
    delete y;
    delete sig;
    delete a;
}


/*******************************************************************************************
 * GetCyliner() calculates the cylinder central axis and its radius based on input points 
 * by using nonlinear least sqaure method
 * Inputs: 
 *        p[][3]; containing a measured point array of [x,y,z];
 *        num_point; number of pints of the input; 
 * Outputs:  
 *         axis normal (nx, ny nz) and the center (x0,y0,z0) of input points that is on axis.
 *         and radius R 
 *         standard devaiation (stdDev) and maxmium deviation (maxDev) 
 *                    
 * return:
 *         void
 * 
 * Notes:
 * Since data fitting is obtained by nonlinear least method that requires initial values
 * To avoid not giving initial value certain procedure has to be followed during the measuremenet
 * Measurement approach
 * 1. measure points around one end of the cylinder. The measurements should be
 *    taken approximately on a section of the cylinder
 * 2. Measure points around the other end of the cylinder. The measuremnts should be
 *    taken approximately on a section of the cylinder
 * The initial [x,y,z] is the central points of the first measurement
 * The initial radius R is the average distance between the initial center and measured points
 * The intial axis normal is determined by the connection of two central points
 * Assume that the number of measurement points for two sections are equal that is num_point/2;
 * 
 * Author: 
 *          10.29.00 
 */

void CStdMath::GetCylinder(int num_point, double  *p[3], double *nx, double *ny, double *nz, double *x0, double *y0, double *z0, double *R, double *stdDev, double *maxDev)
{
  
	// declaration & initialization
	int ndata = num_point;		// number of input points x[i][k]; i = 1,... ndata
	int nd = 3;					// nd: number of given x[i][k],  k = 1,2 ,3;
	int ma = 7;					// number of parameters to be determined after fitting
	int mfit = 7;				// = ma

	int* lista;					// lista[1..ma]
	double ** x;					// x[1..ndata][1..nd]	
	double * y;					// y[1..ndata]
	double * sig;					// sig[1..ndata]
	double * a;					// a[1..ma]

	
	// memory allocation
	lista = new int[ma+1];
	 
	x = new double*[ndata+1];
	for(int i = 1; i <= ndata; i++) {
		x[i] = new double[nd+1];                                                   		 
	}

	y = new double[ndata+1];
	sig = new double[ndata+1];
	a = new double[ma+1];
  
	for(int k = 1; k <= ma; k++) {
		lista[k] = k;			// numbers the parameters a such that 
 						// the first mfit elements correspond to values actually being adjusted
	}
  
	for(int k = 1; k <= ndata; k++) {
		sig[k] = 1;				// standard deviations, this is not experiment data, so set to 1
	    y[k] = 0.0;             // y = 0
	}

    // assign x value
    for(int i = 1; i <= ndata; i++) {
		for (int j = 1; j <= nd; j++){
			x[i][j] = (double) p[i-1][j-1];
		}
    	
	}


    // estimate guess value of NLS fitting
	/* measure two sections around the cylinder one at top othet at bottom
	 * each with assigned points (i.e. 5 points for each)
	 * center of the first section is used as initial (x0,y0,z0) or a[1],a[2],a[3];
	 * average distance from center to measuremnt points are R or a[7]
	 * connection between two section centers are used as initial orientation
	 * or a[4],a[5],a[6];
	 */

	a[1] = 0.0;
	a[2] = 0.0;
	a[3] = 0.0;
	a[4] = 0.0;
	a[5] = 0.0;
	a[6] = 0.0;
	a[7] = 0.0;
    double xc2 = 0.0;
	double yc2 = 0.0;
	double zc2 = 0.0;
	double d1 = 0.0;

    
	// center of first section
	for (int i =1 ; i<= ndata/2; i++){
       a[1] = x[i][1] + a[1];
	   a[2] = x[i][2] + a[2];
	   a[3] = x[i][3] + a[3];
    }
	a[1]= 2*a[1]/ndata;
	a[2]= 2*a[2]/ndata;
	a[3]= 2*a[3]/ndata;

	// average distance to points
	for (int i =1 ; i<= ndata/2; i++){
	double  temp_d = (a[1] - x[i][1])*(a[1]-x[i][1]) + 
				(a[2] - x[i][2])*(a[2]-x[i][2]) +
				(a[3] - x[i][3])*(a[3]-x[i][3]);
	   temp_d = sqrt(temp_d);

       d1 = temp_d + d1;
    }
       d1 = 2*d1/ndata;
	   a[7] = d1;

	// center of sencond section
    for (int i = ndata/2+1 ; i<=ndata; i++){
       xc2 = x[i][1] + xc2;
	   yc2 = x[i][2] + yc2;
	   zc2 = x[i][3] + zc2;
    }
	xc2 = 2*xc2/ndata;
	yc2 = 2*yc2/ndata;
	zc2 = 2*zc2/ndata;
 
    
	// estimate normal
	double temp_l = (xc2-a[1])*(xc2-a[1]) +
		           (yc2-a[2])*(yc2-a[2]) +
				   (zc2-a[3])*(zc2-a[3]);

           temp_l = sqrt(temp_l);

          a[4] = (xc2-a[1])/temp_l;
		  a[5] = (yc2-a[2])/temp_l;
		  a[6] = (zc2-a[3])/temp_l;
    
   // estimate center
		a[1]= (a[1] + xc2)/2.0;
		a[2]= (a[2] + yc2)/2.0;
		a[3]= (a[3] + zc2)/2.0;
	


   // remove sigularity for perfect fitting
   	a[1] = a[1]*0.999999;
	a[2] = a[2]*1.00000001;
	a[3] = a[3]*0.9999999;
	a[4] = a[4]*0.999999;
	a[5] = a[5]*0.999999;
	a[6] = a[6]*0.999999;
	a[7] = a[7]*0.999999;
 
  
	NonLinearLeastSquare(x, y, sig, ndata, a, ma, lista, mfit, 2);
    // HOW TO USE NLLS:
	// (x,y): input point pair y= f(ax); normally y = 0;
	// sig: weighting factors. normally sig = 1;
	// ndata: number of data point.
	// a: coefficients to be determined -- outcome
	// ma: number of a
	// mfit: number of a to be fitted
	// lista: array of a [1 ...ma];



    
	// calculate the fitting error 

	double max_del_d = 0.0;
    double dev = 0;
    double d, tmp_ang, d_sq, del_d;
	
   	double norm_sq = a[4]*a[4] + a[5]*a[5] + a[6]*a[6];

    for (int i = 0; i< num_point; i++){
	  tmp_ang = (p[i][0]-a[1])*a[4]+(p[i][1]-a[2])*a[5]+(p[i][2]-a[3])*a[6];
      
	  d_sq = (p[i][0]-a[1])*(p[i][0]-a[1]) + (p[i][1]-a[2])*(p[i][1]-a[2]) + (p[i][2]-a[3])*(p[i][2]-a[3])
		    - tmp_ang*tmp_ang/norm_sq;
	  // distance between point x[i][] and axis 
	  d = sqrt(d_sq);
	  del_d = fabs(d-a[7]);
	  if (max_del_d < del_d) {
		    max_del_d = del_d;
	  }
      dev = dev + del_d;
	}
	  *stdDev = dev/(double)num_point;
      *maxDev = max_del_d;
   

	  // normalize the vector
       *x0 = a[1];
	   *y0 = a[2];
	   *z0 = a[3];
	   *nx = a[4]/sqrt(norm_sq);
	   *ny = a[5]/sqrt(norm_sq);
	   *nz = a[6]/sqrt(norm_sq);
	   *R = a[7];

    delete lista;
	for(int i = 1; i <= ndata; i++) {
		delete x[i];                                                   		 
	}
    delete x;
    delete y;
    delete sig;
    delete a;


}


/*******************************************************************************************
 * GetCyliner() calculates the cylinder central axis and its radius based on input points 
 * by using nonlinear least sqaure method
 * Inputs: 
 *        p[][3]; containing a measured point array of [x,y,z];
 *        num_point; number of pints of the input; 
 *        int_ori[3]:  initial ori4entation of cylinder axis [Nx, Ny, Nz]
 * Outputs:  
 *         axis normal (nx, ny nz) and the center (x0,y0,z0) of input points that is on axis.
 *         and radius R 
 *         standard devaiation (stdDev) and maxmium deviation (maxDev) 
 *                    
 * return:
 *         void
 * 
 * Notes:
 * Since data fitting is obtained by nonlinear least method that requires initial values
 * 
 * The initial [x,y,z] is the average point of all the measurement data
 * The initial radius R is the average distance between the initial center and measured points
 * 
 * Author: 
 *          05.25.04 
 */

void CStdMath::GetCylinder(int num_point, double  *p[3], double init_ori[3], double *nx, double *ny, double *nz, double *x0, double *y0, double *z0, double *R, double *stdDev, double *maxDev)
{
  
	// declaration & initialization
	int ndata = num_point;		// number of input points x[i][k]; i = 1,... ndata
	int nd = 3;					// nd: number of given x[i][k],  k = 1,2 ,3;
	int ma = 7;					// number of parameters to be determined after fitting
	int mfit = 7;				// = ma

	int* lista;					// lista[1..ma]
	double ** x;					// x[1..ndata][1..nd]	
	double * y;					// y[1..ndata]
	double * sig;					// sig[1..ndata]
	double * a;					// a[1..ma]

	
	// memory allocation
	lista = new int[ma+1];
	 
	x = new double*[ndata+1];
	for(int i = 1; i <= ndata; i++) {
		x[i] = new double[nd+1];                                                   		 
	}

	y = new double[ndata+1];
	sig = new double[ndata+1];
	a = new double[ma+1];
  
	for(int k = 1; k <= ma; k++) {
		lista[k] = k;			// numbers the parameters a such that 
 						// the first mfit elements correspond to values actually being adjusted
	}
  
	for(int k = 1; k <= ndata; k++) {
		sig[k] = 1;				// standard deviations, this is not experiment data, so set to 1
	    y[k] = 0.0;             // y = 0
	}

    // assign x value
    for(int i = 1; i <= ndata; i++) {
		for (int j = 1; j <= nd; j++){
			x[i][j] = (double) p[i-1][j-1];
		}
    	
	}


    // estimate guess value of NLS fitting
	/* measure two sections around the cylinder one at top othet at bottom
	 * each with assigned points (i.e. 5 points for each)
	 * center of the first section is used as initial (x0,y0,z0) or a[1],a[2],a[3];
	 * average distance from center to measuremnt points are R or a[7]
	 * connection between two section centers are used as initial orientation
	 * or a[4],a[5],a[6];
	 */


	a[1] = 0.0;
	a[2] = 0.0;
	a[3] = 0.0;
	a[4] = init_ori[0]+0.05;
	a[5] = init_ori[1]+0.05;
	a[6] = init_ori[2]+0.05;
	a[7] = 0.0;
    double xc2 = 0.0;
	double yc2 = 0.0;
	double zc2 = 0.0;
	double d1 = 0.0;

    
	// center of all measured points
	for (int i =1 ; i<= ndata; i++){
       a[1] = x[i][1] + a[1];
	   a[2] = x[i][2] + a[2];
	   a[3] = x[i][3] + a[3];
    }
	a[1]= a[1]/ndata;
	a[2]= a[2]/ndata;
	a[3]= a[3]/ndata;

	// average distance to points (should min distance)
	for (int i =1 ; i<= ndata/2; i++){
	double  temp_d = (a[1] - x[i][1])*(a[1]-x[i][1]) + 
				(a[2] - x[i][2])*(a[2]-x[i][2]) +
				(a[3] - x[i][3])*(a[3]-x[i][3]);
	   temp_d = sqrt(temp_d);

       d1 = temp_d + d1;
    }
       d1 = 2*d1/ndata;
	   a[7] = d1;


 
    
	NonLinearLeastSquare(x, y, sig, ndata, a, ma, lista, mfit, 2);
    // HOW TO USE NLLS:
	// (x,y): input point pair y= f(ax); normally y = 0;
	// sig: weighting factors. normally sig = 1;
	// ndata: number of data point.
	// a: coefficients to be determined -- outcome
	// ma: number of a
	// mfit: number of a to be fitted
	// lista: array of a [1 ...ma];



    
	// calculate the fitting error 

	double max_del_d = 0.0;
    double dev = 0;
    double d, tmp_ang, d_sq, del_d;
	
   	double norm_sq = a[4]*a[4] + a[5]*a[5] + a[6]*a[6];

    for (int i = 0; i< num_point; i++){
	  tmp_ang = (p[i][0]-a[1])*a[4]+(p[i][1]-a[2])*a[5]+(p[i][2]-a[3])*a[6];
      
	  d_sq = (p[i][0]-a[1])*(p[i][0]-a[1]) + (p[i][1]-a[2])*(p[i][1]-a[2]) + (p[i][2]-a[3])*(p[i][2]-a[3])
		    - tmp_ang*tmp_ang/norm_sq;
	  // distance between point x[i][] and axis 
	  d = sqrt(d_sq);
	  del_d = fabs(d-a[7]);
	  if (max_del_d < del_d) {
		    max_del_d = del_d;
	  }
      dev = dev + del_d;
	}
	  *stdDev = dev/(double)num_point;
      *maxDev = max_del_d;
   

	  // normalize the vector
       *x0 = a[1];
	   *y0 = a[2];
	   *z0 = a[3];
	   *nx = a[4]/sqrt(norm_sq);
	   *ny = a[5]/sqrt(norm_sq);
	   *nz = a[6]/sqrt(norm_sq);
	   *R = a[7];

    delete lista;
	for(int i = 1; i <= ndata; i++) {
		delete x[i];                                                   		 
	}
    delete x;
    delete y;
    delete sig;
    delete a;


}


/*******************************************************************************************
 * GetRotationAxisByPlanes() calculates the rotational and its center based on input points 
 * by using nonlinear least sqaure method
 * Inputs: 
 *        p[][8]; containing a measured point array of [xw1,yw1,zw1, theta, nx, ny, nz,d];
 *        num_point; number of pints of the input; 
 *        init_par[6]:
                   axis normal (nx, ny nz) and the center (x0,y0,z0) of input points that is on axis.
 *      
 * Outputs:  
 *         cal_par[6]: 
 *                    axis normal (nx, ny nz) and the center (x0,y0,z0) of input points that is on axis.
 *      
 *         standard devaiation (stdDev) and maxmium deviation (maxDev) 
 *                    
 * return:
 *         void
 * 
 * Notes:
 * 
 */

void CStdMath::GetRotationAxisByPlanes(int num_point, double  *p[8], double init_par[6], double cal_par[6], double *stdDev, double *maxDev)
{
  
	// declaration & initialization
	int ndata = num_point;		// number of input points x[i][k]; i = 1,... ndata
	int nd = 8;					// nd: number of given x[i][k],  k = 1,2 ,3;
	int ma = 6;					// number of parameters
	int mfit = 2;				// number of parameters to be determined after fitting

	int* lista;					// lista[1..ma]
	double ** x;					// x[1..ndata][1..nd]	
	double * y;					// y[1..ndata]
	double * sig;					// sig[1..ndata]
	double * a;					// a[1..ma]
    int k;
	double chisq, alamda;
	
	// memory allocation
	lista = new int[ma+1];
	 
	x = new double*[ndata+1];
	for(int i = 1; i <= ndata; i++) {
		x[i] = new double[nd+1];                                                   		 
	}

	y = new double[ndata+1];
	sig = new double[ndata+1];
	a = new double[ma+1];
  
//	for(k = 1; k <= ma; k++) {
//		lista[k] = k;			// numbers the parameters a such that 
					// the first mfit elements correspond to values actually being adjusted
//	}


	lista[1] = 1;    //X0
	lista[2] = 2;    //Y0
	lista[3] = 4;    //nx
	lista[4] = 5;    //ny
	lista[5] = 6;    //nz
	lista[6] = 3;;   //Z0

  
	for(k = 1; k <= ndata; k++) {
		sig[k] = 1;				// standard deviations, this is not experiment data, so set to 1
	    y[k] = 0.0;             // y = 0
	}

    // assign x value
    for(int i = 1; i <= ndata; i++) {
		for (int j = 1; j <= nd; j++){
			x[i][j] = (double) p[i-1][j-1];
		}
    	
	}

 
//	a[1] = init_par[0];
//	a[2] = init_par[1];
//	a[3] = init_par[2];
//	a[4] = init_par[3];
//	a[5] = init_par[4];
//	a[6] = init_par[5];
	
   
  
//	NonLinearLeastSquare(x, y, sig, ndata, a, ma, lista, mfit, 11);
	
	double **cvar = new double*[ma + 1];
	for (int i = 0; i < ma+1; i++)
		cvar[i] = new double[ma + 1];

	double **alpha = new double*[ma + 1];
	for (int i = 0; i < mfit+1; i++)
		alpha[i] = new double[ma + 1];

    


    // Use different initial condition to void falling into local maxima	
#define NUM_TRY 10
	int NOISE_AMPLITUDE=3;

    double chisq_min = 1000000;
    int flag;

	for (int ii=0; ii<NUM_TRY; ii++){
		int idum=10;
	    ran0(&idum);  //initialize random number generator
        
		if ( ii==0){ 
		  flag = 0;                 // make sure to try initial value
		} else {
          flag = 1;
		}


		if (ii == 1){                 // best guss
			a[1] = init_par[0] + 2.5;
			a[2] = init_par[1] + 2;
		} else {
			a[1] = init_par[0] + flag*NOISE_AMPLITUDE*(ran0(&idum) - 0.5)*2;
			a[2] = init_par[1] + flag*NOISE_AMPLITUDE*(ran0(&idum) - 0.5)*2;
		}
		a[3] = init_par[2];
		a[4] = init_par[3];
		a[5] = init_par[4];
		a[6] = init_par[5];
		
		chisq = 0;
		alamda = -1;   
		
		NonLinearLeastSquare(x, y, sig, ndata, a, ma, lista, mfit, cvar, alpha, &chisq, &alamda, 11);

		if (chisq_min > chisq)  {
				chisq_min = chisq;
				cal_par[0] = a[1];         // X0
				cal_par[1] = a[2];         // Y0
				cal_par[2] = a[3];         // Z0
				cal_par[3] = a[4];		   // nx
				cal_par[4] = a[5];         // ny
				cal_par[5] = a[6];		   // nz
		}
	}


	
// HOW TO USE NLLS:
// (x,y): input point pair y= f(ax); normally y = 0;
// sig: weighting factors. normally sig = 1;
// ndata: number of data point.
// a: coefficients to be determined -- outcome
// ma: number of a
// mfit: number of a to be fitted
// lista: array of a [1 ...ma];

  	 // normalize the vector
	double den = sqrt(cal_par[3]*cal_par[3]+ cal_par[4]*cal_par[4] + cal_par[5]*cal_par[5]);
	cal_par[3] = cal_par[3]/den;
	cal_par[4] = cal_par[4]/den;
	cal_par[5] = cal_par[5]/den;


//	*stdDev = (double) sqrt(chisq/ndata);
	*stdDev = (double) sqrt(chisq_min/ndata);


   	for(int i = 1; i <= ndata; i++) {
		delete x[i];                                                   		 
	}


	for (int i = 0; i < ma+1; i++){
		delete []cvar[i];
	}
    delete []cvar;

	for (int i = 0; i < mfit+1; i++){
		delete []alpha[i];
	} delete []alpha;

    delete []lista;
	delete []x;
    delete []y;
    delete []sig;
    delete []a;


}


/*******************************************************************************************
 * CircleFitting_2D () calculates the circle parameters based on input points on the circle 
 * Inputs: 
 *        cir[][3]; containing a array of [x,y];
 *        num_point; number of pints of the input; 
 * 
 * Outputs:  
 *           parameters of Circle parameters including 
 *           center position (xc, yc), 
 *           radius R
 *           
 *          
 * return:
 *         void
 * 
 * Notes for approach
 * 1. The general equation for a circle is
 * x^2 + y^2 + a1* x + a2* y + a3 = 0;   -- eq(1)
 * where the center and radius are given by 
 *  xc = -a1/2;  yx= -a2/2; R = sqrt[(a1*a1 + a2*a2)/4 - a3];
 * Eq(1) can be solved by solving LS eqaution
 * AX = B
 * where X[3] = [a1 a2 a3]';
 *       A[M][3] is a Mx3 matrix and its row is made of each [x y 1] 
 *       B [M]  is a Mx1 matrix and its row (one element) is made of [-x^2-y^2]

 * Author: 
 *          6.9.00 
 */
void CStdMath::CircleFitting_2D (int num_point, double *cir[2], double *xc, double *yc, double *R)
{
    double **A;
	double *B;
	double X[3];
	int i;

	A = new double*[num_point];
	B = new double[num_point];

	
    for(int i = 0; i < num_point; i++) {
	  A[i] = new double[3];                                                   
	}
	
	for (int i=0; i<num_point; i++){
		B[i] = 1;
	}

	// Build A , B matrix
	for (int i=0; i<num_point; i++){
		A[i][0] = cir[i][0];
		A[i][1] = cir[i][1];
		A[i][2] = 1.0;
		B[i] = -cir[i][0]*cir[i][0] - cir[i][1]*cir[i][1];
	}

	LinearLeastSquare(num_point, 3, A, B, X);

   // circle parameters
   *xc = -X[0]/2.0; 
   *yc = -X[1]/2.0; 
   *R = sqrt((X[0]*X[0] + X[1]*X[1])/4.0 - X[2]);


   delete []B;
   for(int i = 0; i < num_point; i++) {
	  delete []A[i];                                                   
   }
   delete []A;
}

// Calculate the intersectional point between two lines
// Input: two line parameters
//       double line1_par[6] =  (x1, y1, z1, nx1, ny1, nz1)
//       double line2_par[6] =  (x2, y2, z2, nx2, ny2, nz2)
// Output: intersectional point that has the minimal distance to two lines 
//       double intersect[3] = (x0, y0, z0) 
//       double *err (distance from the intersectional point to the line)
//
// Note 
//     Two line equations are
//     x = x1+nx1*t1               (1)
//     y = y1+ny1*t1
//     z = z1+nz1*t1
//     x = x2+nx2*t2
//     y = y2+ny2*t2
//     z = z2+nz2*t2
// The intersectional points are the solution of the overdetermined equation
//     nx1*t1 - nx2*t2 = x2-x1;
//     ny1*t1 - ny2*t2 = y2-y1;
//     nz1*t1 - nz2*t2 = z2-z1;
//  (t1 and t2) can be solved wby linear LS with the form AX = B
//   where A = [nx1  -nx2      B = [ x2-x1
//              ny1  -ny2            y2-y1
//              nz1  -nz2]           z2-z1]
//
// Then intersect point can be obtainde from equation (1) 
//  02.03.05
void CStdMath::Line2LineIntersect(double line1_par[6], double line2_par[6], double intersect[3], double *err)
{
	double x1, y1, z1, nx1,ny1, nz1;
	double x2, y2, z2, nx2,ny2, nz2;
    double **A;
	double *B;
	double X[2];
	int num_row = 3;
	int i;


	x1 = line1_par[0];
	y1 = line1_par[1];
    z1 = line1_par[2];
	nx1 = line1_par[3];
	ny1 = line1_par[4];
    nz1 = line1_par[5];
	

	x2 = line2_par[0];
	y2 = line2_par[1];
    z2 = line2_par[2];
	nx2 = line2_par[3];
	ny2 = line2_par[4];
    nz2 = line2_par[5];
	
	A = new double*[num_row];
	B = new double[num_row];
    for(int i = 0; i < num_row; i++) {
		A[i] = new double[2];                                                   
	}
	
	// Build A , B matrix
	//   where A = [nx1  -nx2      B = [ x2-x1
	//              ny1  -ny2            y2-y1
	//              nz1  -nz2]           z2-z1]
	A[0][0] =  nx1;
	A[0][1] = -nx2;
	A[1][0] =  ny1;
	A[1][1] = -ny2;
	A[2][0] =  nz1;
	A[2][1] = -nz2;
	
	B[0] = x2-x1;
	B[1] = y2-y1;
	B[2] = z2-z1;

	LinearLeastSquare(num_row, 2, A, B, X);

	intersect[0] = x1 + nx1*X[0];
	intersect[1] = y1 + ny1*X[0];
	intersect[2] = z1 + nz1*X[0];

	double d;

//  PointLineDistance(intersect[0], intersect[1], intersect[2], x1, y1, z1, nx1, ny1, nz1, &d);
   	PointLineDistance(intersect[0], intersect[1], intersect[2], x2, y2, z2, nx2, ny2, nz2, &d); //Dai Liang 2007.1.25
	*err = d;

    // clean up mem
	delete []B;
	for(int i = 0; i < num_row; i++) {
		delete []A[i];                                                   
	}
	delete []A;


}

// Return standard dev and max error
void CStdMath::CircleFitting_2D (int num_point, double *cir[2], double *xc, double *yc, double *R, double *stdDev, double *maxDev)
{
    double **A;
	double *B;
	double X[3];
	int i;

	A = new double*[num_point];
	B = new double[num_point];

	
    for(int i = 0; i < num_point; i++) {
	  A[i] = new double[3];                                                   
	}
	
	for (int i=0; i<num_point; i++){
		B[i] = 1;
	}

	// Build A , B matrix
	for (int i=0; i<num_point; i++){
		A[i][0] = cir[i][0];
		A[i][1] = cir[i][1];
		A[i][2] = 1.0;
		B[i] = -cir[i][0]*cir[i][0] - cir[i][1]*cir[i][1];
	}

	LinearLeastSquare(num_point, 3, A, B, X);

   // circle parameters
   *xc = -X[0]/2.0; 
   *yc = -X[1]/2.0; 
   *R = sqrt((X[0]*X[0] + X[1]*X[1])/4.0 - X[2]);



   	// calculate starndard deveation and max devation
	  
	  double maxError, sum, delta, dis;
	  maxError = -100000;
	  sum = 0;

	  for (	i=0; i< num_point; i++){
		   dis = sqrt((cir[i][0] - *xc)*(cir[i][0]-*xc) + (cir[i][1] - *yc)*(cir[i][1]-*yc));
		   delta = fabs(dis-*R);
		   if (maxError < delta) maxError = delta;
		   sum = sum + delta;
	  }
	  *stdDev = sum/num_point;
	  *maxDev = maxError;



   delete []B;
   for(int i = 0; i < num_point; i++) {
	  delete []A[i];                                                   
   }
   delete []A;
}

/*************************************************	

Function Name: PT_IJK_to_PT_RPY

Functionality:	

Input: Point p1, Direction Vector, Normal Vector 
		 
Return: PT_RPY

***************************************************/

PT_RPY CStdMath::PTIJK_to_PTRPY( double P1[], double d[], double n[])
{
	double xx, xy, xz, zx, zy, zz, yx, yy, yz,  mag;
	double roll, pitch, yaw;
	PT_RPY pt;

	// Vector in X-Z Plane
	
	xx = d[0];
	xy = d[1];
	xz = d[2];
	mag = sqrt(xx*xx+xy*xy+xz*xz);
	xx = xx / mag;
	xy = xy / mag;
	xz = xz / mag;

	// Z-axis
	zx = n[0];
	zy = n[1];
	zz = n[2];
	mag = sqrt(zx*zx+zy*zy+zz*zz);
	zx = zx / mag;
	zy = zy / mag;
	zz = zz / mag;

	// Y-axis
	yx = zy*xz - zz*xy;
	yy = zz*xx - zx*xz;
	yz = zx*xy - zy*xx;
	mag = sqrt(yx*yx+yy*yy+yz*yz);
	yx = yx / mag;
	yy = yy / mag;
	yz = yz / mag;

	// X-axis
	xx = yy*zz - yz*zy;
	xy = yz*zx - yx*zz;
	xz = yx*zy - yy*zx;
	mag = sqrt(xx*xx+xy*xy+xz*xz);
	xx = xx / mag;
	xy = xy / mag;
	xz = xz / mag;

	// Roll, Pitch, and Yaw Computation
	pitch = atan2(-xz,sqrt(xx*xx + xy*xy) );
	
	if ((pitch*180/PI < 90 + DELTA) & (pitch*180/PI > 90 - DELTA))
	{
		roll  = 0;
		pitch = PI/2;
		yaw   = atan2(-yx, yy);
	}
	else if ((pitch*180/PI > -90 - DELTA) & (pitch*180/PI < -90 + DELTA))
	{
		roll  = 0;
		pitch = -PI/2;
		yaw   = -atan2(-yx, yy);
	}
	else 
	{
		roll  = atan2(xy, xx);
		yaw   = atan2(yz, zz);
	}

	pt.roll  = roll * 180/PI;
	pt.pitch = pitch * 180/PI;
	pt.yaw   = yaw * 180/PI;

	// Origin
	pt.x  = P1[0];
	pt.y  = P1[1];
	pt.z  = P1[2];

	/****************************************************
		Roll  Represent Rotation around Z in ABB Robot
		Pitch Represent Rotation around Y in ABB Robot
		Yaw   Represent Rotation around X in ABB Robot
	*****************************************************/

	return pt;
}


/***************************************************************
 * CalculateRotAxisBy2Mat(double m1[3][3], double m1[3][3], double n[3])
 * 
 * Based on the two rotation matrix
 * calculate the rotation axis
 * input: m1[3][3], m2[3][3] rotation matrix
 * output: rotation axis n[3] = {nx, ny, nz}
 * return: rotation angle
 *
 * 6/13/01
 ***************************************************************/

double CStdMath::CalculateRotAxisBy2Mat(double m1[3][3], double m2[3][3], double n[3])
{
	double m3[3][3], t[3][3];
	double **M1 = new double *[3];
	for (int i = 0; i < 3; i++)
	{
		M1[i] = t[i];
		for (int j = 0; j < 3; j++)
			M1[i][j] = m1[i][j];
	}

	MatrixInverse(3, M1);

	MatrixMultiply33(t, m2, m3);

	delete M1;

	CalculateRotAxis(m3, n);

	double cos_beta = (m3[0][0] + m3[1][1] + m3[2][2] - 1) / 2;
	return acos(cos_beta);
	//return 1;
}


/***************************************************************
 * CalculateRotAxis(double m1, double n[3])
 * 
 * Based on the rotation matrix
 * calculate the rotation axis
 * input: m1[3][3] rotation matrix
 * output: rotation axis n[3] = {nx, ny, nz}
 *
 * Note:
 *       This function is based on Eq. (5.34) of book
 *       Computer Graphics using Open GL by F.S. Hill.
 * 
 * 6/9/01
 ***************************************************************/

int CStdMath::CalculateRotAxis(double m[3][3], double n[3])
{
	double cos_beta = (m[0][0] + m[1][1] + m[2][2] - 1) / 2;
	double sin_beta = sin(acos(cos_beta));
	n[0] = (m[2][1] - m[1][2]) / ( 2 * sin_beta);
	n[1] = (m[0][2] - m[2][0]) / ( 2 * sin_beta);
	n[2] = (m[1][0] - m[0][1]) / ( 2 * sin_beta);
	double norm = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	for (int i = 0; i < 3; i++)
		n[i] /= norm;                  // keep directrion toward dow
	return 1;
}


/***************************************************************
 * CalculateRotCenter(double p1[3], double p2[3], double n[3], double c[3])
 * 
 * Based on the position of two points (one is a rotation of another)
 * and rotation axis, calculate the rotation center
 * assuming that the rotation angle is 2*beta
 * input: p1[3], p2[3], n[3]  (two points and one orientation)
 * output: c[3]               (rotation center) 
 * 
 * 
 * 
 * 
 * 6/9/01
 ***************************************************************/

int  CStdMath::CalculateRotCenter(double p1[3], double p2[3], double n[3], double c[3], double beta)
{
  
	double p3[3];   // point between p1 and p2
	double dis;     // distance between p1 and p2
	double nr[3];   // vector toward rotation center
	double nd[3];   // vector formed by p1-p2

	
	p3[0] = (p1[0] + p2[0])/2;
	p3[1] = (p1[1] + p2[1])/2;
	p3[2] = (p1[2] + p2[2])/2;


	dis = sqrt((p1[0] - p2[0])*(p1[0] - p2[0]) + (p1[1] - p2[1])*(p1[1] - p2[1]) + (p1[2] - p2[2])*(p1[2] - p2[2]));

   // p1 is on the left of p2

	nd[0] = (p2[0]-p1[0])/dis;
	nd[1] = (p2[1]-p1[1])/dis;
	nd[2] = (p2[2]-p1[2])/dis;

	// n is pointing up and nr is pointing out
	//nr = nd X n
	crossProd(nd, n, nr);

	c[0] = p3[0] + nr[0] * dis /2 / tan(beta);
	c[1] = p3[1] + nr[1] * dis /2 / tan(beta);
	c[2] = p3[2] + nr[2] * dis /2 / tan(beta);

//  c[0] = p3[0] - nr[0] * dis /2 / tan(beta);
//  c[1] = p3[1] - nr[1] * dis /2 / tan(beta);
//  c[2] = p3[2] - nr[2] * dis /2 / tan(beta);


	return 1;

}


/***************************************************************
 * CalculateCoodinatefrom2Vectors(double vector_z[3], double orig[3], double vector_x[3], t[4][4])
 * 
 * Based on z axis vector_y[3] and origin orig[3], and x ditection vector_x[3] 
 * deteremine the table coodinate system described by a 4x4 transform matrix
 * input: vector_y[3], orig[3]  ( rotation center)
 *        vector_x[3],          ( x direction)
 * output: t[4][4];  
 * 
 * 7/26/01
 ***************************************************************/

int  CStdMath::CalculateCoordinatefrom2Vectors(double vector_y[3], double orig[3], double vector_x[3], double t[4][4])
{
   double y_norm, x_norm;
   double vector_z[3], z_norm;

   // normalize the input
   y_norm = sqrt(vector_y[0]*vector_y[0] + vector_y[1]*vector_y[1] + vector_y[2]*vector_y[2]);
   vector_y[0] = vector_y[0]/y_norm;
   vector_y[1] = vector_y[1]/y_norm;
   vector_y[2] = vector_y[2]/y_norm;
   x_norm = sqrt(vector_x[0]*vector_x[0] + vector_x[1]*vector_x[1] + vector_x[2]*vector_x[2]);
   vector_x[0] = vector_x[0]/x_norm;
   vector_x[1] = vector_x[1]/x_norm;
   vector_x[2] = vector_x[2]/x_norm;

   // calculate z direction  z = vector_x X vector_y
   crossProd(vector_x, vector_y, vector_z);
	
   // normalize y vector
   z_norm = sqrt(vector_z[0]*vector_z[0] + vector_z[1]*vector_z[1] + vector_z[2]*vector_z[2]);
   vector_z[0] = vector_z[0]/z_norm;
   vector_z[1] = vector_z[1]/z_norm;
   vector_z[2] = vector_z[2]/z_norm;


   // create a rotation matrix m[3][3] of table in camera coodinate (table to camera)
   // m[3][3] = [vector_x vector_y vector_z]
   // the relation between table and camera coodinate is 
   // Xc = m*Xt + Orig
  
   double **inv_m;
   int i, j;

   inv_m = new double*[4];
   for (int i=0; i<4; i++){
	   inv_m[i] = new double[4];
   }

    for (int i=0; i<3; i++){
	   inv_m[i][0] = vector_x[i];
	   inv_m[i][1] = vector_y[i];
	   inv_m[i][2] = vector_z[i];
	}

//    for (i=0; i<3; i++){
//	   inv_m[0][i] = vector_x[i];
//	   inv_m[1][i] = vector_y[i];
//	   inv_m[2][i] = vector_z[i];
//	}


	inv_m[0][3] = orig[0];  //x
	inv_m[1][3] = orig[1];  //y
    inv_m[2][3] = orig[2];  //z
	inv_m[3][0] = 0;
	inv_m[3][1] = 0;
	inv_m[3][2] = 0;
	inv_m[3][3] = 1;


   // inverse the matrix
   // MatrixInverse(4, inv_m);

   for (int i=0; i<4; i++){
	   for (j = 0; j<4; j++) {
			t[i][j] = inv_m[i][j];
	  }
   }

   //cleanup up
   for (int i=0; i<4; i++){
	  delete  []inv_m[i];
   }
   delete []inv_m;

  

   return 1;


}



 /*******************************************************************************************
 * LinePlaneIntersection() calculates the intersection point between a line and a plane
 * inputs:  a line passing throught (xl,yl,zl) point and with dirction (nlx,nly,nlz)
 *          a plane (npx,npy,npz, d) with equation npx*x+npy*x + npz*x =d
 * 
 * Outputs:  
 *          Intersection point coordinate (xi,yi,zi)  
 * Notes:
 *         The line equation is
 *         x = xl + nlx*t
 *         y = yl + nly*t
 *         z = z1 + nlz*t                        eq.(1)
 *         The plane equation is
 *         npx*x + npy*y + npz*z = d             eq.(2)
 *         Substituting eq.(1) into eq.(2) results in
 *         t = [d - (xl*npx+yl*npy+zl*npz)]/[nlx*npx + nly*npy + nlz*npz];
 *        Then the intersection point is given by eq.(1)
 *          
 * return:
 * 
 * Author: 
 *          11.15.00 
 */

void CStdMath::LinePlaneIntersection(double xl, double yl, double zl, double nlx, double nly, double nlz, double npx, double npy, double npz, double d, double *xi, double *yi, double *zi)
{
 
	double t, ang;

	ang = nlx*npx + nly*npy + nlz*npz;
	if(ang==0)
	{
		*xi=0;*yi=0;*zi=0;
		return;
	}
		
	t  = (d - (xl*npx+yl*npy+zl*npz))/ang;

	*xi = xl + nlx*t;
	*yi = yl + nly*t;
	*zi = zl + nlz*t;
}

/*******************************************************************************************
 * PointLineProjection() calculates the projected coordinate on a line
 * inputs: (x0,y0,z0) point coordinate that will be projected to a line
 *          a line passing through the point  (xl,yl,zl) and  with direction (nx,ny,nz) 
 *         
 * 
 * Outputs:  
 *          Projected point coordinate (xp,yp,zp)  
 *           
 *          
 * return:void
 *
 * Note:   The projected point (xp, yp, zp) can be obtainde by intersecting the line and 
 *         the plane that passies through (x0,y0,z0) and has normal (nx,ny,nz)
 *         The plane equation is given by
 *         nx*(x-x0) + ny*(y-y0) + nz*(z-z0) = 0
 *         or in its standard form nx*x + ny*y + nz*z = d
 *         where d = nx*x0 + ny*y0m + nz*z0;
 * 
 * Author: 
 *          11.15.00 
 */

void CStdMath::PointLineProjection(double x0, double y0, double z0, double x1, double y1, double z1, double nx, double ny, double nz,double *xp, double *yp, double *zp)
{
   double d;

   d = nx*x0 + ny*y0 + nz*z0;

   // call line plane intersection function
   LinePlaneIntersection(x1,y1,z1,nx,ny,nz, nx,ny,nz,d,xp,yp,zp);
}


/*******************************************************************************************
 * PointLineDistance() calculates the distance between a point to a line
 * inputs: (x0,y0,z0):  reference point coordinate 
 *         (xl, yl, zl, nx, ny, nz):  a line passing through the point  (xl,yl,zl)
 *                                     and  with direction (nx,ny,nz) 
 *         
 * Outputs:  
 *          double d: calculated distance 
 *          
 * return:void
 * 
 * Author: 
 *          10.11.02
 */

void CStdMath::PointLineDistance(double x0, double y0, double z0, double x1, double y1, double z1, double nx, double ny, double nz,double *d)
{
   double xp, yp, zp;
   double d_sq;

   // deteremine the projected point
   PointLineProjection(x0, y0, z0, x1, y1, z1, nx, ny, nz,  &xp, &yp, &zp);

   // calculate distance 
   d_sq = (x0-xp)*(x0-xp)+ (y0-yp)*(y0-yp)+(z0-zp)*(z0-zp);

   *d = sqrt(d_sq);
   
}



/*******************************************************************************************
 * PointLineDistance_2D() calculates the distance between a point  to a line
 * inputs: (x0,y0) point coordinate that will be projected to a line
 *         (a,b) line coefficients
 *         line_mode: if line_mode = 0 --> line is represented as (x = a*y + b) 
 *                    if line_mode = 1 --> line is represented as (y = a*x + b)         
 * Outputs:  
 *          double d: calculated distance 
 *                     
 * return:void
 * Note:   
 *   Distance between (x0, y0) to yhe line ax+by+c = 0 is
 *          (a*x0 + b*y0 +c)/sqrt(a^2 + b^2); 
 * 
 * 
 * Author: 
 *          10.11.02
 */

void CStdMath::PointLineDistance_2D(double x0, double y0, double line_mode, double a, double b, double *d)
{
    // if line_mode = 0 --> x = a*y + b
	if (line_mode == 0){
		*d = fabs((x0-a*y0-b)/sqrt(1+a*a));
	} else {      // y = a*x+b
		*d = fabs((a*x0 - y0 +b)/sqrt(1+a*a));
	}
   
}

/*******************************************************************************************
 * PointProjection() calculates the projected coordinate on the plane
 * inputs: pint coordinate thta will be projected and plane eqation
 *          (x1,y1,z1) and (nx*x + ny*y + nz*z = d)
 *  
 * 
 * Outputs:  
 *          Projected point coordinate (x2,y2,z2)  
 *           
 *          
 * return:
 * 
 * Author: 
 *          6.9.00 
 */

void CStdMath::PointProjection(double x1, double y1, double z1, double nx, double ny, double nz, double d, double *x2, double *y2, double *z2)
{
 
	double t;

	t = d - (nx*x1 + ny*y1 + nz*z1);

	*x2 = x1 + nx*t;
	*y2 = y1 + ny*t;
	*z2 = z1 + nz*t;

}

/*************************************************************************
 * EllipseFitting_2D () calculates the ellipse parameters based on input 
 * points, it does a least square fitting and returns the fitting ellipse
 * parameters
 * Inputs: 
 *        ellipse[][2]; containing a array of [x,y];
 *        num_point; number of pints of the input; 
 * 
 * Outputs:  
 *           parameters of fitted ellipse, including 
 *           center position (xc, yc), 
 *           long axis: longAxi
 *			 short axis: shortAxi
 *           orientation angle: Theta
 *          
 * return:
 *         void
 * 
 * Notes for approach
 * 1. The general equation for a ellipse is
 *  x^2 + a(1)xy + a(2)y^2 + a(3)x + a(4)y + a(5) = 0 -- eq(1)
 *  
 * Eq(1) can be solved by solving LS eqaution
 * AX = B
 * where X[5] = [a1 a2 a3 a4 a5]';
 *       A[M][5] is a Mx3 matrix and its row is made of each [x y 1] 
 *       B [M]  is a Mx1 matrix and its row (one element) is made of [-x^2]

 * Author: ABB
 *         
 */
void CStdMath::EllipseFitting_2D(int		numPoint, 
								double *ellipse[2], 
								int     fittedPoint,
								double *newEllipse[2],
								double *xc, 
								double *yc, 
								double *xAxi, 
								double *yAxi, 
								double *Theta)
{
    double **A;
	double *B;
	double X[5];
	int i;
	double d1[2];
	double dx;
	double x;
	double y;
	double ang;

	double  ct;
	double  st;
	
	double  a1;
	double  t;
	double  tc;

	A = new double*[numPoint];
	B = new double[numPoint];

	// Build A , B matrix	
    for(int i = 0; i < numPoint; i++) 
	{
	  A[i] = new double[5];
	  A[i][0] = ellipse[i][0]*ellipse[i][1];
	  A[i][1] = ellipse[i][1]*ellipse[i][1];
	  A[i][2] = ellipse[i][0];
	  A[i][3] = ellipse[i][1];
	  A[i][4] = 1.0;

	  B[i] = - ellipse[i][0]*ellipse[i][0];
	}
	
  LinearLeastSquare(numPoint, 5 , A, B, X);

  // reallocate memory for ellipse to have fittedPoint memory

   // Ellipse parameters

 //get scaled major/minor axes
   *Theta = atan2(X[0], 1-X[1])/2;
   tc     = tan(*Theta);
   t      = sqrt((2*tc-X[0])/(X[0]*tc*tc+2*tc));

   ct    = cos(*Theta);
   st    = sin(*Theta);
   a1    = ct*ct + t*t*st*st;
	
	*xc = (a1*X[2]*ct+a1*X[3]*st)/(-2.0);
	*yc = (a1*X[2]*st-a1*X[3]*ct)/(2.0*t*t);

    // get major/minor axis radii
     *xAxi  = sqrt((*xc)*(*xc)+t*t*(*yc)*(*yc)-a1*X[4]);
     *yAxi = (*xAxi)/t;
   
	 delete []B;
	 for(int i = 0; i < numPoint; i++) {
		delete []A[i];                                                   
	 }
	 delete []A;

	// replace original ellipse data with fitted data
	 dx = 2.0*PI/fittedPoint;

	 for(int i = 0; i<fittedPoint; i++)
	{
	  ang = (i+1)*dx;
	  x = (*xAxi)*cos(ang)+(*xc);
	  y = (*yAxi)*sin(ang)+(*yc);
	  d1[0] = x*cos(*Theta)-y*sin(*Theta);
	  d1[1] = x*sin(*Theta)+y*cos(*Theta);
	  newEllipse[i][0] = d1[0];
	  newEllipse[i][1] = d1[1];
	}
}

void CStdMath::EllipseFitting_2D(int numPoint, 
								double *ellipse[2], 
								int fittedPoint,
								double *newEllipse[2])
{
	double **A;
	double *B;
	double X[5];
	double d1[2];
	double x;
	double y;
	double ang;
	int i;

	double  ct;
	double  st;
	double  xAxi;
	double  yAxi;
	double  Theta;
	double  xc;
	double  yc;
	double  dx;
	
	double  a1;
	double  t;
	double  tc;

	A = new double*[numPoint];
	B = new double[numPoint];

	// Build A , B matrix	
    for(int i = 0; i < numPoint; i++) 
	{
	  A[i] = new double[5];
	  A[i][0] = ellipse[i][0]*ellipse[i][1];
	  A[i][1] = ellipse[i][1]*ellipse[i][1];
	  A[i][2] = ellipse[i][0];
	  A[i][3] = ellipse[i][1];
	  A[i][4] = 1.0;

	  B[i] = - ellipse[i][0]*ellipse[i][0];
	}
	
  LinearLeastSquare(numPoint, 5 , A, B, X);

   // Ellipse parameters

 //get scaled major/minor axes
   Theta = atan2(X[0], 1-X[1])/2;
   tc     = tan(Theta);
   t      = sqrt((2*tc-X[0])/(X[0]*tc*tc+2*tc));

   ct    = cos(Theta);
   st    = sin(Theta);
   a1    = ct*ct + t*t*st*st;

	xc = (a1*X[2]*ct+a1*X[3]*st)/(-2.0);
	yc = (a1*X[2]*st-a1*X[3]*ct)/(2.0*t*t);


    // get major/minor axis radii
     xAxi  = sqrt(xc*xc+t*t*yc*yc-a1*X[4]);
     yAxi = xAxi/t;
   
	 delete []B;
	 for(int i = 0; i < numPoint; i++) {
		delete []A[i];                                                   
	 }
	 delete []A;

	 // replace original ellipse data with fitted data
	 dx = 2.0*PI/fittedPoint;

	 for(int i = 0; i<fittedPoint; i++)
	{
	  ang = (i+1)*dx;
	  x = xAxi*cos(ang)+xc;
	  y = yAxi*sin(ang)+yc;
	  d1[0] = x*cos(Theta)-y*sin(Theta);
	  d1[1] = x*sin(Theta)+y*cos(Theta);
	  newEllipse[i][0] = d1[0];
	  newEllipse[i][1] = d1[1];
	}
	
}


// Rotate the input points Point3D[]
// into the plane that is in parellel with xy plane (Z value are the same)
// input: points in 3D space (they may lay on a plane): Point3D[]
// output: points after transform: *Point3D[] in which all Z are the same
//         its corresponding 2D coordinate is *Point2D[]
//         note that the rotation angles are roll/pitch/yaw
void CStdMath::Transfer3Dto2D(int numPoint, double *Point3D[], double *Point2D[])
{
	double temp_x,temp_y,temp_z,x1,y1,z1,x2,y2,z2,x,y,z;
 	double dis;
    double real_nx, real_ny, real_nz;
	double p[3], d[3], n[3];
    PT_RPY	Coord;

	double nx,ny,nz;
	GetPlaneNormal(numPoint, Point3D, &nx, &ny, &nz, &dis);
    // normalize the orientation vector in case they are not normalized
	real_nx = nx/sqrt(nx*nx + ny*ny+ nz*nz);
	real_ny = ny/sqrt(nx*nx + ny*ny+ nz*nz);
	real_nz = nz/sqrt(nx*nx + ny*ny+ nz*nz);


	/* 
	 * transform the input points into the plane that is in parellel with xy plane
	 */
	// assume that this vector pass origin
	// since we only consider the orientation
	p[0] = 0.0;
	p[1] = 0.0;
	p[2] = 0.0;

	n[0] = real_nx;
	n[1] = real_ny;
	n[2] = real_nz;

	/* 
	 * Need determine (give) direction vector to create RPY or coordinate
	 * since ijy info is not enogh to determine coordinate
	 */
	 
	// Assume direction vector d[] is normal to n[] or (i,j,k)
	if (real_nx != 0.0){
	  d[0] = -real_ny*0.5/real_nx;
	  d[1] = 0.5;
	  d[2] = 0.0;
	} else if (real_ny != 0.0){
        d[0] = 0.5;
		d[1] = -real_nx*0.5/real_ny;
	    d[2] = 0.0;
    } else { // direction is x 
		d[0] = 1.0;
	    d[1] = 0.0;
		d[2] = 0.0;
	}


	// Define a coordinate for circle based on above normal and direction vector
	Coord = PTIJK_to_PTRPY(p, d, n);


    /*
	 * roll, pitch & yaw indicate how the designed circle coordinate transfers 
	 *into world coordinate. 
	 * From world to circle coordinate reverse the direction and angle
	 */

	double roll,pitch,yaw;
	roll = Coord.roll*PI/180.0;
	pitch = Coord.pitch*PI/180;
	yaw = Coord.yaw*PI/180.0;
	
	for (int i=0; i< numPoint; i++){
	 
	  /* 
	   * Notes:
	   * Rotate it to align with [nx, ny, nz]
	   * rotate along x axis 
	   * y' = cos()*y + sin()*z;
	   * z' = -sin()*y + cos()*z;
	   * rotate along y axis
	   * x' = cos()*x - sin()*z;
	   * z' = sin()*x + cos()*z;
	   * rotate along z axis
	   * x' = cos()*x + sin()*y;
	   * y' = -sin()*x + cos()*y;
	   */
	   
	  /*
	   * Note: 
	   * from circle to  world (forward transform)
	   * a. rotate along z with roll degree
	   * b. rotate along y with pitxh degree
	   * c. rotate along x with yaw degree
	   * From world to circle (inverse trnasform)
	   * a. rotate along x with -yaw degree
	   * b. rotate along y with -pitch degree
	   * c. rotate along z with -roll degree
	   * Here roll pitch, and yaw have been revesed already
	   */


      temp_x = Point3D[i][0];
	  temp_y = Point3D[i][1];
	  temp_z = Point3D[i][2];


	  // first rotate along z axis with roll degree
	  x1 = cos(roll)*temp_x + sin(roll)*temp_y;
	  y1 = -sin(roll)*temp_x + cos(roll)*temp_y;
      z1 = temp_z;
      
	  // rotate along y axis with pitch degree
	  x2 = cos(pitch)*x1 - sin(pitch)*z1;
	  z2 = sin(pitch)*x1 + cos(pitch)*z1;
	  y2 = y1;

      
	  // last rotate along x axis w yaw degree
	  y = cos(yaw)*y2+ sin(yaw)*z2;
	  z = -sin(yaw)*y2 + cos(yaw)*z2;
	  x = x2;
	
      // rotated circle points will be
	  // z should be close to constatnt
	  Point3D[i][0] = x;
	  Point3D[i][1] = y;
	  Point3D[i][2] = z;
	
	  Point2D[i][0] = Point3D[i][0];
	  Point2D[i][1] = Point3D[i][1];
   
	}
}

void CStdMath::Transfer2Dto3D(int numPoint, double *Point2D[], double *Point3D[])
{
	double temp_x,temp_y,temp_z,x1,y1,z1,x2,y2,z2,x,y,z;

	// revese the rotation direction 
	
	roll  = -roll;
	pitch = -pitch;
	yaw   = -yaw;
     

	for(int i =0; i< numPoint; i++)
	{
	/************** 
	  * rotate the circle center back to reflect the result in the world coordinate
	  ***************/

      temp_x = Point2D[i][0];
	  temp_y = Point2D[i][1];
	  temp_z = Point3D[i][2];   // all z should be very close   -- better to take average   

	  // first rotate along x axis (yaw)
	  y1 = cos(yaw)*temp_y + sin(yaw)*temp_z;
	  z1 = -sin(yaw)*temp_y + cos(yaw)*temp_z;
	  x1 = temp_x;

	  // rotate along y axis  (pitch)
	  x2 = cos(pitch)*x1- sin(pitch)*z1;
	  z2 = sin(pitch)*x1 + cos(pitch)*z1;
	  y2 = y1;

	  // rotate along z axis (roll)
      x = cos(roll)*x2 + sin(roll)*y2;
	  y = -sin(roll)*x2 + cos(roll)*y2;
      z = z2;
	  	  	  
      // final result after rotaion back
	  Point3D[i][0] = x;
	  Point3D[i][1] = y;
	  Point3D[i][2] = z;
	}
}


// Rotate the input points Point3D[i][3]
// into the plane that is in parellel with xy plane (Z value are the same)
// input: points in 3D space (they may lay on a plane): Point3D[]
// output: points after transform: *Point3D[] in which all Z are the same
//         its rotation angles are *rpy[3]
void CStdMath::Rotate3DinZ(int numPoint, double *Point3D[], double rpy[3])
{
	double temp_x,temp_y,temp_z,x1,y1,z1,x2,y2,z2,x,y,z;
 	double nx, ny, nz,dis;
    double real_nx, real_ny, real_nz;
	double p[3], d[3], n[3];
    PT_RPY	Coord;

  
	GetPlaneNormal(numPoint, Point3D, &nx, &ny, &nz, &dis);
    // normalize the orientation vector in case they are not normalized
	real_nx = nx/sqrt(nx*nx + ny*ny+ nz*nz);
	real_ny = ny/sqrt(nx*nx + ny*ny+ nz*nz);
	real_nz = nz/sqrt(nx*nx + ny*ny+ nz*nz);


	/* 
	 * transform the input points into the plane that is in parellel with xy plane
	 */
	// assume that this vector pass origin
	// since we only consider the orientation
	p[0] = 0.0;
	p[1] = 0.0;
	p[2] = 0.0;

	n[0] = real_nx;
	n[1] = real_ny;
	n[2] = real_nz;

	/* 
	 * Need determine (give) direction vector to create RPY or coordinate
	 * since ijy info is not enogh to determine coordinate
	 */
	 
	// Assume direction vector d[] is normal to n[] or (i,j,k)
	if (real_nx != 0.0){
	  d[0] = -real_ny*0.5/real_nx;
	  d[1] = 0.5;
	  d[2] = 0.0;
	} else if (real_ny != 0.0){
        d[0] = 0.5;
		d[1] = -real_nx*0.5/real_ny;
	    d[2] = 0.0;
    } else { // direction is x 
		d[0] = 1.0;
	    d[1] = 0.0;
		d[2] = 0.0;
	}


	// Define a coordinate for circle based on above normal and direction vector
	Coord = PTIJK_to_PTRPY(p, d, n);


    /*
	 * roll, pitch & yaw indicate how the designed circle coordinate transfers 
	 *into world coordinate. 
	 * From world to circle coordinate reverse the direction and angle
	 */

	rpy[0] = Coord.roll*PI/180.0;
	rpy[1] = Coord.pitch*PI/180;
	rpy[2] = Coord.yaw*PI/180.0;
	
	
	for (int i=0; i< numPoint; i++){
	  

      temp_x = Point3D[i][0];
	  temp_y = Point3D[i][1];
	  temp_z = Point3D[i][2];


	  // first rotate along z axis with roll degree

	  x1 = cos(roll)*temp_x + sin(roll)*temp_y;
	  y1 = -sin(roll)*temp_x + cos(roll)*temp_y;
      z1 = temp_z;
      
	  // rotate along y axis with pitch degree
	  x2 = cos(pitch)*x1 - sin(pitch)*z1;
	  z2 = sin(pitch)*x1 + cos(pitch)*z1;
	  y2 = y1;

      
	  // last rotate along x axis w yaw degree
	  y = cos(yaw)*y2+ sin(yaw)*z2;
	  z = -sin(yaw)*y2 + cos(yaw)*z2;
	  x = x2;
	
      // rotated circle points will be
	  // z should be close to constatnt
	  Point3D[i][0] = x;
	  Point3D[i][1] = y;
	  Point3D[i][2] = z;
   
	}
}

// Rotate 3D points Point3D[i][3]
// with angles rpy[3]
// input/out: Point3D[i][3]
void CStdMath::Rotate3D(int numPoint, double *Point3D[], double rpy[3])
{
	double temp_x,temp_y,temp_z,x1,y1,z1,x2,y2,z2,x,y,z;
    double roll, pitch, yaw;

	// revese the rotation direction 
	roll  = -rpy[0];
	pitch = -rpy[1];
	yaw   = -rpy[2];
     

	for(int i =0; i< numPoint; i++)
	{
	/************** 
	  * rotate the circle center back to reflect the result in the world coordinate
	  ***************/

      temp_x = Point3D[i][0];
	  temp_y = Point3D[i][1];
	  temp_z = Point3D[i][2];   // all z should be very close   -- better to take average   

	  // first rotate along x axis (yaw)
	  y1 = cos(yaw)*temp_y + sin(yaw)*temp_z;
	  z1 = -sin(yaw)*temp_y + cos(yaw)*temp_z;
	  x1 = temp_x;

	  // rotate along y axis  (pitch)
	  x2 = cos(pitch)*x1- sin(pitch)*z1;
	  z2 = sin(pitch)*x1 + cos(pitch)*z1;
	  y2 = y1;

	  // rotate along z axis (roll)
      x = cos(roll)*x2 + sin(roll)*y2;
	  y = -sin(roll)*x2 + cos(roll)*y2;
      z = z2;
	  	  	  
      // final result after rotaion back
	  Point3D[i][0] = x;
	  Point3D[i][1] = y;
	  Point3D[i][2] = z;
	}
}




// EllipseFitting3D() will fit a group of(numPoint) spacial points into an 
// ellipse, it then replace these points with fittedPoint spacial points. 
// To avoid costly transformation,
//  first:  the plane where the spacial points fits is calculated;
//  second: all points are projected to XY plane to get 2D ellipse points;
//  third:  2D points are then fitted into an ellipse using LSR;
//  forth:  2D points are then projected to the plane calculated in step 1 to 
//          generate 3D ellipse points. 
//  Modification: XY plane may not be a good option if the plane calculated
//                from step one is perpendicular to XY plane  

void CStdMath::EllipseFitting3D(int     numPoint, 
							   double *ellipse[3], 
							   int     fittedPoint,
							   double *newEllipse[3])
{
  double  **ellipse2D;
  double  **newEllipse2D;
  int     i;
  int     planeId =0;

  ellipse2D = new double*[numPoint];
  for(int i =0; i<numPoint; i++)
  {
	ellipse2D[i] = new double [2];
  }

  newEllipse2D = new double*[fittedPoint];
  for(int i=0; i< fittedPoint; i++)
  {
	newEllipse2D[i] = new double[2];
  }
 
  // first:  the plane where the spacial points fits is calculated;
  double nx,ny,nz,d;
  GetPlaneNormal(numPoint, ellipse, &nx, &ny, &nz, &d);
  

  // determine the projected plane
  
  if(nz >= nx && nz >= ny)
	 planeId = 0;              // XY plane is the projected plane

  if(ny > nx && ny > nz)
	 planeId = 1;              // XZ plane is the projected plane

  if(nx > ny && nx > nz)
	 planeId = 2;              // YZ plane is the projected plane


  // second: all points are projected to the project plane to get 
  //         2D ellipse points;
 switch(planeId)
 {
   case 0:
	  for(int i=0; i< numPoint; i++)
	  {
		ellipse2D[i][0] = ellipse[i][0];
		ellipse2D[i][1] = ellipse[i][1];
	  }
	  break;
  case 1:
	  for(int i=0; i< numPoint; i++)
	  {
		ellipse2D[i][0] = ellipse[i][0];
		ellipse2D[i][1] = ellipse[i][2];
	  }
	  break;
  case 2:
	  for(int i=0; i< numPoint; i++)
	  {
		ellipse2D[i][0] = ellipse[i][1];
		ellipse2D[i][1] = ellipse[i][2];
	  }
	  break;
 }
  // third, fit the ellipse2D points, and generate fittedPoint number
  // of ellipse2D points

  EllipseFitting_2D(numPoint, ellipse2D, fittedPoint, newEllipse2D);

  // forth, transfer ellipse2D point into ellipse3D point

     switch(planeId)
 {
   case 0:
	  for(int i =0; i<fittedPoint; i++)
	  {
		newEllipse[i][0] = newEllipse2D[i][0];
		newEllipse[i][1] = newEllipse2D[i][1];
		newEllipse[i][2] = (d-nx*newEllipse[i][0]-ny*newEllipse[i][1])/nz;
	  }
	  break;
   case 1:
	  for(int i =0; i<fittedPoint; i++)
	  {
		newEllipse[i][0] = newEllipse2D[i][0];
		newEllipse[i][2] = newEllipse2D[i][1];
		newEllipse[i][1] = (d-nx*newEllipse[i][0]-nz*newEllipse[i][2])/ny;
	  }
	  break;
   case 2:
	  for(int i =0; i<fittedPoint; i++)
	  {
		newEllipse[i][1] = newEllipse2D[i][0];
		newEllipse[i][2] = newEllipse2D[i][1];
		newEllipse[i][0] = (d-nz*newEllipse[i][2]-ny*newEllipse[i][1])/nx;
	  }
	  break;
 }
  // finally, free memory of ellipse2D

  for(int i=0; i<numPoint; i++)
  {
	delete[] ellipse2D[i];
  }
  delete[] ellipse2D;

  for(int i=0; i<fittedPoint; i++){
  	delete[]newEllipse2D[i];
  }
  delete[]newEllipse2D;

}

// Similar with above function
// Output:	parameter of ellipse 
//           center position (cx, cy, cz), 
//           long axis: xAxi
//			 short axis: yAxi
//           orientation angle: Theta

void CStdMath::EllipseFitting3D(int     numPoint, 
							   double *ellipse[3], 
							   int     fittedPoint,
							   double *newEllipse[3],
							   double *cx,
							   double *cy,
							   double *cz,
							   double *xAxi, 
							   double *yAxi, 
							   double *Theta)
{
  double  **ellipse2D;
  double  **newEllipse2D;
  int     i;
  int     planeId =0;

  ellipse2D = new double*[numPoint];
  for(int i =0; i<numPoint; i++)
  {
	ellipse2D[i] = new double [2];
  }

  newEllipse2D = new double*[fittedPoint];
  for(int i=0; i< fittedPoint; i++)
  {
	newEllipse2D[i] = new double[2];
  }
  // first:  the plane where the spacial points fits is calculated;

  GetPlaneNormal(numPoint, ellipse, &nx, &ny, &nz, &d);
  

  // determine the projected plane
  
  if(nz >= nx && nz >= ny)
	 planeId = 0;              // XY plane is the projected plane

  if(ny > nx && ny > nz)
	 planeId = 1;              // XZ plane is the projected plane

  if(nx > ny && nx > nz)
	 planeId = 2;              // YZ plane is the projected plane


  // second: all points are projected to the project plane to get 
  //         2D ellipse points;
 switch(planeId)
 {
   case 0:
	  for(int i=0; i< numPoint; i++)
	  {
		ellipse2D[i][0] = ellipse[i][0];
		ellipse2D[i][1] = ellipse[i][1];
	  }
	  break;
  case 1:
	  for(int i=0; i< numPoint; i++)
	  {
		ellipse2D[i][0] = ellipse[i][0];
		ellipse2D[i][1] = ellipse[i][2];
	  }
	  break;
  case 2:
	  for(int i=0; i< numPoint; i++)
	  {
		ellipse2D[i][0] = ellipse[i][1];
		ellipse2D[i][1] = ellipse[i][2];
	  }
	  break;
 }
  // third, fit the ellipse2D points, and generate fittedPoint number
  // of ellipse2D points

  double xc2D, yc2D;
  EllipseFitting_2D(numPoint, ellipse2D, fittedPoint, newEllipse2D, &xc2D, &yc2D, xAxi, yAxi, Theta);

  // forth, transfer ellipse2D point into ellipse3D point

     switch(planeId)
 {
   case 0:
		*cx = xc2D;
		*cy = yc2D;
		*cz = (d - nx*(*cx) - ny*(*cy)) / nz;
	  for(int i =0; i<fittedPoint; i++)
	  {
		newEllipse[i][0] = newEllipse2D[i][0];
		newEllipse[i][1] = newEllipse2D[i][1];
		newEllipse[i][2] = (d-nx*newEllipse[i][0]-ny*newEllipse[i][1])/nz;
	  }
	  break;
   case 1:
	   *cx = xc2D;
	   *cz = yc2D;
	   *cy = (d - nx*(*cx) - nz*(*cz)) / ny;
	  for(int i =0; i<fittedPoint; i++)
	  {
		newEllipse[i][0] = newEllipse2D[i][0];
		newEllipse[i][2] = newEllipse2D[i][1];
		newEllipse[i][1] = (d-nx*newEllipse[i][0]-nz*newEllipse[i][2])/ny;
	  }
	  break;
   case 2:
	   *cy = xc2D;
	   *cz = yc2D;
	   *cx = (d - ny*(*cy) - nz*(*cz)) / nx;
	  for(int i =0; i<fittedPoint; i++)
	  {
		newEllipse[i][1] = newEllipse2D[i][0];
		newEllipse[i][2] = newEllipse2D[i][1];
		newEllipse[i][0] = (d-nz*newEllipse[i][2]-ny*newEllipse[i][1])/nx;
	  }
	  break;
 }
  // finally, free memory of ellipse2D

  for(int i=0; i<numPoint; i++)
  {
	delete[] ellipse2D[i];
  }
  delete[] ellipse2D;
  for(int i=0; i<fittedPoint; i++){
  	delete[]newEllipse2D[i];
  }
  delete[]newEllipse2D;
}



bool CStdMath::Normalize(double &nx,double &ny, double &nz)
{
  double temp;

  temp = sqrt(nx*nx + ny*ny + nz*nz);

  nx = nx/temp;
  ny = ny/temp;
  nz = nz/temp;

  return true;
}



/*******************************************************************************************
 * GetXShiftParams() calculates the X-shift table cofficients (k1,k2,k3)  
 * by using nonlinear least sqaure method
 * (k1,k2,k3) is deteremined such that when motor moves L distance
 * the laser imaging system has (k1*L, k2*L, k3*L) translation relative to world coodinate
 *
 * Inputs: 
 *        p[][5]; containing a measured point array of [x,y,z] on the a fixed table coodinate 
 *                and motor move distance L and sphere radius R
 *               p[][5] = {x,y,z,L,R} where R is constant
 *        num_point; number of pints of the input; 
 *        int_par[3]; intial paramters (init_k1, init_k2, init_k3);
 * Outputs:  
 *         (k1,k2,k3)
 *         standard devaiation (stdDev) and maxmium deviation (maxDev) of fitting 
 *                    
 * return:
 *         void
 * 
 * Notes:
 * 	 x[1] = xi, x[2] = yi, x[3] = zi; x[4] = Li; x[5] = R;
 *	 a[1] = k1, a[2] = k2, a[3] = k3;
 *		// The fitting constraint is
 *		// (xi+k1*Li)^2 + (yi+k2*Li)^2 + (zi+k3*Li)^2 = R^2;
 *		// where (xi,yi,zi) are laser point position in table when x-shift is not considered
 *		// Li is the motor shift distance (encoder reading)
 *		// when x-shift moves Li distance the laser imaging system has translation of 
 *		// (k1*Li, k2*Li, k3*Li) 
 *		// (ki, k2,k3) need to be calibrated
 * 
 * 
 *          12.27.01 
 */

void CStdMath::GetXShiftParams(int num_point, double  *p[5], double init_par[3], double *k1, double *k2, double *k3, double *stdDev, double *maxDev)
{
  
	// declaration & initialization
	int ndata = num_point;		// number of input points x[i][k]; i = 1,... ndata
	int nd = 5;					// nd: number of given x[i][k],  k = 1,2 ,3;
	int ma = 3;					// number of parameters to be determined after fitting
	int mfit = 3;				// = ma

	int* lista;					// lista[1..ma]
	double ** x;					// x[1..ndata][1..nd]	
	double * y;					// y[1..ndata]
	double * sig;					// sig[1..ndata]
	double * a;					// a[1..ma]


	// memory allocation
	lista = new int[ma+1];
	 
	x = new double*[ndata+1];
	for(int i = 1; i <= ndata; i++) {
		x[i] = new double[nd+1];                                                   		 
	}

	y = new double[ndata+1];
	sig = new double[ndata+1];
	a = new double[ma+1];
  
	for(int k = 1; k <= ma; k++) {
		lista[k] = k;			// numbers the parameters a such that 
 						// the first mfit elements correspond to values actually being adjusted
	}
  
	for(int k = 1; k <= ndata; k++) {
		sig[k] = 1;				// standard deviations, this is not experiment data, so set to 1
	    y[k] = 0.0;             // y = 0
	}

    // assign x value
    for(int i = 1; i <= ndata; i++) {
		for (int j = 1; j <= nd; j++){
			x[i][j] = (double) p[i-1][j-1];
		}
    	
	}


   // initial value of parameters
	a[1] = init_par[0];
	a[2] = init_par[1];
	a[3] = init_par[2];
 
 	NonLinearLeastSquare(x, y, sig, ndata, a, ma, lista, mfit, 7);
    // HOW TO USE NLLS:
	// (x,y): input point pair y= f(ax); normally y = 0;
	// sig: weighting factors. normally sig = 1;
	// ndata: number of data point.
	// a: coefficients to be determined -- outcome
	// ma: number of a
	// mfit: number of a to be fitted
	// lista: array of a [1 ...ma];

    
	// calculate the fitting error 
	double max_delta = -10000.0;
    double dd, delta, delta_sq, delta_sum;
	
	delta_sum = 0.0;

   
    for (int i = 0; i< num_point; i++){
	    
            // fitting error to sphare center
			dd = sqrt((p[i][0]+a[1]*p[i][3])*(p[i][0]+a[1]*p[i][3])+(p[i][1]+a[2]*p[i][3])*(p[i][1]+a[2]*p[i][3]) + (p[i][2]+a[3]*p[i][3])*(p[i][2]+a[3]*p[i][3])); 
			delta = fabs(dd - p[i][4]);
			delta_sq = (dd - p[i][4])*(dd - p[i][4]);
            delta_sum = delta_sum + delta_sq;
			if (max_delta < delta) max_delta = delta;


	}
	  *stdDev = sqrt(delta_sum/(double)num_point);
      *maxDev = max_delta;
   

	  // normalize the vector
       *k1 = a[1];
	   *k2 = a[2];
	   *k3 = a[3];

    delete lista;
	for(int i = 1; i <= ndata; i++) {
		delete x[i];                                                   		 
	}
    delete x;
    delete y;
    delete sig;
    delete a;

}

// The same as above except for inclusion of sphere center (x0.y0,z0)
// 
/*******************************************************************************************
 * GetXShiftParams() calculates the X-shift table cofficients (k1,k2,k3)  
 * by using nonlinear least sqaure method
 * (k1,k2,k3) is deteremined such that when motor moves L distance
 * the laser imaging system has (k1*L, k2*L, k3*L) translation relative to world coodinate
 *
 * Inputs: 
 *        p[][5]; containing a measured point array of [x,y,z] on the a fixed table coodinate 
 *                and motor move distance L and sphere radius R
 *               p[][5] = {x,y,z,L,R} where R is constant
 *        num_point; number of pints of the input; 
 *        int_par[6]; intial paramters (init_k1, init_k2, init_k3, init_x0, init_y0, init_z0);
 * Outputs:  
 *         (k1,k2,k3, x0, y0, z0)
 *         standard devaiation (stdDev) and maxmium deviation (maxDev) of fitting 
 *                    
 * return:
 *         void
 * 
 * Notes:
 * 	 x[1] = xi, x[2] = yi, x[3] = zi; x[4] = Li; x[5] = R;
 *	 a[1] = k1, a[2] = k2, a[3] = k3; a[4] = x0, a[5] = y0; a[6] = z0;
 *		// The fitting constraint is
 *		// (xi+k1*Li - x0)^2 + (yi+k2*Li - y0)^2 + (zi+k3*Li - z0)^2 = R^2;
 *		// where (xi,yi,zi) are laser point position in table when x-shift is not considered
 *		// Li is the motor shift distance (encoder reading)
 *		// when x-shift moves Li distance the laser imaging system has translation of 
 *		// (k1*Li, k2*Li, k3*Li) 
 *		// (ki, k2,k3) need to be calibrated
 * 
 * 
 *          12.27.01 
 */


void CStdMath::GetXShiftParams(int num_point, double  *p[5], double init_par[6], double *k1, double *k2, double *k3, double *x0, double *y0, double *z0, double *stdDev, double *maxDev)
{
  
	// declaration & initialization
	int ndata = num_point;		// number of input points x[i][k]; i = 1,... ndata
	int nd = 5;					// nd: number of given x[i][k],  k = 1,2 ,3;
	int ma = 6;					// number of parameters to be determined after fitting
	int mfit = 6;				// = ma

	int* lista;					// lista[1..ma]
	double ** x;					// x[1..ndata][1..nd]	
	double * y;					// y[1..ndata]
	double * sig;					// sig[1..ndata]
	double * a;					// a[1..ma]


	// memory allocation
	lista = new int[ma+1];
	 
	x = new double*[ndata+1];
	for(int i = 1; i <= ndata; i++) {
		x[i] = new double[nd+1];                                                   		 
	}

	y = new double[ndata+1];
	sig = new double[ndata+1];
	a = new double[ma+1];
  
	for(int k = 1; k <= ma; k++) {
		lista[k] = k;			// numbers the parameters a such that 
 						// the first mfit elements correspond to values actually being adjusted
	}
  
	for(int k = 1; k <= ndata; k++) {
		sig[k] = 1;				// standard deviations, this is not experiment data, so set to 1
	    y[k] = 0.0;             // y = 0
	}

    // assign x value
    for(int i = 1; i <= ndata; i++) {
		for (int j = 1; j <= nd; j++){
			x[i][j] = (double) p[i-1][j-1];
		}
    	
	}


   // initial value of parameters
	a[1] = init_par[0];
	a[2] = init_par[1];
	a[3] = init_par[2];
	a[4] = init_par[3];
	a[5] = init_par[4];
	a[6] = init_par[5];
 
 	NonLinearLeastSquare(x, y, sig, ndata, a, ma, lista, mfit, 8);
    // HOW TO USE NLLS:
	// (x,y): input point pair y= f(ax); normally y = 0;
	// sig: weighting factors. normally sig = 1;
	// ndata: number of data point.
	// a: coefficients to be determined -- outcome
	// ma: number of a
	// mfit: number of a to be fitted
	// lista: array of a [1 ...ma];

    
	// calculate the fitting error 
	double max_delta = -10000.0;
    double dd, delta, delta_sq, delta_sum;
	
	delta_sum = 0.0;

   
    for (int i = 0; i< num_point; i++){
	    
            // fitting error to sphare center
			dd = sqrt((p[i][0]+a[1]*p[i][3]-a[4])*(p[i][0]+a[1]*p[i][3]-a[4])+(p[i][1]+a[2]*p[i][3]-a[5])*(p[i][1]+a[2]*p[i][3]-a[5]) + (p[i][2]+a[3]*p[i][3]-a[6])*(p[i][2]+a[3]*p[i][3]-a[6])); 
			delta = fabs(dd - p[i][4]);
			delta_sq = (dd - p[i][4])*(dd - p[i][4]);
            delta_sum = delta_sum + delta_sq;
			if (max_delta < delta) max_delta = delta;


	}
	  *stdDev = sqrt(delta_sum/(double)num_point);
      *maxDev = max_delta;
   

	  // normalize the vector
       *k1 = a[1];
	   *k2 = a[2];
	   *k3 = a[3];

	   *x0 = a[4];
	   *y0 = a[5];
	   *z0 = a[6];

    delete lista;
	for(int i = 1; i <= ndata; i++) {
		delete x[i];                                                   		 
	}
    delete x;
    delete y;
    delete sig;
    delete a;

}


// The same as above except that (x0, y0,z0) is not adjusted
// 
/*******************************************************************************************
 * GetXShiftParams() calculates the X-shift table cofficients (k1,k2,k3)  
 * by using nonlinear least sqaure method
 * (k1,k2,k3) is deteremined such that when motor moves L distance
 * the laser imaging system has (k1*L, k2*L, k3*L) translation relative to world coodinate
 *
 * Inputs: 
 *        p[][5]; containing a measured point array of [x,y,z] on the a fixed table coodinate 
 *                and motor move distance L and sphere radius R
 *               p[][5] = {x,y,z,L,R} where R is constant
 *        num_point; number of pints of the input; 
 *        int_par[6]; intial paramters (init_k1, init_k2, init_k3, init_x0, init_y0, init_z0);
 * Outputs:  
 *         (k1,k2,k3, x0, y0, z0)
 *         standard devaiation (stdDev) and maxmium deviation (maxDev) of fitting 
 *                    
 * return:
 *         void
 * 
 * Notes:
 * 	 x[1] = xi, x[2] = yi, x[3] = zi; x[4] = Li; x[5] = R;
 *	 a[1] = k1, a[2] = k2, a[3] = k3; a[4] = x0, a[5] = y0; a[6] = z0;
 *		// The fitting constraint is
 *		// (xi+k1*Li - x0)^2 + (yi+k2*Li - y0)^2 + (zi+k3*Li - z0)^2 = R^2;
 *		// where (xi,yi,zi) are laser point position in table when x-shift is not considered
 *		// Li is the motor shift distance (encoder reading)
 *		// when x-shift moves Li distance the laser imaging system has translation of 
 *	
 *    // init: {k1, k2, k3, x0,y0, z0}
 *    // {x0, y0, z0} is not adjusted 
 *    // {k1, k2, k3} will be adjusted
 * 
 *          1.30.02 
 */


void CStdMath::GetXShiftParams2(int num_point, double  *p[5], double init_par[6], double *k1, double *k2, double *k3, double *x0, double *y0, double *z0, double *stdDev, double *maxDev)
{
  
	// declaration & initialization
	int ndata = num_point;		// number of input points x[i][k]; i = 1,... ndata
	int nd = 5;					// nd: number of given x[i][k],  k = 1,2 ,3;
	int ma = 6;					// number of parameters to be determined after fitting
	int mfit = 3;				// = ma

	int* lista;					// lista[1..ma]
	double ** x;					// x[1..ndata][1..nd]	
	double * y;					// y[1..ndata]
	double * sig;					// sig[1..ndata]
	double * a;					// a[1..ma]


	// memory allocation
	lista = new int[ma+1];
	 
	x = new double*[ndata+1];
	for(int i = 1; i <= ndata; i++) {
		x[i] = new double[nd+1];                                                   		 
	}

	y = new double[ndata+1];
	sig = new double[ndata+1];
	a = new double[ma+1];
  
	for(int k = 1; k <= ma; k++) {
		lista[k] = k;			// numbers the parameters a such that 
 						// the first mfit elements correspond to values actually being adjusted
	}
  
	for(int k = 1; k <= ndata; k++) {
		sig[k] = 1;				// standard deviations, this is not experiment data, so set to 1
	    y[k] = 0.0;             // y = 0
	}

    // assign x value
    for(int i = 1; i <= ndata; i++) {
		for (int j = 1; j <= nd; j++){
			x[i][j] = (double) p[i-1][j-1];
		}
    	
	}


   // initial value of parameters
	a[1] = init_par[0];
	a[2] = init_par[1];
	a[3] = init_par[2];
	a[4] = init_par[3];
	a[5] = init_par[4];
	a[6] = init_par[5];
 
 	NonLinearLeastSquare(x, y, sig, ndata, a, ma, lista, mfit, 8);
    // HOW TO USE NLLS:
	// (x,y): input point pair y= f(ax); normally y = 0;
	// sig: weighting factors. normally sig = 1;
	// ndata: number of data point.
	// a: coefficients to be determined -- outcome
	// ma: number of a
	// mfit: number of a to be fitted
	// lista: array of a [1 ...ma];

    
	// calculate the fitting error 
	double max_delta = -10000.0;
    double dd, delta, delta_sq, delta_sum;
	
	delta_sum = 0.0;

   
    for (int i = 0; i< num_point; i++){
	    
            // fitting error to sphare center
			dd = sqrt((p[i][0]+a[1]*p[i][3]-a[4])*(p[i][0]+a[1]*p[i][3]-a[4])+(p[i][1]+a[2]*p[i][3]-a[5])*(p[i][1]+a[2]*p[i][3]-a[5]) + (p[i][2]+a[3]*p[i][3]-a[6])*(p[i][2]+a[3]*p[i][3]-a[6])); 
			delta = fabs(dd - p[i][4]);
			delta_sq = (dd - p[i][4])*(dd - p[i][4]);
            delta_sum = delta_sum + delta_sq;
			if (max_delta < delta) max_delta = delta;


	}
	  *stdDev = sqrt(delta_sum/(double)num_point);
      *maxDev = max_delta;
   

	  // normalize the vector
       *k1 = a[1];
	   *k2 = a[2];
	   *k3 = a[3];

	   *x0 = a[4];
	   *y0 = a[5];
	   *z0 = a[6];

    delete lista;
	for(int i = 1; i <= ndata; i++) {
		delete x[i];                                                   		 
	}
    delete x;
    delete y;
    delete sig;
    delete a;

}


/*******************************************************************************************
 * GetXShiftAndRotateParams() calculates the X-shift table cofficients (k1,k2,k3) and rotation
 * center (rot_x, 0, rot_z);
 * by using nonlinear least sqaure method
 * (k1,k2,k3) is deteremined such that when motor moves L distance
 * the laser imaging system has (k1*L, k2*L, k3*L) translation relative to world coodinate
 * 
 * Inputs: 
 *        p[][6]; containing a measured point array of [x,y,z] on the vision coodinate (a fixed table) 
 *                and rotation angle theta, linear motor move distance L and sphere radius R
 *               p[][5] = {x,y,z,theta,L,R} where R is constant
 *        num_point; number of pints of the input; 
 *        int_par[8]; intial paramters (init_k1, init_k2, init_k3, x0, y0, z0, rot_x, rot_z);
 *                    where (x0,y0,z0) is the position of calibration sphere
 * Outputs:  
 *         (k1,k2,k3, x0,y0, z0, rot_x, rot_z)
 *         standard devaiation (stdDev) and maxmium deviation (maxDev) of fitting 
 *                    
 * return:
 *         void
 * 
 * Notes:
 *		x[1] = xi, x[2] = yi, x[3] = zi; x[4] = thetai, x[5] = Li; x[6] = R;
 *		a[1] = k1, a[2] = k2, a[3] = k3; a[4] = x0, a[5]=y0, a[6] = z0, a[7] = rot_x, a[8] = rot_z;
 *		The fitting constraint is
 *		x' = cos(thetai)(xi+k1*Li - rot_x) - sin(thetai)(zi + k3*Li-rot_z) + rot_x;
 *		y' = yi + k2*Li;
 *		z' = sin(thetai)(xi+k1*Li - rot_x) + cos(thetai)(zi + k3*Li-rot_z) + rot_z;
 *
 *		(x' -x0)^2 + (y'-y0)^2 + (z'-z0)^2 = R^2;                          (1)
 *                                                
 *		where (xi,yi,zi) are laser point position in table when x-shift is not considered
 *		Li is the motor shift distance (encoder reading)
 *		thetai is the rotation angle
 *		when x-shift moves Li distance the laser imaging system has translation of 
 *		(k1*Li, k2*Li, k3*Li) 
 * 
 * 
 *      01/09/02 
 * 
 *      Other approach to derive the equation:	
 *		Convert local coodinate (xi,yi,zi) into rotation table coodinate
 *		x' = cos(thetai)(xi+k1*Li - rot_x) - sin(thetai)(zi + k3*Li-rot_z);
 *		y' = yi + k2*Li;
 *		z' = sin(thetai)(xi+k1*Li - rot_x) + cos(thetai)(zi + k3*Li-rot_z);
 *		
 *      Convert into sphere coordinate 
 *      (assume the only different between rotation and sphere coodinates are translation:
 *      (rot_x - x0), (rot_y- y0), (rot_z-z0);
 *     
 *      x" = x' + (rot_x - x0);
 *      y" = y' + (rot_y - y0);
 *      z" = z' + (rot_z - z0);
 *     
 *  	(x")^2 + (y")^2 + (z")^2 = R^2;                              (2)
 *    
 *       Eq. (2) is the same as eq. (1)
 */



void CStdMath::GetXShiftAndRotateParams(int num_point, double  *p[6], double init_par[8], double *k1, double *k2, double *k3, double *x0, double *y0, double *z0, double *rot_x, double *rot_z, double *stdDev, double *maxDev)
{
  
	// declaration & initialization
	int ndata = num_point;		// number of input points x[i][k]; i = 1,... ndata
	int nd = 6;					// nd: number of given x[i][k],  k = 1,2 ,3;
	int ma = 8;					// number of parameters to be determined after fitting
	int mfit = 8;				// = ma

	int* lista;					// lista[1..ma]
	double ** x;					// x[1..ndata][1..nd]	
	double * y;					// y[1..ndata]
	double * sig;					// sig[1..ndata]
	double * a;					// a[1..ma]


	// memory allocation
	lista = new int[ma+1];
	 
	x = new double*[ndata+1];
	for(int i = 1; i <= ndata; i++) {
		x[i] = new double[nd+1];                                                   		 
	}

	y = new double[ndata+1];
	sig = new double[ndata+1];
	a = new double[ma+1];
  
	for(int k = 1; k <= ma; k++) {
		lista[k] = k;			// numbers the parameters a such that 
 						// the first mfit elements correspond to values actually being adjusted
	}
  
	for(int k = 1; k <= ndata; k++) {
		sig[k] = 1;				// standard deviations, this is not experiment data, so set to 1
	    y[k] = 0.0;             // y = 0
	}

    // assign x value
    for(int i = 1; i <= ndata; i++) {
		for (int j = 1; j <= nd; j++){
			x[i][j] = (double) p[i-1][j-1];
		}
    	
	}


   // initial value of parameters
	a[1] = init_par[0];
	a[2] = init_par[1];
	a[3] = init_par[2];
	a[4] = init_par[3];
	a[5] = init_par[4];
	a[6] = init_par[5];
	a[7] = init_par[6];
	a[8] = init_par[7];
 
 	NonLinearLeastSquare(x, y, sig, ndata, a, ma, lista, mfit, 9);
    // HOW TO USE NLLS:
	// (x,y): input point pair y= f(ax); normally y = 0;
	// sig: weighting factors. normally sig = 1;
	// ndata: number of data point.
	// a: coefficients to be determined -- outcome
	// ma: number of a
	// mfit: number of a to be fitted
	// lista: array of a [1 ...ma];

    
     
       *k1 = a[1];
	   *k2 = a[2];
	   *k3 = a[3];

	   *x0 = a[4];
	   *y0 = a[5];
	   *z0 = a[6];

	   *rot_x = a[7];
	   *rot_z = a[8];


	// calculate the fitting error 
	double max_delta = -10000.0;
    double dd, delta, delta_sq, delta_sum;
	double x_prime, y_prime, z_prime;

	delta_sum = 0.0;


   
    for (int i = 0; i< num_point; i++){
	    
            // fitting error to sphare center
		    double xi = p[i][0];
			double yi = p[i][1];
			double zi = p[i][2];
			double thetai = p[i][3];
			double Li = p[i][4];
			double R = p[i][5];
			double kk1 = a[1];
			double kk2 = a[2];
			double kk3 = a[3];
			double xx0 = a[4];
			double yy0 = a[5];
			double zz0 = a[6];
			double rotx = a[7];
			double rotz = a[8];
		     
		    
		    x_prime = cos(thetai)*(xi + kk1*Li-rotx)-sin(thetai)*(zi + kk3*Li - rotz) + rotx;
			y_prime = yi + kk2*Li;
			z_prime = sin(thetai)*(xi + kk1*Li-rotx)+cos(thetai)*(zi + kk3*Li - rotz) + rotz;
		    
			dd = sqrt((x_prime-xx0)*(x_prime-xx0)+(y_prime - yy0)*(y_prime - yy0) + (z_prime - zz0)*(z_prime - zz0)); 
			delta = fabs(dd - R);
			delta_sq = (dd - R)*(dd - R);
            delta_sum = delta_sum + delta_sq;
			if (max_delta < delta) max_delta = delta;


	}
	  *stdDev = sqrt(delta_sum/(double)num_point);
      *maxDev = max_delta;
    delete lista;
	for(int i = 1; i <= ndata; i++) {
		delete x[i];                                                   		 
	}
    delete x;
    delete y;
    delete sig;
    delete a;
}

/*    
 *      The same as GetXShiftAndRotateParams()
 *      except for that 
 *      a[1] = k1, a[2] = k2, a[3] = k3; a[4] = x0, a[5]=y0, a[6] = z0, a[7] = rot_x, a[8] = rot_z;
 *      {k1, k2, k3) are fixed the inital value will not be changed after fitting
 *      { x0, y0, z0, rot_x, rot_z } will be adjusted 
 * 
 *      01/15/02 
 */

void CStdMath::GetXShiftAndRotateParams2(int num_point, double  *p[6], double init_par[8], double *k1, double *k2, double *k3, double *x0, double *y0, double *z0, double *rot_x, double *rot_z, double *stdDev, double *maxDev)
{
  
	// declaration & initialization
	int ndata = num_point;		// number of input points x[i][k]; i = 1,... ndata
	int nd = 6;					// nd: number of given x[i][k],  k = 1,2 ,3;
	int ma = 8;					// number of parameters to be determined after fitting
	int mfit = 5;				// only adjust 5 parameters from 8 initial values

	int* lista;					// lista[1..ma]
	double ** x;					// x[1..ndata][1..nd]	
	double * y;					// y[1..ndata]
	double * sig;					// sig[1..ndata]
	double * a;					// a[1..ma]


	// memory allocation
	lista = new int[ma+1];
	 
	x = new double*[ndata+1];
	for(int i = 1; i <= ndata; i++) {
		x[i] = new double[nd+1];                                                   		 
	}

	y = new double[ndata+1];
	sig = new double[ndata+1];
	a = new double[ma+1];
  
	for(int k = 1; k <= ma; k++) {
	//	lista[k] = k;			// numbers the parameters a[1..ma] such that 
 						// the first mfit elements correspond to values actually being adjusted
		// move the first three element in a[] array to the last three 
		if (k<6){
            lista[k] = k + 3;
		} else {
			lista[k] = k - 5;
		}
	}
  
	for(int k = 1; k <= ndata; k++) {
		sig[k] = 1;				// standard deviations, this is not experiment data, so set to 1
	    y[k] = 0.0;             // y = 0
	}

    // assign x value
    for(int i = 1; i <= ndata; i++) {
		for (int j = 1; j <= nd; j++){
			x[i][j] = (double) p[i-1][j-1];
		}
    	
	}


   // initial value of parameters
	a[1] = init_par[0];
	a[2] = init_par[1];
	a[3] = init_par[2];
	a[4] = init_par[3];
	a[5] = init_par[4];
	a[6] = init_par[5];
	a[7] = init_par[6];
	a[8] = init_par[7];
 
 	NonLinearLeastSquare(x, y, sig, ndata, a, ma, lista, mfit, 9);
    // HOW TO USE NLLS:
	// (x,y): input point pair y= f(ax); normally y = 0;
	// sig: weighting factors. normally sig = 1;
	// ndata: number of data point.
	// a: coefficients to be determined -- outcome
	// ma: number of a
	// mfit: number of a to be fitted
	// lista: array of a [1 ...ma];

    
     
       *k1 = a[1];
	   *k2 = a[2];
	   *k3 = a[3];

	   *x0 = a[4];
	   *y0 = a[5];
	   *z0 = a[6];

	   *rot_x = a[7];
	   *rot_z = a[8];


	// calculate the fitting error 
	double max_delta = -10000.0;
    double dd, delta, delta_sq, delta_sum;
	double x_prime, y_prime, z_prime;

	delta_sum = 0.0;


   
    for (int i = 0; i< num_point; i++){
	    
            // fitting error to sphare center
		    double xi = p[i][0];
			double yi = p[i][1];
			double zi = p[i][2];
			double thetai = p[i][3];
			double Li = p[i][4];
			double R = p[i][5];
			double kk1 = a[1];
			double kk2 = a[2];
			double kk3 = a[3];
			double xx0 = a[4];
			double yy0 = a[5];
			double zz0 = a[6];
			double rotx = a[7];
			double rotz = a[8];
		     
		    
		    x_prime = cos(thetai)*(xi + kk1*Li-rotx)-sin(thetai)*(zi + kk3*Li - rotz) + rotx;
			y_prime = yi + kk2*Li;
			z_prime = sin(thetai)*(xi + kk1*Li-rotx)+cos(thetai)*(zi + kk3*Li - rotz) + rotz;
		    
			dd = sqrt((x_prime-xx0)*(x_prime-xx0)+(y_prime - yy0)*(y_prime - yy0) + (z_prime - zz0)*(z_prime - zz0)); 
			delta = fabs(dd - R);
			delta_sq = (dd - R)*(dd - R);
            delta_sum = delta_sum + delta_sq;
			if (max_delta < delta) max_delta = delta;


	}
	*stdDev = sqrt(delta_sum/(double)num_point);
    *maxDev = max_delta;

    delete lista;
	for(int i = 1; i <= ndata; i++) {
		delete x[i];                                                   		 
	}
    delete x;
    delete y;
    delete sig;
    delete a;
}


/*    
 *      The same as GetXShiftAndRotateParams()
 *      except for that 
 *      a[1] = k1, a[2] = k2, a[3] = k3; a[4] = x0, a[5]=y0, a[6] = z0, a[7] = rot_x, a[8] = rot_z;
 *      {x0, y0, z0, rot_x, rot_z) are fixed the inital value will not be changed after fitting
 *      --- only adjust {k1, k2, k3}
 * 
 *      01/15/02 
 */

void CStdMath::GetXShiftAndRotateParams3(int num_point, double  *p[6], double init_par[8], double *k1, double *k2, double *k3, double *x0, double *y0, double *z0, double *rot_x, double *rot_z, double *stdDev, double *maxDev)
{
  
	// declaration & initialization
	int ndata = num_point;		// number of input points x[i][k]; i = 1,... ndata
	int nd = 6;					// nd: number of given x[i][k],  k = 1,2 ,3;
	int ma = 8;					// number of parameters to be determined after fitting
	int mfit = 3;				// only adjust 3 parameters from 8 initial values

	int* lista;					// lista[1..ma]
	double ** x;					// x[1..ndata][1..nd]	
	double * y;					// y[1..ndata]
	double * sig;					// sig[1..ndata]
	double * a;					// a[1..ma]


	// memory allocation
	lista = new int[ma+1];
	 
	x = new double*[ndata+1];
	for(int i = 1; i <= ndata; i++) {
		x[i] = new double[nd+1];                                                   		 
	}

	y = new double[ndata+1];
	sig = new double[ndata+1];
	a = new double[ma+1];
  
	for(int k = 1; k <= ma; k++) {
	    lista[k] = k;	// numbers the parameters a[1..ma] such that 
 						// the first mfit elements correspond to values actually being adjusted
	
	
	}
  
	for(int k = 1; k <= ndata; k++) {
		sig[k] = 1;				// standard deviations, this is not experiment data, so set to 1
	    y[k] = 0.0;             // y = 0
	}

    // assign x value
    for(int i = 1; i <= ndata; i++) {
		for (int j = 1; j <= nd; j++){
			x[i][j] = (double) p[i-1][j-1];
		}
    	
	}


   // initial value of parameters
	a[1] = init_par[0];
	a[2] = init_par[1];
	a[3] = init_par[2];
	a[4] = init_par[3];
	a[5] = init_par[4];
	a[6] = init_par[5];
	a[7] = init_par[6];
	a[8] = init_par[7];
 
 	NonLinearLeastSquare(x, y, sig, ndata, a, ma, lista, mfit, 9);
    // HOW TO USE NLLS:
	// (x,y): input point pair y= f(ax); normally y = 0;
	// sig: weighting factors. normally sig = 1;
	// ndata: number of data point.
	// a: coefficients to be determined -- outcome
	// ma: number of a
	// mfit: number of a to be fitted
	// lista: array of a [1 ...ma];

    
     
       *k1 = a[1];
	   *k2 = a[2];
	   *k3 = a[3];

	   *x0 = a[4];
	   *y0 = a[5];
	   *z0 = a[6];

	   *rot_x = a[7];
	   *rot_z = a[8];


	// calculate the fitting error 
	double max_delta = -10000.0;
    double dd, delta, delta_sq, delta_sum;
	double x_prime, y_prime, z_prime;

	delta_sum = 0.0;


   
    for (int i = 0; i< num_point; i++){
	    
            // fitting error to sphare center
		    double xi = p[i][0];
			double yi = p[i][1];
			double zi = p[i][2];
			double thetai = p[i][3];
			double Li = p[i][4];
			double R = p[i][5];
			double kk1 = a[1];
			double kk2 = a[2];
			double kk3 = a[3];
			double xx0 = a[4];
			double yy0 = a[5];
			double zz0 = a[6];
			double rotx = a[7];
			double rotz = a[8];
		     
		    
		    x_prime = cos(thetai)*(xi + kk1*Li-rotx)-sin(thetai)*(zi + kk3*Li - rotz) + rotx;
			y_prime = yi + kk2*Li;
			z_prime = sin(thetai)*(xi + kk1*Li-rotx)+cos(thetai)*(zi + kk3*Li - rotz) + rotz;
		    
			dd = sqrt((x_prime-xx0)*(x_prime-xx0)+(y_prime - yy0)*(y_prime - yy0) + (z_prime - zz0)*(z_prime - zz0)); 
			delta = fabs(dd - R);
			delta_sq = (dd - R)*(dd - R);
            delta_sum = delta_sum + delta_sq;
			if (max_delta < delta) max_delta = delta;


	}
	*stdDev = sqrt(delta_sum/(double)num_point);
    *maxDev = max_delta;

    delete lista;
	for(int i = 1; i <= ndata; i++) {
		delete x[i];                                                   		 
	}
    delete x;
    delete y;
    delete sig;
    delete a;
}


/*    
 *      Similar to GetXShiftAndRotateParams()
 *      except for that 
 *      a[1] = k1, a[2] = k2, a[3] = k3; a[4] = x0, a[5]=y0, a[6] = z0, a[7] = rot_x, a[8] = rot_z; a[9] = del_theta;
 *      {k1, k2, k3) are fixed the inital value will not be changed after fitting
 * 
 *      [x0,y0,z0,rot_z, rot_z, del_theta] will be adjusted
 *
 *      01/29/02 
 */

void CStdMath::GetXShiftAndRotateParams4(int num_point, double  *p[6], double init_par[9], double *k1, double *k2, double *k3, double *x0, double *y0, double *z0, double *rot_x, double *rot_z, double *del_theta, double *stdDev, double *maxDev)
{
  
	// declaration & initialization
	int ndata = num_point;		// number of input points x[i][k]; i = 1,... ndata
	int nd = 6;					// nd: number of given x[i][k],  k = 1,2 ,3;
	int ma = 9;					// number of parameters to be determined after fitting
	int mfit = 6;				// only adjust 5 parameters from 8 initial values

	int* lista;					// lista[1..ma]
	double ** x;					// x[1..ndata][1..nd]	
	double * y;					// y[1..ndata]
	double * sig;					// sig[1..ndata]
	double * a;					// a[1..ma]


	// memory allocation
	lista = new int[ma+1];
	 
	x = new double*[ndata+1];
	for(int i = 1; i <= ndata; i++) {
		x[i] = new double[nd+1];                                                   		 
	}

	y = new double[ndata+1];
	sig = new double[ndata+1];
	a = new double[ma+1];
  
	for(int k = 1; k <= ma; k++) {
	//	lista[k] = k;			// numbers the parameters a[1..ma] such that 
 						// the first mfit elements correspond to values actually being adjusted
		// move the first three element in a[] array to the last three 
		if (k<7){
            lista[k] = k + 3;
		} else {
			lista[k] = k - 6;
		}
	}
  
	for(int k = 1; k <= ndata; k++) {
		sig[k] = 1;				// standard deviations, this is not experiment data, so set to 1
	    y[k] = 0.0;             // y = 0
	}

    // assign x value
    for(int i = 1; i <= ndata; i++) {
		for (int j = 1; j <= nd; j++){
			x[i][j] = (double) p[i-1][j-1];
		}
    	
	}


   // initial value of parameters
	a[1] = init_par[0];
	a[2] = init_par[1];
	a[3] = init_par[2];
	a[4] = init_par[3];
	a[5] = init_par[4];
	a[6] = init_par[5];
	a[7] = init_par[6];
	a[8] = init_par[7];
	a[9] = init_par[8];
 
 	NonLinearLeastSquare(x, y, sig, ndata, a, ma, lista, mfit, 10);
    // HOW TO USE NLLS:
	// (x,y): input point pair y= f(ax); normally y = 0;
	// sig: weighting factors. normally sig = 1;
	// ndata: number of data point.
	// a: coefficients to be determined -- outcome
	// ma: number of a
	// mfit: number of a to be fitted
	// lista: array of a [1 ...ma];

    
     
       *k1 = a[1];
	   *k2 = a[2];
	   *k3 = a[3];

	   *x0 = a[4];
	   *y0 = a[5];
	   *z0 = a[6];

	   *rot_x = a[7];
	   *rot_z = a[8];

	   *del_theta = a[9];


	// calculate the fitting error 
	double max_delta = -10000.0;
    double dd, delta, delta_sq, delta_sum;
	double x_prime, y_prime, z_prime;

	delta_sum = 0.0;


   
    for (int i = 0; i< num_point; i++){
	    
            // fitting error to sphare center
		    double xi = p[i][0];
			double yi = p[i][1];
			double zi = p[i][2];
			double thetai = p[i][3];
			double Li = p[i][4];
			double R = p[i][5];
			double kk1 = a[1];
			double kk2 = a[2];
			double kk3 = a[3];
			double xx0 = a[4];
			double yy0 = a[5];
			double zz0 = a[6];
			double rotx = a[7];
			double rotz = a[8];
			double deltheta = a[9];
		    double thetai_prime;
		    
			thetai_prime = thetai + deltheta;

		    x_prime = cos(thetai_prime)*(xi + kk1*Li-rotx)-sin(thetai_prime)*(zi + kk3*Li - rotz) + rotx;
			y_prime = yi + kk2*Li;
			z_prime = sin(thetai_prime)*(xi + kk1*Li-rotx)+cos(thetai_prime)*(zi + kk3*Li - rotz) + rotz;
		    
			dd = sqrt((x_prime-xx0)*(x_prime-xx0)+(y_prime - yy0)*(y_prime - yy0) + (z_prime - zz0)*(z_prime - zz0)); 
			delta = fabs(dd - R);
			delta_sq = (dd - R)*(dd - R);
            delta_sum = delta_sum + delta_sq;
			if (max_delta < delta) max_delta = delta;


	}
	*stdDev = sqrt(delta_sum/(double)num_point);
    *maxDev = max_delta;

    delete lista;
	for(int i = 1; i <= ndata; i++) {
		delete x[i];                                                   		 
	}
    delete x;
    delete y;
    delete sig;
    delete a;
}



/*    
 *      The same as GetXShiftAndRotateParams()
 *      except for that 
 *      a[1] = k1, a[2] = k2, a[3] = k3; a[4] = x0, a[5]=y0, a[6] = z0, a[7] = rot_x, a[8] = rot_z;
 *      {k1, k2, k3, x0, y0, z0} are fixed the inital value will not be changed after fitting
 *      {rot_x, rot_z } will be adjusted 
 * 
 *      01/30/02 
 */

void CStdMath::GetXShiftAndRotateParams5(int num_point, double  *p[6], double init_par[8], double *k1, double *k2, double *k3, double *x0, double *y0, double *z0, double *rot_x, double *rot_z, double *stdDev, double *maxDev)
{
  
	// declaration & initialization
	int ndata = num_point;		// number of input points x[i][k]; i = 1,... ndata
	int nd = 6;					// nd: number of given x[i][k],  k = 1,2 ,3;
	int ma = 8;					// number of parameters to be determined after fitting
	int mfit = 2;				// only adjust 2 parameters from 8 initial values

	int* lista;					// lista[1..ma]
	double ** x;					// x[1..ndata][1..nd]	
	double * y;					// y[1..ndata]
	double * sig;					// sig[1..ndata]
	double * a;					// a[1..ma]


	// memory allocation
	lista = new int[ma+1];
	 
	x = new double*[ndata+1];
	for(int i = 1; i <= ndata; i++) {
		x[i] = new double[nd+1];                                                   		 
	}

	y = new double[ndata+1];
	sig = new double[ndata+1];
	a = new double[ma+1];
  
	for(int k = 1; k <= ma; k++) {
	//	lista[k] = k;			// numbers the parameters a[1..ma] such that 
 						// the first mfit elements correspond to values actually being adjusted
		// move the first three element in a[] array to the last three 
		if (k<3){
            lista[k] = k + 6;
		} else {
			lista[k] = k - 2;
		}
	}
  
	for(int k = 1; k <= ndata; k++) {
		sig[k] = 1;				// standard deviations, this is not experiment data, so set to 1
	    y[k] = 0.0;             // y = 0
	}

    // assign x value
    for(int i = 1; i <= ndata; i++) {
		for (int j = 1; j <= nd; j++){
			x[i][j] = (double) p[i-1][j-1];
		}
    	
	}


   // initial value of parameters
	a[1] = init_par[0];
	a[2] = init_par[1];
	a[3] = init_par[2];
	a[4] = init_par[3];
	a[5] = init_par[4];
	a[6] = init_par[5];
	a[7] = init_par[6];
	a[8] = init_par[7];
 
 	NonLinearLeastSquare(x, y, sig, ndata, a, ma, lista, mfit, 9);
    // HOW TO USE NLLS:
	// (x,y): input point pair y= f(ax); normally y = 0;
	// sig: weighting factors. normally sig = 1;
	// ndata: number of data point.
	// a: coefficients to be determined -- outcome
	// ma: number of a
	// mfit: number of a to be fitted
	// lista: array of a [1 ...ma];

    
     
       *k1 = a[1];
	   *k2 = a[2];
	   *k3 = a[3];

	   *x0 = a[4];
	   *y0 = a[5];
	   *z0 = a[6];

	   *rot_x = a[7];
	   *rot_z = a[8];


	// calculate the fitting error 
	double max_delta = -10000.0;
    double dd, delta, delta_sq, delta_sum;
	double x_prime, y_prime, z_prime;

	delta_sum = 0.0;


   
    for (int i = 0; i< num_point; i++){
	    
            // fitting error to sphare center
		    double xi = p[i][0];
			double yi = p[i][1];
			double zi = p[i][2];
			double thetai = p[i][3];
			double Li = p[i][4];
			double R = p[i][5];
			double kk1 = a[1];
			double kk2 = a[2];
			double kk3 = a[3];
			double xx0 = a[4];
			double yy0 = a[5];
			double zz0 = a[6];
			double rotx = a[7];
			double rotz = a[8];
		     
		    
		    x_prime = cos(thetai)*(xi + kk1*Li-rotx)-sin(thetai)*(zi + kk3*Li - rotz) + rotx;
			y_prime = yi + kk2*Li;
			z_prime = sin(thetai)*(xi + kk1*Li-rotx)+cos(thetai)*(zi + kk3*Li - rotz) + rotz;
		    
			dd = sqrt((x_prime-xx0)*(x_prime-xx0)+(y_prime - yy0)*(y_prime - yy0) + (z_prime - zz0)*(z_prime - zz0)); 
			delta = fabs(dd - R);
			delta_sq = (dd - R)*(dd - R);
            delta_sum = delta_sum + delta_sq;
			if (max_delta < delta) max_delta = delta;


	}
	*stdDev = sqrt(delta_sum/(double)num_point);
    *maxDev = max_delta;

    delete lista;
	for(int i = 1; i <= ndata; i++) {
		delete x[i];                                                   		 
	}
    delete x;
    delete y;
    delete sig;
    delete a;
}

/*  Input:  a[][3], an array of points (x, y, z).
*           bum_point, number of points;
*   Output: (kx, ky, kz), defining line direction.
*           (x0, y0, z0), a point in the line.
*           Both outputs determine a line.
*  Notes for Approach: The parameters is solved by minimizing the sum of the square of 
*                      distance between every point to the line. It can be proved that 
*                       the optimal line direction is determined by the covariance matrix
*                       of the data (x, y, z). And the line must pass the geometry center
*                       of these points.
*  Chunming Li, 02/04/02  */
void CStdMath::LineFitting3D(int num_point, double *a[], double *kx, double *ky, double *kz, double *x0, double *y0, double *z0, double *std, double *maxErr)
{
    double *mean = new double[3];
    double **cov = new double*[3];
    int i, j, k;
    for (k = 0; k < 3; k++) {
        cov[k] = new double[3];
    }
    for (int i = 0; i < 3; i++) {
        mean[i] = 0;
        for (j = 0; j < 3; j++) {
            cov[i][j] = 0;            
        }
    }

    for (k = 0; k < num_point; k++) {
        for (int i = 0; i < 3; i++) {
            mean[i] += a[k][i];
            for (j = 0; j< 3; j++) {
                cov[i][j] += a[k][i] * a[k][j];
            }
        }       
    }
    for (int i = 0; i < 3; i++) {
        mean[i] = mean[i] / num_point;
        for (j = 0; j < 3; j++) {
            cov[i][j] = cov[i][j] / num_point;            
        }
    }

    for (int i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            cov[i][j] = cov[i][j] - mean[i] * mean[j];            
        }
    }
    *x0 = mean[0];
    *y0 = mean[1];
    *z0 = mean[2];
    double lineDirection[3];
    MaxEigenvector(cov, 3, lineDirection);
    double sum = 0, dist;
    *maxErr = 0;
	//*x0 = pointInLine[0], *y0 = pointInLine[1], *z0 = pointInLine[2];
    for (k = 0; k < num_point; k++) {
        dist = DistanceFromPointToLine(a[k][0], a[k][1], a[k][2], 
                                       lineDirection, *x0, *y0, *z0);
        sum += dist * dist;
        if (*maxErr < dist) *maxErr = dist;
    }
    *std = sqrt(sum / num_point);
    *kx = lineDirection[0];
    *ky = lineDirection[1];
    *kz = lineDirection[2];

	double test;
	double proj_on_line;
	proj_on_line = ((*x0) * (*kx) + (*y0) * (*ky) + (*z0) * (*kz))
		           / ((*kx) * (*kx) + (*ky) * (*ky) + (*kz) * (*kz));
	test = sqrt((*x0) * (*x0) + (*y0) * (*y0) + (*z0) * (*z0));
	*x0 -= proj_on_line * (*kx);
	*y0 -= proj_on_line * (*ky);
	*z0 -= proj_on_line * (*kz);
	test = sqrt((*x0) * (*x0) + (*y0) * (*y0) + (*z0) * (*z0));
	

    delete[] mean;
    for (k = 0; k < 3; k++) {
        delete[] cov[k];
    }
    delete[] cov;
}

/// Use svd to get eigenvector associated with the maximum eigenvalue of a semi-positive matrix.
void CStdMath::MaxEigenvector(double **A, int nrow, double *eigenVector)
{
	int rl, rh, cl, ch;
 	double** u;
	double*	d;
	double** v;

	rl = 1;
	rh = nrow;
	cl = 1;
	ch = nrow;

	u = matrix(rl, rh, cl, ch);
	d = vector(cl, ch);
	v = matrix(rl, rh, cl, ch);

	int k, j;
	for(k = rl; k <= rh; k++) {
		for(j = cl; j <= ch; j++) {
			u[k][j] = (double)A[k-1][j-1];    
		}
	}
	svdcmp(u, rh, ch, d, v);
	int ind_max = cl;
	double sMax = d[cl], sMin = d[ch];
	for (k = 0; k < ch; k++) {
		if(sMax < d[k+1])
		{
			sMax = d[k+1];
			ind_max = k+1;
		}	
	}
	for (j = 0; j < nrow; j++) {
		eigenVector[j] = (double)u[j+1][ind_max];
		double test = eigenVector[j];
	}
	free_matrix(u,rl,rh,cl,ch);
	free_matrix(v,rl,rh,cl,ch);
	free_vector(d, cl,ch);
}

// Used to calculated standard deviation in LineFitting3D()
// Chunming Li, 02/20/02
double CStdMath::DistanceFromPointToLine(double x, double y, double z, double *lineDirection, double x0, double y0, double z0)
{
    double dist, dx, dy, dz;
    double proj;
    proj = (x - x0) * lineDirection[0] 
           + (y - y0) * lineDirection[1]
           + (z - z0) * lineDirection[2];  /// inner product (X-X0).*D
    dx = x - x0 - proj * lineDirection[0];
    dy = y - y0 - proj * lineDirection[1];
    dz = z - z0 - proj * lineDirection[2];
    dist = sqrt(dx * dx + dy * dy + dz * dz);
    return dist;
}

/*  Input:  a[][3], an array of points (x, y, z).
*           bum_point, number of points;
*   Output: (kx, ky, kz), defining line direction.
*  Notes for Approach: The parameters is solved by minimizing the sum of the square of 
*                      distance between every point to the line. It can be proved that 
*                       the optimal line direction is determined by the covariance matrix
*                       of the data (x, y, z). And the line must pass the geometry center
*                       of these points.
*  Chunming Li, 02/04/02  */

void CStdMath::LineFitting3D(int num_point, double *a[], 
                            double *kx, double *ky, double *kz, 
                            double *std, double *maxErr)
{
    double *mean = new double[3];
    double **cov = new double*[3];
    int i, j, k;
    for (k = 0; k < 3; k++) {
        cov[k] = new double[3];
    }
    for (int i = 0; i < 3; i++) {
        mean[i] = 0;
        for (j = 0; j < 3; j++) {
            cov[i][j] = 0;            
        }
    }

    for (k = 0; k < num_point; k++) {
        for (int i = 0; i < 3; i++) {
            for (j = 0; j< 3; j++) {
                cov[i][j] += a[k][i] * a[k][j];
            }
        }       
    }
    for (int i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            cov[i][j] = cov[i][j] / num_point;            
        }
    }

    for (int i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            cov[i][j] = cov[i][j];          
        }
    }
    double lineDirection[3];
    MaxEigenvector(cov, 3, lineDirection);
    double sum = 0, dist;
    *maxErr = 0;
    for (k = 0; k < num_point; k++) {
        dist = DistanceFromPointToLine(a[k][0], a[k][1], a[k][2], 
                                       lineDirection, 0, 0, 0);
        sum += dist * dist;
        if (*maxErr < dist) *maxErr = dist;
    }
    *std = sqrt(sum / num_point);
    *kx = lineDirection[0];
    *ky = lineDirection[1];
    *kz = lineDirection[2];
    delete[] mean;
    for (k = 0; k < 3; k++) {
        delete[] cov[k];
    }
    delete[] cov;

}


// Used in line and plane fitting
// Chunming Li, 02/16/02
void CStdMath::SampleCovarianceAndMean(int numSample, int dim, double *sampleData[], double *mean, double **cov)
{  
    int i, j, k;
    for (int i = 0; i < dim; i++) {
        mean[i] = 0;
    }
    for (int i = 0; i < dim; i++) {
        for (j = 0; j < dim;  j++) {
            cov[i][j] = 0;
        }
    }
    for (k = 0; k < numSample; k++) {
        for (int i = 0; i < dim; i++) {
            mean[i] += sampleData[k][i];
        }
    }
    for (int i = 0; i < dim; i++) {
        mean[i] = mean[i] / numSample;
    }
    for (k = 0; k < numSample; k++) {        
        for (int i = 0; i < dim; i++) {
            for (j = 0; j < dim;  j++) {
                cov[i][j] += sampleData[k][i] * sampleData[k][j];
            }                       
        }
    }
    for (int i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            cov[i][j] = cov[i][j] / numSample - mean[i]*mean[j];
        }
    }
}

/*  Input:  a[][3], an array of points (x, y, z).
*           bum_point, number of points;
*   Output: nx, ny, nz, and d, defining a plane.
*  Notes for Approach: The parameters is solved by minimizing the sum of the square of 
*                      distance between every point to the plane. It can be proved that 
*                       the optimal plane has the normal being the minimum eigenvalue eigenvector
*                       of the covariance matrix of the group of (x, y, z). And the plane must pass 
*                       the geometry center of these points.
*  Chunming Li, 02/20/02  */
void CStdMath::PlaneFitting(int num_point, double *a[], double *nx, double *ny, double *nz, double *d, double *stdDev, double *maxDev)
{
    double *center = new double[3];
    double **cov = new double*[3];
    int i, j, k;
    for (k = 0; k < 3; k++) {
        cov[k] = new double[3];
    }
    SampleCovarianceAndMean(num_point, 3, a, center, cov);
    double *normal = new double[3];
    MinEigenvector(cov, 3, normal);
    *nx = normal[0];
    *ny = normal[1];
    *nz = normal[2];
    *d = normal[0] * center[0] + normal[1] * center[1] + normal[2] * center[2];
    double dis;
	double max_d = 0.0;
	double dev = 0.0;
    double *X = new double[3];
    double sum = 0;
    for (int i = 0; i < 3; i++) {
        sum += normal[i] * center[i];
    }
    for (int i = 0; i < 3; i++) {
        X[i] = normal[i] / (*d);
    }
    dis = 0;
	for (int i=0; i<num_point; i++){
		  // distance from measured points to the fitted plane
		  dis = fabs((a[i][0]-center[0])*(*nx)+(a[i][1]-center[1])*(*ny)
                     +(a[i][2]-center[2])*(*nz));
		  if (dis > max_d) max_d = dis;
		  dev = dev + dis * dis;
    }
    dev = sqrt(dev/double(num_point));
    *stdDev = dev;
    *maxDev = max_d;
    delete[] normal;
    delete[] center;
	delete[] X;
    for (k = 0; k < 3; k++) {
        delete[] cov[k];
    }
    delete[] cov;

}
// Used to find the normal of the fitted plane in PlaneFitting()
// Chunming Li, 02/20/02
void CStdMath::MinEigenvector(double **A, int nrow, double *eigenVector)
{
    int rl, rh, cl, ch;
 	double** u;
	double*	d;
	double** v;

	rl = 1;
	rh = nrow;
	cl = 1;
	ch = nrow;

	u = matrix(rl, rh, cl, ch);
	d = vector(cl, ch);
	v = matrix(rl, rh, cl, ch);

	int k, j;
	for(k = rl; k <= rh; k++) {
		for(j = cl; j <= ch; j++) {
			u[k][j] = (double)A[k-1][j-1];    
		}
	}
	svdcmp(u, rh, ch, d, v);
	int ind_min = cl;
	double sMin = d[cl];
	for (k = 0; k < ch; k++) {
		if(sMin > d[k+1])
		{
			sMin = d[k+1];
			ind_min = k+1;
		}	
	}
    double sum = 0;
	for (j = 0; j < nrow; j++) {
		eigenVector[j] = (double)u[j+1][ind_min];
        sum += eigenVector[j] * eigenVector[j];
		double test = eigenVector[j];
	}
    for (j = 0; j < nrow; j++) {
		eigenVector[j] = eigenVector[j] / sqrt(sum);
	}
	free_matrix(u,rl,rh,cl,ch);
	free_matrix(v,rl,rh,cl,ch);
	free_vector(d, cl,ch);
}

// Rotate a point around an axis
/* Input:  xyz[3], a point to be rotated,
*          rot_ang, the angle of rotation,
*          axis[6], axis[0], axis[1], and axis[2] give the direction of line,
*                   axis[3], axis[4], and axis[5] are the coordinates of a point on the line
*  Output: xyz_out[3], point after rotation.
*/                                  /* Chunming Li, 02/20/02 */
void CStdMath::AxisRotation(double xyz[], double rot_ang, double axis[], double xyz_out[])
{
    double c, s, ux, uy, uz, x0, y0, z0;
    c = cos(rot_ang);
    s = sin(rot_ang);
    ux = axis[0];
    uy = axis[1];
    uz = axis[2];
    x0 = axis[3];
    y0 = axis[4];
    z0 = axis[5];    

    double xyz_shift[3], xyz_shift_out[3];
    xyz_shift[0] = xyz[0] - x0;
    xyz_shift[1] = xyz[1] - y0;
    xyz_shift[2] = xyz[2] - z0;
    AxisRotation(xyz_shift, rot_ang, ux, uy, uz, xyz_shift_out);
    xyz_out[0] = xyz_shift_out[0] + x0;
    xyz_out[1] = xyz_shift_out[1] + y0;
    xyz_out[2] = xyz_shift_out[2] + z0;
}

// Rotate a point around an axis
/* Input:  xyz[3], a point to be rotated,
*          rot_ang, the angle of rotation,
*          ux, uy, uz, give the direction of the axis passing the origin
*  Output: xyz_out[3], point after rotation.
*/                                  /* Chunming Li, 02/20/02 */
void CStdMath::AxisRotation(double xyz[], double rot_ang, double ux, double uy, double uz, double xyz_out[])
{
    double c, s, x0, y0, z0;
    c = cos(rot_ang);
    s = sin(rot_ang);
    
    double R[3][3];
    R[0][0] = c + (1 - c) * ux * ux;
    R[0][1] = (1 - c) * uy * ux - s * uz;
    R[0][2] = (1 - c) * uz * ux + s * uy;
    R[1][0] = (1 - c) * ux * uy + s * uz;
    R[1][1] = c + (1 - c) * uy * uy;
    R[1][2] = (1 - c) * uz * uy - s * ux;
    R[2][0] = (1 - c) * ux * uz - s * uy;
    R[2][1] = (1 - c) * uy * uz + s * ux;
    R[2][2] = c + (1 - c) * uz * uz;

    xyz_out[0] = R[0][0] * xyz[0] + R[0][1] * xyz[1] + R[0][2] * xyz[2];
    xyz_out[1] = R[1][0] * xyz[0] + R[1][1] * xyz[1] + R[1][2] * xyz[2]; 
    xyz_out[2] = R[2][0] * xyz[0] + R[2][1] * xyz[1] + R[2][2] * xyz[2];
}

/**
 *   LineFitting2D() is written for the special case in current project, where
 *   all fitted lines are almost vertical. One may expect bigger than expected
 *   error when use this function to fit horizental lines.
 *   line function in the format of X = A*Y + B is used, where Y value is used 
 *   to calculate X value for the fitted output
 *   Input:
 *   @numPnt         number of poits used for fitting
 *   @*coefficients  A, B value
 *   @*X             X coordinates
 *   @*Y             Y coordinates
 *   @*stdev         standard square root deviation
 */

void CStdMath::LineFitting2D(int     numPnt, 
                            double *coefficients, 
                            double *x, 
                            double *y, 
                            double *stdev)                            
{
    double **A;
    double *B;
    int      i;
    double   maxError;

    if(numPnt < 1)
    {
        return;
    }

    A = new double *[numPnt];

    for(int i = 0; i<numPnt; i++)
    {
        A[i] = new double[2];
        A[i][0]=A[i][1]=0;
    }
    for(int i = 0; i<numPnt; i++)
    {
        A[i][0] = y[i];
        A[i][1] = 1;
     }

    LinearLeastSquare(numPnt, 2, A, x, coefficients, stdev, &maxError );

    // deallocate memory
    for(int i = 0; i<numPnt; i++)
    {
        delete A[i];
    }
        delete A;
}


// Calclulate transform matrix M such that
// M*T1 = T2;
// where M is 4x4 matrix
// T1 and T2 are 4xN matrxi
// M = T2*T1'*inv(T1*T1')
// input: double T1[]:       pointer to 4xN matrix in the form of T11,... Tn1,... T41 .. T4n
//		  double T2[]:       pointer to 4xN matrix in the form of ....
//        int	iNum:		 number of point 
// output: M[16];            Transform matrix in the form of M11 ....M44
bool CStdMath::CalculateTransform(int iNum, double T1[], double T2[], double M[16])
{
 
	// Calculate transpose of T1
	double *Transp_T1;
    double temp[16];
	double temp2[16];
	double **inv_temp;
	
	Transp_T1 = new double[4*iNum];
  
	inv_temp = new double*[4];

	for (int ii=0; ii<4; ii++){
		inv_temp[ii] = new double[4];
	}


	TransposeMatrix(T1, 4, iNum, Transp_T1);

	// create T1*T1' -->temp[16]
	MultiplyMatrix(T1, 4, iNum, Transp_T1, 4, temp);

	// inverse temp --> temp;
	for(int kk=0; kk<4; kk++){
		for(int jj=0; jj<4; jj++){
			inv_temp[kk][jj] = temp[kk*4+jj];
		}
	}

	MatrixInverse(4, inv_temp);
	
	for(int kk=0; kk<4; kk++){
		for(int jj=0; jj<4; jj++){
			temp[kk*4+jj] = inv_temp[kk][jj];
		}
	}

	// Calculate T2*T1' --> temp2
	MultiplyMatrix(T2, 4, iNum, Transp_T1, 4, temp2);

   // calculate M  = T2*T1'*inv(T1*T1')
   //           M = temp2 *temp

	MultiplyMatrix(temp2, 4, 4, temp, 4, M);

	delete 	[]Transp_T1;
	
	for (int ii=0; ii<4; ii++){
		delete []inv_temp[ii];
	} delete []inv_temp;
   
	
	return true;
}

// Calclulate transform matrix M such that
// T1*M = T2; or M = pinv(T1)*T2;
// where M is 4x4 matrix
// T1 and T2 are 4xN matrxi
// M = inv(T1'*T1)*(T1'*T2)
// input: double T1[]:       pointer to 4xN matrix in the form of T11,... Tn1,... T41 .. T4n
//		  double T2[]:       pointer to 4xN matrix in the form of ....
//        int	iNum:		 number of point 
// output: M[16];            Transform matrix in the form of M11 ....M44
bool CStdMath::CalculateTransform2(int iNum, double T1[], double T2[], double M[16])
{
 
	// Calculate transpose of T1
	double *Transp_T1;
    double temp[16];
	double temp2[16];
	double **inv_temp;
	
	Transp_T1 = new double[4*iNum];
  
	inv_temp = new double*[4];

	for (int ii=0; ii<4; ii++){
		inv_temp[ii] = new double[4];
	}


	TransposeMatrix(T1, 4, iNum, Transp_T1);

	// create T1*T1' -->temp[16]
	// MultiplyMatrix(T1, 4, iNum, Transp_T1, 4, temp);
	// create T1'*T1 -->temp[16]
    MultiplyMatrix(Transp_T1, 4, iNum, T1, 4, temp);

	// inverse temp --> temp;
	for(int kk=0; kk<4; kk++){
		for(int jj=0; jj<4; jj++){
			inv_temp[kk][jj] = temp[kk*4+jj];
		}
	}

	MatrixInverse(4, inv_temp);
	
	for(int kk=0; kk<4; kk++){
		for(int jj=0; jj<4; jj++){
			temp[kk*4+jj] = inv_temp[kk][jj];
		}
	}

	// Calculate T2*T1' --> temp2
	// MultiplyMatrix(T2, 4, iNum, Transp_T1, 4, temp2);
	// Calculate T1'*T2 --> temp2
	MultiplyMatrix(Transp_T1, 4, iNum, T2, 4, temp2);


   // calculate M  = inv(T1'*T1)*(T1'*T2)
   //           M = temp*temp2

	MultiplyMatrix(temp, 4, 4, temp2, 4, M);

	delete 	[]Transp_T1;
	
	for (int ii=0; ii<4; ii++){
		delete []inv_temp[ii];
	} delete []inv_temp;
   
	
	return true;
}


// Calclulate transform matrix M such that
// T1*M = T2; or M = pinv(T1)*T2;
// where M is 4x4 matrix
// T1 and T2 are 4x4 matrxi
// 
// input: double T1[]:       pointer to 4x4 matrix in the form of T11,...T14; .. T41 .. 
//		  double T2[]:       pointer to 4x4 matrix in the form of ....
// output: M[16];            Transform matrix in the form of M11 ....M44
bool CStdMath::CalculateTransform3(double T1[], double T2[], double M[16])
{
 
    double temp[16];
	double **inv_temp;
 
	inv_temp = new double*[4];
	for (int ii=0; ii<4; ii++){
		inv_temp[ii] = new double[4];
	}


	// inverse T1 --> temp;
	for(int kk=0; kk<4; kk++){
		for(int jj=0; jj<4; jj++){
			inv_temp[kk][jj] = T1[kk*4+jj];
		}
	}

	MatrixInverse(4, inv_temp);
	
	for(int kk=0; kk<4; kk++){
		for(int jj=0; jj<4; jj++){
			temp[kk*4+jj] = inv_temp[kk][jj];
		}
	}

	// Calculate inv(T1)*T2 --> temp*T2
	MultiplyMatrix(temp, 4, 4, T2, 4, M);


	for (int ii=0; ii<4; ii++){
		delete []inv_temp[ii];
	} delete []inv_temp;
   
	
	return true;
}


// Functionality:	transpose matrix by row to col, col to row
// Input:
//      double[][] m1 ---matrix 1
//		int row ---- row of matrix 1
//		int col ---- col of  matrix 1 
//		double[][] m2 ---matrix 2
bool CStdMath::TransposeMatrix(double  m1[], int row, int col, double m2[])
{
	for( int i=0; i<row; i++)
		for(int j=0; j<col; j++)
			m2[(row*j)+i]=m1[(col*i)+j];
	return true;
}


/*This function returns a uniform random variable between 0.0 and 1.0 using a system-supplied
routine rand().  Set idum to any negative value to initialize or reinitialize the sequence*/
double CStdMath::ran0(int *idum)
{
	static double y, maxran, v[98];  //the exact number 98 is unimportant
	double dum;
	static int iff=0;
	int j;
	

	if(*idum<0 || iff==0){
		iff=1;
		maxran=(double)(RAND_MAX+1.0);
		srand((unsigned int)(*idum));
		*idum=1;
		for(j=1;j<=97; j++) dum=rand();  //exercise the system routine, especially important 
		//if the system's multiplier is small
		for(j=1; j<=97; j++) v[j]=rand();  //then same 97 values
		y=rand();
	}
	j=(int)(1+97.0*y/maxran);  //this is where we start if not initializing. Use the previously
	//saved random number y to get an index j between 1 and 97.
	if(j>97 || j<1) nrerror("RAN0: This cannot happen.");
	y=v[j];
	v[j]=rand();
	return y/maxran;
}

/*This function returns a normally distributed random variable with zero mean and unit
variance, using ran0(idum) as the source of uniform random variable*/
double CStdMath::gaussdev(int *idum)
{
	static int iset=0;
	static double gset;
	double fac, r, v1, v2;

	if(iset==0){
		do{
			v1=2.0*ran0(idum)-1.0;    //pick two uniform numbers in the square extending
			v2=2.0*ran0(idum)-1.0;    //from -1 to +1 in each direction, see if they are
			r=v1*v1+v2*v2;            // in the unit circle, and if they are not, try again
		}while(r>=1.0 || r<=0);  
		fac=sqrt(-2.0*log(r)/r);
		gset=v1*fac;    //Box-Muller transformation
		iset=1;  //set flag
		return v2*fac;
	}else{  //we have an extra gaussian ready, so unset the flag, and return it
		iset=0;
		return gset;
	}
}

	
// The same as NonLinearLeastSquare() 
// excep that 
// pass two arguments
// int num_line, double **line_pars
// Add 5% random perturd for initial value 
// and compare chisquare, final optimization is deteremined based on
// minimization chisquare
// 
//  02/16/04
void CStdMath::NonLinearLeastSquare_modelfit(double** x, double y[], double sig[], int ndata, 
			double a[], double fiterror[2], int ma, int lista[], int mfit, int num_line, double **line_pars)
{

 	double** covar;				// covar[1..ma][1..ma]
	double** alpha;				// alpha[1..ma][1..ma]
	double* chisq;
	double* chisq_init;	
	double* alamda;
	int maxLoop;
	double  chiTol;
	double chisq_bk;
	double old_chi;
	int j;

  	covar = new double*[ma+1];
	alpha = new double*[ma+1];

 	for(int i = 1; i <= ma; i++) {
		covar[i] = new double[ma+1];                                                   
		alpha[i] = new double[ma+1];                                                   
	}

	chisq = new double();
	chisq_init = new double();
	alamda = new double();

	// initialization
	*alamda = -100.0f;			// must be initialized to less than zero
	*chisq = 0.0f;
	*chisq_init = 0.0f;

	/////////////////////////////////////////////////////////////////////////
	// Evaulate the initial --new added 02/18/04 -QT
	static double *beta;
	beta=vector(1,ma);
	mrqcof_modelfit(x,y,sig,ndata,a,ma,lista,mfit,alpha,beta,chisq_init,num_line, line_pars);
    //////////////////////////////////////////////////////////////////////////	

	// on the first call provide an initial guess for the parameters a, and set alamda < 0
	mrqmin_modelfit(x, y, sig, ndata, a, ma, lista, mfit, covar, alpha, chisq, alamda, num_line, line_pars);

	int count = 0;	// use count to control max loop
	
	
      maxLoop = 500;
	  chiTol = (double)0.1;

    chisq_bk = 10000.0;

    // get out loop when the residual is less than the threhold
    while(fabs(*chisq) > chiTol && count < maxLoop) {
	    old_chi = *chisq;
       	mrqmin_modelfit(x, y, sig, ndata, a, ma, lista, mfit, covar, alpha, chisq, alamda, num_line, line_pars);
		if((fabs(*chisq-old_chi)<0.00001)&&(*alamda>1.0e12))  //stop condition, numbers should be customized
		break;
	   
		count++;
	}

 //	*fiterror = *chisq;
	fiterror[0] = *chisq;
	fiterror[1] = *chisq_init;





    // clean up memory -- QT 6/5/01
 	for(int i = 1; i <= ma; i++) {
		delete []covar[i];                                                   
		delete []alpha[i];                                                   
	}

	delete []covar;
	delete []alpha;

	delete chisq;
	delete chisq_init;
	delete alamda;

}



// INPUT: 
// (alphe, beta, gamma, tx, ty, tz) describe a rigid body transform
//  where aplphe, beta, and gamma are in degrees
// OUTPUT: a 4x4 matrix
//
void CStdMath::ConvertAngle2Matrix(double alfa, double beta, double gama, double tx, double ty, double tz, double t[4][4])
{
#define pi 3.1415926

	alfa = alfa*pi/180;
	beta = beta*pi/180;
	gama = gama*pi/180;
	
	t[0][0] = cos(alfa)*cos(beta);
	t[0][1] = cos(alfa)*sin(beta)*sin(gama)-sin(alfa)*cos(gama);
	t[0][2] = cos(alfa)*sin(beta)*cos(gama)+sin(alfa)*sin(gama);
	t[0][3] = tx;
	t[1][0] = sin(alfa)*cos(beta);
	t[1][1] = sin(alfa)*sin(beta)*sin(gama) + cos(alfa)*cos(gama);
	t[1][2] = sin(alfa)*sin(beta)*cos(gama) - cos(alfa)*sin(gama);
	t[1][3] = ty;
	t[2][0] = -sin(beta);
	t[2][1] = cos(beta)*sin(gama);
	t[2][2] = cos(beta)*cos(gama);    
    t[2][3] = tz;
	t[3][0] = 0;
	t[3][1] = 0;
	t[3][2] = 0;
	t[3][3] = 1;

}



/*******************************************************************************************
 * create circle points based on its parameters
 * Inputs: 
 *         center position (xc, yc,zc), 
 *         radius R
 *         orientation (nx, ny, nz) in space
 *         number of points (N)
 *         out put format (iyk or rpy); rpy_flag
 * Outputs:
 *         cir[][6]   containing a array of [x,y,z,nx,ny,nz] or [x,y,z r,p,y];
 * return:
 *         void
 * Author: 
 *          5.17.00 
 */
void CStdMath::CreateCirclePoints(double xc, double yc, double zc, double R, double nx, double ny, double nz, int N, int rpy_flag, double cir[][6])
{
  
    double temp_x,temp_y,temp_z,x1,y1,z1,x2,y2,z2,x,y,z;
	double angle;
    double real_nx, real_ny, real_nz, roll, pitch, yaw;

    double p[3], d[3], n[3];
    PT_RPY	Coord;


    // normalize the orientation vector
	real_nx = nx/sqrt(nx*nx + ny*ny+ nz*nz);
	real_ny = ny/sqrt(nx*nx + ny*ny+ nz*nz);
	real_nz = nz/sqrt(nx*nx + ny*ny+ nz*nz);
    

	p[0] = xc;
	p[1] = yc;
	p[2] = zc;

	n[0] = real_nx;
	n[1] = real_ny;
	n[2] = real_nz;

	/* 
	 * Need determine (give) direction vector to create RPY or coordinate
	 * since ijy info is not enogh to determine coordinate
	 */
	 
	// Assume direction vector d[] is normal to n[] or (i,j,k)
	
	if (real_nx != 0.0){
	  d[0] = -real_ny*0.5/real_nx;
	  d[1] = 0.5;
	  d[2] = 0.0;
	} else if (real_ny != 0.0){
        d[0] = 0.5;
		d[1] = -real_nx*0.5/real_ny;
	    d[2] = 0.0;
    } else { // direction is x 
		d[0] = 1.0;
	    d[1] = 0.0;
		d[2] = 0.0;
	}


	// Define a coordinate for circle based on above normal and direction vector
	Coord = PTIJK_to_PTRPY(p, d, n);


    /*
	 * roll, pitch & yaw indicate how the designed circle coordinate transfers 
	 *into world coordinate. 
	 * From world to circle coordinate reverse the direction and angle
	 */

	roll = -Coord.roll*PI/180.0;;
	pitch = -Coord.pitch*PI/180.0;
	yaw = -Coord.yaw*PI/180.0;
		
	for (int i=0; i< N; i++){
	 // create circle in origin and plane xy;
	  angle = (double) i/(double) N*2.0 *3.1415926;
	  temp_x = R*cos(angle);
	  temp_y = R*sin(angle);
	  temp_z = 0;
     
	  /* 
	   * Notes:
	   * Rotate it to align with [nx, ny, nz]
	   * rotate along x axis 
	   * y' = cos()*y + sin()*z;
	   * z' = -sin()*y + cos()*z;
	   * where cos() = nx and sin() = sqrt(1-sqr(nx));
	   * rotate along y axis
	   * x' = cos()*x - sin()*z;
	   * z' = sin()*x + cos()*z;
	   * where cos() = ny and sin() = sqrt(1-sqr(ny));
	   * rotate along z axis
	   * x' = cos()*x + sin()*y;
	   * y' = -sin()*x + cos()*y;
	   */
	   
	  /*
	   * Note: 
	   * from circle to  world
	   * a. rotate along z with roll degree
	   * b. rotate along y with pitxh degree
	   * c. rotate along x with yaw degree
	   * From world to circle
	   * a. rotate along x with -yaw degree
	   * b. rotate along y with -pitch degree
	   * c. rotate along z with -roll degree
	   * Here roll pitch, and yaw have been revesed already
	   */



	  // first rotate along x axis (yaw)
	  y1 = cos(yaw)*temp_y;
	  z1 = -sin(yaw)*temp_y;
	  x1 = temp_x;

	  // rotate along y axis  (pitch)
	  x2 = cos(pitch)*temp_x - sin(pitch)*z1;
	  z2 = sin(pitch)*temp_x + cos(pitch)*z1;
	  y2 = y1;

	  // rotate along z axis (roll)
      x = cos(roll)*x2 + sin(roll)*y2;
	  y = -sin(roll)*x2 + cos(roll)*y2;
      z = z2;

	  	  	  
      // translate to the designed center
	  cir[i][0] = x + xc;
	  cir[i][1] = y + yc;
	  cir[i][2] = z + zc;
	 
	  if (rpy_flag == 1){
		cir[i][3] = Coord.roll;
		cir[i][4] = Coord.pitch;
		cir[i][5] = Coord.yaw;
	  } else  {
		cir[i][3] = real_nx;
	    cir[i][4] = real_ny;
	    cir[i][5] = real_nz;
      }



  	}
}




/*******************************************************************************************
 * create arc points based on its parameters
 * Inputs: 
 *         center position (xc, yc,zc), 
 *         radius R
 *         orientation (nx, ny, nz) in space
 *         number of points (N)
 *         boundary points bp[9] {x1y1z1 x2y2z2 x3y3z3} three points are in the arc 
 *         out put format (iyk or rpy); rpy_flag
 * Outputs:
 *         arc[][6]   containing a array of [x,y,z,nx,ny,nz] or [x,y,z r,p,y];
 * return:
 *         void
 * Author: 
 *          08.24.2004
 */
 /*
void CStdMath::CreateArcPoints(double xc, double yc, double zc, double R, double nx, double ny, double nz, int N, double bp[9], int rpy_flag, double arc[][6])
{
  
    double temp_x,temp_y,temp_z,x1,y1,z1,x2,y2,z2,x,y,z;
	double angle;
    double real_nx, real_ny, real_nz, roll, pitch, yaw;

    double p[3], d[3], n[3];
    PT_RPY	Coord;


    // normalize the orientation vector
	real_nx = nx/sqrt(nx*nx + ny*ny+ nz*nz);
	real_ny = ny/sqrt(nx*nx + ny*ny+ nz*nz);
	real_nz = nz/sqrt(nx*nx + ny*ny+ nz*nz);
    

	p[0] = xc;
	p[1] = yc;
	p[2] = zc;

	n[0] = real_nx;
	n[1] = real_ny;
	n[2] = real_nz;

	 
	 // Need determine (give) direction vector to create RPY or coordinate
	 //since ijy info is not enogh to determine coordinate
	 //
	 
	// Assume direction vector d[] is normal to n[] or (i,j,k)
	
	if (real_nx != 0.0){
	  d[0] = -real_ny*0.5/real_nx;
	  d[1] = 0.5;
	  d[2] = 0.0;
	} else if (real_ny != 0.0){
        d[0] = 0.5;
		d[1] = -real_nx*0.5/real_ny;
	    d[2] = 0.0;
    } else { // direction is x 
		d[0] = 1.0;
	    d[1] = 0.0;
		d[2] = 0.0;
	}


	// Define a coordinate for circle based on above normal and direction vector
	Coord = PTIJK_to_PTRPY(p, d, n);

    
	 // roll, pitch & yaw indicate how the designed circle coordinate transfers 
	 //into world coordinate. 
	 // From world to circle coordinate reverse the direction and angle
	 //

	roll = -Coord.roll*PI/180.0;;
	pitch = -Coord.pitch*PI/180.0;
	yaw = -Coord.yaw*PI/180.0;
		
	for (int i=0; i< N; i++){
	 // create circle in origin and plane xy;
	  angle = (double) i/(double) N*2.0 *3.1415926;
	  temp_x = R*cos(angle);
	  temp_y = R*sin(angle);
	  temp_z = 0;
     
	  // first rotate along x axis (yaw)
	  y1 = cos(yaw)*temp_y;
	  z1 = -sin(yaw)*temp_y;
	  x1 = temp_x;

	  // rotate along y axis  (pitch)
	  x2 = cos(pitch)*temp_x - sin(pitch)*z1;
	  z2 = sin(pitch)*temp_x + cos(pitch)*z1;
	  y2 = y1;

	  // rotate along z axis (roll)
      x = cos(roll)*x2 + sin(roll)*y2;
	  y = -sin(roll)*x2 + cos(roll)*y2;
      z = z2;

 	  
      // translate to the designed center
	  arc[i][0] = x + xc;
	  arc[i][1] = y + yc;
	  arc[i][2] = z + zc;
	 
	  if (rpy_flag == TRUE){
		arc[i][3] = Coord.roll;
		arc[i][4] = Coord.pitch;
		arc[i][5] = Coord.yaw;
	  } else  {
		arc[i][3] = real_nx;
	    arc[i][4] = real_ny;
	    arc[i][5] = real_nz;
      }


  	}
}

*/


//Given a circle parameters (xc,yc,zc,R,nx,ny,nz)
// create Arc points arc[N][6]
// arc portion is defined by three sample point samplePt[9][3];
void CStdMath::CreateArcPoints(double xc, double yc, double zc, double R, double nx, double ny, double nz, int N, double samplePt[9], int rpy_flag, double arc[][6])
{
  
    double temp_x,temp_y,temp_z,x1,y1,z1,x2,y2,z2,x3,y3,z3,x,y,z;
	double angle;
	double angle1, angle2, anglemin, anglemax, anglemid;
	int N1;
    double real_nx, real_ny, real_nz, roll, pitch, yaw;
	double x0, y0, z0, xm, ym, zm, xt, yt, zt;

    double p[3], d[3], n[3];
    PT_RPY	Coord;


	x0 = samplePt[0];
	y0 = samplePt[1];
	z0 = samplePt[2];
	xm = samplePt[3];
	ym = samplePt[4];
	zm = samplePt[5];
	xt = samplePt[6];
	yt = samplePt[7];
	zt = samplePt[8];


    // normalize the orientation vector
	real_nx = nx/sqrt(nx*nx + ny*ny+ nz*nz);
	real_ny = ny/sqrt(nx*nx + ny*ny+ nz*nz);
	real_nz = nz/sqrt(nx*nx + ny*ny+ nz*nz);
    

	p[0] = xc;
	p[1] = yc;
	p[2] = zc;

	n[0] = real_nx;
	n[1] = real_ny;
	n[2] = real_nz;

	/* 
	 * Need determine (give) direction vector to create RPY or coordinate
	 * since ijy info is not enogh to determine coordinate
	 */
	 
	// Assume direction vector d[] is normal to n[] or (i,j,k)
	
	if (real_nx != 0.0){
	  d[0] = -real_ny*0.5/real_nx;
	  d[1] = 0.5;
	  d[2] = 0.0;
	} else if (real_ny != 0.0){
        d[0] = 0.5;
		d[1] = -real_nx*0.5/real_ny;
	    d[2] = 0.0;
    } else { // direction is x 
		d[0] = 1.0;
	    d[1] = 0.0;
		d[2] = 0.0;
	}


	// Define a coordinate for circle based on above normal and direction vector
	Coord = PTIJK_to_PTRPY(p, d, n);

     /*
	 * roll, pitch & yaw indicate how the designed circle coordinate transfers 
	 *into world coordinate. 
	 * From world to circle coordinate reverse the direction and angle
	 */

	roll = -Coord.roll*PI/180.0;;
	pitch = -Coord.pitch*PI/180.0;
	yaw = -Coord.yaw*PI/180.0;


// rotate
 // translate to the designed center
	  x1 = x0 - xc;
	  y1 = y0 - yc;
	  z1 = z0 - zc;
// rotate along z axis (-roll)
      x2= cos(roll)*x1 - sin(roll)*y1;
	  y2 = sin(roll)*x1 + cos(roll)*y1;
      z2 = z1;
 // rotate along y axis  (-pitch)
	  x3 = cos(pitch)*x2 + sin(pitch)*z2;
	  z3 = -sin(pitch)*x2 + cos(pitch)*z2;
	  y3 = y2;
// first rotate along x axis (-yaw)
	  y = cos(yaw)*y3 - sin(yaw)*z3;
	  z = sin(yaw)*y3 + cos(yaw)*z3;
	  x = x3;
// angle1
	  angle1=atan2(y,x);
	  if (angle1<0){
		  angle1=2.0*3.1415926+angle1;
	  }


// rotate
 // translate to the designed center
	  x1 = xt - xc;
	  y1 = yt - yc;
	  z1 = zt -zc;
// rotate along z axis (-roll)
      x2= cos(roll)*x1 - sin(roll)*y1;
	  y2 = sin(roll)*x1 + cos(roll)*y1;
      z2 = z1;
 // rotate along y axis  (-pitch)
	  x3 = cos(pitch)*x2 + sin(pitch)*z2;
	  z3 = -sin(pitch)*x2 + cos(pitch)*z2;
	  y3 = y2;
// first rotate along x axis (-yaw)
	  y = cos(yaw)*y3 - sin(yaw)*z3;
	  z = sin(yaw)*y3 + cos(yaw)*z3;
	  x = x3;
// angle2
	  angle2=atan2(y,x);
	  if (angle2<0){
		  angle2=angle2+2.0*3.1415926;
	  }
// rotate
 // translate to the designed center
	  x1 = xm - xc;
	  y1 = ym - yc;
	  z1 = zm -zc;
// rotate along z axis (-roll)
      x2= cos(roll)*x1 - sin(roll)*y1;
	  y2 = sin(roll)*x1 + cos(roll)*y1;
      z2 = z1;
 // rotate along y axis  (-pitch)
	  x3 = cos(pitch)*x2 + sin(pitch)*z2;
	  z3 = -sin(pitch)*x2 + cos(pitch)*z2;
	  y3 = y2;
// first rotate along x axis (-yaw)
	  y = cos(yaw)*y3 - sin(yaw)*z3;
	  z = sin(yaw)*y3 + cos(yaw)*z3;
	  x = x3;
// angle3
	  anglemid=atan2(y,x);
	  if (anglemid<0){
		  anglemid=anglemid+2.0*3.1415926;
	  }
//
	  if (angle1>angle2){
		  anglemin=angle2;
		  anglemax=angle1;
	  }else{
		  anglemin=angle1;
		  anglemax=angle2;
	  }
//
	  if((anglemin<anglemid)&&(anglemax>anglemid)){
		  for (int i=0; i< N; i++){
			 // create circle in origin and plane xy;
			  angle = (double) i/(double) N*(anglemax-anglemin)+anglemin;
			  temp_x = R*cos(angle);
			  temp_y = R*sin(angle);
			  temp_z = 0;
     
			  // first rotate along x axis (yaw)
			  y1 = cos(yaw)*temp_y;
			  z1 = -sin(yaw)*temp_y;
			  x1 = temp_x;

			  // rotate along y axis  (pitch)
			  x2 = cos(pitch)*temp_x - sin(pitch)*z1;
			  z2 = sin(pitch)*temp_x + cos(pitch)*z1;
			  y2 = y1;

			  // rotate along z axis (roll)
			  x = cos(roll)*x2 + sin(roll)*y2;
			  y = -sin(roll)*x2 + cos(roll)*y2;
			  z = z2;

	  	  			  
			  // translate to the designed center
			  arc[i][0] = x + xc;
			  arc[i][1] = y + yc;
			  arc[i][2] = z + zc;
			 
			  if (rpy_flag == 1){
				arc[i][3] = Coord.roll;
				arc[i][4] = Coord.pitch;
				arc[i][5] = Coord.yaw;
			  } else  {
				arc[i][3] = real_nx;
				arc[i][4] = real_ny;
				arc[i][5] = real_nz;
			  }
		  }   // end for
	  }else{
		  N1=int(N*anglemin/(anglemin+2.0*3.1415926-anglemax));
		  for (int i=0; i< N1;i++){
			  // create circle in origin and plane xy;
			  angle = (double) i/(double) N1*anglemin;
			  temp_x = R*cos(angle);
			  temp_y = R*sin(angle);
			  temp_z = 0;
     
			  // first rotate along x axis (yaw)
			  y1 = cos(yaw)*temp_y;
			  z1 = -sin(yaw)*temp_y;
			  x1 = temp_x;

			  // rotate along y axis  (pitch)
			  x2 = cos(pitch)*temp_x - sin(pitch)*z1;
			  z2 = sin(pitch)*temp_x + cos(pitch)*z1;
			  y2 = y1;

			  // rotate along z axis (roll)
			  x = cos(roll)*x2 + sin(roll)*y2;
			  y = -sin(roll)*x2 + cos(roll)*y2;
			  z = z2;

	  	  			  
			  // translate to the designed center
			  arc[i][0] = x + xc;
			  arc[i][1] = y + yc;
			  arc[i][2] = z + zc;
			 
			  if (rpy_flag == 1){
				arc[i][3] = Coord.roll;
				arc[i][4] = Coord.pitch;
				arc[i][5] = Coord.yaw;
			  } else  {
				arc[i][3] = real_nx;
				arc[i][4] = real_ny;
				arc[i][5] = real_nz;
			  }
		  }

		  for(int i=N1+1; i<N;i++){
			  // create circle in origin and plane xy;
			  angle = (double) (i-N-1)/(double) (N-N1)*(2.0*3.1415926-anglemax)+anglemax;
			  temp_x = R*cos(angle);
			  temp_y = R*sin(angle);
			  temp_z = 0;
     
			
			  // first rotate along x axis (yaw)
			  y1 = cos(yaw)*temp_y;
			  z1 = -sin(yaw)*temp_y;
			  x1 = temp_x;

			  // rotate along y axis  (pitch)
			  x2 = cos(pitch)*temp_x - sin(pitch)*z1;
			  z2 = sin(pitch)*temp_x + cos(pitch)*z1;
			  y2 = y1;

			  // rotate along z axis (roll)
			  x = cos(roll)*x2 + sin(roll)*y2;
			  y = -sin(roll)*x2 + cos(roll)*y2;
			  z = z2;

	  	  			  
			  // translate to the designed center
			  arc[i][0] = x + xc;
			  arc[i][1] = y + yc;
			  arc[i][2] = z + zc;
			 
			  if (rpy_flag == 1){
				arc[i][3] = Coord.roll;
				arc[i][4] = Coord.pitch;
				arc[i][5] = Coord.yaw;
			  } else  {
				arc[i][3] = real_nx;
				arc[i][4] = real_ny;
				arc[i][5] = real_nz;
			  }
		  }
	 }
}

/* Register one set of points to another set of points
		Quaternion-based algorithm is used
   Input:	PtsFrom		PtsFrom[numPts][3] coordinates of points
			PtsTo		PtsTo[numPts][3] coordiantes of another points
			numPts		number of points
   Output:	matrix		tranformation matrix from PtsFrom to PtsTo  
			ms			mean square point matching error
	Huagen Liu 11/5/2004		*/
void CStdMath::RegisterationPoints2Points(double **PtsFrom, double **PtsTo, long numPts, double **matrix, double *ms)
{
	Vector UFrom;		// average vector of PtsFrom
	Vector UTo;			// average vector of PtsTo
	Matrix MatrixCC(3,3);	// cross-covariance matrix of sets PtsFromand PtsTo
	double trace;
	double eigen[4];
	
	double sumX, sumY, sumZ;
	long i,j;

	// Calculate average vectors of PtsFrom and PtsTo
	sumX = sumY = sumZ = 0.0;
	for (int i=0; i<numPts; i++) {
		sumX += PtsFrom[i][0];
		sumY += PtsFrom[i][1];
		sumZ += PtsFrom[i][2];
	}
	UFrom.value[0] = sumX / numPts;
	UFrom.value[1] = sumY / numPts;
	UFrom.value[2] = sumZ / numPts;

	sumX = sumY = sumZ = 0.0;
	for (int i=0; i<numPts; i++) {
		sumX += PtsTo[i][0];
		sumY += PtsTo[i][1];
		sumZ += PtsTo[i][2];
	}
	UTo.value[0] = sumX / numPts;
	UTo.value[1] = sumY / numPts;
	UTo.value[2] = sumZ / numPts;

	// Calculate cross covariance matrix
	Matrix tmpFrom(3,1), tmpTo(1,3);
	Matrix tmp;

	MatrixCC.Zero();
	for (int i=0; i<numPts; i++) {
		tmpFrom.SetCol(0, PtsFrom[i][0]-UFrom.value[0], PtsFrom[i][1]-UFrom.value[1], PtsFrom[i][2]-UFrom.value[2]);
		tmpTo.SetRow(0, PtsTo[i][0]-UTo.value[0], PtsTo[i][1]-UTo.value[1], PtsTo[i][2]-UTo.value[2]);
		tmp = tmpFrom * tmpTo;
		MatrixCC = MatrixCC + tmp;
	}

	for (int i=0; i<3; i++)
		for (j=0; j<3; j++)
			MatrixCC.value[i][j] = MatrixCC.value[i][j] / numPts;

	trace = MatrixCC.Trace();

	Matrix tmp1(3,3);
	tmp1.Identity ();
	tmp1 = tmp1 * trace;

	tmp = MatrixCC.Transpose ();
	tmp = MatrixCC + tmp;
	tmp = tmp - tmp1;

	Vector delta(3);
	delta.value[0] = MatrixCC.value[1][2] - MatrixCC.value[2][1];
	delta.value[1] = MatrixCC.value[2][0] - MatrixCC.value[0][2];
	delta.value[2] = MatrixCC.value[0][1] - MatrixCC.value[1][0];

	Matrix Qr(4,4);
	Qr.SetRow(0, trace, delta.value[0], delta.value[1], delta.value[2]);
	Qr.SetRow(1, delta.value[0], tmp.value[0][0], tmp.value[0][1], tmp.value[0][2]);
	Qr.SetRow(2, delta.value[1], tmp.value[1][0], tmp.value[1][1], tmp.value[1][2]);
	Qr.SetRow(3, delta.value[2], tmp.value[2][0], tmp.value[2][1], tmp.value[2][2]);

	//char msg[256];
	//sprintf(msg,"%s\n", Qr.toString().c_str());
	//AfxMessageBox(msg);

	MaxEigenvector(Qr.value, 4, eigen);

	//sprintf(msg,"Eigenvalues: %.4f %.4f %.4f %.4f\n", eigen[0], eigen[1], eigen[2], eigen[3]);
	//AfxMessageBox(msg);

	Matrix R(3,3);
	R.value[0][0] = eigen[0] * eigen[0] + eigen[1] * eigen[1] - eigen[2] * eigen[2] - eigen[3] * eigen[3];
	R.value[0][1] = 2.0 * (eigen[1] * eigen[2] - eigen[0] * eigen[3]);
	R.value[0][2] = 2.0 * (eigen[1] * eigen[3] + eigen[0] * eigen[2]);
	R.value[1][0] = 2.0 * (eigen[1] * eigen[2] + eigen[0] * eigen[3]);
	R.value[1][1] = eigen[0] * eigen[0] + eigen[2] * eigen[2] - eigen[1] * eigen[1] - eigen[3] * eigen[3];
	R.value[1][2] = 2.0 * (eigen[2] * eigen[3] - eigen[0] * eigen[1]);
	R.value[2][0] = 2.0 * (eigen[1] * eigen[3] - eigen[0] * eigen[2]);
	R.value[2][1] =  2.0 * (eigen[2] * eigen[3] + eigen[0] * eigen[1]);
	R.value[2][2] = eigen[0] * eigen[0] + eigen[3] * eigen[3] - eigen[2] * eigen[2] - eigen[1] * eigen[1];

	Vector Qt;
	Qt = R * UFrom;
	Qt = UTo - Qt;

	Matrix Q(4,4);
	Q.Identity ();

	for (int i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			Q.value[i][j] = R.value[i][j];
			matrix[i][j] = R.value[i][j];
		}
		Q.value[i][3] = matrix[i][3] = Qt.value[i];
		matrix[3][i] = 0;
	}
	matrix [3][3] = 1;

	//sprintf(msg,"%s\n",Q.toString().c_str());
	//AfxMessageBox(msg);

	double sumErr = 0;
	Vector fromPt(4);
	Vector tVector;
	
	for (int i=0; i<numPts; i++) {
		fromPt.value[0] = PtsFrom[i][0];
		fromPt.value[1] = PtsFrom[i][1];
		fromPt.value[2] = PtsFrom[i][2];
		fromPt.value[3] = 1.0;
		tVector = Q * fromPt;
		sumErr +=  (tVector.value[0] - PtsTo[i][0]) * (tVector.value[0] - PtsTo[i][0]) + 
						(tVector.value[1] - PtsTo[i][1]) * (tVector.value[1] - PtsTo[i][1]) +
						(tVector.value[2] - PtsTo[i][2]) * (tVector.value[2] - PtsTo[i][2]) ;
	}
	*ms = sqrt(sumErr / numPts);
}


bool CStdMath::IsLeagleTriangle(double* pdTriangle, double dMinDistance)
{
	double dDist12 = (pdTriangle[0] - pdTriangle[3])*(pdTriangle[0] - pdTriangle[3])
		            +(pdTriangle[1] - pdTriangle[4])*(pdTriangle[1] - pdTriangle[4]) 
		            +(pdTriangle[2] - pdTriangle[5])*(pdTriangle[2] - pdTriangle[5]);
	double dDist13 = (pdTriangle[0] - pdTriangle[6])*(pdTriangle[0] - pdTriangle[6])
		            +(pdTriangle[1] - pdTriangle[7])*(pdTriangle[1] - pdTriangle[7]) 
		            +(pdTriangle[2] - pdTriangle[8])*(pdTriangle[2] - pdTriangle[8]);
	double dDist23 = (pdTriangle[3] - pdTriangle[6])*(pdTriangle[3] - pdTriangle[6])
		            +(pdTriangle[4] - pdTriangle[7])*(pdTriangle[4] - pdTriangle[7]) 
		            +(pdTriangle[5] - pdTriangle[8])*(pdTriangle[5] - pdTriangle[8]);
	double dSqrD = dMinDistance*dMinDistance;
	if ((dDist12<dSqrD)||(dDist13<dSqrD)||(dDist23<dSqrD))
		return false ;
	return true ;
}

/*************************************************************************
Function: Calculate square distance of a point to a triangle by
input three points coordinate.
return : double as square distance.
input  : double* pdPt, the point coordinate are a array of double[3]
         respectively as x ,y and z.
		 double* pdTriangle, three point coordinate of the triangle,
		 are a array of double[9] respectively as x1 ,y1, z1,
		 x2 ,y2, z2, x3 ,y3, z3.
output : double* pdS point S parameter of matched point in the triangle.
         double* pdT point T parameter of matched point in the triangle.

note   : if the triangle become a line or a point, the T parameter is
         negative, and by S parameter we can calculate the matched point.
**************************************************************************/
double CStdMath::gCalcSqrDistPtToTriangle(double* pdPt, double* pdTriangle, double *pdS, double* pdT)
{
	double dDiffx,dDiffy,dDiffz ;
	double dEdge0x = (pdTriangle[3] - pdTriangle[0]) ;
	double dEdge0y = (pdTriangle[4] - pdTriangle[1]) ;
	double dEdge0z = (pdTriangle[5] - pdTriangle[2]) ;
	double dEdge1x = (pdTriangle[6] - pdTriangle[0]) ;
	double dEdge1y = (pdTriangle[7] - pdTriangle[1]) ;
	double dEdge1z = (pdTriangle[8] - pdTriangle[2]) ;

	dDiffx = pdTriangle[0] - pdPt[0] ;
	dDiffy = pdTriangle[1] - pdPt[1] ;
	dDiffz = pdTriangle[2] - pdPt[2] ;

	double dA00 = dEdge0x*dEdge0x + dEdge0y*dEdge0y + dEdge0z*dEdge0z ;
	double dA01 = dEdge0x*dEdge1x + dEdge0y*dEdge1y + dEdge0z*dEdge1z ;
	double dA11 = dEdge1x*dEdge1x + dEdge1y*dEdge1y + dEdge1z*dEdge1z ;
    
	double dB0 = dDiffx*dEdge0x + dDiffy*dEdge0y + dDiffz*dEdge0z ;
	double dB1 = dDiffx*dEdge1x + dDiffy*dEdge1y + dDiffz*dEdge1z ;

	double dC  = dDiffx*dDiffx + dDiffy*dDiffy + dDiffz*dDiffz ;
	double dDet = dA00*dA11 - dA01*dA01 ;
	
	if (dDet < 0)
		dDet = -dDet ;
    	
	//calculate the distance of point to line.
	if (dDet < EPSLON)                  //this triangle has become a line or a point. 
	{
		double dS0,dS1,dS2;
		double dDistE0 = 0.0;
		double dDistE1 = 0.0;
		double dDistE2 = 0.0;
		double dA22 = (pdTriangle[6] - pdTriangle[3])*(pdTriangle[6] - pdTriangle[3]) 
			         +(pdTriangle[7] - pdTriangle[4])*(pdTriangle[7] - pdTriangle[4]) 
					 +(pdTriangle[8] - pdTriangle[5])*(pdTriangle[8] - pdTriangle[5]) ;
		//calculate distance to first edge.
		if (dA00 < EPSLON)
		{
			dDistE0 = dC   ;
			dS0     = 0.0  ;
		}
		else
		{
            dS0  =  -dB0 / dA00 ;
			if (dS0 <= 0.0)
			{
				dDistE0 = dC   ;
				dS0     = 0.0  ;
			}
			else if (dS0 >= 1.0)
			{
				dDistE0 = (pdTriangle[3] - pdPt[0])*(pdTriangle[3] - pdPt[0])
					     +(pdTriangle[4] - pdPt[1])*(pdTriangle[4] - pdPt[1])
						 +(pdTriangle[5] - pdPt[2])*(pdTriangle[5] - pdPt[2]) ;
				dS0     = 1.0  ;
			}
			else
			{
				dDistE0 = dC - dB0*dB0/dA00 ;
			}
		}

		//calculate distance to second edge.
		if (dA11 < EPSLON)
		{
			dDistE1 = dC ;
		}
		else
		{
			dS1  =  -dB1 / dA11 ;
			if (dS1 <= 0.0)
			{
				dDistE1 = dC   ;
				dS1     = 0.0  ;
			}
			else if (dS1 >= 1.0)
			{
				dDistE1 = (pdTriangle[6] - pdPt[0])*(pdTriangle[6] - pdPt[0])
					     +(pdTriangle[7] - pdPt[1])*(pdTriangle[7] - pdPt[1])
						 +(pdTriangle[8] - pdPt[2])*(pdTriangle[8] - pdPt[2]) ;
				dS1     = 1.0  ;
			}
			else
			{
				dDistE1 = dC - dB1*dB1/dA11 ;
			}
		}

        //calculate distance to third edge.
		dDiffx = pdTriangle[3] - pdPt[0] ;
	    dDiffy = pdTriangle[4] - pdPt[1] ;
	    dDiffz = pdTriangle[5] - pdPt[2] ;
		dC = dDiffx*dDiffx + dDiffy*dDiffy + dDiffz*dDiffz ;
		if (dA22 < EPSLON)
		{
			dDistE2 = dC ;
			dS2     = 0.0;
		}
		else
		{
			double dB2 = dDiffx*(pdTriangle[6] - pdTriangle[3])
				       + dDiffy*(pdTriangle[7] - pdTriangle[4])
					   + dDiffz*(pdTriangle[8] - pdTriangle[5]) ;
			dS2  =  -dB2 / dA22 ;
			if (dS2 <= 0.0)
			{
				dDistE2 = dC   ;
				dS2     = 0.0  ;
			}
			else if (dS2 >= 1.0)
			{
				dDistE2 = (pdTriangle[6] - pdPt[0])*(pdTriangle[6] - pdPt[0])
					     +(pdTriangle[7] - pdPt[1])*(pdTriangle[7] - pdPt[1])
						 +(pdTriangle[8] - pdPt[2])*(pdTriangle[8] - pdPt[2]) ;
				dS2     = 1.0  ;
			}
			else
			{
				dDistE2 = dC - dB2*dB2/dA22 ;
			}
		}

		if (pdT != NULL)
		{
			*pdT = -1.0  ;
		}

		if (dDistE0 > dDistE1)
		{
			dDistE0 = dDistE1 ;
			dS0     = dS1     ;
			if (pdT != NULL)
			{
				*pdT = -2.0  ;
			}
		}
		if (dDistE0 > dDistE2)
		{
			dDistE0 = dDistE2 ;
			dS0     = dS2     ;
			if (pdT != NULL)
			{
				*pdT = -3.0  ;
			}
		}
		if (pdS != NULL)
		{
		    *pdS  = dS0 ;
		}
		return dDistE0;
	}  //end if (dDet < EPSLON) . When this triangle has become a line or a point.

	double dS = dA01*dB1 - dA11*dB0;
    double dT = dA01*dB0 - dA00*dB1;

    double dSqrDist = 0.0;

    if ( dS + dT <= dDet )
    {
        if ( dS < 0.0 )
        {
            if ( dT < 0.0 )  // region 4
            {
                if ( dB0 < 0.0 )
                {
                    dT = 0.0;
                    if ( -dB0 >= dA00 )
                    {
                        dS = 1.0;
                        dSqrDist = dA00 + 2.0*dB0 + dC;
                    }
                    else
                    {
                        dS = -dB0/dA00;
                        dSqrDist = dB0*dS + dC;
                    }
                }
                else
                {
                    dS = 0.0;
                    if ( dB1 >= 0.0 )
                    {
                        dT = 0.0;
                        dSqrDist = dC;
                    }
                    else if ( -dB1 >= dA11 )
                    {
                        dT = 1.0;
                        dSqrDist = dA11+ 2.0*dB1 + dC;
                    }
                    else
                    {
                        dT = -dB1/dA11;
                        dSqrDist = dB1*dT + dC;
                    }
                }
            }
            else  // region 3
            {
                dS = 0.0;
                if ( dB1 >= 0.0 )
                {
                    dT = 0.0;
                    dSqrDist = dC;
                }
                else if ( -dB1 >= dA11 )
                {
                    dT = 1.0;
                    dSqrDist = dA11 + 2.0*dB1 + dC;
                }
                else
                {
                    dT = -dB1/dA11;
                    dSqrDist = dB1*dT + dC;
                }
            }
        }
        else if ( dT < 0.0 )  // region 5
        {
            dT = 0.0;
            if ( dB0 >= 0.0 )
            {
                dS = 0.0;
                dSqrDist = dC;
            }
            else if ( -dB0 >= dA00 )
            {
                dS = 1.0;
                dSqrDist = dA00 + 2.0*dB0 + dC;
            }
            else
            {
                dS = -dB0/dA00;
                dSqrDist = dB0*dS + dC;
            }
        }
        else  // region 0
        {
            // minimum at interior point
            double dInvDet = 1.0/dDet;
            dS *= dInvDet;
            dT *= dInvDet;
            dSqrDist = dS*(dA00*dS + dA01*dT + 2.0*dB0) + dT*(dA01*dS + dA11*dT + 2.0*dB1) + dC;
        }
    }
    else
    {
        double dTmp0, dTmp1, dNumer, dDenom;

        if ( dS < 0.0 )  // region 2
        {
            dTmp0 = dA01 + dB0;
            dTmp1 = dA11 + dB1;
            if ( dTmp1 > dTmp0 )
            {
                dNumer = dTmp1 - dTmp0;
                dDenom = dA00 - 2.0*dA01 + dA11;
                if ( dNumer >= dDenom )
                {
                    dS = 1.0;
                    dT = 0.0;
                    dSqrDist = dA00 + 2.0*dB0 + dC;
                }
                else
                {
                    dS = dNumer/dDenom;
                    dT = 1.0 - dS;
                    dSqrDist = dS*(dA00*dS + dA01*dT + 2.0*dB0) + dT*(dA01*dS + dA11*dT + 2.0*dB1) + dC;
                }
            }
            else
            {
                dS = 0.0;
                if ( dTmp1 <= 0.0 )
                {
                    dT = 1.0;
                    dSqrDist = dA11 + 2.0*dB1 + dC;
                }
                else if ( dB1 >= 0.0 )
                {
                    dT = 0.0;
                    dSqrDist = dC;
                }
                else
                {
                    dT = -dB1/dA11;
                    dSqrDist = dB1*dT + dC;
                }
            }
        }
        else if ( dT < 0.0 )  // region 6
        {
            dTmp0 = dA01 + dB1;
            dTmp1 = dA00 + dB0;
            if ( dTmp1 > dTmp0 )
            {
                dNumer = dTmp1 - dTmp0;
                dDenom = dA00 - 2.0*dA01 + dA11;
                if ( dNumer >= dDenom )
                {
                    dT = 1.0;
                    dS = 0.0;
                    dSqrDist = dA11 + 2.0*dB1 + dC;
                }
                else
                {
                    dT = dNumer/dDenom;
                    dS = 1.0 - dT;
                    dSqrDist = dS*(dA00*dS + dA01*dT + 2.0*dB0) + dT*(dA01*dS + dA11*dT + 2.0*dB1) + dC;
                }
            }
            else
            {
                dT = 0.0;
                if ( dTmp1 <= 0.0 )
                {
                    dS = 1.0;
                    dSqrDist = dA00 + 2.0*dB0 + dC;
                }
                else if ( dB0 >= 0.0 )
                {
                    dS = 0.0;
                    dSqrDist = dC;
                }
                else
                {
                    dS = -dB0/dA00;
                    dSqrDist = dB0*dS+dC;
                }
            }
        }
        else  // region 1
        {
            dNumer = dA11 + dB1 - dA01 - dB0;
            if ( dNumer <= 0.0 )
            {
                dS = 0.0;
                dT = 1.0;
                dSqrDist = dA11 + 2.0*dB1 + dC;
            }
            else
            {
                dDenom = dA00 - 2.0*dA01 + dA11;
                if ( dNumer >= dDenom )
                {
                    dS = 1.0;
                    dT = 0.0;
                    dSqrDist = dA00 + 2.0*dB0 + dC;
                }
                else
                {
                    dS = dNumer/dDenom;
                    dT = 1.0 - dS;
                    dSqrDist = dS*(dA00*dS + dA01*dT + 2.0*dB0) + dT*(dA01*dS + dA11*dT + 2.0*dB1) + dC;
                }
            }
        }
    }

    if ( pdS )
        *pdS = dS;

    if ( pdT )
        *pdT = dT;

    if (dSqrDist < 0)
		dSqrDist = -dSqrDist;

	return dSqrDist;
}

/*************************************************************************
Function: Calculate square distance of a point to a triangle by
input base point coordinate, edge0 and edge1 coordinate.
return : double as square distance.
input  : double* pdPt:
         the point coordinate are a array of double[3] respectively as x ,y 
		 and z.
		 double* pdTriangle:
		 base point coordinate of the triangle,and Edge0 ,Edge1 coordinate 
		 are a array of double[9] respectively as x1 ,y1, z1, x2-x1 ,y2-y1, 
		 z2-z1, x3-x1 ,y3-y1, z3-z1.
output : double* pdS :
         point S parameter of matched point in the triangle.
         double* pdT:
		 point T parameter of matched point in the triangle.
**************************************************************************/
double CStdMath::gCalcSqrDistPtToTriangle1(double* pdPt, double* pdTriangle, double *pdS, double* pdT)
{
	double dDiffx,dDiffy,dDiffz ;
	double dEdge0x = pdTriangle[3] ;
	double dEdge0y = pdTriangle[4] ;
	double dEdge0z = pdTriangle[5] ;
	double dEdge1x = pdTriangle[6] ;
	double dEdge1y = pdTriangle[7] ;
	double dEdge1z = pdTriangle[8] ;

	dDiffx = pdTriangle[0] - pdPt[0] ;
	dDiffy = pdTriangle[1] - pdPt[1] ;
	dDiffz = pdTriangle[2] - pdPt[2] ;

	double dA00 = dEdge0x*dEdge0x + dEdge0y*dEdge0y + dEdge0z*dEdge0z ;
	double dA01 = dEdge0x*dEdge1x + dEdge0y*dEdge1y + dEdge0z*dEdge1z ;
	double dA11 = dEdge1x*dEdge1x + dEdge1y*dEdge1y + dEdge1z*dEdge1z ;
    
	double dB0 = dDiffx*dEdge0x + dDiffy*dEdge0y + dDiffz*dEdge0z ;
	double dB1 = dDiffx*dEdge1x + dDiffy*dEdge1y + dDiffz*dEdge1z ;

	double dC  = dDiffx*dDiffx + dDiffy*dDiffy + dDiffz*dDiffz ;
	double dDet = dA00*dA11 - dA01*dA01 ;
	if (dDet < 0)
		dDet = -dDet ;

	double dS = dA01*dB1 - dA11*dB0;
    double dT = dA01*dB0 - dA00*dB1;

    double dSqrDist;

    if ( dS + dT <= dDet )
    {
        if ( dS < 0.0 )
        {
            if ( dT < 0.0 )  // region 4
            {
                if ( dB0 < 0.0 )
                {
                    dT = 0.0;
                    if ( -dB0 >= dA00 )
                    {
                        dS = 1.0;
                        dSqrDist = dA00 + 2.0*dB0 + dC;
                    }
                    else
                    {
                        dS = -dB0/dA00;
                        dSqrDist = dB0*dS + dC;
                    }
                }
                else
                {
                    dS = 0.0;
                    if ( dB1 >= 0.0 )
                    {
                        dT = 0.0;
                        dSqrDist = dC;
                    }
                    else if ( -dB1 >= dA11 )
                    {
                        dT = 1.0;
                        dSqrDist = dA11+ 2.0*dB1 + dC;
                    }
                    else
                    {
                        dT = -dB1/dA11;
                        dSqrDist = dB1*dT + dC;
                    }
                }
            }
            else  // region 3
            {
                dS = 0.0;
                if ( dB1 >= 0.0 )
                {
                    dT = 0.0;
                    dSqrDist = dC;
                }
                else if ( -dB1 >= dA11 )
                {
                    dT = 1.0;
                    dSqrDist = dA11 + 2.0*dB1 + dC;
                }
                else
                {
                    dT = -dB1/dA11;
                    dSqrDist = dB1*dT + dC;
                }
            }
        }
        else if ( dT < 0.0 )  // region 5
        {
            dT = 0.0;
            if ( dB0 >= 0.0 )
            {
                dS = 0.0;
                dSqrDist = dC;
            }
            else if ( -dB0 >= dA00 )
            {
                dS = 1.0;
                dSqrDist = dA00 + 2.0*dB0 + dC;
            }
            else
            {
                dS = -dB0/dA00;
                dSqrDist = dB0*dS + dC;
            }
        }
        else  // region 0
        {
            // minimum at interior point
            double dInvDet = 1.0/dDet;
            dS *= dInvDet;
            dT *= dInvDet;
            dSqrDist = dS*(dA00*dS + dA01*dT + 2.0*dB0) + dT*(dA01*dS + dA11*dT + 2.0*dB1) + dC;
        }
    }//end of if ( dS + dT <= dDet )
    else
    {
        double dTmp0, dTmp1, dNumer, dDenom;

        if ( dS < 0.0 )  // region 2
        {
            dTmp0 = dA01 + dB0;
            dTmp1 = dA11 + dB1;
            if ( dTmp1 > dTmp0 )
            {
                dNumer = dTmp1 - dTmp0;
                dDenom = dA00 - 2.0*dA01 + dA11;
                if ( dNumer >= dDenom )
                {
                    dS = 1.0;
                    dT = 0.0;
                    dSqrDist = dA00 + 2.0*dB0 + dC;
                }
                else
                {
                    dS = dNumer/dDenom;
                    dT = 1.0 - dS;
                    dSqrDist = dS*(dA00*dS + dA01*dT + 2.0*dB0) + dT*(dA01*dS + dA11*dT + 2.0*dB1) + dC;
                }
            }
            else
            {
                dS = 0.0;
                if ( dTmp1 <= 0.0 )
                {
                    dT = 1.0;
                    dSqrDist = dA11 + 2.0*dB1 + dC;
                }
                else if ( dB1 >= 0.0 )
                {
                    dT = 0.0;
                    dSqrDist = dC;
                }
                else
                {
                    dT = -dB1/dA11;
                    dSqrDist = dB1*dT + dC;
                }
            }
        }//end of  if ( dS < 0.0 )  // region 2
        else if ( dT < 0.0 )  // region 6
        {
            dTmp0 = dA01 + dB1;
            dTmp1 = dA00 + dB0;
            if ( dTmp1 > dTmp0 )
            {
                dNumer = dTmp1 - dTmp0;
                dDenom = dA00 - 2.0*dA01 + dA11;
                if ( dNumer >= dDenom )
                {
                    dT = 1.0;
                    dS = 0.0;
                    dSqrDist = dA11 + 2.0*dB1 + dC;
                }
                else
                {
                    dT = dNumer/dDenom;
                    dS = 1.0 - dT;
                    dSqrDist = dS*(dA00*dS + dA01*dT + 2.0*dB0) + dT*(dA01*dS + dA11*dT + 2.0*dB1) + dC;
                }
            }
            else
            {
                dT = 0.0;
                if ( dTmp1 <= 0.0 )
                {
                    dS = 1.0;
                    dSqrDist = dA00 + 2.0*dB0 + dC;
                }
                else if ( dB0 >= 0.0 )
                {
                    dS = 0.0;
                    dSqrDist = dC;
                }
                else
                {
                    dS = -dB0/dA00;
                    dSqrDist = dB0*dS+dC;
                }
            }
        } //end of else if ( dT < 0.0 )  // region 6
        else  // region 1
        {
            dNumer = dA11 + dB1 - dA01 - dB0;
            if ( dNumer <= 0.0 )
            {
                dS = 0.0;
                dT = 1.0;
                dSqrDist = dA11 + 2.0*dB1 + dC;
            }
            else
            {
                dDenom = dA00 - 2.0*dA01 + dA11;
                if ( dNumer >= dDenom )
                {
                    dS = 1.0;
                    dT = 0.0;
                    dSqrDist = dA00 + 2.0*dB0 + dC;
                }
                else
                {
                    dS = dNumer/dDenom;
                    dT = 1.0 - dS;
                    dSqrDist = dS*(dA00*dS + dA01*dT + 2.0*dB0) + dT*(dA01*dS + dA11*dT + 2.0*dB1) + dC;
                }
            }
        }//end of else  // region 1
    }//end of else //if ( dS + dT <= dDet )

    if ( pdS )
        *pdS = dS;

    if ( pdT )
        *pdT = dT;

    if (dSqrDist < 0)
		dSqrDist = -dSqrDist;

	return dSqrDist;
}

/**********************************************************************
Function: Calculate distance of a point to a triangle.
return  : double as distance.
input   : double* pdPt:
          the point coordinate are a array of double[3] respectively as
		  x ,y and z.
		  double* pdTriangle:
		  three point coordinate of the triangle,are a array of double[9]
		  respectively as x1 ,y1, z1,x2 ,y2, z2, x3 ,y3, z3.
***********************************************************************/
double CStdMath::gCalcDistPtToTriangle(double* pdPt, double* pdTriangle)
{
	double dDist = 0;
	dDist = gCalcSqrDistPtToTriangle(pdPt, pdTriangle, NULL, NULL);
	dDist = sqrt(dDist);
	return dDist;
}

/***************************************************************************
Function: Calculate distance of a point to a triangle and out the 
          matched point's coordinate in the triangle.
return  : double as distance.
input   : double* pdPt: 
          the point coordinate are a array of double[3]
          respectively as x ,y and z.
		  double* pdTriangle:
		  three point coordinate of the triangle,
		  are a array of double[9] respectively as x1 ,y1, z1,
		  x2 ,y2, z2, x3 ,y3, z3.
output  : double* pdMp:
          point coordinate of the matched point in the triangle,
          are a array of double[3].
*****************************************************************************/
double CStdMath::gCalcDistPtToTriangleWithMp(double* pdPt, double* pdTriangle, double *pdMp)
{
	double dDist = 0;
	double dS,dT ;
	dDist = gCalcSqrDistPtToTriangle(pdPt, pdTriangle, &dS, &dT);
	pdMp[0] = pdTriangle[0] + dS*(pdTriangle[3] - pdTriangle[0]) 
		    + dT*(pdTriangle[6] - pdTriangle[0]) ;
	pdMp[1] = pdTriangle[1] + dS*(pdTriangle[4] - pdTriangle[1]) 
		    + dT*(pdTriangle[7] - pdTriangle[1]) ;
	pdMp[2] = pdTriangle[2] + dS*(pdTriangle[5] - pdTriangle[2]) 
		    + dT*(pdTriangle[8] - pdTriangle[2]) ;
	dDist = sqrt(dDist);
	return dDist;
}

/*****************************************************************************
Function: Calculate max-edge length of all input triangle.
return  : double . the max length of all triangle' edges.
input   : int iTriNum:
          Triangle's number.
		  double* pdTriangle:
		  point to data of triangles vertex coordinate.
*****************************************************************************/
double CStdMath::gMaxLengthOfAllTriangle(int iTriNum, double* pdTriangle,double *pdAveLen)
{
	int i ;
	double dDist ,dAveD = 0.0;
	double dMaxLen = -1.0e305 ;
	double* pdCur = pdTriangle ;
	for (int i=0; i<iTriNum; i++)
	{
		dDist = (pdCur[3] - pdCur[0])*(pdCur[3] - pdCur[0])
			   +(pdCur[4] - pdCur[1])*(pdCur[4] - pdCur[1])
			   +(pdCur[5] - pdCur[2])*(pdCur[5] - pdCur[2]) ;
		dAveD += sqrt(dDist) ;
		if (dDist > dMaxLen)
		{
			dMaxLen = dDist ;
		}

		dDist = (pdCur[6] - pdCur[0])*(pdCur[6] - pdCur[0])
			   +(pdCur[7] - pdCur[1])*(pdCur[7] - pdCur[1])
			   +(pdCur[8] - pdCur[2])*(pdCur[8] - pdCur[2]) ;
		dAveD += sqrt(dDist) ;
		if (dDist > dMaxLen)
		{
			dMaxLen = dDist ;
		}

		dDist = (pdCur[3] - pdCur[6])*(pdCur[3] - pdCur[6])
			   +(pdCur[4] - pdCur[7])*(pdCur[4] - pdCur[7])
			   +(pdCur[5] - pdCur[8])*(pdCur[5] - pdCur[8]) ;
		dAveD += sqrt(dDist) ;
		if (dDist > dMaxLen)
		{
			dMaxLen = dDist ;
		}

		pdCur += 9;
	}

	dAveD /= 3.0*iTriNum ;
	*pdAveLen = dAveD ;
	return sqrt(dMaxLen) ;
}

/*****************************************************************************
Function: Calculate a point is located at which side of a triangle.
return  : double. The dotproduct of the vector which start from the base point 
          of the triangle to the point and  the vector of the triangle's normal.
		  if >0, the point is at the same side,otherwise at another side.
input   : double *pdPt:
          the point coordinate are a array of double[3] respectively as x ,y 
		  and z.
          double* pdTriangle:
		  three vertex coordinate of the triangle,are a array of double[9]
		  respectively as x1 ,y1, z1,x2 ,y2, z2, x3 ,y3, z3.
*****************************************************************************/
double CStdMath::gPtWithTriangle(double *pdPt,double* pdTriangle)
{
	double daE0[3],daE1[3],daEp[3], daCrossProd[3],dDotProd;
	daE0[0] = pdTriangle[3] - pdTriangle[0];
	daE0[1] = pdTriangle[4] - pdTriangle[1];
	daE0[2] = pdTriangle[5] - pdTriangle[2];
	daE1[0] = pdTriangle[6] - pdTriangle[0];
	daE1[1] = pdTriangle[7] - pdTriangle[1];
	daE1[2] = pdTriangle[8] - pdTriangle[2];

	daEp[0] = pdPt[0] - pdTriangle[0];
	daEp[1] = pdPt[1] - pdTriangle[1];
	daEp[2] = pdPt[2] - pdTriangle[2];

    daCrossProd[0] = daE0[1]*daE1[2] - daE0[2]*daE1[1] ;
	daCrossProd[1] = daE0[2]*daE1[0] - daE0[0]*daE1[2] ;
	daCrossProd[2] = daE0[0]*daE1[1] - daE0[1]*daE1[0] ;
   
    dDotProd = daCrossProd[0]*daEp[0] + daCrossProd[1]*daEp[1] + daCrossProd[2]*daEp[2] ;

	return dDotProd ;
}

/*****************************************************************************
Function: Calculate the errors of all input point to the surface 
which is composed of all input triangles.
return : bool.if memory is not enough ,return false.
input  : int iPtNum: 
         points number.
         int iTriNum: 
		 triangles number.
		 double* pdTriTab:
		 point to data of triangles vertex coordinate.
         double dMaxErr: 
		 Max Error of all point.
		 int iHitValue: a parameter to determine the catch's size.
input and output:
		 double *pdPtTab:
		 point to the point's data, and return error values.
*****************************************************************************/
bool CStdMath::gCalcError(int iPtNum, int iTriNum, double *pdPtTab,double* pdTriTab,
				double dMaxErr,double dHitValue)
{
	int i ,j, iPos ,iCatchNum ,iCircleNum ,iHitCount;
	double dSqrDist, dErr ,daHitPt[3],daPt[3],dMaxLen, dAbsSum;
	double dMaxCalcDist ,dMaxHitLen, dMaxHitLen1 ;
	double *pCur;
	double *pdTriTab1;
	bool   bHit = NULL ;
    
	//Get the catch.
	pdTriTab1 = new double[9*MAX_TRINUM];
	if (!pdTriTab1)
		return false;

	//calculating max edge length in all triangle.
	double dAveLen ;
    dMaxLen =  gMaxLengthOfAllTriangle(iTriNum, pdTriTab, &dAveLen);

	dMaxCalcDist = dMaxErr + dMaxLen ;           //if the x,y,z absolute distances exceed this value
	                                             //can not be calculate .
	dMaxHitLen   = dMaxErr + dMaxLen*dHitValue ; //the criterion to find catched triangles.
	if (dHitValue > 5.0)
        dMaxHitLen1  = dMaxLen*(dHitValue - 2.0) ;     //the criterion to determine if the point hit the 
	                                                   //catched triangles.
	else
		dMaxHitLen1  = dMaxHitLen*3/5.0 ;
	
	//Initial hit catch.
	iHitCount = 0 ;
	dErr = 1.0e305 ;
	iPos = 0;
	pCur = pdTriTab;
	iCircleNum = iTriNum ;
	daHitPt[0] = pdPtTab[0];
	daHitPt[1] = pdPtTab[1];
	daHitPt[2] = pdPtTab[2];	
	for (int i=0; i<iCircleNum; i++)
	{
		dAbsSum = fabs(pCur[0]-daHitPt[0])+fabs(pCur[1]-daHitPt[1])
			     +fabs(pCur[2]-daHitPt[2]) ;
		if (dAbsSum < dMaxHitLen)
		{
			pdTriTab1[iHitCount++] = pCur[0] ;
			pdTriTab1[iHitCount++] = pCur[1] ;
			pdTriTab1[iHitCount++] = pCur[2] ;
			pdTriTab1[iHitCount++] = pCur[3] ;
			pdTriTab1[iHitCount++] = pCur[4] ;
			pdTriTab1[iHitCount++] = pCur[5] ;
			pdTriTab1[iHitCount++] = pCur[6] ;
			pdTriTab1[iHitCount++] = pCur[7] ;
			pdTriTab1[iHitCount++] = pCur[8] ;
			if (dAbsSum < dMaxCalcDist)
			{
				dSqrDist = gCalcSqrDistPtToTriangle(daHitPt, pCur, NULL, NULL);
				
				if ((dSqrDist < dErr)&&(dSqrDist >= 0))
				{
					dErr = dSqrDist ; 
					iPos = i;
				}
			}
		}
		pCur += 9 ;
	}
    iCatchNum = iHitCount/9 ;
	pdPtTab[0] = sqrt(dErr);
    pCur = pdTriTab + 9*iPos;
	if (gPtWithTriangle(daHitPt,pCur) < 0)
	{
		pdPtTab[0] = -pdPtTab[0] ;
	}
	
	//calculating errors
	for (j=1; j<iPtNum; j++)
	{
		daPt[0] = pdPtTab[3*j] ;
		daPt[1] = pdPtTab[3*j+1] ;
		daPt[2] = pdPtTab[3*j+2] ;

		dAbsSum = fabs(daHitPt[0]-daPt[0])+fabs(daHitPt[1]-daPt[1])+fabs(daHitPt[2]-daPt[2]) ;
		if (dAbsSum < dMaxHitLen1)
		{
			bHit = true ;
			pCur = pdTriTab1;
			iCircleNum = iCatchNum;
		}
		else
		{
			bHit = false ;
			pCur = pdTriTab;
			iCircleNum = iTriNum ;
		}
		
		dErr = 1.0e305 ;
		iPos = 0 ;
		if (bHit)
		{
			for (int i=0; i<iCircleNum; i++)
			{
				if (fabs(pCur[0]-daPt[0])+fabs(pCur[1]-daPt[1])+fabs(pCur[2]-daPt[2])
					< dMaxCalcDist)
				{
					dSqrDist = gCalcSqrDistPtToTriangle(daPt, pCur, NULL, NULL);
					if ((dSqrDist < dErr)&&(dSqrDist >= 0))
					{
						dErr = dSqrDist ; 
						iPos = i;
					}
				}
				pCur += 9 ;
			}//end of for (i=0; i<iCircleNum; i++)
			pdPtTab[j] = sqrt(dErr);
			pCur = pdTriTab1 + 9*iPos;
			if (gPtWithTriangle(daPt,pCur) < 0)
			{
				pdPtTab[j] = -pdPtTab[j] ;
			}
		}   //end of if (bHit)
		else
		{
			daHitPt[0] = daPt[0];
			daHitPt[1] = daPt[1];
			daHitPt[2] = daPt[2];
			iHitCount = 0;
			for (int i=0; i<iCircleNum; i++)
			{
				dAbsSum = fabs(pCur[0]-daHitPt[0])+fabs(pCur[1]-daHitPt[1])
						 +fabs(pCur[2]-daHitPt[2]) ;
				if (dAbsSum < dMaxHitLen)
				{
					pdTriTab1[iHitCount++] = pCur[0] ;
					pdTriTab1[iHitCount++] = pCur[1] ;
					pdTriTab1[iHitCount++] = pCur[2] ;
					pdTriTab1[iHitCount++] = pCur[3] ;
					pdTriTab1[iHitCount++] = pCur[4] ;
					pdTriTab1[iHitCount++] = pCur[5] ;
					pdTriTab1[iHitCount++] = pCur[6] ;
					pdTriTab1[iHitCount++] = pCur[7] ;
					pdTriTab1[iHitCount++] = pCur[8] ;
					if (dAbsSum< dMaxCalcDist)
					{
						dSqrDist = gCalcSqrDistPtToTriangle(daPt, pCur, NULL, NULL);
						if ((dSqrDist < dErr)&&(dSqrDist >= 0))
						{
							dErr = dSqrDist ; 
							iPos = i;
						}
					}
				}
				pCur += 9 ;
			}//end of for (i=0; i<iCircleNum; i++)
			iCatchNum = iHitCount/9 ;
			pdPtTab[j] = sqrt(dErr);
			pCur = pdTriTab + 9*iPos;
			if (gPtWithTriangle(daPt,pCur) < 0)
			{
				pdPtTab[j] = -pdPtTab[j] ;
			}
		} // end of else  //if (bHit)
    }//end of for (j=1; j<iPtNum; j++)

	delete[] pdTriTab1 ;
	return true ;
}
namespace stdMath {

unsigned char GetR(unsigned int v)
{
    return v>>16;
}
unsigned char GetG(unsigned int  v)
{
    v = v<<16;
    v = v>> 24;
    return v;
}
unsigned char GetB(unsigned int v)
{
    v = v<<24;
    v = v>> 24;
    return v;
}
#ifdef WIN32
#else
    unsigned char GetRValue(unsigned int v)
    {
        return GetR(v);
    }
    unsigned char GetGValue(unsigned int  v)
    {
        return GetG(v);
    }
    unsigned char GetBValue(unsigned int v)
    {
        return GetB(v);
    }
#endif
    unsigned int RGBA(unsigned char r,unsigned char g,unsigned char b)
{
    unsigned int cc=r;
    cc = cc<<8;
    cc= cc + g;
    cc = cc<<8;
    cc= cc + b;
    
    return cc;
}

}


