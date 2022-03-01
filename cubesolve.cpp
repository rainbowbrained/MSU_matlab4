#include "mex.h"
#include <iostream>
#include <complex>
#include <vector>
#include <limits>  
using namespace std;
void cubesolve(vector<complex<double>>& A, vector<complex<double>>& B, vector<complex<double>>& C, double *out_x1_real, double *out_x1_imag, double *out_x2_real, double *out_x2_imag, double *out_x3_real, double *out_x3_imag)
{
    for (size_t i = 0; i < A.size(); i++)
        {
        const complex<double> a = A[i], b = B[i], c = C[i];
        complex<double> x1, x2, x3, D, cube1, cube2, y1, y2, y3;
        double r1, r2, phi1, phi2, epsilon;
        const double pi = 3.14115926535;
        if (abs(a) > epsilon)
        {
            complex<double> p = b/a;
            complex<double> q = c/a;
            D = pow(q/double(2), 2) - pow(p/double(3), 3);
            cube1 = -q/double(2)+sqrt(D);
            cube2 = -q/double(2)-sqrt(D);
            r1 = pow(abs(cube1), double(1)/3);
            phi1 = arg(cube1)/3;
            y1 = polar(r1, phi1);
            y2 = polar(r1, phi1+2*pi/3);
            y3 = polar(r1, phi1+4*pi/3);
            x1 = y1 - p/double(3)/y1;
            x2 = y2 - p/double(3)/y2;
            x3 = y3 - p/double(3)/y3;
        }
        else if (abs(b) > epsilon)
        {
            x1 = -c/b;
            x2 = -c/b;
            x3 = -c/b;
        }
        else
        {
            x1 = numeric_limits<double>::quiet_NaN();
            x2 = numeric_limits<double>::quiet_NaN();
            x3 = numeric_limits<double>::quiet_NaN();
        }

        out_x1_real[i] = x1.real();
        out_x1_imag[i] = x1.imag();
        out_x2_real[i] = x2.real();
        out_x2_imag[i] = x2.imag();
        out_x3_real[i] = x3.real();
        out_x3_imag[i] = x3.imag();
        }
}
/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs,  const mxArray *prhs[])
{
    size_t nrows, ncols;                   /* size of matrices */
    /* input matrices */
    vector<complex<double>> A, B, C;
    /* output matrices */
    double *out_x1_real, *out_x1_imag, *out_x2_real, *out_x2_imag, *out_x3_real, *out_x3_imag;
    /* check for proper number of arguments */
    if(nrhs != 3) {
        mexErrMsgIdAndTxt("MyToolbox:cubesolve:nrhs","Three inputs required.");
    }
    if(nlhs < 2) {
        mexErrMsgIdAndTxt("MyToolbox:cubesolve:nlhs","At least two outputs required.");
    }
    if (nlhs > 3){
        mexErrMsgIdAndTxt("MyToolbox:cubesolve:nlhs","Too many outputs.");
    }
    if( (mxGetM(prhs[0])!= mxGetM(prhs[1])) || (mxGetM(prhs[1])!= mxGetM(prhs[2]))
    ||(mxGetN(prhs[0])!= mxGetN(prhs[1])) || (mxGetN(prhs[1])!= mxGetN(prhs[2])) ) {
        mexErrMsgIdAndTxt("MyToolbox:cubesolve:nrhs","Input matrices must be the same size.");
    }
    
    /* get dimensions of the input matrix */
    nrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    
    double *A_real = (double *)mxGetPr(prhs[0]);
    double *A_imag = (double *)mxGetPi(prhs[0]);
    double *B_real = (double *)mxGetPr(prhs[1]);
    double *B_imag = (double *)mxGetPi(prhs[1]);
    double *C_real = (double *)mxGetPr(prhs[2]);
    double *C_imag = (double *)mxGetPi(prhs[2]);
    for (size_t i = 0; i < nrows; i++)
      for (size_t j = 0; j < ncols; j++)
        {
            A.push_back(complex<double>(*A_real, *A_imag));
            B.push_back(complex<double>(*B_real, *B_imag));
            C.push_back(complex<double>(*C_real, *C_imag));
            A_real++;
            A_imag++;
            B_real++;
            B_imag++;
            C_real++;
            C_imag++;
        }
    
    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix((mwSize) nrows, (mwSize) ncols, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix((mwSize) nrows, (mwSize) ncols, mxCOMPLEX);
    plhs[2] = mxCreateDoubleMatrix((mwSize) nrows, (mwSize) ncols, mxCOMPLEX);
    out_x1_real = mxGetPr(plhs[0]);
    out_x1_imag = mxGetPi(plhs[0]);
    out_x2_real = mxGetPr(plhs[1]);
    out_x2_imag = mxGetPi(plhs[1]);
    out_x3_real = mxGetPr(plhs[2]);
    out_x3_imag = mxGetPi(plhs[2]);
    
    /* call the computational routine */
    cubesolve(A, B, C, out_x1_real, out_x1_imag, out_x2_real, out_x2_imag, out_x3_real, out_x3_imag);
}