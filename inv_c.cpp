#include "mex.h"
#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <limits>  
using namespace std;


#define A_IN prhs[0]
#define B_OUT plhs[0]

void DoComputation(double *B, double *A, int M, int N, double p)
{
    double colnorm;
    int m, n;
    for(n = 0; n < N; n++)
    {
    /* Compute the norm of the nth column */
        for(m = 0, colnorm = 0.0; m < M; m++)
        colnorm += pow(A[m + M*n], p);
        colnorm = pow(fabs(colnorm), 1.0/p);
        /* Fill the nth column of B */
        for(m = 0; m < M; m++)
            B[m + M*n] = A[m + M*n]/colnorm;
    }
}


void myprint(double *A, size_t n) {
    mexPrintf("matrix\n");
    char buffer[20];
    for (size_t i = 0; i < n; i++) {
        size_t k=sprintf (buffer, "%f |", A[i]);
        mexPrintf(buffer, k);
    }
    mexPrintf("\n");
    return;
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs,  const mxArray *prhs[])
{
    size_t nrows, ncols;                   /* size of matrices */
    if(nrhs!=1) { 
        mexErrMsgIdAndTxt("MyToolbox:inv_c:nrhs","Three inputs required.");
    }
    if(nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:inv_c:nlhs","Must have 1 output arguments.");
    }
    nrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
   
    if(ncols!=nrows) { 
        mexErrMsgIdAndTxt("MyToolbox:inv_c:nrhs","Matrix must be square.");
    }
    
    double * A1 = (double *)mxGetPr(A_IN);
    vector<double> A;
     
    
    //char buffer [100];
    //size_t k=sprintf (buffer, "%d  = ncols, %d  = nrows, %d. A[1] = %f\n", ncols, nrows, ncols*nrows,   *A1);
    //mexPrintf(buffer, k);

    for (size_t i = 0; i < ncols*nrows; i++) {
        //char buffer1 [100];
        //size_t k=sprintf (buffer1, "%f  = A[%d],\n", *A1, i);
        //mexPrintf(buffer1, k);
        A.push_back(double(*A1));
        A1++;
    }
    
    B_OUT = mxCreateDoubleMatrix((mwSize) nrows, (mwSize) ncols, mxREAL);
    //ncols = mxGetN(A_IN);
    
   
    //The elements of A are stored contiguously in memory in column-major format,
    //A[m + M*n] corresponds to A(m+1,n+1)
    double * E = mxGetPr(B_OUT);
    for(size_t n = 0; n < ncols; n++)
    {
        E[n*ncols+n] = 1;
    }
    //myprint(E, ncols*ncols);

    for(size_t n = 0; n < ncols; n++)
    {
        double tmp1 = A[n*ncols+n];
        if (abs(tmp1) < numeric_limits<double>::epsilon()) {
            mexWarnMsgIdAndTxt("MyToolbox:inv_c","Bad matrix");
        }
        for (size_t k = 0; k < ncols; k++) {
            A[n*ncols+k] = A[n*ncols+k]/tmp1;
            E[n*ncols+k] = E[n*ncols+k]/tmp1;
        }

        for (size_t k = n+1; k < ncols; k++) {
            double tmp2 = A[k*ncols+n];
            for (size_t k1 = 0; k1 < ncols; k1++) {
                A[k*ncols+k1] = A[k*ncols+k1] - A[n*ncols+k1]*tmp2;
                E[k*ncols+k1] = E[k*ncols+k1] - E[n*ncols+k1]*tmp2;
            }
        }
    }
    
    //myprint(E, ncols*ncols);
    
    for(int n = ncols - 1; n >= 0; n--)
    {
        for(int i = n-1; i >= 0; i--)
        { 
            double tmp3 = A[i*ncols+n];
            //char buffer1 [100];
            //size_t k=sprintf (buffer1, "%d = i, %d = n, tmp3 = %f\n", i, n, tmp3);
            //mexPrintf(buffer1, k);
            for(size_t n1 = 0; n1 < ncols; n1++)
            {
                A[i*ncols+n1] = A[i*ncols+n1] - A[n*ncols+n1]*tmp3;
                E[i*ncols+n1] = E[i*ncols+n1] - E[n*ncols+n1]*tmp3;
            }
        }
    }

return;
}
