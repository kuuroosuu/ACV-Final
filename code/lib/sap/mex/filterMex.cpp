// How to compile:
// mex  filterMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp
#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	B	prhs[0]
#define	A	prhs[1]
#define	X	prhs[2]
/* Output Arguments */
#define	OUT	plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double	*b, *a, *x, *out;
	int bLen, aLen, xLen, m, n;
	char message[200];
	cAudio myAudio;

	/* Check for proper number of arguments */
	if (nrhs!=3){
		sprintf(message, "Usage: out = %s(b, a, x)", mexFunctionName());
		mexErrMsgTxt(message);
	}

	/* Dimensions of the input matrix */
	m = mxGetM(X);
	n = mxGetN(X);

	/* Create a matrix for the return argument */
	OUT = mxCreateDoubleMatrix(m, n, mxREAL);

	/* Assign pointers to the various parameters */
	out = mxGetPr(OUT);
	b = mxGetPr(B); bLen=(int)mxGetM(B)*mxGetN(B);
	a = mxGetPr(A);	aLen=(int)mxGetM(A)*mxGetN(A);
	x = mxGetPr(X);

	/* Do the column-wise computations */
	for (int j=0; j<n; j++) {
		memcpy(out+j*m, x+j*m, m*sizeof(double));
		myAudio.filter(b, bLen, a, aLen, out+j*m, m);
	}
}
