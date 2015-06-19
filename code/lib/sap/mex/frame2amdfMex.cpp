// Frame 2 AMDF, double version
// How to compile:
// MATLAB 7.1: mex frame2amdfMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp -output frame2amdfMex.dll
// Others: mex frame2amdfMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp
#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	FRAME		prhs[0]
#define	MAXSHIFT	prhs[1]
#define	METHOD		prhs[2]

/* Output Arguments */
#define	AMDF	plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *frame, *amdf;
	int i, frameSize, maxShift, method;

	/* Check for proper number of arguments */
	if (nrhs<1) {
		char message[200];
		strcpy(message, mexFunctionName());
		strcat(message, " requires 1 input arguments.\n");
		strcat(message, "Usage: amdf = ");
		strcat(message, mexFunctionName());
		strcat(message, "(frame, maxShift, method)\n");
		strcat(message, "\tThis function is synchronized with frame2amdf.m.\n");
		mexErrMsgTxt(message);
	}
	
	/* Assign pointers to the various parameters */
	frame = mxGetPr(FRAME);
	frameSize = mxGetM(FRAME)*mxGetN(FRAME);
	maxShift=frameSize;
	method=1;
	if (nrhs>=2)
		maxShift=mxGetScalar(MAXSHIFT);
	if (nrhs>=3)
		method=mxGetScalar(METHOD);

	/* Create a matrix for the return argument */
	AMDF = mxCreateDoubleMatrix(maxShift, 1, mxREAL);
	amdf = mxGetPr(AMDF);

	/* Do the actual computation */
	frame2amdf(frame, frameSize, amdf, maxShift, method);
}
