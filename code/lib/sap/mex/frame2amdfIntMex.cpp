// Frame 2 AMDF, integer version
// How to compile:
// MATLAB 7.1: mex frame2amdfIntMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp -output frame2amdfIntMex.dll
// Others: mex frame2amdfIntMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp
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
	int *intFrame, *intAmdf, i, frameSize, maxShift, method;

	/* Check for proper number of arguments */
	if (nrhs<1) {
		char message[200];
		strcpy(message, mexFunctionName());
		strcat(message, " requires 1 input arguments.\n");
		strcat(message, "Usage: amdf = ");
		strcat(message, mexFunctionName());
		strcat(message, "(frame, method)\n");
		strcat(message, "\tThis function is synchronized with frame2amdf.m.\n");
		strcat(message, "\tNote that the input frame must be a vector of integers!\n");
		mexErrMsgTxt(message);
	}
	
	/* Assign pointers to the various parameters */
	frame = mxGetPr(FRAME);
	frameSize = mxGetM(FRAME)*mxGetN(FRAME);
	if (vecMax(frame, frameSize)<=1)
		mexPrintf("Input frame signals are all <= 1! Perhaps you should use wavReadInt.m instead?");
	maxShift=frameSize;
	method=1;
	if (nrhs>=2)
		maxShift=mxGetScalar(MAXSHIFT);
	if (nrhs>=3)
		method=mxGetScalar(METHOD);

	intFrame = (int *)malloc(frameSize*sizeof(int));
	double2intVec(frame, frameSize, intFrame);
	intAmdf = (int *)malloc(sizeof(int)*maxShift);

	/* Do the actual computation */
	frame2amdf(intFrame, frameSize, intAmdf, maxShift, method);

	/* Create a matrix for the return argument */
	AMDF = mxCreateDoubleMatrix(maxShift, 1, mxREAL);
	amdf = mxGetPr(AMDF);
	int2doubleVec(intAmdf, maxShift, amdf);

	free(intFrame);
	free(intAmdf);
}
