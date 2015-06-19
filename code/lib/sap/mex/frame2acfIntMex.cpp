// Frame 2 ACF, integer version
// How to compile:
// MATLAB 7.1: mex frame2acfMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp -output frame2acfMex.dll
// Others: mex frame2acfMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp
#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	FRAME		prhs[0]
#define	MAXSHIFT	prhs[1]
#define	METHOD		prhs[2]
#define	USEMEANSUBTRACTION	prhs[2]

/* Output Arguments */
#define	ACF	plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *frame, *acf;
	int i, frameSize, maxShift, method, useMeanSubtraction;
	cAudio myAudio;

	/* Check for proper number of arguments */
	if (nrhs<1) {
		char message[200];
		strcpy(message, mexFunctionName());
		strcat(message, " requires 1 input arguments.\n");
		strcat(message, "Usage: acf = ");
		strcat(message, mexFunctionName());
		strcat(message, "(frame, maxShift, method, useMeanSubtraction)\n");
		strcat(message, "\tThis function is synchronized with frame2acf.m.\n");
		strcat(message, "\tNote that the input frame must be a vector of integers!\n");
		mexErrMsgTxt(message);
	}
	
	/* Assign pointers to the various parameters */
	frame = mxGetPr(FRAME);
	frameSize = mxGetM(FRAME)*mxGetN(FRAME);
	if (vecMax(frame, frameSize)<=1)
		mexPrintf("Input frame signals are all <=1 in frame2acfIntMex()! Perhaps you should use wavReadInt.m instead?\n");
	// Default parameters
	maxShift=frameSize;
	method=1;
	useMeanSubtraction=1;
	if (nrhs>=2)
		maxShift=mxGetScalar(MAXSHIFT);
	if (nrhs>=3)
		method=mxGetScalar(METHOD);
	if (nrhs>=4)
		useMeanSubtraction=mxGetScalar(USEMEANSUBTRACTION);

	int *frameInt = (int *)malloc(frameSize*sizeof(int));
	double2intVec(frame, frameSize, frameInt);

	/* Create a matrix for the return argument */
	ACF = mxCreateDoubleMatrix(maxShift, 1, mxREAL);
	acf = mxGetPr(ACF);
	
	
	/* Do the actual computation */
	int *acfInt = (int *)malloc(maxShift*sizeof(int));
	frame2acf(frameInt, frameSize, acfInt, maxShift, method, useMeanSubtraction);
	int2doubleVec(acfInt, maxShift, acf);

	free(frameInt);
	free(acfInt);
}
