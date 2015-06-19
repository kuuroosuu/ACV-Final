// Frame 2 ACF, double version
// How to compile:
// MATLAB 7.1: mex frame2acfMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp -output frame2acfMex.dll
// Others: mex frame2acfMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp
#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	FRAMEMAT	prhs[0]
#define USEHAMMINGWIN	prhs[1]

/* Output Arguments */
#define	MAGSPEC		plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *frameMat, *magSpec;
	int i, frameSize, frameNum, useHammingWin;
	cAudio myAudio;

	/* Check for proper number of arguments */
	if (nrhs<2){
		char message[200];
		strcpy(message, mexFunctionName());
		strcat(message, " requires 2 input arguments.\n");
		strcat(message, "Usage: magSpec = ");
		strcat(message, mexFunctionName());
		strcat(message, "(frameMat, useHammingWin)\n");
		strcat(message, "\tThe length of a frame must be equal to 2^n.\n");
		mexErrMsgTxt(message);
	}
	
	/* Assign pointers to the various parameters */
	frameMat = mxGetPr(FRAMEMAT);
	useHammingWin = (int)mxGetScalar(USEHAMMINGWIN);
	frameSize = mxGetM(FRAMEMAT);
	frameNum = mxGetN(FRAMEMAT);

	/* Create a matrix for the return argument */
	MAGSPEC = mxCreateDoubleMatrix(frameSize, frameNum, mxREAL);
	magSpec = mxGetPr(MAGSPEC);
	/* Do the actual computation */
	myAudio.absFft4frameMat(frameMat, frameSize, frameNum, useHammingWin, magSpec);
}
