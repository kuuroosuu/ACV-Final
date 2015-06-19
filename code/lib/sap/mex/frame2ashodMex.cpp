// Frame 2 asdf
// How to compile:
// MATLAB 7.1: mex frame2ashodMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp -output frame2ashodMex.dll
// Others: mex frame2ashodMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp
#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	FRAME		prhs[0]
#define	DIFFORDER	prhs[1]

/* Output Arguments */
#define	ASHOD	plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *frame, ashod;
	int frameSize, diffOrder;
	if (nrhs<2){
		char message[200];
		sprintf(message, "Usage: ashod = %s(vec, diffOrder)", mexFunctionName());
		mexErrMsgTxt(message);
	}
	
	/* Assign pointers to the various parameters */
	frame=mxGetPr(FRAME);
	frameSize=mxGetM(FRAME)*mxGetN(FRAME);
	diffOrder=mxGetScalar(DIFFORDER);
	ashod=frame2ashod(frame, frameSize, diffOrder);

	/* Create a matrix for the return argument */
	ASHOD = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxGetPr(ASHOD)[0]=ashod;
}

