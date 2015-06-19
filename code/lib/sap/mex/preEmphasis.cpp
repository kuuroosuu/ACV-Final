// Pre-emphasis
// How to compile:
// MATLAB 7.1: mex preEmphasis.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp -output wave2amdfMatMex.dll
// Others: mex preEmphasis.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp

#include <string.h>
#include <math.h>
#include <stdio.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	INWAV  prhs[0]
/* Output Arguments */
#define	OUTWAV plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *inwav, *outwav;
	double b[2]={1, -0.95};
	double a[1]={1};
	int i,index,wavLen;
	cAudio myAudio;
	
	/* Check for proper number of arguments */
	if (nrhs!=1) {
		char message[200];
		strcpy(message, mexFunctionName());
		strcat(message, " requires one input arguments.");
		mexErrMsgTxt(message);
	}

	/* Assign pointers to the various parameters */
	inwav = mxGetPr(INWAV);
	
	/* Create a matrix for the return argument */
	wavLen = mxGetM(INWAV)*mxGetN(INWAV);
	OUTWAV    = mxCreateDoubleMatrix(mxGetM(INWAV), mxGetN(INWAV), mxREAL);
	outwav    = mxGetPr(OUTWAV);
	memcpy(outwav, inwav, wavLen);
	
	//MATLAB code : Inwav = Inwav - 0.95*[0 ;Inwav(1:end-1)];
	myAudio.filter(b, 2, a, 1, outwav, wavLen);	
}
