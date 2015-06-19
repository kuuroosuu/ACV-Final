// How to compile:
// mex zcRateMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp
// mex zcRateMex.cpp d:/users/jang/c/lib/audio.lib d:/users/jang/c/lib/utility.cpp
// mex zcRateMex.cpp d:/users/jang/c/lib/audio.lib utility.lib === xxx!

//#pragma comment(lib, "audio.lib")	// easy compiling: mex zcRateMex.cpp utility.cpp
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	WAVE		prhs[0]	// Wave data
/* Output Arguments */
#define	OUT	plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double	*wave, *out;
	int waveLen;
	
	/* Check for proper number of arguments */
	if (nrhs!=1) {
		char message[200];
		strcpy(message, mexFunctionName());
		strcat(message, " requires 1 input arguments.\n");
		strcat(message, "Usage: output = ");
		strcat(message, mexFunctionName());
		strcat(message, "(wave)");
		mexErrMsgTxt(message);
	}

	/* Dimensions of the input matrix */
	waveLen = mxGetM(WAVE)*mxGetN(WAVE);

	/* Create a matrix for the return argument */
	OUT = mxCreateDoubleMatrix(1, 1, mxREAL);

	/* Assign pointers to the various parameters */
	wave = mxGetPr(WAVE);
	out = mxGetPr(OUT);
	out[0]=vecZcr(wave, waveLen);
}

