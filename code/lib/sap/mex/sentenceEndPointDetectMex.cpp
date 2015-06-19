#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	WAVE		prhs[0]	// Wave data
#define	FS		prhs[1]	// Sample rate
#define	NBITS		prhs[2]	// No. of bits
/* Output Arguments */
#define	EP	plhs[0]
//#define	VOLUME	plhs[1]
//#define	ZCR	plhs[2]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *wave, *ep, *volume, *zcr;
	int m, n, i, j, fs, nbits, waveLen, beginIndex, endIndex;
	cAudio myAudio;

	/* Check for proper number of arguments */
	if (nrhs!=3) {
		char message[200];
		strcpy(message, mexFunctionName());
		strcat(message, " requires 3 input arguments.\n");
		strcat(message, "Usage: [startIndex, endIndex] = ");
		strcat(message, mexFunctionName());
		strcat(message, "(wave, fs, nbits)");
		mexErrMsgTxt(message);
	}

	/* Dimensions of the input matrix */
	waveLen = mxGetM(WAVE)*mxGetN(WAVE);

	/* Create a matrix for the return argument */
	EP = mxCreateDoubleMatrix(1, 2, mxREAL);

	/* Assign pointers to the various parameters */
	wave = mxGetPr(WAVE);
	fs = (int) mxGetScalar(FS);
	nbits = (int) mxGetScalar(NBITS);
	ep = mxGetPr(EP);

	myAudio.sentenceEndPointDetect(wave, waveLen, fs, &beginIndex, &endIndex);
	ep[0]=beginIndex+1;	// MATLAB indexing
	ep[1]=endIndex+1;	// MATLAB indexing
}

