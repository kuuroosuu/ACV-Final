#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	WAVE	prhs[0]
#define	FS	prhs[1]
#define	DELAY	prhs[2]
#define	GAIN	prhs[2]

/* Output Arguments */
#define	WAVEOUT	plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *wave, *waveOut, delay, gain;
	int waveLen, fs;
	cAudio myAudio;
	
	/* Check for proper number of arguments */
	if (nrhs!=4){
		char message[200];
		sprintf(message, "Usage: out = %s(wave, fs, delay, gain)", mexFunctionName());
		mexErrMsgTxt(message);
	}

	/* Assign pointers to the various parameters */
	wave = mxGetPr(WAVE);
	waveLen = mxGetM(WAVE)*mxGetN(WAVE);
	fs = (int)mxGetPr(FS)[0];
	delay = mxGetPr(DELAY)[0];
	gain = mxGetPr(GAIN)[0];

	/* Create a matrix for the return argument */
	WAVEOUT = mxCreateDoubleMatrix(waveLen, 1, mxREAL);
	waveOut = mxGetPr(WAVEOUT);
	
	/* Do the actual computation */
	myAudio.waveEcho(wave, waveLen, fs, delay, gain, waveOut);
}
