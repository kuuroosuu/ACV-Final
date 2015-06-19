// How to compile:
// MATLAB 7.1: mex epdByVolZcrMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp -output epdByVolZcrMex.dll
// Others: mex epdByVolZcrMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp

#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	WAVE		prhs[0]	// Wave data
#define	FS		prhs[1]	// Sample rate
/* Output Arguments */
#define	EP	plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *wave, *ep, *volume, *zcr;
	int m, n, i, j, fs, waveLen, index1, index2;
	char message[200];
	cAudio myAudio;

	/* Check for proper number of arguments */
	if (nrhs!=3) {
		strcpy(message, mexFunctionName());
		strcat(message, " requires 3 input arguments.\n");
		strcat(message, "Usage: ep = ");
		strcat(message, mexFunctionName());
		strcat(message, "(wave, fs)");
		mexErrMsgTxt(message);
	}

	/* Dimensions of the input matrix */
	waveLen = mxGetM(WAVE)*mxGetN(WAVE);
	int frameSize=320, overlap=240;		// 注意：這部分必須和 endPointDetect() 一致！
	int frameNum = (waveLen-overlap)/(frameSize-overlap);

	/* Create a matrix for the return argument */
	EP = mxCreateDoubleMatrix(1, 2, mxREAL);

	/* Assign pointers to the various parameters */
	wave = mxGetPr(WAVE);
	fs = (int) mxGetScalar(FS);
	ep = mxGetPr(EP);
	
	if (vecMax(wave, waveLen)<=1)
		mexErrMsgTxt("Input wave signals are not likely to be integers since all are within the range [-1 1]! Perhaps you should use wavReadInt.m instead?");

	int *waveInt = (int *)malloc(waveLen*sizeof(int));
	double2intVec(wave, waveLen, waveInt);
	myAudio.epdByVolZcr(waveInt, waveLen, fs, &index1, &index2);
	free(waveInt);
	
	ep[0]=index1+1;	// MATLAB index
	ep[1]=index2+1;	// MATLAB index
}
