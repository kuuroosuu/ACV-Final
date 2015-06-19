// Abs of spectrogram
// How to compile:
// MATLAB 7.1: mex frame2acfMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp -output frame2acfMex.dll
// Others: mex frame2acfMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp
#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	WAVE		prhs[0]
#define	FRAMESIZE	prhs[1]
#define	OVERLAP		prhs[2]
#define	PADDEDSIZE	prhs[3]
#define USEHAMMINGWIN	prhs[4]

/* Output Arguments */
#define	ABSSPECGRAM	plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *wave, *absSpecgram;
	int waveLen, frameSize, overlap, paddedSize, frameNum, useHammingWin;
	cAudio myAudio;

	/* Check for proper number of arguments */
	if (nrhs<2){
		char message[200];
		strcpy(message, mexFunctionName());
		strcat(message, " requires 5 input arguments.\n");
		strcat(message, "Usage: absSpecgram = ");
		strcat(message, mexFunctionName());
		strcat(message, "(wave, frameSize, overlap, paddedSize, useHammingWin)\n");
		strcat(message, "\tpaddedSize must be equal to 2^n.\n");
		mexErrMsgTxt(message);
	}
	
	/* Assign pointers to the various parameters */
	wave = mxGetPr(WAVE);
	waveLen=mxGetM(WAVE)*mxGetN(WAVE);
	frameSize = (int)mxGetScalar(FRAMESIZE);
	overlap = (int)mxGetScalar(OVERLAP);
	paddedSize = (int)mxGetScalar(PADDEDSIZE);
	useHammingWin = (int)mxGetScalar(USEHAMMINGWIN);
	frameNum=myAudio.getFrameNum2(waveLen, frameSize, overlap);

	/* Create a matrix for the return argument */
	ABSSPECGRAM = mxCreateDoubleMatrix(1+paddedSize/2, frameNum, mxREAL);
	absSpecgram = mxGetPr(ABSSPECGRAM);
	/* Do the actual computation */
	myAudio.absSpectrogram(wave, waveLen, frameSize, overlap, paddedSize, useHammingWin, absSpecgram);
}
