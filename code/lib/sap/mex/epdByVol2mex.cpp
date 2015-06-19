// How to compile:
// MATLAB 7.1: mex epdByVol2mex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp -output epdByVol2mex.dll
// Others: mex epdByVol2mex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp

#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	WAVE		prhs[0]	// Wave data
#define	FS		prhs[1]	// Sample rate
/* Output Arguments */
#define	EPINDEX		plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	char message[200];
	cAudio myAudio;

	/* Check for proper number of arguments */
	if (nrhs!=2){
		strcpy(message, mexFunctionName());
		strcat(message, " requires 2 input arguments.\n");
		strcat(message, "Usage: epIndex = ");
		strcat(message, mexFunctionName());
		strcat(message, "(wave, fs)");
		mexErrMsgTxt(message);
	}

	/* Dimensions of the input matrix */
	int waveLen = mxGetM(WAVE)*mxGetN(WAVE);
	int frameSize=480, overlap=320;		// 注意：這部分必須和 endPointDetect03() 一致！
	int frameNum=(waveLen-overlap)/(frameSize-overlap);

	/* Assign pointers to the various parameters */
	double *wave = mxGetPr(WAVE);
	int fs = (int) mxGetScalar(FS);
	
	if (vecMax(wave, waveLen)-vecMin(wave, waveLen)<=2)
		mexErrMsgTxt("The range of the given wave vector is too small! Perhaps it should be converted to integer format...");

	int epIndexLen=20;	// 最多只能有 10 個 voicedSegment
	int *epIndex=(int *)malloc(epIndexLen*sizeof(int));
	int *waveInt=(int *)malloc(waveLen*sizeof(int));
	double2intVec(wave, waveLen, waveInt);
	int realEpIndexLen=myAudio.epdByVol2(waveInt, waveLen, fs, epIndex, epIndexLen);
	
	/* Create a matrix for the return argument */
	EPINDEX = mxCreateDoubleMatrix(1, realEpIndexLen, mxREAL);
	for (int i=0; i<realEpIndexLen; i++)
		mxGetPr(EPINDEX)[i]=epIndex[i]+1;		// MATLAB index

	free(epIndex);
	free(waveInt);
}
