// Wave to pf matrix conversion
// How to compile:
// mex -Id:/users/jang/c/lib wave2pfMatMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp -output wave2pfMatMex
#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	WAVE 		prhs[0]
#define	FRAMESIZE	prhs[1]
#define	OVERLAP		prhs[2]
#define	PFLEN		prhs[3]
#define	PFTYPE		prhs[4]
#define	PFMETHOD	prhs[5]
/* Output Arguments */
#define	PFMAT	plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double	*wave, *pfMat;
	int frameSize, overlap, pfLen, waveLen, frameNum, pfMethod, *waveInt, *pfMatInt;
	PF_TYPE pfType;
	cAudio myAudio;

	/* Check for proper number of arguments */
	if (nrhs!=6){
		char message[200];
		sprintf(message, "Usage: out = %s(wave, frameSize, overlap, pfLen, pfType, pfMethod)\n\t pfType = 0 for AMDF, 1 for ACF\n", mexFunctionName());
		mexErrMsgTxt(message);
	}

	waveLen = mxGetM(WAVE)*mxGetN(WAVE);
	wave = mxGetPr(WAVE);
	frameSize = (int)mxGetScalar(FRAMESIZE);
	overlap = (int)mxGetScalar(OVERLAP);
	pfLen = (int)mxGetScalar(PFLEN);
	pfType = (PF_TYPE)(int)mxGetScalar(PFTYPE);			// 0 for AMDF, 1 for ACF
	pfMethod = (int)mxGetScalar(PFMETHOD);
	frameNum=myAudio.getFrameNum2(waveLen, frameSize, overlap);	// 音框個數
	
//	waveInt = (int *)malloc(waveLen*sizeof(int));
//	double2intVec(wave, waveLen, waveInt);

	// 輸出矩陣
	PFMAT = mxCreateDoubleMatrix(pfLen, frameNum, mxREAL);
	pfMat = mxGetPr(PFMAT);

// Integer version, which cause ACF to overflow
/*
	pfMatInt = (int *)malloc(pfLen*frameNum*sizeof(int));
	// 實際計算
	wave2pfMat(waveInt, waveLen, frameSize, overlap, pfLen, pfType, pfMethod, pfMatInt);	
	int2doubleVec(pfMatInt, pfLen*frameNum, pfMat);
	// Free memory
	free(pfMatInt);
*/

// Floating version
//	wave2pfMat(waveInt, waveLen, frameSize, overlap, pfType, pfMethod, pfMat);
	wave2pfMat(wave, waveLen, frameSize, overlap, pfLen, pfType, pfMethod, pfMat);

//	free(waveInt);
}
