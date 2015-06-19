// How to compile:
// MATLAB 7.1: mex epdByVolHodMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp -output epdByVolHodMex.dll
// Others: mex -Id:/users/jang/c/lib epdByVolHodMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp

#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	WAVE		prhs[0]	// Wave data
#define	FS		prhs[1]	// Sample rate
#define	NBITS		prhs[2]	// Bit resolution (not used!)
#define	EPD_PARAM	prhs[3]
/* Output Arguments */
#define	EP_IN_SAMPLE	plhs[0]
#define	EP_IN_FRAME	plhs[1]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *wave;
	int m, n, i, j, fs, nbits, waveLen, index1, index2, segmentNum;
	char message[200];
	cAudio myAudio;

	/* Check for proper number of arguments */
	if ((nrhs<2)||(nrhs>4)){
		char message[200];
		sprintf(message, "Usage: mfcc = %s(y, fs)\n\tmfcc = %s(y, fs, nbits, epdParam)", mexFunctionName(), mexFunctionName());
		mexErrMsgTxt(message);
	}

	/* Dimensions of the input matrix */
	waveLen = mxGetM(WAVE)*mxGetN(WAVE);

	/* Create a matrix for the return argument */
	EP_IN_SAMPLE = mxCreateDoubleMatrix(1, 2, mxREAL);
	EP_IN_FRAME  = mxCreateDoubleMatrix(1, 2, mxREAL);

	/* Assign pointers to the various parameters */
	wave = mxGetPr(WAVE);
	fs = (int) mxGetScalar(FS);
	
	if (vecMax(wave, waveLen)<=1)
		mexErrMsgTxt("Input wave signals are not likely to be integers since all are within the range [-1 1]! Perhaps you should use wavReadInt.m instead?");

	int *waveInt = (int *)malloc(waveLen*sizeof(int));
	double2intVec(wave, waveLen, waveInt);
	EPDPARAM epdParam; myAudio.epdParamSet(fs, epdParam);
	if (nrhs==4){
		mxArray *fieldValue;
		fieldValue=mxGetField(EPD_PARAM, 0, "frameSize");		epdParam.frameSize=*mxGetPr(fieldValue);
		fieldValue=mxGetField(EPD_PARAM, 0, "overlap");			epdParam.overlap=*mxGetPr(fieldValue);
		fieldValue=mxGetField(EPD_PARAM, 0, "volRatio");		epdParam.overlap=*mxGetPr(fieldValue);
		fieldValue=mxGetField(EPD_PARAM, 0, "vhRatio");			epdParam.vhRatio=*mxGetPr(fieldValue);
		fieldValue=mxGetField(EPD_PARAM, 0, "diffOrder");		epdParam.diffOrder=*mxGetPr(fieldValue);
		fieldValue=mxGetField(EPD_PARAM, 0, "volWeight");		epdParam.volWeight=*mxGetPr(fieldValue);
		fieldValue=mxGetField(EPD_PARAM, 0, "extendNum");		epdParam.extendNum=*mxGetPr(fieldValue);
		fieldValue=mxGetField(EPD_PARAM, 0, "minSegment");		epdParam.minSegment=*mxGetPr(fieldValue);
		fieldValue=mxGetField(EPD_PARAM, 0, "maxSilBetweenWord");	epdParam.maxSilBetweenWord=*mxGetPr(fieldValue);
		fieldValue=mxGetField(EPD_PARAM, 0, "minLastWordDuration");	epdParam.minLastWordDuration=*mxGetPr(fieldValue);
	}
	SEGMENT *segment=myAudio.epdByVolHod(waveInt, waveLen, fs, epdParam, &segmentNum);
	
	index1 = myAudio.frame2sampleIndex(segment[0].begin, epdParam.frameSize, epdParam.overlap);		// Conversion of frame to sample index 
	index2 = myAudio.frame2sampleIndex(segment[segmentNum-1].end, epdParam.frameSize, epdParam.overlap);	// Conversion of frame to sample index
	mxGetPr(EP_IN_SAMPLE)[0]=index1+1;	// MATLAB index
	mxGetPr(EP_IN_SAMPLE)[1]=index2+1;	// MATLAB index
	mxGetPr(EP_IN_FRAME)[0]=segment[0].begin+1;		// MATLAB index
	mxGetPr(EP_IN_FRAME)[1]=segment[segmentNum-1].end+1;	// MATLAB index

	free(waveInt);
	free(segment);
}
