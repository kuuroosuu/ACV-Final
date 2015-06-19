// How to compile:
// mex -Id:/users/jang/c/lib -Id:/users/jang/c/lib/audio -Id:/users/jang/c/lib/utility -Id:/users/jang/c/lib/wave ptByPfMex.cpp d:/users/jang/c/lib/audio/audio.cpp d:/users/jang/c/lib/utility/utility.cpp d:/users/jang/c/lib/wave/waveRead4pda.cpp -output ptByPfMex.dll

#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	WAVE			prhs[0]
#define	FS			prhs[1]
#define	NBITS			prhs[2]
#define	PT_PARAM		prhs[3]
#define	DATA_MODE		prhs[4]		// 0 for int, 1 for double

/* Output Arguments */
#define	PITCH	plhs[0]
#define	CLARITY	plhs[1]
#define	PFMAT	plhs[2]

#define	printf mexPrintf

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double	*wave, *pfMat, *pitch, *clarity;
	int i, waveLen, *waveInt, debug=1;
	cAudio myAudio;

	/* Check for proper number of arguments */
	if ((nrhs<3) || (nrhs>5)){
		char message[200];
		sprintf(message, "Error!\n");
		sprintf(message, "%sUsage 1: [pitch clarity pfMat] = %s(wave, fs, nbits)\n", message, mexFunctionName());
		sprintf(message, "%sUsage 2: [pitch clarity pfMat] = %s(wave, fs, nbits, ptParam)\n", message, mexFunctionName());
		sprintf(message, "%s\t(ptParam can be obtain via ptParamSet.m.)\n", message);
		sprintf(message, "%sUsage 3: [pitch clarity pfMat] = %s(wave, fs, nbits, ptParam, dataMode)\n", message, mexFunctionName());
		sprintf(message, "%s\t(dataMode = 0 for fixed-point, 1 for float-point operation, which is the default)\n", message);
		sprintf(message, "%s(If dataMod=0, pfType=acf, fs=16bits, then you run into the risk of overflow!)\n", message);
		mexErrMsgTxt(message);
	}

	waveLen = mxGetM(WAVE)*mxGetN(WAVE);
	wave = mxGetPr(WAVE);
	if ((-1<=vecMin(wave, waveLen)) && (vecMax(wave, waveLen)<=1)){
		printf("Warning: The input wave should be transformed into integers in the range [-128, 127] or [-32768, 32767].\n");
	//	sprintf(message, "Error: The input wave should be transformed into integers in the range [-128, 127] or [-32768, 32767].\n");
	//	mexErrMsgTxt(message);
	}
	
	int fs = (int)mxGetScalar(FS);
	int nbits = (int)mxGetScalar(NBITS);
	// Create parameters for pitch tracking
	PTPARAM ptParam; myAudio.ptParamSet(&ptParam, fs, nbits);
//	if (debug)
myAudio.ptParamPrint(&ptParam);
	if (nrhs>=4){
		mxArray *fieldValue;
		fieldValue=mxGetField(PT_PARAM, 0, "frameSize");		ptParam.frameSize=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "overlap");			ptParam.overlap=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "pfType");			ptParam.pfType=(PF_TYPE)(int)*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "pfWeight");			ptParam.pfWeight=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "indexDiffWeight");		ptParam.indexDiffWeight=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "dpMethod");			ptParam.dpMethod=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "dpPosDiffOrder");		ptParam.dpPosDiffOrder=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "dpUseLocalOptim");		ptParam.dpUseLocalOptim=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "minPitch");			ptParam.minPitch=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "maxPitch");			ptParam.maxPitch=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "pitchRangeMax");		ptParam.pitchRangeMax=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "pitchDiffMax");		ptParam.pitchDiffMax=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "useVolThreshold");		ptParam.useVolThreshold=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "cutLeadingTrailingZero");	ptParam.cutLeadingTrailingZero=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "medianFilterOrder");	ptParam.medianFilterOrder=*mxGetPr(fieldValue);
	}
//mexPrintf("ptParam.pfWeight=%d\n", ptParam.pfWeight);
	myAudio.ptParamWrite(&ptParam, (char *)"ptParam.output");
	int dataMode = (nrhs<5)?1:(int)mxGetScalar(DATA_MODE);		// The default value is 1 for floating-point operation
//	mexPrintf("nrhs=%d, dataMode=%d\n", nrhs, dataMode);
	int frameSize = ptParam.frameSize;
	int overlap = ptParam.overlap;
	
	// ====== Determine frameNum
	int frameNum=myAudio.getFrameNum2(waveLen, ptParam.frameSize, ptParam.overlap);
	int pitchLen, pfLen=frameSize/2;
	waveInt = (int *)malloc(waveLen*sizeof(int));
	double2intVec(wave, waveLen, waveInt);
vecWrite(waveInt, waveLen, "wave.txt");
	// Output matrix
	PITCH = mxCreateDoubleMatrix(1, frameNum, mxREAL);
	pitch = mxGetPr(PITCH);
	CLARITY = mxCreateDoubleMatrix(1, frameNum, mxREAL);
	clarity = mxGetPr(CLARITY);
	PFMAT = mxCreateDoubleMatrix(pfLen, frameNum, mxREAL);
	pfMat = mxGetPr(PFMAT);
	
	if (dataMode==0){		// data type is int
		int *pitchInt=(int *)malloc(frameNum*sizeof(int));
		int *clarityInt=(int *)malloc(frameNum*sizeof(int));
		int *pfMatInt = (int *)malloc(pfLen*frameNum*sizeof(int));
		pitchLen=ptByPf(waveInt, waveLen, ptParam, pitchInt, clarityInt, pfMatInt);
		int2doubleVec(pitchInt, frameNum, pitch);
		int2doubleVec(clarityInt, frameNum, clarity);
		int2doubleVec(pfMatInt, pfLen*frameNum, pfMat);
	//	for (i=0; i<frameNum; i++){
	//		pitch[i]/=1024;		// Scale down if you can detect type of T
	//		clarity[i]/=1024;	// Scale down if you can detect type of T
	//	}
		free(pitchInt);
		free(clarityInt);
		free(pfMatInt);
	} else if (dataMode==1){	// data type is double
		pitchLen=ptByPf(waveInt, waveLen, ptParam, pitch, clarity, pfMat);
	} else {
		char message[200];
		sprintf(message, "Unknown data mode = %d", dataMode);
		mexErrMsgTxt(message);
	}
	if (pitchLen!=frameNum){	// Pitch length is different from frameNum due to removal of leading/trailing zeros, etc
		double *temp=(double *)malloc(pitchLen*sizeof(double));
		// For pitch
		for (i=0; i<pitchLen; i++)
			temp[i]=pitch[i];
		mxDestroyArray(PITCH);		// Delete obsolete output
		PITCH = mxCreateDoubleMatrix(1, pitchLen, mxREAL);	// Create a new output
		pitch = mxGetPr(PITCH);
		for (i=0; i<pitchLen; i++)
			pitch[i]=temp[i];
		// For clarity
		for (i=0; i<pitchLen; i++)
			temp[i]=clarity[i];
		mxDestroyArray(CLARITY);		// Delete obsolete output
		CLARITY = mxCreateDoubleMatrix(1, pitchLen, mxREAL);	// Create a new output
		clarity = mxGetPr(CLARITY);
		for (i=0; i<pitchLen; i++)
			clarity[i]=temp[i];
		// Free memory
		free(temp);
	}

// Why we cannot allocate memory inside switch?
/*
	switch(dataMode){
		case 0:	// ====== Pitch tracking by int PF min/max picking
			int *pitchInt=(int *)malloc(frameNum*sizeof(int));
			int *clarityInt=(int *)malloc(frameNum*sizeof(int));
			int *pfMatInt = (int *)malloc(pfLen*frameNum*sizeof(int));
			ptByPf(waveInt, waveLen, ptParam, pitchInt, clarityInt, pfMatInt);
			int2doubleVec(pitchInt, frameNum, pitch);
			int2doubleVec(clarityInt, frameNum, clarity);
			int2doubleVec(pfMatInt, pfLen*frameNum, pfMat);
			for (i=0; i<frameNum; i++){
				pitch[i]/=1024;		// Scale down
				clarity[i]/=1024;	// Scale down
			}
			free(pitchInt);
			free(clarityInt);
			free(pfMatInt);
			break;
		case 1:	// ====== Pitch tracking by float PF min/max picking
			ptByPf(waveInt, waveLen, ptParam, pitch, clarity, pfMat);
			break;
		default:
			myExit("Unknown mode in searchSylVs()!");
	}
*/
	free(waveInt);
}
