// How to compile:
// mex -I/users/jang/c/lib -I/users/jang/c/lib/audio -I/users/jang/c/lib/utility ptByDpOverPfMex.cpp /users/jang/c/lib/audio/audio.cpp /users/jang/c/lib/utility/utility.cpp

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
#define DPPATH	plhs[3]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double	*wave, *pfMat, *origPitch, *pitch, *origClarity, *clarity, *dpPath;
	int i, waveLen, pitchLen, *waveInt, debug=0;
	cAudio myAudio;
	char message[200];

	/* Check for proper number of arguments */
	if ((nrhs<3) || (nrhs>4)){
		sprintf(message, "Format error!\nUsage 1: out = %s(wave, fs, nbits)\nUsage 2: out = %s(wave, fs, nbits, ptParam)", mexFunctionName(), mexFunctionName());
		mexErrMsgTxt(message);
	}

	waveLen = mxGetM(WAVE)*mxGetN(WAVE);
	wave = mxGetPr(WAVE);
	
	if ((-1<=vecMin(wave, waveLen)) && (vecMax(wave, waveLen)<=1)){
		sprintf(message, "Error: The input wave should be transformed into integers in the range [-128, 127] or [-32768, 32767].\n");
		mexErrMsgTxt(message);
	}
	
	int fs = (int)mxGetScalar(FS);
	int nbits = (int)mxGetScalar(NBITS);

	// Create parameters for pitch tracking
//	PRINT(fs); PRINT(nbits);	// Why this doesn't work?
//	printf("fs=%d, nbits=%d\n", fs, nbits);
	PTPARAM ptParam; myAudio.ptParamSet(&ptParam, fs, nbits, ACF);	// The last argument is AMDF or ACF
	if (debug) myAudio.ptParamPrint(&ptParam);

	if (nrhs==4){
	//	char* tempArray[] = {"frameSize", "overlap", "pfLen", "pfType"};
	//	int len=sizeof(tempArray)/sizeof(tempArray[0]);
	//	printf("Array length=%d\n", len);
	//	const vector<string> fieldName(tempArray, tempArray+len);	// Copy to a string vector in C++
	//	for(int i=0; i<fieldName.size(); i++){
	//	}
		
		mxArray *fieldValue;
		fieldValue=mxGetField(PT_PARAM, 0, "frameSize");		if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.frameSize=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "overlap");			if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.overlap=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "pfLen");			if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.pfLen=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "pfType");			if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.pfType=(PF_TYPE)(int)*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "pfWeight");			if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.pfWeight=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "indexDiffWeight");		if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.indexDiffWeight=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "dpMethod");			if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.dpMethod=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "dpPosDiffOrder");		if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.dpPosDiffOrder=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "dpUseLocalOptim");		if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.dpUseLocalOptim=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "minPitch");			if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.minPitch=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "maxPitch");			if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.maxPitch=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "pitchRangeMax");		if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.pitchRangeMax=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "pitchDiffMax");		if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.pitchDiffMax=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "useVolThreshold");		if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.useVolThreshold=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "useClarityThreshold");	if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.useClarityThreshold=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "cutLeadingTrailingZero");	if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.cutLeadingTrailingZero=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "medianFilterOrder");	if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.medianFilterOrder=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "clarityRatio");		if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.clarityRatio=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "pfWeightedByClarity");	if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.pfWeightedByClarity=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "useParabolicFit");		if (fieldValue==NULL) mexErrMsgTxt("Field not exist!"); ptParam.useParabolaFit=*mxGetPr(fieldValue);
	}
	if (debug) myAudio.ptParamPrint(&ptParam);
//	mexPrintf("ptParam.frameSize=%d\n", ptParam.frameSize);
//	mexPrintf("ptParam.overlap=%d\n", ptParam.overlap);
//	mexPrintf("ptParam.indexDiffWeight=%d\n", ptParam.indexDiffWeight);
//	mexPrintf("ptParam.pfWeight=%d\n", ptParam.pfWeight);

	int frameSize=ptParam.frameSize;
	int overlap=ptParam.overlap;
	int pfLen=ptParam.pfLen;
	
	// ====== Determine frameNum
	int frameNum=myAudio.getFrameNum2(waveLen, ptParam.frameSize, ptParam.overlap);

	// Using int wave. This is required by epdByVol
	waveInt=(int *)malloc(waveLen*sizeof(int));
	double2intVec(wave, waveLen, waveInt);
	
	origPitch=(double *)malloc(frameNum*sizeof(double));
	origClarity=(double *)malloc(frameNum*sizeof(double));
	
	// ====== Pitch tracking by DP using int pfMat
//	int *pfMatInt = (int *)malloc(pfLen*frameNum*sizeof(int));
//	myAudio.ptByDpOverPf(waveInt, waveLen, ptParam, pitch, pfMatInt);
//	int2doubleVec(pfMatInt, pfLen*frameNum, pfMat);
//	free(pfMatInt);

	// Output matrix for pfMat
	PFMAT = mxCreateDoubleMatrix(pfLen, frameNum, mxREAL);
	pfMat = mxGetPr(PFMAT);
	
	DPPATH = mxCreateDoubleMatrix(1, frameNum, mxREAL);
	dpPath = mxGetPr(DPPATH);

	int *dpPathInt=(int *)malloc(frameNum*sizeof(int));	// ptByDpOverPf() takes int data type.
	
	// ====== Pitch tracking by DP using double pfMat
	pitchLen=ptByDpOverPf(waveInt, waveLen, ptParam, origPitch, origClarity, pfMat, dpPathInt);
//	vecWrite(pfMat, pfLen*frameNum, "pfMat.txt");
//	ptByDpOverPf(wave, waveLen, ptParam, pitch, pfMat);

	// Output matrix for pitch. It's placed here since the pitch length might be changed due to leading/trailing zero removal.
	PITCH = mxCreateDoubleMatrix(1, pitchLen, mxREAL);
	pitch = mxGetPr(PITCH);
	memcpy(pitch, origPitch, pitchLen*sizeof(double));
	
	CLARITY = mxCreateDoubleMatrix(1, pitchLen, mxREAL);
	clarity = mxGetPr(CLARITY);
	memcpy(clarity, origClarity, pitchLen*sizeof(double));
	
	for (i=0; i<frameNum; i++)
		dpPath[i]=dpPathInt[i]++;		// For MATLAB indexing

//	for (i=0; i<frameNum; i++)
//		pitch[i]=69.0+12.0*log(ptParam.fs/(double)dpPathInt[i]/440.0)/LOG2;		// log(x)=ln(x)

	// Free memory
	free(waveInt);
	free(origPitch);
	free(origClarity);
	free(dpPathInt);
}
