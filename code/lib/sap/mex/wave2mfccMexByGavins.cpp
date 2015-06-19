// Wave to mfcc conversion
// How to compile:
// MATLAB 7.1: mex wave2mfccMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/fft.cpp d:/users/jang/c/lib/myMfcc.cpp d:/users/jang/c/lib/utility.cpp -output wave2amdfMatMex.dll
// Others: mex wave2mfccMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/fft.cpp d:/users/jang/c/lib/myMfcc.cpp d:/users/jang/c/lib/utility.cpp

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"
#include "audio.hpp"
#include "mfccByGavins.h"

/* Input Arguments */
#define INPUTWAVE	prhs[0]
#define FS		prhs[1]
#define MP		prhs[2]
/* Output Arguments */
#define PARAM plhs[0]

typedef struct {
	int frameSize;
	int overlap;
	int tbfNum;
	int cepsNum;
	int useEnergy;
	int useDelta;
} MPTYPE;

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *inputWave, *parameter, *wave, *out, *frameMat, *energy, *legalFrame, *logEnergy, *wFrameMat , *magSpecMat , *tbfCoef, *cepsMat, *mfcc0;
	double sum, normalize_coef;
	int i, j, k, featureNum, frameCount, frameCount2, legalFrameCount;
	double b[2] = {1, -0.95};
	double a[1] = {1};
	cAudio myAudio;

	/* Check for proper number of arguments */
	if (nrhs<3){
		char message[200];
		sprintf(message, "Usage: mfcc = %s(y, fs, mp)", mexFunctionName());
		mexErrMsgTxt(message);
	}
	
	/* Assign pointers to the various parameters */
	inputWave = mxGetPr(INPUTWAVE);
	int waveLen= mxGetM(INPUTWAVE)*mxGetN(INPUTWAVE);
	int fs=(int)mxGetScalar(FS);
	MPTYPE mp;
	mxArray *fieldValue;
	fieldValue=mxGetField(MP, 0, "frameSize");	mp.frameSize=*mxGetPr(fieldValue);
	fieldValue=mxGetField(MP, 0, "overlap");	mp.overlap=*mxGetPr(fieldValue);
	fieldValue=mxGetField(MP, 0, "tbfNum");		mp.tbfNum=*mxGetPr(fieldValue);
	fieldValue=mxGetField(MP, 0, "cepsNum");	mp.cepsNum=*mxGetPr(fieldValue);
	fieldValue=mxGetField(MP, 0, "useEnergy");	mp.useEnergy=*mxGetPr(fieldValue);
	fieldValue=mxGetField(MP, 0, "useDelta");	mp.useDelta=*mxGetPr(fieldValue);

//mexPrintf("mp.frameSize=%d\n", mp.frameSize);
//mexPrintf("mp.overlap=%d\n", mp.overlap);
//mexPrintf("mp.tbfNum=%d\n", mp.tbfNum);
//mexPrintf("mp.cepsNum=%d\n", mp.cepsNum);
//mexPrintf("mp.useEnergy=%d\n", mp.useEnergy);
//mexPrintf("mp.useDelta=%d\n", mp.useDelta);
	
	// Get input wave
	wave = (double *)malloc(waveLen*sizeof(double));
	memcpy(wave, inputWave, waveLen*sizeof(double));
	/* ====== Step 1: Pre-emphasis */
	myAudio.filter(b, 2, a, 1, wave, waveLen);
	/* ====== Step 2: Frame blocking */
	/* Get frame , energy , legalFrameCount */
	frameCount = myAudio.getFrameNum(waveLen, mp.frameSize, mp.overlap);
//mexPrintf("waveLen=%d, frameCount=%d\n", waveLen, frameCount);
	frameMat = (double *)malloc(mp.frameSize*frameCount*sizeof(double));
	buffer(wave, waveLen, frameMat, mp.frameSize, mp.overlap);

/*
int frameNum = getFrameNum2(waveLen, mp.frameSize, mp.overlap);
int *framePos = new int[frameNum];
buffer2(waveLen, mp.frameSize, mp.overlap, framePos);
mexPrintf("waveLen=%d, frameNum=%d\n", waveLen, frameNum);

	energy = (double *)malloc(frameCount*sizeof(double));
	getEnergy(frameMat, mp.frameSize, frameCount, energy);
	vecWrite(energy, frameCount, "energy.txt");

	// Compute legalFrameCount
	legalFrameCount = 0;
	for (i=0; i<frameCount; i++)
		if (energy[i]>=threshold)
			legalFrameCount++;
printf("frameCount=%d, legalFrameCount=%d\n", frameCount, legalFrameCount);

	// ====== Step 3: Edge Detection process (Get legalFrame)
	legalFrame = (double *)malloc(mp.frameSize*legalFrameCount*sizeof(double));
	edgeDetect(frameMat, legalFrame, energy, mp.frameSize, frameCount);
	// Compute log energy
*/

legalFrameCount=frameCount;
legalFrame=frameMat;

	logEnergy = (double *)malloc(legalFrameCount*sizeof(double));
	getLogEnergy(legalFrame, mp.frameSize, legalFrameCount, logEnergy);
	// ====== Step 4: Hamming windowing
	wFrameMat = (double *)malloc(mp.frameSize*legalFrameCount*sizeof(double));
	hammingWin(wFrameMat, legalFrame, mp.frameSize, legalFrameCount);
	// ====== Step 5: FFT process
	/* Get magSpecMat */
	magSpecMat = (double *)malloc(mp.frameSize*legalFrameCount*sizeof(double));
	myAudio.absFft4frameMat(wFrameMat, mp.frameSize, legalFrameCount, 1, magSpecMat);
	// ====== Step 6: Triangular Bandpass Filter process
	/* Get tbfCoef */
	tbfCoef = (double *)malloc(mp.tbfNum*legalFrameCount*sizeof(double));
	getTbfCoef(magSpecMat, mp.frameSize, legalFrameCount, tbfCoef, mp.tbfNum);
	// ===== Step 7: Mel Cepstrum process (Get cepstral coefficients)
	cepsMat = (double *)malloc(mp.cepsNum*legalFrameCount*sizeof(double));
	getCEPS(tbfCoef, mp.tbfNum, legalFrameCount, cepsMat, mp.cepsNum);
	// ====== Combine cepstrum with log energy to get mfcc0, if necessary
	featureNum=mp.cepsNum+mp.useEnergy;
	mfcc0=(double *)malloc(featureNum*legalFrameCount*sizeof(double));
	for (i=0; i<legalFrameCount; i++){
		memcpy(mfcc0+i*featureNum, cepsMat+i*mp.cepsNum, mp.cepsNum*sizeof(double));
		if (mp.useEnergy)
			mfcc0[i*featureNum+mp.cepsNum]=logEnergy[i];
	}
	// ====== Step 8: Delta function process
	switch (mp.useDelta){
		case 0:
			PARAM = mxCreateDoubleMatrix(featureNum, legalFrameCount, mxREAL);
			parameter = mxGetPr(PARAM);
			memcpy(parameter, mfcc0, featureNum*legalFrameCount*sizeof(double));
			break;
		case 1:
			PARAM = mxCreateDoubleMatrix(featureNum*2,legalFrameCount, mxREAL);
			parameter = mxGetPr(PARAM);
			ComputeDELTAFCN(mfcc0, parameter, featureNum, legalFrameCount);
			break;
		default:
			mexErrMsgTxt("Unsupported mp.useDelta!\n");
	}

	free(wave);
//	free(framePos);
	free(mfcc0);
	free(cepsMat);
	free(tbfCoef);
	free(magSpecMat);
	free(wFrameMat);
	free(logEnergy);
	free(legalFrame);
	free(energy);
	free(frameMat);
}
