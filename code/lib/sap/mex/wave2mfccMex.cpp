// Wave to mfcc conversion
// How to compile:
// MATLAB 7.1: mex wave2mfccMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/fft.cpp d:/users/jang/c/lib/myMfcc.cpp d:/users/jang/c/lib/utility.cpp -output wave2amdfMatMex.dll
// Others: mex wave2mfccMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp d:/users/jang/c/lib\mfcc\CMEL.cpp d:/users/jang/c/lib\mfcc\CFEABUF.cpp d:/users/jang/c/lib\mfcc\CSigP.cpp

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <CONIO.H>
#include <iostream>
#include "mex.h"
#include "audio.hpp"
#include "CMel.h"
#include "CWavePCM.h"

/* Input Arguments */
#define INPUTWAVE		prhs[0]
#define FS			prhs[1]
#define NBITS			prhs[2]
#define CFGFILE			prhs[3]
/* Output Arguments */
#define MFCC plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *inputWave, *mfcc;
	int i, frameNum, feaDim;
	char message[200], cfgFile[200];
	
	if (nrhs<4){	
		sprintf(message, "Usage: mfcc = %s(y, fs, nbits, cfgFile)", mexFunctionName());
		mexErrMsgTxt(message);
	}
	
	/* Assign pointers to the various parameters */
	inputWave = mxGetPr(INPUTWAVE);
	int waveLen= mxGetM(INPUTWAVE)*mxGetN(INPUTWAVE);
	int fs=(int)mxGetScalar(FS);
	int nbits=(int)mxGetScalar(NBITS);
	
	int status = mxGetString(CFGFILE, cfgFile, 200);
	if(status != 0) 
		mexWarnMsgTxt("Not enough space. String is truncated.");

	if (vecMax(inputWave, waveLen)<=1)
		mexPrintf("Input wave signals are not likely to be integers since all are within the range [-1 1]! Perhaps you should use wavReadInt.m instead?\n");
	
	fstream file;
	file.open(cfgFile, ios::in);
	if (!file){
		file.close();
		sprintf(message, "Cannot open config file %s!", cfgFile);
		mexErrMsgTxt(message);
	}

	int *wave = new int[waveLen];
	for (i=0; i<waveLen; i++)
		wave[i]=inputWave[i];

	CMel myMel;
	myMel.wave2mfcc(wave, fs, nbits, waveLen, cfgFile, &frameNum, &feaDim);
	float *feature=new float[frameNum*feaDim];
	memcpy(feature, myMel.m_feature, frameNum*feaDim*sizeof(float));
	
	delete [] wave;
	
//	mexPrintf("fs=%d\n", fs);
//	mexPrintf("nbits=%d\n", nbits);
//	mexPrintf("waveLen=%d\n", waveLen);
//	mexPrintf("cfgFile=%s\n", cfgFile);
//	mexPrintf("frameNum=%d\n", frameNum);
//	mexPrintf("feaDim=%d\n", feaDim);
	
	// Create output
	MFCC = mxCreateDoubleMatrix(feaDim, frameNum, mxREAL);
	mfcc = mxGetPr(MFCC);
	for (i=0; i<feaDim*frameNum; i++)
		mfcc[i]=feature[i];
}
