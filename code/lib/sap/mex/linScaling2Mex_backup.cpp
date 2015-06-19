// Linear scaling
// How to compile:
// mex -I/users/jang/c/lib/audio -I/users/jang/c/lib/utility -I/users/jang/c/lib/wave linScalingMex.cpp d:/users/jang/c/lib/wave/waveRead4pda.cpp d:/users/jang/c/lib/audio/audio.cpp d:/users/jang/c/lib/utility/utility.cpp -output linScalingMex.dll

#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define SCALEDPITCHVEC		prhs[0]		// 使用者的輸入音高(after linear scaling to create scaled versions)
#define SCALEDPITCHLEN		prhs[1]		// Length of each scaled vectors
#define STDPITCH		prhs[2]		// 標準答案的音高
#define DISTANCETYPE		prhs[3]		// 1: L1 norm (shift to median), 2: L2 norm (shift to mean)
#define DISTBOUND		prhs[4]		// Distance bound (for partial computation)
/* Output Arguments */
#define	MINDIST			plhs[0]
#define	MININDEX		plhs[1]
#define	ALLDIST			plhs[2]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	/* Default values */
	int i, distanceType=1;	// L1 norm
	double distBound=INF;
	
	// Error checking
	if ((nrhs<2)||(nrhs>7)||(nlhs>3)) {
		char message[200];
		strcpy(message, mexFunctionName());
		strcat(message, " requires 2~6 input arguments and 0~3 output arguments.\n");
		strcat(message, "Usage: [minDist, minIndex, allDistance] = ");
		strcat(message, mexFunctionName());
		strcat(message, "(scaledPitch, stdPitch, distType, distBound)");
		mexErrMsgTxt(message);
	}
	if (nrhs>=4) distanceType = (int)mxGetScalar(DISTANCETYPE);
	if (nrhs>=5) distBound=mxGetScalar(DISTBOUND);
	
	double *scaledPitchVec=mxGetPr(SCALEDPITCHVEC);
	double *scaledPitchLen=mxGetPr(SCALEDPITCHLEN);
	double *stdPitch = mxGetPr(STDPITCH);
	int stdPitchLen = mxGetN(STDPITCH)*mxGetM(STDPITCH);
	int maxPitchLen=mxGetM(SCALEDPITCHVEC);
	int resolution=mxGetN(SCALEDPITCHVEC);

	double *distVec = new double [resolution];
	double *muVec = new double [resolution];
	double *newPitch = new double[maxPitchLen];		// Resampled pitch with the max. length
	double *tempPitch = new double[maxPitchLen];
	int minIndex;

	SCALEDVECSET scaledVecSet;
	scaledVecSet.maxPitchLen=maxPitchLen;
//	scaledVecSet.lowerRatio=lowerRatio;
//	scaledVecSet.upperRatio=upperRatio;
	scaledVecSet.resolution=resolution;
	scaledVecSet.scaledVec=scaledPitchVec;
	scaledVecSet.scaledVecLen=(int *)malloc(resolution*sizeof(int));
	for (i=0; i<resolution; i++)
		scaledVecSet.scaledVecLen[i]=scaledPitchLen[i];

//printf("Before calling linearScaling...\n");
//	linearScaledVecCreate(pitch, pitchLen, lowerRatio, upperRatio, resolution, &scaledVecSet);
	double distanceMin=linearScaling2(&scaledVecSet, stdPitch, stdPitchLen, distanceType, distVec, muVec, newPitch, tempPitch, minIndex, distBound);
//printf("After calling linearScaling...\n");
//	printf("minIndex=%d\n", minIndex);
//	printf("bestRatio=%f\n", lowerRatio+minIndex*(upperRatio-lowerRatio)/(resolution-1));

	//開始準備要輸出的資料
	MINDIST = mxCreateDoubleMatrix(1,1, mxREAL);
	mxGetPr(MINDIST)[0] = distanceMin;			//最小距離

	if (nlhs>1){	// 若需要內插及平移後的 pitch
		MININDEX = mxCreateDoubleMatrix(1, 1, mxREAL);
		mxGetPr(MININDEX)[0]=minIndex+1;	// Plus one to get MATLAB indexing
	}

	if (nlhs>2){	// 若需要 all distance
		ALLDIST = mxCreateDoubleMatrix(1, resolution, mxREAL);
		for (int i=0; i<resolution; i++)
			mxGetPr(ALLDIST)[i]=distVec[i];
	}
	delete [] distVec;
	delete [] muVec;
	delete [] newPitch;
	delete [] tempPitch;
	free(scaledVecSet.scaledVecLen);
}
