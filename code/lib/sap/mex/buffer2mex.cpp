// To compile: mex buffer2mex.cpp d:\users\jang\c\lib\audio.cpp d:\users\jang\c\lib\utility.cpp
#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	Y		prhs[0]
#define	FRAMESIZE	prhs[1]
#define	OVERLAP		prhs[2]

/* Output Arguments */
#define	FRAMEMAT		plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *y, *frameMat;
	int frameSize, overlap, frameNum, yLen, i;
	char message[200];
	cAudio myAudio;

	/* Check for proper number of arguments */
	if (nrhs!=3) {
		strcpy(message, mexFunctionName());
		strcat(message, " requires 3 input arguments.\n");
		strcat(message, "Usage: frameMat = ");
		strcat(message, mexFunctionName());
		strcat(message, "(y, frameSize, overlap)");
		mexErrMsgTxt(message);
	}

	/* Assign pointers to the various parameters */
	y = mxGetPr(Y);
	yLen = mxGetM(Y)*mxGetN(Y);
	frameSize = (int)mxGetPr(FRAMESIZE)[0];
	overlap = (int)mxGetPr(OVERLAP)[0];

	if (overlap >= frameSize)
		mexErrMsgTxt("Overlap must be less than frame size!");
	
	/* Create a matrix for the return argument */
	frameNum = myAudio.getFrameNum2(yLen, frameSize, overlap);
	FRAMEMAT = mxCreateDoubleMatrix(frameSize, frameNum, mxREAL);
	frameMat = mxGetPr(FRAMEMAT);

	int *framePos = new int[frameNum];
	int *yInt = new int[yLen];
	for (i=0; i<yLen; i++)
		yInt[i]=y[i];
	myAudio.buffer2(yLen, frameSize, overlap, framePos);	// Generate correct framePos
	for (i=0; i<frameNum; i++)
		memcpy(frameMat+i*frameSize, y+framePos[i], frameSize*sizeof(double));
	delete [] yInt;
	delete [] framePos;
}
