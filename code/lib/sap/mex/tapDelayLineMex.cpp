// How to compile:
// mex -I/users/jang/c/lib/audio tapDelayLineMex.cpp /users/jang/c/lib/audio/audioEffect.cpp

#include <string.h>
#include <math.h>
#include "mex.h"
#include "audioEffect.hpp"

/* Input Arguments */
#define	X		prhs[0]
#define	PRM		prhs[1]

/* Output Arguments */
#define	Y		plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	/* Check for proper number of arguments */
	if (nrhs!=2){
		char message[200];
		sprintf(message, "Format error!\nUsage: y = %s(x, prm)\n", mexFunctionName());
		mexErrMsgTxt(message);
	}

	int xLen=mxGetM(X)*mxGetN(X);
	double *x=mxGetPr(X);

	TDL_PRM prm;
	mxArray *fieldValue;
	fieldValue=mxGetField(PRM, 0, "delay");
	prm.delayNum=mxGetM(fieldValue)*mxGetN(fieldValue);
	prm.delay=(int *)malloc(prm.delayNum*sizeof(int));
	for (int i=0; i<prm.delayNum; i++)
		prm.delay[i]=(int)(mxGetPr(fieldValue)[i]);
	fieldValue=mxGetField(PRM, 0, "gain");
	prm.gain=(double *)malloc(prm.delayNum*sizeof(double));
	for (int i=0; i<prm.delayNum; i++)
		prm.gain[i]=mxGetPr(fieldValue)[i];

//	mexPrintf("prm.delayNum=%d\n", prm.delayNum);
//	for (int i=0; i<prm.delayNum; i++) mexPrintf("prm.delay[%d]=%d\n", i, prm.delay[i]);
//	for (int i=0; i<prm.delayNum; i++) mexPrintf("prm.gain[%d]=%f\n", i, prm.gain[i]);

	// Output matrix for pfMat
	Y = mxCreateDoubleMatrix(mxGetM(X), mxGetN(X), mxREAL);
	double *y=mxGetPr(Y);
	
	tapDelayLine(x, xLen, y, prm);
	
	// Free memory
	free(prm.delay);
	free(prm.gain);
}

