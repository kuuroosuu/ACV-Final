// How to compile:
// mex -I/users/jang/c/lib/audio lpCombMex.cpp /users/jang/c/lib/audio/audioEffect.cpp

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

	int xLen=(int)(mxGetM(X)*mxGetN(X));
	double *x=mxGetPr(X);

	LPCOMB_PRM prm;
	mxArray *fieldValue;
	fieldValue=mxGetField(PRM, 0, "delay");	prm.delay=(int)*mxGetPr(fieldValue);
	fieldValue=mxGetField(PRM, 0, "a");	prm.a=*mxGetPr(fieldValue);
	fieldValue=mxGetField(PRM, 0, "g");	prm.g=*mxGetPr(fieldValue);

//	mexPrintf("prm.delay=%d\n", prm.delay);
//	mexPrintf("prm.bl=%f\n", prm.bl);
//	mexPrintf("prm.fb=%f\n", prm.fb);
//	mexPrintf("prm.ff=%f\n", prm.ff);

	// Output matrix for pfMat
	Y = mxCreateDoubleMatrix((int)mxGetM(X), (int)mxGetN(X), mxREAL);
	double *y=mxGetPr(Y);
	
	lpComb(x, xLen, y, prm);
}

