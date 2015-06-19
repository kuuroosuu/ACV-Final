#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	VEC1	prhs[0]

/* Output Arguments */
#define	VEC2	plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *vec1, *vec2;
	int i, *intVec, size;
	cAudio myAudio;
	
	/* Check for proper number of arguments */
	if (nrhs!=1) {
		char message[200];
		sprintf(message, "Usage: vec2 = %s(vec1)", mexFunctionName());
		mexErrMsgTxt(message);
	}
	
	/* Assign pointers to the various parameters */
	vec1 = mxGetPr(VEC1);
	size=mxGetM(VEC1)*mxGetN(VEC1);
	intVec = (int *)malloc(size*sizeof(int));
	double2intVec(vec1, size, intVec);
	
	myAudio.localAverage(intVec, size);
	
	/* Create a matrix for the return argument */
	VEC2 = mxCreateDoubleMatrix(mxGetM(VEC1), mxGetN(VEC1), mxREAL);
	vec2 = mxGetPr(VEC2);
	
	for (i=0; i<size; i++)
		vec2[i]=intVec[i];

	free(intVec);
}
