#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	VOLUME		prhs[0]
#define	TESTNUM		prhs[0]
#define	DIVISOR		prhs[0]
/* Output Arguments */
#define	THRESHOLD	plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	int *intVolume, volumeLen, testNum, divisor;
	double *volume, *threshold;
	cAudio myAudio;
	
	/* Check for proper number of arguments */
	if (nrhs!=1) {
		char message[200];
		strcpy(message, mexFunctionName());
		strcat(message, " requires 1 input arguments.\n");
		strcat(message, "Usage: ");
		strcat(message, mexFunctionName());
		strcat(message, " volume");
		mexErrMsgTxt(message);
	}

	/* Assign pointers to the various parameters */
	volumeLen = mxGetM(VOLUME)*mxGetN(VOLUME);
	volume = mxGetPr(VOLUME);
	testNum=5;	// Default value
	divisor=4;	// Default value
	if (nrhs>1) testNum = mxGetScalar(TESTNUM);
	if (nrhs>2) divisor = mxGetScalar(DIVISOR);
	
	intVolume = (int *)malloc(volumeLen*sizeof(int));
	double2intVec(volume, volumeLen, intVolume);
	
	/* Create a matrix for the return argument */
	THRESHOLD = mxCreateDoubleMatrix(1, 1, mxREAL);
	threshold = mxGetPr(THRESHOLD);
	threshold[0]=myAudio.getVolumeThreshold(intVolume, volumeLen, testNum, divisor);
	
	free(intVolume);
}
