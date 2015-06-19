#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	Y1	prhs[0]
#define	FS1	prhs[1]
#define	FS2	prhs[2]

/* Output Arguments */
#define	Y2	plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *y1, *y2;
	int *intY1, *intY2;
	int fs1, fs2, len1, len2, i;
	cAudio myAudio;

	/* Check for proper number of arguments */
	if (nrhs!=3) {
		char message[200];
		strcpy(message, mexFunctionName());
		strcat(message, " requires 3 input arguments.\n");
		strcat(message, "Usage: y2 = ");
		strcat(message, mexFunctionName());
		strcat(message, "(y1, fs1, fs2)");
		mexErrMsgTxt(message);
	}

	/* Assign pointers to the various parameters */
	
	fs1 = (int)mxGetPr(FS1)[0];
	fs2 = (int)mxGetPr(FS2)[0];
	y1 = mxGetPr(Y1);
	len1 = mxGetM(Y1)*mxGetN(Y1);
	intY1 = (int *)malloc(len1*sizeof(int));
	double2intVec(y1, len1, intY1);
	
	len2 = (len1-1)*fs2/fs1+1;
	intY2 = (int *)malloc(len2*sizeof(int));

	/* Do the actual computation */
	myAudio.waveResample(intY1, len1, fs1, fs2, intY2, len2);

	/* Create a matrix for the return argument */
	Y2 = mxCreateDoubleMatrix(1, len2, mxREAL);
	y2 = mxGetPr(Y2);
	int2doubleVec(intY2, len2, y2);

	free(intY1);
	free(intY2);
}
