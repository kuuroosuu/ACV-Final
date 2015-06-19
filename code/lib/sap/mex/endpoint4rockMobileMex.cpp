#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	WAVE		prhs[0]	// Wave data
#define	FS		prhs[1]	// Sample rate
#define	NBITS		prhs[2]	// No. of bits
/* Output Arguments */
#define	OUT	plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double	*wave, *out, sum;
	int m, n, i, j, fs, nbits, waveLen;
	char message[200];
	cAudio myAudio;

	/* Check for proper number of arguments */
	if (nrhs!=3) {
		strcpy(message, mexFunctionName());
		strcat(message, " requires 3 input arguments.\n");
		strcat(message, "Usage: [startIndex, endIndex] = ");
		strcat(message, mexFunctionName());
		strcat(message, "(wave, fs, nbits)");
		mexErrMsgTxt(message);
	}

	/* Dimensions of the input matrix */
	waveLen = mxGetM(WAVE)*mxGetN(WAVE);

	/* Create a matrix for the return argument */
	OUT = mxCreateDoubleMatrix(1, 2, mxREAL);

	/* Assign pointers to the various parameters */
	wave = mxGetPr(WAVE);
	fs = (int) mxGetScalar(FS);
	nbits = (int) mxGetScalar(NBITS);
	out = mxGetPr(OUT);

	myAudio.endPointDetect4rockMobile(wave, waveLen, fs, nbits, out);
}
