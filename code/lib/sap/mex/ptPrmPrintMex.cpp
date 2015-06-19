// Print pt (pitch tracking) parameters
// How to compile:
// MATLAB 7.1: mex ptParamPrintMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp -output ptParamPrintMex.dll
// Others: mex ptParamPrintMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp
#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	FS	prhs[0]
#define	NBITS	prhs[1]
/* Output Arguments */

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double	*wave, *amdfMat, *pitch;
	int fs, nbits;
	PTPARAM ptParam;
	cAudio myAudio;

	/* Check for proper number of arguments */
	if (nrhs!=2){
		char message[200];
		sprintf(message, "Usage: out = %s(fs, nbits)", mexFunctionName());
		mexErrMsgTxt(message);
	}

	fs    = (int)mxGetScalar(FS);
	nbits = (int)mxGetScalar(NBITS);
	myAudio.ptParamSet(&ptParam, fs, nbits, AMDF);
	myAudio.ptParamPrint(&ptParam);
}
