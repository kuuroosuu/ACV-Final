#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	NOTEPITCH	prhs[0]
#define NOTEDURATION	prhs[1]
#define FS		prhs[2]
/* Output Arguments */
#define	WAVE		plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double	*notePitch, *noteDuration, *wave;
	int fs, noteLen, waveLen;
	cAudio myAudio;

	/* Check for proper number of arguments */
	if (nrhs!=3){
		char message[200];
		strcpy(message, mexFunctionName());
		strcat(message, " requires 3 input arguments: notePitch, noteDuration, fs.");
		mexErrMsgTxt(message);
	}

	/* Dimensions of the input matrix */
	noteLen = mxGetM(NOTEPITCH)*mxGetN(NOTEPITCH);
	
	/* Assign pointers to the various parameters */
	notePitch = mxGetPr(NOTEPITCH);
	noteDuration = mxGetPr(NOTEDURATION);
	fs = (int)mxGetScalar(FS);
	
	/* Create a matrix for the return argument */
	waveLen=(int)(fs*vecSum(noteDuration, noteLen));
//	printf("sum=%f\n", vecSum(noteDuration, noteLen));
//	printf("waveLen=%d\n", waveLen);
	WAVE = mxCreateDoubleMatrix(waveLen, 1, mxREAL);
	wave = mxGetPr(WAVE);
	
	myAudio.note2wave(notePitch, noteLen, noteDuration, fs, wave, waveLen);
}
