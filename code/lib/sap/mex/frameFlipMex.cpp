#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	FRAME	prhs[0]
/* Output Arguments */
#define	FRAME2	plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double	*frame, *frame2, *amdfMat, *pitch;
	int i, j, k, m, n, frameSize, overlap, amdfLen, waveLen, frameNum, fs, nbits, *frameInt, *amdf, *waveInt, *framePos, *amdfMatInt;
	cAudio myAudio;

	/* Check for proper number of arguments */
	if (nrhs!=1){
		char message[200];
		sprintf(message, "Usage: out = %s(frame)", mexFunctionName());
		mexErrMsgTxt(message);
	}

	// Input
	m=mxGetM(FRAME);
	n=mxGetN(FRAME);
	frameSize=m*n;
	frame=mxGetPr(FRAME);
	
	// Output
	FRAME2=mxCreateDoubleMatrix(m, n, mxREAL);
	frame2=mxGetPr(FRAME2);
	memcpy(frame2, frame, m*n*sizeof(double));
	
	// Flip it!
	myAudio.frameFlip(frame2, frameSize);
}
