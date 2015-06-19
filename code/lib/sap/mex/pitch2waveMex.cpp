#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	PITCH		prhs[0]
#define PITCHRATE	prhs[1]
#define FS		prhs[2]
/* Output Arguments */
#define	WAVE		plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double	*pitch, *wave, pitchRate;
	int fs, pitchLen, waveLen;
	cAudio myAudio;

	/* Check for proper number of arguments */
	if (nrhs!=3){
		char message[200];
		strcpy(message, mexFunctionName());
		strcat(message, " requires 3 input arguments: pitch, pitchRate, fs.");
		mexErrMsgTxt(message);
	}

	/* Dimensions of the input matrix */
	pitchLen = mxGetM(PITCH)*mxGetN(PITCH);
	
	/* Assign pointers to the various parameters */
	pitch = mxGetPr(PITCH);
	pitchRate = mxGetScalar(PITCHRATE);
	fs = (int)mxGetScalar(FS);
	
	/* If some of the elements are nan/inf, change them to 0 */
	double *pitch2=(double *)malloc(pitchLen*sizeof(double));
	memcpy(pitch2, pitch, pitchLen*sizeof(double));
	for (int i=0; i<pitchLen; i++)
		if (mxIsNaN(pitch[i])||mxIsInf(pitch[i]))
			pitch2[i]=0.0;
	
	/* Create a matrix for the return argument */
	waveLen=(int)floor(fs*pitchLen/pitchRate+0.5);
	WAVE = mxCreateDoubleMatrix(waveLen, 1, mxREAL);
	wave = mxGetPr(WAVE);

	myAudio.pitch2wave(pitch2, pitchLen, pitchRate, fs, wave, waveLen);
	free(pitch2);
}
