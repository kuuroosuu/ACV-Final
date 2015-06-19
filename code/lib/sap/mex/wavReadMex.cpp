// How to compile:
// mex -Id:/users/jang/c/lib -Id:/users/jang/c/lib/wave -Id:/users/jang/c/lib/utility wavReadMex.cpp d:/users/jang/c/lib/wave/CWavePCM.cpp d:/users/jang/c/lib/utility/utility.cpp -output wavReadMex.dll

#include <string.h>
#include <math.h>
#include "mex.h"
#include "CWavePCM.h"
#include <iostream>

using namespace std;

/* Input Arguments */
#define	WAVEFILENAME		prhs[0]

/* Output Arguments */
#define	Y	plhs[0]
#define	FS	plhs[1]
#define	NBITS	plhs[2]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	/* Check for proper number of arguments */
	if ((nrhs<1) || (nrhs>1)){
		char message[200];
		sprintf(message, "Format error!\nUsage 1: [y, fs, nbits] = %s(waveFile)", mexFunctionName(), mexFunctionName());
		mexErrMsgTxt(message);
	}
	
	if (mxIsChar(WAVEFILENAME)!=1)
		mexErrMsgTxt("Input must be a string.");

	/* copy the string data from prhs[0] into a C string waveFileName. */
	char *waveFileName = mxArrayToString(WAVEFILENAME);
	if (waveFileName==NULL) 
		mexErrMsgTxt("Could not convert input to string.");

	CWavePCM wvWave;
	if (wvWave.ReadFromWaveFile(waveFileName)==false){
		cout << "Cannot read " << waveFileName << "!" << endl;
		exit(1);
	}
	int fs=wvWave.QuerySamplesPerSec();
	int nbits=wvWave.QueryBitsPerSample();
	int waveLen=wvWave.QueryNumOfWaveSamples();
	
	// 把 wave 中的 sample 轉成 int 取出來
	int *wave = new int[waveLen];
	wvWave.GetWaveData(wave, waveLen);
	
	// Output matrix
	Y = mxCreateDoubleMatrix(waveLen, 1, mxREAL);
	int maximum=pow(2.0, nbits-1);
	int shift=(nbits==8)?128:0;
	for (int i=0; i<waveLen; i++)
		mxGetPr(Y)[i]=(wave[i]-shift)/(double)maximum;
	FS = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxGetPr(FS)[0]=fs;
	NBITS = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxGetPr(NBITS)[0]=nbits;
	
	// Free memory
	delete [] wave;
	mxFree(waveFileName);
}
