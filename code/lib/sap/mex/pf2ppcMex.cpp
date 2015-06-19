// SMDF to PPC
// How to compile:
// MATLAB 7.1: mex pf2ppcMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp -output pf2ppcMex.dll
// Others: mex pf2ppcMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp
#include <string.h>
#include <math.h>
#include "mex.h"
#include "audio.hpp"

/* Input Arguments */
#define	PF		prhs[0]
#define	FS		prhs[1]
#define	NBITS		prhs[2]
#define	PT_PARAM	prhs[2]
/* Output Arguments */
#define	LOCALOPTINDEX	plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double	*pf, *out;
	int i, j, fs, nbits, pfLen, localMinCount, *localOptIndex, *pfInt;
	cAudio myAudio;
	
	/* Check for proper number of arguments */
	if ((nrhs<3) || (nrhs>4)){
		char message[200];
		sprintf(message, "Format error!\nUsage 1: out = %s(pf, fs, nbits)\nUsage 2: out = %s(pf, fs, nbits, ptParam)", mexFunctionName(), mexFunctionName());
		mexErrMsgTxt(message);
	}

	/* Dimensions of the input matrix */
	pfLen = mxGetM(PF)*mxGetN(PF);
	
	/* Assign pointers to the various parameters */
	pf = mxGetPr(PF);
	fs = (int)mxGetScalar(FS);
	pfInt = (int *)malloc(pfLen*sizeof(int));
	double2intVec(pf, pfLen, pfInt);
	localOptIndex=(int *)malloc(pfLen*sizeof(int));
	for (i=0; i<pfLen; i++)
		localOptIndex[i]=0;

	PTPARAM ptParam; myAudio.ptParamSet(&ptParam, fs, nbits, AMDF);
	if (nrhs==4){
		mxArray *fieldValue;
		fieldValue=mxGetField(PT_PARAM, 0, "frameSize");	ptParam.frameSize=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "overlap");		ptParam.overlap=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "pfType");		ptParam.pfType=(PF_TYPE)(int)*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "pfWeight");		ptParam.pfWeight=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "indexDiffWeight");	ptParam.indexDiffWeight=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "dpMethod");		ptParam.dpMethod=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "dpPosDiffOrder");	ptParam.dpPosDiffOrder=*mxGetPr(fieldValue);
		fieldValue=mxGetField(PT_PARAM, 0, "dpUseLocalOptim");	ptParam.dpUseLocalOptim=*mxGetPr(fieldValue);
	}
//	mexPrintf("ptParam.ppcNum = %d\n", ptParam.ppcNum);

//	localMinCount=myAudio.pf2ppc(pfInt, pfLen, ptParam, localOptIndex);
	localMinCount=pf2ppc(pfInt, pfLen, ptParam, localOptIndex);

	/* Create a matrix for the return argument */
	LOCALOPTINDEX = mxCreateDoubleMatrix(localMinCount, 1, mxREAL);
	out = mxGetPr(LOCALOPTINDEX);
	for (i=0; i<localMinCount; i++)
		out[i]=(double)(localOptIndex[i]+1);	// MATLAB index
	
	free(pfInt);
	free(localOptIndex);
}
