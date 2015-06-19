#include <string.h>
#include <math.h>
#include "mex.h"
#include "utility.hpp"

/* Input Arguments */
#define	IN	prhs[0]
/* Output Arguments */
#define	OUT	plhs[0]

void mexFunction(
	int nlhs,	mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double	*in, *out;
	int m, n, i, j;
	double sum;
	char message[200];

	/* Check for proper number of arguments */
	if (nrhs==0) {
		strcpy(message, mexFunctionName());
		strcat(message, " requires one input arguments.");
		mexErrMsgTxt(message);
	}

	/* Dimensions of the input matrix */
	m = mxGetM(IN);
	n = mxGetN(IN);

	/* Create a matrix for the return argument */
	OUT = mxCreateDoubleMatrix(m, n, mxREAL);

	/* Assign pointers to the various parameters */
	out = mxGetPr(OUT);
	in = mxGetPr(IN);

	memcpy(out, in, m*n*sizeof(double));

	/* Do the column-wise computations */
	for (j=0; j<n; j++)
		vecZeroMean(out+j*m, m);
}
