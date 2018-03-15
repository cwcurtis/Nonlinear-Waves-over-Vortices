/*********************************************************************
 * qf_comp.cpp
 *
 * This file computes the associated potentials in the FMM method.
 * 
 * Keep in mind:
 * <> Use 0-based indexing as always in C or C++
 * <> Indexing is column-based as in Matlab (not row-based as in C)
 * <> Use linear indexing.  [x*dimy+y] instead of [x][y]
 *
 *
 ********************************************************************/
#include <matrix.h>
#include <mex.h>   
#include <complex.h>

/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//declare variables
    mxArray *zinr, *zini, *kinr, *kini, *qoutr, *qouti;
    const mwSize *dims, *dimsp;
    double *zvalsr, *zvalsi, *kcr, *kci, *qvalsr, *qvalsi;
    double zvr, zvi, qvr, qvi;
    int dimz, dimp, ndim, ndimp;
    int i,j;

//associate inputs
    zinr = mxDuplicateArray(prhs[0]);
    zini = mxDuplicateArray(prhs[1]);
    kinr = mxDuplicateArray(prhs[2]);
    kini = mxDuplicateArray(prhs[3]);

//figure out dimensions
    ndim = mxGetNumberOfDimensions(prhs[0]);
    ndimp = mxGetNumberOfDimensions(prhs[2]);
    dims = mxGetDimensions(prhs[0]);
    dimsp = mxGetDimensions(prhs[2]);
    dimz = (int)dims[0]; dimp = (int)dimsp[0];

//associate outputs
    qoutr = plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    qouti = plhs[1] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    
//associate pointers
    zvalsr = mxGetPr(zinr);
    zvalsi = mxGetPr(zini);
    kcr = mxGetPr(kinr);
    kci = mxGetPr(kini);
    qvalsr = mxGetPr(qoutr);
    qvalsi = mxGetPr(qouti);    
//do something
    
    for(i=0;i<dimz;i++)
    {
        zvr = zvalsr[i];
        zvi = zvalsi[i];
        qvalsr[i] = kcr[dimp-2] + kcr[dimp-1]*zvr - kci[dimp-1]*zvi;
        qvalsi[i] = kci[dimp-2] + kcr[dimp-1]*zvi + kci[dimp-1]*zvr;
               
        for(j=2;j<dimp;j++)
        {
            qvr = qvalsr[i];
            qvi = qvalsi[i];
            qvalsr[i] = kcr[dimp-1-j] + zvr*qvr-zvi*qvi;             
            qvalsi[i] = kci[dimp-1-j] + zvr*qvi+zvi*qvr;             
        }
        qvr = qvalsr[i];
        qvi = qvalsi[i];
        qvalsr[i] = zvr*qvr - zvi*qvi;
        qvalsi[i] = zvr*qvi + zvi*qvr;
    }

    return;
}