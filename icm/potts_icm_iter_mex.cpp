#include "mex.h"
#include <string.h>
#include <math.h>
#include <vector>
/*
 * Usage:
 *   l = potts_icm_iter_mex(U, W, l);
 *
 * Solve iteratively for 
 *   E(l) = sum_i U(i,l_i) + sum_ij -w_ij [l_i = l_j]
 *
 * Where [.] is 1 if the argument is true
 *
 * Use Gauss - Seidel iterations (i.e., use updated values in current iteration)
 *
 * U is double NxL real matrix (unary term)
 * W is sparse NxN real matrix 
 * l is double Nx1 vector of labels 
 *
 *
 * compile using:
 * >> mex -O -largeArrayDims potts_icm_iter_mex.cpp
 */

#line   __LINE__  "potts_icm_iter_mex"

#define     STR(s)      #s  
#define     ERR_CODE(a,b)   a ":" "line_" STR(b)


// INPUTS
enum {
    UIN = 0,
    WIN,
    LIN,
    NIN 
};

// OUTPUTS
enum {
    LOUT = 0,
    NOUT
};
// template<typename T>
// inline
// T min(const T& a, const T& b)
// {
//     return (a<b)?a:b;
// }

template<typename T>
inline
mwIndex AdjustCol(const typename std::vector<T>& buffer)        
{
    // what is the maximal entry?
    mwIndex mi = 1; // 1-based matlab index
    T mn = buffer[0];
    for (mwIndex k = 1 ; k < buffer.size() ; k++ ) {
        if ( buffer[k] < mn ) {
            mn = buffer[k];
            mi = k+1; // 1-based matlab index/label
        }
    }  
    return mi;
}

void
mexFunction(
    int nout,
    mxArray* pout[],
    int nin,
    const mxArray* pin[])
{
    if ( nin != NIN )
         mexErrMsgIdAndTxt(ERR_CODE(__FILE__, __LINE__),"must have %d inputs", NIN);
    if (nout==0)
        return;
    if (nout != NOUT )
         mexErrMsgIdAndTxt(ERR_CODE(__FILE__, __LINE__),"must have exactly %d output", NOUT);
    
    if ( mxIsComplex(pin[UIN]) || mxIsSparse(pin[UIN]) || 
            mxGetClassID(pin[UIN]) != mxDOUBLE_CLASS || 
            mxGetNumberOfDimensions(pin[UIN]) != 2 )
        mexErrMsgIdAndTxt(ERR_CODE(__FILE__, __LINE__),"U must be real full double matrix");

    if ( mxIsComplex(pin[WIN]) || !mxIsSparse(pin[WIN]) )
        mexErrMsgIdAndTxt(ERR_CODE(__FILE__, __LINE__),"W must be real sparse matrix");
    
    if ( mxIsComplex(pin[LIN]) || mxIsSparse(pin[LIN]) || mxGetClassID(pin[LIN]) != mxDOUBLE_CLASS )
        mexErrMsgIdAndTxt(ERR_CODE(__FILE__, __LINE__),"l must be real full double vector");
    
    // sizes
    mwSize N = mxGetN(pin[WIN]); // number of variables
    mwSize L = mxGetN(pin[UIN]); // number of labels
    
    if ( mxGetN(pin[WIN]) != mxGetM(pin[WIN]) || 
            mxGetN(pin[WIN]) != N || mxGetNumberOfDimensions(pin[WIN])!=2 ||
            mxGetNumberOfElements(pin[LIN]) != N || 
            mxGetM(pin[UIN]) != N )
        mexErrMsgIdAndTxt(ERR_CODE(__FILE__, __LINE__),"matrix - vector dimensions mismatch");

    // allocate output vector   
    
    pout[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    double* pl = mxGetPr(pout[0]);
    
    double* pr = mxGetPr(pin[LIN]);
    memcpy(pl, pr, N*sizeof(double));
    
   
    
    /* computation starts */
    double *pu = mxGetPr(pin[UIN]);
    
    pr = mxGetPr(pin[WIN]);
    mwIndex* pir = mxGetIr(pin[WIN]);
    mwIndex* pjc = mxGetJc(pin[WIN]);
 
    // compute each column of u' and then make a decision about it 
    // before moving on to the next columns
    for (mwSize col=0; col< N; col++)  {
        
        std::vector<double> buffer(L); // start a vector the size of L
        for ( mwIndex l(0) ; l < L ; l++ )
            buffer[l] = pu[col + N*l];        
            
        
        // perform sparse multiplication
        for (mwIndex ri = pjc[col] ; // starting row index
        ri < pjc[col+1]  ; // stopping row index
        ri++)  {
            if ( col != pir[ri] ) {
                // only off-diagonal elements are participating
                buffer[ pl[ pir[ri] ]-1 ] -= pr[ri];
            }
            // pir[ri] -> current row
            // col -> current col
            // pr[ri] -> W[pir[ri], col]
            
        }
        pl[col] = AdjustCol(buffer); // make a decision for this pixel
        
    }
}

