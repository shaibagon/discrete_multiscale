#include "mex.h"
#include <string.h>
#include <math.h>
#include <vector>
/*
 * Usage:
 *   l = gen_eng_iter_mex(U, V, W, l);
 *
 * Solve iteratively for 
 *   E(l) = min_l sum_i U(i,l_i) + sum_ij w_ij V(l_i, l_j)
 *
 *
 * Use Gauss - Seidel iterations (i.e., use updated values in current iteration)
 *
 * U is double LxN real matrix (unary term)
 * V is double LxL real matrix (invariant smothness), assumed symmteric
 * W is sparse NxN real matrix 
 * l is double Nx1 vector of init labels 
 *
 *
 * compile using:
 * >> mex -O -largeArrayDims gen_eng_icm_iter_mex.cpp
 */

#line   __LINE__  "gen_eng_icm_iter_mex"

#define     STR(s)      #s  
#define     ERR_CODE(a,b)   a ":" "line_" STR(b)


// INPUTS
enum {
    UIN = 0,
    VIN,
    WIN,
    LIN,
    NIN 
};

// OUTPUTS
enum {
    LOUT = 0,
    NOUT
};

template<typename T>
inline
void plus_times(T* buffer, T wij, const T* vlj, mwSize L)
{
    
    for (mwSize ii(0) ; ii < L ; ii++ ) {
        buffer[ii] += wij*vlj[ii];
    }
}
template<typename T>
inline
mwIndex AdjustCol(T* buffer, mwSize L)        
{
    // what is the maximal entry?
    mwIndex mi = 1; // 1-based matlab index
    T mn = buffer[0];
    
    for (mwSize ii(1) ; ii < L ; ii++ ) {
        if ( buffer[ii] < mn ) {
            mn = buffer[ii];
            mi = ii+1; // 1-based matlab index/label
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

    if ( mxIsComplex(pin[VIN]) || mxIsSparse(pin[VIN]) ||
            mxGetClassID(pin[VIN]) != mxDOUBLE_CLASS ||
            mxGetNumberOfDimensions(pin[VIN]) != 2 )
        mexErrMsgIdAndTxt(ERR_CODE(__FILE__, __LINE__), "V must be real full double matrix");

    if ( mxIsComplex(pin[WIN]) || !mxIsSparse(pin[WIN]) )
        mexErrMsgIdAndTxt(ERR_CODE(__FILE__, __LINE__),"W must be real sparse matrix");
    
    if ( mxIsComplex(pin[LIN]) || mxIsSparse(pin[LIN]) || mxGetClassID(pin[LIN]) != mxDOUBLE_CLASS )
        mexErrMsgIdAndTxt(ERR_CODE(__FILE__, __LINE__),"l must be real full double vector");
    
    // sizes
    mwSize N = mxGetN(pin[WIN]); // number of variables
    mwSize L = mxGetM(pin[UIN]); // number of labels
    
    if ( mxGetN(pin[WIN]) != mxGetM(pin[WIN]) || 
            mxGetN(pin[WIN]) != N || mxGetNumberOfDimensions(pin[WIN])!=2 ||
            mxGetNumberOfElements(pin[LIN]) != N || 
            mxGetN(pin[UIN]) != N ||
            mxGetM(pin[UIN]) != L ||
            mxGetN(pin[VIN]) != L || 
            mxGetM(pin[VIN]) != L ||
            mxGetNumberOfDimensions(pin[VIN])!=2)
        mexErrMsgIdAndTxt(ERR_CODE(__FILE__, __LINE__),"matrix - vector dimensions mismatch");

    // allocate output vector   
    
    pout[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    double* pl = mxGetPr(pout[0]);
    
    double* pol = mxGetPr(pin[LIN]);
    memcpy(pl, pol, N*sizeof(double));
    
   
    
    /* computation starts */
    double* pu = mxGetPr(pin[UIN]);
    double* pv = mxGetPr(pin[VIN]);
    
    double* pr = mxGetPr(pin[WIN]);
    mwIndex* pir = mxGetIr(pin[WIN]);
    mwIndex* pjc = mxGetJc(pin[WIN]);
 
    
    double* buffer=new double[L]; // start a vector the size of L
    // foreach variavl i (col)
    // compute buffer[l_i] = E(i=l_i) = U_i(l_i) + sum_j w_ij V(l_i, l_j)
    
    // compute each column of u' and then make a decision about it 
    // before moving on to the next columns
    for (mwSize col=0; col< N; col++)  {
        
        // copy unary term to buffer
        memcpy(buffer, pu+col*L, L*sizeof(double));
        
                
        // perform sparse multiplication
        for (mwIndex ri = pjc[col] ; // starting row index
        ri < pjc[col+1]  ; // stopping row index
        ri++)  {
            if ( col != pir[ri] ) {
                // only off-diagonal elements are participating
                if ( pl[ pir[ri] ] >= 1 ) {
                    plus_times(buffer, pr[ri] /*w_ij*/,
                            pv + static_cast<mwSize>(pl[ pir[ri] ]-1/*L is 1-based matlab index */)*L,
                            L);
                }
            }
            // pir[ri] -> current row
            // col -> current col
            // pr[ri] -> W[pir[ri], col]
            
        }
        
        pl[col] = AdjustCol(buffer, L); // make a decision for this pixel
        
    }
    delete[] buffer;
}

