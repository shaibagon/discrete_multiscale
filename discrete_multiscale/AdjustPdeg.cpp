#include "mex.h"
#include <vector>
#include <algorithm>

/*
 *
 * Adjust COLUMNS of P to have only maxD non zero elements
 *
 * aP = AdjustPdeg(P, maxD);
 *
 * Compile
 * >> mex -O -largeArrayDims AdjustPdeg.cpp
 */

#line   __LINE__  "AdjustPdeg"

#define     STR(s)      #s  
#define     ERR_CODE(a,b)   a ":" "line_" STR(b)


// INPUTS
enum {
    PIN = 0,
    MIN,
    NIN 
};

// OUTPUTS
enum {
    POUT = 0,
    NOUT
};

class ValIndPair {
public:
    mwIndex i;
    double  v;
    
    ValIndPair():i(0),v(0) {}
    ValIndPair(mwIndex i_, double v_):i(i_),v(v_) {}
};

bool operator>(const ValIndPair& rhs, const ValIndPair& lhs)
{
    return rhs.v>lhs.v;
}

bool operator>=(const ValIndPair& rhs, const ValIndPair& lhs)
{
    return rhs.v>=lhs.v;
}

bool operator==(const ValIndPair& rhs, const ValIndPair& lhs)
{
    return rhs.v==lhs.v;
}

bool operator!=(const ValIndPair& rhs, const ValIndPair& lhs)
{
    return rhs.v!=lhs.v;
}

bool operator<(const ValIndPair& rhs, const ValIndPair& lhs)
{
    return rhs.v<lhs.v;
}

bool operator<=(const ValIndPair& rhs, const ValIndPair& lhs)
{
    return rhs.v<=lhs.v;
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
    
    if ( mxIsComplex(pin[PIN]) || !mxIsSparse(pin[PIN]) || 
            mxGetClassID(pin[PIN]) != mxDOUBLE_CLASS || 
            mxGetNumberOfDimensions(pin[PIN]) != 2 )
        mexErrMsgIdAndTxt(ERR_CODE(__FILE__, __LINE__),"P must be real sparse double matrix");
    
    if ( mxGetNumberOfElements(pin[MIN]) != 1 || mxIsComplex(pin[MIN]) )
        mexErrMsgIdAndTxt(ERR_CODE(__FILE__, __LINE__),"maxD must be a real scalar");
    
    
    // max number of enelemnts
    mwIndex maxD = static_cast<mwIndex>(mxGetScalar(pin[MIN]));
    
    // input P dimentions
    mwSize C = mxGetN(pin[PIN]);
    mwSize R = mxGetM(pin[PIN]);
    
    double* pr = mxGetPr(pin[PIN]);
    mwIndex* pir = mxGetIr(pin[PIN]);
    mwIndex* pjc = mxGetJc(pin[PIN]);
 
    
    // allocate output
    pout[POUT] = mxCreateSparse(R, C, maxD*C /*only maxD nonzeros per column*/, mxREAL);
    double* opr = mxGetPr(pout[POUT]);
    mwIndex* opir = mxGetIr(pout[POUT]);
    mwIndex* opjc = mxGetJc(pout[POUT]);
 
    // process matrix
    mwIndex ori = 0;
    for ( mwIndex col(0); col < C ; col++ ) {
        std::vector<ValIndPair> vector(0); 
        
        // push into heap
        for (mwIndex ri = pjc[col] ; // starting row index
            ri < pjc[col+1]  ; // stopping row index
            ri++)  {
                vector.push_back( ValIndPair(pir[ri], -pr[ri]) ); /* v in minus because partial sort is in ascending order and not descending... */
            }
            // pir[ri] -> current row
            // col -> current col
            // pr[ri] -> W[pir[ri], col]           
               
//        mexPrintf("col %d\tGot %d elements\n", col, vector.size());
        
        opjc[col] = ori;
        mwIndex M = vector.size();
        if ( vector.size() > maxD ) {
            M = maxD;
            partial_sort(vector.begin(), vector.begin() + M, vector.end() );
        }
        
        // pop only the maxD maximal elements maximal
        for ( mwIndex mi(0) ; mi < M ; mi++ ) {            
            ValIndPair curr = vector[mi];
            
//            mexPrintf("\t\t poped v=%f, i=%d\n", -curr.v, curr.i); /* v in minus because partial sort is in ascending order and not descending... */
            
            opr[ori] = -curr.v; /* v in minus because partial sort is in ascending order and not descending... */
            opir[ori] = curr.i;
            
            ori++;
        }
    }
    opjc[C] = ori;    
    
        
}
    
    
