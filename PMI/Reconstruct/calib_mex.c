//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%  $Author: dboas $
//%
//%  $Date: 2000/05/25 13:14:47 $
//%
//%  $Revision: 1.1.1.1 $
//%
//%  $Log: calib_mex.c,v $
//%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
//%  initial
//%
//%  Revision 1.2  2000/01/10 00:14:14  dboas
//%  Storing the source and detector lists for use by other functions
//%
//%  Revision 1.1  1999/12/03 14:03:02  dboas
//%  The MEX file to handle the case when measurements are not
//%  made between every source and every detector.
//%
//%  Revision 1.0  1999/12/03 14:00:28  dboas
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Generic Outline for Mex files 
// Place your c function in the top of this file after include statements
// Then call your function with the appropriet function call inside the 
// 	generic mexFunction function
//

#include "mex.h"
#include "matrix.h"
#include <stdio.h>

// Add any include files need for the c-code here

// Add you C-Function here.....
void MyFunction(double *g, double *g_std, double *f, long nSrcs, long nDets, double *coef, double *lsq )
{
	int i, j, k, l, ii, jj;
	int nD, nS;
	double As, Ad;
	double foo;
	double *L;
	double numer, denom;
	double *g_std2;

	L = (double *)malloc(nDets*nSrcs * sizeof(double));
	g_std2 = (double *)malloc(nDets*nSrcs * sizeof(double));

	*lsq = 0;
	for( l=0; l<nSrcs; l++ ) {
		for( k=0; k<nDets; k++ ) {

			coef[l+k*nSrcs] = 0;
			if( f[l+k*nSrcs]!=0 ) {
				for( i=0; i<nSrcs; i++ ) {
					for( j=0; j<nDets; j++ ) {
					
						// Calculate A_s^{ik}
						As = 0;
						nD = 0;
						for( jj=0; jj<nDets; jj++ ) {
							if( f[i+jj*nSrcs]!=0 ) {
								foo = 0;
								nS = 0;
								nD++;
								for( ii=0; ii<nSrcs; ii++ ) {
									if( f[ii+jj*nSrcs]!=0 && f[ii+k*nSrcs]!=0 ) {
										nS++;
										foo += g[ii+jj*nSrcs] * f[ii+k*nSrcs] / (f[ii+jj*nSrcs] * g[ii+k*nSrcs]);
									}
								}
								As += nS * g[i+jj*nSrcs] / (f[i+jj*nSrcs] * foo);
							}
						}
						As /= nD;

						// Calculate A_d^{jl}
						Ad = 0;
						nS = 0;
						for( ii=0; ii<nSrcs; ii++ ) {
							if( f[ii+j*nSrcs]!=0 ) {
								foo = 0;
								nD = 0;
								nS++;
								for( jj=0; jj<nDets; jj++ ) {
									if( f[ii+jj*nSrcs]!=0 && f[l+jj*nSrcs]!=0 ) {
										nD++;
										foo += g[ii+jj*nSrcs] * f[l+jj*nSrcs] / (f[ii+jj*nSrcs] * g[l+jj*nSrcs]);
									}
								}
								Ad += nD * g[ii+j*nSrcs] / (f[ii+j*nSrcs] * foo);
							}
						}
						Ad /= nS;
	
						L[i+j*nSrcs] = As * Ad;
					}   // end loop on j
				}	// end loop on i

				for( i=0; i<nSrcs; i++ ) {
					for( j=0; j<nDets; j++ ) {
						if( f[i+j*nSrcs]!=0 ) {
							g_std2[i+j*nSrcs] = f[i+j*nSrcs] * f[i+j*nSrcs];
//g_std[i+j*nSrcs] * g_std[i+j*nSrcs];
						}
					}
				}

				numer = 0;
				denom = 0;
				for( i=0; i<nSrcs; i++ ) {
					for( j=0; j<nDets; j++ ) {
						if( f[i+j*nSrcs]!=0 ) {
							numer += g[i+j*nSrcs] * L[i+j*nSrcs] * f[i+j*nSrcs] / g_std2[i+j*nSrcs];
							denom += L[i+j*nSrcs] * L[i+j*nSrcs] * f[i+j*nSrcs] * f[i+j*nSrcs] / g_std2[i+j*nSrcs];
						}
					}
				}
				coef[l+k*nSrcs] = denom / numer;

				for( i=0; i<nSrcs; i++ ) {
					for( j=0; j<nDets; j++ ) {
						if( f[i+j*nSrcs]!=0 ) {
							*lsq += (L[i+j*nSrcs]*f[i+j*nSrcs]/coef[l+k*nSrcs] - g[i+j*nSrcs])*(L[i+j*nSrcs]*f[i+j*nSrcs]/coef[l+k*nSrcs] - g[i+j*nSrcs]) / g_std2[i+j*nSrcs];
						}
					}
				}				
			}

		} // end loop on k
	} // end loop on l

	free(L);
	free(g_std2);
}




// Use this if you need to find out what is the class of the variable passed by
//	Matlab.
mxClassID determine_class(mxArray  *array_ptr){ mxClassID  category;
 /* Display the class name of each input mxArray. */ 
   category = mxGetClassID(array_ptr);   switch (category) {
     case mxCHAR_CLASS: mexPrintf("String ");       break;
     case mxSTRUCT_CLASS:mexPrintf("Structure ");    break;
     case mxSPARSE_CLASS:mexPrintf("Sparse ");       break; 
     case mxCELL_CLASS: mexPrintf("Cell ");         break; 
     case mxUNKNOWN_CLASS:mexWarnMsgTxt("Unknown Class.");          break;
     default:           mexPrintf("Full Numeric 'Double' "); break;    } }

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
   long nDets, nSrcs;
   long nDetsf, nSrcsf;
   double *g, *g_std, *f;
  	mxArray *answer_ptr;
  	mxArray *answer_ptr2;
	double *coef, *lsq;
	int status;
 
  if( nrhs!=3 ) {
    mexErrMsgTxt( "Usage: calib_mex( g, g_std, f )\n\twhere g is a matrix of measurements and f is a matrix of calculated fluences.\n" );
  }
 
 
/* Start of Input Processing */	
	nDets = mxGetN(prhs[0]);
   nSrcs = mxGetM(prhs[0]);
   g = (double *) mxGetPr(prhs[0]);
   
   g_std = (double *) mxGetPr(prhs[1]);

	nDetsf = mxGetN(prhs[2]);
   nSrcsf = mxGetM(prhs[2]);
	if( nDetsf!=nDets || nSrcsf!=nSrcs ) {
		mexErrMsgTxt( "The dimensions of matrix g and matrix f must match.\n" );
	}
   f = (double *) mxGetPr(prhs[2]);


// 	Setup where we will return results
   answer_ptr = mxCreateDoubleMatrix(nSrcs,nDets,mxREAL );
   coef=(double *)mxCalloc(nDets*nSrcs,mxGetElementSize(answer_ptr));
   answer_ptr2 = mxCreateDoubleMatrix(1,1,mxREAL );
   lsq=(double *)mxCalloc(1,mxGetElementSize(answer_ptr2));

//	Call C-Function
	MyFunction( g, g_std, f, nSrcs, nDets, coef, lsq );

//	Return results in either the left-hand arg. or standard 'ans' variable
/* set phi to the output given by fdtd c-function */
   	mxSetM(answer_ptr, nSrcs);
   	mxSetN(answer_ptr, nDets);
   	mxSetPr(answer_ptr, coef);
   	mxSetPi(answer_ptr, NULL);

   	mxSetM(answer_ptr2, 1);
   	mxSetN(answer_ptr2, 1);
   	mxSetPr(answer_ptr2, lsq);
   	mxSetPi(answer_ptr2, NULL);

/* Examine output (left-hand-side) arguments. */
//If no output arg used in call we want to create the default 'ans' variable;
	if(nlhs<=0)
    	{
   		mxSetName(answer_ptr, "ans");
		/* Put "ans" into the MATLAB workspace or who ever called us*/
   		status = mexPutArray(answer_ptr, "caller");
   		if(status!=0)
			mexErrMsgTxt("Problem returning answer");

 		/* Clean up temporary array and exit */
   		mxSetPr(answer_ptr, NULL);
   		mxDestroyArray(answer_ptr);
   }
   else
   {
    		/* Return the answer */
     	plhs[0] = answer_ptr; 
     	plhs[1] = answer_ptr2; 
   }
 
}
/* END */
