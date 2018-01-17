/*
//mmex INCLUDES = { fileparts( which( 'vec.m' ) ) };
*/

#include "myMEX.h"
#include "distances.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  int    NSD, k;

  double *A0 = NULL;
  double *A1 = NULL;
  double *B0 = NULL;
  double *B1 = NULL;
  
  double a[1], b[1], d2;
  
  int    aa, bb, nVA, nVB, nFA, nFB;
  double *VA, *VB, *FA, *FB;
  double *out_d, *out_a, *out_b;
  
  mxArray  *mxVA, *mxVB, *mxFA, *mxFB;
  int AisClosed = 0;
  int BisClosed = 0;
  int A0id, A1id, B0id, B1id;

  if( nrhs != 2 ){
    myErrMsgTxt("2 input arguments were expected.");
  }

  
  mxFA = NULL; AisClosed = 0;
  if( mxIsCell( prhs[0] ) ){
    if( mxGetNumberOfElements( prhs[0] ) == 1 ){
      mxVA = mxGetCell( prhs[0] , 0 );
      AisClosed = 1;
    } else if( mxGetNumberOfElements( prhs[0] ) == 2 ){
      mxVA = mxGetCell( prhs[0] , 0 );
      mxFA = mxGetCell( prhs[0] , 1 );
      if( mxGetNumberOfElements( mxFA ) == 0 ){ mxFA = NULL; AisClosed = 1; }
    } else {
      myErrMsgTxt("As a cell a pair vertices - faces was expected for the first argument.");
    }
  } else if( mxIsStruct( prhs[0] ) ){
    if( mxGetFieldNumber( prhs[0] ,"xyz" ) >= 0 ){ mxVA = mxGetField( prhs[0] ,0,"xyz" );
    } else if( mxGetFieldNumber( prhs[0] ,"vertices" ) >= 0 ){ mxVA = mxGetField( prhs[0] ,0,"vertices" );
    } else { myErrMsgTxt("As a struct a \".xyz\" field is compulsory for the first argument.");
    }
    if( mxGetFieldNumber( prhs[0] ,"tri" ) > 0 ){ mxFA = mxGetField( prhs[0] ,0,"tri" );
    } else if( mxGetFieldNumber( prhs[0] ,"faces" ) > 0 ){ mxFA = mxGetField( prhs[0] ,0,"faces" );
    }
  } else {
    mxVA = prhs[0];
  }
  if( mxFA != NULL && mySize( mxFA , 1 ) != 2 ){ mexErrMsgTxt("Segments for the first polygon were expected as an Sx2 matrix."); }


  mxFB = NULL; BisClosed = 0;
  if( mxIsCell( prhs[1] ) ){
    if( mxGetNumberOfElements( prhs[1] ) == 1 ){
      mxVB = mxGetCell( prhs[1] , 0 );
      BisClosed = 1;
    } else if( mxGetNumberOfElements( prhs[1] ) == 2 ){
      mxVB = mxGetCell( prhs[1] , 0 );
      mxFB = mxGetCell( prhs[1] , 1 );
      if( mxGetNumberOfElements( mxFB ) == 0 ){ mxFB = NULL; BisClosed = 1; }
    } else {
      myErrMsgTxt("As a cell a pair vertices - faces was expected for the second argument.");
    }
  } else if( mxIsStruct( prhs[1] ) ){
    if( mxGetFieldNumber( prhs[1] ,"xyz" ) >= 0 ){ mxVB = mxGetField( prhs[1] ,0,"xyz" );
    } else if( mxGetFieldNumber( prhs[1] ,"vertices" ) >= 0 ){ mxVB = mxGetField( prhs[1] ,0,"vertices" );
    } else { myErrMsgTxt("As a struct a \".xyz\" field is compulsory for the second argument.");
    }
    if( mxGetFieldNumber( prhs[1] ,"tri" ) > 0 ){ mxFB = mxGetField( prhs[1] ,0,"tri" );
    } else if( mxGetFieldNumber( prhs[1] ,"faces" ) > 0 ){ mxFB = mxGetField( prhs[1] ,0,"faces" );
    }
  } else {
    mxVB = prhs[1];
  }
  if( mxFB != NULL && mySize( mxFB , 1 ) != 2 ){ mexErrMsgTxt("Segments for the second polygon were expected as an Sx2 matrix."); }

  
  
  NSD = mySize( mxVA , 1 );
  if( mySize( mxVB , 1 ) != NSD ){
    myErrMsgTxt("Number of Spatial Dimensions does not coincide.");
  }
  
  VA = myGetPr( mxVA ); nVA = mySize(  mxVA , 0 );
  if( mxFA == NULL  &&  AisClosed ){          nFA = nVA;
  } else if( mxFA == NULL  &&  !AisClosed ){  nFA = nVA - 1;
  } else {              FA = myGetPr( mxFA ); nFA = mySize(  mxFA , 0 );
  }
  
  
  VB = myGetPr( mxVB ); nVB = mySize(  mxVB , 0 );
  if( mxFB == NULL  &&  BisClosed ){          nFB = nVB;
  } else if( mxFB == NULL  &&  !BisClosed ){  nFB = nVB - 1;
  } else {              FB = myGetPr( mxFB ); nFB = mySize(  mxFB , 0 );
  }
  
  
  
    plhs[0] = mxCreateDoubleMatrix( nFA , nFB , mxREAL ); out_d = mxGetPr( plhs[0] );
  if( nlhs > 1 ){
    plhs[1] = mxCreateDoubleMatrix( nFA , nFB , mxREAL ); out_a = mxGetPr( plhs[1] );
  }
  if( nlhs > 2 ){
    plhs[2] = mxCreateDoubleMatrix( nFA , nFB , mxREAL ); out_b = mxGetPr( plhs[2] );
  }
  
  
  A0 = mxMalloc( NSD * sizeof( double ) );
  A1 = mxMalloc( NSD * sizeof( double ) );
  B0 = mxMalloc( NSD * sizeof( double ) );
  B1 = mxMalloc( NSD * sizeof( double ) );

  for( aa = 0 ; aa < nFA ; aa++ ){
    if( mxFA == NULL && !AisClosed ){        A0id = aa;                          A1id = aa+1;
    } else if( mxFA == NULL && AisClosed ){  A0id = aa;                          A1id = (aa+1) % nVA;
    } else {                                 A0id = ((int)FA[ aa       ])-1;     A1id = ((int)FA[ aa + nFA ])-1;
    }
    if( NSD == 2 ){         A0[0] = VA[ A0id ]; A0[1] = VA[ A0id + nVA ];
                            A1[0] = VA[ A1id ]; A1[1] = VA[ A1id + nVA ];
    } else if( NSD == 3 ){  A0[0] = VA[ A0id ]; A0[1] = VA[ A0id + nVA ]; A0[2] = VA[ A0id + 2*nVA ];
                            A1[0] = VA[ A1id ]; A1[1] = VA[ A1id + nVA ]; A1[2] = VA[ A1id + 2*nVA ];
    } else {                for( k = 0 ; k < NSD ; k++ ){ A0[k] = VA[ A0id + k*nVA ]; A1[k] = VA[ A1id + k*nVA ]; };
    }
    
    
    for( bb = 0 ; bb < nFB ; bb++ ){
      if( mxFB == NULL && !BisClosed ){        B0id = bb;                          B1id = bb+1;
      } else if( mxFB == NULL && BisClosed ){  B0id = bb;                          B1id = (bb+1) % nVB;
      } else {                                 B0id = ((int)FB[ bb       ])-1;     B1id = ((int)FB[ bb + nFB ])-1;
      }
      if( NSD == 2 ){         B0[0] = VB[ B0id ]; B0[1] = VB[ B0id + nVB ];
                              B1[0] = VB[ B1id ]; B1[1] = VB[ B1id + nVB ];
      } else if( NSD == 3 ){  B0[0] = VB[ B0id ]; B0[1] = VB[ B0id + nVB ]; B0[2] = VB[ B0id + 2*nVB ];
                              B1[0] = VB[ B1id ]; B1[1] = VB[ B1id + nVB ]; B1[2] = VB[ B1id + 2*nVB ];
      } else {                for( k = 0 ; k < NSD ; k++ ){ B0[k] = VB[ B0id + k*nVB ]; B1[k] = VB[ B1id + k*nVB ]; };
      }
      
      param_Ray2Ray( A0 , A1 , B0 , B1 , NSD , a , b );
      if( a[0] < 0 ){ a[0] = 0; }
      if( a[0] > 1 ){ a[0] = 1; }
      if( b[0] < 0 ){ b[0] = 0; }
      if( b[0] > 1 ){ b[0] = 1; }
      d2 = d2_RayParam2RayParam( A0 , A1 , a[0] , B0 , B1 , b[0] , NSD );
      
                      out_d[ aa + bb*nFA ] = sqrt(d2);
      if( nlhs > 1 ){ out_a[ aa + bb*nFA ] = a[0];       }
      if( nlhs > 2 ){ out_b[ aa + bb*nFA ] = b[0];       }
    }
  }
  
  EXIT:
    if( A0 != NULL ){ mxFree( A0 ); }
    if( A1 != NULL ){ mxFree( A1 ); }
    if( B0 != NULL ){ mxFree( B0 ); }
    if( B1 != NULL ){ mxFree( B1 ); }
    myFreeALLOCATES();
}

