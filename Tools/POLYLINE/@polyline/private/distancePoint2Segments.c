/*

V = [0 0;1 0;2 2;0 1;0 0];
V = Interp1D( V , 1:size(V,1) , linspace(1,size(V,1),100) );
%V = V + randn(size(V))/50;
F = [ 1:size(V,1)-1 ; 2:size(V,1) ].';
P = rand(2000,2)*3;


close all;
 plot3d( V(F(:,1),:) , V(F(:,2),:) , '.-r' ); axis equal
hplot3d( P , 'ob' )
[s,cp,d,c] = distancePoint2Segments( P , V );
hplot3d( P , cp , '.-k' );

*/
/*
//mmex INCLUDES = { fileparts( which( 'vec.m' ) ) };
*/

#include "myMEX.h"
#include "distances.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  int    NSD, k;

  double *A0 = NULL;
  double *A1 = NULL;
  double *best_A0 = NULL;
  double *best_A1 = NULL;
  double *XYZ = NULL;
  
  double a[1], d2, best_a, this_d2;
  
  int    s, nVA, nFA, nP, pp, best_s;
  double *VA, *FA, *P;
  double *out_s, *out_cp, *out_d, *out_c;
  
  mxArray  *mxVA, *mxFA;
  int AisClosed = 0;
  int A0id, A1id;
  
  

  if( nrhs != 2 ){
    myErrMsgTxt("2 input arguments were expected.");
  }

  mxFA = NULL; AisClosed = 0;
  if( mxIsCell( prhs[1] ) ){
    if( mxGetNumberOfElements( prhs[1] ) == 1 ){
      mxVA = mxGetCell( prhs[1] , 0 );
      AisClosed = 1;
    } else if( mxGetNumberOfElements( prhs[1] ) == 2 ){
      mxVA = mxGetCell( prhs[1] , 0 );
      mxFA = mxGetCell( prhs[1] , 1 );
      if( mxGetNumberOfElements( mxFA ) == 0 ){ mxFA = NULL; AisClosed = 1; }
    } else {
      myErrMsgTxt("As a cell a pair vertices - faces was expected for the first argument.");
    }
  } else if( mxIsStruct( prhs[1] ) ){
    if( mxGetFieldNumber( prhs[1] ,"xyz" ) >= 0 ){ mxVA = mxGetField( prhs[1] ,0,"xyz" );
    } else if( mxGetFieldNumber( prhs[1] ,"vertices" ) >= 0 ){ mxVA = mxGetField( prhs[1] ,0,"vertices" );
    } else { myErrMsgTxt("As a struct a \".xyz\" field is compulsory for the first argument.");
    }
    if( mxGetFieldNumber( prhs[1] ,"tri" ) > 0 ){ mxFA = mxGetField( prhs[1] ,0,"tri" );
    } else if( mxGetFieldNumber( prhs[1] ,"faces" ) > 0 ){ mxFA = mxGetField( prhs[1] ,0,"faces" );
    }
  } else {
    mxVA = prhs[1];
  }
  if( mxFA != NULL && mySize( mxFA , 1 ) != 2 ){ mexErrMsgTxt("Segments for the first polygon were expected as an Sx2 matrix."); }


  NSD = mySize( prhs[0] , 1 );
  if( mySize( mxVA , 1 ) != NSD ){
    myErrMsgTxt("Number of Spatial Dimensions does not coincide.");
  }
  
  
  P  = myGetPr( prhs[0] ); nP = mySize(  prhs[0] , 0 );
  
  VA = myGetPr( mxVA ); nVA = mySize(  mxVA , 0 );
  if( mxFA == NULL  &&  AisClosed ){          nFA = nVA;
  } else if( mxFA == NULL  &&  !AisClosed ){  nFA = nVA - 1;
  } else {              FA = myGetPr( mxFA ); nFA = mySize(  mxFA , 0 );
  }
  

    plhs[0] = mxCreateDoubleMatrix( nP , 1 , mxREAL ); out_s = mxGetPr( plhs[0] );  //segment id
  if( nlhs > 1 ){
    plhs[1] = mxCreateDoubleMatrix( nP , NSD , mxREAL ); out_cp = mxGetPr( plhs[1] ); //closed point (on the segment id)
  }
  if( nlhs > 2 ){
    plhs[2] = mxCreateDoubleMatrix( nP , 1 , mxREAL ); out_d = mxGetPr( plhs[2] );  //distance
  }
  if( nlhs > 3 ){
    plhs[3] = mxCreateDoubleMatrix( nP , 1 , mxREAL ); out_c = mxGetPr( plhs[3] );  //parametric coor within   segment
  }
  
  
  XYZ = mxMalloc( NSD * sizeof( double ) );
  A0  = mxMalloc( NSD * sizeof( double ) );
  A1  = mxMalloc( NSD * sizeof( double ) );
  best_A0  = mxMalloc( NSD * sizeof( double ) );
  best_A1  = mxMalloc( NSD * sizeof( double ) );

  
  for( pp = 0 ; pp < nP ; pp++ ){
    if( NSD == 2 ){         XYZ[0] = P[ pp ]; XYZ[1] = P[ pp + nP ];
    } else if( NSD == 3 ){  XYZ[0] = P[ pp ]; XYZ[1] = P[ pp + nP ]; XYZ[2] = P[ pp + 2*nP ];
    } else {                for( k = 0 ; k < NSD ; k++ ){ XYZ[k] = P[ pp + k*nP ]; };
    }
    
    d2 = 1e308; best_s = -1;
    for( s = 0 ; s < nFA ; s++ ){
      if( mxFA == NULL && !AisClosed ){        A0id = s;                          A1id = s+1;
      } else if( mxFA == NULL && AisClosed ){  A0id = s;                          A1id = (s+1) % nVA;
      } else {                                 A0id = ((int)FA[ s       ])-1;     A1id = ((int)FA[ s + nFA ])-1;
      }
      if( NSD == 2 ){         A0[0] = VA[ A0id ]; A0[1] = VA[ A0id + nVA ];
                              A1[0] = VA[ A1id ]; A1[1] = VA[ A1id + nVA ];
      } else if( NSD == 3 ){  A0[0] = VA[ A0id ]; A0[1] = VA[ A0id + nVA ]; A0[2] = VA[ A0id + 2*nVA ];
                              A1[0] = VA[ A1id ]; A1[1] = VA[ A1id + nVA ]; A1[2] = VA[ A1id + 2*nVA ];
      } else {                for( k = 0 ; k < NSD ; k++ ){ A0[k] = VA[ A0id + k*nVA ]; A1[k] = VA[ A1id + k*nVA ]; };
      }

      param_Point2Ray( XYZ , A0 , A1 , NSD , a );
      if( a[0] < 0 ){        a[0] = 0;
      } else if( a[0] > 1 ){ a[0] = 1;
      }
      
      this_d2 = d2_Point2RayParam( XYZ , A0 , A1 , a[0] , NSD );
      
      if( this_d2 < d2 ){
        d2 = this_d2;
        best_s  = s;
        best_a  = a[0];
        memcpy( best_A0 , A0 , NSD * sizeof( double ) );
        memcpy( best_A1 , A1 , NSD * sizeof( double ) );
      }
      
    }
    
    out_s[ pp ] = best_s + 1;
    if( nlhs > 1 ){
      RayParam2Point( best_A0 , best_A1 , NSD , best_a , XYZ );
      for( k = 0 ; k < NSD ; k++ ){ out_cp[ pp + k*nP ] = XYZ[k]; }
    }
    if( nlhs > 2 ){
      out_d[ pp ] = sqrt(d2);
    }
    if( nlhs > 3 ){
      out_c[ pp ] = best_a;
    }

  }
  
  EXIT:
    if( XYZ != NULL ){ mxFree( XYZ ); }
    if( A0  != NULL ){ mxFree( A0  ); }
    if( A1  != NULL ){ mxFree( A1  ); }
    if( best_A0  != NULL ){ mxFree( best_A0  ); }
    if( best_A1  != NULL ){ mxFree( best_A1  ); }
    myFreeALLOCATES();
}

