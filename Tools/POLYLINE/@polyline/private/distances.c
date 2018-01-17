double Solve2x2( double b1 , double b2 , double m11 , double m21 , double m12 , double m22 , double *x1 , double *x2 );
double d2_Point2Ray( double *P , double *B0 , double *B1 , int NSD );
double d2_Point2RayParam( double *P , double *B0 , double *B1 , double b , int NSD );
double d2_RayParam2RayParam( double *A0 , double *A1 , double a , double *B0 , double *B1 , double b , int NSD );
double d2_Point2Point( double *P , double *Q , int NSD );
void   RayParam2Point( double *A0 , double *A1 , int NSD , double a , double *P );
void   param_Ray2Ray( double *A0 , double *A1 , double *B0 , double *B1 , int NSD , double *a , double *b );
void   param_Point2Ray( double *P , double *B0 , double *B1 , int NSD , double *a );



double Solve2x2( double b1 , double b2 , double m11 , double m21 , double m12 , double m22 , double *x1 , double *x2 ){
  double det;
  
  det = m11*m22 - m12*m21;
  
  *x1 =  -(b2*m12 - b1*m22)/det;
  *x2 =   (b2*m11 - b1*m21)/det;
  
  return( det );
}




void RayParam2Point( double *A0 , double *A1 , int NSD , double a , double *P ){
  int d;
  
  if( NSD == 2 ){
    P[0] = A0[0] + a*( A1[0] - A0[0] );
    P[1] = A0[1] + a*( A1[1] - A0[1] );
  } else if( NSD == 3 ){
    P[0] = A0[0] + a*( A1[0] - A0[0] );
    P[1] = A0[1] + a*( A1[1] - A0[1] );
    P[2] = A0[2] + a*( A1[2] - A0[2] );
  } else {
    for( d = 0 ; d < NSD ; d++ ){
      P[d] = A0[d] + a*( A1[d] - A0[d] );
    }
  }
  
}


double d2_Point2Ray( double *P , double *B0 , double *B1 , int NSD ){
  double D2, b[1];
  param_Point2Ray( P , B0 , B1 , NSD , b );
  D2 = d2_Point2RayParam( P , B0 , B1 , b[0] , NSD );
  return( D2 );
}

double d2_Point2RayParam( double *P , double *B0 , double *B1 , double b , int NSD ){

  int    d;
  double D2, p, q, r;
  
  if( NSD == 2 ){
    p = B0[0] + b*( B1[0] - B0[0] ) - P[0];
    q = B0[1] + b*( B1[1] - B0[1] ) - P[1];
    D2 = p*p + q*q;
  } else if( NSD == 3 ){
    p = B0[0] + b*( B1[0] - B0[0] ) - P[0];
    q = B0[1] + b*( B1[1] - B0[1] ) - P[1];
    q = B0[2] + b*( B1[2] - B0[2] ) - P[2];
    D2 = p*p + q*q + r*r;
  } else {
    D2 = 0;
    for( d = 0 ; d < NSD ; d++ ){
      p = P[d];
      q = B0[d] + b*( B1[d] - B0[d] );
      r = p-q;
      D2 += r*r;
    }
  }
  return( D2 );

}

double d2_RayParam2RayParam( double *A0 , double *A1 , double a , double *B0 , double *B1 , double b , int NSD ){
  int    d;
  double D2, p, q, r;
  
  if( NSD == 2 ){
    p = A0[0] + a*( A1[0] - A0[0] ) - ( B0[0] + b*( B1[0] - B0[0] ) );
    q = A0[1] + a*( A1[1] - A0[1] ) - ( B0[1] + b*( B1[1] - B0[1] ) );
    D2 = p*p + q*q;
  } else if( NSD == 3 ){
    p = A0[0] + a*( A1[0] - A0[0] ) - ( B0[0] + b*( B1[0] - B0[0] ) );
    q = A0[1] + a*( A1[1] - A0[1] ) - ( B0[1] + b*( B1[1] - B0[1] ) );
    r = A0[2] + a*( A1[2] - A0[2] ) - ( B0[2] + b*( B1[2] - B0[2] ) );
    D2 = p*p + q*q + r*r;
  } else {
    D2 = 0;
    for( d = 0 ; d < NSD ; d++ ){
      p = A0[d] + a*( A1[d] - A0[d] );
      q = B0[d] + b*( B1[d] - B0[d] );
      r = p-q;
      D2 += r*r;
    }
  }
  return( D2 );
}




double d2_Point2Point( double *A , double *B , int NSD ){
  int    d;
  double D2, p, q, r;
  
  if( NSD == 2 ){
    p = A[0] - B[0];
    q = A[1] - B[1];
    D2 = p*p + q*q;
  } else if( NSD == 3 ){
    p = A[0] - B[0];
    q = A[1] - B[1];
    r = A[2] - B[2];
    D2 = p*p + q*q + r*r;
  } else {
    D2 = 0;
    for( d = 0 ; d < NSD ; d++ ){
      p = A[d]; q = B[d]; r = p-q;
      D2 += r*r;
    }
  }
  return( D2 );
}

#define  SINGULAR_DET    1e-8
void param_Ray2Ray( double *A0 , double *A1 , double *B0 , double *B1 , int NSD , double *a , double *b ){
  int    d;
  double A0xA0, A0xA1, A0xB0, A0xB1, A1xA1, A1xB0, A1xB1, B0xB0, B0xB1, B1xB1;
  double r, det;
  double D2, thisD2;
  double w[1];

  if( NSD == 3 ){
    #define dot3(x,y) ( x[0]*y[0] + x[1]*y[1] + x[2]*y[2] )
    
    A0xA0 = dot3( A0 , A0 );
    A0xA1 = dot3( A0 , A1 );
    A0xB0 = dot3( A0 , B0 );
    A0xB1 = dot3( A0 , B1 );
    
    A1xA1 = dot3( A1 , A1 );
    A1xB0 = dot3( A1 , B0 );
    A1xB1 = dot3( A1 , B1 );
    
    B0xB0 = dot3( B0 , B0 );
    B0xB1 = dot3( B0 , B1 );
    
    B1xB1 = dot3( B1 , B1 );

  } else if( NSD == 2 ){
    #define dot2(x,y) ( x[0]*y[0] + x[1]*y[1] )
    
    A0xA0 = dot2( A0 , A0 );
    A0xA1 = dot2( A0 , A1 );
    A0xB0 = dot2( A0 , B0 );
    A0xB1 = dot2( A0 , B1 );
    
    A1xA1 = dot2( A1 , A1 );
    A1xB0 = dot2( A1 , B0 );
    A1xB1 = dot2( A1 , B1 );
    
    B0xB0 = dot2( B0 , B0 );
    B0xB1 = dot2( B0 , B1 );
    
    B1xB1 = dot2( B1 , B1 );
    
  } else {
    
    A0xA0 = 0; A0xA1 = 0; A0xB0 = 0; A0xB1 = 0; A1xA1 = 0; A1xB0 = 0; A1xB1 = 0; B0xB0 = 0; B0xB1 = 0; B1xB1 = 0;
    for( d = 0 ; d < NSD ; d++ ){
      A0xA0 += A0[d]*A0[d];
      A0xA1 += A0[d]*A1[d];
      A0xB0 += A0[d]*B0[d];
      A0xB1 += A0[d]*B1[d];

      A1xA1 += A1[d]*A1[d];
      A1xB0 += A1[d]*B0[d];
      A1xB1 += A1[d]*B1[d];

      B0xB0 += B0[d]*B0[d];
      B0xB1 += B0[d]*B1[d];

      B1xB1 += B1[d]*B1[d];
    }
  }
  
  r = - A0xB0 + A1xB0 + A0xB1 - A1xB1;
  det = Solve2x2( + A0xA0 - A0xA1 - A0xB0 + A1xB0 ,
                  - A0xB0 + B0xB0 + A0xB1 - B0xB1 ,
                  + A0xA0 - 2*A0xA1 + A1xA1 , r                ,
                   r               , B0xB0 - 2*B0xB1 + B1xB1   ,
                   a , b );
  if( fabs( det ) < SINGULAR_DET ){  // they look parallel
    
    D2 = 1e308;
    
    //try A0-B0B1 distance
    param_Point2Ray( A0 , B0 , B1 , NSD , w );
    if( w[0] < 0 ){ w[0] = 0; } else if( w[0] > 1 ){ w[0] = 1; }
    thisD2 = d2_Point2RayParam( A0 , B0 , B1 , w[0] , NSD );
//     mexPrintf("thisD2: %g\n",thisD2);
    if( 1 || thisD2 < D2 ){
      D2 = thisD2; *a = 0; *b = w[0];
    }
    
    //try A1-B0B1 distance
    param_Point2Ray( A1 , B0 , B1 , NSD , w );
    if( w[0] < 0 ){ w[0] = 0; } else if( w[0] > 1 ){ w[0] = 1; }
    thisD2 = d2_Point2RayParam( A1 , B0 , B1 , w[0] , NSD );
//     mexPrintf("thisD2: %g\n",thisD2);
    if( thisD2 < D2 ){
      D2 = thisD2; *a = 1; *b = w[0];
    }

    //try B0-A0A1 distance
    param_Point2Ray( B0 , A0 , A1 , NSD , w );
    if( w[0] < 0 ){ w[0] = 0; } else if( w[0] > 1 ){ w[0] = 1; }
    thisD2 = d2_Point2RayParam( B0 , A0 , A1 , w[0] , NSD );
//     mexPrintf("thisD2: %g\n",thisD2);
    if( thisD2 < D2 ){
      D2 = thisD2; *a = w[0]; *b = 0;
    }
    
    //try B1-A0A1 distance
    param_Point2Ray( B1 , A0 , A1 , NSD , w );
    if( w[0] < 0 ){ w[0] = 0; } else if( w[0] > 1 ){ w[0] = 1; }
    thisD2 = d2_Point2RayParam( B1 , A0 , A1 , w[0] , NSD );
//     mexPrintf("thisD2: %g\n",thisD2);
    if( thisD2 < D2 ){
      D2 = thisD2; *a = w[0]; *b = 1;
    }
  }
}


void param_Point2Ray( double *P , double *B0 , double *B1 , int NSD , double *b ){

  int    d;
  double B0xB0, B0xB1, B0xP, B1xB1, B1xP;
  
  if( NSD == 3 ){
    #define dot3(x,y) ( x[0]*y[0] + x[1]*y[1] + x[2]*y[2] )
    
    B0xB0 = dot3( B0 , B0 );
    B0xB1 = dot3( B0 , B1 );
    B0xP  = dot3( B0 , P  );
    B1xB1 = dot3( B1 , B1 );
    B1xP  = dot3( B1 , P  );

  } else if( NSD == 2 ){
    #define dot2(x,y) ( x[0]*y[0] + x[1]*y[1] )
    
    B0xB0 = dot2( B0 , B0 );
    B0xB1 = dot2( B0 , B1 );
    B0xP  = dot2( B0 , P  );
    B1xB1 = dot2( B1 , B1 );
    B1xP  = dot2( B1 , P  );
    
  } else {
    
    B0xB0 = 0; B0xB1 = 0; B0xP = 0; B1xB1 = 0; B1xP = 0;
    for( d = 0 ; d < NSD ; d++ ){
      B0xB0 += B0[d]*B0[d];
      B0xB1 += B0[d]*B1[d];
      B0xP  += B0[d]*P[d];
      B1xB1 += B1[d]*B1[d];
      B1xP  += B1[d]*P[d];
    }
  }  
  
  *b = -( B0xB1 - B0xB0 - B1xP + B0xP ) / ( B1xB1 - 2*B0xB1 + B0xB0 );
}

