function [X,E,H] = Optimize1D_smart( f , P )


  if nargin < 2, P = struct; end
  P = joinfields( P , 'MAX_ITS'         ,  Inf   );       %maximum number of iterations
  P = joinfields( P , 'EPS'             ,  0     );
  
  P = joinfields( P , 'MIN'             , -Inf   );       %minimo valor de X1
  P = joinfields( P , 'MAX'             ,  Inf   );       %maximo valor de X4

  P = joinfields( P , 'VERBOSE'         , []     );       %espera una funcion que recibe:  verbose_fcn( level , str , varargin )
  P = joinfields( P , 'PLOT'            , false  );       %espera un handle donde plotear


  
  BETTERS = P.getBETTERS();
  if      numel( BETTERS ) == 0
    error('no puedo');
  elseif  numel( BETTERS ) == 1
    error('no puedo');
  elseif  numel( BETTERS ) == 2
    f( ( BETTERS{1}{1} + BETTERS{2}{1} )/2 );
    BETTERS = P.getBETTERS();
  end
  
  

  H.PARAMETERS = P;
  H.ITSM = 0;
  

  try, feval( P.VERBOSE , 3 , '    SMART search: '); end

  xs_prev = [];

  while 1
    H.ITSM = H.ITSM + 1;


    xs = cellfun( @(b) b{1} , BETTERS );
    if isequal( xs , xs_prev )      ,  H.STOP = sprintf('stuck'  );                                             break; end
    

    xs_prev = xs;
    
    es = cellfun( @(b) b{2} , BETTERS );
    E  = min( es );

    W = diag( diff( dualVector( xs ) ).^2 );
%     W = diag( diff( dualVector( xs ) ).^1 );
    
    USE_Q = false;
    if numel( xs ) == 3
      USE_Q = true;
    end

    if ~USE_Q
      
      X = [ xs(:).^3 , xs(:).^2 , xs(:) , ones(numel(xs),1) ];

      setWarning('off','MATLAB:rankDeficientMatrix');setWarning('off','MATLAB:nearlySingularMatrix');setWarning('off','MATLAB:singularMatrix');
      ABCD = es(:)'*W*X / ( X'*W*X );
      restoreWarning(  'MATLAB:rankDeficientMatrix');restoreWarning(  'MATLAB:nearlySingularMatrix');restoreWarning(  'MATLAB:singularMatrix');
      
      discriminant = ABCD(2)^2 - 3*ABCD(1)*ABCD(3);

      if abs( ABCD(1) ) < 1e-8
        USE_Q = true;
      elseif discriminant < 0
        USE_Q = true;
      elseif discriminant == 0
        X = -ABCD(2) / ( 3*ABCD(1) );
      else
        discriminant = sqrt( discriminant );
        r1 = ( - ABCD(2) - discriminant )/( 3*ABCD(1) );
        r2 = ( - ABCD(2) + discriminant )/( 3*ABCD(1) );

        if ( r1.^(3:-1:0) )*ABCD(:)  <  ( r2.^(3:-1:0) )*ABCD(:)
          X = r1;
        else
          X = r2;
        end
        
      end
    end
    
    
    if USE_Q
      X = [ xs(:).^2 , xs(:) , ones(numel(xs),1) ];

      setWarning('off','MATLAB:rankDeficientMatrix');setWarning('off','MATLAB:nearlySingularMatrix');setWarning('off','MATLAB:singularMatrix');
  %     ABC = linsolve( X.'*W*X , X'*W*es(:) , struct('SYM',true) );
      ABC = es(:)'*W*X / ( X'*W*X );
      restoreWarning(  'MATLAB:rankDeficientMatrix');restoreWarning(  'MATLAB:nearlySingularMatrix');restoreWarning(  'MATLAB:singularMatrix');

      if ABC(1) <  0                  ,  H.STOP = sprintf('concave' );                                            break; end
      if ABC(1) == 0                  ,  H.STOP = sprintf('line' );                                               break; end

      X = - ABC(2)/( 2 * ABC(1) );
    end
    
    if isnan(X)                     ,  H.STOP = sprintf('X NaN' );            break; end
    if isinf(X)                     ,  H.STOP = sprintf('X Inf' );            break; end
    if X > P.MAX                    ,  H.STOP = sprintf('X > P.MAX' );        break; end
    if X < P.MIN                    ,  H.STOP = sprintf('X < P.MIN' );        break; end

    E = f( X );
    
    if P.PLOT
      line( 'Parent', P.PLOT ,'XData', xs ,'YData', es ,'marker','o','markerfacecolor',[0.5 0.5 0.5],'color',[0.5 0.5 0.5] , 'linestyle','none','markersize',3 );
      
      xx = linspace( min(xs) , max(xs) , 150 );
      
      if USE_Q
        ee = ABC(:).' * [ xx.^2 ; xx ; ones(1,150) ];
        line( 'Parent', P.PLOT ,'XData', xx ,'YData', ee ,'marker','none','color',[0.5 0.5 0.5] , 'linestyle',':' );
      else
        ee = ABCD(:).' * [ xx.^3 ; xx.^2 ; xx ; ones(1,150) ];
        line( 'Parent', P.PLOT ,'XData', xx ,'YData', ee ,'marker','none','color',[0.5 0.1 0.5] , 'linestyle',':' );
      end
    end
    
    
    if USE_Q
      try, feval( P.VERBOSE , 3 , '\{    SMART search: }r  QUADRATIC %.25g   ->  %.25g' , min(es) , E ); end
    else
      try, feval( P.VERBOSE , 3 , '\{    SMART search: }r  CUBIC     %.25g   ->  %.25g' , min(es) , E ); end
    end

    
    if H.ITSM >= P.MAX_ITS          ,  H.STOP = sprintf('H.ITSM(%d) > P.MAX_ITS(%d)' ,H.ITSM,P.MAX_ITS);        break; end
    
    BETTERS = P.getBETTERS();
  end


  try, feval( P.VERBOSE , 3 , '\r    SMART search finish:  %s\n', H.STOP ); end
  
  

end
