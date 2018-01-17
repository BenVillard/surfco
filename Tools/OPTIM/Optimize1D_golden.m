function [X,E,H] = Optimize1D_golden( f , X , P )

  if nargin < 3, P = struct; end
  P = joinfields( P , 'MAX_ITS'         ,  Inf   );       %maximum number of iterations
  P = joinfields( P , 'EPS'             ,  0     );
  
  P = joinfields( P , 'MIN'             , -Inf   );       %minimo valor de X1
  P = joinfields( P , 'MAX'             ,  Inf   );       %maximo valor de X4

  P = joinfields( P , 'VERBOSE'         , []     );       %espera una funcion que recibe:  verbose_fcn( level , str , varargin )
  P = joinfields( P , 'PLOT'            , false  );       %espera un handle donde plotear

  
  X1 = X{1}{1}; E1 = X{1}{2};
  X4 = X{2}{1}; E4 = X{2}{2};

  if isnan(X1) || isnan(X4), error('valor de X incorrecto'); end
  if X1 > X4, error('X1 es mayor que X4!!'); end

  if isnan(E1), E1 = f( X1 ); end
  if isnan(E4), E4 = f( X4 ); end

  
  H.PARAMETERS = P;
  H.initialX   = { {X1,E1},{X4,E4} };
  H.ITGS = 0;

  
  try, feval( P.VERBOSE , 3 , '    GOLDEN search:  { { %g , %g } , {%g , %g } }', X1 , E1 , X4 , E4 ); end
  
  
  GR = (3 - sqrt(5))/2;
  DX = X4-X1;
  X2 = X1 + DX*GR; E2 = f( X2 );
  X3 = X4 - DX*GR; E3 = f( X3 );

  

  E0 = min( [ E1  E2  E3  E4 ] );
  
  PREV_CASE = 0;
  while true

    H.ITGS = H.ITGS+1;


    DX = X4-X1;
    
    [E,id] = min([E1 E2 E3 E4]);
    X = [X1 X2 X3 X4]; X = X(id);

    if E < E0
      E0 = E;
    end
    
    
    if     id == 1  &&  ( X4 - DX/(1-GR) ) >= P.MIN
      idd = -1;
    elseif id == 4  &&  ( X1 + DX/(1-GR) ) <= P.MAX
      idd = -4;
    else
      idd = id;
    end
    
    try, feval( P.VERBOSE , 3 , '\{    GOLDEN search: }r%4d. %.15g + [ %g %g %g %g ] -> ( %g %g %g %g )  [%d]', H.ITGS,X1,0,X2-X1,X3-X1,X4-X1,E1,E2,E3,E4,idd ); end
    
    if P.PLOT
      line('Parent',P.PLOT,'XData',[X1 X2 X3 X4],'YData',[E1 E2 E3 E4],'Marker','o','color',[0 0 1]);
    end

    
    if H.ITGS > P.MAX_ITS           ,  H.STOP = sprintf('H.ITGS(%d) > P.MAX_ITS(%d)' ,H.ITGS,P.MAX_ITS);        break; end
    if abs(DX) < P.EPS              ,  H.STOP = sprintf('X1(%g) and X4(%g) in P.EPS(%d)',X1,X4,P.EPS);          break; end
    
    %if X1     < P.MIN               ,  H.STOP = sprintf('X1(%g) < P.LS_GOLDEN_MIN(%g)'         ,X1,P.MIN);      break; end
    %if X4     > P.MAX               ,  H.STOP = sprintf('X4(%g) > P.LS_GOLDEN_MAX(%g)'         ,X4,P.MAX);      break; end
    
    if X2 <= X1                     ,  H.STOP = sprintf('X2(%g) <= X1(%g)'                     ,X2,X1);         break; end
    if X3 <= X2                     ,  H.STOP = sprintf('X3(%g) <= X2(%g)'                     ,X3,X2);         break; end
    if X4 <= X3                     ,  H.STOP = sprintf('X4(%g) <= X3(%g)'                     ,X4,X3);         break; end
    
    switch idd                                                             % ( b -> the best ) ( n -> the new point ) ( [] the point to remove )
      case -1                                                              %   n           b       |  [|]      |
        X3 = X2; E3 = E2;
        X2 = X1; E2 = E1;
        X1 = X4 - (X4-X1)/(1-GR);   E1 = f( X1 );
      case  1                                                              %               b   n   |   |      [|]
        X4 = X3; E4 = E3;
        X3 = X2; E3 = E2;
        X2 = X1 + (X4-X1)*GR;       E2 = f( X2 );
      case  2                                                              %               |   n   b   |      [|]
        X4 = X3; E4 = E3;
        X3 = X2; E3 = E2;
        X2 = X1 + (X4-X1)*GR;       E2 = f( X2 );
      case  3                                                              %              [|]      |   b   n   |
        X1 = X2; E1 = E2;
        X2 = X3; E2 = E3;
        X3 = X4 - (X4-X1)*GR;       E3 = f( X3 );
      case  4                                                              %              [|]      |   |   n   b
        X1 = X2; E1 = E2;
        X2 = X3; E2 = E3;
        X3 = X4 - (X4-X1)*GR;       E3 = f( X3 );
      case -4                                                              %               |      [|]  |       b           n
        X2 = X3; E2 = E3;
        X3 = X4; E3 = E4;
        X4 = X1 + (X4-X1)/(1-GR);   E4 = f( X4 );
    end
    
    PREV_CASE = idd;
  end

  try, feval( P.VERBOSE , 3 , '\r    GOLDEN search finish:  %s\n', H.STOP ); end
  
end
