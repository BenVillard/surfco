function [X,E,H] = Optimize1D_exhaustive( f , X , P )

  if nargin < 3, P = struct; end
  P = joinfields( P , 'MAX_ITS'         ,  Inf   );      %maximum number of iterations
  P = joinfields( P , 'EPS'             ,  0     );      %min value


  P = joinfields( P , 'N'               ,  10    );      %number of partitions
  P = joinfields( P , 'ALLOWLOG'        ,  true  );      %number of partitions
  
  P = joinfields( P , 'MIN'             , -Inf   );      %minimo valor de X1
  P = joinfields( P , 'MAX'             ,  Inf   );      %maximo valor de X4
  
  P = joinfields( P , 'VERBOSE'         , []     );       %espera una funcion que recibe:  verbose_fcn( level , str , varargin )
  P = joinfields( P , 'PLOT'            , false  );       %espera un handle donde plotear

  X1 = X{1}{1}; E1 = X{1}{2};
  X2 = X{2}{1}; E2 = X{2}{2};
  
  if isnan(X1) || isnan(X2), error('valor de X incorrecto'); end
  if X1 > X2, error('X1 es mayor que X2!!'); end

  if isnan(E1), E1 = f( X1 ); end
  if isnan(E2), E2 = f( X2 ); end

  
  H.PARAMETERS = P;
  H.initialX   = { {X1,E1},{X2,E2} };
  H.ITEX = 0;


  try, feval( P.VERBOSE , 3 , '    EXHAUSTIVE search:  { { %g , %g } , {%g , %g } }', X1 , E1 , X2 , E2 ); end
  
  
  
  
  E0 = min( E1 , E2 );
  
  while true
    H.ITEX = H.ITEX+1;
    

    DX = X2-X1;
    if P.ALLOWLOG && ( ( X1 > 0 && X2 > 0 ) || ( X1 < 0 && X2 < 0 ) )
      Xs = unique( [ geospace( X1 , X2 , P.N ) , linspace( X1 , X2 , P.N ) ] );
    elseif P.ALLOWLOG && X1 == 0 && X2 > 0
      try
        Xs = unique( [ geospace( X2/(2^20) , X2 , P.N ) , linspace( X1 , X2 , P.N ) ] );
      catch
        1
      end
    elseif P.ALLOWLOG && X1 < 0 && X2 == 0
      Xs = unique( [ geospace( X1 , X1/(2^20) , P.N ) , linspace( X1 , X2 , P.N ) ] );
    else
      Xs = unique( linspace( X1 , X2 , P.N ) );
    end
    
    Es = zeros( numel(Xs) , 1 );
    for i = 1:numel(Xs), Es(i) = f( Xs(i) ); end
    
    if P.PLOT
      figure; plot(Xs,Es,'.-');
    end
    
    [E,id] = min( Es ); X = Xs(id);
    
    if E < E0
      E0 = E;
    end


    try, feval( P.VERBOSE , 3 , '\{    EXHAUSTIVE search: }r%4d. ( %.20g : (%d) : %.20g ) -> ( %g , %g )', H.ITEX , X1 , numel(Xs) , X2 ,  E , X ); end

    if P.PLOT
      line('Parent',P.PLOT,'XData',Xs,'YData',Es,'Marker','o','color',[0 0.8 0]);
    end
    
    
    if H.ITEX >= P.MAX_ITS                   ,  H.STOP = sprintf('H.ITEX(%d) > P.MAX_ITS(%d)' ,H.ITEX,P.MAX_ITS);                  break; end

    if abs(DX) < P.EPS                      ,  H.STOP = sprintf('X1(%g) and X2(%g) in P.EPS(%d)' ,X1,X2,P.EPS);                   break; end
    if numel(Xs) < 4                        ,  H.STOP = sprintf('converged' );                                                    break; end


    if      id == 1    &&  X1-DX > P.MIN
      X1 = X1 - DX;
      X2 = X1 + DX;
    elseif  id == 1
      X1 = P.MIN;
      X2 = Xs(2);
    elseif  id == numel(Xs)  &&  X2+DX < P.MAX
      X1 = X2 - DX;
      X2 = X2 + DX;
    elseif  id == numel(Xs)
      X1 = Xs(end-1);
      X2 = P.MAX;
    else
      X1 = Xs(id-1);
      X2 = Xs(id+1);
    end
%     disp([X1 X2])

    if X2 <= X1                             ,  H.STOP = sprintf('X2(%g) <= X1(%g)'                 ,X2,X1);                       break; end
    
  end

  try, feval( P.VERBOSE , 3 , '\r    EXHAUSTIVE search finish:  %s\n', H.STOP ); end

end
