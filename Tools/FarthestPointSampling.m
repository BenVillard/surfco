function [P,Pids] = FarthestPointSampling( P , Pids , minD , maxN , VERBOSE )

  if nargin < 5, VERBOSE = false; end

  if isequal( P(1,:) , P(end,:) ), P(end,:) = []; end
  if nargin < 2 || isempty( Pids )
%     D = ipd( P , [] );
%     D(1:size(D,1)+1:end) = NaN;
%     [~,Pids] = min( min( D , [] , 2 ) );
    Pids = 1;
  end
  
  try
    dPoint2Point([0 0],[1 1]);
    IPD = @(x,y)dPoint2Point(x,y);
  catch
    IPD = @(x,y)ipd(x,y);
  end
    
  if nargin < 3 || isempty( minD ), minD = -Inf; end
  if nargin < 4 || isempty( maxN ), maxN =  Inf; end
  maxN = min( maxN , size(P,1) );
  
  %P( Pids , : ) = NaN;
  D = IPD( P , P( Pids , : ) );
  while 1
    [ d , new ] = max( min( D , [] , 2 ) );
    if d < minD, break; end

    Pids(end+1) = new;

    if numel( Pids ) >= maxN, break; end
    
    D = [ D , IPD( P , P( new , : ) ) ];

    if VERBOSE
      fprintf('%d   %g\n' , numel( Pids ) , d );
    end
  
  end
  P = P( sort(Pids) , : );

end
