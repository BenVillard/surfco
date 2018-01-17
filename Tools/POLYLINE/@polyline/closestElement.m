function [ Bp , By , Bd , Bc ] = closestElement( P , x )
%   
% [ piece , xyz_closest_point , distance , index_coordinates ] = closestElement( POLYLINE , [ point1 ; point2 ; ... ] )
% 
% 

  if size( x ,2) ~= nsd( P ), error('points and polyline must have the same NumberSpatialDimensions'); end
  
  x  = double( x );
  nx = size( x ,1);
  
  Bp = zeros( nx ,1);
  Bd =   Inf( nx ,1);
  By =   NaN( nx , size(x,2) );
  Bc =   NaN( nx , 1);
  for p = 1:np( P )
    if isPoint( P , p )
      y = repmat( P.C{p} , [ nx , 1 ] );
      d = sqrt( sum( bsxfun(@minus, P.C{p} , x).^2 ,2) );
      c = ones( nx , 1 );
    else
      [s,y,d,c] = distancePoint2Segments( x , P.C{p} );
      c = s + c;
    end
    
    w = d < Bd;
    Bp(w)   = p;
    By(w,:) = y(w,:);
    Bd(w)   = d(w);
    Bc(w)   = c(w);
  end

end
