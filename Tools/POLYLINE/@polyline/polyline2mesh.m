function M = polyline2mesh( P )

  if ~isSingle(P),  error('only for single polylines'); end
  if  isPoint(P,1), error('only for non singular polylines'); end
  
  M.xyz = P.C{1};
  
  N = size( M.xyz ,1);
  M.tri = [ 1:N-1 ; 2:N ].';
  if isClosed( P ,1 )
    M.tri = [ M.tri ; N 1 ];
  end
  
end
