function P = circshift( P , k )

  if ~isSingle( P ), error('only single polylines'); end
  if ~isclosed( P )
    error('input should be an all closed polyline; try with close(B).');
  end

  if ischar( k )
    switch k
      case 'minY', [~,k] = min( P.C{1}(:,2) ); k = 1-k;
      case 'maxY', [~,k] = max( P.C{1}(:,2) ); k = 1-k;
      case 'minX', [~,k] = min( P.C{1}(:,1) ); k = 1-k;
      case 'maxX', [~,k] = max( P.C{1}(:,1) ); k = 1-k;
    end
  end
  
  P.C{1}(end,:) = [];
  
  P.C{1} = circshift( P.C{1} , [ k , zeros(1,nsd(P)) ] );
  
  P.C{1}(end,:) = P.C{1}( 1 ,:);
  

end
