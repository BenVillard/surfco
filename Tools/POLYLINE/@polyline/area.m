function [varargout] = area( P )

  if ~is2D( P ), error('only 2d polylines can be ''booleaned''!.'); end
  if isempty( P ), error('input cannot be an empty polyline'); end
  if ~isSingle( P ), error('only single polylines'); end
  if ~isclosed( P )
    error('input should be an all closed polyline; try with close(B).');
  end
  
  
  P = polygon( P.C{1} );
  [varargout{1:nargout}] = area( P );

end
