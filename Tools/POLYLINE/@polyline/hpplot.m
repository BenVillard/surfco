function varargout = hpplot( varargin )

  ax = gca;
  ish = ishold( ax );
  if ~ish
    CLEANUP = onCleanup(@()hold(ax,'off'));
    if numel(get(ax,'Children')), hold(ax,'on'); end
  end
  
  [ varargout{1:nargout} ] = pplot( varargin{:} );
  
end
