function varargout = maximize2minimize( f , varargin )

  [varargout{1:nargout}] = feval( f , varargin{:} );


  for o = 1:min( nargout , 3 )
    
    varargout{o} = -varargout{o};
    
  end

end