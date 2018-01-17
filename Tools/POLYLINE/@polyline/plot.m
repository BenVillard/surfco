function hg_ = plot( P , varargin )
%

  nd = nsd(P);
  if nd < 2 || nd > 3
    error('only 2d or 3d polylines can be ploted');
  end


  [varargin,~,parent] = parseargs( varargin , 'parent','$DEFS$',[]);
  if isempty( parent ), parent = gca; end
  cax = ancestor( parent , 'axes' );
  cax = newplot(cax);

  
  varargin = getLinespec( varargin );
  

  X = double( P );
  switch nd
    case 2
      hg = plot( X(:,1) , X(:,2) , varargin{:} , 'Parent', cax );
    case 3
      hg = plot3( X(:,1) , X(:,2) , X(:,3) , varargin{:} , 'Parent', cax );
  end

  if ~ishold( cax )
    if nd == 3, view(3); end
    axis('equal');
  end
  if nargout, hg_ = hg; end

end
