function hg_ = pplot( P , varargin )
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
  COLOR = [];
  [varargin,~,COLOR ] = parseargs(varargin,'color' ,'$DEFS$',COLOR);
  if isempty( COLOR ), COLOR = 'random'; end
  
  MARKER = [];
  [varargin,~,MARKER] = parseargs(varargin,'marker','$DEFS$',MARKER);
  if isempty( MARKER ), MARKER = 'none'; end

  hg = hggroup('Parent',cax); mZ = 0;
  for p = 1:np(P)
    X = P.C{p}; X(:,end+1:3) = 0; mZ = max( mZ , max( abs(X(:,3)) ) );

    tCOLOR = COLOR;
    if strcmpi( tCOLOR ,'random'), tCOLOR = rand(1,3)*0.5 + 0.2; end
    
    if size( X ,1) == 1 
      switch lower( MARKER )
        case 'none'
          hh = line( X(:,1) , X(:,2) , X(:,3) , varargin{:} , 'Parent',hg,'Color',tCOLOR,'Marker','o','MarkerFaceColor',tCOLOR,'MarkerSize',3);
          line( X(:,1) , X(:,2) , X(:,3) , varargin{:} , 'Parent',hg,'Color',tCOLOR,'Marker','o','MarkerFaceColor','none','MarkerSize',5);
        case '.'
          hh = line( X(:,1) , X(:,2) , X(:,3) , varargin{:} , 'Parent',hg,'Color',tCOLOR,'Marker',MARKER);
          line( X(:,1) , X(:,2) , X(:,3) , varargin{:} , 'Parent',hg,'Color',tCOLOR,'Marker','o','MarkerFaceColor','none','MarkerSize',5);
        otherwise
          hh = line( X(:,1) , X(:,2) , X(:,3) , varargin{:} , 'Parent',hg,'Color',tCOLOR,'Marker',MARKER);
          line( X(:,1) , X(:,2) , X(:,3) , varargin{:} , 'Parent',hg,'Color',tCOLOR,'Marker','o','MarkerFaceColor','none','MarkerSize',get(hh,'MarkerSize')+6);
      end
    elseif ~isequal( X(1,:) , X(end,:) )
      hh = line( X(:,1) , X(:,2) , X(:,3) , varargin{:} , 'Parent',hg,'Color',tCOLOR,'Marker',MARKER);
      line( X( 1 ,1) , X( 1 ,2) , X( 1 ,3) , varargin{:} , 'Parent',hg,'Color',tCOLOR,'Marker','s','MarkerFaceColor',tCOLOR,'MarkerSize',8);
      
      drawOrientation( X , varargin{:} , 'Parent',hg,'Color',tCOLOR,'Marker','s','MarkerFaceColor',tCOLOR );
      
      line( X(end,1) , X(end,2) , X(end,3) , varargin{:} , 'Parent',hg,'Color',tCOLOR,'Marker','x','MarkerFaceColor',tCOLOR,'MarkerSize',8,'LineWidth',2);
    elseif isequal( X(1,:) , X(end,:) )
      hh = line( X(:,1) , X(:,2) , X(:,3) , varargin{:} , 'Parent',hg,'Color',tCOLOR,'Marker',MARKER);
      line( X( 1 ,1) , X( 1 ,2) , X( 1 ,3) , varargin{:} , 'Parent',hg,'Color',tCOLOR,'Marker',' o','MarkerFaceColor',tCOLOR,'MarkerSize',8);
      
      drawOrientation( X , varargin{:} , 'Parent',hg,'Color',tCOLOR,'Marker','o','MarkerFaceColor',tCOLOR );
      
    end
  end
  
  if ~ishold( cax )
    if nd == 3 && mZ > 1e-10, view(3); end
    axis('equal');
  end
  if nargout, hg_ = hg; end

end
function X = interp( X , at )
  L = [ 0 ; cumsum( sqrt( sum( diff(X,1,1).^2 ,2) ) ) ];
  L = L/L(end);
  X = Interp1D( X , L , at , 'linear');
end
function drawOrientation( X , varargin )
  Y = interp( X , 0.05 ); line( Y( 1 ,1) , Y( 1 ,2) , Y( 1 ,3) , varargin{:} ,'MarkerSize',6);
  Y = interp( X , 0.10 ); line( Y( 1 ,1) , Y( 1 ,2) , Y( 1 ,3) , varargin{:} ,'MarkerSize',5);
  Y = interp( X , 0.15 ); line( Y( 1 ,1) , Y( 1 ,2) , Y( 1 ,3) , varargin{:} ,'MarkerSize',4);
  Y = interp( X , 0.20 ); line( Y( 1 ,1) , Y( 1 ,2) , Y( 1 ,3) , varargin{:} ,'MarkerSize',3);
  Y = interp( X , 0.25 ); line( Y( 1 ,1) , Y( 1 ,2) , Y( 1 ,3) , varargin{:} ,'MarkerSize',2);
end