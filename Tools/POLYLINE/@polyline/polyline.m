classdef polyline
  properties ( Access = private , Hidden = true )
    C;
  end

  methods ( Hidden = true )
    function P = polyline( varargin )
      if numel(varargin) == 1 && ischar( varargin{1} ) && strcmp( varargin{1} , 'draw' )
        hFig = figure('Name','Draw your polyline   (close this window to finish)');
        
        X = NaN(1,2);
        hL = line('XData',X(:,1),'YData',X(:,2),'Marker','.','LineStyle','-','Color','r');
        hA = get( hL , 'Parent' );
        set( hA , 'XLim',[0 1],'YLim',[0 1],'DataAspectRatio',[1 1 1]);
        
        set( hFig , 'WindowButtonMotionFcn' , @(h,e)moving() );
        set( hFig , 'WindowButtonDownFcn'     , @(h,e)click()  );
        
        
        waitfor( hL );
        
        P = polyline( X );
        figure; pplot( P );
        
        return;
      end
      function moving()
        cp = mean( get( hA , 'CurrentPoint' ) ,1);
        Y = [ vec( get( hL , 'XData' ) ) , vec( get( hL , 'YData' ) ) ];
        Y(end,:) = cp(1:2);
        set( hL , 'XData', Y(:,1) , 'YData', Y(:,2) );
      end
      function click()
        cp = [NaN NaN];
        if pressedkeys(3) == 1
          cp = mean( get( hA , 'CurrentPoint' ) ,1);
        end
        X  = [ X ; cp(1:2) ];
        set( hL , 'XData', X([1:end end],1) , 'YData', X([1:end end],2) );
      end
      
      
      P.C = cell(0,1);
      nsd = -1;
      p = 0;
      for v = 1:numel( varargin )
        x = varargin{v};
        if ~isnumeric(x), error('Numerical coordinates were expected.'); end
        if  isempty(x)  , error('Non empty coordinates were expected.'); end
        if  ndims(x) > 2, error('A matrix (not an array) was expected.'); end
        
        if nsd < 0, nsd = size( x , 2 ); end
        if size(x,2) ~= nsd, error('%d-d coordinates were expected.',nsd); end
        w = all( x(1:end-1,:) == x(2:end,:) ,2);
        x( w ,:) = [];
%         w = any( ~isfinite( x ) ,2);
%         w = ~( w(1:end-1) == w(2:end) );

        w = [ true ; any( ~isfinite( x ) , 2 ) ; true ];
        while ~all(w)
          i0 = find( ~w , 1 ); w(1:i0) = false; 
          i1 = find(  w , 1 ); w(1:i1) = true;
          
          p = p+1;
          P.C{p,1} = x(i0-1:i1-2,:);
        end
      end
    end
    
    function display( P )
      Pn = inputname(1);
      if isempty( Pn ), Pn = 'ans'; end
      fprintf('%s  =  ',Pn);
      disp( P );
    end
    function disp( P )
      if ~numel( P.C )
        fprintf( 'empty Polyline\n' );
      else
        fprintf( '%d-dimensional Polyline with %d segments\n' , nsd(P) , np( P ) );
        for p = 1:np( P )
          fprintf( '%5d vertices' , size( P.C{p} ,1) );
          if     isPoint( P , p )
            fprintf( '  (single point)');
          elseif isClosed( P , p )
            fprintf( '  (closed)');
          end
          fprintf( '\n');
        end
      end
    end
    function n = numel( P )
%       n = np( P );
      n = 1;
    end
    function X = double( P )
      X = double( P.C{1} );
      for p = 2:np( P )
        X(end+1,:) = NaN;
        X = [ X ; P.C{p} ];
      end
    end
    function out = ispolyline( P )
      if isa(P,'polyline'), out = true; else, out = false; error('why here?'); end
    end
    function SS  = uneval( P )
      S = sprintf( 'polyline(' );
      S = [ S , ' ' , uneval( P.C{1} ) ];
      for p = 2:np( P )
        S = [ S , ' , ' , uneval( P.C{p} ) ];
      end
      S = [ S , ' )' ];
      if nargout > 0, SS = S;
      else
        Sn = inputname(1);
        if isempty( Sn  ), fprintf('%s\n', S );
        else,              fprintf('%s = %s\n', Sn , S );
        end
      end
    end
    function out = isidentical( A , B )
      out = false;
      if ~ispoly(A),   return; end
      if ~ispoly(B),   return; end
      if ~isequal( numel(A.C) , numel(B.C) ),  return; end
      if ~isequal( size(A.C{1},2) , size(B.C{1},2) ),  return; end
      for p = 1:np(A)
        if ~isequal( A.C{p} , B.C{p} ),  return; end
      end
      out = true;
    end
%     function out = isequal( A , B )
%       switch char( 48 + [ isa(A,'polyline') , isa(B,'polyline') ] )
%         case '00', error('why here?');
%         case '01', error('A polyline can only be compared against a polyline.'); out = isequal( A   , B.C );
%         case '10', error('A polyline can only be compared against a polyline.'); out = isequal( A.C , B   );
%         case '11', out = isidentical( A , B );
%       end      
%     end
    
    %%queries: usually results in arrays with a single value per piece
    function sz = size( P , varargin )
      sz = [ np(P) , nsd(P) ];
%       sz = NaN( 1 , np(P) );
%       for p = 1:np(P)
%         sz(p) = size( P.C{p} ,1);
%       end
      sz = sz( varargin{:} );
    end
    

    function P = vertcat( varargin ),   P = cat(1,varargin{:}); end
    function P = horzcat( varargin ),   P = cat(1,varargin{:}); end
    function A = cat( ~ , A , varargin )
      if numel( varargin ) == 0, return; end
      nd = nsd( A );
      for v = 1:numel( varargin )
        B = varargin{v};
        if      ispoly( B )
        elseif  isnumeric( B )
          B = polyline( B );
        end
        if isnan( nd ), nd = nsd(B); end
        if nsd( B ) ~= nd
          error('inconsistent number of spatial dimensions.');
        end
        A.C = vertcat( A.C , B.C );
      end
    end

    function P = uplus( P )
      return;
    end
    function P = uminus( P )
      for p = 1:np(P), P.C{p} = -P.C{p}; end
    end
    function P = times( P , f )
      if ispoly( f ) && isnumeric( P )
        P = times( f , P );
        return;
      end
      if ~isnumeric( f ), error('a numeric factor was expected'); end
      f = double( f );
      if      numel( f ) == 1
        for p = 1:np(P), P.C{p} = f * P.C{p}; end
      elseif  size( f ,1) == 1 && numel( f ) == nsd( P )
        for p = 1:np(P), P.C{p} = bsxfun( @times , f , P.C{p} ); end
      else
        error('invalid multiplication');
      end
    end
    function P = mtimes( T , P )
      if ispoly( T ) && isnumeric( P )
        P = times( P , T );
        return;
      end
      if ~isnumeric( T ), error('a numeric factor was expected'); end
      T = double( T );
      nd = nsd( P );
      if      numel( T ) == 1
        P = times( P , T );
      elseif  size( T ,1) == 1 && numel( T ) == nd
        P = times( P , T );
      elseif  size( T ,1) == size( T ,2) && size( T ,1) == nd
        for p = 1:np(P)
          P.C{p} = P.C{p} * T.';
        end
      elseif  size( T ,1) == size( T ,2) && size( T ,1) == nd+1 && ~any( T(nd+1,1:nd) ) && T(nd+1,nd+1) == 1
        for p = 1:np(P)
          P.C{p} = bsxfun( @plus , P.C{p} * T(1:nd,1:nd).' , T(1:nd,nd+1).' );
        end
      else
        error('invalid multiplication');
      end
    end
    function P = plus( P , s )
      if ispoly( s ) && isnumeric( P )
        P = plus( s , P );
        return;
      end
      if ~isnumeric( s ), error('a numeric summand was expected'); end
      s = double( s );
      if      numel( s ) == 1
        for p = 1:np(P), P.C{p} = s + P.C{p}; end
      elseif  size( s ,1) == 1 && numel( s ) == nsd( P )
        for p = 1:np(P), P.C{p} = bsxfun( @plus , s , P.C{p} ); end
      else
        error('invalid sumation');
      end
    end
    function P = minus( P , r )
      if ~isnumeric( r ), error('a numeric --summand was expected'); end
      r = double( r );
      if      numel( r ) == 1
        for p = 1:np(P), P.C{p} = P.C{p} - r; end
      elseif  size( r ,1) == 1 && numel( r ) == nsd( P )
        for p = 1:np(P), P.C{p} = bsxfun( @minus , P.C{p} , r ); end
      else
        error('invalid sumation');
      end
    end
    
    function P = flip( P , pp )
      if nargin < 2, pp = 1:np(P); end
      for p = 1:np(P);
        try,    P.C{p} = flip(    P.C{p} , 1 );
        catch,  P.C{p} = flipdim( P.C{p} , 1 );
        end
      end
    end
    function P = transform( P , varargin )
      if numel(varargin) == 1 && isnumeric( varargin{1} )
        T = varargin{1};
      else
        T = maketransform( varargin{:} );
      end
      
      nd = nsd(P);
      for p = 1:np(P)
        if nd == 2 && size( T ,1) == 4
          P.C{p}(:,3) = 0;
        end
        P.C{p} = bsxfun( @plus , P.C{p} * T(1:3,1:3).' , T(1:3,4).' );
      end
    end
    function varargout = set( h , varargin )
      if ~strcmp( get( h , 'Type' ) , 'line' )
        error('only lines can be set');
      end
      [varargin,i,P] = parseargs(varargin,'data','$DEFS$',[]);
      if ~i || isempty( P )
        error('only property ''Data'' can be set by a polyline');
      end
      X = double( P );
      try, 
        [ varargout{1:nargout} ] = set(h, varargin{:} ,'XData',X(:,1),'YData',X(:,2),'ZData',X(:,3));
      catch
        set(h, varargin{:} ,'XData',X(:,1),'YData',X(:,2),'ZData',X(:,3));
      end
    end
    
    
    
    
    
    %functions defined in separate files.
    varargout =         subsref( P , varargin );
    varargout =        subsasgn( P , varargin );
    varargout =            plot( P , varargin );
    varargout =          hpplot( P , varargin );
  end

  methods ( Access = public , Hidden = false )
    %%general queries: global values for the whole polyline
    function o = NumberOfSpatialDimensions( P )
      o = nsd(P);
    end
    function o = NumberOfPieces( P )
      o = np( P );
    end
    function L = Length( P )
      if ~np(P)
        L = [];
      else
        L = sum( length(P) );
      end
    end
    function L = ArcLength( P , normalized )
      if nargin < 2, normalized = false; end
      if ~np(P)
        L = [];
      else
        a = arclength( P );
        L = a{1};
        for p = 2:numel(a)
          L = [ L , L(end)+a{p} ];
        end
        if normalized, L = L/L(end); end
      end
    end

    %%piece-wise queries: results in arrays with a single value per piece
    function o = NumberOfNodes( P )
      o = nn( P );
    end
    function o = isclosed( P , varargin )
      o = false( 1 , np(P) );
      for p = 1:np(P)
        o(p) = isClosed( P , p );
      end
      o = o( varargin{:} );
    end
    function o = length( P )
      o = NaN( 1 , np(P) );
      for p = 1:np(P)
        o(p) = sum( sqrt( sum( diff( P.C{p} ,1,1).^2 ,2) ) );
      end
    end
    
    %%methods for single polylines
    function o = coordinates( P )
      o = P.C.';
    end
    function o = arclength( P , normalized )
      if nargin < 2, normalized = false; end
      o = cell( 1 , np(P) );
      for p = 1:np(P)
        o{p} = [ 0 ; cumsum( sqrt( sum( diff( P.C{p} ,1,1).^2 ,2) ) ) ].';
        if normalized && o{p}(end) > 0, o{p} = o{p}/o{p}(end); end
      end
    end
    
    
    
  end
 
  methods ( Access = private , Hidden = true )
    function o  = isPlanar(P)
      o = size( P.C{1} ,2) == 2;
      if o, return; end
      XYZ = vertcat( P.C{:} );
      XYZ = bsxfun( @minus , XYZ , mean(XYZ,1) );
      [~,d] = svd( XYZ ,0);
      o = d(end) < 1e-8;
    end
    
    function o  = is2D(P),   o = size( P.C{1} ,2) == 2; end
    function nd = nsd( P ), nd = NaN; try, nd = size( P.C{1} ,2); end; end  %NumberOfSpatialDimensions
    function n  = np( P ),   n = numel( P.C );          end                 %NumberOfPieces
    function o  = nn( P , p )                                               %NumberOfNodes
      if nargin > 1
        o = size( P.C{p} ,1);
      else
        o = zeros( 1 ,np(P) );
        for p = 1:np(P), o(p) = size( P.C{p} ,1); end
      end
    end

    function X  = uniqueCoordinates( P , p )
      X = P.C{p};
      w = all( ~diff(X,1,1) ,2);
      X(w,:) = [];
      if size(X,1)>1 && isequal( X(1,:) , X(end,:) )
        X(end,:) = [];
      end
    end
    function o  = isSingle( P )
      o = numel( P.C ) == 1;
    end
    function o  = isPoint( P , p )
      o = size( P.C{p} ,1) == 1;
    end
    function o  = isClosed( P , p )
      o = numel( P.C{p} ) > 1 && isequal( P.C{p}(1,:) , P.C{p}(end,:) );
    end
  end
end

function o = ispoly( x ),   o = isa( x , 'polyline' );        end
