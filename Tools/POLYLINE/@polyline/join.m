function P = join( varargin )

    th = Inf;
    if isnumeric( varargin{end} ) && isscalar( varargin{end} )
        th = varargin{end};
        varargin(end) = [];
    end


  nd = NaN;
  for v = 1:numel( varargin )
    V = varargin{v};
    if ~isa( V , 'polyline' ), error('all inputs must be polylines'); end
    if isnan(nd), nd = nsd( V ); end
    if nd ~= nsd(V), error('all input must be of the same NumberSpatialDimensions'); end
    
    A = V.C{1}; V.C(1) = [];
    while ~isempty( V.C )
      d = NaN( 1 , numel(V.C) );
      for p = 1:np(V)
        d(p) = dP2P( A , V.C{p} );
      end
      [~,p] = min( d );
      
      B = V.C{p}; V.C(p) = [];
      A = joinP( A , B , th );
    end
    
    varargin{v} = A;
  end
  
  P = varargin{1};
  for v = 2:numel( varargin )
    P = joinP( P , varargin{v} , th ); 
  end
  P = polyline( P );
  
end


function [d,id] = dP2P( A , B )
  a0 = A( 1 ,:); a1 = A(end,:);
  b0 = B( 1 ,:); b1 = B(end,:);

  d = sum( [  a0 - b0;...
              a0 - b1;...
              a1 - b0;...
              a1 - b1 ].^2,2);
  [d,id] = min(d);
end
function A = joinP( A , B , th )
  [d,id] = dP2P( A , B );
  if d > th
      A(end+1,:) = NaN;
      A = [ A ; B ];
      return;
  end
  
  switch id
    case 1, A = sflip( A ,1);
    case 2, A = sflip( A ,1); B = sflip( B ,1);
    case 3,
    case 4, B = sflip( B ,1);
  end
  if d == 0
    A = [ A ; B( 2:end ,:) ];
  else
    A = [ A ; B ];
  end
end
function x = sflip( x ,d)
  try,   x = flip( x , d );
  catch, x = flipdim( x , d );
  end
end
