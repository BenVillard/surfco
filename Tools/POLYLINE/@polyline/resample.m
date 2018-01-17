function [P,nt] = resample( P , mode , nt )

  if ~isSingle( P ), error('only single polylines can be resamplig'); end

  if ~ischar( mode )
    nt    = mode;
    mode = 'index';
  end
  if     mode( 1 ) == '+', addOriginal = true; mode( 1 ) = [];
  elseif mode(end) == '+', addOriginal = true; mode(end) = [];
  else,                    addOriginal = false;
  end
    
  switch lower( mode )
    case {'index','i'}
      if ~issorted(nt), error('new points should be sorted'); end
      ot = 1:size( P.C{1} ,1);
      
    case {'arclength','al'}
      if ~issorted(nt), error('new points should be sorted'); end
      ot = [ 0 ; cumsum( sqrt( sum( diff( P.C{1} ,1,1).^2 ,2) ) ) ];
      
    case {'normalized','w','norm','n'}
      if ~issorted(nt), error('new points should be sorted'); end
      ot = [ 0 ; cumsum( sqrt( sum( diff( P.C{1} ,1,1).^2 ,2) ) ) ];
      if size( P.C{1} ,1)>1, ot = ot/ot(end); end
      
    case {'every','e'}
      if numel(nt) ~= 1, error('a scalar radius is expected'); end
      if nt <= 0, error('a positive scalar radius is expected'); end
      R2 = nt^2;
      ot = 1:size( P.C{1} ,1);
      
      c = P.C{1}(1,:); nt = 1;
      while 1
        t0 = floor(nt(end));
        t1 = find( sum( bsxfun(@minus, c , P.C{1}( t0+1:end ,:) ).^2 ,2) > R2 ,1);
        if isempty( t1 ), break; end
        t1 = t0 + t1; p1 = P.C{1}(t1,:);
        t0 = t1-1;    p0 = P.C{1}(t0,:);
        
        %tt = fzero( @(tt)R2-sum( ( c - Interp1D( P.C{1}(t0:t1,:) , t0:t1 , max(nt(end),tt),'closest') ).^2 ) , nt(end) );
        tt = circle_segment( c(:) , p0(:) , p1(:) , R2 );
        c  = p0 + tt * ( p1 - p0 );
        nt = [ nt , tt + t0 ];
      end
      
    otherwise
      error('unknown resampling mode');
  end
  
  if nt( 1 ) < ot( 1 ), error('first point before the first point'); end
  if nt(end) > ot(end), error('last point after the last point'); end
  if addOriginal
    nt = sort( [ nt(:) ; ot(:) ] );
  end
  try
    P.C{1} = Interp1D( P.C{1} , ot , nt , 'linear','closest');
  catch
    for d = 1:size( P.C{1} , 2 )
      O(:,d) = interp1( ot(:) , P.C{1}(:,d) , nt(:) , 'linear' );
    end
    P.C{1} = O;
  end


end

function t = circle_segment( c , p , q , R2 )

r  = q - p;
rr = r.' * r;
d  = p - c;

t = ( sqrt( ( d.'*r )^2 + ( R2 - d.'*d )*rr ) - d.'*r )/rr;

end
