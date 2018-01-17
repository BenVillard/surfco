function x = at( P , p , c , mode )
  
  if nargin < 4, mode = 'index'; end


  nc = numel( c );
  if numel(p)==1, p = zeros(nc,1)+p; end
  if numel(p) ~= numel(c), error('pieces and coordinates do not coincide in number'); end
  
  x = NaN( nc , nsd(P) );
  for pp = unique( p(:).' )
    w = p == pp;
    
    switch lower( mode )
      case {'i','index'}
        oc = 1:size( P.C{pp} ,1);
        
      case {'arclength','al','a'}
        oc = [ 0 ; cumsum( sqrt( sum( diff( P.C{pp} ,1,1).^2 ,2) ) ) ];
        
      case {'normalized','w','norm','n'}
        oc = [ 0 ; cumsum( sqrt( sum( diff( P.C{pp} ,1,1).^2 ,2) ) ) ];
        if size( P.C{pp} ,1)>1, oc = oc/oc(end); end
        
    end
    
    y = Interp1D( P.C{pp} , oc , c(w) , 'linear' );

    y( c(w) < oc( 1 ) ,:) = -Inf;
    y( c(w) > oc(end) ,:) = +Inf;
    
    x(w,:) = y;
  end
  
end
