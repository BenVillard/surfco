function [x,pa,pb,ca,cb] = intersection( A , B )


  if ~isa( A , 'polyline' ), error('A must be a polyline'); end
  if isempty( A ), error('A cannot be an empty polyline'); end
    
  if ~isa( B , 'polyline' ), error('B must be a polyline'); end
  if isempty( B ), error('B cannot be an empty polyline'); end

  nd = nsd(A);
  if nd ~= nsd( B ), error('polyline must be of the same NumberSpatialDimensions'); end
  
  x = []; pa = []; pb = []; ca = []; cb = [];
  for a = 1:np(A)
    for b = 1:np(B)
      [d,aa,bb] = distanceSegment2Segment( A.C{a} , B.C{b} );
      [i,j] = ndgrid( 1:size( A.C{a} ,1)-1 , 1:size( B.C{b} ,1)-1 );
      w = d < 1e-10;
      d = d(w); aa = aa(w) + i(w); bb = bb(w) + j(w);
      [~,ord] = sort( aa );
      d = d(ord); aa = aa(ord); bb = bb(ord);
      
      y = ( at( A , a , aa , 'i' ) + at( B , b , bb , 'i' ) )/2;
      
      x  = [ x ; y ];
      pa( end+1:size(x,1) ,1) = a;
      pb( end+1:size(x,1) ,1) = b;
      ca = [ ca ; aa ];
      cb = [ cb ; bb ];
    end
  end
  


end
