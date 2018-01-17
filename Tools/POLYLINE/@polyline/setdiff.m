function C = setdiff( A , B )

  if ~isa( A , 'polyline' ), error('A must be a polyline'); end
  if ~isa( B , 'polyline' ), error('B must be a polyline'); end

  if ~isPlanar( A ), error('only planar polylines can be ''booleaned''!.'); end
  if ~isPlanar( B ), error('only planar polylines can be ''booleaned''!.'); end

  if isempty( A ), error('A cannot be an empty polyline'); end
  if isempty( B ), error('B cannot be an empty polyline'); end
  
%   if ~all( isclosed( A ) ), error('A should be an all closed polyline; try with close(A).'); end
  if ~all( isclosed( B ) ), error('B should be an all closed polyline; try with close(B).'); end

  ndA = nsd(A); if ndA ~= 2  &&  ndA ~= 3, error('A must be 2d or planar 3d polyline.'); end
  ndB = nsd(B); if ndB ~= 2  &&  ndB ~= 3, error('B must be 2d or planar 3d polyline.'); end
  if ndA ~= ndB, error('NumberOfSpatialDimensions must coincide.'); end

  Z = eye(4);
  if ndA == 3
    [Z,iZ] = getPlane( A.C{1} );
    for c = 1:np(A)
      A.C{c} = bsxfun( @plus , A.C{c}*iZ(1:3,1:3).' , iZ(1:3,4).' );
      if max(abs(A.C{c}(:,3))) > 1e-12, error('imposible to convert A to 2d'); end
      A.C{c} = A.C{c}(:,1:2);
    end
    for c = 1:np(B)
      B.C{c} = bsxfun( @plus , B.C{c}*iZ(1:3,1:3).' , iZ(1:3,4).' );
      if max(abs(B.C{c}(:,3))) > 1e-12, error('imposible to convert B to 2d'); end
      B.C{c} = B.C{c}(:,1:2);
    end
  end
  
    
  AA = polygon( A.C{1} );
  for p = 2:np( A ), AA = union( AA , polygon( A.C{p} ) ); end

  BB = polygon( B.C{1} );
  for p = 2:np( B ), BB = union( BB , polygon( B.C{p} ) ); end


  CC = setdiff( AA , BB );
%   hplot( CC , '.:g' )
  
  M = linspace(0,1,9); M = [ M ; fliplr(M) ].';
  C = polyline();
  for p = 1:size( CC ,1)
    P = CC(p).XY; nP = size( P ,1);
    P = struct('xyz',P,'tri', [ ( 1:nP ).' , [2:nP 1].' ]);
    
    for s = 1:size( P.tri ,1)
      [is,id]=ismember( P.xyz( P.tri(s,:) ,:) , A.C{1} , 'rows' );
      if all(is) && abs(diff( id )) == 1, continue; end
      [~,~,d] = closestElement( A , M*P.xyz( P.tri(s,:) ,:) );
      if max(d) > 1e-10
        P.tri(s,:) = 0;
      end
    end
    P.tri( ~any( P.tri ,2) ,:) = [];
    if isempty( P.tri ), continue; end
    PP = tri2Polyline( P.tri );
    
    for pp = 1:numel(PP)
      C = [ C , polyline( P.xyz( PP{pp} ,: ) ) ];
    end
  end

  if ndA == 3
    for c = 1:np(C)
      C.C{c}(:,3) = 0;
    end
  end
  if ~isequal( Z , eye(4) )
    for c = 1:np(C)
      C.C{c} = bsxfun( @plus , C.C{c}*Z(1:3,1:3).' , Z(1:3,4).' );
    end
  end
    
end


function C = tri2Polyline( E )
  C = {};
  while numel(E)
    tC = E(1,:); E(1,:) = [];
    while 1
      w = find( any( E == tC(end) ,2) ,1);
      if isempty(w), break; end
      tC = [ tC , setdiff( E(w,:) , tC(end) ) ]; E(w,:) = [];
    end
    while 1
      w = find( any( E == tC(1) ,2) ,1);
      if isempty(w), break; end
      tC = [ setdiff( E(w,:) , tC(1) ) , tC ]; E(w,:) = [];
    end
    C{end+1} = tC;
  end
end




