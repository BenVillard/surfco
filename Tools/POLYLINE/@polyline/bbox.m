function bb = bbox( P )

  nd = nsd( P );
  bb = [ Inf(1,nd) ; -Inf(1,nd) ];
  
  for p = 1:np(P)
    bb(1,:) = min( bb(1,:) , min( P.C{p} , [] , 1 ) );
    bb(2,:) = max( bb(2,:) , max( P.C{p} , [] , 1 ) );
  end


end
