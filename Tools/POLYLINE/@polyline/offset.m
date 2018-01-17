function O = offset( P , d )

  nd = nsd(P);
  if nd ~= 2  &&  nd ~= 3, error('only 2d or planar 3d polylines'); end
  if ~isPlanar( P ), error('only for 2D polyline'); end
  if ~isSingle( P ), error('only single polylines'); end

  if nd == 3
    [Z,iZ] = getPlane( P.C{1} );
    P.C{1} = bsxfun( @plus , P.C{1}*iZ(1:3,1:3).' , iZ(1:3,4).' );
    P.C{1} = P.C{1}(:,1:2);
  end
  

  if isclosed( P )

    P = polygon( P.C{1} );
    P = offset( P , d , 'ROUNDed' );

    O = polyline();
    for p = 1:size(P,1)
      PP = polyline( P(p).XY );
      PP = close( PP );
      O = [ O ; PP ];
    end
    
  else
    
    S  = P;
    P0 = P.C{1}( 1 ,:);
    P1 = P.C{1}(end,:);
    
    P = polygon( P.C{1} );
    P = offset( P , d , 'ROUNDed' );
    
    mp = @(p) ( p(1:end-1,:) + p(2:end,: ) )/2;
    O = polyline();
    for p = 1:size(P,1)
      PP = polyline( P(p).XY );

      [ ~ , ~ , ~ , c ] = closestElement( PP , [ P0 ; P1 ] );
      c = sort(c);
      A = PP.part( [ 1 , c(1) ] );     [~,~,dA] = closestElement( S , mp( A.C{1} ) ); dA = max( dA );
      B = PP.part( c );                [~,~,dB] = closestElement( S , mp( B.C{1} ) ); dB = max( dB );
      C = PP.part( [ c(2) , Inf ] );   [~,~,dC] = closestElement( S , mp( C.C{1} ) ); dC = max( dC );
      [dd,m] = max( [ dA , dB , dC ] );
      if dd > d * 1.01
        switch m
          case 1, PP = join( [ B ; C ] );
          case 2, PP = join( [ A ; C ] );
          case 3, PP = join( [ A ; B ] );
        end
      end
      O = [ O ; PP ];
    end
    
  end

  if nd == 3
    O.C{1}(:,3) = 0;
    O.C{1} = bsxfun( @plus , O.C{1}*Z(1:3,1:3).' , Z(1:3,4).' );
  end

end
