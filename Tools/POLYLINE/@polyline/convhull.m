function P = convhull( P )

  if ~is2D(P), error('convhull is only for 2d polylines'); end

  for p = 1:np( P )
    X = uniqueCoordinates( P ,p);
    n = size( X ,1);

    if n > 2
      K = convhulln( X );
      P.C{p} = X( [ K(:,1) ; K(end,2) ] ,:);
    end
  end

end
