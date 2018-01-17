function P = orient( P )

  if ~is2D(P), error('only 2d polylines have orientation'); end

  for p = 1:np(P)
    if orientation( P , p ) < 0
      P.C{p} = double( flip( P , p ) );
    end
  end
  
end