function P = close( P , th )

  if nargin < 2, th = Inf; end

  for p = 1:np(P)
    if ~isPoint( P , p ) && ~isClosed( P , p )
      if isinf( th ) ||  sqrt( sum( (P.C{p}(1,:) - P.C{p}(end,:)).^2) ) < th
        P.C{p} = [ P.C{p} ;  P.C{p}(1,:) ];
      end
    end
  end

end
