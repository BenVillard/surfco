function P = polylinefun( fcn , P )

  PP = P; PP.C = {};
  for p = 1:P.np
    PP.C{1} = P.C{p};
    PP = feval( fcn , PP );
    P.C{p} = PP.C{1};
  end

end
