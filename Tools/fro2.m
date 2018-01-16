function n = fro2( x , dim )

  if nargin < 2
    
    n = x(:).'*x(:);
    
  else
    
    n = sum( x.^2 , dim );
    
  end

end
