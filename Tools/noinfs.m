function x = noinfs( x , v )

  if nargin == 1

    x = x( ~isinf(x) );

  else
    
    x( isinf(x) ) = v;
    
  end
    
end
