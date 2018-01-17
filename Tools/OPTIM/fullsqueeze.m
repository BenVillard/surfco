function x = fullsqueeze( x )

    x = squeeze(x);
    if isvector(x)
      x = x(:);
    end
    
end
