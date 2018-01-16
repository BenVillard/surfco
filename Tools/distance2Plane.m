function [d,xy] = distance2Plane( xyz , P , signed )

  if nargin < 3, signed = false; end
  signed = ~~signed;
  
  [P,iP] = getPlane( P );
  
  d = NaN;
  try
    d = xyz * iP(3,1:3).' + iP(3,4);
    d = d(:);
    if isempty(d), d = NaN; end
  end
  if nargin > 1
    xy = transform( xyz , iP );
    xy(:,3) = 0;
    xy = transform( xy , P );
  end

  if ~signed, d = abs( d ); end
  
end
