function o = orientation( P , pp )

  if ~is2D(P), error('only 2d polylines have orientation'); end


  if nargin < 2, pp = 1:np( P ); end

  o = NaN( 1 , np( P ) );
  for p = pp
    X = uniqueCoordinates( P ,p);
    n = size( X ,1);

    if n < 3
      o(p) = 0;
    else

      X = convhull( X );
      X(end) = [];
      X = circshift( X , 1-minInd(X) );
      
      if issorted( X )
        o(p) = 1;
      else
        o(p) = -1;
      end
      
%       [~,i2] = max( X(:,2) );
%       if i2 == 1
%         i1 = n; i3 = 2;
%       elseif i2 == n
%         i1 = n-1; i3 = 1;
%       else
%         i1 = i2-1; i3 = i2+1;
%       end
%       cr = cross( [ X(i2,:)-X(i1,:) , 0 ] , [ X(i3,:)-X(i2,:) , 0 ] );
% 
%       o(p) = sign( cr(:,3) );
    end
    
  end

end
