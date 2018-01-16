function D = ipd( X , Y , d , varargin )

  if nargin < 3, d = 'euclid'; end

  if isempty(Y), Y = X; end
  if ndims( X ) > 2, error('X must be a 2d matrix'); end
  if ndims( Y ) > 2, error('Y must be a 2d matrix'); end
  
  NSD = size(X,2);
  if size(Y,2) ~= NSD, error('number of spatial dimensions must coincide.'); end

  
  if strcmpi( d , 'normal' )
    D = 2*min( asin( ipd( X , Y , 'euclid' ) / 2 ) , asin( ipd( X , -Y , 'euclid' ) / 2 ) );
    return;
  end

  
  szX = size(X,1);
  szY = size(Y,1);
  try,   D = NaN( szX , szY ,'double');
  catch, D = NaN( szX , szY ,'single');
  end
 
  if ~strcmpi( d , 'euclid' )
    X = reshape( X , [ szX , 1 , NSD ] );
  end
  
  ev = min( szY , 5000 );
  while ev >= 1
    try
      for i = 1:ev:szY
        w = ( i ):min( i+ev-1 , szY );
        
        YY = Y(w,1:NSD);
        if ~strcmpi( d , 'euclid' )
          YY = permute( YY , [3 1 2] );
        end

        switch lower(d)
          case 'normal'

          case 'euclid'
            D(:,w) = sqrt( abs( bsxfun( @plus , dot( X , X , 2 ) , dot( YY , YY , 2 ).' ) - 2*( X * YY.' ) ) );
          case 'euclid2'
            D(:,w) = sqrt( sum( bsxfun(@minus,X,YY).^2 ,3) );
          case 'abs'
            D(:,w) =  sum( abs( bsxfun(@minus,X,YY) ) ,3);
          case 'max'
            D(:,w) =  max( abs( bsxfun(@minus,X,YY) ) ,[],3);
          case 'min'
            D(:,w) =  min( abs( bsxfun(@minus,X,YY) ) ,[],3);
          otherwise,
            ev = -2; break;
        end
      end
      break;
    end
    ev = round( ev / 2.001 );
  end
  if ev < 0
    error('invalid distance');
  end
  
end
