function X = ppinv( A , b )

  [m,n] = size( A );

  [U,S,V] = safe_svd( A );

  S = diag( S );

  tol = max(m,n) * eps(max(S));
  r   = sum( S > tol );
  

  if nargin < 2
    
    X = V(:,1:r) * bsxfun( @ldivide , S(1:r) , U(:,1:r).' );
    
  else

%     if r
      X = V(:,1:r) * bsxfun( @ldivide , S(1:r) , U(:,1:r).' * b );
%     else
%       X = zeros(1,size(b,2));
%     end
    
  end



%{
%   if nargin < 2
%     X = ( A.'*A )\A.';
%   else
%     b = A.' * b;
%     X = ( A.'*A )\b;
%   end
%   return;
%   if isempty(A), X = zeros(size(A'),class(A)); return; end

  [m,n] = size(A);
%   try,
%     [U,S,V] = svd(A,0);
%   catch
    [U,S,V] = svdecon(A);
%   end
  if     m > 1,  s = diag(S);
  elseif m == 1, s = S(1);
  else,          s = 0;
  end
  tol = max(m,n) * eps(max(s));
  r = s > tol;
  if sum( r ) < 1
    X = zeros( size(A') , class(A) );
  else
    U = U(:,r);
    s = s( r );
    V = V(:,r);
    if nargin > 1, U = b.' * U; end
    X = V * bsxfun( @ldivide , s(:) , U.' );
  end
  
  
  
  
  function [U,S,V] = svdecon(X)
    [m,n] = size(X);
    
    if  m <= n
      [ U , S ] = eig( X * X.' );
      
      S = diag(S);
      S( S < eps(1) ) = 0;
      [S,order] = sort( S , 'descend' );
      U = U(:,order);
      
      if nargout > 2
        V = X.' * U;
        S = sqrt(S);
        V = bsxfun( @(x,c)x./c , V , S.' );
        S = diag(S);
      end
    else
      C = X'*X;
      [V,S] = eig(C);
      clear C;
      
      [S,order] = sort(abs(diag(S)),'descend');
      V = V(:,order);
      
      U = X*V; % convert evecs from X'*X to X*X'. the evals are the same.
      %s = sqrt(sum(U.^2,1))';
      s = sqrt(S);
      U = bsxfun(@(x,c)x./c, U, s');
      S = diag(s);
    end

  end
%}
end

function [U,S,V] = safe_svd( A )

  OK = false; if OK, return; end
  try,
    [U,S,V] = svd( A , 'econ' );

    OK = true; if OK, return; end
  end
  
  try,
    [V,S,U] = svd( A.' , 'econ' );

    OK = true; if OK, return; end
  end
  
  [m,n] = size(A);
  try
    if  m > n, error('1'); end
    C =  A * A.';  C = ( C + C.' )/2;
    [ U , S ] = eig( C ); U = real( U ); S = real( diag( S ) );

    [ S , order ] = sort( abs(S) , 'descend' );
    w = S > eps(1)*S(1);
    S = S(w);
    U = U( : , order(w) );

    S = sqrt( S );
    if nargout > 2
      V = A.' * U;
%       V = bsxfun( @(x,c)x./c , V , s.' );
      V = bsxfun( @rdivide , V , S.' );
    end
    S = diag( S );

    OK = true; if OK, return; end
  end

  try
    if  m <= n, error('1'); end
    C = A.'*A; C = ( C + C.' )/2;
    [V,S] = eig(C); V = real( V ); S = real( diag( S ) );
    
    [ S , order ] = sort( abs(S) , 'descend' );
    w = S > eps(1)*S(1);
    S = S(w);
    V = V( : , order(w) );
    
    U = A*V; % convert evecs from X'*X to X*X'. the evals are the same.
    S = sqrt(S);
    %U = bsxfun( @(x,c)x./c , U , s' );
    U = bsxfun( @rdivide , U , S.' );
    S = diag(S);

    OK = true; if OK, return; end
  end
  
  matname = fullfile( prefdir , 'no_svd_for_this_matrix.mat' );
  save( matname , 'A' );
  error('cannot compute a descent svd for this matrix.\nCheck it from: "%s"', matname );
  
end
