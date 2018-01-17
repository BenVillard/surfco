function x = nonans( x , v , varargin )

  if nargin == 1

    x = x( ~isnan(x) );

  elseif isnumeric(v)
    
    x( isnan(x) ) = v;
    
  elseif ischar(v)

    switch lower(v)
      case 'mingradient'
        [varargin,RETURN_OP] = parseargs( varargin , 'OPERATOR','$FORCE$',{true,false} );
        [varargin,i,initial] = parseargs( varargin , 'INITial','$DEFS$',[] );

        BM = nan;
        [varargin,BM] = parseargs( varargin , 'symmetric' , '$FORCE$',{'symmetric' ,BM} );
        [varargin,BM] = parseargs( varargin , 'circular'  , '$FORCE$',{'circular'  ,BM} );
        [varargin,BM] = parseargs( varargin , 'replicate' , '$FORCE$',{'replicate' ,BM} );
        
        if isempty(varargin)
          varargin{1} = ndims(x);
        end
        
        L = [];
        for d = 1:varargin{1}
          if size(x,d) < 2, continue; end
          L = [ L ; imfiltermtx( {size(x)}, vec([-1;1],d),BM,'conv','same') ];
        end

        if RETURN_OP,   x = minnorm( x , L , 'OPERATOR');
        else        ,   x = minnorm( x , L , 'INITIAL' , initial );
        end

      case 'minlaplacian'
        [varargin,RETURN_OP] = parseargs( varargin , 'OPERATOR','$FORCE$',{true,false} );
        [varargin,i,initial] = parseargs( varargin , 'INITial','$DEFS$',[] );
        

        BM = nan;
        [varargin,BM] = parseargs( varargin , 'symmetric' , '$FORCE$',{'symmetric' ,BM} );
        [varargin,BM] = parseargs( varargin , 'circular'  , '$FORCE$',{'circular'  ,BM} );
        [varargin,BM] = parseargs( varargin , 'replicate' , '$FORCE$',{'replicate' ,BM} );

        if isempty(varargin)
          varargin{1} = ndims(x);
        end
        
        L = zeros( [ 3*ones(1,ndims(x)) , 1 ] );
        for d = 1:varargin{1}
          if size(x,d) < 3, continue; end
          L = L + resize( vec( [-2 1 1] , d ) , size(L) );
        end
        L = fftshift( L );
        
        L = imfiltermtx( { size(x) } , L , BM , 'conv' , 'same' );

        
        if RETURN_OP,   x = minnorm( x , L , 'OPERATOR');
        else        ,   x = minnorm( x , L , 'INITIAL' , initial );
        end

      case 'minnorm'
        x = minnorm( x , varargin{:} );
      
      case 'dct'
        x = dctinpaint( x , varargin{:} );

      case {'chessboard','cityblock','euclidean','quasi-euclidean'}
        [xx,inds] = crop( x , 1 , isnan(x) );
        x(inds{:}) = bwdistinpaint( xx , lower(v) );

      case {'inpaint'}
        switch ndims( x )
          case {1,2}
            [xx,inds] = crop( x , 2 , isnan(x) );
            x(inds{:}) = inpaint_nans( xx , varargin{:} );
          case 3
            [xx,inds] = crop( x , 2 , isnan(x) );
            x(inds{:}) = inpaint_nans3( xx , varargin{:} );
          otherwise
            error('inpaint only valid for 1, 2 or 3 dims');
        end
        
      otherwise
        error('unknown method');

    end

  end

  function x = minnorm( x , Lop , varargin )
    [varargin,RETURN_OP] = parseargs( varargin , 'OPERATOR','$FORCE$',{true,false} );
    [varargin,i,initial] = parseargs( varargin , 'INITial','$DEFS$',[] );

    fixed = ~isnan( x(:) );
    if all(fixed), return; end
    
    if ~issparse( Lop )
      Lop = imfiltermtx( { size( x ) } , Lop , varargin{:} );
    end
    if size( Lop , 2 ) ~= numel(x), error('incorrect Lop sizes'); end
    
    if ~RETURN_OP
    
      x( ~fixed ) = 0;

  %     C = sparse( find( ~fixed ) , 1:nnz( ~fixed ) , 1 , size(Lop,1) , nnz( ~fixed ) );
  %     x( ~fixed ) = -( Lop * C ) \ ( Lop* x(:) );

  %     x( ~fixed ) = - ( Lop(:,~fixed) ) \ ( Lop* x(:) );

  %     hay que usar lsqr (si esta bien condicionada)  o  qmr (si no is well conditioned )
  
      if ~isempty(initial)
        if ~isequal( size(initial) , size(x) )
          error('initial have to be of X''s size');
        end
        initial = initial( ~fixed );
      end
  
%       [ res ,flag ] = lsqr( Lop(:,~fixed) ,  - ( Lop * x(:) ) , 1e-8 , 1000 ,[],[], initial(:) );
      [ res ,flag ] = lsqr( Lop(:,~fixed) ,  - ( Lop * x(:) ) , 1e-9 , 1e6 ,[],[], initial(:) );
      if any( ~isfinite(res(:)) )
        res = ( Lop(:,~fixed) ) \ (  - ( Lop * x(:) ) );
      else
        if flag, warning('not converged. se necesitan mas iteraciones'); end
      end
      
      x( ~fixed ) = res;
      
    else
      
      A = Lop( :,~fixed);
      
      Lop(:, ~fixed ) = 0;
      
      x = speye( numel(x) );
      try
        try
          pA = pseudoinverse( A );
        catch
          pA = pseudoinverse( A , [] , 'lsqr' );
        end
        x( ~fixed , : ) = - ( pA * Lop );
      catch
        try
          x( ~fixed , : ) = - inverse( A.' * A ) * ( A.' * Lop );
        catch
          x( ~fixed , : ) = - inv( A.' * A ) * ( A.' * Lop );
        end
      end
      
    end

% M = randn(100);
% x = randn(size(M,2),1);
% b = M*x;
% [Q,R,E] = qr(M,0);
% iE = setv( E*0 , E , 1:numel(E) );
% plot(  getv( R\(Q.'*b) , iE )  - x )
    
  end
  
  

  function y = bwdistinpaint( x , mode )
    y = nan(size(x));
    for d = [ 0 , ndims(x):-1:1 ]
      if d > 0
        x = flipdim( x , d ); 
        y = flipdim( y , d ); 
      end
      [D,L] = bwdist( ~isnan( x ) , mode );
      L( L > 16777215 ) = 0;
      y( ~~L ) = x( L( ~~L ) );
      if d > 0
        x = flipdim( x , d ); 
        y = flipdim( y , d ); 
      end
      if ~any( isnan(y(:)) ), break; end
    end
    
    if any( isnan( y(:) ) )
      error('aun quedan nans!!');
    end
  end
  
  
  function y = dctinpaint( x , N , y )
 
    W = ~isnan( x );
    if all( W(:) )
      y = x;
      return;
    end

    x = double(x);
    
    if nargin < 2, N = 100; end

    sizx = size(x);
    d = ndims(x);
    Lambda = zeros(sizx);
    for i = 1:d
      siz0 = ones(1,d);
      siz0(i) = sizx(i);
      Lambda = bsxfun(@plus,Lambda,...
        cos(pi*(reshape(1:sizx(i),siz0)-1)/sizx(i)));
    end
    Lambda = 2*( Lambda - d );

    if nargin < 3
      y = zeros(sizx);
      y( ~W ) = mean( vec( x( W ) ) );
    elseif ischar( y )
      y = nonans( x , y );
    elseif isscalar( y )
      y = zeros(sizx) + y;
    end
    x( ~W ) = 0;
    y( W  ) = x( W );

    % Smoothness parameters: from high to negligible values
    s = logspace( 3 , -3 , N );

    RF = 2; % relaxation factor
    Lambda = Lambda.^2;

    for i = 1:N
      Gamma = 1./( 1 + s(i)*Lambda );
      y = RF * idctn( Gamma .* dctn( W .* (x-y) +y ) ) + ( 1 - RF ) * y;
    end
    y( W ) = x( W );
  
    
  end

  
end
