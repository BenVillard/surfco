function [X,E] = ExhaustiveSearch( F , X , R , D , varargin )
% 
% 
% ExhautiveSearch( fcn , X0 , radius , divisions , 'verbose' )
% 
% 
%


  if nargin < 4 || isempty( D )
    D = 3;
  end

  if ~isscalar( R ), error('radius R has to be an scalar'); end
  R_init = R;

  N = numel(X);
  sz = size(X);
  Rsz = @(x)reshape(x,sz);
  
  if isscalar(D)
    D = ones(N,1)*D;
  end
  if any(D < 2)  , error('number of divisions have to be greater than 1.'); end
  
  shrinkFactor = min(D);
  shrinkFactor = (shrinkFactor-1)/(shrinkFactor);
  
  [varargin,i,MAX_ITERATIONS] = parseargs(varargin,'maxITERATIONS','$DEFS$',Inf);
  [varargin,i,MAX_TIME      ] = parseargs(varargin,'maxTIME'      ,'$DEFS$',Inf);
  IS_VECTORIZED  = true;
  VERBOSE        = false;
  
  [varargin,VERBOSE] = parseargs(varargin,'verbose','$FORCE$',{true,VERBOSE});
  
  FCN = [];
  [varargin,~,FCN] = parseargs(varargin,'fcn','$DEFS$',FCN);
  
  
  IT = 0;
  START_TIME = clock;
  onBoundary = true;
  typical_sz = [];
  while 1

    for c = 1:N
      xx{c} = linspace( X(c) - R , X(c) + R , D(c) );
      xx{c} = unique( [ X(c) xx{c} ] );
    end
      
    Xs = ndmat_mx( xx{:} );
    
    if ~onBoundary %|| rem(IT,10) ~= 0

      sz = cellfun('prodofsize',xx);
        
      if isequal( sz , typical_sz )
        Z = typical_Z;
      else
        Z  = false( [sz,1] );
        
        for c = 1:N
          indexs{c} = unique( ceil( ( sz(c) + [ 0 1 ] )/2 ) );
        end
        for c = 1:N
          Z( indexs{1:c-1} , : , indexs{c+1:end} ) = true;
        end
        if isempty( typical_sz )
          typical_sz = sz;
          typical_Z  = Z;
        end
      end

      Xs = Xs( Z(:) , : );
    end
    Xs = [ Xs ; X(:).' ];
    [~,i] = unique( Xs , 'rows' ,'first' );
    Xs = Xs( sort(i) , : );
    

    if size( Xs , 1 ) == 1
      X = Rsz( Xs(1,:) );
      break;
    end

    if IS_VECTORIZED    
      try                    %vectorized way
        E = F( Xs );
        if numel( E ) ~= size( Xs , 1 )
          error('no vectorized function');
        end
      catch 
	      IS_VECTORIZED = false;
        continue;
      end
    else                  %not vectorized way

      E = zeros( size(Xs,1) , 1 );  for i = 1:size(Xs,1), E(i) = F( Rsz( Xs(i,:) ) ); end

    end
    
    
    d = sum( bsxfun( @minus , Xs , X(:).' ).^2 ,2);

    
    [~,id] = sortrows( [ E , d ] ); id = id(1);
    E = E( id );
    X = Rsz( Xs(id,:) );
    
    IT = IT+1;
    if IT > MAX_ITERATIONS
      break;
    end
    if etime( clock , START_TIME ) > MAX_TIME
      break;
    end

    onBoundary = false;
    for c = 1:N
      xxid  = find( xx{c} == X(c) , 1 ,'first' );

      if numel(xx{c}) > 1  && numel(xxid) &&  ( xxid == 1 || xxid == numel(xx{c}) ) 
        onBoundary = true;
        break;
      end
    end
    
    if onBoundary
      R = R / shrinkFactor;
    else
      R = R * shrinkFactor;
    end
    
    if VERBOSE
      fprintf('R: %30.20g   E: %30.20g   X: ',R,E);
      fprintf('%g ',X(:));
      fprintf('\n');
    end
    if ~isempty( FCN )
      try
        feval( FCN , Rsz(X) , E );
      end
      drawnow('expose');
    end
    
    if ~rem( IT , 1000 ), R = max( R , R_init ); end
    
  end
  
end
