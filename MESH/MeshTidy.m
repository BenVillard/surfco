function [M,pID] = MeshTidy( M , th , REMOVE_COLLAPSED , METRIC )

  %first remove nans
  w = find( any( isnan( M.xyz ) ,2 ) );
  M.tri( any( ismember( M.tri , w ) ,2) ,:) = [];
  M.xyz( w ,:) = -1e+308;

  dots = repmat( {':'} , [1,20]);

  if nargin < 2 || isempty( th )
    EL = [];
    for i = 1:size(M.tri,2)
      for j = i+1:size(M.tri,2)
        EL = [ EL ; M.tri(:,[i,j]) ];
      end
    end
    EL = sort( EL , 2 );
    EL = unique( EL ,'rows' );
    EL = M.xyz( EL(:,1) , : ) - M.xyz( EL(:,2) , : );
    EL = sum( EL.^2 , 2 );
    
    EL( ~EL ) = [];
    th = sqrt( min( EL ) )/2;
  end
  if nargin < 3 || isempty( REMOVE_COLLAPSED ), REMOVE_COLLAPSED = true;  end
  if nargin < 4 || isempty( METRIC ),           METRIC = 1;               end
  METRIC = sqrt(METRIC(:).');
  
  
  Fs = fieldnames( M );
  Fxyz = Fs( strncmp( Fs , 'xyz' , 3 ) ); Fxyz = Fxyz(:).';
  Ftri = Fs( strncmp( Fs , 'tri' , 3 ) ); Ftri = Ftri(:).';
  
  
  
  usedNODES = unique( M.tri(:) , 'sorted' );
  newIDX = zeros( max( usedNODES ) , 1 );
  newIDX( usedNODES ) = 1:numel(usedNODES);
  
  M.tri = reshape( newIDX( M.tri ) , size( M.tri ) );
  for f = Fxyz
    M.(f{1}) = M.(f{1})(usedNODES,dots{:});
  end

  
  
  if th < 0
    %remove the collapsed faces
    if REMOVE_COLLAPSED
      properFaces = all( diff( sort( M.tri , 2 ) , 1 , 2 ) , 2 );
      for f = Ftri
        M.(f{1}) = M.(f{1})( properFaces ,dots{:});
      end
    end    
    return;
  end

  X = [];
  for f = Fxyz
    X = [ X , M.(f{1}) ];
  end
  X = bsxfun(@times, X , METRIC(1:min(end,size(X,2))) );
  X( : , all( isnan( X ) , 1 ) ) = [];

  [~,usedNODES,newIDX] = unique( M.xyz , 'rows' , 'stable' );
  if numel( usedNODES ) ~= size( M.xyz , 1 )
    usedNODES = sort( usedNODES );
    M.tri = newIDX( M.tri );
    for f = Fxyz
      M.(f{1}) = M.(f{1})(usedNODES,dots{:});
    end
    if th > 0, X = X( usedNODES ,:); end
  end  
  
  if th == 0
    %remove the collapsed faces
    if REMOVE_COLLAPSED
      properFaces = all( diff( sort( M.tri , 2 ) , 1 , 2 ) , 2 );
      for f = Ftri
        M.(f{1}) = M.(f{1})( properFaces ,dots{:});
      end
    end    
    return;
  end
  
  
  
  %inter-nodes distance
  D = ipd( X , [] );
  D( D > th ) = Inf;

  %new indexes... on each "cluster" used the node with lower index
  newIDX = double( isfinite( D ) );
  newIDX( ~~newIDX ) = 1:sum( newIDX(:) );
  newIDX(  ~newIDX ) = Inf;
  [~,newIDX] = min( newIDX , [] , 2 );
  
  %remap faces
  M.tri = newIDX( M.tri );
  
  %remove the fullCollapsed faces
  if true
    notCollapsedFaces = ~~var( M.tri , [] , 2 );
    for f = Ftri
      M.(f{1}) = M.(f{1})( notCollapsedFaces ,dots{:});
    end
  end
  
  %remove the collapsed faces
  if REMOVE_COLLAPSED
    properFaces = all( diff( sort( M.tri , 2 ) , 1 , 2 ) , 2 );
    for f = Ftri
      M.(f{1}) = M.(f{1})( properFaces ,dots{:});
    end
  end
  
  usedNODES = unique( M.tri , 'sorted' );
  %average the coordinate and attributes of the "clusters"
  for n = usedNODES(:).'
    w = ~isinf( D(n,:) );
    for f = Fxyz
      M.(f{1})(n,dots{:}) = mean( M.(f{1})(w,dots{:}) , 1 );
    end
  end
  
  %remove unused nodes
  for f = Fxyz
    M.(f{1}) = M.(f{1})(usedNODES,dots{:});
  end
  
  %remap the faces after unused nodes were removed.
  newIDX( usedNODES ) = 1:numel( usedNODES );
  M.tri = newIDX( M.tri );
  
  if nargout > 1
    nX = [];
    for f = Fxyz
      nX = [ nX , M.(f{1}) ];
    end
    nX = bsxfun(@times, nX , METRIC(1:min(end,size(X,2))) );
    nX( : , all( isnan( nX ) , 1 ) ) = [];
    
    %D = bsxfun( @minus , permute( nX , [1 3 2] ) , permute( X , [3 1 2] ) );
    %D = D.^2;
    %D = sum( D , 3 );
    D = ipd( nX , X );

    [~,pID] = min( D , [] , 2 );
  end
  
  
  
if 0
%% find singular non-manifoldness
ED = meshEdges( M.tri );
A = sparse( ED(:,1) , ED(:,2) , true , size( M.xyz,1) , size( M.xyz,1) );
A = A+A.';
arrayfun( @(i)numel( SortChain( find( A(:,i) ) , ED ) ) , 1:size(A,2) ) .* full(~~sum( A ,1))
  
end
  
  
  
  
  
end
