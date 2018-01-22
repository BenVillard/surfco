function [M,pID] = TidyMesh( M , th , REMOVE_COLLAPSED , METRIC )

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
  if nargin < 3 || isempty( REMOVE_COLLAPSED ), REMOVE_COLLAPSED = false; end
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
    M.(f{1}) = M.(f{1})(usedNODES,:);
  end

  
  
  if th < 0, return; end

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
      M.(f{1}) = M.(f{1})(usedNODES,:);
    end
    if th > 0, X = X( usedNODES ,:); end
  end  
  
  
  %remove the collapsed faces
  if REMOVE_COLLAPSED
    properFaces = all( diff( sort( M.tri , 2 ) , 1 , 2 ) , 2 );
    for f = Ftri
      M.(f{1}) = M.(f{1})( properFaces ,:);
    end
  end
  
  
  if th == 0, return; end
  
  
  
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
      M.(f{1}) = M.(f{1})( notCollapsedFaces ,:);
    end
  end
  
  %remove the collapsed faces
  if REMOVE_COLLAPSED
    properFaces = all( diff( sort( M.tri , 2 ) , 1 , 2 ) , 2 );
    for f = Ftri
      M.(f{1}) = M.(f{1})( properFaces ,:);
    end
  end
  
  usedNODES = unique( M.tri , 'sorted' );
  %average the coordinate and attributes of the "clusters"
  for n = usedNODES(:).'
    w = ~isinf( D(n,:) );
    for f = Fxyz
      M.(f{1})(n,:) = mean( M.(f{1})(w,:) , 1 );
    end
  end
  
  %remove unused nodes
  for f = Fxyz
    M.(f{1}) = M.(f{1})(usedNODES,:);
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
end
