function S = meshSlice( M , P , OUTmode )

  if nargin < 3, OUTmode = []; end


  if isnumeric( M )
    M = struct( 'xyz' , M , 'tri' , [ ( 1:size(M,1)-1 ).' , ( 2:size(M,1) ).'  ] ,'celltype',3 );
  end
  
  M.celltype = meshCelltype( M );

  if isempty( OUTmode )
    switch M.celltype
      case 3, OUTmode = 3;
      case 5, OUTmode = 3;
      case 10, OUTmode = 1;
    end
  end

  
  
  
  S = struct('xyz',[],'tri',[]);
  switch M.celltype
    case 3,  S.celltype = 1;
    case 5,  S.celltype = 3;
    case 10, S.celltype = 5;
  end  
  if size( M.tri ,1) > 1
    M.xyz(:,end+1:3) = 0;
    
    [~,iZ] = getPlane( P );
    Z = M.xyz * iZ( 3 , 1:size( M.xyz ,2) ).' + iZ(3,4);

    if ~all( all( Z( M.tri ) > 0 ,2) | all( Z( M.tri ) < 0 ,2) )
      S = MeshZeroContour( M , Z );
    end
  end

  
  
  if ischar( OUTmode )
    switch lower( OUTmode )
      case {'raw'},         OUTmode = 0;
      case {'mesh'},        OUTmode = 1;
      case {'tidy'},        OUTmode = 2;
      case {'line'},        OUTmode = 3;
      case {'cell'},        OUTmode = 4;
      otherwise, error('unknown OUTmode');
    end
  end
  
  switch OUTmode
    case 0
      %return as it is
    case 1
      %sort the elements by their connectivity (the largest piece first)
      S = MeshTidy( S , 0 );
      conn = meshFacesConnectivity( S );
      nconn = accumarray( conn , 1 );
      [~,nconn] = sort( nconn , 'descend' );
      conn = nconn( conn );
      [~,order] = sort( conn );
      for f = fieldnames( S ).', if ~strncmp( f{1} , 'tri' , 3 ), continue; end
        S.(f{1}) = S.(f{1})( order ,:);
      end
    case 2
      %mesh as a in segment order
      S      = MeshTidy( S , 0 );
      [~,S]  = Mesh2Polyline( S );
    case 3
      if S.celltype == 3
        S = MeshTidy( S , 0 );
        try,    S = Mesh2Polyline( S );
        catch,  S = mesh2contours( S );
        end
        
      elseif S.celltype == 1
        S = S.xyz;
      end
    case 4
      S = MeshTidy( S , 0 );
      C = Mesh2Polyline( S.tri );
      for c = 1:numel(C)
        C{c} = S.xyz( C{c}(1,:) ,:);
      end
      S = C;
    otherwise, error('incorrect OUTmode');
  end
  
  return;

  

  M.tri( all( reshape( M.xyz( M.tri , 3 ) , [] , size(M.tri,2) ) > 0 , 2 ),: ) = [];
  M.tri( all( reshape( M.xyz( M.tri , 3 ) , [] , size(M.tri,2) ) < 0 , 2 ),: ) = [];

  if numel( M.tri ) < 1, 
    return; end

  M = TidyMesh( M , -1 );

  if size(M.tri,2) == 2
    
    w = M.xyz( M.tri(:,2) ,3) ./ ( M.xyz( M.tri(:,2) ,3) - M.xyz( M.tri(:,1) ,3) );
    S = bsxfun(@times,M.xyz( M.tri(:,1) ,:),w) + bsxfun(@times,M.xyz( M.tri(:,2) ,:),1-w);

    S = bsxfun( @plus , S * Z(1:3,1:3).' , Z(1:3,4).' );
    
    return;
  end
  
  
  
  
  M = vtkClipPolyData( M , [0 0 0;0 0 1] , 'GenerateClipScalarsOff' , [] );
  M = vtkFeatureEdges( M , 'BoundaryEdgesOn',[],'FeatureEdgesOff',[],'NonManifoldEdgesOff',[],'ManifoldEdgesOff',[],'ColoringOff',[]);

  if isempty(M) || isempty(fieldnames(M))
    return;
  end

  
  
  
  M.tri( M.tri(:,3) == 0 , 3 ) = M.tri( M.tri(:,3) == 0 , 1 );

  M.tri( any( abs( reshape( M.xyz( M.tri , 3 ) , [] , 3 ) ) > 1e-6 , 2 ) , : ) = [];
  M.xyz(:,3) = 0;

  M.xyz = bsxfun( @plus , M.xyz * Z(1:3,1:3).' , Z(1:3,4).' );
  
  if FAST
    S = M.xyz; conn = M.tri;
    return;
  end
  
  
  M = vtkCleanPolyData( M , 'SetAbsoluteTolerance',1e-12,'ToleranceIsAbsoluteOn',[],'ConvertLinesToPointsOn',[],'PointMergingOn',[]);
  if ~isfield( M , 'tri' )
    return;
  end 
  
  M.tri( sum( M.tri == 0 , 2 ) > 1 ,: ) = [];

  if ~isfield(M,'tri') || isempty(M.tri)
    return;
  end
  
   M.tri(:,3) = [];

  
  
  while any( M.tri(:) )
    
    node_head = M.tri( find( M.tri , 1 ) );
    
    node_head  = getHead( M.tri , node_head );
    chain = node_head;
    
    isClosed = false;
    while true
      node_row = find( any( M.tri == chain(end) , 2 ) , 1 , 'first' );
      if isempty( node_row ), break; end
      
      node_new = setdiff( M.tri( node_row , :) , chain(end) );
      M.tri( node_row , : ) = 0;
      
      chain = [ chain , node_new ];

      if node_head == node_new
        isClosed = true;
        break; 
      end
    end
    
    conn = [ conn ; size(S,1) + [ ( 1:numel(chain)-1 ).' ,  ( 2:numel(chain) ).'  ] ];
    if isClosed
      conn(end,2) = size(S,1)+1;
    elseif allClosed
      conn = [ conn ; conn(end,2) size(S,1)+1 ];
    end
      

    S = [ S(:) ; chain(:) ; NaN ];
    
    
  end

  
  S(end) = [];   %remove the last NaN
  
  M.xyz = [ M.xyz ; NaN NaN NaN ];
  S( isnan( S ) ) = size( M.xyz , 1 );
  S = M.xyz( S , : );

  

  
  function h = getHead( G , h0 )
    
    h = h0;
    while true
      h_idx = find( any( G == h , 2 ) , 1 , 'first' );
      if isempty( h_idx ), break; end
      
      h = setdiff( G( h_idx , :) , h );
      
      if h == h0, break; end

      G( h_idx , : ) = 0;
    end
    
  end






%   M = vtkPolyDataConnectivityFilter( M , 'ColorRegionsOn',[],'SetExtractionModeToAllRegions',[]);
% 
%   M.tri( M.tri(:,3) == 0 , 3 ) = M.tri( M.tri(:,3) == 0 , 1 );
% 
%   
%   M.xyzRegionId = M.xyzRegionId + 1;
%   n_vert = accumarray( M.xyzRegionId , 1 );
%   [n_vert , id_vert ] = sort( n_vert , 'descend' );
%   
%   B = zeros(size(M.xyz,1)+numel(id_vert),3);
%   bb = 1;
%   for id = id_vert(:)'
%     MM = DeletePoints( M , find( M.xyzRegionId ~= id ) );
%     MM.tri( MM.tri(:,1) == MM.tri(:,3) | MM.tri(:,1) == MM.tri(:,2) , 3 ) = 0;
%     if ~any( MM.tri(:,3) )
%       MM.tri(:,3) = [];
%     else
%       continue;
%     end
% 
%     MMM = MM.tri;
%     if isempty(MMM), continue; end
%     
%     row = 1;
%     n = MMM(row,2);
%     MMM(row,:)= 0;
%     while any( MMM(:) )
%       row = find( MMM(:,1) == n , 1 , 'last' );
%       if ~isempty( row )
%         n = MMM(row,2);
%         MMM(row,:)= 0;
%         continue;
%       end
%       row = find( MMM(:,2) == n , 1 , 'last' );
%       if ~isempty( row )
%         n = MMM(row,1);
%         MMM(row,:)= 0;
%         continue;
%       end
%       break;
%     end
%     
%     
%     while any( MM.tri(:) );
%       row = find( MM.tri(:,1) == n , 1 , 'last' );
%       if ~isempty( row )
%         B(bb,1:3) = MM.xyz( MM.tri(row,1) , : );  bb = bb+1;
%         n = MM.tri(row,2);
%         MM.tri(row,:) = 0;
%         continue;
%       end
% 
%       row = find( MM.tri(:,2) == n , 1 , 'last' );
%       if ~isempty( row )
%         B(bb,1:3) = MM.xyz( MM.tri(row,2) , : );  bb = bb+1;
%         n = MM.tri(row,1);
%         MM.tri(row,:) = 0;
%         continue;
%       end
%       break;
%     end
%     B(bb,1:3) = NaN; bb = bb+1;
%       
%   end
%   bb = bb-2;
%   B = B(1:bb,:);
  
end
