function [B,conn] = SliceMesh( M , plane , allClosed , FAST )

%   if numel(plane) == 6 && size(plane,1) == 3 && size(plane,2) == 3
%     
%   end
  
  if nargin < 3, allClosed = false; end
  if nargin < 4, FAST = false; end

  B = []; conn = [];

  if isnumeric( M )
    M = struct( 'xyz' , M , 'tri' , [ 1:size(M,1)-1 ; 2:size(M,1) ].' );
    b = find( any( ~isfinite( M.xyz ) ,2) );
    M.tri( any( ismember( M.tri , b ) ,2),:) = [];
  else
    M = MakeMesh( M );
  end
  

  if numel( M.tri ) < 1, 
    return; end

  [Z,iZ] = getPlane( plane );
  M.xyz = bsxfun( @plus , M.xyz * iZ(1:3,1:3).' , iZ(1:3,4).' );

  M.tri( all( reshape( M.xyz( M.tri , 3 ) , [] , size(M.tri,2) ) > 0 , 2 ),: ) = [];
  M.tri( all( reshape( M.xyz( M.tri , 3 ) , [] , size(M.tri,2) ) < 0 , 2 ),: ) = [];

  if numel( M.tri ) < 1, 
    return; end

  M = TidyMesh( M , -1 );

  if size(M.tri,2) == 2
    
    w = M.xyz( M.tri(:,2) ,3) ./ ( M.xyz( M.tri(:,2) ,3) - M.xyz( M.tri(:,1) ,3) );
    B = bsxfun(@times,M.xyz( M.tri(:,1) ,:),w) + bsxfun(@times,M.xyz( M.tri(:,2) ,:),1-w);

    B = bsxfun( @plus , B * Z(1:3,1:3).' , Z(1:3,4).' );
    
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
    B = M.xyz; conn = M.tri;
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
    
    conn = [ conn ; size(B,1) + [ ( 1:numel(chain)-1 ).' ,  ( 2:numel(chain) ).'  ] ];
    if isClosed
      conn(end,2) = size(B,1)+1;
    elseif allClosed
      conn = [ conn ; conn(end,2) size(B,1)+1 ];
    end
      

    B = [ B(:) ; chain(:) ; NaN ];
    
    
  end

  
  B(end) = [];   %remove the last NaN
  
  M.xyz = [ M.xyz ; NaN NaN NaN ];
  B( isnan( B ) ) = size( M.xyz , 1 );
  B = M.xyz( B , : );

  

  
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
