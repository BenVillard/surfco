function M = MakeMesh( varargin )
% 
% mesh= MakeMesh( 'Points', xyz
%                 'Triangles', tri
%
   
  M = struct();
  
  if nargin == 1 && isempty( varargin{1} )
    M = struct('xyz',[],'tri',[]);
    return;
  end

  if numel( varargin ) == 2 && isnumeric( varargin{1} ) && ischar( varargin{2} ) && strcmpi( varargin{2} , 'contour' )
    try, M.xyz = double( varargin{1} ); end
    try, M.tri = [ 1:size(M.xyz,1)-1 ; 2:size(M.xyz,1) ].'; end
    M = MeshRemoveNodes( M , any( isnan( M.xyz ) , 2) );
    
    return;
  end
  
  
  if isa( varargin{1} , 'delaunayTriangulation' )
    M.xyz = varargin{1}.Points;
    M.tri = varargin{1}.ConnectivityList;
    return;
  end
    
  
  if numel( varargin ) == 1 && numel( varargin{1} ) == 1 && ishandle( varargin{1} )
    M.xyz = get( varargin{1} , 'Vertices' );
    M.tri = get( varargin{1} , 'Faces'    );
    return;
  end

  
  if numel( varargin ) == 1 && isnumeric( varargin{1} )
    M.xyz = double( varargin{1} );
    return;
  end
  
  if numel( varargin ) == 1 && isstruct( varargin{1} )
    if isfield( varargin{1} , 'vertices' ) && ~isfield( varargin{1} , 'xyz' )
      varargin{1}.xyz = varargin{1}.vertices;
    end
    if isfield( varargin{1} , 'faces' ) && ~isfield( varargin{1} , 'tri' )
      varargin{1}.tri = varargin{1}.faces;
    end
  end
  
  
  if numel( varargin ) == 1 && isstruct( varargin{1} )
    try, M.xyz = double( varargin{1}.xyz ); end
    try, M.tri = double( varargin{1}.tri ); end
    return;
  end

  if numel( varargin ) == 2 && isnumeric( varargin{1} ) && isnumeric( varargin{2} )
    try, M.xyz = double( varargin{1} ); end
    try, M.tri = double( varargin{2} ); end
    return;
  end
  
  if numel( varargin ) == 2 && isnumeric( varargin{1} ) && isstruct( varargin{2} )
    try, M.xyz = double( varargin{1} ); end
    try, M.tri = double( varargin{2}.tri ); end
    return;
  end
  
  if numel( varargin ) == 2 && isstruct( varargin{1} ) && isnumeric( varargin{2} )
    try, M.xyz = double( varargin{1}.xyz ); end
    try, M.tri = double( varargin{2} ); end
    return;
  end

  if numel( varargin ) == 2 && isstruct( varargin{1} ) && isstruct( varargin{2} )
    try, M.xyz = double( varargin{1}.xyz ); end
    try, M.tri = double( varargin{2}.tri ); end
    return;
  end
  
% zz= parseargs( varargin,'points','p');
% if zz
%   M.xyz= varargin{zz+1};
% end
% if size( M.xyz,2 ) < 3
%   M.xyz(:,3)=0;
% end
% 
% zz= parseargs( varargin,'triangles','tri','t','triangle');
% if zz
%   M.tri= varargin{zz+1};
% else
%   M.tri= delaunayn( M.xyz(:,[1 2 3]) );
% end

end