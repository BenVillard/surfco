function c = meshCelltype( M )

  % EMPTY_CELL                       = 0,
  % VERTEX                           = 1,
  % POLY_VERTEX                      = 2,
  % LINE                             = 3,
  % POLY_LINE                        = 4,
  % TRIANGLE                         = 5,
  % TRIANGLE_STRIP                   = 6,
  % POLYGON                          = 7,
  % PIXEL                            = 8,
  % QUAD                             = 9,
  % TETRA                            = 10,
  
  nsd = 3;  %default NumberOfSpatialDimensions
  if isfield( M , 'xyz' ), nsd = size( M.xyz ,2); end
  nnf = []; %default NumberOfNodesPerFace
  if isfield( M , 'tri' ), nnf = size( M.tri ,2); end
  
  if isfield( M , 'celltype' )
    c = M.celltype;
    if ~isscalar(c) 
      if numel(c) ~= size( M.tri ,1)
        error('incorrect number of celltypes specified.');
      end
      if all( c == c(1) )
        c = c(1);
      end
    end
    if isscalar( c )
      switch c
        case 0,  error('celltype = 0?');
        case 1,  if ~isempty(nnf) &&  nnf ~= 1, error('celltype is 1 (POINT), then NumberOfNodesPerFace should be 1.'); end
        case 2,  error('celltype = 2?');
        case 3,  if ~isempty(nnf) &&  nnf ~= 2, error('celltype is 3 (LINE), then NumberOfNodesPerFace should be 2.'); end
        case 4,  error('celltype = 4?');
        case 5,  if ~isempty(nnf) &&  nnf ~= 3, error('celltype is 5 (TRIANGLE), then NumberOfNodesPerFace should be 3.'); end
        case 6,  error('celltype = 6?');
        case 7,  error('celltype = 7?');
        case 8,  error('celltype = 8?');
        case 9,  error('celltype = 9?');
        case 10, if ~isempty(nnf) &&  nnf ~= 4, error('celltype is 10 (TETRA), then NumberOfNodesPerFace should be 4.'); end
      end
    end
  else
    Mtype = [ nnf , nsd ];
    if     0
    elseif isequal( Mtype , [2 3] ) || isequal( Mtype , [2 2] )
      c = 3;
    elseif isequal( Mtype , [3 3] ) || isequal( Mtype , [3 2] )
      c = 5;
    elseif isequal( Mtype , [4 2] )
      c = 9;
    elseif isequal( Mtype , [4 3] )
      c = 10;
    else
      c = 0;
      %error('unknown type of mesh');
    end
  end

end
