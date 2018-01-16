function M = MeshSubdivide( M , W )

if 0
  M.xyz = randn(2000,3); M.tri = delaunayn( M.xyz ); M.celltype = 10; M.triL = (1:size(M.tri,1)).';
  MM = MeshSubdivide( M , 1:11:size(M.tri,1) );
%   plMESH( MM )
  w = MeshQuality( MM , 'volume' ) < 0; unique( MM.triCASE(w) )

  %%
end

if 0
  
    p1=1;p2=2;p3=3;p4=4;p12=5;p13=6;p14=7;p23=8;p24=9;p34=10;
    T=[];F=[];w=[];

    
    
    
    
    V = struct('tri',T,'xyz',[0,0,0;2,0,0;1,1.6,0;1,0.6,1.4;1,0,0;0.5,0.8,0;0.5,0.3,0.7;1.5,0.8,0;1.5,0.3,0.7;1,1.1,0.7],'celltype',10);
    MeshQuality( V , 'volume' )
    plMESH(V);
  
  %%  
end
if 0
%   M.xyz = rand(20,2); M.tri = delaunayn( M.xyz ); M.triL = (1:size(M.tri,1)).';
  W = [ 1 8 10 20 25 ]; W( W>size(M.tri,1) ) = []; cm = rand(size(M.tri,1),3)/2+0.5;
  subplot(121); plotMESH(  M , 'td','L' ); colormap(cm); colorbar
  MM = MeshSubdivide( M , W );
  subplot(122); plotMESH( MM , 'td','L' ); colormap(cm); colorbar
  hplot3d( getv( FacesCenter( M ) , W , ':' ) ,'*r' )
  hplotMESH( MeshBoundary(MM) , '-','EdgeColor','r','LineWidth',3)
  
  %%
end


  nT  = size( M.tri , 1);   %number of faces
  if nargin < 2, W = 1:nT; end
  if islogical( W )
    if numel(W) ~= nT, error('a number of triangles logical were expected'); end
    W = find(W);
  else
    W = unique( W(:) ,'sorted');
  end
  if isempty( W ), return; end
  if any( mod(W,1) ), error('indexes must be integers or logicals.'); end
  if any( W < 1 ),    error('indexes must be positive integers.');    end
  if any( W > nT ),   error('Index exceeds number of faces.');        end

  T  = M.tri;
  F   = ( 1:nT ).';         %face indexes
  nP  = size( M.xyz , 1);   %number of points
  

  if 0
  elseif size( T ,2) == 2   %segments case

    %middle points on edges
    E   = T( W ,:);
    for f = fieldnames( M ).', if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
      M.(f{1}) = [ M.(f{1}) ; ( M.(f{1})( E(:,1) ,:) + M.(f{1})( E(:,2) ,:) )/2 ];
    end
    P = nP + ( 1:size(E,1) ); P = P(:);    %indexes of the new points


    %old and new faces
    T = [   T ; ... 
          [ T( W ,1) , P        ] ;...
          [ P        , T( W ,2) ] ;...
        ];
    F = [ F ; W ; W ];      %original ids of the new faces
    
  elseif size( T ,2) == 3          %% triangular mesh

    allE = sort( [ T(:,[1 2]) ; T(:,[2 3]) ; T(:,[1 3]) ] ,2);

    if 0                                %% case 1, in case of cuadrilateral faces, split the whole triengle into 4

      while 1
        E = allE( [ W ; W + nT ; W + 2*nT ] , : );
        E = unique( E , 'rows' );
        TtD = reshape( ismember( allE , E , 'rows' ) , nT , 3 );
        Wp = W;
        W = find( sum( TtD ,2) > 1 );
        if isequal( Wp , W ), break; end
      end
      
    else                         %% case 2, try to keep to minimum the number of faces
      
      E = allE( [ W ; W + nT ; W + 2*nT ] , : );
      E = unique( E , 'rows' );
      
    end


    %triangles containing edges
    [ ~ , ET(:,1) ] = ismember( sort( T(:,[1 2]) ,2) , E , 'rows' );
    [ ~ , ET(:,2) ] = ismember( sort( T(:,[2 3]) ,2) , E , 'rows' );
    [ ~ , ET(:,3) ] = ismember( sort( T(:,[1 3]) ,2) , E , 'rows' );
    %original faces to be removed
    W = find( any( ~~ET ,2) );


    %middle points on edges
    for f = fieldnames( M ).', if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
      M.(f{1}) = [ M.(f{1}) ; ( M.(f{1})( E(:,1) ,:) + M.(f{1})( E(:,2) ,:) )/2 ];
    end
    P = nP + ( 1:size(E,1) ); P = P(:);    %indexes of the new points


    %first, triangles to be divided into 4.
    w = find( ~~ET(:,1) & ~~ET(:,2) & ~~ET(:,3) );
    T = [ T ; ...
          [     T(w, 1 )   , P( ET(w, 1 ) ) , P( ET(w, 3 ) ) ] ;...
          [ P( ET(w, 1 ) ) ,     T(w, 2 )   , P( ET(w, 2 ) ) ] ;...
          [ P( ET(w, 3 ) ) , P( ET(w, 2 ) ) ,     T(w, 3 )   ] ;...
          [ P( ET(w, 1 ) ) , P( ET(w, 2 ) ) , P( ET(w, 3 ) ) ] ;...
        ];
    F = [ F ; w ; w ; w ; w ];


    %triangles to be divided at edge 1-2
    w = find( ~~ET(:,1) & ~ET(:,2) & ~ET(:,3) );
    T = [ T ; ...
          [     T(w, 1 )   , P( ET(w, 1 ) ) ,     T(w, 3 )   ] ;...
          [ P( ET(w, 1 ) ) ,     T(w, 2 )   ,     T(w, 3 )   ] ;...
        ];
    F = [ F ; w ; w ];


    %triangles to be divided at edge 2-3
    w = find( ~ET(:,1) & ~~ET(:,2) & ~ET(:,3) );
    T = [ T ; ...
          [     T(w, 1 )   ,     T(w, 2 )   , P( ET(w, 2 ) ) ] ;...
          [     T(w, 1 )   , P( ET(w, 2 ) ) ,     T(w, 3 )   ] ;...
        ];
    F = [ F ; w ; w ];


    %triangles to be divided at edge 1-3
    w = find( ~ET(:,1) & ~ET(:,2) & ~~ET(:,3) );
    T = [ T ; ...
          [     T(w, 1 )   ,     T(w, 2 )   , P( ET(w, 3 ) ) ] ;...
          [     T(w, 2 )   ,     T(w, 3 )   , P( ET(w, 3 ) ) ] ;...
        ];
    F = [ F ; w ; w ];

    
    %triangles to be divided at edge 1-2 && 2-3
    w = find( ~~ET(:,1) & ~~ET(:,2) & ~ET(:,3) );
    T = [ T ; ...
          [ P( ET(w, 1 ) ) ,     T(w, 2 )   , P( ET(w, 2 ) ) ] ;...
          [     T(w, 1 )   , P( ET(w, 1 ) ) ,     T(w, 3 )   ] ;...
          [ P( ET(w, 1 ) ) , P( ET(w, 2 ) ) ,     T(w, 3 )   ] ;...
        ];
    F = [ F ; w ; w ; w ];

    
    %triangles to be divided at edge 1-2 && 1-3
    w = find( ~~ET(:,1) & ~ET(:,2) & ~~ET(:,3) );
    T = [ T ; ...
          [     T(w, 1 )   , P( ET(w, 1 ) ) , P( ET(w, 3 ) ) ] ;...
          [ P( ET(w, 1 ) ) ,     T(w, 3 )   , P( ET(w, 3 ) ) ] ;...
          [ P( ET(w, 1 ) ) ,     T(w, 2 )   ,     T(w, 3 )   ] ;...
        ];
    F = [ F ; w ; w ; w ];
    

    %triangles to be divided at edge 2-3 && 1-3
    w = find( ~ET(:,1) & ~~ET(:,2) & ~~ET(:,3) );
    T = [ T ; ...
          [ P( ET(w, 3 ) ) , P( ET(w, 2 ) ) ,     T(w, 3 )   ] ;...
          [     T(w, 1 )   , P( ET(w, 2 ) ) , P( ET(w, 3 ) ) ] ;...
          [     T(w, 1 )   ,     T(w, 2 )   , P( ET(w, 2 ) ) ] ;...
        ];
    F = [ F ; w ; w ; w ];
    
    
  elseif size( T ,2) == 4  &&  isfield( M , 'celltype' )  &&   M.celltype == 10          %% tetrahedral mesh

    allE = sort( [ T(:,[1 2]) ; T(:,[1 3]) ; T(:,[1 4]) ; T(:,[2 3]) ; T(:,[2 4]) ; T(:,[3 4]) ] ,2);

    while 1
      E = allE( [ W ; W + nT ; W + 2*nT ; W + 3*nT ; W + 4*nT ; W + 5*nT ] , : );
      E = unique( E , 'rows' );
      ET = reshape( ismember( allE , E , 'rows' ) , nT , 6 );
      Wp = W;
      
      W = false;
      W = W | sum( ET ,2) > 3;
      W = W | all( bsxfun(@eq, ET , [0 0 1 0 1 1] ) ,2);
      W = W | all( bsxfun(@eq, ET , [0 0 1 1 0 1] ) ,2);
      W = W | all( bsxfun(@eq, ET , [0 0 1 1 1 0] ) ,2);
      W = W | all( bsxfun(@eq, ET , [0 1 0 0 1 1] ) ,2);
      W = W | all( bsxfun(@eq, ET , [0 1 0 1 0 1] ) ,2);
      W = W | all( bsxfun(@eq, ET , [0 1 0 1 1 0] ) ,2);
      W = W | all( bsxfun(@eq, ET , [0 1 1 0 1 0] ) ,2);
      W = W | all( bsxfun(@eq, ET , [0 1 1 1 0 0] ) ,2);
      W = W | all( bsxfun(@eq, ET , [1 0 0 0 1 1] ) ,2);
      W = W | all( bsxfun(@eq, ET , [1 0 0 1 0 1] ) ,2);
      W = W | all( bsxfun(@eq, ET , [1 0 0 1 1 0] ) ,2);
      W = W | all( bsxfun(@eq, ET , [1 0 1 0 0 1] ) ,2);
      W = W | all( bsxfun(@eq, ET , [1 0 1 1 0 0] ) ,2);
      W = W | all( bsxfun(@eq, ET , [1 1 0 0 0 1] ) ,2);
      W = W | all( bsxfun(@eq, ET , [1 1 0 0 1 0] ) ,2);
      W = W | all( bsxfun(@eq, ET , [1 1 1 0 0 0] ) ,2);

      W = find( W );
      if isequal( Wp , W ), break; end
    end
      
    ET = zeros( nT , 6 );
    %tetras containing edges
    [ ~ , ET(:,1) ] = ismember( sort( T(:,[1 2]) ,2) , E , 'rows' );
    [ ~ , ET(:,2) ] = ismember( sort( T(:,[1 3]) ,2) , E , 'rows' );
    [ ~ , ET(:,3) ] = ismember( sort( T(:,[1 4]) ,2) , E , 'rows' );
    [ ~ , ET(:,4) ] = ismember( sort( T(:,[2 3]) ,2) , E , 'rows' );
    [ ~ , ET(:,5) ] = ismember( sort( T(:,[2 4]) ,2) , E , 'rows' );
    [ ~ , ET(:,6) ] = ismember( sort( T(:,[3 4]) ,2) , E , 'rows' );
    %original faces to be removed
    W = find( any( ~~ET ,2) ); %size(W)
   
    
    %middle points on edges
    for f = fieldnames( M ).', if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
      M.(f{1}) = [ M.(f{1}) ; ( M.(f{1})( E(:,1) ,:) + M.(f{1})( E(:,2) ,:) )/2 ];
    end
    P = nP + ( 1:size(E,1) ); P = P(:);    %indexes of the new points


    %C = zeros( nT , 1);
    
    %first, tetras to be divided into 8.
    w = find( all(  bsxfun(@eq, ~~ET , ~~[1 1 1 1 1 1] ) ,2) );
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p12  ,  p13  ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 8.1 ];
    T = [ T ;    p12  ,  p2   ,  p23  ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 8.2 ];
    T = [ T ;    p13  ,  p23  ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 8.3 ];
    T = [ T ;    p14  ,  p24  ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 8.4 ];
    T = [ T ;    p12  ,  p13  ,  p14  ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 8.5 ];
    T = [ T ;    p12  ,  p13  ,  p24  ,  p23  ]; F = [ F ; w ]; %C = [ C ; w*0 + 8.6 ];
    T = [ T ;    p13  ,  p14  ,  p24  ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 8.7 ];
    T = [ T ;    p13  ,  p23  ,  p34  ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 8.8 ];
    end
  
    %tetras to be divided at edge 1-2
    w = find( all(  bsxfun(@eq, ~~ET , ~~[1 0 0 0 0 0] ) ,2) );
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p12  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 12.1 ];
    T = [ T ;    p12  ,  p2   ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 12.2 ];
    end
    
    %tetras to be divided at edge 1-3
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 1 0 0 0 0] ) ,2) );
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p2   ,  p13  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 13.1 ];
    T = [ T ;    p13  ,  p2   ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 13.2 ];
    end
    
    %tetras to be divided at edge 1-4
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 0 1 0 0 0] ) ,2) );
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p2   ,  p3   ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 14.1 ];
    T = [ T ;    p14  ,  p2   ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 14.2 ];
    end
    
    %tetras to be divided at edge 2-3
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 0 0 1 0 0] ) ,2) );
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p2   ,  p23  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 23.1 ];
    T = [ T ;    p1   ,  p23  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 23.2 ];
    end
    
    %tetras to be divided at edge 2-4
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 0 0 0 1 0] ) ,2) );
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p2   ,  p3   ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 24.1 ];
    T = [ T ;    p1   ,  p24  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 24.2 ];
    end

    %tetras to be divided at edge 3-4
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 0 0 0 0 1] ) ,2) );
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p2   ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 34.1 ];
    T = [ T ;    p1   ,  p2   ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 34.2 ];
    end
    
    %tetras to be divided at edge 1-3 & 2-4
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 1 0 0 1 0] ) ,2) );
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p13  ,  p4   ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1324.1 ];
    T = [ T ;    p13  ,  p3   ,  p4   ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1324.2 ];
    T = [ T ;    p1   ,  p13  ,  p24  ,  p2   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1324.3 ];
    T = [ T ;    p13  ,  p3   ,  p24  ,  p2   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1324.4 ];
    end
    
    %tetras to be divided at edge 1-2 & 3-4
    w = find( all(  bsxfun(@eq, ~~ET , ~~[1 0 0 0 0 1] ) ,2) );
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p12  ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1234.1 ];
    T = [ T ;    p1   ,  p12  ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1234.2 ];
    T = [ T ;    p12  ,  p2   ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1234.3 ];
    T = [ T ;    p12  ,  p2   ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1234.4 ];
    end

    %tetras to be divided at edge 1-4 & 2-3
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 0 1 1 0 0] ) ,2) );
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p2   ,  p23  ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1423.1 ];
    T = [ T ;    p1   ,  p23  ,  p3   ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1423.2 ];
    T = [ T ;    p14  ,  p2   ,  p23  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1423.3 ];
    T = [ T ;    p14  ,  p23  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1423.4 ];
    end
    
    
    %tetras to be divided at edge 1-2 & 1-3 & 2-3
    w = find( all(  bsxfun(@eq, ~~ET , ~~[1 1 0 1 0 0] ) ,2) );
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p12  ,  p13  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 121323.1 ];
    T = [ T ;    p12  ,  p2   ,  p23  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 121323.2 ];
    T = [ T ;    p13  ,  p23  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 121323.3 ];
    T = [ T ;    p12  ,  p23  ,  p13  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 121323.4 ];
    end

    %tetras to be divided at edge 2-3 & 2-4 & 3-4
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 0 0 1 1 1] ) ,2) );
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p2   ,  p23  ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 232434.1 ];
    T = [ T ;    p1   ,  p24  ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 232434.2 ];
    T = [ T ;    p1   ,  p23  ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 232434.3 ];
    T = [ T ;    p1   ,  p24  ,  p23  ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 232434.4 ];
    end
    
    %tetras to be divided at edge 1-3 & 1-4 & 3-4
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 1 1 0 0 1] ) ,2) );
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p2   ,  p13  ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 131434.1 ];
    T = [ T ;    p13  ,  p2   ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 131434.2 ];
    T = [ T ;    p14  ,  p2   ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 131434.3 ];
    T = [ T ;    p14  ,  p2   ,  p13  ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 131434.4 ];
    end
    
    %tetras to be divided at edge 1-2 & 1-4 & 2-4
    w = find( all(  bsxfun(@eq, ~~ET , ~~[1 0 1 0 1 0] ) ,2) );
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p12  ,  p3   ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 121424.1 ];
    T = [ T ;    p12  ,  p2   ,  p3   ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 121424.2 ];
    T = [ T ;    p14  ,  p24  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 121424.3 ];
    T = [ T ;    p12  ,  p24  ,  p3   ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 121424.4 ];
    end
    
    
    
    %tetras to be divided at edge 2-4 & 3-4
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 0 0 0 1 1] ) ,2) ); w( T(w,2) > T(w,3) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p24  ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 2434.1 ];
    T = [ T ;    p1   ,  p2   ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2434.2 ];
    T = [ T ;    p1   ,  p34  ,  p24  ,  p2   ]; F = [ F ; w ]; %C = [ C ; w*0 + 2434.3 ];
    end
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 0 0 0 1 1] ) ,2) ); w( T(w,2) < T(w,3) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p24  ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 2434.4 ];
    T = [ T ;    p1   ,  p3   ,  p24  ,  p2   ]; F = [ F ; w ]; %C = [ C ; w*0 + 2434.5 ];
    T = [ T ;    p1   ,  p24  ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2434.6 ];
    end
    
    %tetras to be divided at edge 2-3 & 2-4
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 0 0 1 1 0] ) ,2) ); w( T(w,3) > T(w,4) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p2   ,  p23  ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2324.1 ];
    T = [ T ;    p1   ,  p3   ,  p4   ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2324.2 ];
    T = [ T ;    p1   ,  p24  ,  p23  ,  p3   ]; F = [ F ; w ]; %C = [ C ; w*0 + 2324.3 ];
    end
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 0 0 1 1 0] ) ,2) ); w( T(w,3) < T(w,4) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p2   ,  p23  ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2324.4 ];
    T = [ T ;    p1   ,  p4   ,  p24  ,  p23  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2324.5 ];
    T = [ T ;    p1   ,  p23  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 2324.6 ];
    end
    
    %tetras to be divided at edge 2-3 & 3-4
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 0 0 1 0 1] ) ,2) ); w( T(w,2) > T(w,4) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p23  ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2334.1 ];
    T = [ T ;    p1   ,  p2   ,  p23  ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2334.2 ];
    T = [ T ;    p1   ,  p34  ,  p4   ,  p2   ]; F = [ F ; w ]; %C = [ C ; w*0 + 2334.3 ];
    end
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 0 0 1 0 1] ) ,2) ); w( T(w,2) < T(w,4) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p23  ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2334.4 ];
    T = [ T ;    p1   ,  p4   ,  p23  ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2334.5 ];
    T = [ T ;    p1   ,  p23  ,  p4   ,  p2   ]; F = [ F ; w ]; %C = [ C ; w*0 + 2334.6 ];
    end

    %tetras to be divided at edge 1-2 & 1-3
    w = find( all(  bsxfun(@eq, ~~ET , ~~[1 1 0 0 0 0] ) ,2) ); w( T(w,2) > T(w,3) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p12  ,  p13  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1213.1 ];
    T = [ T ;    p2   ,  p3   ,  p13  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1213.2 ];
    T = [ T ;    p13  ,  p12  ,  p2   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1213.3 ];
    end
    w = find( all(  bsxfun(@eq, ~~ET , ~~[1 1 0 0 0 0] ) ,2) ); w( T(w,2) < T(w,3) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p12  ,  p13  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1213.4 ];
    T = [ T ;    p3   ,  p13  ,  p12  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1213.5 ];
    T = [ T ;    p12  ,  p2   ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1213.6 ];
    end
    
    %tetras to be divided at edge 1-2 & 1-4
    w = find( all(  bsxfun(@eq, ~~ET , ~~[1 0 1 0 0 0] ) ,2) ); w( T(w,2) > T(w,4) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p12  ,  p3   ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1214.1 ];
    T = [ T ;    p2   ,  p4   ,  p3   ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1214.2 ];
    T = [ T ;    p14  ,  p12  ,  p3   ,  p2   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1214.3 ];
    end
    w = find( all(  bsxfun(@eq, ~~ET , ~~[1 0 1 0 0 0] ) ,2) ); w( T(w,2) < T(w,4) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p12  ,  p3   ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1214.4 ];
    T = [ T ;    p4   ,  p14  ,  p3   ,  p12  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1214.5 ];
    T = [ T ;    p12  ,  p2   ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1214.6 ];
    end
    
    %tetras to be divided at edge 1-2 & 2-4
    w = find( all(  bsxfun(@eq, ~~ET , ~~[1 0 0 0 1 0] ) ,2) ); w( T(w,1) > T(w,4) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p12  ,  p2   ,  p3   ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1224.1 ];
    T = [ T ;    p1   ,  p12  ,  p3   ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1224.2 ];
    T = [ T ;    p24  ,  p4   ,  p3   ,  p1   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1224.3 ];
    end
    w = find( all(  bsxfun(@eq, ~~ET , ~~[1 0 0 0 1 0] ) ,2) ); w( T(w,1) < T(w,4) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p12  ,  p2   ,  p3   ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1224.4 ];
    T = [ T ;    p4   ,  p12  ,  p3   ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1224.5 ];
    T = [ T ;    p12  ,  p4   ,  p3   ,  p1   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1224.6 ];
    end
    
    %tetras to be divided at edge 1-2 & 2-3
    w = find( all(  bsxfun(@eq, ~~ET , ~~[1 0 0 1 0 0] ) ,2) ); w( T(w,1) > T(w,3) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p12  ,  p2   ,  p23  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1223.1 ];
    T = [ T ;    p1   ,  p12  ,  p23  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1223.2 ];
    T = [ T ;    p23  ,  p3   ,  p1   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1223.3 ];
    end
    w = find( all(  bsxfun(@eq, ~~ET , ~~[1 0 0 1 0 0] ) ,2) ); w( T(w,1) < T(w,3) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p12  ,  p2   ,  p23  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1223.4 ];
    T = [ T ;    p3   ,  p12  ,  p23  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1223.5 ];
    T = [ T ;    p12  ,  p3   ,  p1   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1223.6 ];
    end

    %tetras to be divided at edge 1-3 & 1-4
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 1 1 0 0 0] ) ,2) ); w( T(w,3) > T(w,4) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p2   ,  p13  ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1314.1 ];
    T = [ T ;    p3   ,  p2   ,  p14  ,  p13  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1314.2 ];
    T = [ T ;    p14  ,  p2   ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1314.3 ];
    end
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 1 1 0 0 0] ) ,2) ); w( T(w,3) < T(w,4) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p1   ,  p2   ,  p13  ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1314.4 ];
    T = [ T ;    p4   ,  p2   ,  p14  ,  p13  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1314.5 ];
    T = [ T ;    p13  ,  p2   ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1314.6 ];
    end
    
    %tetras to be divided at edge 1-4 & 3-4
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 0 1 0 0 1] ) ,2) ); w( T(w,1) > T(w,3) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p14  ,  p2   ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1434.1 ];
    T = [ T ;    p1   ,  p2   ,  p34  ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1434.2 ];
    T = [ T ;    p34  ,  p2   ,  p1   ,  p3   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1434.3 ];
    end
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 0 1 0 0 1] ) ,2) ); w( T(w,1) < T(w,3) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p14  ,  p2   ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1434.4 ];
    T = [ T ;    p3   ,  p2   ,  p34  ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1434.5 ];
    T = [ T ;    p14  ,  p2   ,  p1   ,  p3   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1434.6 ];
    end
    
    %tetras to be divided at edge 1-4 & 2-4
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 0 1 0 1 0] ) ,2) ); w( T(w,1) > T(w,2) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p14  ,  p24  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1424.1 ];
    T = [ T ;    p1   ,  p24  ,  p3   ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1424.2 ];
    T = [ T ;    p24  ,  p1   ,  p3   ,  p2   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1424.3 ];
    end
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 0 1 0 1 0] ) ,2) ); w( T(w,1) < T(w,2) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p14  ,  p24  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1424.4 ];
    T = [ T ;    p2   ,  p24  ,  p3   ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1424.5 ];
    T = [ T ;    p14  ,  p1   ,  p3   ,  p2   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1424.6 ];
    end
    
    %tetras to be divided at edge 1-3 & 3-4
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 1 0 0 0 1] ) ,2) ); w( T(w,1) > T(w,4) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p13  ,  p2   ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1334.1 ];
    T = [ T ;    p1   ,  p2   ,  p13  ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1334.2 ];
    T = [ T ;    p34  ,  p2   ,  p4   ,  p1   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1334.3 ];
    end
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 1 0 0 0 1] ) ,2) ); w( T(w,1) < T(w,4) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p13  ,  p2   ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1334.4 ];
    T = [ T ;    p4   ,  p2   ,  p13  ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1334.5 ];
    T = [ T ;    p13  ,  p2   ,  p4   ,  p1   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1334.6 ];
    end
    
    %tetras to be divided at edge 1-3 & 2-3
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 1 0 1 0 0] ) ,2) ); w( T(w,1) > T(w,2) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p13  ,  p23  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1323.1 ];
    T = [ T ;    p1   ,  p23  ,  p13  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1323.2 ];
    T = [ T ;    p23  ,  p1   ,  p2   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1323.3 ];
    end
    w = find( all(  bsxfun(@eq, ~~ET , ~~[0 1 0 1 0 0] ) ,2) ); w( T(w,1) < T(w,2) ) = [];
    if ~isempty(w)
    p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
    T = [ T ;    p13  ,  p23  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1323.4 ];
    T = [ T ;    p4   ,  p23  ,  p13  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1323.5 ];
    T = [ T ;    p13  ,  p1   ,  p2   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1323.6 ];
    end
    
    
    ET( ~any( ET ,2) ,:) = [];
    ET = unique( ~~ET , 'rows' );
    if ~isempty( ET )
      warning('ET is not empty');
      ET
    end
  end
  
  M.tri      = T;
  F(W)       = [];
  M.tri(W,:) = [];            %remove the original faces
  %try, C(W,:)= []; end
  
  
  [F,ord] = sort( F );      %reorder the new faces in their "original position"
  M.tri = M.tri(ord,:);
  %try, C     = C(ord); end
  
  for f = fieldnames( M ).'
    if strcmp( f{1} , 'tri' ), continue; end
    if ~strncmp( f{1} , 'tri' , 3 ), continue; end
    M.(f{1}) = M.(f{1})(F,:);
  end
  %try, M.triCASE = C; end
  
end
