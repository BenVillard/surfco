function M = SurFCo( C , varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SurFCo is a tool to compute a surface mesh from biological delineations, contours or segmentations. It has been developped at the Institute of Biomedical Engineering (IBME), at the University of Oxford, under the supervision of Professor Vicente Grau, and Dr. Ernesto Zacur.
% 
% This work has been published in the following publications, as such, if you use the code we would highly appreciate you citing them:
% 
% [1] B. Villard, V. Grau, and E. Zacur, Surface mesh reconstruction from cardiac MRI contours, J. Imaging, vol. 4(1), no. 16, 2018.
% 
% [2] B. Villard, V. Carapella, R. Ariga, V. Grau, and E. Zacur, Cardiac Mesh Reconstruction from Sparse, Heterogeneous Contours. In: Valdés Hernández M., González-Castro V. (Eds.) Medical Image Understanding and Analysis. MIUA 2017. Communications in Computer and Information Science, Vol. 723. Springer, Cham
%
% Please refer to [1] for a detailed workings of the method as well as for parameter choices. 
% 
%
% Authors: Benjamin Villard, Ernesto Zacur <benjamin.villard@eng.ox.ac.uk>
% Copyright © 2018 University of Oxford
% Version: 0.1.0
%
% University of Oxford means the Chancellor, Masters and Scholars of
% the University of Oxford, having an administrative office at
% Wellington Square, Oxford OX1 2JD, UK. 
%
% This file is part of SurFCo.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. The offer of this
% program under the terms of the License is subject to the License
% being interpreted in accordance with English Law and subject to any
% action against the University of Oxford being under the jurisdiction
% of the English Courts.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% Default Parameters

%%%%%%%%%%%%%%% Values obtained from Multi - Param Study %%%%%%%%%%%%%%  
% FARTERPOINTS         = -2    ; % In mm. The minimum distance beetween any point samples.
% FARTHESTP_RESAMPLING = 20    ;
% SMTHDEC_ITER         = 20    ;
% TARGETREDUCTION      = 0.72  ;  %reduction at decimation step
% MAX_DEFORMATION_ITS  = 250   ;
% PERCENTAGE           = 0.1   ;  %percentage of pushing force
% THRES_E              = -Inf  ;
% SMOOTH_STRENGTH      = 200   ;
% LAMBDA               = 100000;
% FAKEBUTTERFLY_ITER   = 0     ;
% SMOOTH_LAMBDA        = 0.1   ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  enableVTK                             ;
  FARTERPOINTS         = []             ;  %number of points to use ( As Neg: breaking distance between points)
  FARTHESTP_RESAMPLING = []             ;  %resample seeds every Nth iteration
  SMTHDEC_ITER         = []             ;  %number of Smoothing + Decimation iteration 
  bLID                 = NaN            ;  %bottom LID
  uLID                 = NaN            ;  %upper LID
  RESOLUTION           = 50             ;  %number of angular partitions in initial mesh
  TARGETREDUCTION      = 0.72            ;  %reduction at decimation step
  MAX_DEFORMATION_ITS  = 250            ;  %Number of Deformation Iterations
  PERCENTAGE           = 0.1            ;  %percentage of pushing force
  VPLOT                = false          ;
  DEFORM               = true           ;
  THRES_E              = 1e-5           ;
  INITIAL_MESH         = []             ;
  LAMBDA               = 100000         ;
  SUBDIV_METHOD        = 'FakeButterfly';  % 'loop'    'butterfly'    'FakeButterfly'    'linear0'    'linear'        'pushing'
                                           %  {'FakeButterfly',0}  perform loop + TPS
                                           %  {'FakeButterfly',100}  perform linear + smooth_100 + TPS
  EAP                  = []             ;
  SMOOTH_STRENGTH      = 200            ;
  SMOOTH_LAMBDA        = 0.1            ;
  SMOOTH_ALG           = 'vtk';  % 'vtk'  'cotan'   'uniform'   'cotan2' (Peyre Implementation)   {'cotan','implicit'}   {'uniform','explicit'}
  
  
  
  [varargin,VPLOT]  = parseargs(varargin,'plot'    ,'$FORCE$',{true,VPLOT}) ;
  [varargin,VPLOT]  = parseargs(varargin,'noplot'  ,'$FORCE$',{false,VPLOT});
  vprintf           = @(varargin)[]                                         ;
  if VPLOT, vprintf = @(varargin)title(sprintf(varargin{:}))                ; end


  [varargin,~,TARGETREDUCTION     ]  = parseargs(varargin,'TARGETREDUCTION'                     ,'$DEFS$',TARGETREDUCTION     );
  [varargin,~,MAX_DEFORMATION_ITS ]  = parseargs(varargin,'MAX_DEFORMATION_ITS','mdi'           ,'$DEFS$',MAX_DEFORMATION_ITS );
  [varargin,~,THRES_E             ]  = parseargs(varargin,'THRES_E'                             ,'$DEFS$',THRES_E             );
  [varargin,~,LAMBDA              ]  = parseargs(varargin,'LAMBDA'                              ,'$DEFS$',LAMBDA              );
  [varargin,~,FARTERPOINTS        ]  = parseargs(varargin,'FARTERPOINTS'                        ,'$DEFS$',FARTERPOINTS        );
  [varargin,~,FARTHESTP_RESAMPLING]  = parseargs(varargin,'FARTHESTP_RESAMPLING','FARTHESTP_RES','$DEFS$',FARTHESTP_RESAMPLING);
  [varargin,~,SMTHDEC_ITER        ]  = parseargs(varargin,'SMTHDEC_ITER'                        ,'$DEFS$',SMTHDEC_ITER        );
  [varargin,~,INITIAL_MESH        ]  = parseargs(varargin,'INITIAL_MESH'                        ,'$DEFS$',INITIAL_MESH        );
  [varargin,~,PERCENTAGE          ]  = parseargs(varargin,'PERCENTAGE'                          ,'$DEFS$',PERCENTAGE          );
  [varargin,~,SUBDIV_METHOD       ]  = parseargs(varargin,'SUBDIV_METHOD'                       ,'$DEFS$',SUBDIV_METHOD       );
  [varargin,~,EAP                 ]  = parseargs(varargin,'ExtraAnchorPoints'                   ,'$DEFS$',EAP                 );
  [varargin,~,SMOOTH_STRENGTH     ]  = parseargs(varargin,'SMOOTH_STRENGTH','STIFFNESS'         ,'$DEFS$',SMOOTH_STRENGTH     );
  [varargin,~,SMOOTH_LAMBDA       ]  = parseargs(varargin,'SMOOTH_LAMBDA'                       ,'$DEFS$',SMOOTH_LAMBDA       );
  [varargin,~,SMOOTH_ALG          ]  = parseargs(varargin,'SMOOTH_ALG'                          ,'$DEFS$',SMOOTH_ALG          );
 
  %% remove empties
  C( all( cellfun( 'isempty' , C ) ,2) ,:) = [];
  
  for c = 1:numel(C)
    thisC = C{c}; if isempty( thisC ), continue; end
    thisC = polyline( thisC );
    for p = 1:thisC.np
      thisC(p) = resample( thisC(p) , '+e' , 1 );
    end
    C{c} = double( thisC );
  end
  
  
  %% join contours in the same row
  for c = 1:size(C,1), C{c,1} = vertcat( C{c,:} ); end
  C = C(:,1);
  
  %% clean contours
  for c = 1:numel(C)
    thisC = C{c};
    thisC( any( ~isfinite(thisC) ,2) ,:) = [];
    d = sqrt( sum( diff( thisC , [] , 1 ).^2 , 2 ) );
    thisC( find( d < 1e-30 )+1 ,:) = [];
    C{c} = thisC;
  end
  C( cellfun( 'isempty' , C ) ) = [];
  
  nC = numel( C );
  allCountourPoints = cell2mat( C(:) );
  
  
  hFig = 0;
  if VPLOT
      hFig = figure( 'IntegerHandle','off','NumberTitle','off','Name','Fitting surface to points','NextPlot','new','color','w','CreateFcn',[]);
      
      for c = 1:nC
        line( C{c}(:,1) , C{c}(:,2) , C{c}(:,3) , 'color','b','linewidth',2,'linestyle','-','marker','.');
      end
      try, line( EAP(:,1) , EAP(:,2) , EAP(:,3) , 'color','r','linewidth',1,'linestyle','none','marker','+'); end
      view(3); axis('equal');
  end
  
  
  
  %% initial MESH
  if isempty( INITIAL_MESH )
    [varargin,~,bLID] = parseargs( varargin , 'BottomLID', 'bLID' , '$DEFS$' , bLID );
    [varargin,~,uLID] = parseargs( varargin , 'UpperLID' , 'uLID' , '$DEFS$' , uLID );
    [varargin,~,RESOLUTION] = parseargs( varargin , 'RESOLUTION' , '$DEFS$' , RESOLUTION );

    M = computeRuledSurface( C , bLID , uLID , RESOLUTION , VPLOT );
  else
    M = INITIAL_MESH;
  end

  
  if VPLOT
    hM = patch('vertices',M.xyz,'Faces',M.tri,'EdgeColor',[0,0.6980,0.8941],'FaceColor',[0.5294,0.8078,0.9216],'FaceAlpha',0.3);
  end

  if MAX_DEFORMATION_ITS <= 0
    return;
  end
  
  if isempty( INITIAL_MESH )
    M = SUBDIVIDE( M , SUBDIV_METHOD );
    M = SMOOTH( M , SMOOTH_STRENGTH  , SMOOTH_ALG , SMOOTH_LAMBDA );
    M = DECIMATE( M , TARGETREDUCTION );
    if VPLOT, set( hM ,'vertices', M.xyz , 'Faces' , M.tri ); drawnow; end
  end    
  
  %%%%%%%%% %END% 
  %% Construction of the MESH
  if isempty( FARTERPOINTS )
    FARTERPOINTS = max( 150 , 25*nC );
  end
  if isempty( FARTHESTP_RESAMPLING )
    FARTHESTP_RESAMPLING = ceil(MAX_DEFORMATION_ITS./10);
  end
  if isempty( SMTHDEC_ITER )
    SMTHDEC_ITER = 20;
  end
  
  SCALE = max( vec( ipd( allCountourPoints , [] ) ) );
  
  M.xyz = M.xyz/SCALE;

  if VPLOT
    set( hM ,'vertices',SCALE*M.xyz,'Faces',M.tri);
    hP = line( NaN , NaN , NaN , 'marker','o','color',[0 0 0],'linestyle','none','markerfacecolor',[0,0.4,0],'markersize',6);
    drawnow;
  end
  
  percentage = @( from , to ) from + ( to - from ) * PERCENTAGE;
  
  Ep = Inf; last_ds = 0;
  for it = 0:MAX_DEFORMATION_ITS
    if it > 0 && ~rem( it, SMTHDEC_ITER ) && it < MAX_DEFORMATION_ITS
      Ep = Inf;

      M = SUBDIVIDE( M , SUBDIV_METHOD );
      M = SMOOTH( M , SMOOTH_STRENGTH  , SMOOTH_ALG , SMOOTH_LAMBDA );
      M = DECIMATE( M , TARGETREDUCTION );
      if VPLOT, set( hM ,'vertices',SCALE*M.xyz,'Faces',M.tri); drawnow; end
    end
    
    % Redo Farthest Point sampling every Nth Iteration
    if ~rem( it, min( 1e300 , FARTHESTP_RESAMPLING ) ) && it < MAX_DEFORMATION_ITS
      Ep = Inf;
      
      [~,~,ds] = vtkClosestElement( MakeMesh( M.xyz*SCALE , M ) , allCountourPoints );
      [~,ds] = sort( ds ); ds = ds(end-1);
      
      if ds ~= last_ds
        if FARTERPOINTS < 0
          P = FarthestPointSampling( allCountourPoints, ds , -FARTERPOINTS , Inf );
        else
          P = FarthestPointSampling( allCountourPoints, ds , -1 , FARTERPOINTS );
        end
        last_ds = ds;
        P = [ P ; EAP ];
        P = P/SCALE;
      end
      if VPLOT
        set( hP ,'XData',P(:,1)*SCALE,'YData',P(:,2)*SCALE,'ZData',P(:,3)*SCALE);
        set( hP ,'markerfacecolor',rand(1,3) );
      end
    end
    
    [~,cp,d] = vtkClosestElement( M , P ); E = max( d ).^2;    
    vprintf('it: %3d  , E = %g\n', it , E );
    if E < THRES_E
      break; end;
%     if E > Ep, M = Mp;
%       break; end;
%     disp(['E = ' num2str(E) '  ThreshE = ' num2str(THRES_E) '   Ep = ' num2str(Ep) ]);

    w = true( numel(d) ,1);
%     [~,w] = sort(d); w = w( 1:end-1 );
%     w = d < max(d)*2;
    Ep = E; Mp = M;
    M.xyz = InterpolatingSplines( cp(w,:) , percentage( cp(w,:) , P(w,:) ) , M.xyz , 'r' ,'LAMBDA' , LAMBDA/size(cp,1) ); %LAMBDA/size(cp,1)
    if VPLOT, set( hM ,'Vertices',SCALE*M.xyz,'Faces',M.tri); drawnow; end
  end
  
  M.xyz = M.xyz*SCALE;
  
  if VPLOT
    for c = 1:nC
      try
        xyz = SliceMesh( M , getPlane( C{c} ) );
        line(xyz(:,1),xyz(:,2),xyz(:,3) ,'color','m','linewidth',2,'linestyle','-');
      end
    end
  end
  
%   if 1, set(0,'DefaultFigureWindowStyle','docked'); end
  
end

function X = orientCurve( X )
  o = convhull( X(:,1:2) );
  o(end) = [];
  o = circshift( o , 1-minInd(o) );

  if ~issorted( o )
    X = flipdim( X , 1);
  end
end
function [M,C,D] = computeRuledSurface( C , bLID , uLID , RESOLUTION , VPLOT )
  %% get the normals of each contour
  nC = numel(C);
  Ns = NaN(nC,3);
  for c = 1:nC
    Ns(c,:) = getPlane( C{c} , '+z' , 'normal' );
  end
  
  %% compute the "mode" of these normals
  dN2N = abs( 2*min( asind( ipd( Ns , Ns )/2 ) , asind( ipd( -Ns , Ns )/2 ) ) );
  
  [~,N] = min( sum( log( 180 + dN2N ) , 1 ) );
  w = dN2N( N ,:) < 20;

  Z = meanNormal( Ns( w ,:) ); Z = Z(:);
  if Z(3) < 0, Z = -Z; end
  
  C = C(w) ; nC = numel( C );
  
  %% sort contours, from bottom to up
  Zs = NaN(nC,1);
  for c = 1:nC
    Zs(c) = mean( C{c} * Z );
  end
  [Zs,ord] = sort( Zs , 'ascend' );
  C  = C( ord );
  P = C{1}(1,:);

  %% rotation pointing upwards
  R = [ null( Z.' ) , Z ];
  if det( R ) < 0, R(:,1) = -R(:,1); end
  if     isnan( bLID )
  elseif isscalar( bLID )  &&   bLID > 0
    C = [     ( mean( C{ 1 } * R , 1 ) - [0 0 bLID] ) * R.' ; C ]; nC = nC + 1;
  elseif isscalar( bLID )  &&   bLID < 0
    error('to be implemented');
%     C = [ bsxfun( @plus , C{end} * R  , [0 0 -bLID] ) * R.' ; C ]; nC = nC + 1;
  elseif numel( bLID ) == 3
    C = [  bLID(:).'  ;  C  ]; nC = nC + 1;
  else
    error('unknown specification of bLID.');
  end
  

  if     isnan( uLID )
  elseif isscalar( uLID )  &&   uLID > 0
    C = [ C ; ( mean( C{end} * R , 1 ) + [0 0 uLID] ) * R.'     ]; nC = nC + 1;
  elseif isscalar( uLID )  &&   uLID < 0
    C = [ C ; bsxfun( @plus , C{end} * R  , [0 0 -uLID] ) * R.'     ]; nC = nC + 1;
  elseif numel( uLID ) == 3
    C = [ C; uLID(:).'  ]; nC = nC + 1;
  else
    error('unknown specification of bLID.');
  end
  
  
  if ~~VPLOT
    for c = 1:nC
      line( C{c}(:,1) , C{c}(:,2) , C{c}(:,3) , 'Marker','.' , 'Color' , [ 1-(c/(2*nC)) , 0 , 0 ],'LineStyle','none' );
    end
  end
  
  
  M.xyz = [];  M.tri = [];
  for c = 1:nC
    thisC = C{c};
    if size( thisC , 1 ) > 1
      thisC = thisC * R;
      if isequal( thisC(1,:) , thisC(end,:) ), thisC( end ,:) = []; end  %open if it is closed
      thisC = orientCurve( thisC );                                      %orient counter-clockwise
      thisC = thisC * R.';
      d2P = ipd( thisC , P );
      [~,P] = min( d2P );
      thisC = circshift( thisC , [ 1-P , 0 ] );
      P = thisC(1,:);
      thisC = thisC( [1:end,1] , :);                                     %close the contour again
      L = [ 0 ; cumsum( sqrt( sum( diff( thisC , [] , 1 ).^2 , 2 ) ) ) ];
      
      w = find( diff(L)==0 ) + 1; thisC(w,:) = []; L(w,:) = [];
      
      
      try
        new_xyz = Interp1D( thisC , L , linspace(0,L(end),RESOLUTION) );
      catch
        new_xyz(:,1) = interp1( L , thisC(:,1) , linspace(0,L(end),RESOLUTION) );
        new_xyz(:,2) = interp1( L , thisC(:,2) , linspace(0,L(end),RESOLUTION) );
        new_xyz(:,3) = interp1( L , thisC(:,3) , linspace(0,L(end),RESOLUTION) );
      end
    else
      new_xyz = ones( RESOLUTION , 1 ) * thisC;
    end
    if VPLOT && 0
      line( new_xyz(:,1) , new_xyz(:,2) , new_xyz(:,3) , 'marker','o','color',[0 0 0],'linestyle','none','markerfacecolor','r','markersize',5);
    end
    new_xyz([1 end],:) = [1;1]*mean( new_xyz([1 end],:) ,1);
    M.xyz = [ M.xyz ; new_xyz ];
    M.tri = [ M.tri ;...
      (c-1)*RESOLUTION + [ 1:RESOLUTION-1 ; 2:RESOLUTION              ; RESOLUTION+1:2*RESOLUTION-1 ].' ;...
      (c-1)*RESOLUTION + [ 2:RESOLUTION   ; RESOLUTION+2:2*RESOLUTION ; RESOLUTION+1:2*RESOLUTION-1 ].' ];
  end
  M.tri( any( M.tri > size(M.xyz,1) , 2 ) , : ) = [];
  
  M = TidyMesh( M , 0 , true );
  M.tri( meshQuality( M , 'area' ) == 0 , :) = [];

end
function Z = meanNormal( N )
  M    = @(ae) [ cos(ae(2)) * cos(ae(1)) , cos(ae(2)) * sin(ae(1)) , sin(ae(2)) ];
  dN2N = @(a,b) 2*min( asin(ipd(a,b)/2) , asin(ipd(a,-b)/2) );
  E    = @(m) sum( abs( dN2N( m , N ) ).^1 );
  
  m = mean(N,1);
  [ m(1),m(2),m(3) ] = cart2sph( m(1),m(2),m(3) );
  m = m([1,2]);
  try
        m = Optimize( @(ae)E(M(ae)) , m , 'methods',{'conjugate','coordinate',1},'ls',{'quadratic','golden','quadratic'},struct('COMPUTE_NUMERICAL_JACOBIAN',{{'a'}}),'noplot','verbose',0);
  catch
        m = fminsearch(@(ae) E(M(ae)),m);
  end
  Z = M(m);
  if 1, set(0,'DefaultFigureWindowStyle','docked'); end
end

function M2Sub = SUBDIVIDE( M2Sub , SUBDIV_METHOD )
  FAKEBUTTERFLY_ITER = 0;
  if iscell( SUBDIV_METHOD )
    FAKEBUTTERFLY_ITER = SUBDIV_METHOD{2};
    SUBDIV_METHOD = SUBDIV_METHOD{1};
  end

  method = SUBDIV_METHOD;
  M2Sub = MakeMesh( M2Sub );
  switch lower( method )
    case 'loop'
      M2Sub = vtkLoopSubdivisionFilter( M2Sub );

    case 'butterfly'
      Ms = vtkButterflySubdivisionFilter( M2Sub  );
      d_check  = vtkClosestElement__o3( M2Sub , Ms.xyz );
      if any( d_check > prctile( d_check , 95 ) * 10 )
        M2Sub = MeshSubdivide( M2Sub );
      else
        M2Sub = Ms;
      end

    case 'fakebutterfly'
      M2Sub = MeshFakeButterfly( M2Sub, FAKEBUTTERFLY_ITER );

    case 'linear0'
      M2Sub = MeshSubdivide( M2Sub );

    case 'linear'
      %linear + smooth
      M2Sub = MeshSubdivide( M2Sub );
      M2Sub = vtkSmoothPolyDataFilter( M2Sub ,'SetNumberOfIterations' , 25 ,...
        'SetFeatureAngle'           , 180  ,...
        'SetEdgeAngle'              , 180  ,...
        'SetConvergence'            , 0    ,...
        'SetFeatureEdgeSmoothing'   , true ,...
        'SetBoundarySmoothing'      , true );

    case 'pushing'
      %linear + smooth + pushing back
      M2Sub = MeshSubdivide( M2Sub );
      M2Sub = MeshSmooth( M2Sub , 25 );
  end

end
function M = SMOOTH( M , its , alg , lambda )
  if ~iscell( alg ), alg = {alg}; end
  M = MeshSmooth( M , its  , alg{:} , 'lambda' , lambda );
end
function M = DECIMATE( M , tr )
  M = vtkQuadricDecimation( M , 'SetTargetReduction', tr );
end
