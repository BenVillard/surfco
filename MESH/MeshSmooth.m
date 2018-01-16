function M = MeshSmooth( M , its , varargin )
%   if nargin < 1
%     vtkSmoothPolyDataFilter();
%     return;
%   end

  PLOT = false;
  [varargin,PLOT] = parseargs( varargin , 'plot' , '$FORCE$',{true,PLOT});
  
  ALG = 'vtk';
  [varargin,ALG] = parseargs( varargin , 'vtk'     , '$FORCE$',{'vtk'    ,ALG});
  [varargin,ALG] = parseargs( varargin , 'uniform' , '$FORCE$',{'uniform',ALG});
  [varargin,ALG] = parseargs( varargin , 'cotan'   , '$FORCE$',{'cotan'  ,ALG});
  [varargin,ALG] = parseargs( varargin , 'cotan2'  , '$FORCE$',{'cotan2' ,ALG});

  SCHEME = 'explicit';
  [varargin,SCHEME] = parseargs(varargin,'explicit','$FORCE$',{'explitic',SCHEME});
  [varargin,SCHEME] = parseargs(varargin,'implicit','$FORCE$',{'implicit',SCHEME});

  LAMBDA = 0.1;
  [varargin,~,LAMBDA] = parseargs(varargin,'LAMBDA','$DEFS$',LAMBDA);

  if PLOT
  	figure;
    plotMESH( M , 'EdgeColor','r','FaceColor','none' );
    h = hplotMESH( M );
  end
  
  M0xyz = double( M.xyz );
  for it = its(:).'
    
    if strcmpi( ALG , 'vtk' )
      
      sM = vtkSmoothPolyDataFilter( MakeMesh(M) , 'SetNumberOfIterations' , ceil(it) , ...
              'SetRelaxationFactor'     , LAMBDA ,...
              'SetFeatureAngle'         , 180    ,...
              'SetEdgeAngle'            , 180    ,...
              'SetFeatureEdgeSmoothing' , true   ,...
              'SetBoundarySmoothing'    , true   ,...
              'SetGenerateErrorScalars' , false  ,...
              'SetGenerateErrorVectors' , false  ,...
              varargin{:} );
      M.xyz = sM.xyz;
      
    elseif strcmpi( ALG , 'cotan' ) || strcmpi( ALG , 'uniform' )
      
      if isempty( which( 'laplacian_smooth' ) )
        gptoolboxPath = fullfile(fileparts(fileparts(mfilename('fullpath'))),'thirdParty','gptoolbox');
        addpath( fullfile( gptoolboxPath , 'mesh' ) );
        addpath( fullfile( gptoolboxPath , 'matrix' ) );
      end
      
      factor = 1;
      while 1
        sM = laplacian_smooth( double(M.xyz) , double(M.tri) , ALG , [] , LAMBDA/factor , SCHEME , M.xyz , ceil( it*factor ) );
        if true || isOK( sM )
          break;
        end
        factor = factor * 1.2;
      end
      M.xyz = sM;
      
    elseif strcmpi( ALG , 'cotan2' )

         M = laplacian_smoothing( M, 'iter', it, 'peyre');
        
    else
        
      error('invalid or unknown ALG.');

    end
      
    M.xyz = transform( M.xyz , MatchPoints( M0xyz , M.xyz , 'Gt' ) );
     
    if PLOT
      set( h , 'Vertices',M.xyz );
      drawnow;
    end
  end

end
