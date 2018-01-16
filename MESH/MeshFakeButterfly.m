function Mb = MeshFakeButterfly( M , varargin )


  Nsmooth = 0;
  if numel( varargin ) && isnumeric( varargin{1} )
    Nsmooth = varargin{1};
    varargin(1) = [];
  end
  
  if ~Nsmooth
  
    Mb = vtkLoopSubdivisionFilter( M );
    
  else

    Mb = MeshSubdivide( M );
    Mb = vtkSmoothPolyDataFilter( Mb ,...
        'SetNumberOfIterations', Nsmooth ,...
        'SetFeatureEdgeSmoothing',true ,...
        'SetFeatureAngle',180 ,...
        'SetEdgeAngle', 180 ,...
        'SetBoundarySmoothing',true ,...
        varargin{:} );
  
  end
    
  from = Mb.xyz( 1:size( M.xyz ,1) ,:);
  to   = M.xyz;
  
  if size( from , 1 ) > 500
    [~,ids] = FarthestPointSampling( from , 1 , 0 , 500 );
    from = from(ids,:);
    to   =   to(ids,:);
  end
  
  
  
  Mb.xyz = InterpolatingSplines( from , to , Mb.xyz , 'r' );


end