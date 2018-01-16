function [M,iM] = getPlane( I , varargin )
% Returns a homogeneous matrix to represent the plane P.
% The plane P can be described by:
% - [ point(:)' ; normal(:)' ] , a 2 x 3 matrix
% - [ p1(:)' ; p2(:)' ; p3(:)' ] , a 3 x 3 matrix where each row is a point lying on the plane
% - [ p1(:)' ; p2(:)' ; ... ; pn(:)' ] , a list of points (as rows) lying on the plane.
%   In this case the plane is computed as the best linear fit to the points (in lsq sense).
% - other representations.... are wellcome.
%   

  MODE = 'matrix';

  R = eye(3);
  P = [0;0;0];

  CHECK_ORIENTATIONS = true;
  if          0
  elseif      isa( I , 'I3D' )
    I = I.centerGrid;

    
    R = I.SpatialTransform(1:3,1:3);
    P = I.center(:);
    P = I.SpatialTransform(1:3,4);
    
  elseif      isstruct( I )  &&  isfield( I , 'ImageOrientationPatient' )  &&  isfield( I , 'ImagePositionPatient' )

    R = reshape( I.ImageOrientationPatient , 3 , 2 );
    R(:,3)= cross( R(:,1), R(:,2) );
    for cc = 1:3, for it = 1:5, R(:,cc) = R(:,cc)/sqrt( R(:,cc).' * R(:,cc) ); end; end

    P = I.ImagePositionPatient(:);

  elseif      isequal( size(I) , [4,4] )

    R = I(1:3,1:3);
    R(:,3) = cross2( R(:,1).' , R(:,2).' );
    for it = 1:5, R(:,1) = R(:,1)/sqrt( R(:,1).'*R(:,1) ); end
    for it = 1:5, R(:,2) = R(:,2)/sqrt( R(:,2).'*R(:,2) ); end
    for it = 1:5, R(:,3) = R(:,3)/sqrt( R(:,3).'*R(:,3) ); end
    
    P = I(1:3,4);
    
  elseif  isequal( size(I) , [2,3] )

    I(2,:) = I(2,:) / sqrt( I(2,:) * I(2,:).' );
    CHECK_ORIENTATIONS = false;
    R = [ null( I(2,:) ) , I(2,:).' ];

    P = I(1,:).';
    
  elseif  isequal( size(I) , [1,4] )

    R = I(1:3);
    [~,id] = max(abs(R));
    P = [ zeros(id-1,1) ; -I(4)/R(id) ; zeros(3-id,1) ];

    for it = 1:5, R = R/sqrt( R(:).'*R(:) ); end
    R = [ null( R(:).' ) , R(:) ];
    
  elseif isnumeric( I ) && size( I,2 ) == 3 && size( I,1 ) >= 3
    
    I( any( ~isfinite(I) , 2 ) ,:) = [];
    
    P = mean( I , 1 );

    [~,d,R] = svd( bsxfun( @minus , I , P ) , 0 );

    P = P(:);
    
  elseif isstruct( I )  &&  isfield( I , 'SpatialTransform' )
    
    R = I.SpatialTransform(1:3,1:3);
    P = I.SpatialTransform(1:3,4);
    
  else
    
    error('not implemented for this input.');
    
  end
  
  if det( R(1:3,1:3) ) < 0
    R(:,[1 2]) = R(:,[2 1]);
%     R(:,3) = -R(:,3);
  end  

  offset = [];
  for v = 1:numel(varargin)
    if ischar( varargin{v} )
      switch lower( varargin{v} )
        case {'matrix','m'},                MODE = 'matrix';
        case {'3p','3points'},              MODE = '3points';
        case {'4p','4points'},              MODE = '4points';
        case {'pn','pointnormal'},          MODE = 'pointnormal';
        case {'n','normal'},                MODE = 'normal';
        case {'abcd','equation'},           MODE = 'abcd';
        case {'+x','x+'}
          if CHECK_ORIENTATIONS && R(1,3) < 0, R = R * diag([1,-1,-1]); end
        case {'+y','y+'}
          if CHECK_ORIENTATIONS && R(2,3) < 0, R = R * diag([1,-1,-1]); end
        case {'+z','z+'}
          if CHECK_ORIENTATIONS && R(3,3) < 0, R = R * diag([1,-1,-1]); end
        case {'-z','z-'}
          if CHECK_ORIENTATIONS && R(3,3) > 0, R = R * diag([1,-1,-1]); end
  %       case {'yvector'}
  %         Y = varargin{v+1};
      end
    elseif isnumeric( varargin{v} ) && isscalar( varargin{v} ) && isempty( offset )
      offset = varargin{v};
    else
      error('unknown argument');
    end
  end
  
  
  M = [ R , P ; 0 , 0 , 0 , 1 ];
  if ~isempty( offset )
    M = M * [ eye(3) , [0;0;offset] ; 0 0 0 1 ];
  end
  
  if nargout > 1, iM = minv( M ); end
  
  transform = @(x,M)bsxfun( @plus , x * M(1:3,1:3).' , M(1:3,4).' );
  switch lower( MODE )
    case 'matrix'
    case '3points',     M = transform( [0 0 0;1 0 0;0 1 0] , M );
    case '4points',     M = transform( [0 0 0;1 0 0;1 1 1;0 1 0] , M );
    case 'pointnormal', M = [ M(1:3,4).' ; M(1:3,3).' ];
    case 'normal',      M = M(1:3,3);
    case 'abcd',        M = [ M(1:3,3).' , -( M(1:3,3).'*M(1:3,4) ) ];
  end
  
end
