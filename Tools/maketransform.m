function [H,dH] = maketransform( varargin )

if false

  tt = @(p) {'center',p(1:3),'cl_sxyzt',p(3+(1:7)),'inv',reshape(p(10+(1:16)),[4 4]),'inv','center',{[5 2 -3]},'ct',reshape(p(26+(1:9)),[3 3]),'inv'};
  p = rand(1,3+7+16+9);

  tt = @(p) {'center',p(1:3),'ct',reshape(p(3+(1:9)),[3 3]),'inv'};
  p = rand(1,3+9);

  tt = @(p) {'center',p(1:3),'rxyx',p(3+(1:3)),'inv','center',{[20 30 50]},'s',p(3+3+1),'rzyz',p(3+3+1+(1:3)),'resetcenter','ryxy',p(3+3+1+3+(1:3)),'t',p(3+3+1+3+3+(1:3)),'center',{[-3 -4 -5]},'q',p(16+(1:3))};
  p = rand(1,19);


  tt = @(p) {'center',p(1:3),'rxyx',p(3+(1:3)),'inv','center',{[20 30 50]},'l_s',p(3+3+1),'rzyz',p(3+3+1+(1:3)),'resetcenter','ryxy',p(3+3+1+3+(1:3)),'t',p(3+3+1+3+3+(1:3)),'center',{[-3 -4 -5]},'q',p(16+(1:3)),'mx'};
  p = rand(1,19);

  %TT = tt(p); maketransform( TT{:} , 'type' )
  

  tt = @(p) {'center',p(1:3),'mx','s',{2},'rxyx',{[10 20 30]}};
  p = rand(1,3);

  
  tt = @(p) {'center',{[10 -20 15]},'mx','s',p,'rxyx',{[10 20 30]}};
  p = 2;

  [H,dH] = maketransform( tt , p );
  
  dHn = NumericalDiff( @(p) maketransform(tt,p) , p );
  
  maxnorm( dH , dHn )
  
  
end

%   if nargin == 2 && isa( varargin{1} , 'function_handle' )
%     args = feval( varargin{1} , varargin{2} );
%     if ~iscell( args ), args = {args}; end
%     if nargout > 1
%       [H,dH] = maketransform( args{:} );
%     else
%       H      = maketransform( args{:} );
%     end
%     return;
%   end

  
  IN_HOMOGENEOUS_COORDINATES = true;
  IN_2D                      = false;
  RETURN_TTYPE = false;

  H     = eye(4);
  TTYPE = 'E';

  compute_dH = false;
  if nargout > 1
    compute_dH = true;
    dH = zeros(16,0);
  end

  if nargin < 1
    rearrangeOutputs();
    return;
  end

  % defaults
  center = [0;0;0]; CENTER_AT_ZERO = true;
   CENTER = eye(4);  CENTER(1:3,4) =  center; 
  ICENTER = eye(4); ICENTER(1:3,4) = -center; 
  if compute_dH
    CENTER_idx = [];
    dC  = zeros(16,0);
    dIC = zeros(16,0);
  end
  
  
  runits = pi/180;     % degree
  

  i = 1;
  while i <= numel(varargin)
    arg = varargin{i};   i=i+1;
    
    if isa( arg , 'function_handle' )
     i = i-1;
     params = varargin{i+1};
     model  = feval( arg , params );
     
     varargin = [ varargin(1:i-1) , model , varargin(i+2:end) ];
     continue;
    end

    if iscell( arg ) || isnumeric( arg )
      i = i-1;
      varargin = [ varargin(1:i-1) , 'h' , varargin(i:end) ];
      continue;
    end
    
    dp  = 1;
    arg = lower( arg );
    switch arg
      case {'l'}
        if nargin ~= 2, error('with ''l'' option, only 2 args expected'); end

        p = varargin{i}; p = double( p(:) );
        switch numel(p)
          case 1             %scale factor only
            H       = diag([p p p]);
            H(4,4)  = 0;
          case 2             %vol_preserving scale
            H = diag([p(1) p(2) -p(1)-p(2)]);
            H(4,4)  = 0;
          case 3             %rotation
            H = tangentmatrix( p(1) , p(2) , p(3) );
          case 4             %isotropic scale + rotation
            H = tangentmatrix( p(2) , p(3) , p(4) );
            H([1 6 11]) = p(1);
          case 6             %rotation + traslation
            H = tangentmatrix( p(1) , p(2) , p(3) );
            H(1:3,4) = p(4:6);
          case 7             %isotropic scale + rotation + traslation
            H = tangentmatrix( p(2) , p(3) , p(4) );
            H([1 6 11]) = p(1);
            H(1:3,4) = p(5:7);
          case 8             %volume preserving affine3x3
            H = zeros(3,3);
            H(1:8) = p;
            H(9)   = -H(1,1) -H(2,2);
            H(4,4) = 0;
          case 9             %affine3x3
            H = zeros(3,3);
            H(1:9) = p;
            H(4,4) = 0;
          case 11            %volume preserving affine3x4
            H = zeros(3,3);
            H(1:8) = p(1:8);
            H(1:3,4) = p(9:11);
            H(3,3) = -H(1,1) -H(2,2);
            H(4,4) = 0;
          case 12            %affine3x4
            H = reshape( p , 3 , 4 );
            H(4,4) = 0;
        end
      
        if compute_dH
          dH = NumericalDiff( @(p) maketransform('l',p) , p , 'i' );
        end
        
        rearrangeOutputs();
        
        return
      case {'resetcenter' 'centeratzero' 'centeraatrigin' 'reset'}
         CENTER = eye(4);
        ICENTER = eye(4);
        
        if compute_dH
          CENTER_idx = [];
          dC  = zeros(16,0);
          dIC = zeros(16,0);
        end
        CENTER_AT_ZERO = true;
      case {'noh','nohomogeneous','3x3'}
        IN_HOMOGENEOUS_COORDINATES = false;
      case {'2d'}
        IN_2D = true;
      case {'type'}
        RETURN_TTYPE = true;
      case {'angleunits' 'units' 'runits' 'rotationunits'}
        runits = varargin{i}; i=i+1;
        if ischar( runits )
          switch lower( runits )
            case {'deg','degrees','degree'},  runits = pi/180;
            case {'rad','radians','radian'},  runits = 1;
          end
        end
      case {'degrees','deg' }, runits = pi/180;
      case {'radians','rad' } , runits = 1;

        
      case {'c' 'center' 'centerat' 'setcenter' 'setcenterat'}
        p = varargin{i}; i = i+1;

        if      iscell( p ) && numel(p) == 2,    dcenter = p{2};                  center = p{1};
        elseif  iscell( p ) && numel(p) == 1,    dcenter = zeros(numel(p{1}),0);  center = p{1};
        else
          dcenter = eye(3);
          center  = double(p);
        end
        if numel( center ) == 2, center = [ center(:) ; 0 ];     end
        if numel( center ) ~= 3, error('3 parameters expected'); end
        
         CENTER = eye(4);  CENTER(1:3,4)=  center;
        ICENTER = eye(4); ICENTER(1:3,4)= -center;
        
        if compute_dH
          CENTER_idx = size( dH , 2 ) + ( 1:size(dcenter,2) );
          dH(:,CENTER_idx) = 0;
          
          dC  = [0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0; 1,0,0;0, 1,0;0,0, 1;0 0 0]*dcenter;
          dIC = [0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;-1,0,0;0,-1,0;0,0,-1;0 0 0]*dcenter;
        end
        
        CENTER_AT_ZERO = false;

        
      case {'inv','inverse','i'}
        H = minv( H );

        try
        if max( abs( H(4,:) - [0 0 0 1] ) ) < 100*eps(1)
          H(4,:) = [0 0 0 1];
        end
        end
        
        if compute_dH
          dH = - kron( H.' , H )*dH;
        end

        
      case {'mx','mirrorx'}
        t = diag([-1 1 1 1]);

        updateCenteredH( t , zeros(16,0) );
        TTYPE = updateTTYPE( 'L' , ~CENTER_AT_ZERO );
        
        
      case {'my','mirrory'}
        t = diag([1 -1 1 1]);

        updateCenteredH( t , zeros(16,0) );
        TTYPE = updateTTYPE( 'L' , ~CENTER_AT_ZERO );

        
      case {'mz','mirrorz'}
        t = diag([1 1 -1 1]);

        updateCenteredH( t , zeros(16,0) );
        TTYPE = updateTTYPE( 'L' , ~CENTER_AT_ZERO );
        
        
      case {'mirror','m'}
        t = diag([-1 -1 -1 1]);

        updateCenteredH( t , zeros(16,0) );
        TTYPE = updateTTYPE( 'L' , ~CENTER_AT_ZERO );
        
        
      otherwise
        p = varargin{i}; i = i+1;
        if      iscell( p ) && numel(p) == 2,    dp = double( p{2} );        p = double( p{1} );
        elseif  iscell( p ) && numel(p) == 1,    dp = zeros(numel(p{1}),0);  p = double( p{1} );
        end
        
        t  = eye(4);
        if isempty( dp ),  dt = zeros(16,0);
        else            ,  dt = 1;
        end

        switch arg
          case {'h'}
            if all( size(p) == [2 2] )
              t(1:2,1:2) = p;

              if compute_dH  && ~isempty(dp)
                dt = [1,0,0,0;0,1,0,0;0,0,0,0;0,0,0,0;0,0,1,0;0,0,0,1;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0]*dp;
              end
              hasT = false;

            elseif all( size(p) == [3 3] )
              t(1:3,1:3) = p;

              if compute_dH  && ~isempty(dp)
                dt = [1,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,1,0,0,0,0,0;0,0,0,0,1,0,0,0,0;0,0,0,0,0,1,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,0,1;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0]*dp;
              end
              hasT = false;

            elseif all( size(p) == [4 4] )
              t  = p;

              if compute_dH  && ~isempty(dp)
                dt = eye(16)*dp;
              end
              hasT = true;

            elseif isempty( p )
              
              hasT = false;
              
            else, error('a 4x4 matrix expected'); end

            updateH( t , dt );
            TTYPE = updateTTYPE( 'G' , hasT );


          case {'ct','centeredtransform'}
            if all( size(p) == [2 2] )
              t(1:2,1:2) = p;

              if compute_dH  && ~isempty(dp)
                dt = [1,0,0,0;0,1,0,0;0,0,0,0;0,0,0,0;0,0,1,0;0,0,0,1;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0]*dp;
              end
              hasT = false;

            elseif all( size(p) == [3 3] )
              t(1:3,1:3) = p;

              if compute_dH  && ~isempty(dp)
                dt = [1,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,1,0,0,0,0,0;0,0,0,0,1,0,0,0,0;0,0,0,0,0,1,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,0,1;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0]*dp;
              end
              hasT = false;

            elseif all( size(p) == [4 4] )
              t  = p;

              if compute_dH  && ~isempty(dp)
                dt = eye(16)*dp;
              end
              hasT = true;

            elseif isempty(p)
              
              hasT = false;
              
            else, error('a 4x4 matrix expected'); end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'G' , hasT );


          case {'general','linear','generallinear'}
            if numel(p) ~= 12, error('12 parameters expected'); end

            t(1:3,1:4) = reshape( p , [3,4] );
            if compute_dH  && ~isempty(dp)
              dt = [1,0,0,0,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,1,0,0,0,0,0,0,0,0;0,0,0,0,1,0,0,0,0,0,0,0;0,0,0,0,0,1,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,1,0,0,0,0,0;0,0,0,0,0,0,0,1,0,0,0,0;0,0,0,0,0,0,0,0,1,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,0,0,0,0,1;0,0,0,0,0,0,0,0,0,0,0,0]*dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'G' , true );


          case {'general9','linear9','generallinear9','general3x3','linear3x3','generallinear3x3'}
            if numel(p) ~= 9, error('9 parameters expected'); end

            t(1:3,1:3) = reshape( p , [3 3] );
            if compute_dH  && ~isempty(dp)
              dt = [1,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,1,0,0,0,0,0;0,0,0,0,1,0,0,0,0;0,0,0,0,0,1,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,0,1;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0]*dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'G' , ~CENTER_AT_ZERO );


          case {'l_sz','cl_sz'}
            if numel(p) ~= 2, error('2 parameters expected'); end

            L = tangentmatrix( 0 , 0 , p(2));
            L([ 1  6]) = p(1);
            t = expm( L ); t(4,:) = [0 0 0 1];

            if compute_dH  && ~isempty(dp)
              dt = d_expm( L )*[1,0;0,1;0 0;0 0;0,-1;1,0;0 0;0 0;0 0;0 0;0 0;0 0;0,0;0,0;0 0;0 0]*dp;
            end

            if strncmp( arg , 'c' , 1 ),  updateCenteredH( t , dt );
            else                       ,          updateH( t , dt );
            end
            TTYPE = updateTTYPE( 'M' , true );

          
          case {'l_ztt','cl_ztt'}
            if numel(p) ~= 3, error('3 parameters expected'); end

            L = tangentmatrix( 0 , 0 , p(1));
            L([13 14]) = p(2:3);
            t = expm( L ); t(4,:) = [0 0 0 1];

            if compute_dH  && ~isempty(dp)
              dt = d_expm( L )*[0 0 0;1,0,0;0 0 0;0 0 0;-1,0,0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0,1,0;0,0,1;0 0 0;0 0 0]*dp;
            end

            if strncmp( arg , 'c' , 1 ),  updateCenteredH( t , dt );
            else                       ,          updateH( t , dt );
            end
            TTYPE = updateTTYPE( 'R' , true );

            
          case {'l_sztt','cl_sztt'}
            if numel(p) ~= 4, error('4 parameters expected'); end

            L = tangentmatrix( 0 , 0 , p(2));
            L([ 1  6]) = p(1);
            L([13 14]) = p(3:4);
            t = expm( L ); t(4,:) = [0 0 0 1];

            if compute_dH  && ~isempty(dp)
              dt = d_expm( L )*[1,0,0,0;0,1,0,0;0 0 0 0;0 0 0 0;0,-1,0,0;1,0,0,0;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;0,0,1,0;0,0,0,1;0 0 0 0;0 0 0 0]*dp;
            end

            if strncmp( arg , 'c' , 1 ),  updateCenteredH( t , dt );
            else                       ,          updateH( t , dt );
            end
            TTYPE = updateTTYPE( 'M' , true );

            
          case {'l_zstt','cl_zstt'}
            if numel(p) ~= 4, error('3 parameters expected'); end

            L = tangentmatrix( 0 , 0 , p(1));
            L([ 1  6]) = p(2);
            L([13 14]) = p(3:4);
            t = expm( L ); t(4,:) = [0 0 0 1];

            if compute_dH  && ~isempty(dp)
              dt = d_expm( L )*[0,1,0,0;1,0,0,0;0 0 0 0;0 0 0 0;-1,0,0,0;0,1,0,0;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;0,0,1,0;0,0,0,1;0 0 0 0;0 0 0 0]*dp;
            end

            if strncmp( arg , 'c' , 1 ),  updateCenteredH( t , dt );
            else                       ,          updateH( t , dt );
            end
            TTYPE = updateTTYPE( 'M' , true );

            
          case {'l_ss','cl_ss'}
            if numel(p)==1
              UNIFORM = 'U';

              t = diag( [ exp([p;p]) ; 1 ; 1  ] );
              if compute_dH  && ~isempty(dp)
                dt = [ t(1:6).' ; zeros(10,1) ]*dp;
              end

            elseif numel(p) == 2
              UNIFORM = 'S';

              t = diag( [ exp(p(:)) ; 1 ; 1 ] );
              if compute_dH  && ~isempty(dp)
                dt = [t(1),0;0 0;0 0;0 0;0 0;0,t(6);0 0;0 0;0 0;0 0;0,0;0 0;0 0;0 0;0 0;0 0]*dp;
              end

            else, error('1 or 2 parameters expected'); end

            if strncmp( arg , 'c' , 1 ),  updateCenteredH( t , dt );
            else                       ,          updateH( t , dt );
            end
            TTYPE = updateTTYPE( UNIFORM , strncmp( arg , 'c' , 1 ) && ~CENTER_AT_ZERO );

            
          case {'l_vp','l_volpreserving','l_volumepreserving','cl_vp','cl_volpreserving','cl_volumepreserving'}
            if numel(p) ~= 11, error('11 parameters expected'); end

            L = zeros(4,4); L([1 2 3 5 6 7 9 10 13 14 15]) = p; L(3,3) = -L(1,1) -L(2,2);
            t = expm( L ); t(4,:) = [0 0 0 1];

            if compute_dH  && ~isempty(dp)
              dt = d_expm( L )*[1,0,0,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0;0,0,0,1,0,0,0,0,0,0,0;0,0,0,0,1,0,0,0,0,0,0;0,0,0,0,0,1,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,1,0,0,0,0;0,0,0,0,0,0,0,1,0,0,0;-1,0,0,0,-1,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,0,0,0,1;0,0,0,0,0,0,0,0,0,0,0]*dp;
            end

            if strncmp( arg , 'c' , 1 ),  updateCenteredH( t , dt );
            else                      ,          updateH( t , dt );
            end
            TTYPE = updateTTYPE( 'V' , true );


          case {'l_vp9','l_volpreserving9','l_volumepreserving9','l_vp3x3','l_volpreserving3x3','l_volumepreserving3x3','cl_vp9','cl_volpreserving9','cl_volumepreserving9','cl_vp3x3','cl_volpreserving3x3','cl_volumepreserving3x3'}
            if numel(p) ~= 8, error('8 parameters expected'); end

            L = zeros(4,4); L([1 2 3 5 6 7 9 10]) = p; L(3,3) = -L(1,1) -L(2,2);
            t = expm( L ); t(4,:) = [0 0 0 1];

            if compute_dH  && ~isempty(dp)
              dt = d_expm( L )*[1,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0;0,0,1,0,0,0,0,0;0,0,0,0,0,0,0,0;0,0,0,1,0,0,0,0;0,0,0,0,1,0,0,0;0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,1;-1,0,0,0,-1,0,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0]*dp;
            end

            if strncmp( arg , 'c' , 1 ),  updateCenteredH( t , dt );
            else                       ,          updateH( t , dt );
            end
            TTYPE = updateTTYPE( 'V' , strncmp( arg , 'c' , 1 ) && ~CENTER_AT_ZERO );


          case {'l_affine','l_affine12','l_affine3x4','cl_affine','cl_affine12','cl_affine3x4'}
            L = [ reshape( p , 3,4 ) ; 0 0 0 0 ];
            t = expm( L ); t(4,:) = [0 0 0 1];

            if compute_dH  && ~isempty(dp)
              dt = d_expm( L )*[1,0,0,0,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,1,0,0,0,0,0,0,0,0;0,0,0,0,1,0,0,0,0,0,0,0;0,0,0,0,0,1,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,1,0,0,0,0,0;0,0,0,0,0,0,0,1,0,0,0,0;0,0,0,0,0,0,0,0,1,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,0,0,0,0,1;0,0,0,0,0,0,0,0,0,0,0,0]*dp;
            end

            if strncmp( arg , 'c' , 1 ),  updateCenteredH( t , dt );
            else                       ,          updateH( t , dt );
            end
            TTYPE = updateTTYPE( 'A' , true );


          case {'l_affine9','l_affine3x3','cl_affine9','cl_affine3x3'}
            if numel(p) ~= 9, error('9 parameters expected'); end

            L = zeros(4,4); L(1:3,1:3) = reshape( p , [3,3] );
            t = expm( L ); t(4,:) = [0 0 0 1];

            if compute_dH  && ~isempty(dp)
              dt = d_expm( L )*[1,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,1,0,0,0,0,0;0,0,0,0,1,0,0,0,0;0,0,0,0,0,1,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,0,1;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0]*dp;
            end

            if strncmp( arg , 'c' , 1 ),  updateCenteredH( t , dt );
            else                       ,          updateH( t , dt );
            end
            TTYPE = updateTTYPE( 'A' , strncmp( arg , 'c' , 1 ) && ~CENTER_AT_ZERO );


          case {'l_sxyzttt','l_sxyzt','cl_sxyzttt','cl_sxyzt'}
            if numel(p) ~= 7, error('7 parameters expected'); end

            L = tangentmatrix( p(2) , p(3) , p(4) );
            L([ 1  6 11]) = p(1);
            L([13 14 15]) = p(5:7);
            t = expm( L ); t(4,:) = [0 0 0 1];

            if compute_dH  && ~isempty(dp)
              dt = d_expm( L )*[1,0,0,0,0,0,0;0,0,0,1,0,0,0;0,0,-1,0,0,0,0;0 0 0 0 0 0 0;0,0,0,-1,0,0,0;1,0,0,0,0,0,0;0,1,0,0,0,0,0;0 0 0 0 0 0 0;0,0,1,0,0,0,0;0,-1,0,0,0,0,0;1,0,0,0,0,0,0;0 0 0 0 0 0 0;0,0,0,0,1,0,0;0,0,0,0,0,1,0;0,0,0,0,0,0,1;0 0 0 0 0 0 0]*dp;
            end

            if strncmp( arg , 'c' , 1 ),  updateCenteredH( t , dt );
            else                       ,          updateH( t , dt );
            end
            TTYPE = updateTTYPE( 'M' , true );


          case {'l_xyzsttt','l_xyzst','cl_xyzsttt','cl_xyzst'}
            if numel(p) ~= 7, error('7 parameters expected'); end

            L = tangentmatrix( p(1) , p(2) , p(3) );
            L([ 1  6 11]) = p(4);
            L([13 14 15]) = p(5:7);
            t = expm( L ); t(4,:) = [0 0 0 1];

            if compute_dH  && ~isempty(dp)
              dt = d_expm( L )*[0,0,0,1,0,0,0;0,0,1,0,0,0,0;0,-1,0,0,0,0,0;0 0 0 0 0 0 0;0,0,-1,0,0,0,0;0,0,0,1,0,0,0;1,0,0,0,0,0,0;0 0 0 0 0 0 0;0,1,0,0,0,0,0;-1,0,0,0,0,0,0;0,0,0,1,0,0,0;0 0 0 0 0 0 0;0,0,0,0,1,0,0;0,0,0,0,0,1,0;0,0,0,0,0,0,1;0 0 0 0 0 0 0]*dp;
            end

            if strncmp( arg , 'c' , 1 ),  updateCenteredH( t , dt );
            else                       ,          updateH( t , dt );
            end
            TTYPE = updateTTYPE( 'M' , true );


          case {'l_xyzttt','l_xyzt','cl_xyzttt','cl_xyzt'}
            if numel(p) ~= 6, error('6 parameters expected'); end

            L = tangentmatrix( p(1) , p(2) , p(3) );
            L([13 14 15]) = p(4:6);
            t = expm( L ); t(4,:) = [0 0 0 1];

            if compute_dH  && ~isempty(dp)
              dt = d_expm( L )*[0 0 0 0 0 0;0,0,1,0,0,0;0,-1,0,0,0,0;0 0 0 0 0 0;0,0,-1,0,0,0;0 0 0 0 0 0;1,0,0,0,0,0;0 0 0 0 0 0;0,1,0,0,0,0;-1,0,0,0,0,0;0 0 0 0 0 0;0 0 0 0 0 0;0,0,0,1,0,0;0,0,0,0,1,0;0,0,0,0,0,1;0 0 0 0 0 0]*dp;
            end

            if strncmp( arg , 'c' , 1 ),  updateCenteredH( t , dt );
            else                       ,          updateH( t , dt );
            end
            TTYPE = updateTTYPE( 'R' , true );

            
          case {'l_xyz','cl_xyz'}
            if numel(p) ~= 3, error('3 parameters expected'); end

            L = tangentmatrix( p(1) , p(2) , p(3) );
            t = expm( L ); t(4,:) = [0 0 0 1];

            if compute_dH  && ~isempty(dp)
              dt = d_expm( L )*[0 0 0;0,0,1;0,-1,0;0 0 0;0,0,-1;0 0 0;1,0,0;0 0 0;0,1,0;-1,0,0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0]*dp;
            end

            if strncmp( arg , 'c' , 1 ),  updateCenteredH( t , dt );
            else                       ,          updateH( t , dt );
            end
            TTYPE = updateTTYPE( 'R' , strncmp( arg , 'c' , 1 ) && ~CENTER_AT_ZERO );


          case {'l_s','l_sss','cl_s','cl_sss'}
            if numel(p)==1
              UNIFORM = 'U';

              t = diag( [ exp([p;p;p]) ; 1  ] );
              if compute_dH  && ~isempty(dp)
                dt = [ t(1:15).' ; 0 ]*dp;
              end

            elseif numel(p) == 3
              UNIFORM = 'S';

              t = diag( [ exp(p(:)) ; 1 ] );
              if compute_dH  && ~isempty(dp)
                dt = [t(1),0,0;0 0 0;0 0 0;0 0 0;0 0 0;0,t(6),0;0 0 0;0 0 0;0 0 0;0 0 0;0,0,t(11);0 0 0;0 0 0;0 0 0;0 0 0;0 0 0]*dp;
              end

            else, error('1 or 3 parameters expected'); end

            if strncmp( arg , 'c' , 1 ),  updateCenteredH( t , dt );
            else                       ,          updateH( t , dt );
            end
            TTYPE = updateTTYPE( UNIFORM , strncmp( arg , 'c' , 1 ) && ~CENTER_AT_ZERO );


          case {'l_ttt','l_t','cl_ttt','cl_t'}
            if numel(p) ~= 3, error('3 parameters expected'); end

            L = zeros(4,4);
            L([13 14 15]) = p;
            t = expm( L ); t(4,:) = [0 0 0 1];

            if compute_dH  && ~isempty(dp)
              dt = d_expm( L )*[0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;1,0,0;0,1,0;0,0,1;0 0 0]*dp;
            end

            updateH( t , dt );
            TTYPE = updateTTYPE( '' , true );


          case {'l_sttt','l_st','cl_sttt','cl_st'}
            if numel(p) ~= 4, error('4 parameters expected'); end

            L = zeros(4,4);
            L([ 1  6 11]) = p(1);
            L([13 14 15]) = p(2:4);
            t = expm( L ); t(4,:) = [0 0 0 1];

            if compute_dH  && ~isempty(dp)
              dt = d_expm( L )*[1,0,0,0;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;1,0,0,0;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;1,0,0,0;0 0 0 0;0,1,0,0;0,0,1,0;0,0,0,1;0 0 0 0]*dp;
            end

            if strncmp( arg , 'c' , 1 ),  updateCenteredH( t , dt );
            else                       ,          updateH( t , dt );
            end
            TTYPE = updateTTYPE( 'S' , true );


          case {'l_sxyz','cl_sxyz'}
            if numel(p) ~= 4, error('4 parameters expected'); end

            L = tangentmatrix( p(2) , p(3) , p(4) );
            L([ 1  6 11 ]) = p(1);
            t = expm( L ); t(4,:) = [0 0 0 1];

            if compute_dH  && ~isempty(dp)
              dt = d_expm( L )*[1,0,0,0;0,0,0,1;0,0,-1,0;0,0,0,0;0,0,0,-1;1,0,0,0;0,1,0,0;0,0,0,0;0,0,1,0;0,-1,0,0;1,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0]*dp;
            end

            if strncmp( arg , 'c' , 1 ),  updateCenteredH( t , dt );
            else                       ,          updateH( t , dt );
            end
            TTYPE = updateTTYPE( 'M' , strncmp( arg , 'c' , 1 ) && ~CENTER_AT_ZERO );


          case {'l_xyzs','cl_xyzs'}
            if numel(p) ~= 4, error('4 parameters expected'); end

            L = tangentmatrix( p(1) , p(2) , p(3) );
            L([ 1  6 11 ]) = p(4);
            t = expm( L ); t(4,:) = [0 0 0 1];

            if compute_dH  && ~isempty(dp)
              dt = d_expm( L )*[0,0,0,1;0,0,1,0;0,-1,0,0;0,0,0,0;0,0,-1,0;0,0,0,1;1,0,0,0;0,0,0,0;0,1,0,0;-1,0,0,0;0,0,0,1;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0]*dp;
            end

            if strncmp( arg , 'c' , 1 ),  updateCenteredH( t , dt );
            else                       ,          updateH( t , dt );
            end
            TTYPE = updateTTYPE( 'M' , strncmp( arg , 'c' , 1 ) && ~CENTER_AT_ZERO );


          case {'q','quat','quaternion'}
            if numel(p) ~= 3, error('3 parameters expected'); end

            if any( abs(p) > 1 )
              H(1) = 0;

              rearrangeOutputs();

              return;  %error('parameters must be lower than 1');
            end

            if compute_dH,    [t,dQ] = quat2mat( p );
            else         ,     t     = quat2mat( p );
            end
            t(4,4) = 1;

            if compute_dH  && ~isempty(dp)
              dt = zeros(16,3);
              dt([1,2,3,5,6,7,9,10,11],:) = dQ;
              dt = dt*dp;
            end

            updateH( t , dt );
            TTYPE = updateTTYPE( 'R' );


          case {'cq','cquat','cquaternion','rq','rquat','rquaternion','relq','relquat','relquaternion','relativeq','relativequat','relativequaternion'}
            if numel(p) ~= 3, error('3 parameters expected'); end

            if any( abs(p) > 1 )
              H(1) = 0;

              rearrangeOutputs();

              return;  %error('parameters must be lower than 1');
            end

            if compute_dH,    [t,dQ] = quat2mat( p );
            else         ,     t     = quat2mat( p );
            end
            t(4,4) = 1;


            if compute_dH  && ~isempty(dp)
              dt = zeros(16,3);
              dt([1,2,3,5,6,7,9,10,11],:) = dQ;
              dt = dt*dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'R' , ~CENTER_AT_ZERO );


          case {'scale','s'}
            if numel(p) == 1
              UNIFORM = 'I';

              t = diag([p p p 1]);

              if compute_dH  && ~isempty(dp)
                dt = [1;0;0;0;0;1;0;0;0;0;1;0;0;0;0;0]*dp;
              end

            elseif numel(p) == 2
              UNIFORM = 'F';

              t = diag([p(1) p(2) 1 1]);
              if compute_dH  && ~isempty(dp)
                dt = [1,0;0 0;0 0;0 0;0 0;0,1;0 0;0 0;0 0;0 0;0 0;0 0;0 0;0 0;0 0;0 0]*dp;
              end

            elseif numel(p) == 3
              UNIFORM = 'F';

              t = diag([p(1) p(2) p(3) 1]);
              if compute_dH  && ~isempty(dp)
                dt = [1,0,0;0 0 0;0 0 0;0 0 0;0 0 0;0,1,0;0 0 0;0 0 0;0 0 0;0 0 0;0,0,1;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0]*dp;
              end

            else, error('1 or 3 parameters expected'); end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( UNIFORM , ~CENTER_AT_ZERO );
  

          case {'scale2','s2'}
            if numel(p) ~= 1, error('1 parameters expected'); end

            t = diag([p p 1 1]);
            if compute_dH  && ~isempty(dp)
              dt = [1;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0]*dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'F' , ~CENTER_AT_ZERO );

           
          case {'sx','xs','xscale','scalex'}
            if numel(p) ~= 1, error('1 parameters expected'); end

            t = diag([p 1 1 1]);

            if compute_dH  && ~isempty(dp)
              dt = [1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]*dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'F' , ~CENTER_AT_ZERO );


          case {'sy','ys','yscale','scaley'}
            if numel(p) ~= 1, error('1 parameters expected'); end

            t = diag([1 p 1 1]);

            if compute_dH  && ~isempty(dp)
              dt = [0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0]*dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'F' , ~CENTER_AT_ZERO );


          case {'sz','zs','zscale','scalez'}
            if numel(p) ~= 1, error('1 parameters expected'); end

            t = diag([1 1 p 1]);

            if compute_dH  && ~isempty(dp)
              dt = [0;0;0;0;0;0;0;0;0;0;1;0;0;0;0;0]*dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'F' , ~CENTER_AT_ZERO );


          case {'rz' 'zr' 'rotatez' 'zrotate' }
            if numel(p) ~= 1, error('1 parameters expected'); end

            t = RZ( p*runits );

            if compute_dH  && ~isempty(dp)
              dt = runits*dRZ( p*runits )*dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'R' , ~CENTER_AT_ZERO );


          case {'rx' 'xr' 'rotatex' 'xrotate' }
            if numel(p) ~= 1, error('1 parameters expected'); end

            t = RX( p*runits );

            if compute_dH  && ~isempty(dp)
              dt = runits*dRX( p*runits )*dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'R' , ~CENTER_AT_ZERO );


          case {'ry' 'yr' 'rotatey' 'yrotate' }
            if numel(p) ~= 1, error('1 parameters expected'); end

            t = RY( p*runits );

            if compute_dH  && ~isempty(dp)
              dt = runits*dRY( p*runits )*dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'R' , ~CENTER_AT_ZERO );


          case {'axisangle' 'rotateaxis' 'raxis' 'aa'}
            if numel(p) ~= 3, error('3 parameters expected'); end
            
            N = p;
            if norm(N), N= N/norm(N);
            else,       N= [0 ; 0 ; 1];    
            end

            a  = varargin{i}; i=i+1;
            a  = a*runits;
            ca = cos(a);   sa= sin(a);
            t  = [  ca+(1-ca)*N(1)*N(1)      (1-ca)*N(1)*N(2)-sa*N(3)  (1-ca)*N(1)*N(3)+sa*N(2)   0  ; ...
                   (1-ca)*N(1)*N(2)+sa*N(3)  ca+(1-ca)*N(2)*N(2)       (1-ca)*N(2)*N(3)-sa*N(1)   0  ; ...
                   (1-ca)*N(1)*N(3)-sa*N(2)  (1-ca)*N(2)*N(3)+sa*N(1)  ca+(1-ca)*N(3)*N(3)        0  ; ...
                    0                         0                         0                         1  ];

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'R' , ~CENTER_AT_ZERO );
              

          case {'rxyz'}

            if numel(p) ~= 3, error('3 parameters expected'); end

            t = RZ( p(3)*runits )*RY( p(2)*runits )*RX( p(1)*runits );

            if compute_dH  && ~isempty(dp)
              dt      = kron(                                    eye(4) , RZ( p(3)*runits )*RY( p(2)*runits ) ) * dRX( p(1)*runits );
              dt(:,2) = kron(                       RX( p(1)*runits ).' , RZ( p(3)*runits )                   ) * dRY( p(2)*runits );
              dt(:,3) = kron( RX( p(1)*runits ).' * RY( p(2)*runits ).' , eye(4)                              ) * dRZ( p(3)*runits );
              dt      = runits * dt * dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'R' , ~CENTER_AT_ZERO );
          case {'rxzy'}
            if numel(p) ~= 3, error('3 parameters expected'); end

            t = RY( p(3)*runits )*RZ( p(2)*runits )*RX( p(1)*runits );

            if compute_dH  && ~isempty(dp)
              dt      = kron(                                    eye(4) , RY( p(3)*runits )*RZ( p(2)*runits ) ) * dRX( p(1)*runits );
              dt(:,2) = kron(                       RX( p(1)*runits ).' , RY( p(3)*runits )                   ) * dRZ( p(2)*runits );
              dt(:,3) = kron( RX( p(1)*runits ).' * RZ( p(2)*runits ).' , eye(4)                              ) * dRY( p(3)*runits );
              dt      = runits * dt * dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'R' , ~CENTER_AT_ZERO );
          case {'ryxz'}
            if numel(p) ~= 3, error('3 parameters expected'); end

            t = RZ( p(3)*runits )*RX( p(2)*runits )*RY( p(1)*runits );

            if compute_dH  && ~isempty(dp)
              dt      = kron(                                    eye(4) , RZ( p(3)*runits )*RX( p(2)*runits ) ) * dRY( p(1)*runits );
              dt(:,2) = kron(                       RY( p(1)*runits ).' , RZ( p(3)*runits )                   ) * dRX( p(2)*runits );
              dt(:,3) = kron( RY( p(1)*runits ).' * RX( p(2)*runits ).' , eye(4)                              ) * dRZ( p(3)*runits );
              dt      = runits * dt * dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'R' , ~CENTER_AT_ZERO );
          case {'ryzx'}
            if numel(p) ~= 3, error('3 parameters expected'); end

            t = RX( p(3)*runits )*RZ( p(2)*runits )*RY( p(1)*runits );

            if compute_dH  && ~isempty(dp)
              dt      = kron(                                    eye(4) , RX( p(3)*runits )*RZ( p(2)*runits ) ) * dRY( p(1)*runits );
              dt(:,2) = kron(                       RY( p(1)*runits ).' , RX( p(3)*runits )                   ) * dRZ( p(2)*runits );
              dt(:,3) = kron( RY( p(1)*runits ).' * RZ( p(2)*runits ).' , eye(4)                              ) * dRX( p(3)*runits );
              dt      = runits * dt * dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'R' , ~CENTER_AT_ZERO );
          case {'rzyx'}
            if numel(p) ~= 3, error('3 parameters expected'); end

            t = RX( p(3)*runits )*RY( p(2)*runits )*RZ( p(1)*runits );

            if compute_dH  && ~isempty(dp)
              dt      = kron(                                    eye(4) , RX( p(3)*runits )*RY( p(2)*runits ) ) * dRZ( p(1)*runits );
              dt(:,2) = kron(                       RZ( p(1)*runits ).' , RX( p(3)*runits )                   ) * dRY( p(2)*runits );
              dt(:,3) = kron( RZ( p(1)*runits ).' * RY( p(2)*runits ).' , eye(4)                              ) * dRX( p(3)*runits );
              dt      = runits * dt * dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'R' , ~CENTER_AT_ZERO );
          case {'rzxy'}
            if numel(p) ~= 3, error('3 parameters expected'); end

            t = RY( p(3)*runits )*RX( p(2)*runits )*RZ( p(1)*runits );

            if compute_dH  && ~isempty(dp)
              dt      = kron(                                    eye(4) , RY( p(3)*runits )*RX( p(2)*runits ) ) * dRZ( p(1)*runits );
              dt(:,2) = kron(                       RZ( p(1)*runits ).' , RY( p(3)*runits )                   ) * dRX( p(2)*runits );
              dt(:,3) = kron( RZ( p(1)*runits ).' * RX( p(2)*runits ).' , eye(4)                              ) * dRY( p(3)*runits );
              dt      = runits * dt * dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'R' , ~CENTER_AT_ZERO );

          case {'rxyx'}
            if numel(p) ~= 3, error('3 parameters expected'); end

            t = RX( p(3)*runits )*RY( p(2)*runits )*RX( p(1)*runits );

            if compute_dH  && ~isempty(dp)
              dt      = kron(                                    eye(4) , RX( p(3)*runits )*RY( p(2)*runits ) ) * dRX( p(1)*runits );
              dt(:,2) = kron(                       RX( p(1)*runits ).' , RX( p(3)*runits )                   ) * dRY( p(2)*runits );
              dt(:,3) = kron( RX( p(1)*runits ).' * RY( p(2)*runits ).' , eye(4)                              ) * dRX( p(3)*runits );
              dt      = runits * dt * dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'R' , ~CENTER_AT_ZERO );
          case {'rxzx'}
            if numel(p) ~= 3, error('3 parameters expected'); end

            t = RX( p(3)*runits )*RZ( p(2)*runits )*RX( p(1)*runits );

            if compute_dH  && ~isempty(dp)
              dt      = kron(                                    eye(4) , RX( p(3)*runits )*RZ( p(2)*runits ) ) * dRX( p(1)*runits );
              dt(:,2) = kron(                       RX( p(1)*runits ).' , RX( p(3)*runits )                   ) * dRZ( p(2)*runits );
              dt(:,3) = kron( RX( p(1)*runits ).' * RZ( p(2)*runits ).' , eye(4)                              ) * dRX( p(3)*runits );
              dt      = runits * dt * dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'R' , ~CENTER_AT_ZERO );
          case {'ryxy'}
            if numel(p) ~= 3, error('3 parameters expected'); end

            t = RY( p(3)*runits )*RX( p(2)*runits )*RY( p(1)*runits );

            if compute_dH  && ~isempty(dp)
              dt      = kron(                                    eye(4) , RY( p(3)*runits )*RX( p(2)*runits ) ) * dRY( p(1)*runits );
              dt(:,2) = kron(                       RY( p(1)*runits ).' , RY( p(3)*runits )                   ) * dRX( p(2)*runits );
              dt(:,3) = kron( RY( p(1)*runits ).' * RX( p(2)*runits ).' , eye(4)                              ) * dRY( p(3)*runits );
              dt      = runits * dt * dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'R' , ~CENTER_AT_ZERO );
          case {'ryzy'}
            if numel(p) ~= 3, error('3 parameters expected'); end

            t = RY( p(3)*runits )*RZ( p(2)*runits )*RY( p(1)*runits );

            if compute_dH  && ~isempty(dp)
              dt      = kron(                                    eye(4) , RY( p(3)*runits )*RZ( p(2)*runits ) ) * dRY( p(1)*runits );
              dt(:,2) = kron(                       RY( p(1)*runits ).' , RY( p(3)*runits )                   ) * dRZ( p(2)*runits );
              dt(:,3) = kron( RY( p(1)*runits ).' * RZ( p(2)*runits ).' , eye(4)                              ) * dRY( p(3)*runits );
              dt      = runits * dt * dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'R' , ~CENTER_AT_ZERO );
          case {'rzxz'}
            if numel(p) ~= 3, error('3 parameters expected'); end

            t = RZ( p(3)*runits )*RX( p(2)*runits )*RZ( p(1)*runits );

            if compute_dH  && ~isempty(dp)
              dt      = kron(                                    eye(4) , RZ( p(3)*runits )*RX( p(2)*runits ) ) * dRZ( p(1)*runits );
              dt(:,2) = kron(                       RZ( p(1)*runits ).' , RZ( p(3)*runits )                   ) * dRX( p(2)*runits );
              dt(:,3) = kron( RZ( p(1)*runits ).' * RX( p(2)*runits ).' , eye(4)                              ) * dRZ( p(3)*runits );
              dt      = runits * dt * dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'R' , ~CENTER_AT_ZERO );
          case {'rzyz'}
            if numel(p) ~= 3, error('3 parameters expected'); end

            t = RZ( p(3)*runits )*RY( p(2)*runits )*RZ( p(1)*runits );

            if compute_dH  && ~isempty(dp)
              dt      = kron(                                    eye(4) , RZ( p(3)*runits )*RY( p(2)*runits ) ) * dRZ( p(1)*runits );
              dt(:,2) = kron(                       RZ( p(1)*runits ).' , RZ( p(3)*runits )                   ) * dRY( p(2)*runits );
              dt(:,3) = kron( RZ( p(1)*runits ).' * RY( p(2)*runits ).' , eye(4)                              ) * dRZ( p(3)*runits );
              dt      = runits * dt * dp;
            end

            updateCenteredH( t , dt );
            TTYPE = updateTTYPE( 'R' , ~CENTER_AT_ZERO );


          case {'t' 'translate' 'move' }
            if numel(p) ~= 3, error('3 parameters expected'); end

            t([13 14 15]) = p;

            if compute_dH  && ~isempty(dp)
              dt = [0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;1,0,0;0,1,0;0,0,1;0 0 0]*dp;
            end

            updateH( t , dt );
            TTYPE = updateTTYPE( '' , true );


          case {'t2','tt' }
            if numel(p) ~= 2, error('2 parameters expected'); end

            if isa(p,'sym'), t=sym(t); end
            t([13 14]) = p;

            if compute_dH  && ~isempty(dp)
              dt = [0 0;0 0;0 0;0 0;0 0;0 0;0 0;0 0;0 0;0 0;0 0;0 0;1,0;0,1;0,0;0 0]*dp;
            end

            updateH( t , dt );
            TTYPE = updateTTYPE( '' , true );


          case {'tx' 'xtranslate' 'xmove'  'translatex' 'movex'}
            if numel(p) ~= 1, error('1 parameters expected'); end

            t( 13 ) = p;

            if compute_dH  && ~isempty(dp)
              dt = [0;0;0;0;0;0;0;0;0;0;0;0;1;0;0;0]*dp;
            end

            updateH( t , dt );
            TTYPE = updateTTYPE( '' , true );


          case {'ty' 'ytranslate' 'ymove'  'translatey' 'movey'}
            if numel(p) ~= 1, error('1 parameters expected'); end

            t( 14 ) = p;

            if compute_dH  && ~isempty(dp)
              dt = [0;0;0;0;0;0;0;0;0;0;0;0;0;1;0;0]*dp;
            end

            updateH( t , dt );
            TTYPE = updateTTYPE( '' , true );


          case {'tz' 'ztranslate' 'zmove'  'translatez' 'movez'}
            if numel(p) ~= 1, error('1 parameters expected'); end

            t( 15 ) = p;

            if compute_dH  && ~isempty(dp)
              dt = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;1;0]*dp;
            end

            updateH( t , dt );
            TTYPE = updateTTYPE( '' , true );


          otherwise
            error('incorrect transform type   %s' , uneval(arg) );
        end
    end
  end

  rearrangeOutputs();
  
  
  if RETURN_TTYPE
    H = TTYPE;
  end
  
  
  function t= RX(t)
    c = cos(t); s = sin(t);
    t = [ 1 0 0 0; 0 c -s 0; 0 s c 0; 0 0 0 1];
  end
  function t= RY(t)
    c = cos(t); s = sin(t);
    t = [ c 0 s 0; 0 1 0 0; -s 0 c 0; 0 0 0 1];
  end
  function t= RZ(t)
    c = cos(t); s = sin(t);
    t = [ c -s 0 0; s c 0 0; 0 0 1 0; 0 0 0 1];
  end

  function t= dRX(t)
    c = cos(t); s = sin(t);
    t = [ 0 0 0 0; 0 -s -c 0; 0 c -s 0; 0 0 0 0];
    t = t(:);
  end
  function t= dRY(t)
    c = cos(t); s = sin(t);
    t = [ -s 0 c 0; 0 0 0 0; -c 0 -s 0; 0 0 0 0];
    t = t(:);
  end
  function t= dRZ(t)
    c = cos(t); s = sin(t);
    t = [ -s -c 0 0; c -s 0 0; 0 0 0 0; 0 0 0 0];
    t = t(:);
  end


  function t = tangentmatrix( x , y , z )
    t = [ 0 -z y 0; z 0 -x 0; -y x 0 0; 0 0 0 0 ];
  end


  function nTTYPE = updateTTYPE( m , t )
    %     types:
    %     E : identity
    %     L : inversion without scale change
    %     U : uniform scale  ( positive , no inversions )
    %     I : uniform scale ( with or without inversions )
    %     S : non_uniform scale ( all positives , no inversions )
    %     F : non uniform scale ( with or without inversions ), allowing flips
    %     R : rotation  ( no inversions )
    %     O : orthogonal ( rotation with or without inversions )
    %     M : similarity ( rotation and uniform scale , no inversions )
    %     N : similarity ( rotation + uniform scale , with or without inversions )
    %     V : volume preserving affine ( no inversions )
    %     W : volume preserving affine ( with or, inversions )
    %     A : affine ( no inversions )
    %     G : general linear ( affine + inversions )
    %     P : perspective ...
    
    
    if nargin < 2, t = false; end

    M = TTYPE(1);
    switch m
      case ''
      case 'L'
        switch M, case 'E', M = 'L';
                  case 'L', M = 'L';
                  case 'U', M = 'I';
                  case 'I', M = 'I';
                  case 'S', M = 'F';
                  case 'F', M = 'F';
                  case 'R', M = 'O';
                  case 'O', M = 'O';
                  case 'M', M = 'N';
                  case 'N', M = 'N';
                  case 'V', M = 'W';
                  case 'W', M = 'W';
                  case 'A', M = 'G';
                  case 'G', M = 'G';
                  case 'P', M = 'P'; end
      case 'U'
        switch M, case 'E', M = 'U';
                  case 'L', M = 'I';
                  case 'U', M = 'U';
                  case 'I', M = 'I';
                  case 'S', M = 'S';
                  case 'F', M = 'F';
                  case 'R', M = 'M';
                  case 'O', M = 'N';
                  case 'M', M = 'M';
                  case 'N', M = 'N';
                  case 'V', M = 'A';
                  case 'W', M = 'G';
                  case 'A', M = 'A';
                  case 'G', M = 'G';
                  case 'P', M = 'P'; end
      case 'I'
        switch M, case 'E', M = 'I';
                  case 'L', M = 'I';
                  case 'U', M = 'I';
                  case 'I', M = 'I';
                  case 'S', M = 'F';
                  case 'F', M = 'F';
                  case 'R', M = 'N';
                  case 'O', M = 'N';
                  case 'M', M = 'N';
                  case 'N', M = 'N';
                  case 'V', M = 'G';
                  case 'W', M = 'G';
                  case 'A', M = 'G';
                  case 'G', M = 'G';
                  case 'P', M = 'P'; end
      case 'S'
        switch M, case 'E', M = 'S';
                  case 'L', M = 'F';
                  case 'U', M = 'S';
                  case 'I', M = 'F';
                  case 'S', M = 'S';
                  case 'F', M = 'F';
                  case 'R', M = 'A';
                  case 'O', M = 'G';
                  case 'M', M = 'A';
                  case 'N', M = 'G';
                  case 'V', M = 'A';
                  case 'W', M = 'G';
                  case 'A', M = 'A';
                  case 'G', M = 'G';
                  case 'P', M = 'P'; end
      case 'F'
        switch M, case 'E', M = 'F';
                  case 'L', M = 'F';
                  case 'U', M = 'F';
                  case 'I', M = 'F';
                  case 'S', M = 'F';
                  case 'F', M = 'F';
                  case 'R', M = 'G';
                  case 'O', M = 'G';
                  case 'M', M = 'G';
                  case 'N', M = 'G';
                  case 'V', M = 'G';
                  case 'W', M = 'G';
                  case 'A', M = 'G';
                  case 'G', M = 'G';
                  case 'P', M = 'P'; end
      case 'R'
        switch M, case 'E', M = 'R';
                  case 'L', M = 'O';
                  case 'U', M = 'M';
                  case 'I', M = 'N';
                  case 'S', M = 'A';
                  case 'F', M = 'G';
                  case 'R', M = 'R';
                  case 'O', M = 'O';
                  case 'M', M = 'M';
                  case 'N', M = 'N';
                  case 'V', M = 'V';
                  case 'W', M = 'W';
                  case 'A', M = 'A';
                  case 'G', M = 'G';
                  case 'P', M = 'P'; end
      case 'O'
        switch M, case 'E', M = 'O';
                  case 'L', M = 'O';
                  case 'U', M = 'N';
                  case 'I', M = 'N';
                  case 'S', M = 'G';
                  case 'F', M = 'G';
                  case 'R', M = 'O';
                  case 'O', M = 'O';
                  case 'M', M = 'N';
                  case 'N', M = 'N';
                  case 'V', M = 'W';
                  case 'W', M = 'W';
                  case 'A', M = 'G';
                  case 'G', M = 'G';
                  case 'P', M = 'P'; end
      case 'M'
        switch M, case 'E', M = 'M';
                  case 'L', M = 'N';
                  case 'U', M = 'M';
                  case 'I', M = 'N';
                  case 'S', M = 'A';
                  case 'F', M = 'G';
                  case 'R', M = 'M';
                  case 'O', M = 'N';
                  case 'M', M = 'M';
                  case 'N', M = 'N';
                  case 'V', M = 'A';
                  case 'W', M = 'G';
                  case 'A', M = 'A';
                  case 'G', M = 'G';
                  case 'P', M = 'P'; end
      case 'N'
        switch M, case 'E', M = 'N';
                  case 'L', M = 'N';
                  case 'U', M = 'N';
                  case 'I', M = 'N';
                  case 'S', M = 'G';
                  case 'F', M = 'G';
                  case 'R', M = 'N';
                  case 'O', M = 'N';
                  case 'M', M = 'N';
                  case 'N', M = 'N';
                  case 'V', M = 'G';
                  case 'W', M = 'G';
                  case 'A', M = 'G';
                  case 'G', M = 'G';
                  case 'P', M = 'P'; end
      case 'V'
        switch M, case 'E', M = 'V';
                  case 'L', M = 'W';
                  case 'U', M = 'A';
                  case 'I', M = 'G';
                  case 'S', M = 'A';
                  case 'F', M = 'G';
                  case 'R', M = 'V';
                  case 'O', M = 'G';
                  case 'M', M = 'A';
                  case 'N', M = 'G';
                  case 'V', M = 'V';
                  case 'W', M = 'W';
                  case 'A', M = 'A';
                  case 'G', M = 'G';
                  case 'P', M = 'P'; end
      case 'W'
        switch M, case 'E', M = 'W';
                  case 'L', M = 'W';
                  case 'U', M = 'G';
                  case 'I', M = 'G';
                  case 'S', M = 'G';
                  case 'F', M = 'G';
                  case 'R', M = 'W';
                  case 'O', M = 'W';
                  case 'M', M = 'G';
                  case 'N', M = 'G';
                  case 'V', M = 'W';
                  case 'W', M = 'W';
                  case 'A', M = 'G';
                  case 'G', M = 'G';
                  case 'P', M = 'P'; end
      case 'A'
        switch M, case 'E', M = 'A';
                  case 'L', M = 'G';
                  case 'U', M = 'A';
                  case 'I', M = 'G';
                  case 'S', M = 'A';
                  case 'F', M = 'G';
                  case 'R', M = 'A';
                  case 'O', M = 'G';
                  case 'M', M = 'A';
                  case 'N', M = 'G';
                  case 'V', M = 'A';
                  case 'W', M = 'G';
                  case 'A', M = 'A';
                  case 'G', M = 'G';
                  case 'P', M = 'P'; end
      case 'G'
        switch M, case 'E', M = 'G';
                  case 'L', M = 'G';
                  case 'U', M = 'G';
                  case 'I', M = 'G';
                  case 'S', M = 'G';
                  case 'F', M = 'G';
                  case 'R', M = 'G';
                  case 'O', M = 'G';
                  case 'M', M = 'G';
                  case 'N', M = 'G';
                  case 'V', M = 'G';
                  case 'W', M = 'G';
                  case 'A', M = 'G';
                  case 'G', M = 'G';
                  case 'P', M = 'P'; end
      case 'P'
        switch M, case 'E', M = 'P';
                  case 'L', M = 'P';
                  case 'U', M = 'P';
                  case 'I', M = 'P';
                  case 'S', M = 'P';
                  case 'F', M = 'P';
                  case 'R', M = 'P';
                  case 'O', M = 'P';
                  case 'M', M = 'P';
                  case 'N', M = 'P';
                  case 'V', M = 'P';
                  case 'W', M = 'P';
                  case 'A', M = 'P';
                  case 'G', M = 'P';
                  case 'P', M = 'P'; end
    end

    if numel( TTYPE ) > 1
      T = TTYPE(2);
    else
      T = '';
    end
    
    if t, T = 'T'; end
    
    nTTYPE = [M T];
  end


  function updateH( t , dt )
    try, if isa(  t ,'sym'),  t = simple(  t ); end; end
    try, if isa( dt ,'sym'), dt = simple( dt ); end; end
    try, if isa(  H ,'sym'),  H = simple(  H ); end; end
    try, if isa( dH ,'sym'), dH = simple( dH ); end; end
    
    if isnumeric(t), t = double( t ); end
    if compute_dH
      if isnumeric(dt), dt = double( dt ); end
      if isempty(dH)
        dH = kron( H.' , eye(4) )*dt;
      else
        dH = [  kron( eye(4) , t )*dH  ,  kron( H.' , eye(4) )*dt ];
      end
    end
    if ~isequal(t,eye(4)), H = t*H; end
  end
  function updateCenteredH( t , dt )
    try, if isa(  t ,'sym'), try,  t = simplify(  t ); catch,  t = simple(  t ); end; end; end
    try, if isa( dt ,'sym'), try, dt = simplify( dt ); catch, dt = simple( dt ); end; end; end
    try, if isa(  H ,'sym'), try,  H = simplify(  H ); catch,  H = simple(  H ); end; end; end
    try, if isa( dH ,'sym'), try, dH = simplify( dH ); catch, dH = simple( dH ); end; end; end

    if isnumeric(t), t = double( t ); end
    if compute_dH
      if isnumeric(dt), dt = double( dt ); end
      if isempty(dH)
        dH = kron( ICENTER*H , CENTER.' ).'*dt;
      else
        dH = [ kron( eye(4) , CENTER*t*ICENTER )*dH   ,  kron( ICENTER*H , CENTER.' ).'*dt ];
      end
      if ~isempty(CENTER_idx)
        dH(:,CENTER_idx) = dH(:,CENTER_idx) + kron( H.' , eye(4) ) * (  kron(    eye(4) , CENTER*t )   * dIC + kron( t*ICENTER , eye(4)   ).' * dC  );
      end
    end
    if ~isnumeric(t) || ~isequal(t,eye(4))
      H = CENTER*t*ICENTER*H;
    end
  end

  function rearrangeOutputs
    if      ~IN_HOMOGENEOUS_COORDINATES  &&  ~IN_2D
      H = H(1:3,1:3);
      if compute_dH, dH = dH([1 2 3 5 6 7 9 10 11],:); end
    elseif  ~IN_HOMOGENEOUS_COORDINATES  &&   IN_2D
      H = H(1:2,1:2);
      if compute_dH, dH = dH([1 2 5 6],:); end
    elseif   IN_HOMOGENEOUS_COORDINATES  &&   IN_2D
      H = H([1 2 4],[1 2 4]);
      if compute_dH, dH = dH([1 2 4 5 6 8 13 14 16],:); end
    elseif   IN_HOMOGENEOUS_COORDINATES  &&  ~IN_2D
    end
  end

end
