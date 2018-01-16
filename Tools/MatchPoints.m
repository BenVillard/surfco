function [T , T_M , err ] = MatchPoints( F , M , type )
%{

N = 3;
norm2=@(x) x(:).'*x(:);
F = [0 2   ;0   0   ;1   0   ];
M = [0 2.2 ;1.1 0.1 ;0.1 0.1 ];
F=bsxfun(@minus,F,mean(F,1))
N=bsxfun(@minus,M,mean(M,1))

H = ( F.' * M )/( M.' * M )
det(ans)


F = [0 2   1 ;0   0   1 ;1   0   1 ; 0 0 0];
M = [0 1.2 1 ;1.1 0.1 1 ;0.1 0.1 1 ; 0 0 0];


%}

  if iscell( M )
    T   = cell( size(M) );
    T_M = cell( size(M) );
    err = cell( size(M) );
    for m = 1:numel(M)
      [ T{m} , T_M{m} , err{m} ] = MatchPoints( F , M{m} , type );
    end
    return;
  end


  if ~isequal( size(F) , size(M) ), error('M and F have to be equal size'); end
  NSD = size( F , 2 );
  if NSD ~= 3  &&  NSD ~= 2       , error(' Nx3   or Nx2  matrices expected.'); end
  N = size(F,1);


  if isa( type , 'function_handle' )

    if NSD == 2
      M(end,3) = 0;
      F(end,3) = 0;
    end
      
    p = invmaketransform( [] , type );

    sM = sum( M , 1 );
    sF = sum( F , 1 );

    MM =   [  M.' * M  ,  sM.' ;  sM , N ];   %MM  -   Mh.'*Mh
    kMMI    =  kron( MM , eye( 4 ) );

%     MF = 2*[  M.' * F  ,  sM.' ;  sF , N ];   %MF  - 2*Mh.'*Fh
%     mf  = vect( MF.' );

    FM = 2*[  F.' * M  ,  sF.' ;  sM , N ];   %MF  - 2*Mh.'*Fh
    fm  = vect( FM );

%     FF =   [  F.' * F  ,  sF.' ;  sF , N ];   %FF  -   Fh.'*Fh
%     tFF = trace(FF);

    FF =   F.' * F;
    tFF = trace(FF) + N;


    p = Optimize( @(p) ener(p) , p , 'methods' ,{'quasinewton',250,'conjugate',50,'descendneg',1,'coordinate',1},'ls',{'quadratic','golden','quadratic'},'verbose',0,struct('MIN_ENERGY',0),'noplot');
    %disp( ener(p) );

    T = maketransform( type , p );

    if NSD == 2, T = T([1 2 4],[1 2 4]); end
    return;
  end

  if NSD == 3
    T = eye( NSD+1 ); delements = [1 6 11]; telements = [13 14 15];
  else
    T = eye( NSD+1 ); delements = [1 5   ]; telements = [ 7  8   ];
  end
  switch sprintf('%s_%d', upper(type) , NSD )
    %  *   E : identity
    %     L : inversion without scale change
    %  *   U : uniform scale  ( positive , no inversions )
    %  *   I : uniform scale ( with or without inversions )
    %  *   S : non_uniform scale ( all positives , no inversions )
    %  *   F : non uniform scale ( with or without inversions ), allowing flips
    %  *   R : rotation  ( no inversions )
    %  *   O : orthogonal ( rotation with or without inversions )
    %  *   M : similarity ( rotation and uniform scale , no inversions )
    %  *   N : similarity ( rotation + uniform scale , with or without inversions )
    %     V : volume preserving affine ( no inversions )
    %     W : volume preserving affine ( with or, inversions )
    %  +   A : affine ( no inversions )
    %  *   G : general linear ( affine + inversions )
    %     P : perspective ...
    case {'E_2','E_3'}                   % E : identity
      T = eye( NSD + 1 );
    case {'ET_3','T_3','ET_2','T_2'}     % E : identity
      T(telements) = ( sum( F , 1 ) - sum( M , 1 ) )/N;
    case 'L_2'                           % L : inversion without scale change
      error('falta por resolver');
    case 'LT_2'                          % L : inversion without scale change
      error('falta por resolver');
    case 'L_3'                           % L : inversion without scale change
      error('falta por resolver');
    case 'LT_3'                          % L : inversion without scale change
      error('falta por resolver');
    case {'I_3','I_2'}                   % I : uniform scale ( with or without inversions )
      T(delements) = trace( F.' * M )/trace( M.' * M );
    case {'IT_3','IT_2'}                 % I : uniform scale ( with or without inversions )
      Fc = bsxfun( @minus , F , mean(F,1) );
      Mc = bsxfun( @minus , M , mean(M,1) );
      
      s = trace( Fc.' * Mc )/trace( Mc.' * Mc );
      T(delements) = s;
      
      T(telements) = ( sum( F , 1 ) - s * sum( M , 1 ) )/N;
    case {'U_3','U_2'}                   % U : uniform scale  ( positive , no inversions )
      T(delements) = max( 0 , trace( F.' * M )/trace( M.' * M ) );
    case {'UT_3','UT_2'}                 % U : uniform scale  ( positive , no inversions )
      Fc = bsxfun( @minus , F , mean(F,1) );
      Mc = bsxfun( @minus , M , mean(M,1) );
      
      s = max( 0 , trace( Fc.' * Mc )/trace( Mc.' * Mc ) );
      T(delements) = s;
      
      T(telements) = ( sum( F , 1 ) - s * sum( M , 1 ) )/N;
    case 'F_3'                           % F : non uniform scale ( with or without inversions ), allowing flips
      MM = M.' * M;
      FM = F.' * M;
      
      T(delements) = [ FM(1,1)/MM(1,1) ; FM(2,2)/MM(2,2) ; FM(3,3)/MM(3,3) ];
    case 'F_2'                           % F : non uniform scale ( with or without inversions ), allowing flips
      MM = M.' * M;
      FM = F.' * M;
      
      T(delements) = [ FM(1,1)/MM(1,1) ; FM(2,2)/MM(2,2) ];
    case 'FT_3'                          % F : non uniform scale ( with or without inversions ), allowing flips
      Fc = bsxfun( @minus , F , mean(F,1) );
      Mc = bsxfun( @minus , M , mean(M,1) );
      
      MM = Mc.' * Mc;
      FM = Fc.' * Mc;
      
      T(delements) = [ FM(1,1)/MM(1,1) ; FM(2,2)/MM(2,2) ; FM(3,3)/MM(3,3) ];
      
      T(telements) = ( sum( F , 1 ) - sum( M * T(1:NSD,1:NSD).' , 1 ) )/N;
    case 'FT_2'                          % F : non uniform scale ( with or without inversions ), allowing flips
      Fc = bsxfun( @minus , F , mean(F,1) );
      Mc = bsxfun( @minus , M , mean(M,1) );
      
      MM = Mc.' * Mc;
      FM = Fc.' * Mc;
      
      T(delements) = [ FM(1,1)/MM(1,1) ; FM(2,2)/MM(2,2) ];
      
      T(telements) = ( sum( F , 1 ) - sum( M * T(1:NSD,1:NSD).' , 1 ) )/N;
    case 'S_3'                           % S : non_uniform scale ( all positives , no inversions )
      MM = M.' * M;
      FM = F.' * M;
      
      T(delements) = max( 0 , [ FM(1,1)/MM(1,1) ; FM(2,2)/MM(2,2) ; FM(3,3)/MM(3,3) ] );
    case 'S_2'                           % S : non_uniform scale ( all positives , no inversions )
      MM = M.' * M;
      FM = F.' * M;
      
      T(delements) = max( 0 , [ FM(1,1)/MM(1,1) ; FM(2,2)/MM(2,2) ] );
    case 'ST_3'                          % S : non_uniform scale ( all positives , no inversions )
      Fc = bsxfun( @minus , F , mean(F,1) );
      Mc = bsxfun( @minus , M , mean(M,1) );
      
      MM = Mc.' * Mc;
      FM = Fc.' * Mc;
      
      T(delements) = max( 0 , [ FM(1,1)/MM(1,1) ; FM(2,2)/MM(2,2) ; FM(3,3)/MM(3,3) ] );
      
      T(telements) = ( sum( F , 1 ) - sum( M * T(1:NSD,1:NSD).' , 1 ) )/N;
    case 'ST_2'                          % S : non_uniform scale ( all positives , no inversions )
      Fc = bsxfun( @minus , F , mean(F,1) );
      Mc = bsxfun( @minus , M , mean(M,1) );
      
      MM = Mc.' * Mc;
      FM = Fc.' * Mc;
      
      T(delements) = max( 0 , [ FM(1,1)/MM(1,1) ; FM(2,2)/MM(2,2) ] );
      
      T(telements) = ( sum( F , 1 ) - sum( M * T(1:NSD,1:NSD).' , 1 ) )/N;
    case {'O_3','O_2'}                   % O : orthogonal ( rotation with or without inversions )
      [U,S,V] = svd( M.' * F );
      
      T(1:NSD,1:NSD) = V*U.';
    case {'OT_3','OT_2'}                 % O : orthogonal ( rotation with or without inversions )
      Fc = bsxfun( @minus , F , mean(F,1) );
      Mc = bsxfun( @minus , M , mean(M,1) );
      
      [U,S,V] = svd( Mc.' * Fc );
      
      T(1:NSD,1:NSD) = V*U.';

      T(telements) = ( sum( F , 1 ) - sum( M * T(1:NSD,1:NSD).' , 1 ) )/N;
    case 'R_3'                           % R : rotation  ( no inversions )
      [U,S,V] = svd( M.' * F );
      
      T(1:NSD,1:NSD) = V*diag([1,1,det(V*U.')])*U.';
    case 'R_2'                           % R : rotation  ( no inversions )
      [U,S,V] = svd( M.' * F );
      
      T(1:NSD,1:NSD) = V*diag([1,det(V*U.')])*U.';
    case 'RT_3'                          % R : rotation  ( no inversions )
      Fc = bsxfun( @minus , F , mean(F,1) );
      Mc = bsxfun( @minus , M , mean(M,1) );
      
      [U,S,V] = svd( Mc.' * Fc );
      
      T(1:NSD,1:NSD) = V*diag([1,1,det(V*U.')])*U.';

      T(telements) = ( sum( F , 1 ) - sum( M * T(1:NSD,1:NSD).' , 1 ) )/N;
    case 'RT_2'                          % R : rotation  ( no inversions )
      Fc = bsxfun( @minus , F , mean(F,1) );
      Mc = bsxfun( @minus , M , mean(M,1) );
      
      [U,S,V] = svd( Mc.' * Fc );
      
      T(1:NSD,1:NSD) = V*diag([1,det(V*U.')])*U.';

      T(telements) = ( sum( F , 1 ) - sum( M * T(1:NSD,1:NSD).' , 1 ) )/N;
    case 'N_3'                           % N : similarity ( rotation + uniform scale , with or without inversions )
      [U,S,V] = svd( M.' * F );
      
      bestE = Inf;
      bestR = eye(NSD);
      trMM  = trace( M.' * M );
      trFF  = trace( F.' * F );
      MF    = M.' * F;
      for sign_change = 0:8
        switch sign_change
          case 0, d = [  0 ,  0 ,  0 ];
          case 1, d = [  1 ,  1 ,  1 ];
          case 2, d = [ -1 ,  1 ,  1 ];
          case 3, d = [  1 , -1 ,  1 ];
          case 4, d = [ -1 , -1 ,  1 ];
          case 5, d = [  1 ,  1 , -1 ];
          case 6, d = [ -1 ,  1 , -1 ];
          case 7, d = [  1 , -1 , -1 ];
          case 8, d = [ -1 , -1 , -1 ];
        end
        R = V * diag(d) * U.';
        s = trace( R * MF ) / trMM;
        R = s*R;
        
        %E = trace( ( M*R.' - F ).' * ( M*R.' - F ) );
        E = s*s*trMM - 2*trace( R * MF ) + trFF;
        if E < bestE, bestE = E;  bestR = R; end
        if bestE == 0, break; end
      end
      
      T(1:NSD,1:NSD) = bestR;
    case 'N_2'                           % N : similarity ( rotation + uniform scale , with or without inversions )
      [U,S,V] = svd( M.' * F );
      
      bestE = Inf;
      bestR = eye(NSD);
      trMM  = trace( M.' * M );
      trFF  = trace( F.' * F );
      MF    = M.' * F;
      for sign_change = 0:4
        switch sign_change
          case 0, d = [  0 ,  0 ];
          case 1, d = [  1 ,  1 ];
          case 2, d = [ -1 ,  1 ];
          case 3, d = [  1 , -1 ];
          case 4, d = [ -1 , -1 ];
        end
        R = V * diag(d) * U.';
        s = trace( R * MF ) / trMM;
        R = s*R;
        
        %E = trace( ( M*R.' - F ).' * ( M*R.' - F ) );
        E = s*s*trMM - 2*trace( R * MF ) + trFF;
        if E < bestE, bestE = E;  bestR = R; end
        if bestE == 0, break; end
      end
      
      T(1:NSD,1:NSD) = bestR;
    case {'NT_3','NT_2'}                 % N : similarity ( rotation + uniform scale , with or without inversions )
      Fc = bsxfun( @minus , F , mean(F,1) );
      Mc = bsxfun( @minus , M , mean(M,1) );
      
      T = MatchPoints( Fc , Mc , 'N' );
      
      T(telements) = ( sum( F , 1 ) - sum( M * T(1:NSD,1:NSD).' , 1 ) )/N;
    case 'M_3'                           % M : similarity ( rotation and uniform scale , no inversions )
      [U,S,V] = svd( M.' * F );
      
      bestE = Inf;
      bestR = eye(NSD);
      trMM  = trace( M.' * M );
      trFF  = trace( F.' * F );
      MF    = M.' * F;
      for sign_change = 0:8
        switch sign_change
          case 0, d = [  0 ,  0 ,  0 ];
          case 1, d = [  1 ,  1 ,  1 ];
          case 2, d = [ -1 ,  1 ,  1 ];
          case 3, d = [  1 , -1 ,  1 ];
          case 4, d = [ -1 , -1 ,  1 ];
          case 5, d = [  1 ,  1 , -1 ];
          case 6, d = [ -1 ,  1 , -1 ];
          case 7, d = [  1 , -1 , -1 ];
          case 8, d = [ -1 , -1 , -1 ];
        end
        
        R = V * diag(d) * U.';
        s = trace( R * MF ) / trMM;
        
        R = s*R;
        if det( R ) <= 0, continue; end
        
        %disp( [ sign_change , s , trace( ( M*R.' - F ).' * ( M*R.' - F ) ) , det( V * diag(d) * U.' ) ] );
        
        %E = trace( ( M*R.' - F ).' * ( M*R.' - F ) );
        E = s*s*trMM - 2*trace( R * MF ) + trFF;
        if E < bestE, bestE = E;  bestR = R; end
        if bestE == 0, break; end
      end
      
      T(1:NSD,1:NSD) = bestR;
    case 'M_2'                           % M : similarity ( rotation and uniform scale , no inversions )
      [U,S,V] = svd( M.' * F );
      
      bestE = Inf;
      bestR = eye(NSD);
      trMM  = trace( M.' * M );
      trFF  = trace( F.' * F );
      MF    = M.' * F;
      for sign_change = 0:4
        switch sign_change
          case 0, d = [  0 ,  0 ];
          case 1, d = [  1 ,  1 ];
          case 2, d = [ -1 ,  1 ];
          case 3, d = [  1 , -1 ];
          case 4, d = [ -1 , -1 ];
        end
        R = V * diag(d) * U.';
        s = trace( R * MF ) / trMM;
        R = s*R;
        
        %disp( [ sign_change , s , trace( ( M*R.' - F ).' * ( M*R.' - F ) ) , det( V * diag(d) * U.' ) ] );
        if det( R ) < 0, continue; end
        
        %E = trace( ( M*R.' - F ).' * ( M*R.' - F ) );
        E = s*s*trMM - 2*trace( R * MF ) + trFF;
        if E < bestE, bestE = E;  bestR = R; end
        if bestE == 0, break; end
      end
      
      T(1:NSD,1:NSD) = bestR;
    case {'MT_3','MT_2'}                 % M : similarity ( rotation and uniform scale , no inversions )
      Fc = bsxfun( @minus , F , mean(F,1) );
      Mc = bsxfun( @minus , M , mean(M,1) );
      
      T = MatchPoints( Fc , Mc , 'M' );
      
      T(telements) = ( sum( F , 1 ) - sum( M * T(1:NSD,1:NSD).' , 1 ) )/N;
    case {'G_3','G_2'}                   % G : general linear ( affine + inversions )
      T(1:NSD,1:NSD) = ( F.' * M ) / ( M.' * M );
    case {'GT_3','GT_2'}                 % G : general linear ( affine + inversions )
      Fc = bsxfun( @minus , F , mean(F,1) );
      Mc = bsxfun( @minus , M , mean(M,1) );
      
      T(1:NSD,1:NSD) = ( Fc.' * Mc ) * pinv( Mc.' * Mc );
      
      T(telements) = ( sum( F , 1 ) - sum( M * T(1:NSD,1:NSD).' , 1 ) )/N;
    case 'A_3'                           % A : affine ( no inversions )
      bestA = ( F.' * M ) / ( M.' * M );
      
      if det( bestA ) < 0
        warning( 'todavia no se bien como resolverlo!!!!' );
        norm2 = @(x) x(:).'*x(:);
        
        As = {};
        
        As{end+1} = MatchPoints(F,M, @(p) {'l_xyz',p(1:3),'s',{[0,exp(p(4)),exp(p(5))],[0,0;exp(p(4)),0;0,exp(p(5))]},'l_xyz',p(6:8)});
        %As{end+1} = MatchPoints(F,M, @(p) {'l_xyz',p(1:3),'s',{[0,p(4)^2,p(5)^2],2*[0,0;p(4),0;0,p(5)]},'l_xyz',p(6:8)});
        As{end+1} = MatchPoints(F,M, @(p) {'l_xyz',p(1:3),'s',{[exp(p(4)),0,exp(p(5))],[exp(p(4)),0;0,0;0,exp(p(5))]},'l_xyz',p(6:8)});
        As{end+1} = MatchPoints(F,M, @(p) {'l_xyz',p(1:3),'s',{[exp(p(4)),exp(p(5)),0],[exp(p(4)),0;0,exp(p(5));0,0]},'l_xyz',p(6:8)});
        
        As{end+1} = MatchPoints(F,M, @(p) {'l_xyz',p(1:3),'s',{[exp(p(4)),0,0],[exp(p(4));0;0]},'l_xyz',p(5:7)});
        As{end+1} = MatchPoints(F,M, @(p) {'l_xyz',p(1:3),'s',{[0,exp(p(4)),0],[0;exp(p(4));0]},'l_xyz',p(5:7)});
        As{end+1} = MatchPoints(F,M, @(p) {'l_xyz',p(1:3),'s',{[0,0,exp(p(4))],[0;0;exp(p(4))]},'l_xyz',p(5:7)});
        
        As{end+1} = zeros(NSD);
  
        %As{end+1} = MatchPoints(F,M, @(p) {'l_xyz',p(1:3),'s',{p(4:6).^2,2*diag(p(4:6))},'l_xyz',p(7:9)});
        As{end+1} = MatchPoints(F,M, @(p) {'l_xyz',p(1:3),'l_s',p(4:6),'l_xyz',p(7:9)});

        bestE = Inf;
        for i = 1:numel(As)
          E = norm2( M*As{i}(1:NSD,1:NSD).' - F );
          if E < bestE, bestE = E; bestA = As{i}(1:NSD,1:NSD); end
        end
      end

      T(1:NSD,1:NSD) = bestA;
    case 'A_2'                           % A : affine ( no inversions )
      bestA = ( F.' * M ) / ( M.' * M );
      
      if det( bestA ) < 0
        warning( 'todavia no se bien como resolverlo!!!!' );
        norm2 = @(x) x(:).'*x(:);
        
        As = {};
        
        As{end+1} = MatchPoints(F,M, @(p) {'rz',p(1),'s',{[exp(p(2)),0,1],[exp(p(2));0;0]},'rz',p(3)});
        As{end+1} = MatchPoints(F,M, @(p) {'rz',p(1),'s',{[0,exp(p(2)),1],[0;exp(p(2));0]},'rz',p(3)});
        
        As{end+1} = zeros(NSD);
  
        As{end+1} = MatchPoints(F,M, @(p) {'rz',p(1),'l_ss',p(2:3),'rz',p(4)});

        bestE = Inf;
        for i = 1:numel(As)
          E = norm2( M*As{i}(1:NSD,1:NSD).' - F );
          if E < bestE, bestE = E; bestA = As{i}(1:NSD,1:NSD); end
        end
      end

      T(1:NSD,1:NSD) = bestA;
    case {'AT_3','AT_2'}                 % A : affine ( no inversions )
      Fc = bsxfun( @minus , F , mean(F,1) );
      Mc = bsxfun( @minus , M , mean(M,1) );
      
      T = MatchPoints( Fc , Mc , 'A' );
      
      T(telements) = ( sum( F , 1 ) - sum( M * T(1:NSD,1:NSD).' , 1 ) )/N;
    case {'V_2','V_3'}                   % V : volume preserving affine ( no inversions )
      error('todavia no se como resolverlo!!!!');
    case {'VT_2','VT_3'}                 % V : volume preserving affine ( no inversions )
      error('todavia no se como resolverlo!!!!');
    case {'W_2','W_3'}                   % W : volume preserving affine ( with or, inversions )
      error('todavia no se como resolverlo!!!!');
    case {'WT_2','WT_3'}                 % W : volume preserving affine ( with or, inversions )
      error('todavia no se como resolverlo!!!!');
    case {'P_2','P_3'}                   % P : perspective ...
      error('todavia no se como resolverlo!!!!');
    case {'PT_2','PT_3'}                 % P : perspective ...
      error('todavia no se como resolverlo!!!!');
  end

  if nargout > 1
    T_M = transform( M , T );
  end
  if nargout > 2
    err = F - T_M;
    err = err(:).' * err(:);
  end
  
  function [E,dE] = ener( p )
    if nargout < 2
    
      H = maketransform( type , p );
      h = H(:);
      
      E = h.' * kMMI * h  -  fm * h + tFF;
      
    else
      
      [H,dH] = maketransform( type , p );
      
      h = H(:);

       E =   h.' * kMMI * h  -  fm * h + tFF;
      dE = ( 2 * h.' * kMMI  -  fm ) * dH;
      %NumericalDiff(@(p) ener(p) , p , 'i' ) - dE
    end
  end
  
end
