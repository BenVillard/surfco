function [BLID, ULID] = getLids ( C )

C(cellfun('isempty', C)) = []; 

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
  
  L = C(~w); 
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
  L = arrayfun(@(c) L{c}*R,1:length(L),'un',0)';
  
  if ~isempty( L )
      %%% Bottom Lid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
dd = []; for t = 1:numel(L), dd = vertcat(dd,L{t}(:,3)); end
          bLID = min(sort(dd));
          
          BLID = mean( C{ 1 } * R ,1); BLID(3) = bLID;
          BLID = BLID * R';
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
 

      %%% Upper Lid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    dd = []; for t = 1:numel(L), dd = vertcat(dd,L{t}(:,3)); end
          bLID = max(sort(dd));
          
          ULID = mean( C{ end } * R ,1); ULID(3) = bLID;
          ULID = ULID * R';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  elseif isempty( L )
           
      BLID = mean( C{ 1 } * R ,1); BLID = BLID * R';
      ULID = mean( C{ end } * R ,1); ULID = ULID * R';
  
  end
  
function Z = meanNormal( N )
  M    = @(ae) [ cos(ae(2)) * cos(ae(1)) , cos(ae(2)) * sin(ae(1)) , sin(ae(2)) ];
  dN2N = @(a,b) 2*min( asin(ipd(a,b)/2) , asin(ipd(a,-b)/2) );
  E    = @(m) sum( abs( dN2N( m , N ) ).^1 );
  
  m = mean(N,1);
  [ m(1),m(2),m(3) ] = cart2sph( m(1),m(2),m(3) );
  m = m([1,2]);
  try,
        m = Optimize( @(ae)E(M(ae)) , m , 'methods',{'conjugate','coordinate',1},'ls',{'quadratic','golden','quadratic'},struct('COMPUTE_NUMERICAL_JACOBIAN',{{'a'}}),'noplot','verbose',0);
  catch
        m = fminsearch(@(ae) E(M(ae)),m);
  end
  Z = M(m);
end









end