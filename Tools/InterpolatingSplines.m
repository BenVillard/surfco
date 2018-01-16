function [uW,IS] = InterpolatingSplines( M , F , W , varargin )
% 
% 
%                                        L
%                                       ---  ||                || 2
% Solves:    min  reg( u )  + LAMBDA * \     ||  u( M ) -  F   ||
%             u                        /     ||      l      l  ||
%                                       ---  ||                ||
%                                       l=1
%
% and return u(W)
%
% Instead of specify reg(.) you have to specify the green function of reg.
%
%
%

if 0
  
figure
% clf
G = [         ndmat( [ linspace(-1,2,100) , NaN ] , -1:0.25:2 )   ;...
      fliplr( ndmat( [ linspace(-1,2,100) , NaN ] , -1:0.25:2 ) ) ];

F = rand(5,2);
F = ndmat( linspace(0,1,2) , linspace(0,1,3) );
% F = ndmat( linspace(0,1,2) , 0 );

getXY = @(h)[ vec(get(h,'XData')) , vec(get(h,'YData')) ];
setXY = @(h,xy)set_(h,'XData',xy(:,1),'YData',xy(:,2));

hF = line( F(:,1) , F(:,2) ,'color',[0 0 0],'linestyle','none','marker','s','markerfacecolor','r','markersize',9,'hittest','off');
hM = arrayfun( @(m)line( F(m,1) , F(m,2) ,'color',[0 0 0],'linestyle','none','marker','o','markerfacecolor','b'),1:size(F,1)).';
hG = line( G(:,1) , G(:,2) ,'color',[1 1 1]*0.7,'hittest','off');
hL = line( NaN , NaN , 'color','k','marker','.','hittest','off');
LAMBDA = eEntry( 'range', [-3 5],'ivalue', 5 ,'step', 0.1 ,'normal','position',[2 , 2 , 0 , 0] );

set(gca,'DataAspectRatio',[1 1 1],'XLim',[-1.5 2.5],'YLim',[-1.5 2.5]);
drawnow;


update_      = @(d){ setXY( hG , d(G) ) , setXY( hL , reshape(permute(cat(3, F , d(F) , NaN(size(F)) ),[3 1 2]),[],2) ) };
update_plus_ = @(d){ setXY( hG , G+d(G) ) , setXY( hL , reshape(permute(cat(3, F , F+d(F) , NaN(size(F)) ),[3 1 2]),[],2) ) };

TPS_fcn = {'rlogr'};
TPS_fcn = {'r'};
TPS_fcn = { @(x,y)-1/( 8*pi )*ipd(x,y) , @(x)[ones(size(x,1),1)] };

update = @()update_plus_( InterpolatingSplines0( getXY(hF) , -getXY(hF) + cell2mat(arrayfun(@(h)getXY(h),hM,'un',0)) ,[] , TPS_fcn{:} ,'LAMBDA',pow10(LAMBDA.v)) );

for m = 1:size(F,1),dragableobject( hM(m) , @(h)update() );end
LAMBDA.continuous = true; LAMBDA.callback_fcn = @(x)update();

%%
end
if 0
from = [zeros(1,2);1,0;0,0.5;1,0.5;0,1;ones(1,2)];
to   = [1.705882352941176,-0.22408963585434138;2.526610644257703,-0.21288515406162434;-0.1859943977591037,0.5784313725490198;2.6610644257703084,0.8921568627450975;.5042016806722687,1.1008403361344539;3.042016806722689,1.1008403361344534];
  
  
G = [         ndmat( [ linspace(-1,2,501) , NaN ] , -1:0.1:2 )   ;...
      fliplr( ndmat( [ linspace(-1,2,501) , NaN ] , -1:0.1:2 ) ) ];
  
LAMBDA = bestLAMBDA;

u = InterpolatingSplines( from , to ,[] ,'r','LAMBDA',pow10(LAMBDA) )


 plot3d( u(G) , '-' ,'Color',[1 1 1]*0.8,'eq');
text( from(:,1)+0.1 , from(:,2) , arrayfun( @num2str,1:size(from,1),'un',0) ,'HorizontalAlignment','l','VerticalAlignment','m','FontSize',12,'BackgroundColor','none','EdgeColor','none','Margin',1 )
hplot3d( from , 'sb','MarkerFaceColor','b');
hplot3d( to   , 'or','MarkerFaceColor','r');
hplot3d( from , to , '.:k');
hplot3d( from , u(from) , '.-m','LineWidth',2,'Marker','o','MarkerSize',6,'MarkerFaceColor','m');


d0 = fro( from-to , 2 );
d1 = fro( u(from)-to , 2 );
[ d0 , d1 , d0 < d1 , d0 - d1 ]

fro2( min( d0 - d1 ,0 ) )

%%

figure;
discrepENER = @(L)( min( d0 - fro( InterpolatingSplines( from , to ,from ,'r','LAMBDA',pow10(L) ) -to ,2) ) );
plotfun( @(L)discrepENER(L) , linspace( -2 , 4 ,101) ,'.-' );hline(0);

bestLAMBDA = fzero( @(L)discrepENER(L) , 0 )


%%
end




  if numel( varargin ) == 0, error('G_fcn and g_fcn should be provided. try with ''rlogr'' or ''r''.'); end
  if ischar( varargin{1} )
    switch lower( varargin{1} )
      case {'rlogr','r2logr'}
        rbf   = @(r) 1/( 8*pi ) * nonans(noinfs( r.^2.*log(r) ,0),0);                       %the Radial Basis Function
        G_fcn = @(x,y) rbf( ipd(x,y) );                                                     %Green's function of the bihamonic
        g_fcn = @(x) [ x , ones(size(x,1),1) ];                                             %null space of the bihamonic
        %g_fcn  = @(x) [];
        
      case 'r'
        G_fcn = @(x,y) -1/( 8*pi ) * nonans(noinfs( ipd(x,y) ,0),0);
        g_fcn = @(x) [ x , ones(size(x,1),1) ];
        
      case 'r3'
        G_fcn = @(x,y) -1/( 8*pi ) * nonans(noinfs( ipd(x,y).^3 ,0),0);
        g_fcn = @(x) [ x , ones(size(x,1),1) ];
        
      otherwise
        error('unknown Green''s function');
        
    end
    varargin(1) = [];
  elseif isa( varargin{1} , 'function_handle' )
    if numel( varargin ) < 2 || ~isa( varargin{2} , 'function_handle' )
      error('Null space function (g_fcn) must be provided.');
    end
    G_fcn = varargin{1};
    g_fcn = varargin{2};
    varargin(1:2) = [];
  end

  [varargin,i,LAMBDA] = parseargs(varargin, 'lambda','$DEFS$',Inf);
  if ~isempty(LAMBDA) && ~isscalar(LAMBDA), error('a scalar is expected in LAMBDA.'); end
  if ~isempty(LAMBDA) && LAMBDA <= 0,       error('LAMBDA must be greater than 0.');  end

  L   = size(M,1);
  NSD = size(M,2);
  nC  = size(F,2);

  %using MASK
  if ~isempty(F) && size(F,1) ~= L, error('Number of points M and fixed values F must coincide.'); end

  [varargin,i,MASK] = parseargs(varargin, 'mask','$DEFS$',[]);
  if isscalar( MASK ) && MASK >= L
    MASK = [];
  elseif isscalar( MASK )
    MASK = unique( round( linspace( 1 , L , MASK ) ) );
  end
  try
    if ~isempty( MASK )
                      M = M( MASK ,1:NSD);
      if ~isempty(F), F = F( MASK ,1:nC ); end
    end
  catch
    error('some error evaluating MASK');
  end
  L = size(M,1);
  %%%%
  

  if 0
  elseif   isempty(LAMBDA)  &&  ~isempty( F )
    error('not implemented yet');
  elseif   isempty(LAMBDA)  &&   isempty( F )
    error('not implemented yet');
  elseif  ~isinf( LAMBDA )  &&  ~isempty( F )
  
  elseif  ~isinf( LAMBDA )  &&   isempty( F )
    

  elseif   isinf( LAMBDA )

  else
    error('invalid calling');
  end
  

  if 0
  elseif  ( ~isempty( F )  &&  ~isempty(LAMBDA)  &&   isinf(LAMBDA) ) ||...
          (  isempty( F )  &&  ~isempty(LAMBDA)  &&   isinf(LAMBDA) )
    G = G_fcn( M , M );  %size should be: L x L
    if max(abs(vec( G - G.'))) > 1e-10
      G_is_sym = false;
      warning('G isn''t symmetric  %g' ,  max(abs(vec( G - G.')))  );
    else, G_is_sym = true;
    end

    g = g_fcn( M );      %size should be: L x K
    K = size( g , 2 );
    if ~K, g = zeros( L , 0 ); end
    
    if ~isempty( F )
      if G_is_sym
        A = [ G   , g           ;...
              g.' , zeros(K,K) ];
        B = [ F            ;...
              zeros(K,nC) ];
      else
        A = [ (G+G.')/2  , zeros(L,K) , G.'         ;...
              zeros(K,L) , zeros(K,K) , g.'         ;...
              G          , g          , zeros(L,L) ];
        B = [ zeros( L , nC )  ;...
              zeros( K , nC )  ;...
              F               ];
      end
      ba = ppinv( A , B );
      beta   = ba(     1:L     ,1:nC);
      alpha  = ba( (L+1):(L+K) ,1:nC);
    end
    
    if nargout > 1 || isempty( W ) || isempty( F )
      inv_G = ppinv(G);
      
      if G_is_sym, S = inv_G;
      else,        S = double( inv_G ); S = ( S + S.' )/2;
      end
      
%       if iscolumn(g) && all( g == 0 )
%         ALPHA = zeros(1,NSD);
%         BETA  = inv_G;
%       else
        ALPHA = ppinv( g.' * S * g  ,  g.' * S );
        BETA  = inv_G * ( eye(L,L) - g * ALPHA );
%       end
      
    end
    if nargout > 1
      IS.F      = F;
      IS.LAMBDA = LAMBDA;
      IS.G      = G;
      IS.ALPHA  = ALPHA;
      IS.BETA   = BETA;

      R         = BETA.' * G * BETA;
      IS.R      = R;

      if isempty(F)
        IS.alpha  = @(F) ALPHA*F;
        IS.beta   = @(F) BETA*F;
        IS.u      = @(x,F) ( g_fcn(x)*ALPHA + G_fcn(x,M)*BETA )*F;
        IS.reg_c  = @(F) diag( F.' * R * F ).';
        IS.reg    = @(F) trace(F.' * R * F);
      else
        IS.alpha  = alpha;
        IS.beta   = beta;
        IS.u      = @(x) g_fcn(x)*alpha + G_fcn(x,M)*beta;
        IS.reg_c  = arrayfun( @(c) beta(:,c).' * IS.G * beta(:,c) , 1:nC );
        IS.reg    = sum( IS.reg_c );
      end
    end
    
  elseif  ( ~isempty( F )  &&  ~isempty(LAMBDA)  &&  ~isinf(LAMBDA) ) ||...
          (  isempty( F )  &&  ~isempty(LAMBDA)  &&  ~isinf(LAMBDA) )
    [~,IS] = InterpolatingSplines( M , [] , [] , G_fcn , g_fcn , 'lambda' , Inf );
    if ~isempty( F )
      %consider also to solve as the solution of the problem
      %[ G + eye(size(G))/LAMBDA , g ; g.' , zeros(K,K) ] * [ beta ; alpha ] = [ F ; zeros(K,nC) ]
      %that is:
      % ba = [ G + eye(size(G))/LAMBDA , g ; g.' , zeros(K,K) ] \ [ F ; zeros(K,nC) ]
      X = LAMBDA * ppinv( IS.R + LAMBDA*eye(L,L) , F );
      beta  = IS.BETA  * X;
      alpha = IS.ALPHA * X;
    end
    
    if nargout > 1 || isempty( W ) || isempty( F )
      Xopt = LAMBDA * ppinv( IS.R + LAMBDA*eye(L,L) );
      ALPHA = IS.ALPHA * Xopt;
      BETA  = IS.BETA  * Xopt;
    end    
    if nargout > 1
      IS.F      = F;
      IS.LAMBDA = LAMBDA;
      
      IS.ALPHA = ALPHA;
      IS.BETA  = BETA;

      Rinf = IS.R;
      R = Xopt.' *Rinf * Xopt;
      IS.R     = R;

      IS.Xopt = Xopt;
      
      if isempty(F)
        IS.alpha  = @(F) ALPHA*F;
        IS.beta   = @(F) BETA*F;
        IS.u      = @(x,F) ( g_fcn(x)*ALPHA + G_fcn(x,M)*BETA )*F;
        IS.reg_c  = @(F) diag( F.' * R * F ).';
        IS.reg    = @(F) trace(F.' * R * F);
        %IS.X      = @(F) Xopt*F;
        IS.X      = @(F) LAMBDA * ( ( Rinf + LAMBDA*eye(L,L) ) \ F );
        IS.ener   = @(F) trace(F.' * ( R + LAMBDA * ( eye(L,L) - Xopt ).' * ( eye(L,L) - Xopt ) ) * F);
      else
        IS.alpha  = alpha;
        IS.beta   = beta;
        IS.u      = @(x) g_fcn(x)*alpha + G_fcn(x,M)*beta;
        IS.reg_c  = arrayfun( @(c) beta(:,c).' * IS.G * beta(:,c) , 1:nC );
        IS.reg    = sum( IS.reg_c );
        IS.X      = X;
        IS.ener   = IS.reg + LAMBDA * fro2( F - X );
      end
    end
    
  elseif  ( ~isempty( F )  &&   isempty(LAMBDA)                     )
    error('not implemented yet');
  elseif  (  isempty( F )  &&   isempty(LAMBDA)                     )
    error('not implemented yet');
  else,   error('invalid calling');
  end
  
  %now evaluate u in W
  if 0
  elseif  ( ~isempty( W )  &&  ~isempty( F )  &&  ~isempty(LAMBDA)  &&   isinf(LAMBDA) ) ||...
          ( ~isempty( W )  &&  ~isempty( F )  &&  ~isempty(LAMBDA)  &&  ~isinf(LAMBDA) )
    nW = size( W , 1 ); uW = NaN( [ nW , nC ] );

    ev = nW; while ev >= 1; try
      for i = 1:ev:nW
        if i > 1, fprintf('%9d of %d\n',i,nW); end
        w = ( i ):min( i+ev-1 , nW ); XYZ = W( w ,1:NSD);
        uW(w, 1:nC ) = G_fcn( XYZ , M ) * beta + g_fcn( XYZ ) * alpha;
      end; break; end;
      ev = round( ev / 2.001 );
    end
    
  elseif  ( ~isempty( W )  &&   isempty( F )  &&  ~isempty(LAMBDA)  &&   isinf(LAMBDA) ) ||...
          ( ~isempty( W )  &&   isempty( F )  &&  ~isempty(LAMBDA)  &&  ~isinf(LAMBDA) )
    nW = size( W , 1 ); linearOp = NaN( [ nW , L ] );

    ev = nW; while ev >= 1; try
      for i = 1:ev:nW
        if i > 1, fprintf('%9d of %d\n',i,nW); end
        w = ( i ):min( i+ev-1 , nW ); XYZ = W( w ,1:NSD);
        linearOp(w, 1:L ) = G_fcn( XYZ , M ) * BETA + g_fcn( XYZ ) * ALPHA;
      end; break; end;
      ev = round( ev / 2.001 );
    end
    uW = @(F) linearOp * F;

  elseif  (  isempty( W )  &&  ~isempty( F )  &&  ~isempty(LAMBDA)  &&   isinf(LAMBDA) ) ||...
          (  isempty( W )  &&  ~isempty( F )  &&  ~isempty(LAMBDA)  &&  ~isinf(LAMBDA) )
    uW = @(x) g_fcn(x)*alpha + G_fcn(x,M)*beta;
        
  elseif  (  isempty( W )  &&   isempty( F )  &&  ~isempty(LAMBDA)  &&   isinf(LAMBDA) ) ||...
          (  isempty( W )  &&   isempty( F )  &&  ~isempty(LAMBDA)  &&  ~isinf(LAMBDA) )
    uW = @(x,F) (g_fcn(x)*ALPHA + G_fcn(x,M)*BETA)*F;

  elseif  ( ~isempty( W )  &&  ~isempty( F )  &&   isempty(LAMBDA)                     )
    error('to do...');
  elseif  ( ~isempty( W )  &&   isempty( F )  &&   isempty(LAMBDA)                     )
    error('to do...');
  elseif  (  isempty( W )  &&  ~isempty( F )  &&   isempty(LAMBDA)                     )
    error('to do...');
  elseif  (  isempty( W )  &&   isempty( F )  &&   isempty(LAMBDA)                     )
    error('to do...');
  end
  
  if isa( uW , 'function_handle'), try, uW = CleanFH( uW ); end; end
  
  if nargout > 1
    for f = fieldnames(IS).'
      if isa( IS.(f{1}) , 'function_handle' )
        try, IS.(f{1}) = CleanFH( IS.(f{1}) ); end
      end
    end
  end

end
