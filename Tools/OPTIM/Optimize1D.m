function [bestX,bestE,H,stored_values] = Optimize1D( f , X , P , varargin )
%
%
%
%
%
%



  if nargin > 2 && ischar( P )
    varargin = [ P , varargin ];
    P = [];
  end
  if nargin < 3 || isempty(P), P = struct; end
  P = setv( P , varargin{:} ); P = touppercasefields( P );
    
  %%defaults
  P = joinfields( P , 'METHODS'           , {'golden'}        );
  P = joinfields( P , 'VERBOSE'           , []                );   %espera una funcion que recibe:  verbose_fcn( level , str , varargin )
  P = joinfields( P , 'PLOT'              , false             );   %espera un handle donde plotear

  P = joinfields( P , 'FORCE_BRACKETING'  , false             );   %fueza el bracketing y usa el ultimo metodo especificado primero
  P = joinfields( P , 'EPS'               , 0                 );   %si  abs(Xi-Xk) < EPS,  Xi es lo mismo que Xj
  P = joinfields( P , 'TOL'               , 0                 );   %tolerancia de la smart function

  P = joinfields( P , 'MIN'               , minnum('single')  );   
  P = joinfields( P , 'MAX'               , maxnum('single')  );   
  

  P = joinfields( P , 'DO_BRACKETING_POS' , []                );   %hace bracketing hacia +Inf
  P = joinfields( P , 'DO_BRACKETING_NEG' , []                );   %hace bracketing hacia -Inf
  P = joinfields( P , 'BRACKETING_RHO'    , 2                 );   %factor de bracketing
  P = joinfields( P , 'BRACKETING_ITS'    , Inf               );   %max iters en bracketing
  P = joinfields( P , 'BRACKETING_INIT'   , eps(single(1))    );   %minimum step in bracketing
  P = joinfields( P , 'BRACKETING_STOP'   , maxnum('single')  );   %maximum step in bracketing
  P.TOL = max( P.TOL , P.EPS/10 );
  %%END defaults
  
  
  if P.PLOT
    if islogical( P.PLOT )  && P.PLOT
      P.PLOT = figure;
    end
    if ~ishandle( P.PLOT )
      error( 'PLOT has to be a handle!!!' );
    end  
    
    switch get( P.PLOT , 'type' )
      case 'figure'
        P.PLOT = axes('Parent',P.PLOT );
      case { 'axes','hggroup','hgtransform' }
        
      otherwise
        error( 'PLOT has to be a figure, axes, group handle ' );
    end
  end
  
  
  
  if ~iscell( P.METHODS ), P.METHODS = { P.METHODS }; end
  
  
  %%PARSE X
  %  []                 -> { {nan,nan,nan} , {nan,nan,nan} }
  %  Xi                 -> { {Xi ,nan,nan} , {nan,nan,nan} }
  %  [Xi Xs]            -> { {Xi,nan,nan} , {Xs,nan,nan} }
  %  [Xi Xs kk]         -> error
  %  {}                 -> { {nan,nan,nan} , {nan,nan,nan} }
  %  { Xi }             -> { {Xi,nan,nan} , {nan,nan,nan} }
  %  { Xi  Xs }         -> { {Xi,nan,nan} , {Xs,nan,nan} }
  %  { {Xi,Ei}   }      -> { {Xi,Ei,nan} , {nan,nan,nan} }
  %  { {Xi,Ei,Ji}   }   -> { {Xi,Ei,Ji} , {nan,nan,nan} }
  %  { [] , Xs   }      -> { {nan,nan,nan} , {Xs,nan,nan} }
  %  { [] , { Xs Es } } -> { {nan,nan,nan} , {Xs,Es,nan} }
  Xi = NaN;  Ei = NaN;  Ji = NaN;
  Xs = NaN;  Es = NaN;  Js = NaN;

  if iscell( X )
    if numel( X ) > 2, error('valor de X incorrecto'); end
    if numel( X ) > 0
      if iscell( X{1} )
        if numel( X{1} ) > 3,   error('valor de X incorrecto'); end
        if numel( X{1} ) > 0,   Xi = X{1}{1};   end
        if numel( X{1} ) > 1,   Ei = X{1}{2};   end
        if numel( X{1} ) > 2,   Ji = X{1}{3};   end
      else
        if numel( X{1} ) > 3,   error('valor de X incorrecto'); end
        if numel( X{1} ) > 0,   Xi = X{1}(1);   end
        if numel( X{1} ) > 1,   Ei = X{1}(2);   end
        if numel( X{1} ) > 2,   Ji = X{1}(3);   end
      end
    end
    if numel( X ) > 1
      if iscell( X{2} )
        if numel( X{2} ) > 3,   error('valor de X incorrecto'); end
        if numel( X{2} ) > 1,   Xs = X{2}{1};   end
        if numel( X{2} ) > 1,   Es = X{2}{2};   end
        if numel( X{2} ) > 2,   Js = X{2}{3};   end
      else
        if numel( X{2} ) > 3,   error('valor de X incorrecto'); end
        if numel( X{2} ) > 0,   Xs = X{2}(1);   end
        if numel( X{2} ) > 1,   Es = X{2}(2);   end
        if numel( X{2} ) > 2,   Js = X{2}(3);   end
      end
    end
  elseif ischar( X )
    error('X no puede ser char');
  else
    if numel( X ) > 0,    Xi = X(1);   end
    if numel( X ) > 1,    Xs = X(2);   end
    if numel( X ) > 2,    error('valor de X incorrecto'); end
  end
  
  if isnan( Xi )  &&  ( ~isnan( Ei )  ||  ~isnan( Ji ) ), error('valor de X incorrecto'); end
  if isnan( Xs )  &&  ( ~isnan( Es )  ||  ~isnan( Js ) ), error('valor de X incorrecto'); end
  if                     isnan( Ei )  &&  ~isnan( Ji )  , error('valor de X incorrecto'); end
  if                     isnan( Es )  &&  ~isnan( Js )  , error('valor de X incorrecto'); end
  if Xi > Xs                                            , error('valor de X incorrecto'); end
  %%END PARSE X

  
  try, feval( P.VERBOSE , 2 , '  Initializing Optimize1D:  { { %g , %g , %g }  ,  { %g , %g , %g } }\n', Xi , Ei , Ji , Xs , Es , Js ); end
  
  
  bestE = Inf;
  bestX = NaN;
  stored_values = zeros(100,3); last_stored = 0;
  if ~isnan( Ei )
    last_stored = last_stored+1;
    stored_values( last_stored , : ) = [ Xi Ei Ji];
    if Ei < bestE
      bestE = Ei;
      bestX = Xi;
    end    
  end
  if ~isnan( Es )
    last_stored = last_stored+1;
    stored_values( last_stored , : ) = [ Xs Es Ji]; 
    if Es < bestE
      bestE = Es;
      bestX = Xs;
    end    
  end
  
  
  hLINE_UPDATE = false;
  if P.PLOT
    hLINE = line('Parent',P.PLOT,'marker','x','color',[1 0 0],'linestyle','-','XData',[Xi Xs],'YData',[Ei Es]);
    hLINE_UPDATE = true;
  end
  
  function ee = F( x )
    a_idx = val2ind( stored_values( 1:last_stored , 1 ) , x );

    if ~isempty( a_idx )  &&  abs( stored_values(a_idx,1) - x ) <= P.TOL
      ee = stored_values(a_idx,2);
    else
      ee = f( x ); jj = NaN;
      
      last_stored = last_stored+1;
      if size( stored_values , 1 ) < last_stored
        stored_values( 2*last_stored , 2 ) = 0;
      end
      stored_values(last_stored,:) = [ x  ee  jj ];
      
      if ee < bestE
        bestE = ee;
        bestX = x;
        
%         if last_stored > 10
%           error('optimize1D:bestEnergyFound','bestEnergyFound');
%         end
        
      end
      
      if hLINE_UPDATE
        xs = stored_values( 1:last_stored ,1);
        [xs,ids] = sort( xs );
        es = stored_values(  ids          ,2);
        
        set( hLINE , 'XData',xs,'YData',es );
      end
      
    end
  end



%   try

    if       isnan(Xi)  && isnan(Xs)
      Xi = 0;
      Xs = 0;
      P.DO_BRACKETING_POS = true;
      P.DO_BRACKETING_NEG = true;
    elseif  ~isnan(Xi)  &&   isnan(Xs)
      Xs = Xi;
      P.DO_BRACKETING_POS = true;
    elseif   isnan(Xi)  &&  ~isnan(Xs)
      Xi = Xs;
      P.DO_BRACKETING_NEG = true;
    elseif isequal( Xi , Xs )
      P.DO_BRACKETING_POS = true;
      P.DO_BRACKETING_NEG = true;
    end

    if isnan( Ei ), Ei = F( Xi ); end
    if isnan( Es ), Es = F( Xs ); end


    H.PARAMETERS = P;
    H.initialX   = { {Xi,Ei,Ji},{Xs,Es,Js} };


    lsorder = 1:numel( P.METHODS );
    if P.FORCE_BRACKETING
      lsorder = lsorder([end 1:end-1]);
      if isempty( P.DO_BRACKETING_POS ),  P.DO_BRACKETING_POS = true; end
      if isempty( P.DO_BRACKETING_NEG ),  P.DO_BRACKETING_NEG = true; end
    end

    Xsoriginal = Xs;
    if ~isempty( P.DO_BRACKETING_POS )  && P.DO_BRACKETING_POS

      try, feval( P.VERBOSE , 2 , '    Performing BRACKETING POS:' ); end

      B  = Es;
      brack_it = 0;

      DX = max( [ P.EPS , P.BRACKETING_INIT , 10*eps(Xs) , abs(Xsoriginal-Xi) ] ) / P.BRACKETING_RHO;
      while 1
        DX = DX * P.BRACKETING_RHO;
        
        Bn = F( min( Xs + DX , P.MAX ) );
        brack_it = brack_it + 1;

        try, feval( P.VERBOSE , 2 , '\{    Performing BRACKETING POS:}r %g -> ( %g )', min( Xs + DX , P.MAX ) , Bn ); end
        
        if P.PLOT
          line('Parent',P.PLOT,'XData',min( Xs + DX , P.MAX ),'YData',Bn,'Marker','>','color',[1 1 1]*0.6);
        end
        
        

        if Bn       >   B                 ,  break;  end
        if DX       >=  P.BRACKETING_STOP ,  break;  end
        if Xs + DX  >=  P.MAX             ,  break;  end
        if brack_it >   P.BRACKETING_ITS  ,  break;  end

        B  = Bn;
      end
      Xs = min( Xs + DX , P.MAX );  Es = Bn;

      try, feval( P.VERBOSE , Inf ); end

    end


    if ~isempty( P.DO_BRACKETING_NEG )  && P.DO_BRACKETING_NEG

      try, feval( P.VERBOSE , 2 , '    Performing BRACKETING NEG:' ); end

      B  = Es;
      brack_it = 0;

      DX = max( [ P.EPS , P.BRACKETING_INIT , 10*eps(Xs) , abs(Xsoriginal-Xi) ] ) / P.BRACKETING_RHO;
      while 1
        DX = DX * P.BRACKETING_RHO;
        
        Bn = F( max( Xs - DX , P.MIN ) );
        brack_it = brack_it + 1;

        try, feval( P.VERBOSE , 2 , '\{    Performing BRACKETING MIN:}r %g -> ( %g )', max( Xs - DX , P.MIN ) , Bn ); end

        if P.PLOT
          line('Parent',P.PLOT,'XData',min( Xs + DX , P.MAX ),'YData',Bn,'Marker','<','color',[1 1 1]*0.6);
        end


        if Bn       >   B                 ,  break;  end
        if DX       >=  P.BRACKETING_STOP ,  break;  end
        if Xs - DX  >=  P.MIN             ,  break;  end
        if brack_it >   P.BRACKETING_ITS  ,  break;  end

        B  = Bn;
      end
      Xi = max( Xs - DX , P.MIN );  Ei = Bn;

      try, feval( P.VERBOSE , Inf ); end

    end

    

    H.afterBractX   = { {Xi,Ei,Ji},{Xs,Es,Js} };


    E0 = min( Ei,Es );

    %%do the Optimization
    sid = 0;
    for lsm = vect( uniquens( P.METHODS(lsorder) ) )
      sid = sid + 1;
      SID = sprintf('SEARCH_%02d',sid);
      H.(SID).METHOD = lsm{1};
      if isfield( P , upper( H.(SID).METHOD ) )
        PP = P.( upper( H.(SID).METHOD ) ); 
      else
        PP = struct;
      end
      if ~isfield( PP , 'VERBOSE'        ), PP.VERBOSE        = P.VERBOSE;        end
      if ~isfield( PP , 'PLOT'           ), PP.PLOT           = P.PLOT;           end
      if ~isfield( PP , 'EPS'            ), PP.EPS            = P.EPS;            end
      if ~isfield( PP , 'MIN'            ), PP.MIN            = P.MIN;            end
      if ~isfield( PP , 'MAX'            ), PP.MAX            = P.MAX;            end


      switch lower( H.(SID).METHOD )
        case 'golden'
          [X,E,HH] = Optimize1D_golden( @(x) F(x) , BETTERS(2) , PP );

        case 'exhaustive'
          [X,E,HH] = Optimize1D_exhaustive( @(x) F(x) , BETTERS(2) , PP );
          
        case 'smart'
          PP.getBETTERS = @() BETTERS( 0 );
          [X,E,HH] = Optimize1D_smart( @(x) F(x) , PP );
          


  %       case 'fixed'
  %         switch nargin
  %           case {0,1},  x = xf;       E = FLS( x );
  %           case 2    ,  x = a_fixed;  E = FLS( x );
  %           case 3    ,  x = a_fixed;  E = M_fixed;
  %         end
  %       case 'backtracking'
  %         [x,E] = backtracking( P , xf );
  %       case 'quadratic'
  %         [x,E] = quadratic( P , xf );

      end
      H.(SID) = joinfields( H.(SID) , HH );

      if E < E0
        try, feval( P.VERBOSE , 2 , '  LineSearch done with METHOD %s\n', H.(SID).METHOD ); end
        break;
      end
    end

    stored_values = stored_values( 1:last_stored ,:);
    
    plotFUN();

%   catch LE
% 
%
%     if isequal( LE.identifier , 'optimize1D:bestEnergyFound' )
%       try, feval( P.VERBOSE , Inf ); feval( P.VERBOSE , 2 , '  better energy found in %d evaluations\n' , last_stored ); end
%       return;
%     else
%       plotFUN();
%       rethrow( LE );
%     end
% 
%   end

  
  function plotFUN
    if ~P.PLOT, return; end
    
    hLINE_UPDATE = false;
    xs = stored_values( 1:last_stored , 1 );
    
    xs = [ min( xs ) , max( xs ) ];
    xs = ( xs-mean(xs) )*1.1 + mean( xs );
    
    xs = unique( [ linspace( xs(1) , xs(2) , 100 ) , stored_values( 1:last_stored , 1 ).' ] );
    
    es = zeros( numel(xs) , 1 );
    for i = 1:numel(xs)
      es(i) = F( xs(i) );
    end
    
    %hl = line( 'Parent',P.PLOT,'XData',xs,'YData',es,'color',[1 0 0],'linestyle','-','marker','none','linewidth',2 );
    set( hLINE , 'XData' , xs , 'YData', es ,'marker','none','linewidth',2 );
    hLINE_UPDATE = true;
    
    hb = line( 'Parent',P.PLOT,'XData',bestX,'YData',bestE,'color',[1 0 0],'linestyle','none','marker','o','markerfacecolor',[0 1 0],'markersize',10,'linewidth',2 );
  end

  
  function B = BETTERS( n )
    if last_stored < 2
      error('it only work with last_stored > 1');
    end
    
    
    xs = stored_values( 1:last_stored , 1 );
    [xs , ids ] = sort( xs );
    es = stored_values( ids , 2 );
    js = stored_values( ids , 3 );

      
      
    if numel( ids ) == n
      
      Bid = 1:n;
      
      
    elseif n == 0
      
      [be,Bid] = min( es );
      id = Bid(1);
      while id < numel( ids )  &&  es(id+1) >= es( id )
        id = id + 1;
        Bid = [ Bid , id ];
      end

      id = Bid(1);
      while id > 1  &&  es(id-1) >= es( id )
        id = id - 1;
        Bid = [ id , Bid ];
      end
      
    elseif n == 2
      
      [be,Bid] = min( es );
      if Bid == 1
        Bid = [1 2];
      elseif Bid == numel( ids )
        Bid = [ Bid-1 Bid ];
      else
        if es(Bid-1) < es(Bid+1)
          Bid = [ Bid-1  Bid   ];
        else
          Bid = [ Bid    Bid+1 ];
        end        
        
      end
        
    elseif n == 3  && numel( ids ) > 2
      
      [be,Bid] = min( es );
      if Bid == 1
        Bid = [1 2 3];
      elseif  Bid == numel(ids)
        Bid = [ Bid-2 Bid-1 Bid ];
      else
        Bid = [ Bid-1 Bid Bid+1 ];
      end
      
    elseif n == 3 && numel(ids) == 2
      
      xs = [ xs(1)   ( xs(1) + xs(2) )/2   xs(2) ];
      [ee,jj] = F( xs(2) );
      es = [ es(1)        ee               es(2) ];
      js = [ js(1)        jj               js(2) ];
      
    else
      
%       while numel( Bid ) < 1
% 
%       end

    end
    
    
    try, feval( P.VERBOSE , -8 , '    bests (%s) of %d\n' , sprintf('%d ', ids(Bid)) , last_stored ); end
    
    
    B = cell(1,n);
    for nn = 1:numel( Bid )
      B{nn} = { xs(Bid(nn))  es(Bid(nn))  js(Bid(nn)) };
    end
    
  end

end
