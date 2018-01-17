function [x,HISTORY,X_PATH,ALL_EVALS] = Optimize( F , x , varargin )
% 
%    x_optim                      = Optimize( F , x_init , ... 
%   [x_optim , history          ] = Optimize( F , x_init , ... 
%   [x_optim , history , x_path ] = Optimize( F , x_init , ... 
% 
%   x_optim : the optim x_found
%   history : a resume of the evaluation
%   x_path  : a cell with the path done by x
% 
%   
%   OPTIONS:
%     'methods' -> MET
%               -> { MET1 , MET2 , ... }
%               -> { MET1 , number_of_iterations , 
%                    MET2 , - time_per_method , 
%                    MET3 , number_of_iterations , - time_per_method , 
%                    MET4 , PARAMETERS_struct , ... 
%                    MET5 , {'param1',p1,'param2',p2} 
%                  }
% 
%       valids MET: 
%         only energy needed:
%                     'coordinate'
%                     'fminsearch'
%
%         jacobian needed:
%                     'descend'
%                     'conjugate'
%                     'quasinewton'
%                     'lmbfgs'
%                     'descendneg'
%                     'g-n.l-m'
%                     'ode'
% 
%         hessian needed: 
%                     'newton'
%                     'levenberg-marquardt'
% 
% 
% 
% GENERAL_OPTIONS: (these can't be modified for every method)
%   MAX_ITERATIONS                    (  Inf  )
%   MAX_TIME                          (  Inf  )
%   MIN_ENERGY                        ( -Inf  )
%   VERBOSE                           (  2    )
%   SAVE_IN_BASE                      ( false )
%   PLOT                              ( true  )
%   LS_PLOT                           ( false )
%   COMPUTE_NUMERICAL_JACOBIAN        ( true  )
%   COMPUTE_NUMERICAL_HESSIAN         ( true  )
%   USE_SMART_FUNCTIONS               ( true  )
% 
% 
% PARAMETERS: (for each method you can define it)
%   P.MAX_ITERATIONS_PER_METHOD     = Inf;
%   P.MAX_TIME_PER_METHOD           = Inf;
%   P.MIN_NORM_GRADIENT             = 0;
%  
%   P.
% 
%   P.LINE_SEARCH                   = {'golden'}
%   P.INITIAL_STEP                      = 0.5;   %if P.INITIAL_STEP is 0, force bracketing
%                                           %negative for factor of previous alpha
%                                           %positive value: a linesearch performed between [0,P.INITIAL_STEP]
% 
%   P.SHAMANSKII_P                  = 5;
% 
%   P.RESTART_CONJUGATE             = ceil( P.RESTART_CONJUGATE*numel(x) );
%   P.RESTART_QUASINEWTON           = ceil( P.RESTART_QUASINEWTON*numel(x) );
%   P.LMBFGS_ORDER                  = 15;
% 
%   P.PERFORM_FIRST_SEARCH          = 0;
% 
%   P.LS_WOLFE_DELTA                = 1; %actually it is an Armijo coeficient
%   P.LS_WOLFE_SIGMA                = 0.9;
%   P.LS_WOLFE_GAMMA                = 2/3;
%   P.LS_WOLFE_MIN_ALPHA            = 1e-20;
% 
%   P.LS_GOLDEN_MIN_DIFF            = 0;
%   P.LS_GOLDEN_MAX_ITS             = 20;
%   P.LS_GOLDEN_MIN_ALPHA           = 1e-20;
% 
%   P.LS_EXHAUSTIVE_RECUR           = 10;
%   P.LS_EXHAUSTIVE_N               = 50;
% 
%   P.LS_BACKTRACKING_FACTOR        = 0.66;
%   P.LS_BACKTRACKING_MIN_ALPHA     = 1e-20;  
% 
%   P.LS_QUADRATIC_IGNORE_SLOPE0    = 1;
%   P.LS_QUADRATIC_RECURSIONS       = 1;
% 
%   P.ALPHA_MINIMUM                 = 1e-8;
%   P.BRACKETING_RHO                = 5;
% 

  %%Inicializaciones
  % GO   -> global options
  GO.MAX_ITERATIONS                =  Inf;
  GO.MAX_TIME                      =  Inf;
  GO.MIN_ENERGY                    = -Inf;
%   GO.FIXED_COORDINATES             =  [];
  GO.VERBOSE                       =  1;
  GO.SAVE_IN_BASE                  =  true;  %store x in getappdata( 0 , 'OPTIMIZE_X')
  GO.PLOT                          =  true;
  GO.LS_PLOT                       =  false;
  GO.COMPUTE_NUMERICAL_JACOBIAN    = { '5' };  % false para no calcularlo
  GO.COMPUTE_NUMERICAL_HESSIAN     = { 'cc' , 'hessian' };
  GO.CHECK_INFS                    = [ minnum('single') , maxnum('single') ]; % false para no chequear
  GO.USE_SMART_FUNCTIONS           =  true;
  GO.VERBOSE_FILE                  =  ''; %'Optimize.log';
  GO.CLEAN_VERBOSE_FILE            =  true;
  GO.CLOSE_PLOT_AT_END             =  true;
  GO.SAVE_X_IN_FILE                =  []; %'Optimize_X.mat';      % [] -> no salva
  GO.SAVE_X_EVERY                  =  Inf;  % positive number -> cada n iteraciones
                                            % negative number -> cada -n seconds
                                            % Inf -> salva solamente al final

  %DP  -> default parameters
  DP.MAX_ITERATIONS_PER_METHOD     = Inf;
  DP.MAX_TIME_PER_METHOD           = Inf;
  DP.FILTER                        = [];
  DP.MIN_NORM_GRADIENT             = 0;
  DP.MIN_ENERGY_DELTA              = 0;
  
  DP.CONJUGATE_METHOD              = 'prp'; %{'hs','fr','prp','cd','ls','dy','hz','mprp','h-1','h-2','h-3','h-4','h-5','h-6'}
  DP.CONJUGATE_RESTART             = -4;    %negative means a factor of numel(x)
  
  DP.NEWTON_HESSIAN_MODIFICATION   = 'eigs'; %{'add','eye','eigs','none'}
  DP.NEWTON_SHAMANSKII_STEP        = 1;
  
  DP.QUASINEWTON_METHOD            = 'bfgs'; %{'bfgs','dfp','sr1'}
  DP.QUASINEWTON_RESTART           = -4;     %negative means a factor of numel(x)
  
  DP.LMBFGS_ORDER                  = 15;

  DP.DAMPED_INITIAL_LAMBDA         = 1/100;
  DP.DAMPED_USE_DIAG_H             = true;
  DP.DAMPED_REDUCE_FACTOR          = 1/10;
  DP.DAMPED_INCREASE_FACTOR        = 5;
  DP.DAMPED_MAX_LAMBDA             = 1e+10;
  
  DP.COORDINATES_ORDER             = @(N) randperm(N);   % 1:N
  
  
  DP.LINE_SEARCH                   = {'golden'}; %{'golden','backtracking','backtracking0','none','exhaustive','wolfe','derivative','fast','quadratic'}
  DP.MIN_STEP                      = 0;     %si vale 0 previene x == x + a*Dn
  DP.INITIAL_STEP                  = -1.5;  %if P.INITIAL_STEP is 0, force bracketing
                                            %negative for factor of
                                            %previous alpha, si alpha = NaN || 0, fuerza bracketing
                                            %positive value: a linesearch performed between [0,P.INITIAL_STEP]

  DP.SMART_LINE_SEARCH_ITS         = 2;
  DP.PERFORM_FIRST_SEARCH          = 1;  %si vale N, las primeras N iteraciones de cada metodo utiliza como LINE_SEARCH el ultimo de la lista (no lo repite si esta en otra posicion tambien)
                                         %ademas en estas primeras N iteraciones fuerza el bracketing
  DP.BRACKETING_RHO                = 2;
  DP.BRACKETING_MIN                = 0.1;   %el bracketing se hace hasta que se encuentra algun minimo local
  DP.BRACKETING_MAX                = 1e+6;
  DP.BRACKETING_ITS                = 3;

  
  DP.LS_GOLDEN_MIN_DIFF            = 0;
  DP.LS_GOLDEN_MAX_ITS             = 100;
  DP.LS_GOLDEN_MIN_ALPHA           = 1e-30;
  DP.LS_GOLDEN_MAX_ALPHA           = 1e+6;

  
  DP.LS_BACKTRACKING_MAX_ITS       = 100;
  DP.LS_BACKTRACKING_MIN_ALPHA     = 1e-30;
  DP.LS_BACKTRACKING_REDUCE_FACTOR = 0.85;
  DP.LS_BACKTRACKING_EXTRA_TESTS   = 5;

  DP.LS_QUADRATIC_MAX_ITS          = 3;
  DP.LS_QUADRATIC_MAX_ALPHA        = 1e+6;
  DP.LS_QUADRATIC_TOL              = 1/20;
  
  DP.LS_EXHAUSTIVE_ITERATIONS      = 10;
  DP.LS_EXHAUSTIVE_N               = 50;

  DP.ODE_EVOLUTION_TIME            = 1;     %tiempo de evolucion
  DP.ODE_TOLERANCE                 = 1e-3;
  DP.ODE_MIN_STEP                  = 0;     %if 0 automatically set a proper value
  DP.ODE_INITIAL_STEP              = 1e-3;
  DP.ODE_USE_ENERGY_VALUE          = true;
  DP.ODE_NORMALIZE_J               = false;

  

  
  GRADIENT_METHODS                  = {'lmbfgs','quasinewton','conjugate','descend','descendneg','ode'};
  HESSIAN_METHODS                   = {'newton','dampednewton'};
  FREE_METHODS                      = {'coordinate','partan','powell'};  %,'fminsearch'
  ALL_METHODS                       = [ GRADIENT_METHODS  HESSIAN_METHODS  FREE_METHODS ];
  
  ALL_CONJUGATE_METHODS             = {'hs','fr','prp','cd','ls','dy','hz','mprp','h-1','h-2','h-3','h-4','h-5','h-6'};
  ALL_NEWTON_HESSIAN_MODIFICATIONS  = {'add','eye','eigs','none'};
  ALL_QUASINEWTON_METHODS           = {'bfgs','dfp','sr1'};
  ALLS_LINE_SEARCHS                 = {'golden','backtracking','quadratic','exhaustive','fixed','smart'};  %,'none','wolfe','derivative','fast'
 
  %%END Inicializaciones
  
  
  if nargin == 0
    x = GO;
    HISTORY = DP;
    X_PATH = struct( 'GRADIENT_METHODS'                 , { GRADIENT_METHODS }                 ,...
                     'HESSIAN_METHODS'                  , { HESSIAN_METHODS  }                 ,...
                     'FREE_METHODS'                     , { FREE_METHODS     }                 ,...
                     'ALL_METHODS'                      , { ALL_METHODS      }                 ,...
                     'ALL_CONJUGATE_METHODS'            , { ALL_CONJUGATE_METHODS }            ,...
                     'ALL_NEWTON_HESSIAN_MODIFICATIONS' , { ALL_NEWTON_HESSIAN_MODIFICATIONS } ,...
                     'ALL_QUASINEWTON_METHODS'          , { ALL_QUASINEWTON_METHODS          } ,...
                     'ALLS_LINE_SEARCHS'                , { ALLS_LINE_SEARCHS                } );

    if nargout <= 1
      x = struct('OPTIONS' , x , 'PARAMETERS' , HISTORY , 'METHODS' , X_PATH );
    end
    return;
  end
  
  E = NaN;
  J = NaN;
  H = NaN;


  if ~isa( F , 'function_handle' )   , error('F have to be a ''function_handle''.');   end
  if ~isfloat( x )                   , error('X_init have to be a float.');            end
  

  if nargout > 1
    SAVE_HISTORY = true; HISTORY = [];
  else
    SAVE_HISTORY = false;
  end
  if nargout > 2
    SAVE_X_PATH = true;  X_PATH = {};
  else
    SAVE_X_PATH = false;
  end
  if nargout > 3
    SAVE_ALL_EVALS = true;
  else
    SAVE_ALL_EVALS = false;
  end


  %%metodos con los que optimiza
  [varargin,i,METHODS ] = parseargs( varargin , 'methods', 'Method','$DEFS$','descend');
  if ~iscell( METHODS ) , METHODS = {METHODS}; end

  %%chequea que ningun metodo de error
  for meth_id = 1:numel( METHODS )
    if ~ischar( METHODS{meth_id} ), continue; end
    if isequal( METHODS{meth_id} , '$BREAK$' ), continue; end
    P = parse_METHODS( METHODS , meth_id );
  end

  %%funcion que recibe METHODS y el indice y devuelve una estructura de
  %%parametros con los que correr el metodo
  function P = parse_METHODS( METHODS , meth_id )
    if ~ismember( strtrim( lower( METHODS{meth_id} )) , ALL_METHODS )
      error('Invalid method  '' %s '' .' , strtrim( lower( METHODS{meth_id} )) );
    end
    
    P = DP;  %assigning the default parameters
    if     numel(METHODS) > meth_id+1  &&  isnumeric( METHODS{meth_id+1} )  &&  isscalar( METHODS{meth_id+1} )  &&  isnumeric( METHODS{meth_id+2} )  &&  isscalar( METHODS{meth_id+2} )

      if     METHODS{meth_id+1} >= 0 ,  P.MAX_ITERATIONS_PER_METHOD =  METHODS{meth_id+1};
      elseif METHODS{meth_id+1} <  0 ,  P.MAX_TIME_PER_METHOD       = -METHODS{meth_id+1};
      end
      if     METHODS{meth_id+2} >= 0 ,  P.MAX_ITERATIONS_PER_METHOD =  METHODS{meth_id+2};
      elseif METHODS{meth_id+2} <  0 ,  P.MAX_TIME_PER_METHOD       = -METHODS{meth_id+2};
      end

    elseif numel(METHODS) > meth_id  &&  isnumeric( METHODS{meth_id+1} )  &&  isscalar( METHODS{meth_id+1} ) && ( numel(METHODS)==meth_id+1 || ischar( METHODS{meth_id+2} ) )

      if     METHODS{meth_id+1} >= 0 ,  P.MAX_ITERATIONS_PER_METHOD =  METHODS{meth_id+1};
      elseif METHODS{meth_id+1} <  0 ,  P.MAX_TIME_PER_METHOD       = -METHODS{meth_id+1};
      end

    elseif numel(METHODS) > meth_id  &&  iscell( METHODS{meth_id+1} )  &&  ( numel(METHODS)==meth_id+1 || ischar( METHODS{meth_id+2} ) )

      params = METHODS{meth_id+1};
      all_params_names = fieldnames( P );
      for p = 1:numel(all_params_names)
        [params,i,val] = parseargs( params , lower(all_params_names{p}) );
        if i ~= 0, P.(all_params_names{p}) = val; end
      end

      if ~isempty( params ),
        fprintf('unknow parameters!!!\n\n');
        fprintf('%s', params{1} );
        fprintf('\ntry remove it\n\n');
        error('INVALID PARSING OF METHODS.'); 
      end
    
    elseif numel(METHODS) > meth_id  &&  isstruct( METHODS{meth_id+1} )  &&  ( numel(METHODS)==meth_id+1 || ischar( METHODS{meth_id+2} ) )
      
      pnames = fieldnames( METHODS{meth_id+1} );
      newM  = cell( 1 , numel( pnames )*2 );
      for p = 1:numel( pnames )
        newM{2*p-1} = pnames{p};
        newM{2*p  } = METHODS{meth_id+1}.(pnames{p});
      end
      METHODS{meth_id+1} = newM;
      P = parse_METHODS( METHODS , meth_id );
    
    elseif numel(METHODS) == meth_id || ischar(METHODS{meth_id+1})
    else
      error('INVALID PARSING OF METHODS.');
    end
    

    %check some parameters
    try,
      feval( P.COORDINATES_ORDER , numel(x) );
    catch
      error('error in   COORDINATES_ORDER   function');
    end      
    
    if ~ischar( P.CONJUGATE_METHOD ) || ~any( strcmpi( P.CONJUGATE_METHOD , ALL_CONJUGATE_METHODS ) )
      error('incorrect CONJUGATE_METHOD');
    end
    if P.CONJUGATE_RESTART   < 0, P.CONJUGATE_RESTART    = ceil( -P.CONJUGATE_RESTART  * numel(x) ); end

    if ~ischar( P.NEWTON_HESSIAN_MODIFICATION ) || ~any( strcmpi( P.NEWTON_HESSIAN_MODIFICATION , ALL_NEWTON_HESSIAN_MODIFICATIONS ) )
      error('incorrect NEWTON_HESSIAN_MODIFICATION');
    end
    
    if ~ischar( P.QUASINEWTON_METHOD ) || ~any( strcmpi( P.QUASINEWTON_METHOD , ALL_QUASINEWTON_METHODS ) )
      error('incorrect QUASINEWTON_METHOD');
    end
    if P.QUASINEWTON_RESTART < 0, P.QUASINEWTON_RESTART  = ceil( -P.QUASINEWTON_RESTART * numel(x) ); end
    
    if ~iscell( P.LINE_SEARCH ), P.LINE_SEARCH = { P.LINE_SEARCH }; end
    if ~all( cellfun( @(x) any( strcmpi( x , ALLS_LINE_SEARCHS ) ) , P.LINE_SEARCH ) )
      error('incorrect LINE_SEARCH');
    end
    
    P.METHOD = strtrim( lower( METHODS{meth_id} ));
  end
  %%END parsing methods


  %%parseo otras cosas....
  %outputfcn
  [varargin,i,OUTPUT_FCN ] = parseargs( varargin , 'OUTPUTfcn' );
  
  %projection
  [varargin,i,PROJECT_FCN ] = parseargs( varargin , 'PROJECTIONfcn' );

  %inner product
  [varargin,i,INNER ] = parseargs( varargin , 'inner','$DEFS$', @(x,y) x(:)'*y(:) );
  NORM  = @(x) sqrt( INNER(x,x) );
  NORM2 = @(x) INNER(x,x);

  %filter descend direction
  [varargin,i, DP.FILTER ] = parseargs( varargin , 'filter' );

  %plot
  [varargin,GO.PLOT ] = parseargs( varargin , 'plot'  ,'$FORCE$', { true  , GO.PLOT } );
  [varargin,GO.PLOT ] = parseargs( varargin , 'noplot','$FORCE$', { false , GO.PLOT } );

  %lsplot
  [varargin,GO.LS_PLOT ] = parseargs( varargin , 'lsplot'    ,'$FORCE$', { true  , GO.LS_PLOT } );
  [varargin,GO.LS_PLOT ] = parseargs( varargin , 'nolsplot'  ,'$FORCE$', { false , GO.LS_PLOT } );
  
  %verbose
  [varargin,i,GO.VERBOSE ] = parseargs( varargin , 'verbose','$DEFS$', GO.VERBOSE );

  
  %line_search
  [varargin,i,DP.LINE_SEARCH ] = parseargs( varargin , 'LineSearch' ,'$DEFS$', DP.LINE_SEARCH );
  if ~iscell( DP.LINE_SEARCH ), P.LINE_SEARCH = { DP.LINE_SEARCH }; end
  if ~all( cellfun( @(x) any( strcmpi( x , ALLS_LINE_SEARCHS ) ) , DP.LINE_SEARCH ) )
    error('incorrect LINE_SEARCH');
  end

  %conjugate_methods
  [varargin,i,DP.CONJUGATE_METHOD ] = parseargs( varargin , 'conjugatemethod','ConjuGate','$DEFS$', DP.CONJUGATE_METHOD );
  if ~ischar( DP.CONJUGATE_METHOD ) || ~any( strcmpi( DP.CONJUGATE_METHOD , ALL_CONJUGATE_METHODS ) )
    error('incorrect CONJUGATE_METHOD');
  end

  %si queda algo en varargin deben ser GO y DP
  if numel( varargin )
    GO_ = varargin{1};     varargin(1) = [];
    fnames = fieldnames( GO_ );
    for p = 1:numel(fnames)
      if isfield( GO , upper( fnames{p} ) )
        GO.(upper( fnames{p} )) = GO_.(fnames{p});
      else
        warning(' cuidado!!! no entiendo esta opcion  %s',fnames{p});
      end
    end
  end

  if numel( varargin )
    DP_ = varargin{1};    varargin(1) = [];

    fnames = fieldnames( DP_ );
    for p = 1:numel(fnames)
      if isfield( DP , upper( fnames{p} ) )
        DP.(upper( fnames{p} )) = DP_.(fnames{p});
      else
        warning(' cuidado!!! no entiendo este parametro  %s',fnames{p});
      end
    end
  end
  %%END creo que ya termine de parsear toda la entrada...
  
  CLEANUP = onCleanup( @() CLEANUP_FUNCTION(GO) );
  function CLEANUP_FUNCTION(GO)
    if ~isempty( getappdata(0,'OPTIMIZE_START') )   &&  getappdata(0,'OPTIMIZE_START')
      fileprintf( GO.VERBOSE_FILE , '\n\nUSER INTERRUP ????:    %s \n\n-=*=--=*=--=*=--=*=--=*=--=*=--=*=--=*=--=*=--=*=-\n\n' , datestr(now, ' HH:MM:SS.FFF       ( mmmm dd )')  );
      rmappdata( 0 , 'OPTIMIZE_START' );
    end
  end
  
  printMARGIN = 0;
  
  %%preparando la funcion
  N  = numel( x ); NN = N*N;
  
%   if ~isempty( GO.FIXED_COORDINATES ) && max( GO.FIXED_COORDINATES(:) ) > N
%     error( 'error en FIXED_COORDINATES');
%   end
  
  sz = size( x );
  x  = x(:);
  NOF_F_evals       = uint32(0);
  NOF_J_evals       = uint32(0);
  NOF_H_evals       = uint32(0);
  NOF_Jn_evals      = uint32(0);
  NOF_Hn_evals      = uint32(0);
  NOF_F_smarts      = uint32(0);
  NOF_F_no_compute  = uint32(0);
  NOF_J_smarts      = uint32(0);
  NOF_H_smarts      = uint32(0);
  NOF_LS_smarts     = uint32(0);
  try
    [E,J,H] = F( R(x) );
    J = reshape( J , N , 1 );
    H = reshape( H , N , N );    
    F_RETURN_J = true;
    F_RETURN_H = true;
    NOF_F_evals = NOF_F_evals + 1;
    NOF_J_evals = NOF_J_evals + 1;
    NOF_H_evals = NOF_H_evals + 1;
    Vprintf(-4,'NOF: first evaluation(f,j,h): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);
  catch LE
    if isequal( LE.identifier , 'MATLAB:TooManyOutputs' )  ||  isequal( LE.identifier , 'MATLAB:maxlhs' )  || isequal( LE.identifier , 'deliver:TooManyOutputs' )  ||   isequal( LE.identifier , 'MATLAB:too_few_values_for_assignment' )
      F_RETURN_H = false;
    else
      error('La primera evaluacion de F , dio error retornando [F,J,H]');
    end
    try
      [E,J] = F( R(x) );
      J = reshape( J , N , 1 );

      F_RETURN_J = true;
      NOF_F_evals = NOF_F_evals + 1;
      NOF_J_evals = NOF_J_evals + 1;
      Vprintf(-4,'NOF: first evaluation(f,j): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);
    catch LE
      if isequal( LE.identifier , 'MATLAB:TooManyOutputs' )   ||  isequal( LE.identifier , 'MATLAB:maxlhs' )  ||   isequal( LE.identifier , 'MATLAB:too_few_values_for_assignment' )
        F_RETURN_J = false;
      else
        disperror(LE)
        error('La primera evaluacion de F , dio error retornando [F,J]');
      end
      try
        E = F( R(x) );
        NOF_F_evals = NOF_F_evals + 1;
        Vprintf(-4,'NOF: first evaluation(f): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);
      catch
        disperror(LE)
        error('Error evaluating function at frist callback!');
      end
    end
  end
  

%   J( GO.FIXED_COORDINATES ) = 0;

  if ~isscalar( E ) || isnan(E) || isinf(E) || ~isreal(E)
    error('Invalid function value (at X_init). ( isscalar: %d  -  isnan: %d  -  isinf: %d  -  isreal: %d )' ,...
        isscalar( E ) , isnan(E) , isinf(E) , isreal(E) );
  end
  if F_RETURN_J && (  numel(J) ~= N  ||  any(isnan(J(:)))  ||  any(isinf(J(:))) ||  any(~isreal(J(:))) )
    error('Invalid jacobian (at X_init). ( numel==N: %d  -  isnan: %d  -  isinf: %d  -  isreal: %d )' ,...
      numel(J) == N  ,  any(isnan(J(:))) ,  any(isinf(J(:))) ,  any(isreal(J(:))) );
  end
  if F_RETURN_H && (  numel(H) ~= NN ||  any(isnan(H(:)))  ||  any(isinf(H(:))) ||  any(~isreal(H(:))) )
    error('Invalid hessian (at X_init). ( numel==NN: %d  -  isnan: %d  -  isinf: %d  -  isreal: %d )' ,...
      numel(H) == N  ,  any(isnan(H(:))) ,  any(isinf(H(:))) ,  any(isreal(H(:))) );
  end
  
  if SAVE_ALL_EVALS, ALL_EVALS = { x(:).' , E }; end

  if GO.USE_SMART_FUNCTIONS
    % SF  -> smart_function
    SF.F.v = E;
    SF.F.x = R(x);
    if F_RETURN_J
      SF.J.v = J;
      SF.J.x = R(x);
    else
      SF.J.v = zeros(N,1); SF.J.v(:) = NaN;
      SF.J.x = R(x);       SF.J.x(:) = NaN;
    end
    if F_RETURN_H
      SF.H.v = H;
      SF.H.x = R(x);
    else
      try
      SF.H.v = zeros(N,N); SF.H.v(:) = NaN;
      SF.H.x = R(x);       SF.H.x(:) = NaN;
      end
    end
  else
    SF.F.x = NaN;
    SF.J.x = NaN;
    SF.H.x = NaN;
  end
  
  F = @(x, varargin )   FJH( F , R(x) , varargin{:} );

  function [f,j,h] = FJH( ff , x , ComputeF )
    fid = fopen( 'optimize_aux','a');
    fwrite( fid , x ,'double');
    fclose(fid);
      
      
    if any(isnan(x(:))) || any(isinf(x(:)))
      error('BAD  x!!!  %d (nans)  -  %d (inf)',sum(isnan(x(:))),sum(isinf(x(:))));
    end
    
    if nargin < 3
      ComputeF = true; 
    elseif ~islogical( ComputeF )
      error('ComputeE is expected to be logical !!');
    end
    
    if nargout <= 1

      if isequal( SF.F.x , x )    &&   ComputeF
        f = SF.F.v;   NOF_F_smarts = NOF_F_smarts + 1;
        Vprintf(-4,'NOF: (sf): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);
      elseif ComputeF
        try
        f = ff(x);
        catch LE
          setappdata(0,'OPTIMIZE_ERROR',x);
          rethrow( LE );
        end
        NOF_F_evals = NOF_F_evals + 1;
        Vprintf(-4,'NOF: (f): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);
      else
        f = SF.F.v;  NOF_F_no_compute = NOF_F_no_compute + 1;
        Vprintf(-4,'NOF: (nf): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);
      end

    elseif nargout == 2
      
      if isequal( SF.F.x , x ) && isequal( SF.J.x , x )

        j = SF.J.v;   NOF_J_smarts = NOF_J_smarts + 1;
        if ComputeF
          f = SF.F.v;  NOF_F_smarts = NOF_F_smarts + 1;
          Vprintf(-4,'NOF: (sf,sj): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);
        else
          f = SF.F.v;  NOF_F_no_compute = NOF_F_no_compute + 1;
          Vprintf(-4,'NOF: (nf,sj): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);
        end
        

      elseif   ~isequal( SF.F.x , x )   &&   isequal( SF.J.x , x )
        
        j = SF.J.v;   NOF_J_smarts = NOF_J_smarts + 1;
        if ComputeF
          try
          f = ff(x);
          catch LE
            setappdata(0,'OPTIMIZE_ERROR',x);
            rethrow( LE );
          end
          NOF_F_evals  = NOF_F_evals + 1;
          
          Vprintf(-4,'NOF: (f,sj): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);
        else
          f = SF.F.v;  NOF_F_no_compute = NOF_F_no_compute + 1;
          Vprintf(-4,'NOF: (nf,sj): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);
        end

      elseif   ~isequal( SF.J.x , x )   &&   F_RETURN_J

        if ComputeF
          try
          [f,j] = ff(x); %j(GO.FIXED_COORDINATES) = 0;
          catch LE
            setappdata(0,'OPTIMIZE_ERROR',x);
            rethrow( LE );
          end
          
          NOF_F_evals = NOF_F_evals + 1;
          NOF_J_evals = NOF_J_evals + 1;
          Vprintf(-4,'NOF: (f,j): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);

        else
          
          setappdata( 0 , 'OPTIMIZE_NO_COMPUTE_F' , 1 );
          try
          [f,j] = ff(x);  %j(GO.FIXED_COORDINATES) = 0;
          catch LE
            setappdata(0,'OPTIMIZE_ERROR',x);
            rethrow( LE );
          end
          
          rmappdata( 0 , 'OPTIMIZE_NO_COMPUTE_F' );

          f = SF.F.v;
          NOF_F_no_compute = NOF_F_no_compute + 1;
          NOF_J_evals      = NOF_J_evals + 1;
          Vprintf(-4,'NOF: (nf,j): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);

        end


      elseif   ~isequal( SF.J.x , x )   &&   ~F_RETURN_J
        
        if ComputeF
        
          try
          f = ff(x);
          catch LE
            setappdata(0,'OPTIMIZE_ERROR',x);
            rethrow( LE );
          end

          NOF_F_evals = NOF_F_evals + 1;

          if isfalse( GO.COMPUTE_NUMERICAL_JACOBIAN )   ,  error('The function has no jacobian.');  end
          try
            j = NumericalDiff( @(z) ff(z) , x , GO.COMPUTE_NUMERICAL_JACOBIAN{:} );
            %j(GO.FIXED_COORDINATES) = 0;
          catch
            setappdata(0,'OPTIMIZE_ERROR_COMPUTING_JACOBIAN',x);
            error('error computing numerical jacobian.');
          end
          NOF_Jn_evals = NOF_Jn_evals + 1;
          Vprintf(-4,'NOF: (f,jn): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);
        
        else
          
          f = SF.F.v;   NOF_F_no_compute = NOF_F_no_compute + 1;

          if isfalse( GO.COMPUTE_NUMERICAL_JACOBIAN )   ,  error('The function has no jacobian.');  end
          try
            j = NumericalDiff( @(z) ff(z) , x , GO.COMPUTE_NUMERICAL_JACOBIAN{:} );
            %j(GO.FIXED_COORDINATES) = 0;
          catch
            setappdata(0,'OPTIMIZE_ERROR_COMPUTING_JACOBIAN',x);
            error('error computing numerical jacobian.');
          end
          NOF_Jn_evals = NOF_Jn_evals + 1;
          Vprintf(-4,'NOF: (nf,jn): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);
          
        end
          
      else

        error('no deberia estar aqui!!!!!!!!!');

      end
      
    elseif nargout == 3
      
      if isequal( SF.F.x , x ) && isequal( SF.J.x , x ) && isequal( SF.H.x , x )
        f = SF.F.v;   NOF_F_smarts = NOF_F_smarts + 1;
        j = SF.J.v;   NOF_J_smarts = NOF_J_smarts + 1;
        h = SF.H.v;   NOF_H_smarts = NOF_H_smarts + 1;
        Vprintf(-4,'NOF: (sf,sj,sh): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);
      elseif isequal( SF.J.x , x ) && isequal( SF.H.x , x )
        try
          f = ff(x);
        catch LE
          setappdata(0,'OPTIMIZE_ERROR',x);
          rethrow( LE );
        end
        NOF_F_evals = NOF_F_evals + 1;
        j = SF.J.v;   NOF_J_smarts = NOF_J_smarts + 1;
        h = SF.H.v;   NOF_H_smarts = NOF_H_smarts + 1;
        Vprintf(-4,'NOF: (f,sj,sh): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);
      elseif isequal( SF.H.x , x )
        h = SF.H.v;   NOF_H_smarts = NOF_H_smarts + 1;
        if F_RETURN_J
          try
          [f,j] = ff(x); %j(GO.FIXED_COORDINATES) = 0;
          catch LE
            setappdata(0,'OPTIMIZE_ERROR',x);
            rethrow( LE );
          end
          
          NOF_F_evals = NOF_F_evals + 1;
          NOF_J_evals = NOF_J_evals + 1;
          Vprintf(-4,'NOF: (f,j,sh): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);
        else
          if isfalse( GO.COMPUTE_NUMERICAL_JACOBIAN )   ,  error('The function has no jacobian.');  end
          try
            j = NumericalDiff( @(z) ff(z) , x , GO.COMPUTE_NUMERICAL_JACOBIAN{:} );
            %j(GO.FIXED_COORDINATES) = 0;
          catch
            setappdata(0,'OPTIMIZE_ERROR_COMPUTING_JACOBIAN',x);
            
            error('error computing numerical jacobian.');
          end
          NOF_Jn_evals = NOF_Jn_evals + 1;
          try
            f = ff(x);
          catch LE
            setappdata(0,'OPTIMIZE_ERROR',x);
            rethrow( LE );
          end
          
          NOF_F_evals = NOF_F_evals + 1;
          Vprintf(-4,'NOF: (f,jn,sh): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);
        end
      else
        if       F_RETURN_J  &&  F_RETURN_H

          try
            [f,j,h] = ff(x); %j(GO.FIXED_COORDINATES) = 0;
          catch LE
            setappdata(0,'OPTIMIZE_ERROR',x);
            rethrow( LE );
          end
          
          NOF_F_evals = NOF_F_evals + 1;
          NOF_J_evals = NOF_J_evals + 1;
          NOF_H_evals = NOF_H_evals + 1;
          Vprintf(-4,'NOF: (f,j,h): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);

        elseif  F_RETURN_J  && ~F_RETURN_H

          if isfalse( GO.COMPUTE_NUMERICAL_HESSIAN )    ,  error('The function has no hessian.');   end
          try
            h = NumericalDiff( @(z) ff(z) , x , GO.COMPUTE_NUMERICAL_HESSIAN{:} );
          catch
            setappdata(0,'OPTIMIZE_ERROR_COMPUTING_HESSIAN',x);
            error('error computing numerical hessian.');
          end            
          NOF_Hn_evals = NOF_Hn_evals + 1;
          
          try
            [f,j] = ff(x);
          catch LE
            setappdata(0,'OPTIMIZE_ERROR',x);
            rethrow( LE );
          end
          
          %j(GO.FIXED_COORDINATES) = 0;
          NOF_F_evals = NOF_F_evals + 1;
          NOF_J_evals = NOF_J_evals + 1;
          Vprintf(-4,'NOF: (f,j,hn): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);

        elseif ~F_RETURN_J  && ~F_RETURN_H
          
          if isfalse( GO.COMPUTE_NUMERICAL_HESSIAN )    ,  error('The function has no hessian.');   end
          try
            h = NumericalDiff( @(z) ff(z) , x , GO.COMPUTE_NUMERICAL_HESSIAN{:} );
          catch
            setappdata(0,'OPTIMIZE_ERROR_COMPUTING_HESSIAN',x);
            error('error computing numerical hessian.');
          end            
          NOF_Hn_evals = NOF_Hn_evals + 1;

          if isfalse( GO.COMPUTE_NUMERICAL_JACOBIAN )   ,  error('The function has no jacobian.');  end
          try
            j = NumericalDiff( @(z) ff(z) , x , GO.COMPUTE_NUMERICAL_JACOBIAN{:} );
            %j(GO.FIXED_COORDINATES) = 0;
          catch
            setappdata(0,'OPTIMIZE_ERROR_COMPUTING_JACOBIAN',x);
            error('error computing numerical jacobian.');
          end
          NOF_Jn_evals = NOF_Jn_evals + 1;

          try
            f = ff(x);
          catch LE
            setappdata(0,'OPTIMIZE_ERROR',x);
            rethrow( LE );
          end
          
          NOF_F_evals = NOF_F_evals + 1;
          Vprintf(-4,'NOF: (f,jn,hn): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);
          
        end
        
      end

    end

    
    if ~isscalar(f)
      setappdata(0,'OPTIMIZE_LAST_X',x);
      error('Invalid function value (it is not a scalar!!)...   x_error = getappdata(0,''OPTIMIZE_LAST_X'');');
    end
    if ~isreal(f)
      setappdata(0,'OPTIMIZE_LAST_X',x);
      error('Invalid function value (it is not a real!!)...   x_error = getappdata(0,''OPTIMIZE_LAST_X'');');
    end
    if isnan(f)
      setappdata(0,'OPTIMIZE_LAST_X',x);
      error('Invalid function value (it is NaN!!)...   x_error = getappdata(0,''OPTIMIZE_LAST_X'');');
    end
    if isinf(f)
      if isfalse( GO.CHECK_INFS )
        setappdata(0,'OPTIMIZE_LAST_X',x);
        error('Invalid function value (it is inf!!)...   x_error = getappdata(0,''OPTIMIZE_LAST_X'');');
      else
        Vwarning('OPTIMIZE:InfinityValue','Invalid function value (it is inf!!)...');
        if f < 0
          f = GO.CHECK_INFS(1); 
        else
          f = GO.CHECK_INFS(2); 
        end
      end
    end
    
    if GO.USE_SMART_FUNCTIONS,  SF.F.v = f; SF.F.x = x;  end
    
    if SAVE_ALL_EVALS
      if ~any( all( bsxfun( @eq , ALL_EVALS{1} , x(:).' ) , 2 ) )
        ALL_EVALS{1} = [ ALL_EVALS{1} ; x(:).' ];
        ALL_EVALS{2} = [ ALL_EVALS{2} ;  f     ];
      end
    end

    if nargout > 1
      if numel(j) ~= N
        setappdata(0,'OPTIMIZE_LAST_X',x);
        error('Invalid jacobian (it has not N components)...   x_error = getappdata(0,''OPTIMIZE_LAST_X'');');
      end
      if any( ~isreal(j(:)) )
        setappdata(0,'OPTIMIZE_LAST_X',x);
        error('Invalid jacobian (it has complex components)...   x_error = getappdata(0,''OPTIMIZE_LAST_X'');');
      end
      if any( isnan(j(:)) )
        setappdata(0,'OPTIMIZE_LAST_X',x);
        error('Invalid jacobian (it has NaNs components)...   x_error = getappdata(0,''OPTIMIZE_LAST_X'');');
      end
      if any( isinf(j(:)) )
        if isfalse( GO.CHECK_INFS )
          setappdata(0,'OPTIMIZE_LAST_X',x);
          error('Invalid jacobian (it has infs!!)...   x_error = getappdata(0,''OPTIMIZE_LAST_X'');');
        else
          Vwarning('OPTIMIZE:InfinityJacobian','Invalid jacobian (it has infs!!)...');
          j( isinf(j) & j < 0 ) = GO.CHECK_INFS(1);
          j( isinf(j) & j > 0 ) = GO.CHECK_INFS(2);
        end
      end

      j = reshape( j , N , 1 );
      if GO.USE_SMART_FUNCTIONS,  SF.J.v = j; SF.J.x = x;  end
    end


    if nargout > 2
      if numel(h) ~= NN
        setappdata(0,'OPTIMIZE_LAST_X',x);
        error('Invalid hessian (it has not N components)...   x_error = getappdata(0,''OPTIMIZE_LAST_X'');');
      end
      if any( ~isreal(h(:)) )
        setappdata(0,'OPTIMIZE_LAST_X',x);
        error('Invalid hessian (it has complex components)...   x_error = getappdata(0,''OPTIMIZE_LAST_X'');');
      end
      if any( isnan(h(:)) )
        setappdata(0,'OPTIMIZE_LAST_X',x);
        error('Invalid hessian (it has NaNs components)...   x_error = getappdata(0,''OPTIMIZE_LAST_X'');');
      end
      if any( isinf(h(:)) )
        if isfalse( GO.CHECK_INFS )
          setappdata(0,'OPTIMIZE_LAST_X',x);
          error('Invalid hessian (it has infs!!)...   x_error = getappdata(0,''OPTIMIZE_LAST_X'');');
        else
          Vwarning('OPTIMIZE:InfinityHessian','Invalid hessian (it has infs!!)...');
          h( isinf(j) & h < 0 ) = GO.CHECK_INFS(1);
          h( isinf(j) & h > 0 ) = GO.CHECK_INFS(2);
        end
      end

      h = reshape( h , N , N );
      if GO.USE_SMART_FUNCTIONS,  SF.H.v = h; SF.H.x = x;  end
    end

  end
  %%END la funcion F ya esta lista!!!
  

  %%ya que tengo el valor de E ... pruebo que OUTPUT_FCN y PROJECT_FCN no den error!!!
  if ~isempty( OUTPUT_FCN )
    try,   OUTPUT_FCN( x , E , 0 ); drawnow;
    catch LE
      disperror(LE);
      error(' error  evaluating OUTPUT_FCN ' );
    end
  end

  if ~isempty( PROJECT_FCN )
    try,   x = PROJECT_FCN( x , 0 );
    catch LE
      disperror(L )
      error(' error  evaluating PROJECT_FCN ' ); 
    end
  end
  %%END pruebo que OUTPUT_FCN y PROJECT_FCN no den error
  
  
  
  %%preparando plots, history, etc...
  if GO.PLOT

    PLOT.hFig = figure( 'Units','pixels'                                ,...
                        'IntegerHandle','off'                           ,...
                        'Toolbar','none'                                ,...
                        'NextPlot','new'                                ,...
                        'MenuBar','none'                                ,...
                        'MenuBar','none'                                ,...
                        'DockControls','on'                             ,...
                        'DoubleBuffer','on'                             ,...
                        'Renderer','OpenGL'                             ,...
                        'RendererMode','Manual'                         ,...
                        'IntegerHandle','off'                           ,...
                        'DockControls','off'                            ,...
                        'Interruptible','on'                            ,...
                        'BusyAction','cancel'                           ,...
                        'Position',[100 500 600 300]                    ,...
                        'NumberTitle','off'                             ,...
                        'Name','Optimizing' );
    AlwaysOnTopButton( PLOT.hFig );

    PLOT.hAxeTime   = axes('Parent', PLOT.hFig , 'Position',[0.1 0.05 0.8 0.9],'yAxisLocation','right','xtick',[] );
    PLOT.hLineTime  = line('Parent', PLOT.hAxeTime,'XData',0,'YData',0,'Color',[0 0 1],'Marker','none','linestyle','-');

    PLOT.hAxe       = axes('Parent', PLOT.hFig , 'Position',[0.1 0.05 0.8 0.9],'Color','none');
    PLOT.hLine      = line('Parent', PLOT.hAxe, 'XData',0,'YData',E,'Color',[1 0 0],'Marker','.','linestyle','-');

    uicontrol( 'Parent',PLOT.hFig , 'style','togglebutton','value',0,'String','Log','callback',@(h,e) PLOT_Y_LOG_SCALE( get(h,'Value') ),'Position',[2 2 30 14]);

  else

    GO.LS_PLOT = false;

  end
  if GO.LS_PLOT
    PLOT.hAxeLS = axes('Parent',PLOT.hFig,'Position',[0.4 0.6 0.4 0.4],'box','on');
    PLOT.LS_control = uicontrol('Parent',PLOT.hFig,'style','slider','callback',@(h,e) ShowLS( get(h,'Value') ),'Position',[2 20 60 14]);
    PLOT.hGLS = [];
    PLOT.LS_X = 0;
  end
  
  
  if ~isempty( GO.VERBOSE_FILE )
    try
      if GO.CLEAN_VERBOSE_FILE,  fid = fopen( GO.VERBOSE_FILE , 'w' );
      else,                      fid = fopen( GO.VERBOSE_FILE , 'a' ); fprintf( fid , '\n\n\n');
      end
      fprintf( fid , '**************************************************************\n');
      fprintf( fid , 'STARTING OPTIMIZATION: ' );
      fprintf( fid , datestr(now, ' HH:MM:SS.FFF       mmmm dd, yyyy  ( dddd )\n') );
      if isunix
        fprintf( fid , '          RUNNING IN :  %s' , HostName( ) );
      elseif ispc
        fprintf( fid , '          RUNNING IN :  WINDOWS' );
      end
      fprintf( fid , '   ( pid:  .%d.  )\n' , feature('getpid') );
      
      fprintf( fid , '                 pwd :  %s\n\n' , pwd );

      fprintf( fid , '                    N:  %d\n' , N );
      fprintf( fid , '    Jacobian Computed:  ' );
      if F_RETURN_J
        fprintf( fid , 'in function' );
      elseif isfalse( GO.COMPUTE_NUMERICAL_JACOBIAN )
        fprintf( fid , 'not computed' );
      else
        fprintf( fid , '%s', uneval( GO.COMPUTE_NUMERICAL_JACOBIAN ) );
      end
      fprintf( fid , '\n' );
      

      fprintf( fid , '     Hessian Computed:  ' );
      if F_RETURN_H
        fprintf( fid , 'in function' );
      elseif isfalse( GO.COMPUTE_NUMERICAL_HESSIAN )
        fprintf( fid , 'not computed' );
      else
        fprintf( fid , '%s', uneval( GO.COMPUTE_NUMERICAL_HESSIAN ) );
      end
      
      fprintf( fid , '\n       GLOBAL OPTIONS:  %s\n' , strrep( strrep( strrep( evalc( 'GO' ) , ' ' , '' ) , char(10) , '  ;  ' ) , ':' , ': ') );
      fprintf( fid , '\n**************************************************************\n\n');
      fclose(fid);
    catch
      try, fclose(fid); end
      warning('\n error open verbose file:  %s \n' , GO.VERBOSE_FILE );
      GO.VERBOSE_FILE = [];
    end
  end

  Vprintf(-4,'');
  Vprintf(-4,'');
  Vprintf(1, '                E_Init: %30.20g\n' , E );
  %%END preparando plots...
  
  
  
  %%algunas variables que necesitare
  METHOD_STARTTIME = clock;
  METHOD_IT = 0;
  alpha = NaN;
  
  try
  %%procedo con la optimizacion
  setappdata( 0 , 'OPTIMIZE_START' , true );
  
  FLS_a = [ ]; FLS_e = [ ];
  IT = 0; E0 = E; Ep = E; STARTTIME = clock; LAST_SAVE_TIME = STARTTIME; LAST_SAVE_IT = 0; LS = []; dE = Inf; small_dE_times = 0;
  while checkUSER_CANCEL() && check_IT()  &&  check_STARTTIME()  &&  check_E()
    %EXPERIMENTAL
    if small_dE_times > numel(METHODS) * 5
      Vprintf( 0 , 'deberia cortar...  parece que ya no avanza mas o desciende muy lento\n\n' );
      
      HISTORY.END_CONDITION = 'Small dE many times';
      break;
    end
    %

    x_old = x;
    for meth_id = 1:numel(METHODS)
      if ~checkUSER_CANCEL() || ~check_IT()  ||  ~check_STARTTIME()  ||  ~check_E(), break; end

      if isequal( METHODS{meth_id} , '$BREAK$' ), break; end
      
      if ~ischar( METHODS{meth_id} ), continue; end
      P = parse_METHODS( METHODS , meth_id );
      
      beta    = NaN;
      lambda  = NaN;
      switch P.METHOD

        case 'ode', START('ODE');
          H = NaN;
          if is1nan(J) || is1nan( E )
            [E,J] = F(x); 
          end

          if P.ODE_NORMALIZE_J, J = J/norm(J); end
          J = FILTER( J );

          alphamin = P.ODE_MIN_STEP;
          alpha = P.ODE_INITIAL_STEP;
          t = 0;
          while checkUSER_CANCEL() && check_METHOD_IT() && check_METHOD_STARTTIME() &&  check_IT()  &&  check_STARTTIME()  &&  check_E()  &&  check_dE()  &&  check_MIN_NORM_GRADIENT()

            if ( P.ODE_EVOLUTION_TIME - t ) < alpha, alpha = P.ODE_EVOLUTION_TIME - t; end

            if ~P.ODE_MIN_STEP
              alphamin = eps( min( abs( x./J ) ) ) * 0.05;   %   x + alphamin*Dn == x
            end

            Vprintf(2,' ');
            while 1
              Vprintf( 2 , '\btry to iter with alpha = %g' , alpha );

              k1 = - alpha * J;

              [ kk , k2 ] = F( x + 0.5*k1           , false ); if P.ODE_NORMALIZE_J, k2 = k2/norm(k2); end; k2 = - alpha * FILTER( k2 );
              [ kk , k3 ] = F( x + 0.5*k1 + 0.75*k2 , false ); if P.ODE_NORMALIZE_J, k3 = k3/norm(k3); end; k3 = - alpha * FILTER( k3 );

              xnew = x + 0.2222222222222222*k1 + 0.3333333333333333*k2 + 0.4444444444444444*k3 ;

              [ Ealpha , Jalpha  ] = F( xnew );
              if P.ODE_NORMALIZE_J, Jalpha = Jalpha/norm(Jalpha); end
              Jalpha = FILTER( Jalpha );

              if alpha <= alphamin
                Vprintf( 2 , '\balpha too small, no allow more reductions' );
                err = -Inf;
                break;
              end
              
              k4 = - alpha * Jalpha;

              err = -0.06944444444444445*k1 +  0.08333333333333333*k2 + 0.1111111111111111*k3 - 0.125*k4;
              err = max( abs( err ) );

              if err > P.ODE_TOLERANCE*2
                alpha = alpha/( 1.2 * realpow( err/P.ODE_TOLERANCE , 0.3333333333333333 ) );
                Vprintf( 2 , '\bDecreasing alpha because err too large ( err: %g  ,  newalpha:  %g )' , err , alpha );
              elseif P.ODE_USE_ENERGY_VALUE
                if Ealpha > E
                  alpha = alpha * 0.5;
                  Vprintf( 2 , '\bDecreasing alpha because Ealpha > E' );
                else
                  break;
                end
              else
                break;
              end
              if alpha < alphamin
                alpha = alphamin;
                Vprintf( 2 , '\balpha too small, setting it to alphamin ( alphamin = %g )' , alphamin );
              end
            end
            Vprintf( 2 , '\b' );
            
            
            if P.ODE_USE_ENERGY_VALUE   &&   Ealpha >= E
              Vprintf( 1 , 'Imposible to descend\n' );
              break;
            end

            Dn = 0.2222222222222222*k1 + 0.3333333333333333*k2 + 0.4444444444444444*k3;
            LINE_SEARCH( struct('LINE_SEARCH',{{'fixed'}},'PERFORM_FIRST_SEARCH',0,'INITIAL_STEP',1,'MIN_STEP',1),1,Ealpha );

            t = t + alpha;
            x = xnew;
            
            E = Ealpha; J = Jalpha;
            updateX( 'alpha:' , alpha , 'time:' , t , 'err:' , err );

            if t >= P.ODE_EVOLUTION_TIME
              Vprintf( 1 , 'EVOLUTION_TIME reached\n' );
              break;
            end
            
            if ~isinf( err )
              alpha  = alpha/max( 0.02 , 1.2 * realpow( err/P.ODE_TOLERANCE , 0.3333333333333333 ) );
              Vprintf( -2 , 'new value of alpha:  %g ', alpha );
            end
            
          end

        case 'descend', START('DESCEND');
          H = NaN;
          if is1nan(J) || is1nan( E )
            [E,J] = F(x); 
          end
          while checkUSER_CANCEL() && check_METHOD_IT() && check_METHOD_STARTTIME() &&  check_IT()  &&  check_STARTTIME()  &&  check_E()  &&  check_dE()  &&  check_MIN_NORM_GRADIENT()
            Dn     = - J/norm(J);
            Dn     = FILTER( Dn );
            SLOPE0 = INNER( J , Dn );
            
            [ alpha , M ] = LINE_SEARCH();
            if  M < E  &&  alpha > 0
              x = x + Dn*alpha;
              [E,J] = F(x);
              if ~isequal(E,M), Vwarning('New Energy value different to the minimum computed at LINESEARCH'); end
              updateX( 'alpha:' , alpha , 'norm(J):',norm(J) );
            else
              Vprintf( 1 , 'no descent direction\n' );
              break;
            end
          end

        case 'descendneg', START('DESCEND POSITIVE and NEGATIVE');
          H = NaN;
          if is1nan(J) || is1nan( E )
            [E,J] = F(x); 
          end
          while checkUSER_CANCEL() && check_METHOD_IT() && check_METHOD_STARTTIME() &&  check_IT()  &&  check_STARTTIME()  &&  check_E()  &&  check_dE()  &&  check_MIN_NORM_GRADIENT()
            Dn     = - J/norm(J);
            Dn     = FILTER( Dn );
            SLOPE0 = INNER( J , Dn );
            
            [ alpha , M ] = DOUBLE_LINE_SEARCH();
            if  M < E
              x = x + Dn*alpha;
              [E,J] = F(x);
              if ~isequal(E,M), Vwarning('New Energy value different to the minimum computed at LINESEARCH'); end
              updateX( 'alpha:' , alpha , 'norm(J):',norm(J) );
              continue;
            else
              Vprintf( 1 , 'no descent direction\n' );
              break;
            end
              
          end

        case 'conjugate'
          switch P.CONJUGATE_METHOD
            case 'hs'  , START('CONJUGATE GRADIENT (HESTENES-STIEFEL)');
            case 'fr'  , START('CONJUGATE GRADIENT (FLETCHER-REEVES)');
            case 'prp' , START('CONJUGATE GRADIENT (POLAK-RIBIERE-POLYAK)');
            case 'cd'  , START('CONJUGATE GRADIENT (CONJUGATE DESCEND)');
            case 'ls'  , START('CONJUGATE GRADIENT (LIU-STOREY)');
            case 'dy'  , START('CONJUGATE GRADIENT (DAI-YUAN)');
            case 'hz'  , START('CONJUGATE GRADIENT (HAGER-ZHANG)');
            case 'mprp', START('CONJUGATE GRADIENT (MODIFIED PRP)');
            case 'h-1' , START('CONJUGATE GRADIENT (HYBRID 1 [PRP-FR])');
            case 'h-2' , START('CONJUGATE GRADIENT (HYBRID 2 [PRP-FR])');
            case 'h-3' , START('CONJUGATE GRADIENT (HYBRID 3 [PRP-FR])');
            case 'h-4' , START('CONJUGATE GRADIENT (HYBRID 4 [HS-DY])');
            case 'h-5' , START('CONJUGATE GRADIENT (HYBRID 5 [HS-DY])');
            case 'h-6' , START('CONJUGATE GRADIENT (HYBRID 6)');
          end
          H = NaN;
          if is1nan(J) || is1nan( E )
            [E,J] = F(x); 
          end
          beta  = 0;
          D     = J;  %solo para inicializar y reservar ya la memoria (podria valer 0)
          J0    = J;  %solo para inicializar y reservar ya la memoria (podria no existir)

          while checkUSER_CANCEL() && check_METHOD_IT() && check_METHOD_STARTTIME() &&  check_IT()  &&  check_STARTTIME()  &&  check_E()  &&  check_dE()  &&  check_MIN_NORM_GRADIENT()
            D      = - J + beta*D;
            Dn     = D/norm(D);
            Dn     = FILTER( Dn );
            SLOPE0 = INNER( J , Dn );
            
            [ alpha , M ] = LINE_SEARCH();
            if  M < E  &&  alpha > 0
              x     = x + Dn*alpha;
              J0    = J;
              [E,J] = F(x);
              if ~isequal(E,M), 
                Vwarning('New Energy value different to the minimum computed at LINESEARCH'); end
              updateX( 'alpha:' , alpha , 'norm(J):',norm(J) ,'beta: ', beta );
            else
              Vprintf( 1 , 'no descent direction\n' );
              break;
            end

            if ~isinf( P.CONJUGATE_RESTART ) && ~mod( METHOD_IT , P.CONJUGATE_RESTART )
              Vprintf( -2 ,'Restarting Conjugate');
              beta = 0;
            else
              switch P.CONJUGATE_METHOD
                case 'hs'
                  beta = INNER( J , ( J - J0 ) ) / INNER( D , ( J - J0 ) );
                  beta = max( beta , 0 );
                case 'fr'
                  beta = NORM2( J ) / NORM2( J0 );
                  beta = max( beta , 0 );
                case 'prp'
                  beta = INNER( J , ( J - J0 ) ) / NORM2( J0 );
                  beta = max( beta , 0 );
                case 'cd'
                  beta = - NORM2( J ) / INNER( D , J0 );
                  beta = max( beta , 0 );
                case 'ls'
                  beta = - INNER( J , ( J - J0 ) ) / INNER( D , J0 );
                  beta = max( beta , 0 );
                case 'dy'
                  beta = NORM2( J ) / INNER( D , ( J - J0 ) );
                  beta = max( beta , 0 );
                case 'hz'
                  y = ( J - J0 );
                  beta = INNER( ( y - 2*D*NORM2( y )/INNER( D , y ) ) , J ) / INNER( D , y );
                  beta = max( beta , -1/( NORM(D) * min( 0.01 , NORM(J0) ) ) );
                case 'mprp'
                  y = ( J - J0 );
                  beta = INNER(J,y) / NORM2(J0);
                  beta = beta - min( beta , 1/2 * NORM2(y)/NORM2(J0)^2 * INNER(J,D) );
                case 'h-1'
                  beta_prp = INNER( J,( J - J0 ) ) / NORM2( J0 );
                  beta_fr  = NORM2(J) / NORM2(J0);
                  if beta_prp > 0 && beta_prp < beta_fr
                    beta = beta_prp;
                  else
                    beta = beta_fr;
                  end
                case 'h-2'
                  beta_prp = INNER( J,( J - J0 ) ) / NORM2( J0 );
                  beta_fr  = NORM2(J) / NORM2(J0);
                  beta = max( 0 , min( beta_prp , beta_fr ) );
                case 'h-3'
                  beta_prp = INNER( J,( J - J0 ) ) / NORM2( J0 );
                  beta_fr  = NORM2(J) / NORM2(J0);
                  beta = max( -beta_fr , min( beta_prp , beta_fr ) );
                case 'h-4'
                  beta_hs  = INNER( J , ( J - J0 ) ) / INNER( D , ( J - J0 ) );
                  beta_dy  = NORM2( J ) / INNER( D , ( J - J0 ) );
                  beta = max(-(1-sigma)/(1+sigma)*beta_dy , min( beta_hs , beta_dy ) );
                case 'h-5'
                  beta_hs  = INNER( J , ( J - J0 ) ) / INNER( D , ( J - J0 ) );
                  beta_dy  = NORM2( J ) / INNER( D , ( J - J0 ) );
                  beta = max( 0 , min( beta_hs , beta_dy ) );
                case 'h-6'
                  beta = NORM2( J ) / max( INNER(D,(J-J0)) , -INNER(J0,D) );
              end
            end
            if isnan(beta) || isinf(beta) || beta > 1e+6
              Vprintf( -2 ,'NaN , Inf or too large beta, reset to beta = 0 ');
              beta = 0;
            end
          end

        case 'quasinewton'
          switch P.QUASINEWTON_METHOD
            case 'bfgs'   , START('QUASI NEWTON (BFGS)');
            case 'dfp'    , START('QUASI NEWTON (DFP)');
            case 'sr1'    , START('QUASI NEWTON (SR1)'); H = eye( N );
          end
          
          H = NaN;
          if is1nan(J) || is1nan( E )
            [E,J] = F(x); 
          end
          while checkUSER_CANCEL() && check_METHOD_IT() && check_METHOD_STARTTIME() &&  check_IT()  &&  check_STARTTIME()  &&  check_E()  &&  check_dE()  &&  check_MIN_NORM_GRADIENT()
            if ~METHOD_IT || ( ~isinf( P.QUASINEWTON_RESTART ) &&  ~mod( METHOD_IT , P.QUASINEWTON_RESTART ) )
              Vprintf( -2 ,'Restarting QuasiNewton');
              inv_H = eye( N );
            else
              y = J - J0;
              if METHOD_IT == 1         %%pag 143 Nocedal & Wright
                inv_H = alpha* y.'*Dn/( y.' * y )*inv_H;
              end

              switch P.QUASINEWTON_METHOD
                case 'bfgs'
                  H = eye( N ) - ( y * Dn.' )/( y.' * Dn );
                  inv_H = H.' * inv_H * H + alpha*( Dn * Dn.' )/( y.' * Dn );
                case 'dfp'
                  inv_H = inv_H + alpha*( Dn * Dn.' )/( y.' * Dn ) - ( inv_H * y * y.' * inv_H )/( y.' * inv_H * y );
                case 'sr1'
                  H = H + ( ( y - alpha * H * Dn ) * ( y - alpha * H * Dn ).' )/( alpha * ( y - alpha * H * Dn ).' *  Dn );
                  inv_H = EnsurePositiveDefinite( H );
              end
            end

            if any( isnan(inv_H(:)) ) || any( isinf(inv_H(:)) )
              Vprintf( -2 ,'NaN or Inf in inv_H, reset to Id');
              inv_H = eye( N ); 
            end

            Dn = - inv_H * J;
            Dn = Dn/norm( Dn );
            Dn = FILTER( Dn );
            SLOPE0 = INNER( J , Dn );

            [ alpha , M ] = LINE_SEARCH();
            if  M < E  &&  alpha > 0
              x = x + Dn*alpha;
              J0 = J;
              [E,J] = F(x);
              if ~isequal(E,M), Vwarning('New Energy value different to the minimum computed at LINESEARCH'); end
              updateX( 'alpha:' , alpha , 'norm(J):',norm(J) );
            else
              Vprintf( 1 , 'no descent direction\n' );
              break;
            end
          end

        case 'lmbfgs',  START('LIMITED-MEMORY QUASI NEWTON BFGS');
          H = NaN;
          if is1nan(J) || is1nan( E )
            [E,J] = F(x); 
          end
          while checkUSER_CANCEL() && check_METHOD_IT() && check_METHOD_STARTTIME() &&  check_IT()  &&  check_STARTTIME()  &&  check_E()  &&  check_dE()  &&  check_MIN_NORM_GRADIENT()
            
            if ~METHOD_IT
              y = {}; s = {}; rho = [];
              D = -J;
            else
              D = J; a = zeros( numel(y) , 1 );
              for l = numel(y):-1:1
                a(l) = INNER( s{l},D ) * rho(l);
                D = D - a(l) * y{l};
              end
%               D = D* ( s{end}'*y{end} )/( y{end}'*y{end} );
              for l = 1:numel(y)
                D = D + s{l}*( a(l) - rho(l) * INNER( y{l},D ) );
              end
              D = -D;
            end
            if any( isnan(D) ) || any( isinf(D) )
              Vprintf( -2 ,'Nan or Inf in D, break!!');
              break; 
            end

            Dn     = D/norm(D);
            Dn     = FILTER( Dn );
            SLOPE0 = INNER( J , Dn );
            
            [ alpha , M ] = LINE_SEARCH();
            if  M < E  &&  alpha > 0
              x = x + Dn*alpha;
              J0     = J;
              [E,J] = F(x);
              try
                if numel(y) >= P.LMBFGS_ORDER
                  y(1) = []; s(1) = []; rho(1) = [];
                  Vprintf(-2,' reducing the stored gradients and secants to %d',numel(y) );
                end
                y{end+1} = J - J0;
                s{end+1} = Dn*alpha;
                rho(end+1) = 1 / INNER( y{end} , s{end} );
                Vprintf(-2,'updating stored gradients and secants ( %d )',numel(y) );
              catch LE
                if strcmp( LE.identifier , 'MATLAB:nomem' )
                  Vprintf(-2 ,'Insuficient Memory to store more steps... ');
                  P.LMBFGS_ORDER = numel(y)-1;
                end
                
                y(1) = []; s(1) = []; rho(1) = [];
                Vprintf(-2 ,'reducing the stored gradients and secants to %d', numel(y) );
                try
                  y{end+1} = J - J0;
                  s{end+1} = Dn*alpha;
                  rho(end+1) = 1 / INNER( y{end},s{end} );
                  Vprintf(-2,'updating stored gradients and secants ( %d )',numel(y) );
                end
              end

              if ~isequal(E,M), Vwarning('New Energy value different to the minimum computed at LINESEARCH'); end
              updateX( 'alpha:' , alpha , 'norm(J):',norm(J) );
            else
              Vprintf( 1 , 'no descent direction\n' );
              break;
            end
          end

        case 'newton'
          switch P.NEWTON_HESSIAN_MODIFICATION
            case 'none',  START('NEWTON without HESSIAN MODIFICATION');
            case 'eigs',  START('NEWTON with EIGENVALUE MODIFICATION');
            case 'eye' ,  START('NEWTON with EYE MODIFICATION');
            case 'add' ,  START('NEWTON with ADD MODIFICATION');
          end
          if is1nan(J) || is1nan( E ) || is1nan( H )
            [E,J,H] = F(x); 
          end
          hessianComputed = true;
          while checkUSER_CANCEL() && check_METHOD_IT() && check_METHOD_STARTTIME() &&  check_IT()  &&  check_STARTTIME()  &&  check_E()  &&  check_dE()  &&  check_MIN_NORM_GRADIENT()
            if hessianComputed
              inv_H  = EnsurePositiveDefinite( H );
              hessianComputed = false;
            end
            
            setWarning('off','MATLAB:nearlySingularMatrix');
            Dn     = - ( inv_H * J );
            restoreWarning(  'MATLAB:nearlySingularMatrix');
            
            
            Dn     = FILTER( Dn );
            SLOPE0 = INNER( J , Dn );
            
            [ alpha , M ] = LINE_SEARCH();
            if  M < E  &&  alpha > 0
              x = x + Dn*alpha;
              if  ( P.NEWTON_SHAMANSKII_STEP > 0  &&  ~mod( METHOD_IT , P.NEWTON_SHAMANSKII_STEP ) ) || P.NEWTON_SHAMANSKII_STEP <= 0
                [E,J,H]  = F(x);
                hessianComputed = true;
              else
                [E,J]  = F(x);
              end
              if ~isequal(E,M), Vwarning('New Energy value different to the minimum computed at LINESEARCH'); end
              updateX( 'alpha:' , alpha , 'norm(J):',norm(J) );
            else
              Vprintf( 1 , 'no descent direction\n' );
              break;
            end

          end

        case 'partan',  START('PARTAN SEARCH');
          J = NaN; H = NaN;
          if is1nan( E )
            E = F(x); 
          end
          while checkUSER_CANCEL() && check_METHOD_IT() && check_METHOD_STARTTIME() &&  check_IT()  &&  check_STARTTIME()  &&  check_E()  &&  check_dE()
%             x_before = x;
            
            x_2 = x - [1;x(1:end-1)*0];
            x_1 = x - [x(1:end-1);1e10];
            
            
            Dn = x - x_2;
            Dn = Dn/norm(Dn);
            
            [ alpha , M ] = LINE_SEARCH();
            if  M < E  &&  alpha > 0
              x0 = x;
              x  = x + Dn*alpha;
              E = M;
            else
              Dn(c) = -1;
             [ alpha , M ] = LINE_SEARCH();
              if  M < E  &&  alpha > 0
                x0 = x;
                x = x + Dn*alpha;
                E = M;
              else
                break;
              end
            end


            Dn = - Dn;

            Vprintf(2,'dir:  -%6d  ',c);
            [ alpha , M ] = LINE_SEARCH();
            Vprintf(2,'\b');
            if  M < E  &&  alpha > 0
              x = x + Dn*alpha;
              E = F(x);
              if ~isequal(E,M), Vwarning('New Energy value different to the minimum computed at LINESEARCH'); end
              updateX( 'negative','alpha:',alpha , sprintf('dir: -%d',c) );
              break;
            end

            if isequal( x_before , x )
              Vprintf(1,'no descent directions\n');
              break;
            end

          end

        case 'coordinate',  START('COORDINATE SEARCH');
          J = NaN; H = NaN;
          if is1nan( E ) 
            E = F(x); 
          end
          while checkUSER_CANCEL() && check_METHOD_IT() && check_METHOD_STARTTIME() &&  check_IT()  &&  check_STARTTIME()  &&  check_E()  &&  check_dE()  
            x_before = x;
            for c = vec( P.COORDINATES_ORDER( N ) ).'
%               if any( GO.FIXED_COORDINATES == c ), continue; end
              Dn     = zeros( N , 1 );
              
              Dn(c)  = 1;
              Dn     = FILTER( Dn );
              Vprintf(2,'dir:  %6d  ',c);
              [ alpha , M ] = DOUBLE_LINE_SEARCH();
              Vprintf(2,'\b');
              if  M < E  &&  alpha ~= 0
                x = x + Dn*alpha;
                E = F(x);
                if ~isequal(E,M), Vwarning('New Energy value different to the minimum computed at LINESEARCH'); end
                updateX( 'alpha:', alpha , sprintf('dir: %d',c) );
                alpha = abs(alpha);
                break;
              end
            end

            if isequal( x_before , x )
              Vprintf(1,'no descent directions\n');
              break;
            end

          end

        case 'powell',  START('POWELL SEARCH');
          %%aca faltan tener en cuenta   GO.FIXED_COORDINATES!!
          J = NaN; H = NaN;
          if is1nan( E ) 
            E = F(x); 
          end
          Ds = eye(N);
          Ds = Ds(:,[end 1:end-1]);
          Dn = Ds(:,1);
          while checkUSER_CANCEL() && check_METHOD_IT() && check_METHOD_STARTTIME() &&  check_IT()  &&  check_STARTTIME()  &&  check_E()  &&  check_dE()
            x_before = x;
            
            Ds = [ Ds(:,2:end) Dn ];
            for c = 1:N
              Dn     = Ds(:,c);
              Dn     = FILTER( Dn );
              
              Vprintf(2,'dir:  %6d  ',c);
              [ alpha , M ] = DOUBLE_LINE_SEARCH();
              Vprintf(2,'\b');
              if  M < E  &&  alpha ~= 0
                x = x + Dn*alpha;
                E = F(x);
                if ~isequal(E,M), Vwarning('New Energy value different to the minimum computed at LINESEARCH'); end
                updateX( 'alpha:', alpha , sprintf('dir: %d',c) );
                alpha = abs(alpha);
              end
            end

            if isequal( x_before , x )
              Vprintf(1,'no descent directions\n');
              break;
            end
            
            Dn = x - x_before;
            Dn = Dn(:)/norm(Dn(:));
            Dn = FILTER( Dn );
            
            Vprintf(2,'dir:  powell' );
            [ alpha , M ] = DOUBLE_LINE_SEARCH();
            Vprintf(2,'\b');
            if  M < E  &&  alpha ~= 0
              x = x + Dn*alpha;
              E = F(x);
              if ~isequal(E,M), Vwarning('New Energy value different to the minimum computed at LINESEARCH'); end
              updateX( 'alpha:', alpha , sprintf('dir: %d',c) );
              alpha = abs(alpha);
            else
              Vprintf(1,'no descent direction\n');
              break;
            end
            
          end

        case 'dampednewton', START('DAMPED NEWTON');
          if is1nan(J) || is1nan( E ) || is1nan( H )
            [E,J,H] = F(x); 
          end
          lambda = max( abs( diag(H) ) ) * P.DAMPED_INITIAL_LAMBDA;
          while checkUSER_CANCEL() && check_METHOD_IT() && check_METHOD_STARTTIME() &&  check_IT()  &&  check_STARTTIME()  &&  check_E()  &&  check_dE()  &&  check_MIN_NORM_GRADIENT()
            
              if P.DAMPED_USE_DIAG_H
                Dn = - ( H + lambda * diag(diag(H)) )\J;
              else
                Dn = - ( H + lambda * eye(N)        )\J;
              end
              Dn     = FILTER( Dn );

              M = F( x + Dn );
              if  M < E
                rho = -( E - M )/( Dn'*( J +  H*Dn /2 ) );  %gain factor
                
                alpha = 1;
                LINE_SEARCH( struct('LINE_SEARCH',{{'fixed'}},'PERFORM_FIRST_SEARCH',0,'INITIAL_STEP',1,'MIN_STEP',1) , 1 , M );

                x = x + Dn;

                if  ( P.NEWTON_SHAMANSKII_STEP > 0  &&  ~mod( METHOD_IT , P.NEWTON_SHAMANSKII_STEP ) ) || P.NEWTON_SHAMANSKII_STEP <= 0
                  [E,J,H]  = F(x);
                else
                  [E,J]  = F(x);
                end

                updateX( 'lambda:' , lambda , 'norm(J):',norm(J) );

%                 lambda = lambda * P.DAMPED_REDUCE_FACTOR;
                if rho < 1
                  lambda =  lambda * ( P.DAMPED_INCREASE_FACTOR  ...
                    - 4*rho*  (      P.DAMPED_REDUCE_FACTOR +     P.DAMPED_INCREASE_FACTOR -  2 ) ...
                    - 2*rho^3*( 3  * P.DAMPED_REDUCE_FACTOR +     P.DAMPED_INCREASE_FACTOR -  4 ) ...
                    + rho^2*  ( 11 * P.DAMPED_REDUCE_FACTOR + 5 * P.DAMPED_INCREASE_FACTOR - 16 ) );
                else
                  lambda = lambda * P.DAMPED_REDUCE_FACTOR;
                end

                Vprintf(-2,'New lambda ( %g )', lambda );

              else

                lambda = lambda * P.DAMPED_INCREASE_FACTOR;
                Vprintf(-2,'Increasing lambda ( %g )', lambda );

              end

              if lambda > P.DAMPED_MAX_LAMBDA

                Vprintf(-2,'lambda ( %g ) too large, break!!', lambda );
                break;

              end

          end
          
      end
      
    end

    if isequal( METHODS{meth_id} , '$BREAK$' )
      HISTORY.END_CONDITION = 'BREAKED';
      break;
    end
    
    if isequal( x , x_old )
      HISTORY.END_CONDITION = 'CONVERGED';
      break;
    end
    
  end
  
  Vprintf(0.1, 'FIN:  %s\n', HISTORY.END_CONDITION );

  x = R(x);
  
  %%ya termino... guardando, etc...
  SAVE_FILE();

  if ~isempty( GO.VERBOSE_FILE )
    try
      fid = fopen( GO.VERBOSE_FILE , 'a' );
      fprintf( fid , '\n**************************\n');
      fprintf( fid , '\nFIN OPTIMIZATION: ' );
      fprintf( fid , datestr(now, ' HH:MM:SS.FFF       mmmm dd, yyyy  ( dddd )') );
      [t,hdms] = time( STARTTIME );
      if hdms(1)
        fprintf( fid , '   (elapsed time:  %.2f secs    [ %d days . %02d:%02d:%02d ])\n' , t , hdms(1) , hdms(2) , hdms(3) , floor( hdms(4) ) );
      else
        fprintf( fid , '   (elapsed time:  %.2f secs    [ %02d:%02d:%02d ])\n' , t , hdms(2) , hdms(3) , floor( hdms(4) ) );
      end
      fprintf( fid , '      ITERATIONS:  %d\n' , IT );
      fprintf( fid , '    FINAL ENERGY:  %s         (initial:  %s ) ( (E-E0)/E0 %%:  %s )  \n' , strtrim(sprintf('%30.20g',E)) , strtrim(sprintf('%30.20g',E0)) , strtrim(sprintf('%30.20g',(E-E0)/E0*100)) );
      fprintf( fid , '\n');
      fprintf( fid , 'Number of Function Evals :  %9u  (no compute: %9d )  [ smarts: %9u ] [ ls_smarts: %9u ]\n' , NOF_F_evals , NOF_F_no_compute , NOF_F_smarts , NOF_LS_smarts );
      fprintf( fid , 'Number of Jacobian Evals :  %9u  (numericals: %9d )  [ smarts: %9u ]\n' , NOF_J_evals , NOF_Jn_evals , NOF_J_smarts );
      fprintf( fid , 'Number of Hessian  Evals :  %9u  (numericals: %9d )  [ smarts: %9u ]\n' , NOF_H_evals , NOF_Hn_evals , NOF_H_smarts );
      fprintf( fid , '\n');
      fprintf( fid , '--------------------------------------------------------------\n');
      fprintf( fid , '**************************************************************\n');
      fprintf( fid , '\n\n');
      fclose(fid);
    catch, try, fclose(fid); end; end
  end

  
  if GO.PLOT
    xl = get( PLOT.hLine, 'XData' ); xl = xl(end)+0.2;
    line( 'Parent', PLOT.hAxe , 'YData',[-1e10 1e-30 1e10] , 'XData',[xl xl xl] ,'YLimInclude','off','Color',[0 1 1],'linewidth',3 );

    text( 'Parent', PLOT.hAxe ,  ...
      'BackgroundColor',  [0.8 0.8 1] ,  ...
      'Clipping',  'off' ,  ...
      'Color',  [0 0 0] ,  ...
      'EdgeColor', 'none' ,... %[1 0 0] ,  ...
      'FontAngle',  'normal' ,  ...
      'FontName',  'Helvetica' ,  ...
      'FontSize',  10 ,  ...
      'FontUnits',  'pixels' ,  ...
      'FontWeight',  'bold' ,  ...
      'HorizontalAlignment',  'left' ,  ...
      'LineStyle',  '-' ,  ...
      'LineWidth',  0.1 ,  ...
      'Margin',  1e-8 ,  ...
      'Rotation',  0 ,  ...
      'VerticalAlignment',  'bottom' ,  ...
      'Visible',  'on' ,  ...
      'XLimInclude',  'off' ,  ...
      'YLimInclude',  'on' ,  ...
      'rotation',90 ,...
      'Position',  [xl E 0] ,  ...
      'String',  sprintf( 'FINISH: %s' , HISTORY.END_CONDITION ) );
  end
  if GO.CLOSE_PLOT_AT_END
    try, delete( PLOT.hFig ); end
  else
    try, linkaxes( [ PLOT.hAxe  PLOT.hAxeTime ] , 'x' ); end
  end   
  try,   rmappdata( 0 , 'OPTIMIZE_START' ); end
  

  %%Experimental!!!
  GO.PLOT_LANDSCAPE_AT_END = 0;
  if GO.PLOT_LANDSCAPE_AT_END
    PLOT.hFig = figure;
    PLOT.hAxe = axes;
    l = line( 'Parent' , PLOT.hAxe , 'XData' , 0 , 'YData' , E ,'marker' , '+' , 'Color' ,[1 0 1] );
    
    D = 1; LSx = [];
    for kk = 1:5, LSx = [ LSx , linspace( -D , D , 21 ) ]; D = D/2; end
    LSx = unique(LSx); LSy = LSx*0;
    
    [E,J] = F(x);
    if norm(J) == 0, J = zeros(N,1); J(1) = 1; end 
    J = J/norm(J);
    J = [ J  null( J' ) ];
    for d = 1:size(J,2)
      Dn = J(:,d)*1e-14;
      for i = 1:numel(LSx), LSy(i) = F( x(:) + LSx(i)*Dn ); end
      l = line( 'Parent' , PLOT.hAxe , 'XData' , LSx , 'YData' , LSy ,'userdata' , Dn );
      if d == 1, set( l , 'Color' , [1 0 0] , 'linewidth' , 2 ); end
    end
  end
  %%Experimental!!!
  
  %%END ya termino!!!
  
  catch  LE

    Vprintf( 0 , '\n\nERROR :    %s \n' , datestr(now, ' HH:MM:SS.FFF       mmmm dd, yyyy  ( dddd )')  );
    
    Vprintf( 2 , '\n\n%s\n%s\n\n' , LE.message , LE.identifier );
    stack = LE.stack;
    for i = 1:numel(stack)
      Vprintf( 2 , '%s\n'  , strrep( evalc( 'stack(i)' ) , '\' , '\\' ) );
    end

    Vprintf( -2 , '\n\n\n%s\n\n\n\n'  , disperror( LE ) );
    
    Vprintf( 2 , '\n-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n\n');

    try,   rmappdata( 0 , 'OPTIMIZE_START' ); end
    
    if  GO.CLOSE_PLOT_AT_END
      try, delete( PLOT.hFig ); end
    else
      try, linkaxes( [ PLOT.hAxe  PLOT.hAxeTime ] , 'x' ); end
    end    
    
    rethrow( LE );
  end
  
  
  


  function y = R(y),   y = reshape( y , sz ); end
  function y = vec(y), y = y(:); end
  function r = is1nan(y), r = numel(y) == 1 && isnan(y); end
  function r = isfalse(y), r = islogical(y) && numel(y) == 1 && ~y; end


  function PLOT_Y_LOG_SCALE( v )
    lostfocus( PLOT.hFig );
    if v
      set( PLOT.hAxe,'YScale','log');
    else
      set( PLOT.hAxe,'YScale','linear');
    end
    drawnow('expose');
  end
  function ShowLS( v )
    lostfocus( PLOT.hFig );
    v = round(v);
    
    set( get( PLOT.hAxeLS , 'children' ) , 'Visible','off' );

    g = findall( PLOT.hAxeLS , 'Tag' , sprintf('it: %d' , v ) );
        
    if isempty(g)
      return;
    end
    
    dispstruct( get( g , 'UserData' ) );
    
    set( g , 'Visible' ,'on' );
%     lims = objbounds( get(g,'children') );
    gc = get(g,'children');
    if iscell(gc), gc = cell2mat(gc(:)); end
    
    lims = objbounds( gc );
    
    limsx = lims(1:2); limsx = (limsx - mean(limsx) )*1.1 + mean(limsx);
    limsy = lims(3:4); limsy = (limsy - mean(limsy) )*1.1 + mean(limsy);
    
    set( PLOT.hAxeLS , 'XLim' , limsx , 'YLim',limsy );
  end
  function h = DrawInGLSline( type , varargin )
    switch lower( type )
      case 'line',  h = line( 'Parent',PLOT.hGLS,varargin{:} );
    end
    
    lims = objbounds( get(PLOT.hGLS,'children') );
    if isempty( lims ) , return; end
    
    limsx = lims(1:2); limsx = (limsx - mean(limsx) )*1.5 + mean(limsx);
    limsy = lims(3:4); limsy = (limsy - mean(limsy) )*1.5 + mean(limsy);

    if diff( limsx ) <= 2e-6, limsx = mean(limsx) + [ -1e-6 1e-6 ]; end
    if diff( limsy ) <= 2e-6, limsy = mean(limsy) + [ -1e-6 1e-6 ]; end
    
    try
      set( PLOT.hAxeLS , 'xlim', limsx ,'ylim' , limsy );
    end
    drawnow;
  end


  function updateX( varargin )
    IT  = IT+1;
    METHOD_IT = METHOD_IT + 1;
    
    ellapsed_time =  time(STARTTIME);

    
    dE = Ep-E;
    str = sprintf('[%7d.%4ds].Energy: %30.15g ( %.2e )   ' , IT , floor( ellapsed_time ), E , dE );
    
    Ep = E;
    while numel(varargin)
      switch class( varargin{1} )
        case 'char',                str = [ str  sprintf(' %s '    , varargin{1} ) ];
        case {'double','single'},   str = [ str  sprintf(' %+.3e  ', varargin{1} ) ];
      end
      varargin(1) = [];
    end
    fileprintf( GO.VERBOSE_FILE , ' ' );
    Vprintf(1,[ str '\n' ]);


    if ~isempty( PROJECT_FCN )
      x = PROJECT_FCN( x , IT );
    end

    if ~isempty( OUTPUT_FCN )
      OUTPUT_FCN( x , E , IT );
    end


    if ~GO.PLOT && ~isempty( OUTPUT_FCN )
      drawnow;
    end

    if GO.SAVE_IN_BASE, setappdata( 0 , 'OPTIMIZE_X' , x ); end
    
    if SAVE_X_PATH
      X_PATH{end+1} = x;
      if GO.SAVE_IN_BASE, setappdata( 0 , 'OPTIMIZE_X_PATH' , X_PATH ); end
    end


    if GO.SAVE_IN_BASE && SAVE_ALL_EVALS , setappdata( 0 , 'OPTIMIZE_ALL_EVALS' , ALL_EVALS ); end

    
    if SAVE_HISTORY
      try,   HISTORY.it(IT) = IT;
      catch, HISTORY.it(IT) = NaN;
      end

      try,   HISTORY.energy(IT) = E;
      catch, HISTORY.energy(IT) = NaN;
      end
    
      try,   HISTORY.method{IT} = P.METHOD;
      catch, HISTORY.method{IT} = '';
      end
    
      try,   HISTORY.beta(IT) = beta;
      catch, HISTORY.beta(IT) = NaN;
      end

      try,   HISTORY.lambda(IT) = lambda;
      catch, HISTORY.lambda(IT) = NaN;
      end

      try,   HISTORY.LineSearch{IT} = LS;
      catch, HISTORY.LineSearch{IT} = [];
      end
      
      try,   HISTORY.normG(IT) = norm( J(:) );
      catch, HISTORY.normG(IT) = NaN;
      end

      try,   HISTORY.maxnormG(IT) = max(abs( J(:) ));
      catch, HISTORY.maxnormG(IT) = NaN;
      end

      try,   HISTORY.time(IT) = ellapsed_time;
      catch, HISTORY.time(IT) = NaN;
      end

      try
        HISTORY.evals.NOF_F_evals(IT)   = NOF_F_evals;
        HISTORY.evals.NOF_J_evals(IT)   = NOF_J_evals;
        HISTORY.evals.NOF_Jn_evals(IT)  = NOF_Jn_evals;
        HISTORY.evals.NOF_H_evals(IT)   = NOF_H_evals;
        HISTORY.evals.NOF_Hn_evals(IT)  = NOF_Hn_evals;

        HISTORY.evals.NOF_F_smarts(IT)  = NOF_F_smarts;
        HISTORY.evals.NOF_J_smarts(IT)  = NOF_J_smarts;
        HISTORY.evals.NOF_H_smarts(IT)  = NOF_H_smarts;
        HISTORY.evals.NOF_LS_smarts(IT) = NOF_LS_smarts;

        HISTORY.evals.NOF_F_no_compute(IT) = NOF_F_no_compute;
      end
      
      if GO.SAVE_IN_BASE, setappdata( 0 , 'OPTIMIZE_HISTORY' , HISTORY ); end
    end    
  
    if GO.SAVE_X_EVERY >= 0  &&  IT - LAST_SAVE_IT >= GO.SAVE_X_EVERY
      SAVE_FILE();
    elseif  GO.SAVE_X_EVERY < 0  &&  time( LAST_SAVE_TIME ) >= -GO.SAVE_X_EVERY
      SAVE_FILE();
    end

    if GO.PLOT
      xl = get( PLOT.hLine , 'XData' );
      yl = get( PLOT.hLine , 'YData' );
      set( PLOT.hLine , 'XData', [ xl(:) ; xl(end)+1 ] , 'YData', [ yl(:) ; E ] );
      
      yl = get( PLOT.hLineTime , 'YData');
      set( PLOT.hLineTime , 'XData', [ xl(:) ; xl(end)+1 ] , 'YData', [ yl(:) ; ellapsed_time ] );
      set( [ PLOT.hAxe PLOT.hAxeTime ] , 'XLim' , [0 xl(end)+2] );
      drawnow;
    end

  end
  function Vwarning( varargin )
    if ~isempty( GO.VERBOSE_FILE )
      try
        fid = fopen( GO.VERBOSE_FILE , 'a' );
        fprintf( fid , '*******Warning:  ' );
        fprintf( fid , varargin{:} );
        fprintf( fid , '\n\n' );
        fclose(fid);
      catch
        try fclose(fid); end
      end
      try fclose(fid); end
    end
    warning( varargin{:} );
  end
  function Vprintf(varargin)
    if isnumeric( varargin{1} )
      level = varargin{1};
      varargin(1) = [];
    else
      level = 1;
    end

    
    str = sprintf( varargin{:} );

    while ~isempty(str) && str(1) == 8
      str(1) = [];
      if level > 0   &&   abs( GO.VERBOSE ) >= level
        if GO.VERBOSE > 0  
          fprintf( repmat( '\b' , 1 , printMARGIN(end) ) );
        end
        printMARGIN(end) = [];
      end        
    end
    
    
    if ~isempty( GO.VERBOSE_FILE )  &&  ceil( abs( GO.VERBOSE ) ) >= abs( level )  && ~isempty( str )
      str_with_spaces = [ repmat( ' ' , 1 , sum( printMARGIN ) )  str ];
      fileprintf( GO.VERBOSE_FILE , str_with_spaces );
    end
    
    if level >= 0  &&  abs( GO.VERBOSE ) >= level
      if ~isempty(str)
        if GO.VERBOSE > 0  ||  level == 0
          fprintf( '%s', str );
        end
        printMARGIN = [ printMARGIN  numel(str) ];
      end
      if ~isempty(str)  &&  str(end) == 10
        printMARGIN = 0;
      end
    end
    
  end
  function SAVE_FILE()
    if ~isempty(  GO.SAVE_X_IN_FILE  )
      if       SAVE_HISTORY  &&  SAVE_X_PATH
        save( '-v7.3' , GO.SAVE_X_IN_FILE , 'x' , 'HISTORY' , 'X_PATH' ); 
      elseif  ~SAVE_HISTORY  &&  SAVE_X_PATH
        save( '-v7.3' , GO.SAVE_X_IN_FILE , 'x' , 'X_PATH' ); 
      elseif   SAVE_HISTORY  && ~SAVE_X_PATH
        save( '-v7.3' , GO.SAVE_X_IN_FILE , 'x' , 'HISTORY' ); 
      elseif  ~SAVE_HISTORY  && ~SAVE_X_PATH
        save( '-v7.3' , GO.SAVE_X_IN_FILE , 'x' ); 
      end
      Vprintf( 0 , 'saving ... %s   [  %s  ]  ( %s )\n', GO.SAVE_X_IN_FILE , datestr(now, ' HH:MM:SS   dd/mmm/yy ') , fixname(GO.SAVE_X_IN_FILE) );
      LAST_SAVE_IT = IT;
      LAST_SAVE_TIME = clock;
    end
  end


  function START(mode)

    Vprintf(1,'Starting %s...\n',mode);
    
    fileprintf(GO.VERBOSE_FILE,'-------------------------------------\nStarting %s\n         PARAMETERS: %s\n-------------------------------------' , mode,strrep( strrep( strrep( evalc( 'P' ) , ' ' , '' ) , char(10) , '  ;  ' ) , ':' , ': ' ) );

    if GO.PLOT
      xl = get( PLOT.hLine, 'XData' ); xl = xl(end)+0.2;
      line( 'Parent', PLOT.hAxe , 'YData',[-1e10 1e-30 1e10] , 'XData',[xl xl xl] ,'YLimInclude','off','Color',[1 0 1] );
      
      if numel( METHODS ) > 1
        text( 'Parent', PLOT.hAxe ,  ...
              'BackgroundColor',  [1 1 0] ,  ...
              'Clipping',  'off' ,  ...
              'Color',  [1 0 0] ,  ...
              'EdgeColor', 'none' ,... %[1 0 0] ,  ...
              'FontAngle',  'normal' ,  ...
              'FontName',  'Helvetica' ,  ...
              'FontSize',  10 ,  ...
              'FontUnits',  'pixels' ,  ...
              'FontWeight',  'bold' ,  ...
              'HorizontalAlignment',  'left' ,  ...
              'LineStyle',  '-' ,  ...
              'LineWidth',  0.1 ,  ...
              'Margin',  1e-8 ,  ...
              'Rotation',  0 ,  ...
              'VerticalAlignment',  'bottom' ,  ...
              'Visible',  'on' ,  ...
              'XLimInclude',  'off' ,  ...
              'YLimInclude',  'on' ,  ...
              'Position',  [xl E 0] ,  ...
              'String',  mode              );
      end
      
      drawnow;
    end

    METHOD_STARTTIME = clock;
    METHOD_IT = 0;
    dE        = Inf;
  end
  function [t,dhms] = time( INIT_TIME )
    if nargin < 1
      INIT_TIME = 0;
    end
    if isscalar( INIT_TIME )
      if INIT_TIME == 0
        INIT_TIME = STARTTIME;
      elseif INIT_TIME == 1
        INIT_TIME = METHOD_STARTTIME;
      else
        warning(' INTERNAL!!!! , error in  INIT_TIME ' );
      end
    end
    t = etime( clock , INIT_TIME );
    if nargout > 1
      s = t;
      d = floor( s/3600/24 );
      s = s - d*3600*24;
      h = floor( s/3600 );
      s = s - h*3600;
      m = floor( s/60 );
      s = s - m*60;
      dhms=[d h m s];
    end
  end
  function v = checkUSER_CANCEL( )
    v = 1;
  end
  function v = check_METHOD_IT( )
    if METHOD_IT < P.MAX_ITERATIONS_PER_METHOD
      v = true;
    else
      v = false;
      Vprintf( 1 , 'MAX_ITERATIONS_PER_METHOD ( %d )  reached\n' , P.MAX_ITERATIONS_PER_METHOD );
    end
  end
  function v = check_METHOD_STARTTIME( )
    if time(METHOD_STARTTIME) < P.MAX_TIME_PER_METHOD
      v = true;
    else
      v = false;
      Vprintf( 1 , 'MAX_TIME_PER_METHOD ( %f secs )  reached\n' , P.MAX_TIME_PER_METHOD );
    end
  end
  function v = check_dE( )
    if dE > P.MIN_ENERGY_DELTA  ||  dE < 0
      v = true;
      if isinf( dE ), return; end
      small_dE_times = 0;
    else
      v = false;
      small_dE_times = small_dE_times + 1;
      Vprintf( 1 , 'MIN_ENERGY_DELTA ( %g )  reached ( %d )\n' , P.MIN_ENERGY_DELTA , small_dE_times );
      dE = Inf;
    end
  end
  function v = check_IT( )
    if IT < GO.MAX_ITERATIONS
      v = true;
    else
      v = false;
      if GO.MAX_ITERATIONS > 0
        Vprintf( 1 , 'MAX_ITERATIONS ( %d )  reached\n' , GO.MAX_ITERATIONS );
        HISTORY.END_CONDITION = sprintf('MAX_ITERATIONS (%d) reached',GO.MAX_ITERATIONS);
      else
        GO.MAX_ITERATIONS = -1;
      end
    end
  end
  function v = check_STARTTIME( )
    if time( STARTTIME ) < GO.MAX_TIME
      v = true;
    else
      v = false;
      if GO.MAX_TIME > 0
        Vprintf( 1 , 'MAX_TIME ( %f secs )  reached\n' , GO.MAX_TIME );
        HISTORY.END_CONDITION = sprintf('MAX_TIME ( %g ) seconds reached',GO.MAX_TIME);
      else
        GO.MAX_TIME = -1;
      end
    end
  end
  function v = check_E( )
    if E > GO.MIN_ENERGY
      v = true;
    else
      v = false;
      Vprintf( 1 , 'MIN_ENERGY ( %g )  reached\n' , GO.MIN_ENERGY );
      HISTORY.END_CONDITION = sprintf('MIN_ENERGY ( %g ) reached',GO.MIN_ENERGY);
    end
  end
  function v = check_MIN_NORM_GRADIENT( )
    if norm(J) > P.MIN_NORM_GRADIENT
      v = true;
    else
      v = false;
      Vprintf( 1 ,'MIN_GRADIENT_NORM : %g ( < %g )\n' , norm(J) , P.MIN_NORM_GRADIENT );
      HISTORY.END_CONDITION = sprintf('MIN_GRADIENT_NORM : %g ( < %g )\n' , norm(J) , P.MIN_NORM_GRADIENT );
    end
  end


  function hn = HostName()
  %   hn = getenv('HOSTNAME');
    [a,hn] = system('hostname -s');
  %   hn = hn( 1:find(hn=='.',1,'first')-1 );
    hn = strtrim(hn);
  end

  function d = FILTER(d)
    if ~isempty( P.FILTER )
      try
        Dn = P.FILTER( Dn );
      catch
        Vwarning('error at filtering'); 
      end
    end
  end

  function inv_H = EnsurePositiveDefinite( H )
    switch P.NEWTON_HESSIAN_MODIFICATION
      case 'none'
        inv_H = inv(H);
      case 'add'
        PD = 0; try, chol(H); PD = 1; end
        while ~PD
          H = H + abs(max(H(:)))*eye(N);
          try, chol(H); PD = 1; end
        end
        inv_H = inv(H);
      case 'eigs'
        try
          chol(H);
          try
            inv_H = inverse(H);
          catch
            
            setWarning('off','MATLAB:nearlySingularMatrix');
            inv_H = inv(H);
            restoreWarning(  'MATLAB:nearlySingularMatrix');
            
          end
        catch
          [eig_vecs,eig_vals] = eig(H);
          eig_vals = diag( eig_vals );
          min_eig_val = min( 1e-8 , eps(max(eig_vals))*1e+5 );
          eig_vals( eig_vals < min_eig_val ) = min_eig_val;
          inv_H = eig_vecs*diag( 1./eig_vals )*eig_vecs';
        end
      case 'eye'
        try
          chol(H);

          setWarning('off','MATLAB:nearlySingularMatrix');
          inv_H = inv(H);
          restoreWarning(  'MATLAB:nearlySingularMatrix');

        catch
          eig_vals = eig(H);
          min_eig_val = min( 1e-8 , eps(max(eig_vals))*1e+5 );

          setWarning('off','MATLAB:nearlySingularMatrix');
          inv_H = inv( H + max( 0 ,min_eig_val - min(eig_vals) )*eye(size(H)) );
          restoreWarning(  'MATLAB:nearlySingularMatrix');
        end
    end

    try,
      if any( isnan( inv_H(:) ) ) || any( isinf( inv_H(:) ) )
        inv_H = eye( N );
      end
    end
  end



  function els = FLS( xx )
    
    xxid = find( FLS_a == xx , 1 );
    if isempty( xxid )
    
      els   = F( x + xx*Dn );

      FLS_e = [ FLS_e   els  ];
      FLS_a = [ FLS_a   xx   ];

    else
      
      els   = FLS_e( xxid );
      NOF_LS_smarts = NOF_LS_smarts + 1;
      Vprintf(-4,'NOF: (lsf): %d , %d , %d , %d , %d , %d , %d , %d , %d , %d',NOF_F_evals,NOF_J_evals,NOF_H_evals,NOF_Jn_evals,NOF_Hn_evals,NOF_F_smarts,NOF_F_no_compute,NOF_J_smarts,NOF_H_smarts,NOF_LS_smarts);

    end
    
  end
  function [a,M] = LINE_SEARCH( LSP , a_fixed , M_fixed , varargin )
    if ~GO.LS_PLOT && nargin == 3
      a = a_fixed;
      M = M_fixed;
      return;
    end
    
    FLS_a = 0;
    FLS_e = E;
    
    if nargin == 0, LSP = P; end
    
    if GO.LS_PLOT
%       allLS = setdiff( findall( PLOT.hAxeLS , 'type','hggroup') , findall( PLOT.hAxeLS , 'Tag' , sprintf('it: %d' , IT-1 ) ) );
      allLS = findall( PLOT.hAxeLS , 'type','hggroup');
      set( allLS , 'Visible','off' );
      PLOT.hGLS = hggroup('Parent', PLOT.hAxeLS , 'Tag' , sprintf('it: %d',IT) , 'UserData' , [] );

      set( PLOT.LS_control , 'Min', 0 ,'Max', IT + 1 , 'SliderStep' , [1 1]/(IT+1) , 'value', IT );

      DrawInGLSline('line', 'XData' , PLOT.LS_X + [0 0 0], 'YData' , [-1e10 E 1e10] , 'line',':'   ,'YLimInclude','off' ,'color',[0  1  0] , 'color',[0 0.6 0] ,'marker','o','markerfacecolor',[0 1 0] );
    end
    
    DO_BRACKETING = false;
    
    lsorder = 1:numel(LSP.LINE_SEARCH);
    if METHOD_IT < LSP.PERFORM_FIRST_SEARCH
      Vprintf(-2,'PERFORM_FIRST_SEARCH (  METHOD_IT: %d  < %d )', METHOD_IT , LSP.PERFORM_FIRST_SEARCH );
      lsorder = lsorder([end 1:end-1]);
      DO_BRACKETING = true;
    end
    
    if LSP.INITIAL_STEP > 0
      xf = LSP.INITIAL_STEP;
    elseif LSP.INITIAL_STEP == 0
      xf = max( LSP.MIN_STEP , LSP.BRACKETING_MIN );
      DO_BRACKETING = true;
    else
      if isnan( alpha ) || alpha <= 0
        xf = max( LSP.MIN_STEP , LSP.BRACKETING_MIN );
        DO_BRACKETING = true;
      else
        xf = abs( LSP.INITIAL_STEP )*alpha;
      end
    end
    
    if LSP.MIN_STEP == 0
      LSP.MIN_STEP = eps( min( abs( x./Dn ) ) ) / 4;
    end
    
    if DO_BRACKETING
      Vprintf(2,'bracketing... ' );
      Vprintf( -2 , '             xf (             E/M )');
      Vprintf(2,'%15g ( %15g )', 0 , E );

      if GO.LS_PLOT, DrawInGLSline('line','XData',0 + PLOT.LS_X , 'YData' , E , 'line','none','marker','>','markerFaceColor',[0 0 1] ); end
      
      M  = E;
      brack_it = 0;
      while 1
        brack_it = brack_it + 1;
        
        Mb = FLS( xf );
        Vprintf(2,'\b%15g ( %15g )',xf,Mb );

        if GO.LS_PLOT, DrawInGLSline('line', 'XData' , xf + PLOT.LS_X , 'YData' , Mb , 'line','none','marker','>','markerFaceColor',[0 0 1] ); end
        
        if Mb > M, break; end
        M = Mb;
        if brack_it > LSP.BRACKETING_ITS, break; end

        xf = xf * LSP.BRACKETING_RHO;
        if  xf >= LSP.BRACKETING_MAX,     break; end
      end
      Vprintf( 2,'\b\b')
      Vprintf(-2,'... end bracketing');

      if GO.LS_PLOT, DrawInGLSline('line','XData',xf + PLOT.LS_X , 'YData' , Mb , 'line','none','marker','>','markerFaceColor',[1 0 1] ); end
    end
    
    
    xf = max( xf , LSP.MIN_STEP );
    
    if GO.LS_PLOT
      DrawInGLSline('line', 'XData' , xf + PLOT.LS_X + [0 0 0], 'YData' , [-1e10 E 1e10] , 'line',':','YLimInclude','off' ,'color',[0 0 1] ,'marker','o','markerfacecolor',[0 0 1]);
    end
    
    
    smartLS_already_executed = false;
    for lsm = vec( uniquens( LSP.LINE_SEARCH(lsorder) ) )'
      LS = [];
      LS.IT = IT;
      LS.METHOD = lsm{1};
      LS.PARAMETERS = LSP;
      Vprintf(2,'ls: %s ...',LS.METHOD);
      LS.xf = xf;
      
      switch lower( LS.METHOD )
        case 'smart'
          [a,M] = smartLS( LSP );
          smartLS_already_executed = true;
        case 'fixed'
          switch nargin
            case {0,1},  a = xf;       M = FLS( a );
            case 2    ,  a = a_fixed;  M = FLS( a );
            case 3    ,  a = a_fixed;  M = M_fixed;
          end

        case 'golden'
          [a,M] = goldensearch( LSP , 0 , xf );
        case 'backtracking'
          [a,M] = backtracking( LSP , xf );
        case 'quadratic'
          [a,M] = quadratic( LSP , xf );
        case 'exhaustive'
          [a,M] = exhaustive( LSP , xf );
        otherwise
          error( 'error doing LINE_SEARCH with %s', LS.METHOD );
      end
      
      LS.alpha = a;
      LS.energy = M;
      Vprintf(2,'\b... end %s' , LS.METHOD );
      Vprintf(2,'\b' );
      
      Vprintf(-3,'%s', strrep( strrep( strrep( evalc( 'LS' ) , ' ' , '' ) , char(10) , '  ;  ' ) , ':' , ': ') );
      
      if M < E, break; end
    end
    Vprintf(-2,'LS: %s  ... a:  %g   ...  M: %g', LS.METHOD , a , M );
    
    if ~smartLS_already_executed
      Vprintf(2,'smart LS');
      smartLS( LSP );
      Vprintf(2,'\b');
    end
    
    [mm,mmid] = min( FLS_e ); 
    if mm < M
      Vprintf(-2,'STORED LS values give a smaller result, mm: %g  at a: %g   instead of M: %g at a: %g', mm , FLS_a(mmid) , M , a );
      M = mm;
      a = FLS_a(mmid);
    end
    
    
    if GO.LS_PLOT
      LS.FLS_a = FLS_a;
      LS.FLS_e = FLS_e;

      DrawInGLSline('line', 'XData' ,  a + PLOT.LS_X , 'YData' , M , 'line','none','marker','o','markerFaceColor',[1 0 0] ,'markersize', 10 );

      LSx = [0 max(a,xf) ]; LSx = 1.05*(LSx - mean(LSx) ) + mean(LSx);
      LSd = LSx(2) - LSx(1);
      LSx = [ 0  linspace( LSx(1) , LSx(2) , 50 ) ];
      for kk = 1:6
        LSd = LSd/2;
        LSx = unique( [ LSx  linspace(  a - LSd , a + LSd , 21 )  FLS_a ] );
      end
      
      LSy = LSx*0;
      for i = 1:numel( LSx )
        LSy(i) = FLS( LSx(i) );
      end
      
      DrawInGLSline('line', 'XData' , PLOT.LS_X + LSx , 'YData' , LSy , 'line','-' , 'linewidth',2 , 'color',[1 0 0]);

      set( PLOT.hGLS , 'UserData' , LS )
      PLOT.LS_X = PLOT.LS_X + a;
    end

  end
  function [a,M] = DOUBLE_LINE_SEARCH( varargin )
    Vprintf(2,'search on positive:');
    [a,M] = LINE_SEARCH( varargin{:} );
    Vprintf(2,'\b');
    if M < E
      return;
    end
    
    Vprintf(2,'search on negative:');
    Dn = -Dn;
    [a,M] = LINE_SEARCH( varargin{:} );
    Vprintf(2,'\b');
    
    Dn = -Dn;
    a  = -a;
  end
  function [a,M] = goldensearch( LSP , x1 , x4 )
    if LSP.LS_GOLDEN_MIN_DIFF < eps(0) * 100 , LSP.LS_GOLDEN_MIN_DIFF = eps(0) * 100; end
    
    LS.ITGS = 0;
    Vprintf(-2 ,'ITGS - [              x1                 x4 ] (     x4-x1 )' );
    Vprintf( 2 ,'%4d - [ %15g    %15g ] ( %.2e )',LS.ITGS,x1,x4,x4-x1);
    
    e1 = FLS( x1 );
    
    e4 = FLS( x4 );
    GR = (3 - sqrt(5))/2;
    DX = x4-x1;
    x2 = x1 + DX*GR; e2 = FLS( x2 );
    x3 = x4 - DX*GR; e3 = FLS( x3 );

    [M,id] = min([e1 e2 e3 e4]);
    a = [x1 x2 x3 x4]; a = a(id);

    while true
      
      if GO.LS_PLOT, 
        DrawInGLSline('line' ...
        , 'XData' , PLOT.LS_X + [x1 x1 x2 x3 x4 x4]  ...
        , 'YData' , [e1 (min([e1 e2 e3 e4])+max([e1 e2 e3 e4]))/2*[1 1 1 1] e4] ...
        , 'linestyle','-','marker','+' , 'color' , [.3 .3 .4] ); end
      
      if LS.ITGS > LSP.LS_GOLDEN_MAX_ITS  ,  LS.STOP = sprintf('LS.ITGS(%d) > P.LS_GOLDEN_MAX_ITS(%d)',LS.ITGS,LSP.LS_GOLDEN_MAX_ITS);   break; end
      if DX < LSP.LS_GOLDEN_MIN_DIFF      ,  LS.STOP = sprintf('DX(%g) < P.LS_GOLDEN_MIN_DIFF(%g)'    ,DX,LSP.LS_GOLDEN_MIN_DIFF);       break; end
      if x4 <= LSP.MIN_STEP               ,  LS.STOP = sprintf('x4(%g) <= P.MIN_STEP(%g)'             ,x4,LSP.MIN_STEP);                 break; end
      if x4 < LSP.LS_GOLDEN_MIN_ALPHA     ,  LS.STOP = sprintf('x4(%g) < P.LS_GOLDEN_MIN_ALPHA(%g)'   ,x4,LSP.LS_GOLDEN_MIN_ALPHA);      break; end
      if x4 > LSP.LS_GOLDEN_MAX_ALPHA     ,  LS.STOP = sprintf('x4(%g) > P.LS_GOLDEN_MAX_ALPHA(%g)'   ,x4,LSP.LS_GOLDEN_MAX_ALPHA);      break; end        
      if DX < eps(x1)*100                 ,  LS.STOP = sprintf('DX(%g) < eps(x1)*100'                 ,DX);                              break; end
      if isequal(x+x1*Dn , x+x4*Dn)       ,  LS.STOP = sprintf('isequal( x+x1(%g)*Dn , x+x4(%g)*Dn ),',x1,x4);                           break; end
      
      LS.ITGS = LS.ITGS+1;
      [M,id] = min([e1 e2 e3 e4]);
      a = [x1 x2 x3 x4]; a = a(id);
      if id == 1 || id == 2
        x4 = x3; e4 = e3;
        x3 = x2; e3 = e2;
        DX = x4-x1;
        x2 = x1 + DX*GR; 
        if isequal( x + x1*Dn , x + x2*Dn ) || isequal( x + x2*Dn , x + x3*Dn )
          LS.STOP = sprintf('new x2(%g)  equal to x1(%g) or x3(%g)',x2,x1,x3);
          break;
        end
        e2 = FLS( x2 );
      elseif id == 3  %|| id == 4
        x1 = x2; e1 = e2;
        x2 = x3; e2 = e3;
        DX = x4-x1;
        x3 = x4 - DX*GR; 
        if isequal( x + x2*Dn , x + x3*Dn ) || isequal( x + x3*Dn , x + x4*Dn )
          LS.STOP = sprintf('new x3(%g)  equal to x2(%g) or x4(%g)',x3,x2,x4);
          break;
        end
        e3 = FLS( x3 );
      elseif id == 4
        Vprintf(-2,'minimun in x4, increasing range');
        x2 = x3; e2 = e3;
        x3 = x4; e3 = e4;
        x4 = x1 + (x4-x1)/(1-GR); 
        if isequal( x + x3*Dn , x + x4*Dn ),
          LS.STOP = sprintf('new x4(%g)  equal to x3(%g)',x4,x3);
          break;
        end
        e4 = FLS( x4 );
        DX = x4-x1;
      end
      Vprintf( 2 ,'\b%4d - [ %15g    %15g ] ( %.2e )',LS.ITGS,x1,x4,DX);
    end
    Vprintf( 2 , '\b');
    Vprintf(-2 , 'LS.STOP:  %s', LS.STOP );
  end
  function [a,M] = backtracking( LSP , a )
    if isequal( x , x+a*Dn )
      LS.STOP = sprintf('x == x+a(%g)*Dn',a);
      M = E;
      Vprintf(-2 , 'LS.STOP:  %s', LS.STOP );
      return;  
    end
    
    M = FLS( a );
    LS.BTIT = 0;
    Vprintf(-2,'LS.BTIT -               a  (               M )');
    Vprintf(2,'   %4d - %15g  ( %15g )',LS.BTIT,a,M);
    if GO.LS_PLOT, DrawInGLSline('line', 'XData' , a + PLOT.LS_X , 'YData' , M , 'line','none','marker','<','markersize',8 ); end
    
    
    while true

      if  M < E                                   , LS.STOP = sprintf('M(%g) < E(%g)'                              , M , E );                              break;  end
      if LS.BTIT > LSP.LS_BACKTRACKING_MAX_ITS    , LS.STOP = sprintf('LS.BTIT(%d) > P.LS_BACKTRACKING_MAX_ITS(%d)',LS.BTIT,LSP.LS_BACKTRACKING_MAX_ITS);  break;  end
      if a < LSP.LS_BACKTRACKING_MIN_ALPHA        , LS.STOP = sprintf('a(%g) < P.LS_BACKTRACKING_MIN_ALPHA(%g)'    ,a,LSP.LS_BACKTRACKING_MIN_ALPHA);      break;  end
      if a <= LSP.MIN_STEP                        , LS.STOP = sprintf('a(%g) <= P.MIN_STEP(%g)'                    ,a,LSP.MIN_STEP);                       break;  end
      
      LS.BTIT = LS.BTIT + 1;

      ap = a;
      a = a * LSP.LS_BACKTRACKING_REDUCE_FACTOR;
      if a <= LSP.MIN_STEP                        , LS.STOP = sprintf('a(%g) <= P.MIN_STEP(%g)',a,LSP.MIN_STEP); a=ap; break;  end

      M = FLS( a );

      Vprintf(2,'\b   %4d - %15g  ( %15g )',LS.BTIT,a,M);
      if GO.LS_PLOT, DrawInGLSline('line', 'XData' , a + PLOT.LS_X , 'YData' , M , 'line','none','marker','<','markersize',8 ); end

    end
    Vprintf(2,'\b');

    if M<E && LSP.LS_BACKTRACKING_EXTRA_TESTS > 0

      Vprintf(-2 , 'first bactracking stop:  %s', LS.STOP );
      Vprintf(2 , 'trying %d extra test', LSP.LS_BACKTRACKING_EXTRA_TESTS );

      for ex = 1:LSP.LS_BACKTRACKING_EXTRA_TESTS
        LS.BTIT = LS.BTIT + 1;

        ap = a;
        a = a * LSP.LS_BACKTRACKING_REDUCE_FACTOR;
        if a <= LSP.MIN_STEP                        , LS.STOP = sprintf('a(%g) <= P.MIN_STEP(%g)',a,LSP.MIN_STEP); a=ap; break;  end

        Me = FLS( a );

        Vprintf(2,'\b   %4d - %15g  ( %15g ) (extra reducing)',LS.BTIT,a,Me);
        if GO.LS_PLOT, DrawInGLSline('line', 'XData' , a + PLOT.LS_X , 'YData' , Me , 'line','none','marker','<','markersize',8 , 'markerfacecolor',[.5 .5 .5]); end

        if    Me >= M , LS.STOP = sprintf('Me(%g) >= M(%g)', Me , M ); a = ap; break;
        else          , M = Me;
        end
      end
      
      
      if ex == 1 && LSP.LS_BACKTRACKING_EXTRA_TESTS > 1 && LS.BTIT <= 1
        while ex < LSP.LS_BACKTRACKING_EXTRA_TESTS
          ex = ex + 1;
          LS.BTIT = LS.BTIT + 1;

          ap = a;
          a = a / LSP.LS_BACKTRACKING_REDUCE_FACTOR;

          Me = FLS( a );

          Vprintf(2,'\b   %4d - %15g  ( %15g ) (extra enlarge)',LS.BTIT,a,Me);
          if GO.LS_PLOT, DrawInGLSline('line', 'XData' , a + PLOT.LS_X , 'YData' , Me , 'line','none','marker','>','markersize',8 , 'markerfacecolor',[.5 .5 .9]); end

          if    Me >= M , LS.STOP = sprintf('Me(%g) >= M(%g)', Me , M ); a = ap; break;
          else          , M = Me;
          end

        end
      end
    
      Vprintf( 2 , '\b');
    end
      
    Vprintf(-2 , 'LS.STOP:  %s', LS.STOP );
    
  end
  function [a,M] = quadratic( LSP , x2 )
    function [q,A,B] = FindQ(x0,x1,x2,f0,f1,f2)
      A = ( (f0-f2)/(x0-x2) + (f2-f1)/(x1-x2) )/( x0 - x1 );
      B = ( f2*(x1^2-x0^2) + f1*(x0^2-x2^2) + f0*(x2^2-x1^2) ) / ( (x0 - x1)*(x0 - x2)*(x1 - x2) );
      q = - 0.5 * B / A;
    end

    x0 = 0;
    e0 = E;
    e2 = FLS( x2 );
    
    x1 = (x2+x0)*0.5;
    e1 = FLS( x1 );
    
    [M,id] = min([e0 e1 e2]);
    a = [x0 x1 x2]; a = a(id);
    
    LS.ITGS = 0;
    Vprintf( 2,'(       x0        x1        x2 ).(       e0        e1        e2 ) ->  a');
%     Vprintf( 2,'( %8g  %8g  %8g ).( %8g  %8g  %8g ) -> ',x0,x1,x2,e0,e1,e2);
    while 1
      if LS.ITGS >= LSP.LS_QUADRATIC_MAX_ITS ,  LS.STOP = sprintf('LS.ITGS(%d) >= P.LS_QUADRATIC_MAX_ITS(%d)' ,LS.ITGS,LSP.LS_QUADRATIC_MAX_ITS);   break; end
      if x2 <= LSP.MIN_STEP                  ,  LS.STOP = sprintf('x2(%g) <= P.MIN_STEP(%g)'                 ,x2,LSP.MIN_STEP);                    break; end
      if x2 > LSP.LS_QUADRATIC_MAX_ALPHA     ,  LS.STOP = sprintf('x2(%g) > P.LS_QUADRATIC_MAX_ALPHA(%g)'    ,x2,LSP.LS_QUADRATIC_MAX_ALPHA);      break; end        
      if isequal( x , x+x2*Dn )              ,  LS.STOP = sprintf('isequal( x , x+x2(%g)*Dn )'               ,x2);                                 break; end
      if isequal( x0 ,x2 )                   ,  LS.STOP = sprintf('isequal( x0(%g) , x2(%g) )'               ,x0,x2);                              break; end
      if isequal( x0 ,x1 )                   ,  LS.STOP = sprintf('isequal( x0(%g) , x1(%g) )'               ,x0,x1);                              break; end
      if isequal( x1 ,x2 )                   ,  LS.STOP = sprintf('isequal( x1(%g) , x2(%g) )'               ,x1,x2);                              break; end
      if x2-x0 < eps(x1)*100                 ,  LS.STOP = sprintf('x2(%g)-x0(%g) < eps(x1)*100'              ,x2,x0);                              break; end
      
      LS.ITGS = LS.ITGS + 1;

      if GO.LS_PLOT
        DrawInGLSline('line','XData',PLOT.LS_X+[x0 x1 x2],'YData',[e0 e1 e2], 'line','none','marker','s','markersize',8,'markerfacecolor',[0.5 0.5 0.5] );
        xx = linspace( x0 , x2 , 100 );
        yy = ( (xx - x1).*( e2*(xx - x0)*(x0 - x1) + e0*(xx - x2)*(x1 - x2) ) - e1*(xx - x0).*(xx - x2)*(x0 - x2) )/( (x0 - x1)*(x0 - x2)*(x1 - x2) );

        hl = DrawInGLSline('line','XData',PLOT.LS_X+xx,'YData',yy, 'line',':','marker','none','color',[0.5 0.5 0.5] );
      end

      ap = a;
      [a,A,B] = FindQ( x0 , x1 , x2 , e0 , e1 , e2 );
      Vprintf(2,'\b( %8g  %8g  %8g ).( %8g  %8g  %8g ) -> %8g   ( A: %g  B: %g )',x0,x1,x2,e0,e1,e2,a,A,B );
      
      if a < 0 %&& ap < 0
        LS.STOP = 'minimun on negative.';
        break;
      end

      if A > 0
        
        Ma = FLS( a );
        if GO.LS_PLOT
          DrawInGLSline('line','XData',PLOT.LS_X+a,'YData',Ma, 'line','none','marker','s','markersize',9,'markerfacecolor',[1 0 0] );
        end

        if Ma < M
          M = Ma;
        end

        if      a < x0
          DX = ( x2 - x0 )/10;
          
          x2 = x1; e2 = e1;
          x1 = x0; e1 = e0;
          x0 = a;  e0 = Ma;
        elseif  a < x1
          DX = x1-x0;
          
          x2 = x1; e2 = e1;
          x1 = a ; e1 = Ma;
        elseif  a < x2
          DX = x2-x1;
          
          x0 = x1; e0 = e1;
          x1 = a ; e1 = Ma;
        else
          DX = ( x2 - x0 )/10;

          x0 = x1; e0 = e1;
          x1 = x2; e1 = e2;
          x2 = a;  e2 = Ma;
        end
        if LS.ITGS > 1 && abs( ( a-ap )/DX ) < LSP.LS_QUADRATIC_TOL ,  LS.STOP = 'QUADRATIC converged.';   break; end
        
      elseif A < 0

        if e2 < E 
          Vprintf(-2,'concave ... increase ..');

          x1 = x2; e1 = e2;
          x2 = x0 + 2*(x2-x0);
          e2 = FLS( x2 );
          if e2 < M , M = e2; end
        else
          Vprintf(-2,'concave ... reduce ..');
          
          x2 = x1; e2 = e1;
          x1 = 0.5*(x0+x1);
          e1 = FLS( x1 );
          if e1 < M , M = e1; end
        end          
      
      elseif A == 0
        
        Vprintf(-2,'QUADRATIC flat!!');
        
        if B >= 0

          Vprintf(-2,'reduce ..');
          
          x2 = x1; e2 = e1;
          x1 = 0.5*(x0+x1);
          e1 = FLS( x1 );
          if e1 < M , M = e1; end
          
        else

          Vprintf(-2,'increase ..');
          
          x1 = x2; e1 = e2;
          x2 = x0 + 2*(x2-x0);
          e2 = FLS( x2 );
          if e2 < M , M = e2; end

        end
        
      end
      
    end
    
    try
      set( hl , 'linewidth', 2 ,'linestyle', '-' ,'color',[1 0 1] );
    end
    
    [M,id] = min([e0 e1 e2]);
    a = [x0 x1 x2]; a = a(id);

    Vprintf( 2,'\b');
    Vprintf(-2,'LS.STOP: %s', LS.STOP );
  end
  function [a,M] = exhaustive( LSP , xf )
    x0 = 0;
    x1 = xf;
    ee = [ FLS(x0)  FLS(x1) ];
    Vprintf( 2 ,'exhaustive search it: %4d  ( x0: %g   x1: %g ) (%e)', 0 , x0 , x1 , x1-x0 );
    for it = 1:LSP.LS_EXHAUSTIVE_ITERATIONS
      if isequal( x0 , x1 )
        Vprintf(2,'\nInterval too small');
        break;
      end
      if x1 < LSP.MIN_STEP
        Vprintf(2,'\nx1(%g) < MIN_STEP(%g)',x1,LSP.MIN_STEP);
        break;
      end
      
      xx = unique( linspace( x0 , x1 , LSP.LS_EXHAUSTIVE_N ) );
      ee = xx*0;
      for ii = 1:numel(ee)
        ee(ii) = FLS(xx(ii));
      end

      if GO.LS_PLOT
        DrawInGLSline('line','XData',PLOT.LS_X+xx,'YData',ee,'linestyle','none','marker','.','color', [ (it-1)/(max(LSP.LS_EXHAUSTIVE_ITERATIONS,2)-1)  0 1] );
      end

      [M,eeid] = min( ee );
      if eeid == 1
        x0 = xx(1); x1 = xx(2);
      elseif eeid == numel(xx)
        x0 = (xx(1)+xx(end))/2;
        x1 = xx(end) + (xx(1)+xx(end))/2;
      else
        x0 = xx(eeid-1); x1 = xx(eeid+1);
      end

      Vprintf( 2 ,'\bexhaustive search it: %4d  ( x0: %g   x1: %g ) (%e)', it , x0 , x1 , x1-x0 );
    end
    
    Vprintf(2,'\b');
    [M,eeid] = min( ee );
    a = xx(eeid);
    
  end
  function [a,M] = smartLS( LSP )
    function a = fit_quadratic( aa , ee )

%       while numel(aa) > 3
%         [mm,mmid] = max(ee);
%         ee(mmid)=[];
%         aa(mmid)=[];
%       end
      
      polyW = diag( diff( dualVector(aa) ) );
      polyM = [ ones(numel(aa),1) aa(:) aa(:).^2 ];

      setWarning('off','MATLAB:rankDeficientMatrix');
      setWarning('off','MATLAB:nearlySingularMatrix');
      setWarning('off','MATLAB:singularMatrix');
      ABC = linsolve( polyM.' * polyW * polyM ,  polyM.' * polyW * ee(:) );
      restoreWarning(  'MATLAB:rankDeficientMatrix');
      restoreWarning(  'MATLAB:nearlySingularMatrix');
      restoreWarning(  'MATLAB:singularMatrix');

      ABC(4) = 0;
      
      if ABC(3) <= 0
        Vprintf(-2, 'concave fit  %g,%g,%g', ABC(1),ABC(2),ABC(3) );
        a = Inf;
        return;
      end

      a = -0.5 * ABC(2)/ABC(3);
    end
    
    
    function a = fit_cubic( aa , ee )
      l1 = min(aa);
      l2 = max(aa);
      
%       while numel(aa) > 4
%         [mm,mmid] = max(ee);
%         ee(mmid)=[];
%         aa(mmid)=[];
%       end
      
      polyW = diag( diff( dualVector(aa) ) );
      polyM = [ ones(numel(aa),1) aa(:) aa(:).^2 aa(:).^3 ];

      setWarning('off','MATLAB:rankDeficientMatrix');
      setWarning('off','MATLAB:nearlySingularMatrix');
      setWarning('off','MATLAB:singularMatrix');
      ABC = linsolve( polyM.' * polyW * polyM ,  polyM.' * polyW * ee(:) );
      restoreWarning(  'MATLAB:rankDeficientMatrix');
      restoreWarning(  'MATLAB:nearlySingularMatrix');
      restoreWarning(  'MATLAB:singularMatrix');

      
      discri = ABC(3)^2 - 3*ABC(2)*ABC(4);
      if discri < 0 || abs( ABC(4) ) < 1e-10
        a = fit_quadratic(aa,ee);
        return;
      end
      
      discri = sqrt( discri );
      r1 = ( -ABC(3) - discri ) / ( 3 * ABC(4) );
      r2 = ( -ABC(3) + discri ) / ( 3 * ABC(4) );
      
      if isinf( r1 ) || isnan(r1) || isinf( r2 ) || isnan(r2)
        a = fit_quadratic(aa,ee);
        return;
      end
      
      if l1 < r1  && l2 > r1  &&  l1 < r2  && l2 > r2
        a = fit_quadratic(aa,ee);
      elseif l1 < r1  && l2 > r1
        a = r1;
      elseif l1 < r2  && l2 > r2
        a = r2;
      elseif   ABC(2) + 3*ABC(4)*r1  > 0
        a = r1;
      elseif   ABC(2) + 3*ABC(4)*r2  > 0
        a = r2;
      else
        a = fit_quadratic(aa,ee);
      end
      
    end
    
    
    
    
    
    [M,mmid] = min( FLS_e );
    a = FLS_a(mmid);

    if numel( FLS_e ) < 3
      Vprintf(22,'SMART_LS no se puede hacer porque hay muy pocos datos');
    end

    Vprintf( 2 , '%4d - %2d ( %.15g )' , 0 , a ,  M );
    for it = 1:LSP.SMART_LINE_SEARCH_ITS

      aa = FLS_a;
      ee = FLS_e;
      [aa,ord] = sort( aa );
       ee      = ee( ord );
       
      %seleccionar los puntos  para hacer el fiteo.
      [mm,mmid] = min(ee);
      mmid1 = mmid;
      while mmid1 > 1
        if ee(mmid1-1) > ee(mmid1), mmid1 = mmid1-1;
        else,  break; end
      end
      mmid2 = mmid;
      while mmid2 < numel(aa)
        if ee(mmid2+1) > ee(mmid2), mmid2 = mmid2+1;
        else,  break; end
      end
      aa = aa( mmid1:mmid2 );
      ee = ee( mmid1:mmid2 ); 
      
      if numel(aa) < 3
        Vprintf(-2,'Too few points');
        break;
      end
      
      ABC = [];
      if numel(aa) > 3
        a = fit_cubic( aa , ee );
      else
        a = fit_quadratic( aa , ee );
      end
        
      
      if 0
        figure;
        line('XData',FLS_a,'YData',FLS_e,'linestyle','none','marker','o');
        line('XData', aa  ,'YData', ee  ,'linestyle','none','marker','o','markerfacecolor',[1 0 0]);
        
        xxx = centerscale( [ min(FLS_a) , max(FLS_a) ] , 2 );
        xxx = linspace( xxx(1) , xxx(2) , 1000 );
        yyy = ABC(1) + ABC(2)*xxx + ABC(3)*xxx.^2 + ABC(4)*xxx.^3;
        line('XData',xxx,'YData',yyy,'linestyle',':','marker','none');
        
        aaa = FLS_a;
        eee = FLS_e;
        
        line('XData', [a a] ,'YData', [ ABC(1) + ABC(2)*a + ABC(3)*a.^2 + ABC(4)*a.^3 , FLS(a) ] ,'linestyle','none','marker','+','linewidth',2,'markersize',12);
        
        xxx = centerscale( [ min(FLS_a) , max(FLS_a) ] , 2 );
        xxx = unique( [ FLS_a  ,  linspace( xxx(1) , xxx(2) , 100 ) ] );
        yyy = xxx*0;
        for iii = 1:numel(xxx)
          yyy(iii) = FLS(xxx(iii));
        end
        
        line('XData',xxx,'YData',yyy,'linestyle','-','marker','none','color',[1 0 0]);

        FLS_a = aaa;
        FLS_e = eee;
      end
      
      if a < 0
        Vprintf(-2, 'negative a  %g', a );
        break;
      end
      if isnan( a ) || isinf( a ) || a > 1e10
        Vprintf(-2, 'a  invalid  %g', a );
        break;
      end
      if min( abs( FLS_a - a ) )/( max(FLS_a) - min(FLS_a) )  < 1e-10
        Vprintf(-2, 'smart converged');
        break;
      end

      M = FLS( a );

      if GO.LS_PLOT
        DrawInGLSline('line', 'XData' , PLOT.LS_X + a , 'YData' , M , 'line','none' ,'marker','+','markersize',25 , 'color',[0 0.8 0] ,'linewidth' ,2);
      end
      
%       [M,mmid] = min( FLS_e );
%       a = FLS_a(mmid);

      Vprintf( 2 , '\b%4d - %2d ( %.15g )' , it , a ,  M );
    end
    Vprintf(2,'\b');

    [M,mmid] = min( FLS_e );
    a = FLS_a(mmid);
  end































%   function a = derivativesearch( Fu , xi , xf )
%     h = xf/20; h =[-h h];
%     a = 0;
%     if xi == 0
%       f = E; dxi = SLOPE0;
%     else
%       [f,g] = Fu(xi); dxi = INNER( g(:),Dn );
%     end
%     if LS_PLOT > 1
%       line('Parent',hGroupLS,'Color',[0.4 0.4 0.4], ...
%         'linewidth',3,...
%         'XData', xi + LS_PLOT_X + h , ...
%         'YData', f + dxi*h ); 
%     end
%     
%     [f,g] = Fu(xf); dxf = INNER( g(:),Dn );
%     if LS_PLOT > 1
%       line('Parent',hGroupLS,'Color',[0.4 0.4 0.4], ...
%         'linewidth',3,...
%         'XData', xf + LS_PLOT_X + h , ...
%         'YData', f + dxf*h );
%     end
% 
%     while dxf < 0
%       xf      = P.BRACKETING_RHO*xf;
%       [f,g] = Fu(xf); dxf = INNER( g(:),Dn );
%       if LS_PLOT > 1
%         line('Parent',hGroupLS,'Color',[0.4 0.4 0.4], ...
%           'linewidth',3,...
%           'XData', xf + LS_PLOT_X + h , ...
%           'YData', f + dxf*h );
%       end
%     end
%     
%     while dxi < 0 && dxf > 0  &&  xf-xi > 10*eps(xi)
%       a = ( xi*dxf - xf*dxi )/( dxf-dxi );
%       if min( abs(a-xi) , abs(a-xf) ) < ( xf-xi )/100
%         return;
%       end
%       [f,g] = Fu(a); d = INNER( g(:),Dn );
%       if LS_PLOT > 1
%         line('Parent',hGroupLS,'Color',[0.4 0.4 0.4], ...
%           'linewidth',3,...
%           'XData', a + LS_PLOT_X + h , ...
%           'YData', f + d*h );
%       end
%       if d < 0
%         xi = a; dxi = d;
%       elseif d > 0
%         xf = a; dxf = d;
%       else
%         break;
%       end
%     end
%   end
%   function [a,M] = fastsearch( Fu , a , factors )
%     M = Fu( a );
%     if LS_PLOT > 1, line('Parent',hGroupLS,'Color',[0.4 0.4 0.4],'marker','h','linestyle','none','XData',a+LS_PLOT_X,'YData',M,'MarkerFaceColor',[1 0 0]); end
%     if ~iscell( factors ), factors = {factors}; end
%     for f = factors
%       xx = a*setdiff( f{1} , 1);
%       yy = vectorialF( Fu , xx );
% 
%       if LS_PLOT > 1, line('Parent',hGroupLS,'Color',[0.4 0.4 0.4],'marker','d','linestyle','none','XData',xx+LS_PLOT_X,'YData',yy,'MarkerFaceColor',rand(1,3)); end
%       
%       
%       [MM,idx] = min( yy );
%       if MM < M
%         a = xx(idx);
%         M = MM;
%       end
%     end
%   end
%   function stop = FMINSEARCH_OutputFcn( x , optimValues , state )
%     stop = 0;
%     switch state
%       case 'iter'
%         if E == optimValues.fval;
%           Vprintf2( '\b (%s) \n',optimValues.procedure )
%           return;
%         end
%         
%         E = optimValues.fval;
%         updateX( x );
%         if ITT >= P.MAX_ITERATIONS_PER_METHOD || IT >= MAX_ITERATIONS
%           stop = 1;
%         end
%       case 'done'
%         E = optimValues.fval;
%         updateX( x );
%         Vprintf('converged... \n');
%     end
%   end


  
%   Vprintf('                E_Init: %30.20g\n',E );
%   if SAVE_HISTORY, HISTORY.E0 = E; end
% 
%   alpha = 0;
%   
%   P0 = P;
%   x_old = x + Inf;
%   IT = 0; E0 = E;
%   while IT < MAX_ITERATIONS
%     x_old = x;
%     for mm = 1:numel(meths)
%       if IT >= MAX_ITERATIONS, break; end
%       PP = [];
%       m = meths{mm};
%       if ~ischar( m ), continue; end
%       if numel(meths) > mm  &&  isstruct( meths{mm+1} )
%         PP = meths{mm+1};
%       elseif numel(meths) > mm  &&  isnumeric( meths{mm+1} )
%         PP.MAX_ITERATIONS_PER_METHOD = meths{mm+1};
%       else
%         PP = P0;
%       end
%       
%       if ~isempty( PP )
%         for fn = fieldnames( PP )
%           P.(fn{1}) = PP.(fn{1});
%         end
%       end
%       
%       
%       
%       
%       G      = NaN;
%       H      = NaN;
%       beta   = NaN;
%       SLOPE0 = NaN;
%       Dn     = NaN;
%       inv_H  = NaN;
%       G0     = NaN;
%       ITT    = 0;
%       
%       switch m
%         case 'fminsearch',  START('FMINSEARCH');
%           FMINSEARCH_OPTIONS              = optimset;
%           FMINSEARCH_OPTIONS.MaxFunEvals  = Inf;
%           FMINSEARCH_OPTIONS.MaxIter      = Inf;
%           FMINSEARCH_OPTIONS.OutputFcn    = @FMINSEARCH_OutputFcn;
%           FMINSEARCH_OPTIONS.Display      = 'off';
%           
%           FMINSEARCH_OPTIONS.TolFun       = 0;
%           FMINSEARCH_OPTIONS.TolX         = 0;
%           
%           x = fminsearch( @(x) F(R(x)) , x , FMINSEARCH_OPTIONS );
%           x = x(:);
%         
%       end
%     end
%     

end


      %{
      [mm,mmid] = min( FLS_e );

      ee = FLS_e >= mm + 100*eps(mm); ee( mmid ) = 1;
      if sum( ee ) > 2*NumberOfPointsForFit +1
        aa = FLS_a( ee );
        ee = FLS_e( ee );
      else
        aa = FLS_a;
        ee = FLS_e;
      end

      [aa,ord] = sort( aa );
      ee = ee(ord);

      [mm,mmid] = min( ee );
      mmid1 = max( 1 , mmid - NumberOfPointsForFit );
      mmid2 = min( numel(ee) , max( mmid1 + 2*NumberOfPointsForFit , mmid + NumberOfPointsForFit ) );
      mmid1 = max( 1 , mmid2 - 2*NumberOfPointsForFit );

      aa = vec( aa( mmid1:mmid2 ) );
      ee = vec( ee( mmid1:mmid2 ) );

      setWarning('off','MATLAB:rankDeficientMatrix');
      ABC = linsolve( [ aa.^2 ,  aa  , ones(numel(aa),1) ] , ee );
%         ABC = linsolve( [ FLS_a(:).^2 , FLS_a(:) , ones(numel(FLS_a(:)),1) ] , FLS_e(:) );
      restoreWarning(  'MATLAB:rankDeficientMatrix');

      new_a = - 0.5 * ABC(2) / ABC(1);
      %}

      %{
      ee = 0; dd = E;
      while sum(ee) < 3
        ee = FLS_e <= dd;
        dd = min( FLS_e( ~ee ) );
      end

      aa = FLS_a(ee);
      ee = FLS_e(ee);
      [aa,ord] = sort( aa );
       ee      = ee( ord );

      polyM = sum( [    2*( aa(2:end)    - aa(1:end-1)    ) ;...
                          ( aa(2:end).^2 - aa(1:end-1).^2 ) ;...
                      2/3*( aa(2:end).^3 - aa(1:end-1).^3 ) ;...
                      1/2*( aa(2:end).^4 - aa(1:end-1).^4 ) ;...
                      2/5*( aa(2:end).^5 - aa(1:end-1).^5 ) ;...
                      1/3*( aa(2:end).^6 - aa(1:end-1).^6 ) ;...
                      2/7*( aa(2:end).^7 - aa(1:end-1).^7 ) ] , 2 );

      polyM = polyM([ 1 2 3 4 ; 2 3 4 5 ; 3 4 5 6 ; 4 5 6 7 ]);

      polyR = sum( [      diff(aa).*(                          ee(1:end-1) + ee(2:end) )                                ; ...
                      1/3*diff(aa).*( aa(1:end-1).*        ( 2*ee(1:end-1) + ee(2:end)) + aa(2:end).*(ee(1:end-1)+2*ee(2:end)) )  ; ...
                      1/6*diff(aa).*( aa(1:end-1).*aa(2:end).*(2*ee(1:end-1)+2*ee(2:end)) + aa(1:end-1).^2.*(3*ee(1:end-1)+ee(2:end)) + aa(2:end).^2.*(ee(1:end-1)+3*ee(2:end)) )  ; ...
                     1/10*diff(aa).*( aa(1:end-1).^3.*(4*ee(1:end-1)+ee(2:end)) + aa(1:end-1).^2.*aa(2:end).*(3*ee(1:end-1)+2*ee(2:end)) + aa(1:end-1).*aa(2:end).^2.*(2*ee(1:end-1)+3*ee(2:end)) +  aa(2:end).^3.*(ee(1:end-1)+4*ee(2:end)) ) ] , 2 );

      ABC = polyM\polyR;        

      if abs( ABC(4) ) > 1e-8

        discri = ABC(3)^2 - 3*ABC(2)*ABC(4);
        if discri > 0
          discri = sqrt( discri );
          if      2*ABC(3) + 6*ABC(4)*( -ABC(3) + discri )/(3*ABC(4)) > 0
            new_a = ( -ABC(3) + discri )/(3*ABC(4));
          elseif  2*ABC(3) + 6*ABC(4)*( -ABC(3) - discri )/(3*ABC(4)) > 0
            new_a = ( -ABC(3) - discri )/(3*ABC(4));
          else
            new_a = Inf;
          end
        else
            new_a = Inf;
        end

      else

        new_a = -0.5 * ABC(2)/ABC(3);

      end
      %}
