function [ varargout ] = parseargs( in , varargin )

if 0
  %     1  2      3  4  5  6      7 8    9       10   11    12 13
  in = {0,'first',1,'s',2,'third',3,3.1,'second',2.2,'yes','N',Inf};
  [var,pos,val] = parseargs(in,'First','$DEFS$',-1);
  [var,pos,val] = parseargs(in,'Second','sec','$DEFS$',-2);
  [var,pos,val1,val2] = parseargs(in,'third','$DEFS$',-3,-3.1);
  [var,pos] = parseargs(in,'Yes','$FORCE$',{1,2})
  [var,pos] = parseargs(in,'No','$FORCE$',{1,2})
  [var,pos] = parseargs(in,'No')
  [var,pos] = parseargs(in,'maybe')
  [var,pos] = parseargs(in,'maybe','$FORCE$',{1,2})
end
  
  try
    [ varargout{1:nargout} ] = parseargs_mex( in , varargin{:} );
    return;
  catch LE
    if ~strcmp( LE.identifier , 'MATLAB:UndefinedFunction' )
      throw(LE)
    end
  end

%   if nargout < 1, error('at least one output is spected'); end
  if ~iscell( in ), error( 'IN should be a cell'); end
  
  Nouts           = max( nargout - 2 , 0 );
  POSITION_output = [];
  POSITION_set    = false;
  DEFAULTS_output = cell(1,Nouts);

  OPTSid = find( strcmp( varargin , '$DEFS$' ) | strcmp( varargin , '$FORCE$' ) ,1);
  if isempty( OPTSid ), OPTSid = numel( varargin ) + 1; end
  OPTS = varargin( OPTSid:end );

  o = 1;
  while o <= numel( OPTS )
    if strcmp( OPTS{o} , '$FORCE$' )
      POSITION_set = true;
      %%%% TRY CATCH ADDED BY BEN 05/12/17
%       try,
          if numel( OPTS ) < o + 1 || ( numel( OPTS{o+1} ) == 1 && ischar( OPTS{o+1} ) && strcmp( OPTS{o+1} , '$DEFS$' ) )
              POSITION_output = true;
          else
              POSITION_output = OPTS{o+1};
          end
%       catch
%           if numel( OPTS ) < o + 1 || strcmp( OPTS{o} , '$DEFS$' )
%               POSITION_output = true;
%           else
%               POSITION_output = OPTS{o+1};
%           end
%       end
          
      o = o + 1 + 1;
      continue;
    end
    if strcmp( OPTS{o} , '$DEFS$' )
      if numel( OPTS ) ~= o + Nouts
        error('specification of %d defaults is expected',Nouts);
      end
      DEFAULTS_output = OPTS(o+(1:Nouts));
      o = o + Nouts + 1;
      continue;
    end
    error('it is not expected to be here!!');
  end
  if POSITION_set
    if ~iscell( POSITION_output )
      POSITION_output = { POSITION_output };
    end
    if numel( POSITION_output ) == 1
      POSITION_output{2} = false;
    end
    if numel( POSITION_output ) > 2
      error('invalid $FORCE$ specification');
    end
  end
  
  
  KEYS = cell( 1 , 2*(OPTSid-1) );
  lk = 0;
  for k = 1:( OPTSid-1 )
    K = varargin{k};
    if ~isstring(K)
      error('only strings are allowed as keys');
    end
    lk = lk + 1; KEYS{ lk } = lower( K );
    
    K( K >= 'a' & K <= 'z' ) = [];
    if isempty( K ), continue; end
    lk = lk + 1; KEYS{ lk } = lower( K );
  end
  KEYS( (lk+1):end ) = [];
  

  toDELETE = [];
  POSITION = 0;
  i = 1;
  while i <= numel( in )
    if ~isstring( in{i} ), i = i+1; continue; end
    if ~any( strcmp( lower( in{i} ) , KEYS ) ), i = i+1; continue; end
  
    if numel( in ) < i+Nouts
      error('after key ''%s'', %d inputs are expected.',KEYS{1},Nouts);
    end
    
    toDELETE = i:(i+Nouts);
    POSITION = i;
    i = i + Nouts + 1;
  end
  
  if POSITION == 0 && nargout < 2

    if POSITION_set
      varargout{1} = POSITION_output{2};
    else
      varargout{1} = 0;
    end
  
  elseif POSITION ~= 0 && nargout < 2
    
    if POSITION_set
      varargout{1} = POSITION_output{1};
    else
      varargout{1} = POSITION;
    end
    
  elseif POSITION == 0

    for i = 1:Nouts
      varargout{2+i} = DEFAULTS_output{ i };
    end
    if POSITION_set
      varargout{2} = POSITION_output{2};
    else
      varargout{2} = 0;
    end
    varargout{1} = in;

  else
    
    for i = 1:Nouts
      varargout{2+i} = in{ POSITION + i };
    end
    if POSITION_set
      varargout{2} = POSITION_output{1};
    else
      varargout{2} = POSITION;
    end
    in( toDELETE ) = [];
    varargout{1} = in;

  end
  
end
function s = isstring( x )
  s = ischar( x ) && ~isempty( x ) && ndims( x ) <= 2 && size( x , 1 ) == 1;
end
