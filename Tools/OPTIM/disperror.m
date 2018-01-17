function E_ = disperror( E , links )

  isC = false;
  
  if nargin == 0
    E = lasterror;
    links = true;
  end
  
  if nargin == 1
    if isnumeric( E ) || islogical( E )
      links = E;
      E = lasterror;
    else
      links = true;
    end
  end
  
  if isempty( E.message )
    fprintf('No error message stored\n');
    return;
  end

  returnE = false;
  if nargout > 0, 
    returnE = true;
    E_ = ''; 
  end
    
  nprintf( 'error_identifier : %s\n' , E.identifier );
  nprintf( 'error_message    : %s\n' , E.message    );

  if isa( E , 'MException' )
    
    setWarning('off','MATLAB:structOnObject');
    EE = struct( E );
    restoreWarning(  'MATLAB:structOnObject');
    
    nprintf( 'error_cause      : %s\n' , uneval( EE.cause ) );
    nprintf( 'error_type       : %s\n' , uneval( EE.type  ) );
  end
    
  nprintf('\n');
  for i=1:numel( E.stack )
    
    s = E.stack(i);
    ff = which( s.file );
    [ p , N , ex ] = fileparts( ff );
    %N = [ N ex ];
    n = s.name;
    
    if links
      href = sprintf('matlab:opentoline(''%s'',%d)',ff,s.line);

      if strcmp( n , N )
        nprintf( '   <a href="%s">%s ( line: %d )</a>\n' , href , ff , s.line );
      else
        nprintf( '   <a href="%s">%s -> %s ( line: %d )</a>\n' , href , ff , n , s.line );
      end
    
    else
      
      if strcmp( n , N )
        nprintf( '   %s ( line: %d )\n' ,  ff , s.line );
      else
        nprintf( '   %s -> %s ( line: %d )\n' , ff , n , s.line );
      end
      
    end
    
    nprintf('\n');
    
  end
  

  function nprintf( varargin )
    if isC, fprintf( 2 , varargin{:} ); end
    fprintf( varargin{:} );
    if returnE
      E_ = [ E_ , sprintf( varargin{:} ) ];
    end
  end
end
