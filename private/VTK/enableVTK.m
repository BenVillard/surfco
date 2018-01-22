function enableVTK

  try
    evalc( 'vtkPolyDataReader()' );
    return;
  catch
    switch computer
      case 'PCWIN64'
        
        while 1
          d = fullfile( fileparts(which('enableVTK')) , 'w64' );
          if isdir( d ), setenv( 'path' , [ getenv('path') , ';' , d ] ); break; end
          
          d = fullfile( fileparts(which('enableVTK')), 'w64_mt' );
          if isdir( d ), setenv( 'path' , [ getenv('path') , ';' , d ] ); break; end
          
          break;
        end
      
      case 'PCWIN32'
        setenv( 'path' , [ getenv('path') , ';' , fullfile( fileparts(which('enableVTK')) , 'vtk' , 'w32' ) ] );

       case 'MACI64'
        setenv( 'path' , [ getenv('path') , ';' , fullfile( fileparts(which('enableVTK')) , 'vtk' , 'maci64' ) ] );
        
      otherwise
        error('cannot enable VTK within a matlab runtime.')
    end
  end

  try
    evalc( 'vtkPolyDataReader()' );
  catch
    error('VTK couldn''t be enabled.');
  end
  
end
