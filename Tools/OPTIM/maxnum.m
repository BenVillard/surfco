function m = maxnum( type )

  if nargin < 1, type = 'double'; end
  if ~ischar( type ), type = class( type ); end
  
  switch lower( type )
    case 'logical',       m = double( 1 );
      
    case 'uint8',         m = double( intmax('uint8') );
    case 'int8',          m = double( intmax('int8') );

    case 'uint16',        m = double( intmax('uint16') );
    case 'int16',         m = double( intmax('int16') );

    case 'uint32',        m = double( intmax('uint32') );
    case 'int32',         m = double( intmax('int32') );

    case 'uint64',        m = double( intmax('uint64') );
    case 'int64',         m = double( intmax('int64') );

    case 'single',        m = 3.4028235677973362e038;    % 3.402823567797336427480734639795e38
    case 'double',        m = 1.7976931348623157e308;    % 1.7976931348623158079372897140530341507993413271003782693e308
      
  end

 
end