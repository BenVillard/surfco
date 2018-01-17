function [P,nt] = part( P , mode , ft )

  if ~isSingle( P ), error('only single polylines can be resamplig'); end

  if ~ischar( mode )
    ft   = mode;
    mode = 'index';
  end

  if numel( ft ) ~= 2, error('a starting coordinate and a endding coordinate was expected in from_to'); end
  
  FLIP_at_END = false;
  if ~issorted( ft )
    FLIP_at_END = true;
    ft = ft([2 1]);
  end
  
  switch lower( mode )
    case {'index','i'}
      if isinf( ft(1) ) && ft(1) < 0
        ft(1) = 1;
      end
      if isinf( ft(2) ) && ft(1) > 0
        ft(2) = nn(P,1);
      end
      nt = unique( [ ft(1) , ceil(ft(1)):floor(ft(2)) , ft(2) ] );
      P = resample( P , 'index' , nt );
      
    otherwise
      error('not implemented yet');
  end

  
  if FLIP_at_END
    try,    P.C{1} = flip(    P.C{1} , 1 );
    catch,  P.C{1} = flipdim( P.C{1} , 1 );
    end
  end
  
end