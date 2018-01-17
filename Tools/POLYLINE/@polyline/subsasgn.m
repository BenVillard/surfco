function P = subsasgn( P , s , in )

  if numel( s ) ~= 1, error('too complex subasgn'); end
  
  nd = nsd( P );
  
  switch s(1).type
    case '()'
      if isempty( in )
        P.C( s(1).subs{:} ) = in;
        return;
      end
      
      if ~isa( in , 'polyline' )
        in = polyline( in );
      end
      if ~isnan( nd ) && ~isequal( nd , nsd(in) )
        error('inconsistent NumberSpatialDimensions');
      end
        
      P.C( s(1).subs{:} ) = in.C;
        
      
    case '.'
      switch s(1).subs
        case 'nsd'
          if numel(in) ~= 1 || in < 1 || rem( in ,1)
            error('invalied new NumberSpatialDimensions');
          end
          if      in == nd
          elseif  in <  nd
            for p = 1:np( P )
              P.C{p} = P.C{p}(:,1:in);
            end
          elseif  in >  nd
            for p = 1:np( P )
              P.C{p}(:,end+1:in) = 0;
            end
          end
          
        otherwise, error('no allowed assign');
      end
      
  end


end
