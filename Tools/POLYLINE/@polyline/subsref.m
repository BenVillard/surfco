function [ P ] = subsref( P , s )

  for ss = 1:numel( s )
    if ~isa( P , 'polyline' )
      P = subsref( P , s(ss:end) );
      return;
    end
    
    switch s(ss).type
      case { '{}' }
        error('no curly subsref allowed for POLYLINE');
      
      case { '()' }
        if numel( s(ss).subs ) > 1, error('only global indexing is allowed for POLYLINE'); end
        try
          P.C = subsref( P.C , s(ss) );
        catch LE
          throw( LE );
        end

      case { '.' }
        if ismethod( P , s(ss).subs )
          if numel(s) == ss+1 && strcmp( s(ss+1).type , '()' )
            try
              P = feval( s(ss).subs , P , s(ss+1).subs{:} );
              return;
            end
          end
          P = feval( s(ss).subs , P );
        else
          switch s(ss).subs
            case {'nsd'},                   P = nsd(P);
            case {'np'},                    P = np(P);
            case {'n','nn'},                P = NumberOfNodes( P );
            case {'data','coords'},         P = coordinates( P );
            otherwise, error('invalid calling .%s for a polyline object',s(ss).subs);
          end
        end
    end
  end

end
