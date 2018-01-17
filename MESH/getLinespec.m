function v = getLinespec( v )

  if ~iscell(v), v = {v}; end

  V_lstyle = ''; V_color = ''; V_marker = '';

  if numel(v) && ischar( v{1} )
    [V_lstyle,V_color,V_marker,msg] = colstyle( v{1} );
    if isempty( msg ), v(1) = []; end
  end
  if ~isempty(V_marker), v = [ 'Marker',    V_marker , v ]; end
  if ~isempty(V_lstyle), v = [ 'LineStyle', V_lstyle , v ]; end
  if ~isempty(V_color )
%     switch lower(V_color)
%       case 'r', V_color = [1 0 0];
%       case 'g', V_color = [0 1 0];
%     end
    
    v = [ 'Color',     V_color  , v ]; 
  end


end

