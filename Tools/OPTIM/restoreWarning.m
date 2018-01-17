function restoreWarning( wid )

  ALL_W = getappdata(0,'WARNING_STATES');
  if isempty( ALL_W ), ALL_W = {}; end

  incell_id = find( strcmp( ALL_W(:,1) , wid ) , 1 , 'last' );
  if isempty(incell_id)
    error('no hay estado guardado para este warning_id');
  end
  
  old_state = ALL_W{incell_id, 2 };
  ALL_W(incell_id,:) = [];
  setappdata(0,'WARNING_STATES',ALL_W );
  
  warning( old_state , wid );

  
end
