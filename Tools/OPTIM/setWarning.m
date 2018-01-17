function setWarning( new_state , wid )


  old_state = warning( 'query' , wid );
  warning( new_state , wid );
  
  ALL_W = getappdata(0,'WARNING_STATES');
  if isempty( ALL_W ), ALL_W = {}; end
  ALL_W = [ ALL_W ; { wid   old_state.state } ];

  setappdata(0,'WARNING_STATES',ALL_W );
  
end

