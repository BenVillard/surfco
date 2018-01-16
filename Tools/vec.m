function v = vec( v , dim )
  
  if nargin < 2 , dim = 1; end

%   if iscell( v )
%     for i= numel(v):-1:1
%       if isempty( v{i} )
%         v(i)=[];
%       end
%     end
%     v= cell2mat( v(:) );
%   else
    v= v(:);
%   end
  
  if dim ~= 1
    
    v = reshape( v , [ ones(1,dim-1) , numel(v) ] );
    
  end
  
  
end
