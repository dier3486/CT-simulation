function r = nodestatecheck(nodestatus, field, nodesindex, logic, oprate)
% a pipeline tool to check the logical state of the nodes
%   r = nodestatecheck(nodestatus, field, nodesindex, logic);
% or r = nodestatecheck(nodestatus, field, nodesindex, logic, opt);
% nodestatus:       the stract to be checked, e.g. status.pipeline, status.pipepool;
% field:            the field name to be checked, e.g. 'sleeping';
% nodesindex:       the index of the nodes to be checked, e.g. 0 : status.pipeline.(nodename).index - 1;
% logic:            'all' or 'any'
% oprate:           oprate function, a string, r.g. @(x)x>1.


nodes = fieldnames(nodestatus);
nodesindex = nodesindex(:)';

if nargin  > 4
    % opt
    funopt = str2func(oprate);
else
    funopt = @(x)x;
end

% get the checks of the nodes
% check = true(0);
% for ii = 1 : length(nodesindex)
%     if isfield(nodestatus.(nodes{nodesindex(ii)+1}), field)
%         check = [check funopt(nodestatus.(nodes{nodesindex(ii)+1}).(field))];
%     end
%     % I know the status.pipeline.(nodes{nodesindex(ii)+1}).index == nodesindex(ii), bug if not satisfied.
% end
check = true(0);
for ii = 1:length(nodes)
    index_ii = nodestatus.(nodes{ii}).index;
    if any(index_ii==nodesindex)
        check = [check funopt(nodestatus.(nodes{ii}).(field))];
    end
end

% logical
switch lower(logic)
    case 'all'
        % all
        r = all(check);
        % note that all([]) is true, therefore if you check an unexisting field it wiil be true.
    case 'any'
        % any
        r = any(check);
    % we need not all(~check), any(~check), do we?
    otherwise
        funlogic = str2func(logic);
        r = funlogic(check);
end

end