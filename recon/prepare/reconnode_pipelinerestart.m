function [dataflow, prmflow, status] = reconnode_pipelinerestart(dataflow, prmflow, status)
% pipeline prepare node
% [dataflow, prmflow, status] = reconnode_pipelineprepare(dataflow, prmflow, status);
% Plz run this node after loading calibration tables

% Copyright Dier Zhang
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% echo onoff
echo_onoff = status.debug.echo_onoff;

% backup status.nodename
orignodename = status.nodename;

% prepare of the pipeline nodes
pipenodes = fieldnames(status.pipeline);

flag_restart = false;
% to record the error message from the prepare fucntions 
preperrormsg = '';
preperrorcode = 0;
% loop the nodes
for ii = 1:length(pipenodes)
    % spell the nodes prepare function
    nodename_slip = regexp(pipenodes{ii}, '_', 'split');
    preparename = ['reconnode_' lower(nodename_slip{1}) 'restart'];
    % Do not name a pipe line node in 'pipeline' ;)  
    status.nodename = pipenodes{ii};
    % to call the node's prepare function
    if any(exist(preparename) == [2 5 6]) % can run
        if echo_onoff
            if flag_restart
                fprintf(', ');
            else
                fprintf(' [');
                flag_restart = true;
            end
            fprintf(pipenodes{ii});
        end
        preparefun = str2func(preparename);
        [dataflow, prmflow, status] = preparefun(dataflow, prmflow, status);
%         if status.jobdone == 0
%             return
%         elseif status.errorcode ~=0
%             % record them
%             preperrorcode = status.errorcode;
%             preperrormsg = [preperrormsg status.errormsg ' '];
%             status.errorcode = 0;
%             status.errormsg = '';
%         end
        % set the flag prepared in status.pipeline to true
        status.pipeline.(status.nodename).prepared = true;
    end
end
if echo_onoff && flag_restart
    fprintf(']   done\n');
end

% call back the nodename
status.nodename = orignodename;

% to return the error message
if preperrorcode ~= 0
    status.errorcode = preperrorcode;
    status.errormsg = preperrormsg;
end

% done
status.torestart = false;
status.jobdone = true;

end