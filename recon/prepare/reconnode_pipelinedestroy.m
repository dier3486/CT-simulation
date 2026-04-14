function [dataflow, prmflow, status] = reconnode_pipelinedestroy(dataflow, prmflow, status)
% pipeline destroy node
% [dataflow, prmflow, status] = reconnode_pipelinedestroy(dataflow, prmflow, status);
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

flag_destoryed = false;
% to record the error message from the prepare fucntions 
dstrerrormsg = '';
dstrerrorcode = 0;
% loop the nodes
for ii = 1:length(pipenodes)
    % spell the nodes prepare function
    nodename_slip = regexp(pipenodes{ii}, '_', 'split');
    destroyname = ['reconnode_' lower(nodename_slip{1}) 'destroy'];
    % Do not name a pipe line node in 'pipeline' ;)  
    status.nodename = pipenodes{ii};
    % to call the node's prepare function
    if any(exist(destroyname) == [2 5 6]) % can run
        if echo_onoff
            if flag_destoryed
                fprintf(', ');
            else
                fprintf(' [');
                flag_destoryed = true;
            end
            fprintf(pipenodes{ii});
        end
        destroyfun = str2func(destroyname);
        [dataflow, prmflow, status] = destroyfun(dataflow, prmflow, status);
        if status.jobdone == 0
            return
        elseif status.errorcode ~=0
            % record them
            dstrerrorcode = status.errorcode;
            dstrerrormsg = [dstrerrormsg status.errormsg ' '];
            status.errorcode = 0;
            status.errormsg = '';
        end
        % set the flag destroyed in status.pipeline to true
        status.pipeline.(status.nodename).destroyed = true;
    end
end
if echo_onoff 
    if flag_destoryed
        fprintf(']');
    end
    fprintf('\n');
end

% unload librarys
for ilib = status.libs
    if libisloaded(ilib{1})
        unloadlibrary(ilib{1});
    end
end

if libisloaded('CUpipeline')
    unloadlibrary('CUpipeline');
end

% recover the nodename
status.nodename = orignodename;

% to return the error message
if dstrerrorcode ~= 0
    status.errorcode = dstrerrorcode;
    status.errormsg = dstrerrormsg;
end

% done
status.jobdone = true;

end
