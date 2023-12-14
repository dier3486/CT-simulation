function [dataflow, prmflow, status] = reconnode_pipelineprepare(dataflow, prmflow, status)
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

% backup status.nodename
orignodename = status.nodename;

% % prepare of read rawdata
% prmflow = readrawdataprepare(prmflow, status);
% status.pipeline.loadrawdata.prepared = true;

% prepare of the pipeline nodes
pipenodes = fieldnames(status.pipeline);
flag_prepared = false;

for ii = 1:length(pipenodes)
    % nodes prepare function
    nodename_slip = regexp(pipenodes{ii}, '_', 'split');
    preparename = ['reconnode_' lower(nodename_slip{1}) 'prepare'];
    % Do not name a pipe line node in 'pipeline' ;)  
    status.nodename = pipenodes{ii};
    % to call node prepare function
    if any(exist(preparename) == [2 5 6]) % can run
        if status.echo_onoff
            if flag_prepared
                fprintf(', ');
            else
                fprintf(' [');
                flag_prepared = true;
            end
            fprintf(nodename_slip{1});
        end
        preparefun = str2func(preparename);
        [dataflow, prmflow, status] = preparefun(dataflow, prmflow, status);
        % set the flag prepared in status.pipeline to true
        status.pipeline.(status.nodename).prepared = true;
    elseif status.pipeline.(status.nodename).pipeline_onoff
        % default prepare of pipepool
        dataflow.pipepool.(pipenodes{ii}) = status.defaultpool;
        % I know the defaultpool is rawdata and rawhead
        % default prepare of private buffer
        dataflow.buffer.(pipenodes{ii}) = struct();
    end
end
if status.echo_onoff && flag_prepared
    fprintf(']');
end
1;

% done

status.nodename = orignodename;
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end


% function prmflow = readrawdataprepare(prmflow, status)
% % prepare of read rawdata
% prmraw = struct();
% 
% % shots
% if isfield(prmflow, 'protocol')
%     prmraw.Nshot = prmflow.protocol.shotnumber;
%     if isfield(prmflow.protocol, 'startshot')
%         prmraw.startshot = prmflow.protocol.startshot;
%     else
%         prmraw.startshot = 1;
%     end
%     if isfield(prmflow.protocol, 'endshot')
%         prmraw.endshot = prmflow.protocol.endshot;
%     else
%         prmraw.endshot = prmraw.Nshot;
%     end
%     prmraw.viewpershot = prmflow.protocol.viewnumber;
% else
%     prmraw.startshot = 1;
%     prmraw.endshot = 1;
%     prmraw.Nshot = 1;
% end
% % recount shotnum
% prmraw.Nshot = min(prmraw.endshot-prmraw.startshot+1, prmraw.Nshot);
% % copy to prmflow
% prmflow.raw = prmraw;
% 
% % .raw from protocol
% if isfield(prmflow, 'protocol')
% %     viewnumber = reconcfg.protocol.viewnumber;
%     prmflow.raw.Nviewprot = prmflow.protocol.viewperrot;
%     % viewnumber
%     prmflow.raw.Nview = prmflow.protocol.viewnumber * prmflow.raw.Nshot;
%     % I know the prmflow.protocol.viewnumber is the view number per shot for axial, and for helical only one shot once.
%     % scan
%     prmflow.raw.scan = lower(prmflow.protocol.scan);
%     % tilt
%     prmflow.raw.gantrytilt =  prmflow.protocol.gantrytilt*(pi/180);
%     % explain focal spot
%     focalspot_0x = focalspot20x(prmflow.protocol.focalspot);
%     spots = fliplr(dec2bin(focalspot_0x)=='1');
%     prmflow.raw.Nfocal = sum(spots);
%     prmflow.raw.focalspot = find(spots);
%     % NOTE: prmflow.protocol.focalspot is the name of the focalspot mode,
%     %       prmflow.raw.focalspot is the index of the focalspot(s).
% end
% 
% % .raw from calibration tables
% if isfield(prmflow.system, 'detector')
%     prmflow.raw.Nslice = prmflow.system.detector.Nmergedslice;
%     prmflow.raw.Npixel = double(prmflow.system.detector.Npixel);
% end
% 
% % initial
% status.pipepool.loadrawdata = struct();
% % no input pool fields for the node loadrawdata
% 
% end