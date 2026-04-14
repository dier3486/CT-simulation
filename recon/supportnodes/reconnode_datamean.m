function [dataflow, prmflow, status] = reconnode_datamean(dataflow, prmflow, status)
% support node, mean of rawdata
% [dataflow, prmflow, status] = reconnode_datamean(dataflow, prmflow, status);

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

% not prepared?
if ~status.pipeline.(status.nodename).prepared
    [dataflow, prmflow, status] = reconnode_datameanprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);
Nshot = nodeprm.Nshotbefore;

% pipeline_onoff
pipeline_onoff = status.currentjob.pipeline_onoff;

% prio
if pipeline_onoff
    % node prio-step
    [dataflow, prmflow, status] = nodepriostep(dataflow, prmflow, status);
    if status.currentjob.topass
        % error or pass
        return;
    end
    % private prio-step, the node type is H-A.1.N
    nextnode = status.currentjob.nextnode;
    carrynode = status.currentjob.carrynode;
    % re-initial the next pool while shot start
    if status.currentjob.pipeline.isshotstart
        % to pass?
        if status.currentjob.pipeline.readnumber <= nodeprm.viewskip
            % to pass
            status.currentjob.jobdone = 3;
            status.currentjob.topass = true;
            return;
        end
        % 1st shot start
        % dataflow.pipepool.(nextnode).WriteStart = 1;
        % dataflow.pipepool.(nextnode).ReadStart = 1;
        % dataflow.pipepool.(nextnode).WritePoint = 1;
        % dataflow.pipepool.(nextnode).ReadPoint = 1;
        dataflow.pipepool.(nextnode).WriteEnd = inf;
        dataflow.pipepool.(nextnode).ReadEnd = inf;

        1;
        status.currentjob.pipeline.Index_in(1) = status.currentjob.pipeline.Index_in(1) + nodeprm.viewskip;
        status.currentjob.pipeline.Index_out(2) = status.currentjob.pipeline.Index_out(2) - nodeprm.viewskip;
        status.currentjob.pipeline.writenumber = status.currentjob.pipeline.writenumber - nodeprm.viewskip;
    end
    % re-check the shot end
    Ishotend1 = dataflow.pipepool.(nodename).WritePoint == dataflow.pipepool.(nodename).WriteEnd + 1;
    status.currentjob.pipeline.isshotend = Ishotend1 && (~nodeprm.overshots || (dataflow.buffer.(nodename).ishot==Nshot));
    % only jobdone while it is 'really' shotend
    if status.currentjob.pipeline.isshotend
        % job done
        dataflow.pipepool.(nextnode).WriteEnd = dataflow.pipepool.(nextnode).WritePoint + ...
            status.currentjob.pipeline.writenumber - 1;
        status.currentjob.pipeline.newAvail = nodeprm.nominalnumber;
        dataflow.pipepool.(nextnode).ReadEnd = nodeprm.nominalnumber;
        % done
        status.jobdone = 1;
    else
        status.currentjob.pipeline.newAvail = 0;
        % techly pass
        status.jobdone = 4;
    end
    % ishot++
    dataflow.buffer.(nodename).ishot = dataflow.buffer.(nodename).ishot + Ishotend1;

end

% main
if pipeline_onoff
    [dataflow.pipepool.(carrynode)(1).data, dataflow.buffer.(nodename)] = ...
        datameanKernelfunction(dataflow.pipepool.(carrynode)(1).data, dataflow.pipepool.(nodename).data, ...
        dataflow.buffer.(nodename), dataflow.pipepool.(nodename), dataflow.pipepool.(nextnode));
else
    dataflow = datameanKernelfunction(dataflow, dataflow);
end

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done

% Kernel funtion
    function [dataOut, buffer] = datameanKernelfunction(dataOut, dataIn, buffer, currpool, nextpool)
        % The anonymous function is static
        debug = [];
        
        % nodeprm = prmflow.pipe.(nodename);
        Nnomi = nodeprm.nominalnumber;       

        if pipeline_onoff
            % ishot = buffer.ishot;
            plconsol = status.currentjob.pipeline;
            
            index_in = poolindex(currpool, plconsol.Index_in);
            % NviewIn = length(index_in);
            index_out = poolindex(nextpool, plconsol.Index_out);
            % NviewOut = length(index_out);

            Nstart = mod(nextpool.WritePoint - 1, Nnomi);
            dataOut.rawdata = rawdatacumsum(dataOut.rawdata, dataIn.rawdata, Nnomi, Nstart, index_out, index_in);

            % rawhead   
            if nextpool.WritePoint <= Nnomi
                IndexOut_head = min(plconsol.Index_out, Nnomi);
                index_out_head = poolindex(nextpool, IndexOut_head);
                Nhead = length(index_out_head);
                index_in_head = index_in(1:Nhead);
                dataOut.rawhead = meanskiphead(dataOut.rawhead, dataIn.rawhead, index_out_head, index_in_head);
            end

            % in shotend
            if plconsol.isshotend
                Nview = nextpool.WritePoint - nextpool.WriteStart + plconsol.writenumber;
                Nsum = ones(1, Nnomi) .* floor(Nview/Nnomi);
                Nsum(1 : mod(Nview, Nnomi)) = Nsum(1 : mod(Nview, Nnomi)) + 1;
                dataOut.rawdata = dataOut.rawdata./Nsum;
            end
            
        else % pipeline off
            if nodeprm.overshots
                dataOut.rawdata = rawdatacumsum([], dataIn.rawdata(:, 1+nodeprm.viewskip : end), Nnomi);
                Nview = size(dataIn.rawdata, 2) - nodeprm.viewskip;
                Nsum = ones(1, Nnomi) .* floor(Nview/Nnomi);
                Nsum(1 : mod(Nview, Nnomi)) = Nsum(1 : mod(Nview, Nnomi)) + 1;
                dataOut.rawdata = dataOut.rawdata./Nsum;
                dataOut.rawhead = meanskiphead(dataOut.rawhead, dataIn.rawhead, 1:Nnomi, (1:Nnomi)+nodeprm.viewskip);
            else
                viewpershot = nodeprm.viewpershotbefore;
                rawout = [];
                rawhead = struct();
                for ishot = 1 : Nshot
                    shotindex = nodeprm.startshotbefore + ishot - 1;
                    % index in 
                    startvindex = sum(viewpershot(1 : shotindex-1));
                    index_in = (1+nodeprm.viewskip : viewpershot(shotindex)) + startvindex;
                    % Nsum
                    Nview = viewpershot(shotindex) - nodeprm.viewskip;
                    Nsum = ones(1, Nnomi) .* floor(Nview/Nnomi);
                    Nsum(1 : mod(Nview, Nnomi)) = Nsum(1 : mod(Nview, Nnomi)) + 1;
                    % rawout
                    rawout = cat(2, rawout, rawdatacumsum([], dataIn.rawdata(:, index_in), Nnomi)./Nsum);
                    % rawhead
                    index_in_head = (1:Nnomi) + startvindex + nodeprm.viewskip;
                    index_out_head = (1:Nnomi) + Nnomi * (ishot-1);
                    rawhead = meanskiphead(rawhead, dataIn.rawhead, index_out_head, index_in_head);
                end
                dataOut.rawdata = rawout;
                dataOut.rawhead = rawhead;
            end

            status.jobdone = true;
        end

    end
end

function rawOut = rawdatacumsum(rawOut, rawIn, Nnomi, Nstart, index_out, index_in)

if nargin < 3
    Nnomi = size(rawOut, 2);
end
if nargin < 4
    Nstart = 0;
end
if nargin<6 || isempty(index_in)
    index_in = 1:size(rawIn, 2);
end
NviewIn = length(index_in);
if nargin<5 || isempty(index_out)
    index_out = mod((0 : NviewIn-1) + Nstart, Nnomi) + 1;
end

Np = size(rawIn, 1);
Nend = mod(Nstart + NviewIn, Nnomi);
Nmid = [min(mod(-Nstart, Nnomi), NviewIn), NviewIn - Nend];

if isempty(rawOut)
   rawOut = zeros(Np, Nnomi, 'like', rawIn); 
end

% start to normial
rawOut(:, index_out(1:Nmid(1))) = rawOut(:, index_out(1:Nmid(1))) + ...
    rawIn(:, index_in(1:Nmid(1)));
% mid
if Nmid(2) > Nmid(1)+1
    rawOut = rawOut + ...
        sum(reshape(rawIn(:, index_in(Nmid(1)+1 : Nmid(2))), Np, Nnomi, []), 3);
end
% normial to end
if Nmid(2) >= 0
    rawOut(:, index_out(Nmid(2)+1 : end)) = rawOut(:, index_out(Nmid(2)+1 : end)) + ...
        rawIn(:, index_in(Nmid(2)+1 : end));
end

end


function rawheadOut = meanskiphead(rawheadOut, rawheadIn, index_out, index_in)

for ifield = fieldnames(rawheadIn)'
    if ~isfield(rawheadOut, ifield{1})
        rawheadOut.(ifield{1}) = [];
    end
    rawheadOut.(ifield{1})(:, index_out) = rawheadIn.(ifield{1})(:, index_in);
end

end
