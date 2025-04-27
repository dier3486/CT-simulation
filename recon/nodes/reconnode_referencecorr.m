function [dataflow, prmflow, status] = reconnode_referencecorr(dataflow, prmflow, status)
% recon node, (air) reference correction
% [dataflow, prmflow, status] = reconnode_referencecorr(dataflow, prmflow, status);

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
    [dataflow, prmflow, status] = reconnode_referenceprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% pipeline_onoff
pipeline_onoff = status.currentjob.pipeline_onoff;

% prio
if pipeline_onoff
    % node prio-step
    [dataflow, prmflow, status] = nodepriostep(dataflow, prmflow, status);
    if status.jobdone==0 || status.jobdone>=3
        % error or pass
        return;
    end
end

% main
nodename = status.nodename;
Nshot = prmflow.raw.Nshot;
if pipeline_onoff
    nextnode = status.currentjob.nextnode;
    carrynode = status.currentjob.carrynode;
    if isempty(nextnode) || strcmpi(nextnode, 'NULL')
        % return;
        status.jobdone = 0;
        status.errorcode = -311;
        status.errormsg = 'The referencecorr node cannot end in last node.';
        return
    end
    Nshot = 1;
    % skip looping the shots
end

% loop the shots
for ishot = 1 : Nshot
    if pipeline_onoff
        [dataflow.pipepool.(carrynode).data, dataflow.buffer.(nodename).reflast] = ...
            refcorrKernelfuntion(dataflow.pipepool.(carrynode).data, ishot, dataflow.pipepool.(nextnode));
    else
        dataflow = refcorrKernelfuntion(dataflow, ishot);
    end
end

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done

    % Kernel funtion
    function [data, reflast] = refcorrKernelfuntion(data, ishot, nextpool)
        % The anonymous function is static
        debug = [];

        % parameters from prmflow
        Nfocal = prmflow.raw.Nfocal;
        blockwindow = prmflow.correction.air.blockwindow;
        refpixel = prmflow.correction.air.refpixelindex;
        referrcut = reshape(prmflow.correction.air.referrcut, [], Nfocal);
        sliceindependent = prmflow.correction.air.sliceindependent;
        Npixel = prmflow.raw.Npixel;
        Nslice = prmflow.raw.Nslice;
        startshot = prmflow.raw.startshot;
        viewpershot = prmflow.raw.viewpershot;
        scan = lower(prmflow.raw.scan);
        
        % the console parameters for pipeline or not
        if pipeline_onoff
            plconsol = status.currentjob.pipeline;
            % set Index_corr
            plconsol.Index_corr = [max(nextpool.AvailPoint+1, nextpool.ReadPoint), nextpool.AvailPoint + plconsol.newAvail];
            Index_union = [min(plconsol.Index_in(1), plconsol.Index_corr(1)), max(plconsol.Index_in(2), plconsol.Index_corr(2))];
            % indexes to use
            index_ref = poolindex(nextpool, Index_union);
            index_out = poolindex(nextpool, plconsol.Index_out);
            Ncorr = plconsol.Index_corr(2) - plconsol.Index_corr(1) + 1;
            index_refblk = max(1, 1 + plconsol.Index_in(1) - plconsol.Index_corr(1)) : ...
                plconsol.Index_in(2) + 1 - min(plconsol.Index_in(1), plconsol.Index_corr(1));
            index_rawref = max(1, 1 + plconsol.Index_corr(1) - plconsol.Index_in(1)) : ...
                    plconsol.Index_corr(2) + 1 - min(plconsol.Index_corr(1), plconsol.Index_in(1));
            index_corr = poolindex(nextpool, plconsol.Index_corr);
            Do2i = plconsol.Do2i;
        else
            shotindex = startshot + ishot - 1;
            % indexes to use
            index_ref = (1:viewpershot(shotindex)) + sum(viewpershot(1:shotindex-1));
            switch scan
                case 'axial'
                    index_out = [index_ref(end-blockwindow+1 : end)  index_ref  index_ref(1 : blockwindow)];
                case {'helical', 'static', 'halfaxial'}
                    index_out = index_ref;
                otherwise
                    status.jobdone = false;
                    status.errorcode = 1;
                    status.errormsg = sprintf('Not support scan mode %s in air correcion', prmflow.raw.scan);
                    return;
            end
            Ncorr = viewpershot(shotindex);
            index_refblk = 1 : Ncorr;
            index_rawref = 1 : Ncorr;
            index_corr = index_ref;
            Do2i = 0;
        end
        
        % get the refernce mean and refernce error by the reference pixels
        [rawref, refblk_new] = airreference(data.rawdata(:, index_ref), referrcut, refpixel, Nslice, sliceindependent);

        % reference block
        data.rawhead.refblock = isrefblocked(data.rawhead.refblock, refblk_new(:, index_refblk), ...
            blockwindow, scan, index_out, Do2i);
        
        % reflast
        if pipeline_onoff
            reflast = dataflow.buffer.(nodename).reflast;
            % re-initial the last reference
            if status.currentjob.pipeline.isshotstart
                reflast(:) = {[]};
            end
        else
            reflast = cell(1, Nfocal);
        end

        % reference correction
        if Ncorr > 0
            % first reading number
            reading1 = data.rawhead.Reading_Number(index_corr(1));
            % loop the focals (DFS) to correct the rawdata
            for ifocal = 1:Nfocal
                % the first view is focal 1 or focal 2? depending on the first reading number
                ifocal_index = (mod(reading1(1), 1) ~= (ifocal-1)) + 1;
                index_rawref_ifc = index_rawref(ifocal_index:Nfocal:Ncorr);
                index_corr_ifc = index_corr(ifocal_index:Nfocal:Ncorr);
                % referencecorr
                [rawref_ifc, reflast{ifocal}] = referencecorr(rawref(:, index_rawref_ifc), ...
                        data.rawhead.refblock(:, index_corr_ifc), reflast{ifocal});
                % add back
                if sliceindependent
                    data.rawdata(:, index_corr_ifc) = data.rawdata(:, index_corr_ifc) - repelem(rawref_ifc, Npixel, 1);
                else
                    data.rawdata(:, index_corr_ifc) = data.rawdata(:, index_corr_ifc) - rawref_ifc;
                end
            end 
        end

        if ~pipeline_onoff
            status.jobdone = true;
        end
        % refblockKernelfuntion END
    end

end

