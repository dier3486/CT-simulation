function [dataflow, prmflow, status] = reconnode_Rebin(dataflow, prmflow, status)
% recon node, rebin 
% [dataflow, prmflow, status] = reconnode_Rebin(dataflow, prmflow, status);
% Support (X)DFS, gantry tilt, NO QDO

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
    [dataflow, prmflow, status] = reconnode_rebinprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% pipeline_onoff
pipeline_onoff = status.currentjob.pipeline_onoff;

nodename = status.nodename;

% prio
if pipeline_onoff
    % node prio-step
    [dataflow, prmflow, status] = nodepriostep(dataflow, prmflow, status);
    if status.currentjob.topass
        % error or pass
        return;
    end
    carrynode = status.currentjob.carrynode;
end

% main
Nshot = prmflow.raw.Nshot;
if pipeline_onoff
    nextnode = status.pipeline.(status.nodename).nextnode;
    if isempty(nextnode) || strcmpi(nextnode, 'NULL')
        % return;
        status.jobdone = 0;
        status.errorcode = -311;
        status.errormsg = 'The rebin node cannot end in last node.';
        return
    end
    Nshot = 1;
    % skip looping the shots
else
    % create dataflow_redirect to output results
    dataflow_redirect = struct();
    dataflow_redirect.rawhead = struct();
    rebinbuffer = struct();
end

% loop the shots
for ishot = 1 : Nshot
    if pipeline_onoff
        [dataflow.pipepool.(carrynode)(1).data, dataflow.buffer.(nodename)] = ...
            rebinKernelfuntion(dataflow.pipepool.(carrynode)(1).data, dataflow.pipepool.(nodename).data, ...
            dataflow.buffer.(nodename), ishot, ...
            dataflow.pipepool.(nodename), dataflow.pipepool.(nextnode));
    else
        [dataflow_redirect, rebinbuffer] = rebinKernelfuntion(dataflow_redirect, dataflow, rebinbuffer, ishot);
    end
end
if ~pipeline_onoff
    % copy the dataflow_redirect to dataflow
    dataflow.rawdata = dataflow_redirect.rawdata;
    dataflow.rawhead = dataflow_redirect.rawhead;
end

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done
1;

% Kernel funtion
    function [dataOut, buffer] = rebinKernelfuntion(dataOut, dataIn, buffer, ishot, currpool, nextpool)
        % The anonymous function is static
        debug = [];

        % parameters
        scan = prmflow.raw.scan;
        rebin = prmflow.rebin;
        Nslice = rebin.Nslice;
        Npixel = rebin.Npixel;
        Nfocal = rebin.Nfocal;
        Nreb = rebin.Nreb;
        Nviewprot = rebin.Nviewprot;
        % delta_view = rebin.delta_view;
        % just use the Nviewprot as the viewblock
%         viewblock = prmflow.rebin.Nviewprot;
%         viewblock_focal = viewblock/Nfocal;
        viewpershot = prmflow.raw.viewpershot;
%         viewpershot_out = viewpershot./Nfocal;
        sloperebin_onoff = rebin.issloperebin;
        anglepercode = rebin.anglepercode;
        delta_anglecode = rebin.delta_anglecode; % double
        angulationzero = rebin.angulationzero; % double

        if Nfocal>1
            % prepare for DFS small disturb alignment
            interpshift = floor(rebin.DFSviewinterp);
            interpalpha = rebin.DFSviewinterp - interpshift;
        end

        
        % the console parameters for pipeline on or off
        if pipeline_onoff
            plconsol = status.currentjob.pipeline;
            NviewIn = plconsol.Index_in(2)-plconsol.Index_in(1)+1;
            % no odd view in DFS
            if mod(NviewIn, Nfocal) ~= 0
                plconsol.Index_in(2) = plconsol.Index_in(2) - mod(NviewIn, Nfocal);
                NviewIn = NviewIn - mod(NviewIn, Nfocal);
            end
            index_in = poolindex(currpool, plconsol.Index_in);
            IndexHead = [plconsol.Index_out(1) - plconsol.Do2i ...
                NviewIn/Nfocal + plconsol.Index_out(1) - plconsol.Do2i - 1];
            index_head = poolindex(nextpool, IndexHead);
            index_out = poolindex(nextpool, plconsol.Index_out);
            NviewOut = length(index_out);
            if strcmpi(scan, 'axial')
                Noutlimit = Nviewprot / Nfocal;
            else
                Noutlimit = inf;
            end
            Do2i = plconsol.Do2i;
        else
            shotindex = prmflow.raw.startshot + ishot - 1;
            % index in and out
            index_in = (1:viewpershot(shotindex)) + sum(viewpershot(1:shotindex-1));
            index_head = (1:floor(viewpershot(shotindex)/Nfocal)) + sum(floor(viewpershot(1:shotindex-1)./Nfocal));
            NviewIn = viewpershot(shotindex);
            index_out = index_head;
            NviewOut = length(index_out);
            Noutlimit = inf;
            Do2i = 0;
            if ~isfield(dataOut, 'rawdata') || isempty(dataOut.rawdata)
                % mostly
                dataOut.rawdata = zeros(Nreb*Nslice, rebin.NviewOut, 'like', dataIn.rawdata);
            end
        end
        % record NviewIn
        prmflow.rebin.viewread = prmflow.rebin.viewread + NviewIn;

        % step1 rawhead
        angle0 = double(dataIn.rawhead.Angle_encoder(index_in(1)));     % uint32 -> int
        % align due to anglecode-distort
        if rebin.viewanglealign && isfield(buffer, 'startanglecode')
            angledis = mod(angle0 - buffer.startanglecode, delta_anglecode);
            if angledis > delta_anglecode/2
                angledis = angledis-delta_anglecode;
            end
        else
            angledis = 0;
        end
        viewdis = -double(angledis)*anglepercode + rebin.DFSviewshift;
        % the rawhead after rebin
        dataOut.rawhead = rawheadofrebin(dataOut.rawhead, dataIn.rawhead, index_head, index_in, Nfocal, angledis, viewdis);
        1;
        
        % view angles shall use the rawhead.Angle_encoder
        if ~isfield(buffer, 'startanglecode')
            if pipeline_onoff && plconsol.isshotstart
                buffer.startanglecode = angle0;
            elseif strcmpi(scan, 'axial') && ishot==1
                buffer.startanglecode = angle0;
            end
        end
        
        % step2 Ratial
        NviewIn_focal = NviewIn/Nfocal + (Nfocal-1)*2;
        viewindex = cast(1:NviewIn_focal, 'like', dataIn.rawdata);
        dataRatial = zeros(Nreb, Nslice, NviewIn_focal+2, 'like', dataIn.rawdata);
        % slope ratial
        if sloperebin_onoff
            % viewangle to be used in Y-shift
            viewangle_head = cast(double(dataIn.rawhead.Angle_encoder(index_in(1:Nfocal:end))).*anglepercode + ...
                angulationzero, 'like', dataIn.rawdata);
            if Nfocal>1
                viewangle_head = ...
                    [viewangle_head(1)*2 - viewangle_head(2)  viewangle_head  viewangle_head(end)*2 - viewangle_head(end-1)];
            end
            viewangle_head = viewangle_head + rebin.DFSviewshift;
            % 1/dfan to use
            dfanun1 = 1/rebin.dfan;
        end

        for islice = 1:Nslice
            index_is = (1:Npixel) + Npixel*(islice-1);
            % get the ith-slice raw data
            if Nfocal==1
                data_islice = dataIn.rawdata(index_is, index_in);
            else
                % DFS small disturb alignment
                data_islice = zeros(Npixel, NviewIn+4, 'like', dataIn.rawdata);
                data_islice(:, 3:end-2) = dataIn.rawdata(index_is, index_in);
                for ifocal = 1:Nfocal
                    data_islice(:, (ifocal - interpshift(ifocal)*Nfocal) : Nfocal : end-(interpshift(ifocal)+1)*Nfocal) = ...
                        data_islice(:, ifocal:Nfocal:end-2).*(1-interpalpha(ifocal)) + ...
                        data_islice(:, Nfocal+ifocal:Nfocal:end).*interpalpha(ifocal);
                end
                if ~strcmpi(scan, 'axial')
                    % to add a boundary for helical
                    % not so necessary
                end
            end
            % reorder by focals
            data_islice = reshape(permute(reshape(data_islice, Npixel, Nfocal, []), [2 1 3]), Npixel*Nfocal, []);
            if sloperebin_onoff
                % interp to equal fanangle (sloperebin on)
                data_islice = interp1(data_islice, rebin.faninterpkern(:, islice), 'linear', 'extrap');
                % interp to ideal fanangle (equal distance)
                phi3 = rebin.idealphi + asin( sin(viewangle_head+rebin.idealphi).*rebin.Yshift(islice) );
                phi3 = phi3.*dfanun1 + rebin.midchannel;
                1;
                dataRatial(:, islice, 2:end-1) = interp2(data_islice, repmat(viewindex, Nreb, 1), phi3, 'linear', 0);
            else
                % I know the faninterpkern is ideal fanangle while sloperebin off
                dataRatial(:, islice, 2:end-1) = interp1(data_islice, rebin.faninterpkern(:, islice), 'linear', 0);
            end
        end
        if sloperebin_onoff
            % interp Z (sloperebin on)
            Zgrid = repmat(rebin.Zgrid, 1, 1,NviewIn_focal);
            Xgrid = repmat(cast((1:Nreb)', 'like', dataIn.rawdata), 1, Nslice, NviewIn_focal);
            viewindex_grid = repmat(reshape(viewindex, 1, 1, []), Nreb, Nslice, 1);
            dataRatial(:,:, 2:end-1) = interp3(dataRatial(:,:, 2:end-1), Zgrid, Xgrid, viewindex_grid, 'linear', 0);
        end
        
        1;
        % Azi-rebin
        pixelindex = cast((1:Nreb)', 'like', dataIn.rawdata);
        viewindex_out = (1 : NviewOut) + Do2i + Nfocal;
        Azidis = single(angledis)/single(delta_anglecode)/Nfocal;
        if isfield(rebin, 'lowaccuracyshift')
            Azidis = round(Azidis.*2^rebin.lowaccuracyshift).*2^(-rebin.lowaccuracyshift);
        end
        fAzi = rebin.idealfAzi - Azidis + cast(viewindex_out, 'like', dataIn.rawdata);
        if ~pipeline_onoff && strcmpi(scan, 'axial')
            % boundary condition of axial
            if Nfocal == 1
                dataRatial(:,:, end) = dataRatial(:,:, 2);
            else
                dataRatial(:,:, 3) = dataRatial(:,:, 3) + dataRatial(:,:, end-1);
                dataRatial(:,:, end-2) = dataRatial(:,:, end-2) + dataRatial(:,:, 2);
                dataRatial(:,:, end-1) = dataRatial(:,:, 3);
            end
            fAzi = mod(fAzi-1-Nfocal, Nviewprot/Nfocal) + 1+Nfocal;
        end

        for islice = 1:Nslice
            index_is = (1:Nreb) + Nreb*(islice-1);
            if NviewOut <= Noutlimit
                dataOut.rawdata(index_is, index_out) = dataOut.rawdata(index_is, index_out) + ...
                    reshape(interp2(squeeze(dataRatial(:, islice, :)), ...
                    fAzi, repmat(pixelindex, 1, NviewOut), 'linear', 0), Nreb, NviewOut);
            else
                % while the index_out is not unique
                dataOut.rawdata(index_is, index_out(1:Noutlimit)) = dataOut.rawdata(index_is, index_out(1:Noutlimit)) + ...
                    reshape(interp2(squeeze(dataRatial(:, islice, :)), ...
                    fAzi(:, 1:Noutlimit), repmat(pixelindex, 1, Noutlimit), 'linear', 0), Nreb, Noutlimit);
                dataOut.rawdata(index_is, index_out(Noutlimit+1:end)) = dataOut.rawdata(index_is, index_out(Noutlimit+1:end)) + ...
                    reshape(interp2(squeeze(dataRatial(:, islice, :)), ...
                    fAzi(:, Noutlimit+1:end), repmat(pixelindex, 1, NviewOut-Noutlimit), 'linear', 0), Nreb, NviewOut-Noutlimit);
            end
        end

        if ~pipeline_onoff
            status.jobdone = true;
        end
        % rebinKernelfuntion END
    end
end


function rawheadOut = rawheadofrebin(rawheadOut, rawheadIn, index_out, index_in, Nfocal, angledis, viewdis)
% 'copy' the input rawhead to output
% support viewangle/2 for DFS and other fucntions

headfields = fieldnames(rawheadIn);
Nfields = length(headfields);

for ifield = 1:Nfields
    datasize = size(rawheadIn.(headfields{ifield}), 1);
    if ~isfield(rawheadOut, headfields{ifield})
        rawheadOut.(headfields{ifield}) = zeros(datasize, 0 , 'like', rawheadIn.(headfields{ifield}));
    end
    switch headfields{ifield}
        case 'refblock'
            rawheadOut.refblock(:, index_out) = squeeze(any(reshape(rawheadIn.refblock(:, index_in), datasize, Nfocal, []), 2));
        case 'viewangle'
            viewdis = cast(viewdis, 'like', rawheadIn.viewangle);
            rawheadOut.viewangle(:, index_out) = rawheadIn.viewangle(index_in(1:Nfocal:end)) + viewdis;
            % hard code, shift the viewangle by pi/2
            rawheadOut.viewangle(:, index_out) = rawheadOut.viewangle(:, index_out) + pi/2;
            % path dependence operation, but we will abandon in using the viewangle.
        case 'Angle_encoder'
            rawheadOut.Angle_encoder(:, index_out) = rawheadIn.Angle_encoder(:, index_in(1:Nfocal:end)) - angledis;
        case {'mA', 'KV'}
            rawheadOut.(headfields{ifield})(:, index_out) = mean(reshape(rawheadIn.(headfields{ifield})(index_in), Nfocal, []), 1);
        case {'Shot_Number', 'Reading_Number', 'Table_encoder', 'Table_gear'}
            rawheadOut.(headfields{ifield})(:, index_out) = rawheadIn.(headfields{ifield})(:, index_in(1:Nfocal:end));
        case 'Shot_Start'
            rawheadOut.Shot_Start(:, index_out) = sum(reshape(rawheadIn.Shot_Start(index_in), Nfocal, []), 1);
            % Warn: could lost the shot end while the viewnumber is not common with the Nfocal.
        otherwise
            % remove the useless fields
            rawheadOut = rmfield(rawheadOut, headfields{ifield});
    end
end

end
