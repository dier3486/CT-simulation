    function [dataflow, prmflow, status] = reconnode_offfocalcorr(dataflow, prmflow, status)
% recon node, off-focal correction
% [dataflow, prmflow, status] = reconnode_offfocalcorr(dataflow, prmflow, status);
% only for axial

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
    [dataflow, prmflow, status] = reconnode_offfocalprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

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
        status.errormsg = 'The offfocalcorr node cannot end in last node.';
        return
    end
    Nshot = 1;
    % skip looping the shots
end

% loop the shots
for ishot = 1 : Nshot
    if pipeline_onoff
        [dataflow.pipepool.(carrynode).data, dataflow.buffer.(nodename).offspacepool] = ...
            offfocalcorrKernelfuntion(dataflow.pipepool.(carrynode).data, ishot, ...
            dataflow.pipepool.(nextnode), dataflow.buffer.(nodename).offspacepool);
    else
        dataflow = offfocalcorrKernelfuntion(dataflow, ishot);
    end
end

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done

    % off-focall correction Kernelfuntion
    function [data, offspacepool] = offfocalcorrKernelfuntion(data, ishot, nextpool, offspacepool)
        % The anonymous function is static
        debug = [];
        
        % parameters to use in prmflow
        Nslice = prmflow.raw.Nslice;
        Npixel = prmflow.raw.Npixel;
        Nviewprot = prmflow.raw.Nviewprot;
        Nfocal = prmflow.raw.Nfocal;
        startshot = prmflow.raw.startshot;
        viewpershot = prmflow.raw.viewpershot;
        scan = lower(prmflow.raw.scan);
        % prepared parameters
        prmoff = prmflow.correction.offfocal;
        slicemerge = prmoff.slicemerge;
        Nslicemerged = Nslice/slicemerge;
        Noffsample = prmoff.Noffsample;
        viewsparse = prmoff.viewsparse;

        if pipeline_onoff
            plconsol = status.currentjob.pipeline;
            % prio step of inner pool
            % 'copy' the next pool to currpool
            offinputpool = poolmirror(nextpool);
            % the offinputpool is a mirror of nextpool used to input the data from the nextpool to offspacepool
            % hard fix of offinputpool
            offinputpool.ReadPoint = nextpool.WritePoint;
            offinputpool.ReadStart = nextpool.WriteStart;
            offinputpool.ReadEnd = nextpool.WriteEnd;
            offinputpool.WritePoint = nextpool.WritePoint + status.currentjob.pipeline.readnumber;
            % currpool.WriteEnd = dataflow.pipepool.(nextnode).WriteEnd;
            offinputpool.AvailPoint = nextpool.WritePoint + status.currentjob.pipeline.readnumber - 1;
            offinputpool.isshotstart = plconsol.isshotstart;
            % pool-priostep
            jobraw2off = struct();
            [offinputpool, offspacepool, jobraw2off] = ...
                poolpriostep(offinputpool, offspacepool, prmoff.pipeline.raw2off, jobraw2off);
            % jobraw2off.Index_in = status.currentjob.pipeline.Index_in;
            % indexes of the input data
            index_renew = poolindex(offinputpool, jobraw2off.Index_in);
            Nrenew = length(index_renew);
        else
            shotindex = startshot + ishot - 1;
            shotstartview = sum(viewpershot(1:shotindex-1));
            % indexes of the input data
            index_renew = (1:viewpershot(shotindex)) + shotstartview;
            Nrenew = length(index_renew);
        end

        % step1 exp (on 2)
        data.rawdata(:, index_renew) = 2.^(-data.rawdata(:, index_renew));

        % step2, merge slices
        Aoff = zeros(Npixel, Nslicemerged, Nrenew + 2*Nfocal, 'like', data.rawdata);
        if slicemerge>1
            if ~slicezebra
                Aoff(:,:, Nfocal+1:end-Nfocal) = squeeze(mean(reshape(data.rawdata(:, index_renew), ...
                    Npixel, slicemerge, Nslicemerged, Nrenew), 2));
            else
                Aoff(:,:, Nfocal+1:end-Nfocal) = reshape(mean(reshape(data.rawdata(:, index_renew), ...
                    Npixel, 2, slicemerge, Nslicemerged/2, Nrenew), 3), Npixel, Nslicemerged, Nrenew);
            end
        else
            Aoff(:,:, Nfocal+1:end-Nfocal) = reshape(data.rawdata(:, index_renew), Npixel, Nslicemerged, Nrenew);
        end
        % boundary condition of Aoff
        switch scan
            case 'axial'
                if ~pipeline_onoff
                    Aoff(:,:, 1:Nfocal) =  Aoff(:,:, end-Nfocal*2+1:end-Nfocal);
                    Aoff(:,:, end-Nfocal+1:end) =  Aoff(:,:, Nfocal+1:Nfocal*2);
                end
            case {'helical', 'halfaxial'}
                if pipeline_onoff
                    if jobraw2off.isshotstart
                        % helical shot start
                        Aoff(:,:, 1:Nfocal) = Aoff(:,:, Nfocal+1:Nfocal*2);
                    end
                    if jobraw2off.isshotend
                        % helical shot end
                        Aoff(:,:, end-Nfocal+1:end) = Aoff(:,:, end-Nfocal*2+1:end-Nfocal);
                    end
                else
                    Aoff(:,:, 1:Nfocal) = Aoff(:,:, Nfocal+1:Nfocal*2);
                    Aoff(:,:, end-Nfocal+1:end) = Aoff(:,:, end-Nfocal*2+1:end-Nfocal);
                end
            otherwise
                0;
        end

        % step3, resample to off-focal space
        if pipeline_onoff
            % indexes of the data to update in off-focal space
            index_out = poolindex(offspacepool, jobraw2off.Index_out);
            Nviewoff = length(index_out);
            Nlimit = offspacepool.poolsize / Nfocal;
            index_voff = ((0 : Nviewoff-1) + jobraw2off.Do2i).*viewsparse + 1;
        else
            % indexes of the data to update in off-focal space
            index_voff = (0 : prmoff.offendview-prmoff.offstartview).*prmoff.viewsparse + prmoff.offstartview;
            Nviewoff = prmoff.Nviewoff;
            Nlimit = inf;
            index_out = 1 : Nviewoff;
            % ini rawoff
            offspacepool.data = struct();
            offspacepool.data.rawdata = zeros(Noffsample*Nslicemerged, Nviewoff + Nfocal, 'single');
        end
        % flying focal
        for ifocal = 1:Nfocal
            index_voff(ifocal:Nfocal:end) = (index_voff(ifocal:Nfocal:end) - 1)./prmoff.viewsparse + 1 - ifocal;
        end
        index_voff = index_voff./Nfocal.*prmoff.viewsparse + 1;

        % interpolation position (Df)
        if ~pipeline_onoff && strcmpi(scan, 'axial')
            Df = mod(index_voff - prmoff.rawinterp2phi - 1, Nviewprot/Nfocal) + 2;
        else
            Df = index_voff - prmoff.rawinterp2phi + 1;
            Df(Df<1) = 1;   Df(Df>floor(Nrenew/Nfocal)+2) = floor(Nrenew/Nfocal)+2;
        end
        
        for ifocal = 1:Nfocal
            % the number and index of the ith-focal
            index_ifocal = ifocal:Nfocal:Nviewoff;
            index_out_ifocal = index_out(index_ifocal);
            Nviewoff_ifocal = length(index_ifocal);
            % loop the slices to do the interpolations
            for islice = 1:Nslicemerged
                sampleindex = (1:Noffsample) + (islice-1).*Noffsample;
                if Nviewoff_ifocal <= Nlimit
                    offspacepool.data.rawdata(sampleindex, index_out_ifocal) = ...
                        offspacepool.data.rawdata(sampleindex, index_out_ifocal) + ...
                        interp2(squeeze(Aoff(:, islice, ifocal:Nfocal:end)) , Df(:, index_ifocal), ...
                        repmat(prmoff.rawinterp2t(:, islice), 1, Nviewoff_ifocal), 'linear', 0);
                else
                    % the index_out could overlap in axial
                    offspacepool.data.rawdata(sampleindex, index_out_ifocal(1:Nlimit)) = ...
                        offspacepool.data.rawdata(sampleindex, index_out_ifocal(1:Nlimit)) + ...
                        interp2(squeeze(Aoff(:, islice, ifocal:Nfocal:end)) , Df(:, index_ifocal(1:Nlimit)), ...
                        repmat(prmoff.rawinterp2t(:, islice), 1, Nlimit), 'linear', 0);
                    % the overlap part
                    offspacepool.data.rawdata(sampleindex, index_out_ifocal(Nlimit+1:end)) = ...
                        offspacepool.data.rawdata(sampleindex, index_out_ifocal(Nlimit+1:end)) + ...
                        interp2(squeeze(Aoff(:, islice, ifocal:Nfocal:end)) , Df(:, index_ifocal(Nlimit+1:end)), ...
                        repmat(prmoff.rawinterp2t(:, islice), 1, Nviewoff_ifocal-Nlimit), 'linear', 0);
                end
            end
        end

        % step4, convolution
        if pipeline_onoff
            IndexConv = [max(offspacepool.ReadPoint, offspacepool.AvailPoint + 1), ...
                         offspacepool.AvailPoint + jobraw2off.newAvail];
            index_conv = poolindex(offspacepool, IndexConv);
        else
            index_conv = 1 : Nviewoff;
        end
        % do convolution by ifft(fft(rawdata)*offkernel)).
        offspacepool.data.rawdata(:, index_conv) = reshape(ifft(fft(reshape(offspacepool.data.rawdata(:, index_conv), ...
            Noffsample, []), Noffsample).*prmoff.offkernel), Noffsample*Nslicemerged, []);

        if ~pipeline_onoff && strcmpi(scan, 'axial')
            % boundary condition of offspace for axial
            offspacepool.data.rawdata(:, end-Nfocal+1 : end) = offspacepool.data.rawdata(:, 1 : Nfocal);
        end

        % post step of inner pool
        if pipeline_onoff
            % move the points
            [~, offspacepool] = movepointsaftercopy([], offspacepool, [], jobraw2off.writenumber);
            % I know the offinputpool is used out
            % new available views
            offspacepool.AvailPoint = offspacepool.AvailPoint + jobraw2off.newAvail;
        end
        
        % step5, resample back to raw space
        if pipeline_onoff
            % pior step of off to raw
            offoutputpool = poolmirror(nextpool);
            % the offoutputpool is a mirror of nextpool used to output the data to the nextpool from offspacepool
            % hard fix of offoutputpool
            offoutputpool.WritePoint = max(nextpool.AvailPoint + 1, nextpool.ReadPoint);
            offoutputpool.WriteStart = nextpool.ReadStart;
            if offoutputpool.circulatemode
                offoutputpool.WriteEnd = nextpool.ReadEnd;
            else
                offoutputpool.WriteEnd = nextpool.WriteEnd - prmoff.pipeline.off2raw.viewrely_out(2);
                offoutputpool.WriteEnd = floor((offoutputpool.WriteEnd - offoutputpool.WriteStart + 1) / viewsparse) ...
                    * viewsparse - 1 + offoutputpool.WriteStart;
            end
            if ~jobraw2off.isshotstart
                % to forget the ReadStart
                offspacepool.ReadStart = -inf;
            end
            offoutputpool.isshotstart = plconsol.isshotstart;
            % pool-priostep
            joboff2raw = struct();
            [offspacepool, offoutputpool, joboff2raw] = ...
                poolpriostep(offspacepool, offoutputpool, prmoff.pipeline.off2raw, joboff2raw);
            % but need to hard fix somthing,
            % to replace the Index_out
            joboff2raw.Index_out = status.currentjob.pipeline.Index_out;
            % Mostly those Index_out are same, but not while the view number is not common with the viewsparse in shot end.
            % to fix negative read point
            if joboff2raw.Index_in(1) < 1
                % forced to move a negative or 0 Index_in(1) to 1
                joboff2raw.Do2i = joboff2raw.Do2i - (1 - joboff2raw.Index_in(1))*viewsparse;
                joboff2raw.Index_in(1) = 1;
                % The H-H.1.S could draw the ReadPoint to negative while the readnumber is too small in shot start. Not a bug,
                % and we can set a big enough minlimit to avoid that.
            end
            
            if joboff2raw.jobdone==3 || joboff2raw.jobdone==6
                % pass the left steps
                return;
                % The A-A.1.S could lay in here, not a bug just pass is ok.
            end
            
            % indexes of the data to use in off-focal space
            index_in = poolindex(offspacepool, joboff2raw.Index_in);
            Nviewoff = length(index_in);
            % indexes of the data to update in next pool
            index_fix = poolindex(offoutputpool, joboff2raw.Index_out);
            Nfix = length(index_fix);
            index_vraw = ((0 : Nfix-1) + joboff2raw.Do2i)./viewsparse + 1;
        else
            index_in = 1 : Nviewoff + Nfocal;
            Nfix = Nrenew;
            index_vraw = (0:Nfix-1)./viewsparse - prmoff.offstartview + 2;
            index_fix = (1 : Nrenew) + shotstartview;
        end
        % flying focal
        for ifocal = 1:Nfocal
            index_vraw(ifocal:Nfocal:end) = (index_vraw(ifocal:Nfocal:end) - ifocal)./Nfocal + 1;
        end

        % interp off space back to raw space
        Afix = zeros(Npixel, Nslicemerged, Nfix, 'single');
        for ifocal = 1 : Nfocal
            % number and index of ith-focal
            Nviewoff_focal = floor((Nviewoff-ifocal)/Nfocal) + 1;
            index_ifocal = ifocal:Nfocal:Nfix;
            Nfix_ifocal = length(index_ifocal);
            % loop the slices
            for islice = 1:Nslicemerged
                sampleindex = (1:Noffsample) + (islice-1).*Noffsample;
                Dfb = index_vraw + prmoff.tinterp2phi(:, islice);
                if ~pipeline_onoff && strcmpi(scan, 'axial')
                    Dfb = mod(Dfb - 1, Nviewoff/Nfocal) + 1;
                else                   
                    Dfb(Dfb<1) = 1;   Dfb(Dfb>Nviewoff_focal) = Nviewoff_focal;
                end
                Afix(:, islice, index_ifocal) = ...
                    interp2(squeeze(offspacepool.data.rawdata(sampleindex, index_in(ifocal:Nfocal:end))), ...
                    Dfb(:, index_ifocal), repmat(prmoff.tinterp2raw(:, islice), 1, Nfix_ifocal), 'linear', 0);

            end
        end

        if pipeline_onoff
            % post step of off to raw
            [offspacepool, ~] = movepointsaftercopy(offspacepool, [], joboff2raw.readnumber);
            % I know the offoutputpool is used out
            % recycle the offspacepool
            offspacepool = poolrecycle(offspacepool);
        end

        % step6, raw space operators
        % measure scale
        Afix = Afix.*prmoff.Dphiscale;
        Afix = real(Afix) + imag(Afix).*prmoff.Dphiscale_odd;
        % permute Afix to move the slice to dim 1
        Afix = reshape(permute(Afix, [2 1 3]), Nslicemerged, []);
        % inverse slice merge
        if slicemerge>1
            Afix = repelem(Afix, slicemerge, 1);
            if slicezebra
                Afix = reshape(permute(reshape(Afix, slicemerge, 2, []), [2 1 3]), Nslice, []);
            end
        end
        % Z cross
        Afix = prmoff.crsMatrix * Afix;
%         if ~prmoff.slicezebra
%             Afix = prmoff.crsMatrix * Afix;
%         else
%             Afix(1:2:end, :) = prmoff.crsMatrix * Afix(1:2:end, :);
%             Afix(2:2:end, :) = prmoff.crsMatrix * Afix(2:2:end, :);
%         end
        
        1;
        % step7, fix rawdata
        % reshape
        Afix = reshape(permute(reshape(Afix, Nslice, Npixel, Nfix), [2 1 3]), Npixel*Nslice, Nfix);
        % fix to rawdata
        Afix = data.rawdata(:, index_fix) - Afix;
        % minimum limit
        minintensity = prmflow.correction.offfocal.minintensity;
        Afix(Afix<minintensity) = minintensity;
        % log2
        data.rawdata(:, index_fix) = -log2(Afix);
        
        % Done
        if ~pipeline_onoff
            status.jobdone = true;
        end
        % offfocalcorrKernelfuntion END
    end


end

