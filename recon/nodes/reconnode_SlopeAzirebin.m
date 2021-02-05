function [dataflow, prmflow, status] = reconnode_SlopeAzirebin(dataflow, prmflow, status)
% recon node, Axial slope-Azi rebin 
% [dataflow, prmflow, status] = reconnode_SlopeAzirebin(dataflow, prmflow, status);
% Support (X)DFS

Nshot =  prmflow.recon.Nshot;
Nslice = prmflow.recon.Nslice;
rebin = prmflow.rebin;
Nviewprot = rebin.Nviewprot;
Nreb = rebin.Nreb;
% I know it was Nviewprot/Nfocal.
idealphi = rebin.idealphi;
delta_view = pi*2/Nviewprot;
isGPU = ~isempty(status.GPUinfo);
% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Nreb, Nslice, []);
% prepare interp
if isGPU
    fAzi = gpuArray(-idealphi(:)./delta_view);
    viewindex = gpuArray(single(1:Nviewprot+1));
    pixelindex = gpuArray(single(1:Nreb)');
else
    fAzi = -idealphi(:)./delta_view;
    viewindex = single(1:Nviewprot+1);
    pixelindex = single(1:Nreb)';
end
startvindex = mod(gather(max(floor(-fAzi))+1), Nviewprot) + 1;
fAzi = mod(fAzi+viewindex([startvindex:Nviewprot 1:startvindex-1])-1, Nviewprot) + 1;

for ishot = 1:Nshot
    vindex = (1:Nviewprot)+(ishot-1)*Nviewprot;
    for islice = 1:Nslice
        % get data per slice per shot
        if isGPU
            data_islice = zeros(Nreb, Nviewprot+1, class(dataflow.rawdata), 'gpuArray');
        else
            data_islice = zeros(Nreb, Nviewprot+1, class(dataflow.rawdata));
        end
        data_islice(:, 1:Nviewprot) = squeeze(dataflow.rawdata(:, islice, vindex));
        % boundary
        data_islice(:, end) = data_islice(:, 1);
        % interp
        dataflow.rawdata(:, islice, vindex) = reshape(gather(interp2(viewindex, pixelindex, data_islice, ...
            fAzi, repmat(pixelindex, 1, Nviewprot), 'linear')), Nreb, 1, Nviewprot);
    end
end

% viewangle
viewangle = reshape(dataflow.rawhead.viewangle, Nviewprot, Nshot);
startviewangle = viewangle(startvindex, :);
dataflow.rawhead.viewangle = [viewangle(startvindex : Nviewprot, :); viewangle(1 : startvindex-1, :)];
dataflow.rawhead.viewangle = dataflow.rawhead.viewangle(:)';

prmflow.recon.startviewangle = startviewangle;
prmflow.recon.delta_view = delta_view;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end