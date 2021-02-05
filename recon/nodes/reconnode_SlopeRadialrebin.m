function [dataflow, prmflow, status] = reconnode_SlopeRadialrebin(dataflow, prmflow, status)
% recon node, Axial slope-radial rebin 
% [dataflow, prmflow, status] = reconnode_SlopeRadialrebin(dataflow, prmflow, status);
% Support (X)DFS

Nview =  prmflow.recon.Nview;
Nslice = prmflow.recon.Nslice;
Npixel = prmflow.recon.Npixel;
% Nviewprot = prmflow.recon.Nviewprot;
Nfocal = prmflow.recon.Nfocal;
rebin = prmflow.rebin;
viewblock = rebin.viewblock;
Nreb = rebin.Nreb;
isGPU = ~isempty(status.GPUinfo);

dataout = zeros(Nreb, Nslice, Nview/Nfocal, class(dataflow.rawdata));
if isGPU
    % GPU
    idealphi = gpuArray(rebin.idealphi);
    pixelindex = gpuArray(single(1:Npixel*Nfocal));
    faninterpkern = gpuArray(rebin.faninterpkern);
    viewangle = gpuArray(dataflow.rawhead.viewangle);
    Yshift = gpuArray(rebin.Yshift);
    dfanun1 = gpuArray(1/rebin.dfan);
    midU = gpuArray(rebin.midchannel);
    sliceindex = gpuArray(single(1:Nslice));
    viewindex = repmat(reshape(gpuArray(single(1:viewblock)), 1, 1, []), Nreb, Nslice, 1);
    Zgrid = repmat(gpuArray(rebin.Zgrid), 1, 1,viewblock);
    Xgrid = repmat(gpuArray(single(1:Nreb)'), 1, Nslice, viewblock);
else
    idealphi = rebin.idealphi;
    pixelindex = single(1:Npixel*Nfocal);
    faninterpkern = rebin.faninterpkern;
    viewangle = dataflow.rawhead.viewangle;
    Yshift = rebin.Yshift;
    dfanun1 = 1/rebin.dfan;
    midU = rebin.midchannel;
    sliceindex = single(1:Nslice);
    viewindex = repmat(reshape(single(1:viewblock), 1, 1, []), Nreb, Nslice, 1);
    Zgrid = repmat(rebin.Zgrid, 1, 1, viewblock);
    Xgrid = repmat(single(1:Nreb)', 1, Nslice, viewblock);
end

dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nview);

for iblk = 1:Nview/viewblock
    vindex_iblk = (1:viewblock)+(iblk-1)*viewblock;
    if isGPU
        data_blk = gpuArray(dataflow.rawdata(:, vindex_iblk));
    else
        data_blk = dataflow.rawdata(:, vindex_iblk);
    end
    data_blk = reshape(permute(reshape(data_blk, Npixel*Nslice, Nfocal, []), [2, 1, 3]), Npixel*Nfocal, Nslice, []); 
    for islice = 1:Nslice
        data_blk(:, islice, :) = reshape(interp1(pixelindex, reshape(data_blk(:, islice, :), Npixel, viewblock), ...
            faninterpkern(:, islice), 'linear', 'extrap'), Npixel, 1, viewblock);
    end
    
    phi3 = idealphi + reshape(asin(reshape(sin(viewangle(vindex_iblk)+idealphi), [], 1)*Yshift), Nreb, []);
    phi3 = phi3.*dfanun1 + midU;
    phi3 = permute(reshape(phi3, Nreb, viewblock, Nslice), [1 3 2]);
    data2_blk = interp3(data_blk, repmat(sliceindex, Nreb, 1, viewblock), phi3, viewindex, 'linear', 0);
    dataout(:,:,vindex_iblk) = gather(interp3(data2_blk, Zgrid, Xgrid, viewindex, 'linear', 0));
end

% return data
dataflow.rawdata = dataout;
% view angle
dataflow.rawhead.viewangle = dataflow.rawhead.viewangle(1:Nfocal:end);

% prm
prmflow.recon.Npixel = rebin.Nreb;
prmflow.recon.Nviewprot = rebin.Nviewprot;
prmflow.recon.Nview = rebin.Nviewprot*prmflow.recon.Nshot;
prmflow.recon.midchannel = rebin.midU_phi;
prmflow.recon.delta_d = rebin.delta_d;
prmflow.recon.delta_z = rebin.delta_z;
prmflow.recon.SID = rebin.SID;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end