% 'slope' rebin test code
% Atilt transform

load('E:\data\simulation\TM\test\tiltrb_test1.mat');
prmflow.recon.gantrytilt = prmflow.protocol.gantrytilt*(pi/180);

gpuDevice;

% Nviewprot = prmflow.recon.Nviewprot;
focalspot = prmflow.recon.focalspot;
focalposition = prmflow.system.focalposition(focalspot, :);

gantrytilt = prmflow.protocol.gantrytilt*(pi/180);

Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nview = prmflow.recon.Nview;
Nviewprot = prmflow.recon.Nviewprot;
Nshot = prmflow.recon.Nshot;
% detector
detector = prmflow.system.detector;

[fanangles0, focalangle0] = detpos2fanangles(detector.position, focalposition);
fanangles0 = reshape(fanangles0, Npixel, Nslice) - focalangle0;

delta_d = detector.hx_ISO;
dfan = atan(double(detector.hx_ISO)/double(detector.SID));
% dfan = detector.hx_ISO/detector.SID;

[dfan1, dfan2] = idealdfan(double(mean(fanangles0, 2)), double(detector.mid_U));

% idealfan = dfan.*((1:Npixel)'-detector.mid_U);
idealfan = dfan1.*((1:Npixel)'-(detector.mid_U));
% idealfan = single(idealfan);
% idealfan = linspace(dfan2.*(1-detector.mid_U), dfan2.*(Npixel-detector.mid_U), Npixel)';
t0 = mod(detector.mid_U, 1);
t1 = ceil(sin(min(idealfan))*detector.SID/delta_d + t0);
t2 = floor(sin(max(idealfan))*detector.SID/delta_d + t0);
idealphi = asin(((t1:t2)'-t0).*(double(delta_d)/double(detector.SID)));
Nphi = t2-t1+1;
midU_phi = -t1+1+t0;

Zdet = mean(reshape(detector.position(:,3), Npixel, Nslice));
Yshift = repelem(-Zdet(:).*tan(gantrytilt), Npixel, 1);
Yshift1 = -Zdet.*tan(gantrytilt)./detector.SDD;

delta_z = detector.hz_ISO;
Zgrid0 = single(1:Nslice);
Zsec = cos(gantrytilt)./cos(idealphi);
Zsec(Zsec>1) = 1;
Zgrid = Zsec*(-(Nslice-1)/2 : (Nslice-1)/2) + (Nslice+1)/2;
Xgrid0 = 1:Nphi;
Xgrid = repmat(Xgrid0(:), 1, Nslice);



faninterpkern = zeros(Npixel, Nslice, 'single');
for islice = 1:Nslice
    faninterpkern(:, islice) = interp1(fanangles0(:,islice), 1:Npixel, idealfan, 'linear', 'extrap');
end

viewangle = dataflow.rawhead.viewangle;

dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nview);


data0 = zeros(Npixel, Nslice, Nview);
data1 = zeros(Nphi, Nslice, Nview);
data2 = zeros(Nphi, Nslice, Nview);
data3 = zeros(Nphi, Nslice, Nview);
tic;
for iview = 1:Nview    
    phi1 = inverseyshift(idealphi, viewangle(iview), Yshift1);
    phi1 = reshape(phi1, Nphi, Nslice);
%     phi1_index = phi1./dfan + detector.mid_U;
    for islice = 1:Nslice
%         data_tmp = interp1(dataflow.rawdata(:, islice, iview), faninterpkern(:, islice), 'linear', 'extrap');
%         data2(:, islice, iview) = interp1(data_tmp, phi1_index(:, islice), 'linear', 0);
        
        data1(:, islice, iview) = interp1(fanangles0(:, islice), dataflow.rawdata(:, islice, iview), phi1(:, islice), 'linear', 0);
    end

    % z interp
    % ignore the components of Yshift;
%     data2(:, :, iview) = interp2(data2(:, :, iview), Zgrid, Xgrid);
    data1(:, :, iview) = interp2(data1(:, :, iview), Zgrid, Xgrid);
end
toc;

idealphi_gpu = gpuArray(idealphi);
viewangle_gpu = gpuArray(viewangle);
Yshift_gpu = gpuArray(Yshift1);
midU_gpu = gpuArray(detector.mid_U);
sliceindex = gpuArray(single(1:Nslice));
Nslice_gpu = gpuArray(Nslice);
Npixel_gpu = gpuArray(Npixel);
Nphi_gpu = gpuArray(Nphi);
faninterpkern = gpuArray(faninterpkern);
fanangles0_gpu = gpuArray(fanangles0);
pixelindex = gpuArray(single(1:Npixel));
Zgrid_gpu = gpuArray(Zgrid);
Xgrid_gpu = gpuArray(Xgrid);
dfan_gpu = gpuArray(dfan1);
tic;
for iview = 1:Nview
    phi1 = idealphi_gpu + asin(Yshift_gpu.*sin(viewangle_gpu(iview)+idealphi_gpu));
    phi1 = phi1./dfan_gpu + midU_gpu;
    data_iview = gpuArray(dataflow.rawdata(:, :, iview));
    data_tmp = interp2(data_iview, repmat(sliceindex, Npixel_gpu, 1), faninterpkern, 'linear', 0);
    data2_iview = interp2(data_tmp, repmat(sliceindex, Nphi_gpu, 1), phi1, 'linear', 0);
    data2(:,:,iview) = gather(interp2(data2_iview, Zgrid_gpu, Xgrid_gpu));
end
toc;

vblock = 32;
data_blk = zeros(Npixel, Nslice, vblock, 'single', 'gpuArray');
viewindex = reshape(gpuArray(single(1:vblock)), 1, 1, []);
vblock_gpu = gpuArray(single(vblock));
tic;
for iblk = 1:Nview/vblock
    vindex_iblk = (1:vblock)+(iblk-1)*vblock;
    data_blk = gpuArray(dataflow.rawdata(:, :, vindex_iblk));
    for islice = 1:Nslice
        data_blk(:, islice, :) = reshape(interp1(pixelindex, reshape(data_blk(:, islice, :), Npixel_gpu, vblock_gpu), ...
            faninterpkern(:, islice), 'linear', 'extrap'), Npixel_gpu, 1, vblock_gpu);
    end
    
    phi3 = idealphi_gpu + reshape(asin(reshape(sin(viewangle_gpu(vindex_iblk)+idealphi_gpu), [], 1)*Yshift_gpu), Nphi_gpu, []);
    phi3 = phi3./dfan_gpu + midU_gpu;
    phi3 = permute(reshape(phi3, Nphi_gpu, vblock_gpu, Nslice_gpu), [1 3 2]);
    data2_blk = interp3(data_blk, repmat(sliceindex, Nphi_gpu, 1, vblock_gpu), phi3, ...
        repmat(viewindex, Nphi_gpu, Nslice_gpu, 1), 'linear', 0);
    data3(:,:,vindex_iblk) = gather(interp3(data2_blk, repmat(Zgrid_gpu, 1, 1, vblock_gpu), ...
        repmat(Xgrid_gpu, 1, 1, vblock_gpu), repmat(viewindex, Nphi_gpu, Nslice_gpu, 1), 'linear', 0));
end
toc;

% 2
pf2 = prmflow;
df2 = dataflow;

pf2.recon.fanangles = repmat(idealphi, 1, Nslice);
pf2.recon.focalangle = focalangle0;
pf2.recon.Npixel = Nphi;
pf2.system.detector.mid_U = midU_phi;
df2.rawdata = data2;
[pf2, ~] = reconnode_rebinprepare(pf2, status);
[df2, pf2, ~] = reconnode_Azirebin(df2, pf2, status);
% pf2.recon.Npixel = rebin.Nreb;
pf2.recon.Nviewprot = pf2.rebin.Nviewprot;
pf2.recon.Nview = pf2.rebin.Nviewprot*Nshot;
pf2.recon.midchannel = pf2.rebin.midchannel;
pf2.recon.delta_d = pf2.rebin.delta_d;
pf2.recon.delta_z = pf2.rebin.delta_z;
pf2.recon.SID = pf2.rebin.SID;

% % Filter
% status.nodename = 'Filter';
% [dataflow, prmflow, status] = reconnode_Filter(dataflow, prmflow, status);
% 
% % BP prepare
% status.nodename = 'Backprojection';
% [dataflow, prmflow, status] = reconnode_BPprepare(dataflow, prmflow, status);

function phi1 = inverseyshift(phi2, theta, Ys)

% a = tan(phi2).^2+1;
% b = -2.*(cos(theta).*tan(phi2) + sin(theta)).*Ys;
% c = (tan(phi2).*cos(theta) + sin(theta)).^2.*Ys.^2-tan(phi2).^2;
% 
% phi1 = asin((-b+sqrt(b.^2-a.*c.*4))./a./2);

phi1 = phi2 + asin(Ys.*sin(theta+phi2));

end

function [dfan1, dfan2] = idealdfan(fanangles, mid_U)

N = length(fanangles);
index = (1:N) - mid_U;

dfan1 = (index*fanangles(:)) / sum(index.^2);

f2 = fanangles(:) - index(:).*dfan1;
dfan2 = (index*f2) / sum(index.^2);

end
