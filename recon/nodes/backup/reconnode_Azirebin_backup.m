function [dataflow, prmflow, status] = reconnode_Azirebin(dataflow, prmflow, status)
% recon node, Azi rebin for axial
% [dataflow, prmflow, status] = reconnode_Azirebin(dataflow, prmflow, status);
% fly-focal is not supported yet

% parameters to use in prmflow
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nshot = prmflow.recon.Nshot;
Nview = prmflow.recon.Nview;
Nviewprot = prmflow.recon.Nviewprot;
% from .rebin
Npixel_rebin = prmflow.rebin.Npixel;
Nviewprot_rebin = prmflow.rebin.Nviewprot;
% I know Nview=Nviewprot*Nshot, for axial

% rebin is prepared
rebin = prmflow.rebin;

% Azi rebin
dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nview);
for ishot = 1:Nshot
    A = zeros(Npixel_rebin*Nslice, Nviewprot_rebin, 'single');
    viewindex = (1:Nviewprot) + (ishot-1)*Nviewprot;
    A(rebin.vindex1_azi) = ...
        dataflow.rawdata(:, viewindex).*repmat(1-rebin.interalpha_azi, 1, Nviewprot);
    A(rebin.vindex2_azi) = A(rebin.vindex2_azi) + ...
        dataflow.rawdata(:, viewindex).*repmat(rebin.interalpha_azi, 1, Nviewprot);
    % start angle for first rebin view, replace the rawdata
    dataflow.rawdata(:, viewindex) = ...
        [A(:, (rebin.startvindex : Nviewprot)) A(:, (1 : rebin.startvindex-1))];
    
    % ref blocked
    % to find out the views whose references are blocked after rebin
    if isfield(dataflow.rawhead, 'refblock') && any(dataflow.rawhead.refblock(:))
        B = false(2, Nviewprot);
        % ref1
        v1_ref1=ceil(rebin.vindex1_azi(1,:)./(Npixel*Nslice));
        v2_ref1=ceil(rebin.vindex2_azi(1,:)./(Npixel*Nslice));
        B(1, v1_ref1) = dataflow.rawhead.refblock(1, viewindex);
        B(1, v2_ref1) = B(1, v2_ref1) | dataflow.rawhead.refblock(1, viewindex);
        % ref2
        v1_ref2=ceil(rebin.vindex1_azi(end,:)./(Npixel*Nslice));
        v2_ref2=ceil(rebin.vindex2_azi(end,:)./(Npixel*Nslice));
        B(2, v1_ref2) = dataflow.rawhead.refblock(2, viewindex);
        B(2, v2_ref2) = B(2, v2_ref2) | dataflow.rawhead.refblock(2, viewindex);
        % start angle
        dataflow.rawhead.refblock(:, viewindex) = ...
            [B(:, (rebin.startvindex : Nviewprot)) B(:, (1 : rebin.startvindex-1))];
    end
end
% viewangle
viewangle = reshape(dataflow.rawhead.viewangle, Nviewprot, Nshot);
startviewangle = viewangle(rebin.startvindex, :);
dataflow.rawhead.viewangle = [viewangle(rebin.startvindex : Nviewprot, :); viewangle(1 : rebin.startvindex-1, :)];
dataflow.rawhead.viewangle = dataflow.rawhead.viewangle(:)';

% prm to return
prmflow.recon.startviewangle = startviewangle;
prmflow.recon.delta_view = rebin.delta_view;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end