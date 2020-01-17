function [dataflow, prmflow, status] = reconnode_beamhardencorr(dataflow, prmflow, status)
% recon node, beamharden correction
% [dataflow, prmflow, status] = reconnode_beamhardencorr(dataflow, prmflow, status);

% parameters to use in prmflow
Nview = prmflow.recon.Nview;
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;

%--- test ---%
% cross talk (before bh)
corr_crs = load('E:\data\rawdata\bhtest\pcrs_2.mat');
p_crs = reshape(corr_crs.p_crs, Npixel, Nslice);
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nview);
dataflow.rawdata = 2.^(-dataflow.rawdata);
dataflow.rawdata(2:end-1,:,:) = dataflow.rawdata(2:end-1,:,:)+(dataflow.rawdata(3:end,:,:)-dataflow.rawdata(1:end-2,:,:)).*p_crs(2:end-1,:);
dataflow.rawdata(dataflow.rawdata<0) = 1;
dataflow.rawdata = -log2(dataflow.rawdata);
dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, []);
%--- test over ---%

% calibration table
bhcorr = prmflow.corrtable.(status.nodename);
bhorder = bhcorr.order;
bhpoly = reshape(bhcorr.main, [], bhorder);

% beam harden polynomial
dataflow.rawdata = reshape(dataflow.rawdata, [], Nview);
dataflow.rawdata = iterpolyval(bhpoly, dataflow.rawdata);

%--- test ---%
% nonlinear1
corr_nl = load('E:\data\rawdata\bhtest\pnl_2.mat');
p_nl = corr_nl.p_nl;
dataflow.rawdata = reshape(dataflow.rawdata, size(p_nl,1), []);
dataflow.rawdata = iterpolyval(p_nl, dataflow.rawdata);

% nonlinear2
corr_nl = load('E:\data\rawdata\bhtest\pnl2_1a.mat');
p2_nl = corr_nl.p2_nl;
dataflow.rawdata = reshape(dataflow.rawdata, size(p2_nl,1), []);
dataflow.rawdata = iterpolyval(p2_nl, dataflow.rawdata);
%--- test over ---%

% %--- test ---%
% % cross talk (after bh)
% corr_crs = load('E:\data\rawdata\bhtest\cross_1.mat');
% p_crs = reshape(corr_crs.p, Npixel, Nslice);
% p_scale = bhpoly(:,end).*p_nl(:,end);
% dataflow.rawdata = dataflow.rawdata./p_scale;
% dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nview);
% dataflow.rawdata = 2.^(-dataflow.rawdata);
% dataflow.rawdata(2:end-1,:,:) = dataflow.rawdata(2:end-1,:,:)+(dataflow.rawdata(3:end,:,:)-dataflow.rawdata(1:end-2,:,:)).*p_crs(2:end-1,:);
% dataflow.rawdata(dataflow.rawdata<0) = 1;
% dataflow.rawdata = -log2(dataflow.rawdata);
% dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, []).*p_scale;
% %--- test over ---%

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end