function [dataflow, prmflow, status] = reconnode_aircorr(dataflow, prmflow, status)
% recon node, air correction
% [dataflow, prmflow, status] = reconnode_aircorr(dataflow, prmflow, status);

% parameters to use in prmflow
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nview = prmflow.recon.Nview;
Nfocal = prmflow.system.Nfocal;

% calibration table
aircorr = prmflow.corrtable.(status.nodename);
% parameters in corr
Nsect = single(aircorr.Nsection);
refpixel = aircorr.refpixel;
Nref = aircorr.refnumber;

% angles of the air table
sectangle = (pi*2/Nsect);
% airangle = (-1:Nsect).*(pi*2/Nsect);
% airmain & airref
aircorr.main = reshape(aircorr.main, [], Nsect);
airmain = [aircorr.main aircorr.main(:,1)];
aircorr.reference = reshape(aircorr.reference, [], Nsect);
airref = [aircorr.reference aircorr.reference(:,1)];

% interp index and weight
retangle = mod(dataflow.rawhead.viewangle - aircorr.firstangle, pi*2);
intp_index = floor(retangle./sectangle)+1;
intp_alpha = retangle./sectangle - intp_index;

% rawref
rawref = airreference(dataflow.rawdata, refpixel, Npixel, Nslice);
% corr rawref
rawref = rawref - airref(:, intp_index).*(1-intp_alpha) - airref(:, intp_index+1).*intp_alpha;

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nview);

% corr rawdata with ref
for ifocal = 1:Nfocal
    viewindex = ifocal:Nfocal:Nview;
    refindex = (1:Nref) + Nref*(ifocal-1);
    dataflow.rawdata(:, viewindex) = dataflow.rawdata(:, viewindex) - mean(rawref(refindex, :), 1);
end

% corr rawdata with air
for ifocal = 1:Nfocal
    viewindex = ifocal:Nfocal:Nview;
    airindex = (1:Npixel*Nslice) + Npixel*Nslice*(ifocal-1);
    dataflow.rawdata(:, viewindex) = dataflow.rawdata(:, viewindex) - airmain(airindex, intp_index).*(1-intp_alpha(viewindex));
    dataflow.rawdata(:, viewindex) = dataflow.rawdata(:, viewindex) - airmain(airindex, intp_index+1).*intp_alpha(viewindex);
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end