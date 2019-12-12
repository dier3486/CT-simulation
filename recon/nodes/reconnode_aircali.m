function [dataflow, prmflow, status] = reconnode_aircali(dataflow, prmflow, status)
% air calibration
% [dataflow, prmflow, status] = reconnode_aircali(dataflow, prmflow, status)

% parameters to use in prmflow
Npixel = prmflow.system.detector.Npixel;
Nslice = prmflow.system.detector.Nslice;
Nfocal = prmflow.system.Nfocal;
Nview = prmflow.Nview;

% parameters to use
caliprm = status.reconcfg.pipe.(status.nodename);
if isfield(caliprm, 'Nsection')
    Nsection = caliprm.Nsection;
else
    Nsection = 24;
end
if isfield(caliprm, 'refpixel')
    refpixel = caliprm.refpixel;
else
    refpixel = 16;
end
if isfield(caliprm, 'firstangle')
    firstangle = caliprm.firstangle;
else
    firstangle = 0;
end
if isfield(caliprm, 'corrversion')
    corrversion = caliprm.corrversion;
else
    corrversion = 'v1.0';
end

% paramters for aircorr
aircorr = caliprmforcorr(prmflow, corrversion);
aircorr.Nsection = Nsection;
aircorr.firstangle = firstangle;
aircorr.refpixel = refpixel;
aircorr.refnumber = 2;

% shift viewangle
viewangle = dataflow.rawhead.viewangle - firstangle;
% airmain and reference
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nview);
[aircorr.main, aircorr.reference] = aircalibration(dataflow.rawdata, viewangle, refpixel, Nsection, Nfocal);
% mainsize
aircorr.mainsize = length(aircorr.main(:));

end
