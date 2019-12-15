function [dataflow, prmflow, status] = readrawdata(reconcfg, dataflow, prmflow, status)
% read raw data in recon/cali pipe line
% [dataflow, prmflow, status] = readrawdata(status.reconcfg, dataflow, prmflow, status);
% or a quick way
% [dataflow, prmflow] = readrawdata(recon_cfgxml);
% the recon_cfgxml is a strucut read from recon xml configure file, inwhich
% recon_cfgxml.rawdata is the rawdata path
% reconcfg.IOstandard is the data format configure path (If you do not 
% know what it is, try to set reconcfg.IOstandard=[] .)

% load raw data
if iscell(reconcfg.rawdata)
    raw = loaddata(reconcfg.rawdata{status.series_index}, reconcfg.IOstandard);
else
    raw = loaddata(reconcfg.rawdata, reconcfg.IOstandard);
end

% data flow
dataflow.rawhead.Angle_encoder = [raw.Angle_encoder];
dataflow.rawhead.Reading_Number = [raw.Reading_Number];
dataflow.rawhead.Integration_Time = [raw.Integration_Time];
dataflow.rawhead.Time_Stamp = [raw.Time_Stamp];
dataflow.rawhead.mA = single([raw.mA]);
dataflow.rawhead.KV = [raw.KV];
dataflow.rawdata = single([raw.Raw_Data]);
% views
dataflow.rawhead.viewangle = (single(dataflow.rawhead.Angle_encoder) - reconcfg.system.angulationzero)./reconcfg.system.angulationcode.*(pi*2);

% Nview
prmflow.recon.Nview = reconcfg.protocol.viewnumber;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

return


