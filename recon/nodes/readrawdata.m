function [dataflow, prmflow, status] = readrawdata(reconcfg, dataflow, prmflow, status)
% recon node to read raw data in recon/cali pipe line
% [dataflow, prmflow, status] = readrawdata(status.reconcfg, dataflow, prmflow, status);
% or a quick way
% [dataflow, prmflow] = readrawdata(recon_cfgxml);
% the recon_cfgxml is a strucut read from recon xml configure file, inwhich
% recon_cfgxml.rawdata is the rawdata file
% reconcfg.IOstandard is the data format configure path (If you do not 
% know what it is, try to set reconcfg.IOstandard=[] .)

% load raw data
[~, ~, fileEXT] = fileparts(reconcfg.rawdata);
switch lower(fileEXT)
    case {'.raw', '.bin'}
        raw = loaddata(reconcfg.rawdata, reconcfg.IOstandard);
        % data flow
        [dataflow.rawhead, dataflow.rawdata] = raw2dataflow(raw);
    case '.mat'
        raw = load(reconcfg.rawdata);
        if isfield(raw, 'rawhead') && isfield(raw, 'rawdata')
            dataflow = raw;
        else
            tmpfield = fieldnames(raw);
            [dataflow.rawhead, dataflow.rawdata] = raw2dataflow(raw.(tmpfield{1}));
        end
    case '.pd'
        dataflow = CRIS2dataflow(reconcfg.rawdata);
    otherwise
        error('Unknown rawdata ext: %s', fileEXT);
end

% load offset
if isfield(reconcfg, 'offset')
    offset = loaddata(reconcfg.offset, reconcfg.IOstandard);
    % to data flow
    dataflow.offset.rawhead.Reading_Number = [offset.Reading_Number];
    dataflow.offset.rawhead.Integration_Time = [offset.Integration_Time];
    dataflow.offset.rawhead.Time_Stamp = [offset.Time_Stamp];
    dataflow.offset.rawdata = single([offset.Raw_Data]);
end

% other
if isfield(reconcfg, 'system')
    % views
    dataflow.rawhead.viewangle = (single(dataflow.rawhead.Angle_encoder) - reconcfg.system.angulationzero) ...
                                 ./reconcfg.system.angulationcode.*(pi*2);
end

% recon parameters
if isfield(reconcfg, 'protocol')
%     viewnumber = reconcfg.protocol.viewnumber;
    prmflow.recon.Nshot = reconcfg.protocol.shotnumber;
    prmflow.recon.Nviewprot = reconcfg.protocol.viewperrot;
    % for Axial
    prmflow.recon.Nview = prmflow.recon.Nviewprot * prmflow.recon.Nshot;
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end


function [rawhead, rawdata] = raw2dataflow(raw)
% raw to dataflow

rawhead.Angle_encoder = [raw.Angle_encoder];
rawhead.Reading_Number = [raw.Reading_Number];
rawhead.Integration_Time = [raw.Integration_Time];
% rawhead.Time_Stamp = [raw.Time_Stamp];
rawhead.mA = single([raw.mA]);
rawhead.KV = single([raw.KV]);
rawdata = single([raw.Raw_Data]);

end

