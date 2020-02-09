function [dataflow, prmflow, status] = reconnode_readrawdata(dataflow, prmflow, status)
% recon node to read raw data in recon/cali pipe line
% [dataflow, prmflow, status] = reconnode_readrawdata(dataflow, prmflow, status);
% NOTE: no quick start, if you only want to read a rawdata plz call loadrawdata.m
% this function is a recon pipe line node, but not an I/O function.
% 

if ~exist(prmflow.rawdata, 'file')
    status.jobdone = false;
    status.errorcode = 1;
    status.errormsg = '[readrawdata] rawdata file not exist.';
    return;
end

% load raw data
dataflow = structmerge(loadrawdata(prmflow.rawdata, prmflow.IOstandard), dataflow);

% load offset
if isfield(prmflow, 'offset') && ~isempty(prmflow.offset)
    dataflow.offset = loadrawdata(prmflow.offset, prmflow.IOstandard);
end

% other
if isfield(prmflow, 'system')
    % views
    dataflow.rawhead.viewangle = (single(dataflow.rawhead.Angle_encoder) - prmflow.system.angulationzero) ...
                                 ./prmflow.system.angulationcode.*(pi*2);
end

% recon parameters
if isfield(prmflow, 'protocol')
%     viewnumber = reconcfg.protocol.viewnumber;
    prmflow.recon.Nshot = prmflow.protocol.shotnumber;
    prmflow.recon.Nviewprot = prmflow.protocol.viewperrot;
    % for Axial
    prmflow.recon.Nview = prmflow.recon.Nviewprot * prmflow.recon.Nshot;
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end


function dataflow = loadrawdata(filename, IOstandard)
% load rawdata (or offset) from the file filename

[~, ~, fileEXT] = fileparts(filename);
switch lower(fileEXT)
    case {'.raw', '.bin'}
        raw = loaddata(filename, IOstandard);
        % data flow
        [dataflow.rawhead, dataflow.rawdata] = raw2dataflow(raw);
    case '.mat'
        raw = load(filename);
        if isfield(raw, 'rawhead') && isfield(raw, 'rawdata')
            dataflow = raw;
        else
            tmpfield = fieldnames(raw);
            [dataflow.rawhead, dataflow.rawdata] = raw2dataflow(raw.(tmpfield{1}));
        end
    case '.pd'
        % external IO of .pd
        dataflow = CRIS2dataflow(filename);
    otherwise
        error('Unknown rawdata ext: %s', fileEXT);
end

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

