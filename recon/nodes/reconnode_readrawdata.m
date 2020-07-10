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

if isfield(prmflow, 'protocol')
    shotnum = prmflow.protocol.shotnumber;
    if isfield(prmflow.protocol, 'startshot')
        startshot = prmflow.protocol.startshot;
    else
        startshot = 1;
    end
    viewpershot = prmflow.protocol.viewnumber;
else
    startshot = 1;
    shotnum = 1;
end

% load raw data
dataflow = structmerge(loadrawdata(prmflow.rawdata, prmflow.IOstandard, startshot, shotnum, viewpershot), dataflow, 0, 0);

% load offset
if isfield(prmflow, 'offset') && ~isempty(prmflow.offset)
    dataflow.offset = loadrawdata(prmflow.offset, prmflow.IOstandard, 1, 1, 10000);
end

% other
if isfield(prmflow, 'system')
    % viewangle (if not exist)
    if ~isfield(dataflow.rawhead, 'viewangle')
        dataflow.rawhead.viewangle = (single(dataflow.rawhead.Angle_encoder) - prmflow.system.angulationzero) ...
                                 ./prmflow.system.angulationcode.*(pi*2);
    end
end

% copy to recon parameters
% shot number
prmflow.recon.Nshot = shotnum;
% from protocol
if isfield(prmflow, 'protocol')
%     viewnumber = reconcfg.protocol.viewnumber;
    prmflow.recon.Nviewprot = prmflow.protocol.viewperrot;
    % viewnumber
    prmflow.recon.Nview = prmflow.protocol.viewnumber * shotnum;
    % I know the prmflow.protocol.viewnumber is the view number per shot for axial, and for helical only one shot once.
    % scan
    prmflow.recon.scan = lower(prmflow.protocol.scan);
    
    % explain focal spot
    focalspot_0x = focalspot20x(prmflow.protocol.focalspot);
    spots = fliplr(dec2bin(focalspot_0x)=='1');
    prmflow.recon.Nfocal = sum(spots);
    prmflow.recon.focalspot = find(spots);
    % NOTE: prmflow.protocol.focalspot is the name of the focalspot mode,
    %       prmflow.recon.focalspot is the index of the focalspot(s).
    
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end


function [dataflow, protocol] = loadrawdata(filename, IOstandard, startshot, shotnum, viewpershot)
% load rawdata (or offset) from the file filename

protocol = [];
[~, ~, fileEXT] = fileparts(filename);
switch lower(fileEXT)
    case {'.raw', '.bin'}
        % tmp code
        raw = loaddata(filename, IOstandard);
        % .raw should have a single code to read 
        
        % data flow
        [dataflow.rawhead, dataflow.rawdata] = raw2dataflow(raw, startshot, shotnum, viewpershot);
    case '.mat'
        % load mat
        raw = load(filename);
        % data flow
        if isfield(raw, 'rawhead') && isfield(raw, 'rawdata')
            dataflow = raw;
            % but we can not select the shot(s) after it has been merged in dataflow
        else
            tmpfield = fieldnames(raw);
            [dataflow.rawhead, dataflow.rawdata] = raw2dataflow(raw.(tmpfield{1}), startshot, shotnum, viewpershot);
        end
    case '.pd'
        % external IO of .pd
        [dataflow, protocol] = CRIS2dataflow(filename, startshot, shotnum);
    otherwise
        error('Unknown rawdata ext: %s', fileEXT);
end

end


function [rawhead, rawdata] = raw2dataflow(raw, startshot, shotnum, viewpershot)
% raw to dataflow

% current shot(s)
Nraw = size(raw(:),1);
viewstart = (startshot-1)*viewpershot+1;
viewend = min(viewstart-1+shotnum*viewpershot, Nraw);
raw = raw(viewstart:viewend);

rawhead.Angle_encoder = [raw.Angle_encoder];
rawhead.Reading_Number = [raw.Reading_Number];
rawhead.Integration_Time = [raw.Integration_Time];
% rawhead.Time_Stamp = [raw.Time_Stamp];
rawhead.mA = single([raw.mA]);
rawhead.KV = single([raw.KV]);
rawdata = single([raw.Raw_Data]);

end

