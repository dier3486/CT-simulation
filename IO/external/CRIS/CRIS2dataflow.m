function [dataflow, protocol] = CRIS2dataflow(pdfile, CRISpath)
% read .pd rawdata
% [dataflow, protocol] = CRIS2dataflow(pdfile, CRISpath);
% INPUT:
%   pdfile          filename 
%   CRISpath        external IO platform

if nargin<2 || isempty(CRISpath)
    % default CRIS path
    CRISpath = 'E:\matlab\CRIS';
    % I know it is
end

% read data
currpath = pwd;
cd(CRISpath);
raw = readRaw(pdfile, 4);
cd(currpath);

% almost protocol
protocol = raw.header;

% to dataflow
dataflow = struct();
% rawdata
dataflow.rawdata = reshape(raw.active.chanData, raw.header.numView, [])';
% rawdata head
dataflow.rawhead = struct();
% I know KV*5 and mA*4
dataflow.rawhead.mA = raw.active.header.mA(:)'.*4.0;
dataflow.rawhead.KV = raw.active.header.kv(:)'.*5.0;
N0 = 69120;
viewcodes = (0:raw.header.numView-1).*(N0 / raw.header.trigFrequency * raw.header.gantrySpeed);
startcode = round(raw.header.startAngle*N0/360);
viewcodes = mod(viewcodes + startcode, N0);
dataflow.rawhead.Angle_encoder = viewcodes;
dataflow.rawhead.Reading_Number = raw.active.header.AdcViewNum(:)';
T0 = round(1.0e9/raw.header.trigFrequency/8);
dataflow.rawhead.Integration_Time = repmat(T0, 1, raw.header.numView);
% views
dataflow.rawhead.viewangle = (single(dataflow.rawhead.Angle_encoder) - 0)./N0.*(pi*2);

% offset
Nvoff = size(raw.offset.chanData, 1);
dataflow.offset = struct();
dataflow.offset.rawdata = reshape(raw.offset.chanData, Nvoff, [])';
dataflow.offset.rawhead = struct();
dataflow.offset.rawhead.Reading_Number = raw.offset.header.AdcViewNum(:)';
dataflow.offset.rawhead.Integration_Time = repmat(T0, 1, Nvoff);

end