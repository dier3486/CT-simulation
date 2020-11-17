function [dataflow, protocol] = CRIS2dataflow(pdfile, shotindex, shotnumber, protocol, CRISpath)
% read .pd rawdata
% [dataflow, protocol] = CRIS2dataflow(pdfile, shotIndex, shotnumber CRISpath);
% or [dataflow, protocol] = CRIS2dataflow(pdfile);
% INPUT:
%   pdfile          filename 
%   shotindex       start which shot
%   shotnumber      to read how many shots
%   CRISpath        external IO platform

if nargin<2
    % from 1st shot
    shotindex = 1;
end
if nargin<3
    % read all shots
    shotnumber = inf;
end
if nargin<4
    protocol = [];
end
if nargin<5 || isempty(CRISpath)
    % default CRIS path
    CRISpath = 'E:\matlab\CRIS_v1.0';
    % I know it is
end

% path
currpath = pwd;
cd(CRISpath);
addpath(genpath('.'));

pdfile = char(pdfile);
% data xml
if isempty(protocol)
    protocol = CRISgetrawxml(pdfile);
end

% read data
raw = CRISreadraw(pdfile, shotindex, shotnumber, protocol);

% path go back
cd(currpath);

% % I know
Nview = size(raw.data, 2);
% Nps = raw.raw_size(1)*raw.raw_size(2);
% Nview = raw.raw_size(3);

% to dataflow
dataflow = struct();
% rawdata
% dataflow.rawdata = reshape(raw.data, [], Nview);
dataflow.rawdata = raw.data;
% rawdata head
dataflow.rawhead = struct();
% I know KV*5 and mA*4
dataflow.rawhead.mA = raw.header.mA(:)'.*4.0;
dataflow.rawhead.KV = raw.header.kV(:)'.*5.0;

% viewangle
dataflow.rawhead.viewangle = (raw.header.viewAngle(:)' - 0).*(pi/180);

% viewcodes (no viewcodes in raw now)
N0 = 69120;
dataflow.rawhead.Angle_encoder = round(raw.header.viewAngle(:)'.*(N0/360));

% no reading number yet
% dataflow.rawhead.Reading_Number = raw.active.header.AdcViewNum(:)';

% the integration time should be
% dataflow.rawhead.Integration_Time = round(raw.header.Trigger_Interval(:)'*1e3/8);
% but
switch raw.header.Data_Format(1)
    case 1
        % 16 bit
        dataflow.rawhead.Integration_Time = ones(1, Nview);
        raw.offset(:) = 0;
    otherwise
        T0 = round(1.0e9/protocol.AcquisitionParameter.trigFrequency/8);
        dataflow.rawhead.Integration_Time = repmat(T0, 1, Nview);
end

% offset
dataflow.offset = struct();
dataflow.offset.rawdata = raw.offset;
dataflow.offset.rawhead = struct();
% dataflow.offset.rawhead.Reading_Number = raw.offset.header.AdcViewNum(:)';
% dataflow.offset.rawhead.Integration_Time = repmat(T0, 1, Nvoff);

end


function raw = CRISreadraw(Rawfilepath, shotIndex, shotNumber, protocol)
% call RawDataManager to read rawdata

struct_global_parameter = struct();
struct_global_parameter.RawData.RawDataPath = Rawfilepath;
struct_global_parameter.ReconParameters.ECGFileName = '';

rawInstance = RawDataManager.GetInstance();
rawInstance.Init(struct_global_parameter);

shotNumber = min(shotNumber, protocol.AcquisitionParameter.ShotNumber-shotIndex+1);
viewNumber_pershot = protocol.AcquisitionParameter.numView;
% viewNumber = viewNumber_pershot*shotNumber;

% ini
raw.data = [];
raw.offset = [];
raw.header = struct();

for ishot = 1:shotNumber
    % view offset
    nViewOffset = (shotIndex + ishot - 2)*viewNumber_pershot + 1;
    % read raw
    [ret, raw_ii] = rawInstance.read_PDType_raw(viewNumber_pershot, nViewOffset);
    if ret
        error('Error in reading .pd file %s.', Rawfilepath) ;
    end
    raw.data = [raw.data  reshape(raw_ii.data, [], viewNumber_pershot)];
    raw.offset = reshape(raw_ii.offset, size(raw.data, 1), []);
    [ret, header_ii] = rawInstance.ReadRawDataHeaders(viewNumber_pershot, nViewOffset);
    if ret
        error('Error in reading data head from .pd file %s.', Rawfilepath);
    end
    raw.header = rawheadermerge(raw.header, header_ii);
end

% fix PX ref
switch raw.header.MachineType(:, 1)'
    case {'INSITUM_16s_338', 'INSITUM_338'}
        % PX
        Nslice = raw.header.CollimatedSliceNum(1);
        Nview_data = size(raw.data, 2);
        Nview_offset = size(raw.offset, 2);
        % skip the ref
        raw.data = reshape(raw.data, [], Nslice, Nview_data);
        raw.data = reshape(raw.data(17:end, :, :), [], Nview_data);
        raw.offset = reshape(raw.offset, [], Nslice, Nview_offset);
        raw.offset = reshape(raw.offset(17:end, :, :), [], Nview_offset);
    otherwise
        % do nothing
        1;
end

end


function head1 = rawheadermerge(head1, head2)

headfields = fieldnames(head2);
Nhf = size(headfields(:),1);
for ihf = 1:Nhf
    hfield_ii = headfields{ihf};
    if isfield(head1, hfield_ii)
        head1.(hfield_ii) = [head1.(hfield_ii) head2.(hfield_ii)'];
    else
        head1.(hfield_ii) = head2.(hfield_ii)';
    end
end

end