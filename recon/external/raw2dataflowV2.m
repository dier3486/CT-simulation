function [rawhead, rawdata, headprm] = raw2dataflowV2(raw, startview, viewnum)
% raw to dataflow structure format, v2.* for deployments

% raw shall be in size 1*Nview
raw = raw(:)';

% current shot(s)
if nargin > 2
    endview = min(startview + viewnum - 1, size(raw, 2));
    raw = raw(startview : endview);
elseif nargin > 1
    raw = raw(startview : end);
end
Nview = size(raw, 2);

% rawhead
% AngleEncoder := angle
rawhead.AngleEncoder = mod([raw.angle], raw(1).angle_count_per_circle);
% ReadingNumber := sample_num
rawhead.ReadingNumber = [raw.sample_num];
% IntegrationTime := integration_time/8, (/8 is the compatibility to orinigal 8ns unit)
rawhead.IntegrationTime = [raw.integration_time]./8;
% ShotNumber := sample_signal_count
rawhead.ShotNumber = [raw.sample_signal_count];
% TimeStamp := timestamp
% rawhead.TimeStamp = [raw.timestamp];
% mA := xray_current/1000
rawhead.mA = single([raw.xray_current])./1000;
% kV := xray_voltage/1000
rawhead.KV = single([raw.xray_voltage])./1000;
% TableEncoder := uint32(corrected_position)
TableEncoder = typecast([raw.corrected_position], 'uint32');
rawhead.TableEncoder = TableEncoder(1:2:end);
% TableGear := running_state
rawhead.TableGear = [raw.running_state];

% We may have a list to configure in reading those tags to rawhead, which could be depending on protocol.

% rawdata
rawdata = single([raw.RawData]);
% I know all the fields of rawhead and rawdata are in size n*Nview

if ~isempty(raw)
    headprm.PackageVersion = [raw(1).major_version raw(1).minor_version];
    headprm.SeriesNumber = 1;
    headprm.StartSlice = 1;
    headprm.EndSlice = raw(1).img_height;
    headprm.SliceMergescale = 1;
    headprm.SliceNumber = raw(1).img_height;
    headprm.RawDataSize = raw(1).img_width*raw(1).img_height;
else
    headprm = [];
end

end