function raw = simurawdatastructV2(SYS, Data, rawdataversion, iw)
% v1
% data cell to structure

raw = struct();
Nview = length(Data.viewangle(:));
startreading = Data.startreading;
raw(Nview) = struct();

% major_version
[raw(:).major_version] = deal(rawdataversion(1));
% minor_version
[raw(:).minor_version] = deal(rawdataversion(2));
% header_size
[raw(:).header_size] = deal(160);
% img_width
[raw(:).img_width] = deal(SYS.detector.Npixel);
% img_height
[raw(:).img_height] = deal(SYS.detector.Nslice);
% bytes_per_pixel
if SYS.datacollector.islog2
    % 16 bit
    [raw(:).bytes_per_pixel] = deal(2);
else
    % 24 bit
    [raw(:).bytes_per_pixel] = deal(3);
end
% pixel_byte_order
[raw(:).pixel_byte_order] = deal(0);
% low_high_order
if isfield(SYS.detector, 'ZebraOrder')
    [raw(:).low_high_order] = deal(SYS.detector.ZebraOrder);
else
    [raw(:).low_high_order] = deal(0);
end
% sample_num (readingnumber)
readingnumber = num2cell((0:Nview-1)+startreading, 1);
[raw(:).sample_num] = readingnumber{:};
% integration_time
integrationtime = round(SYS.datacollector.integrationtime*1000);
[raw(:).integration_time] = deal(integrationtime);
% position (TableEncoder)
movercode = Data.couch(:,3)'./SYS.datacollector.moverlength;
movercode = round(movercode.*SYS.datacollector.moveruppersample);
movercode = num2cell(int64(movercode), 1);
[raw(:).position] = movercode{:};
% corrected_position
[raw(:).corrected_position] = movercode{:};
% correct_position_state
[raw(:).correct_position_state] = deal(1);
% running_state (TableGear)
if isfield(Data, 'couchgear')
    couchgear = num2cell(int8(Data.couchgear), 1);
    [raw(:).running_state] = couchgear{:};
else
    [raw(:).running_state] = deal(0);
end
% angle (angulation)
angcode = SYS.datacollector.angulationcode;
angzero = SYS.datacollector.angulationzero;
angleencoder0 = mod(round(Data.viewangle./(pi*2/angcode))+angzero, angcode);
angleencoder = num2cell(angleencoder0, 1);
[raw(:).angle] = angleencoder{:};
% angle_count_per_circle
[raw(:).angle_count_per_circle] = deal(angcode);
% view_count_per_circle
[raw(:).view_count_per_circle] = deal(SYS.protocol.viewperrot);
% is_xray_on
[raw(:).is_xray_on] = deal(SYS.source.mA{iw}>0);
% xray_index (focalspot)
Nfocal = SYS.source.focalnumber;
focalspots = SYS.protocol.focalspot(:)' - 1;  % fix the start from 0
focalspots = num2cell(repmat(focalspots, 1, Nview/Nfocal));
[raw(:).xray_index] = focalspots{:};
% xray_voltage (KV)
[raw(:).xray_voltage] = deal(SYS.source.KV{iw}*1000);
% xray_current (mA)
[raw(:).xray_current] = deal(SYS.source.mA{iw}*1000);
% timestamp
TimeStamp = zeros(1, Nview, 'uint64');
stamp0 = uint64(0);
stampperview = SYS.protocol.rotationspeed*1e6 / SYS.protocol.viewperrot;
stamp_rot0 = stamp0;
for ii = 1 : ceil(SYS.protocol.rotationnumber)
    N0 = (ii-1)*SYS.protocol.viewperrot;
    Nv = min(SYS.protocol.viewperrot, Nview - N0);
    TimeStamp(N0+(1:Nv)) = uint64(round((0:Nv-1).*stampperview)) + stamp_rot0;
    stamp_rot0 = TimeStamp(N0+Nv);
end
TimeStamp = num2cell(TimeStamp, 1);
[raw(:).timestamp] = TimeStamp{:};
% timestamp_when_recv_pos
[raw(:).timestamp_when_recv_pos] = TimeStamp{:};
% timestamp_when_recv_xray
[raw(:).timestamp_when_recv_xray] = TimeStamp{:};
% is_emergency_stop
[raw(:).is_emergency_stop] = deal(0);
% is_reconstruction
[raw(:).is_reconstruction] = deal(1);
% loss_view_count
[raw(:).loss_view_count] = deal(0);
% err_view_count
[raw(:).err_view_count] = deal(0);
% total_view_count
[raw(:).total_view_count] = readingnumber{:};
% board_count (Nmodule)
if isfield(SYS.detector, 'Nmodule')
    Nmodule = SYS.detector.Nmodule;
else
    Nmodule = ceil(SYS.detector.Npixel / 16);
end
[raw(:).board_count] = deal(Nmodule);
% one_board_info_byte_size
[raw(:).one_board_info_byte_size] = deal(9);    % to be fixed
% das_data_state
[raw(:).das_data_state] = deal(typecast(uint8([Nmodule 0]), 'uint16')); 
% view_crc32 (not support yet)
[raw(:).view_crc32] = deal(0); 
% zero_signal_count (not support yet, we might to have a rotationindex in Data.)
[raw(:).zero_signal_count] = deal(0);
% sample_signal_count (ShotNumber)
shotnumber = num2cell(Data.shotindex, 1);
[raw(:).sample_signal_count] = shotnumber{:};

% (one) board_info
board_info(Nmodule) = struct();
% board_info.board_num
board_num = num2cell(1:Nmodule);
[board_info(:).board_num] = board_num{:};
% board_info.board_is_recved
[board_info(:).board_is_recved] = deal(1);
% board_info.board_state
[board_info(:).board_state] = deal(0);
% board_info.low_board_temp (not support yet)
[board_info(:).low_board_temp] = deal(0);
% board_info.high_board_temp (not support yet)
[board_info(:).high_board_temp] = deal(0);
% board_info
[raw(:).board_info] = deal(board_info);

% raw data
RawData = num2cell(Data.P{iw}, 1);
[raw(:).RawData] = RawData{:};
    
end