function reconxml = simuresultsoutput(SYS, Data)
% outpu the simualtion results (rawdata, calibration tables and reconxml)
% reconxml = simuresultsoutput(SYS, Data);

% to check is the output path exist
if ~exist(SYS.output.path, 'dir')
    mkdir(SYS.output.path);
end

% output the rawdata and air (no offset?)
rawdataoutput(SYS, Data);

if Data.iblock == 1
    % output calibration tables
    corrtableoutput(SYS, Data);

    % recon xml
    [~, reconxml] = reconxmloutput(SYS);
else
    reconxml = [];
end

end