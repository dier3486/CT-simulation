function simuresultsoutput(SYS, Data)
% simuresultsoutput(SYS, Data)

% output the rawdata and air (no offset?)
rawdataoutput(SYS, Data);

% output calibration tables
corrtableoutput(SYS, Data);

% recon xml
reconxmloutput(SYS);

end