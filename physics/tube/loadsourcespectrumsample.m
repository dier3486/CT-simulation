% READ ME, how to get the spectrum of a KVp

% Fist: addpath(genpath(<rootpath>));
% we wana 120KVp, play this:
KVp = 120;
samplekeV = 10:0.5:140;
tubedata = loaddata('tube_spectrumdata_60-160_v1.0.corr');
tubedata.main = reshape(tubedata.main, [], tubedata.KVnumber);
spectdata = reshape(tubedata.main(:, (tubedata.KVtag == KVp)), [], 2);
spectrum = interp1(spectdata(:,1), spectdata(:,2), samplekeV, 'linear', 0);
% or spectrum = spectrumresample(spectdata(:,1), spectdata(:,2), samplekeV);

% DONE!