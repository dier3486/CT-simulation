% READ ME, how to get a material

% Fist: addpath(genpath(<rootpath>));
% we wana 93WNiFe, play this:
samplekeV = 10:140;
material_93WNiFe = materialdefine(loadmaterial('93WNiFe'), samplekeV);
% DONE!