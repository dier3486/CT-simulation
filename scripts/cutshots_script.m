% to cut a shot from multi shots axial

% data after rebin
load E:\matlab\Data\Lcone\test0307\data0.mat

selectshot = 2;

recon = prmflow.recon;
Nshot = recon.Nshot;
Nviewprot = recon.Nviewprot;
Nview = recon.Nview;

viewindex = (1:Nviewprot) + (selectshot-1)*Nviewprot;

% cut rawdata
dataflow.rawdata = reshape(dataflow.rawdata, [], Nview);
dataflow.rawdata = dataflow.rawdata(:, viewindex);
headfields = fieldnames(dataflow.rawhead);
for ii = 1:length(headfields)
    field_ii = headfields{ii};
    try
        dataflow.rawhead.(field_ii) = reshape(dataflow.rawhead.(field_ii), [], Nview);
    catch ME
        continue;
    end
    dataflow.rawhead.(field_ii) = dataflow.rawhead.(field_ii)(:, viewindex);
end

% protocol
prmflow.protocal.shotnumber = 1;
% recon
prmflow.recon.Nshot = 1;
prmflow.recon.Nview = Nviewprot;
prmflow.recon.startviewangle = prmflow.recon.startviewangle(selectshot);

