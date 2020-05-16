path0 = 'E:\matlab\CTsimulation\physics\tube';

tubespect =loaddata(fullfile(path0, 'tube_spectrumdata_v1.0.corr'));
tubespect.main = reshape(tubespect.main, tubespect.Nsample, []);


Nsample = tubespect.Nsample;
samplekeV = tubespect.main(:, 1);
p0 = tubespect.main(:, 2:2:end);
KVtag = tubespect.KVtag;

p70_1 = zeros(Nsample, 1);
for ii = 1:Nsample
   p70_1(ii) = interp1( KVtag, p0(ii,:), 70, 'spline', 'extrap');
end
% I know
peakindex = [58 59 67];
p70_1(peakindex) = [];
samplekeV_1 = samplekeV;
samplekeV_1(peakindex) = [];

p70_2 = interp1(samplekeV_1, p70_1, samplekeV);

p70_3 = p70_2 - p70_2(71).*samplekeV./71;

p70_3(70:end) = 0;
p70_3(p70_3<0) = 0;

tubespect2 = tubespect;
tubespect2.KVtag = [70; KVtag];
tubespect2.KVnumber = 5;
tubespect2.main = [samplekeV p70_3 tubespect.main];

file2 = fullfile(path0, 'tube_spectrumdata2_v1.0.corr');
bincfg = readcfgfile(cfgmatchrule(file2));
packstruct(tubespect2, bincfg, file2);