function r = equivalentthickness(A, B, KV, elementspath, tubecorr)
% to return the equivalent thickness rate of A to B, A = r*B

if nargin<3
    % default reference matrial
    B = 'metalAl';
end

if nargin<3
    % default KV
    KV = 120;
end

if nargin<4
    % default elements path
    mypath = fileparts(which(mfilename));
    elementspath = fullfile(mypath, '\matter\elements\');
end

if nargin<5
    % default tube corr
    tubecorr = fullfile(mypath, '\tube\tube_spectrumdata2_v1.0.corr');
end

tube = loaddata(tubecorr);
tube.main = reshape(tube.main, tube.Nsample, []);
iKV = find(tube.KVtag == KV, 1);
samplekeV = tube.main(5:end, iKV*2-1);
soucespect = tube.main(5:end, iKV*2);
samplekeV = samplekeV(soucespect>0);
soucespect = soucespect(soucespect>0);

material_A = materialdefine(loadmaterial(A), samplekeV, elementspath);
material_B = materialdefine(loadmaterial(B), samplekeV, elementspath);

r0 = double((material_A.mu_total'*(soucespect.*samplekeV))/(material_B.mu_total'*(soucespect.*samplekeV)));
r = fzero(@(x) (exp(-material_A.mu_total')-exp(-material_B.mu_total'.*x))*(soucespect.*samplekeV), r0);

end

