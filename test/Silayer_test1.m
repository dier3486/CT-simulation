T1 = loaddata('D:\matlab\CTsimulation\physics\tube\tube_spectrumdata_60-150_v1.0.corr');

T1.main = reshape(T1.main, T1.Nsample, []);
samplekeV = T1.main(10:end, 1);
sp120 = T1.main(10:end, 13*2);
M1 = materialdefine(loadmaterial('crystalSi'), samplekeV);

thick = 0:50;


