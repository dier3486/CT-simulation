resp0=load('D:\matlab\ct\BCT16\BHtest\response_1219.mat');
samplekeV = resp0.samplekeV;

Npixel = 912;
Nslice = 24;
Nps = Npixel*Nslice;
v = [0 50 100 150];
Nv = length(v);
Rv = rand(Nps, Nv);

R = zeros(Nps, length(samplekeV));
for ii = 1:Nps
    R(ii, :) = spline(v, [0 Rv(ii, :) 0], samplekeV);
end

resp1 = resp0;
resp1.response = R.*0.04+resp0.response.*0.98;