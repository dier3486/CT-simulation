function bonecorr = loadboneharden(corrtable, FPchannelpos)
% reload Boneharden from bone-beamharden calibrationtble
% bonecorr = loadboneharden(corrtable, prmflow.recon.FPchannelpos);

% bone curve
bonecorr.bonecurve = reshape(corrtable.bonecurve, corrtable.bonecurvelength, []);

% and get these values
bonecorr.refrencemu = corrtable.refrencemu;
bonecorr.refrencebonemu = corrtable.refrencebonemu;
HCscale = 1000;
bonecorr.Dscale = 1/HCscale/corrtable.curvescale(1);
bonecorr.bonescale = 1/HCscale/corrtable.curvescale(2);
bonecorr.efffilter = interp1(corrtable.beamposition, corrtable.effbeamfilter, FPchannelpos, 'linear', 'extrap') ...
                       ./corrtable.curvescale(3);
% BBH table resampling
mubonmuw = single(corrtable.refrencebonemu/corrtable.refrencemu);
Nw = 400; Nb = 200; Nf = 20;
gridW = single(linspace(0, 0.5, Nw));
gridB = single(linspace(0, 2.0, Nb));
gridF = single(linspace(0.1, 1.0, Nf));
[gBB, gWW, gFF] = meshgrid(gridB, gridW, gridF);
bonecorr.order = reshape(corrtable.order, 1, []);
bonecorr.curvematrix = reshape(corrtable.curvematrix, bonecorr.order);
bonecorr.BBHmatrix = polyval3dm(bonecorr.curvematrix, gWW, gBB, gFF).*mubonmuw-1;
bonecorr.Wmesh = gWW;
bonecorr.Bmesh = gBB;
bonecorr.Fmesh = gFF;

end