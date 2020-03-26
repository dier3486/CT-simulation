function [dataflow, prmflow, status] = reconnode_beamhardencali(dataflow, prmflow, status)
% recon node, beamharden calibration
% [dataflow, prmflow, status] = reconnode_beamhardencali(dataflow, prmflow, status);

% parameters set in pipe
bhcaliprm = prmflow.pipe.(status.nodename);
if isfield(bhcaliprm, 'polyorder')
    bhpolyorder = bhcaliprm.polyorder;
else
    % default polyorder is 3
    bhpolyorder = 3;
end
% format version of calibration table
if isfield(bhcaliprm, 'corrversion')
    corrversion = bhcaliprm.corrversion;
else
    corrversion = 'v1.0';
end
% to plot 
if isfield(bhcaliprm, 'toplot')
    toplot = bhcaliprm.toplot;
else
    toplot = false;
end

% parameters to use in prmflow
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
focalpos = prmflow.system.focalposition;
detpos = reshape(prmflow.system.detector.position, Npixel, Nslice, 3);
   
% material of the bowtie to fix
mu_1 = prmflow.SYS.collimation.bowtie{1}.material.mu_total(:);
samplekeV = prmflow.SYS.world.samplekeV(:);

% I know the simulation results are
Pempty = dataflow.simudata_bk1.Pair{1};
Pbowtie = dataflow.simudata_bk2.Pair{1};

% the simulated effective empty bowtie
Dempty = log(Pempty*samplekeV);
% the experiment effective bowtie thickness
Dexp = (dataflow.rawdata_bk2 - dataflow.rawdata_bk1).*log(2);
% try to fit the thickness fix 'dfit' to satisfy Dbowtie(dfit)-Dempty = Dexp.
dfit = zeros(Npixel, Nslice);
for ipixel = 1:Npixel*Nslice
    Pbow_ip = Pbowtie(ipixel, :);
    Dexp_ip = Dexp(ipixel)-Dempty(ipixel);
    if isfinite(Dexp_ip)
        dfit(ipixel) = fzero(@(x) -log(Pbow_ip*(exp(-x.*mu_1).*samplekeV))-Dexp_ip, 0);
    end
end
% smooth
for islice = 1:Nslice
    xx=detpos(:,islice,1)-focalpos(1);
    yy=detpos(:,islice,2)-focalpos(2);
    XYangle = atan2(yy, xx) - pi/2;
    dfit(:, islice) = smooth(XYangle, dfit(:,islice), 0.1, 'loess');
end

% plot
if toplot
    figure;
    plot(dfit);
    drawnow;
end

% effective filter's material
bhcaliprm = materialconfigure(bhcaliprm, samplekeV, prmflow.SYS.world.elementsdata, prmflow.SYS.world.materialdata);
if isfield(bhcaliprm, 'material')
    filtermaterial = bhcaliprm.material;
else
    % use the material of bowtie
    filtermaterial = prmflow.SYS.collimation.bowtie{1}.material;
end

% add effective filter
Nfilt = length(prmflow.SYS.collimation.filter);
prmflow.SYS.collimation.filter{Nfilt+1} = struct();
prmflow.SYS.collimation.filter{Nfilt+1}.effect = true;
prmflow.SYS.collimation.filter{Nfilt+1}.thickness = dfit(:);
prmflow.SYS.collimation.filter{Nfilt+1}.material = filtermaterial;
% merge detector slices
prmflow.SYS.detector = mergedetector(prmflow.SYS.detector);
% I know the prmflow.SYS is updated

% simu BH cali
BHcorr = simuBHcali(prmflow.SYS, bhpolyorder, corrversion);
beamhardencorr = BHcorr{1};

% air rate
% I know
airrate = reshape(Dexp./log(2), Npixel, Nslice);
% to smooth
for islice = 1:Nslice
    airrate(:, islice) = smooth(airrate(:, islice), 0.01);
end
beamhardencorr.airrate = airrate(:);

% to return
dataflow.beamhardencorr = beamhardencorr;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end