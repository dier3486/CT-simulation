function aircorr = simuAircali(SYS, Dataflow, Nsection, corrversion)
% simulation of air calibration
% aircorr = simuAircali(SYS, Dataflow, Nsection)

if nargin < 3
    Nsection = 1;
end
if nargin < 4
    corrversion =  'v1.0';
end
Nw = SYS.source.Wnumber;

% paramters to put in table
corrprm = parameterforcorr(SYS, corrversion);

% % corr table baseline (need not)
% aircorr_basefile = [SYS.path.IOstandard, 'air_sample_v1.0.corr'];
% if exist(aircorr_basefile, 'file')
%     aircorr_base = loaddata(aircorr_basefile, SYS.path.IOstandard);
% else
%     % empty baseline
%     aircorr_base = struct();
% end

% initial
aircorr = cell(1, Nw);
aircorr(:) = {struct()};
% loop Nw
for iw = 1:Nw
    % v1.x 
    % values to put in struct
    aircorr{iw}.ID = corrprm.ID;
    aircorr{iw}.Npixel = corrprm.Npixel;
    aircorr{iw}.Nslice = corrprm.slicenumber;
    aircorr{iw}.startslice = corrprm.startslice;
    aircorr{iw}.endslice = corrprm.endslice;
    aircorr{iw}.mergescale = corrprm.mergescale;
    aircorr{iw}.focalspot = corrprm.focalspot;
    aircorr{iw}.focalsize = corrprm.focalsize;
    aircorr{iw}.KV = corrprm.KV{iw};
    aircorr{iw}.mA = corrprm.mA_air{iw};
    aircorr{iw}.bowtie = corrprm.bowtie;
    aircorr{iw}.rotationspeed = corrprm.rotationspeed;
    aircorr{iw}.focalnumber = corrprm.focalnumber;
    aircorr{iw}.refnumber = 2;
    refpixel = 16;
    aircorr{iw}.refpixel = refpixel;
    aircorr{iw}.Nsection = Nsection;
    aircorr{iw}.firstangle = 0;
    aircorr{iw}.mainsize = length(Dataflow.Pair{iw}(:))*Nsection;
    aircorr{iw}.referenceKVmA = repmat(-log2(aircorr{iw}.KV*aircorr{iw}.mA), corrprm.focalnumber, 1);
    aircorr{iw}.referrcut = repmat([0.01; 0.01], 1, corrprm.focalnumber);
    
    % reference
    airref = airreference(Dataflow.Pair{iw}, refpixel, corrprm.Nprange, corrprm.slicenumber);
    aircorr{iw}.reference = repmat(single(airref), 1, Nsection);
    % main
    aircorr{iw}.main = repmat(single(Dataflow.Pair{iw}(:)), 1, Nsection);
    % I know the Dataflow.Pair is the air data
    
    % NOTE: those codes shall be upgraded to aircalibration2.
end

end