function corrtableoutput(SYS, Data)
% output the rawdata and air calibration table
% [raw, aircorr]= rawdataoutput(SYS, Data);

Nw = SYS.source.Wnumber;

% corr names
corrnames = fieldnames(SYS.output.files);
for ii = 1:length(corrnames)
    corr = corrnames{ii};
    switch corr
        case 'air'
            corrversion = SYS.output.corrversion.air;
            corrdata = simuAircali(SYS, Data, 1, corrversion);
        case 'beamharden'
            corrversion = SYS.output.corrversion.beamharden;
            corrdata = simuBHcali(SYS, 4, corrversion);
        otherwise
            continue
    end
    % output
    for iw = 1:Nw
        % output air corr table
        corrfile = fullfile(SYS.output.path, [SYS.output.files.(corr){iw} '.corr']);
        cfgfile = cfgmatchrule(corrfile, SYS.path.IOstandard, corrversion);
        corrcfg = readcfgfile(cfgfile);
        packstruct(corrdata{iw}, corrcfg, corrfile);
    end
end

end