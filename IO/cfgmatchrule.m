function [cfgfile, versionflag] = cfgmatchrule(targetfile, cfgpath, versionflag)
% return the bin format configure file 
% [cfgfile, versionflag] = cfgmatchrule(targetfile, cfgpath, versionflag);
%  or cfgfile = cfgmatchrule(targetfile, cfgpath);

if nargin<2
    cfgpath = '';
end

if nargin<3
    versionflag = '';
end

[filePATH, fileNAME, fileEXT] = fileparts(targetfile);
filekeys = regexp(fileNAME, '_', 'split');

cfgext = {'.json', '.xml', '.mat'};

if isempty(versionflag)
    if filekeys{end}(1) == 'v' && length(filekeys{end})>=2 && isempty(regexp(filekeys{end}(2:end), '[^\.\d]', 'once'))
        % while the fileNAME is *_v1.1.0, the versionflag is v1.1.0
        versionflag = filekeys{end};
    end
    % else the versionflag is ''.
end

switch fileEXT
    case '.raw'
        % rawdata
        if isempty(versionflag)
            % try to get data version from head
            cfgfile_base = fullfile(cfgpath, 'rawdatabase.xml');
            rawdata0 = loadbindata(targetfile, cfgfile_base);
            versionflag = ['v' num2str(rawdata0.PackageVersion(1)) '.' num2str(rawdata0.PackageVersion(2))];
        end
        rawfilename = ['rawdata_', versionflag];
        cfgfile = fullfile(cfgpath, rawfilename);
    case '.corr'
        % calibration table
        if isempty(versionflag)
            versionflag = 'v1.0';
        end
        corrfilename = [filekeys{1}, '_corr_', versionflag];
        cfgfile = fullfile(cfgpath, corrfilename);
    case '.bin'
        % other bin files
        if isempty(versionflag)
            versionflag = 'v1.0';
        end
        switch lower(filekeys{1})
            case 'dataflow'
                % dataflow
                binfilename = ['dataflow_', filekeys{2}, '_', versionflag];
            otherwise
                binfilename = [filekeys{1}, '_', versionflag];
        end
        cfgfile = fullfile(cfgpath, binfilename);
    otherwise
        cfgfile = targetfile;
end
% fill up ext
cfgfile = fillupext(cfgfile, cfgext);
% read configure file (skip)
% cfg = readcfgfile(cfgfile);

end

function filename_out = fillupext(filename_in, cfgext)
% fill up the EXT of a cfg file
filename_out = '';
for iext = 1:length(cfgext)
    if exist([filename_in, cfgext{iext}], 'file')
        filename_out = [filename_in, cfgext{iext}];
        % cfg = readcfgfile([cfgfile, cfgext{iext}]);
        break;
    end
end
end



