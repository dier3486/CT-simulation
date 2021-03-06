function cfgfile = cfgmatchrule(targetfile, cfgpath, versionflag)
% return the bin format configure file 

if nargin<2
    cfgpath = '';
end

if nargin<3
    versionflag = [];
end

[filePATH, fileNAME, fileEXT] = fileparts(targetfile);
filekeys = eval( ['{''', strrep(fileNAME, '_', ''','''), '''}'] );

cfgext = {'.xml', '.json', '.mat'};

switch fileEXT
    case '.raw'
        % rawdata
        if ~isempty(versionflag)
            rawfilename = ['rawdata_', versionflag];
        elseif filekeys{end}(1) == 'v' && ~isempty(str2double(filekeys{end}(2:end)))
            rawfilename = ['rawdata_', filekeys{end}];
        else
            % try to get data version from head
            cfgfile_base = fullfile(cfgpath, 'rawdatabase.xml');
            rawdata0 = loadbindata(targetfile, cfgfile_base);
            versionflag = ['v' num2str(rawdata0.Package_Version(1)) '.' num2str(rawdata0.Package_Version(2))];
            rawfilename = ['rawdata_', versionflag];
        end
        cfgfile = fullfile(cfgpath, rawfilename);
    case '.corr'
        % calibration table
        if ~isempty(versionflag)
            corrfilename = [filekeys{1}, '_corr_', versionflag];
        elseif filekeys{end}(1) == 'v' && ~isempty(str2double(filekeys{end}(2:end)))
            corrfilename = [filekeys{1}, '_corr_', filekeys{end}];
        else
            % default version is 1.0
            corrfilename = [filekeys{1}, '_corr_v1.0'];
        end
        cfgfile = fullfile(cfgpath, corrfilename);
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



