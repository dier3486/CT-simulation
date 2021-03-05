function datastruct = loaddata(datafile, cfgpath, varargin)
% fast data reading interface
% datastruct = loaddata(datafile, cfgpath)

if nargin<2
    cfgpath = '';
end

datafile= char(datafile);

if ~exist(datafile, 'file')
    error(['Not found file ' datafile]); 
end

[~, ~, fileEXT] = fileparts(datafile);

switch lower(fileEXT)
    case '.mat'
        % mat file
        datastruct = load(datafile);
    case {'.corr', '.raw', '.bin'}
        % bin file (calibration table, rawdata, binarry data)
        datastruct = loadbindata(datafile, cfgmatchrule(datafile, cfgpath));
    case {'.xml', '.json'}
        % configure file
        datastruct = readcfgfile(datafile);
    case '.csv'
        % .csv data
        datastruct.data = csvread(datafile);
        % it is not a good idea in using this function to read .csv, 
        % and no options supported
    case '.txt'
        % ascii data
        datastruct.data = load(datafile, '-acsii');
    case {'.dcm', '.dicom'}
        % dicom info
        datastruct = dicominfo(datafile);
        % not dicomread
%     case '.pd'
%         % external CRIS .pd rawdata file
%         % TBC, plz use CRIS2dataflow.m to read .pd files
    case '.ct'
        % external CRIS .ct file 
        datastruct = CRISgetcorr(datafile, varargin{:});
    otherwise
        warning(['Uknown file ext ''' fileEXT ''' to read file.'])
        datastruct = struct();
end
