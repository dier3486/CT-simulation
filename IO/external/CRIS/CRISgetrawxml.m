function rawxml = CRISgetrawxml(datafile)
% create or load .pd data's .xml

% .xml filename
[filepath, filename, fileext] = fileparts(datafile);
filexml = fullfile(filepath, [filename, '.xml']);

if exist(filexml, 'file')
    % load xml
%     tmp = myxml2struct(filexml);
    rawxml = Global_Parameter(filexml);
elseif strcmpi(fileext, '.pd')
    % to create .xml file
    struct_global_parameter = struct();
    struct_global_parameter.RawData.RawDataPath = datafile;
    struct_global_parameter.ReconParameters.ECGFileName = '';
    rawInstance = RawDataManager.GetInstance();
    rawInstance.Init(struct_global_parameter);
    xmlname = rawInstance.GenXmlFromPD(datafile);
%     tmp = myxml2struct(xmlname);
    rawxml = Global_Parameter(xmlname);
else
    % not .pd, return empty 
    rawxml = [];
    % do nothing
end