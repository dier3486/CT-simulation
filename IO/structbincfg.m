function bincfg = structbincfg(S, bincfg)
% to return the bin configure for a structure
%   bincfg = structbincfg(S);
% It is a tool in helping user to setup a configure file of packing a structure to a binary file.
% Suppose we have a data structure S needs to be packed to a binary file, e.g. a calibration table(.corr) or a rawdata(.raw).
% We need to set up a .xml configure file to define the data format, you can find some samples in ~\IO\standard\.
% And to edit the configure file could be tedious and error-prone. If that this function will be helpful,
%   1. setup a sample S in workspace, make use the size and class of all the fields are OK;
%   2. >> bincfg = structbincfg(S);
%   3. >> root = struct(); root.datanameofS = bincfg;
%   4. >> struct2xml(root, xmlfilename);
%   5. check and edit the xmlfilename;
%   6. use this xmlfilename to pack and sparse your data.
% Have a good use.

if nargin<2
    bincfg = struct();
    bincfg.offset = 0;
    bincfg.class = 'Struct';
    bincfg.size = 0;
end

if ~isstruct(S)
    return
end

% number (length)
bincfg.number = length(S);

% initial offset
offset_ini = 0;
% current offset
offset_cur = offset_ini;

% fields
sfields = fieldnames(S);
for ifield = 1: length(sfields)
    % field information
    field_ii = sfields{ifield};
    s_ii = S(1).(field_ii);
    class_ii = class(s_ii);
    
    % to configure
    bincfg.(field_ii).offset = offset_cur;
    bincfg.(field_ii).class = class_ii;
    % to avoid xml eval bug, set the first char in upper case
    bincfg.(field_ii).class(1) = upper(bincfg.(field_ii).class(1));
    switch class_ii
        case 'struct'
            bincfg.(field_ii).size = 0;
            % recurse
            bincfg.(field_ii) = structbincfg(s_ii, bincfg.(field_ii));
            offset_cur = offset_cur + bincfg.(field_ii).size * bincfg.(field_ii).number;
        case 'cell'
            % cell is not in supports
            bincfg.(field_ii).number = 0;
            bincfg.(field_ii).size = 0;
            offset_cur = offset_cur + 0;
        otherwise
            bincfg.(field_ii).number = length(s_ii(:));
            bincfg.(field_ii).size =  classsize(class_ii);
            offset_cur = offset_cur + bincfg.(field_ii).size * bincfg.(field_ii).number;
    end
end
bincfg.size = offset_cur - offset_ini;
    
return
