function [bincfg, offset_align] = structbincfg(S, bincfg, alignpack)
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

if nargin<2 || isempty(bincfg)
    bincfg = struct();
    bincfg.offset = 0;
    bincfg.class = 'Struct';
    bincfg.size = 0;
end
if nargin<3
    alignpack = inf;
elseif ~(alignpack > 0)
    alignpack = 1;
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

% alignment
offset_align = 1;
toalign = true;

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
            [bincfg.(field_ii), offset_align_ii] = structbincfg(s_ii, bincfg.(field_ii), alignpack);
            offset_align = max(offset_align, offset_align_ii);
            % alignment
            if toalign && bincfg.(field_ii).size > 0
                offset_ini = offset_cur - bincfg.(field_ii).offset;
                bincfg.offset = bincfg.(field_ii).offset;
                toalign = false;
            end
        case 'cell'
            % cell is not in supports
            bincfg.(field_ii).number = 0;
            bincfg.(field_ii).size = 0;
            continue;
        otherwise
            bincfg.(field_ii).number = length(s_ii(:));
            bincfg.(field_ii).size =  classsize(class_ii);
            if isfloat(s_ii) && ~isreal(s_ii)
                % complex (e.g. float2)
                bincfg.(field_ii).size = bincfg.(field_ii).size*2;
            end
            offset_align_ii = min(bincfg.(field_ii).size, alignpack);
            % alignment
            if toalign && bincfg.(field_ii).size > 0
                offset_cur = ceil(bincfg.offset / offset_align_ii) * offset_align_ii;
                offset_ini = offset_cur - bincfg.offset;
                bincfg.offset = offset_cur;
                toalign = false;
            end
            bincfg.(field_ii).offset = ceil(offset_cur / offset_align_ii) * offset_align_ii;
            offset_align = max(offset_align, offset_align_ii);
    end
    offset_cur = bincfg.(field_ii).offset + bincfg.(field_ii).size * bincfg.(field_ii).number;
end
bincfg.size = ceil((offset_cur - offset_ini) / offset_align) * offset_align;

return
