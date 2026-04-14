function [data, bincfg] = packstruct(S, bincfg, outputfile, wora)
% transform a struct to bin data and/or write in a file
% [data, bincfg] = packstruct(S, bincfg, outputfile, wora);
% or, packstruct(S, bincfg, outputfile);
% INPUT:
%   S               the structure to pack
%   bincfg          the configure of the data format
%   outputfile      output file name (path)
%   wora            'w', 'a', 'W' or 'A', method in writing the file, default is 'w'.
% OUTPUT:
%   data            the bin data to (or did) save in file
%   bincfg          the returned configure of the data format
%   write the outputfile
% NOTE: If you don't know how to get the bincfg, try this: bincfg = readcfgfile(cfgmatchrule(outputfile, cfgpath));
% where the cfgpath is the folder of the format configure files, ~/IO/standard/;
% But if the bincfg file not exists, you may need to set up it which is used to configre the data format. You may look up the
% function structbincfg.m for more information.

if nargin<2 || isempty(bincfg)
    bincfg = structbincfg(S);
else
    if regexp(class(bincfg), '^lib\.')
        % bincfg is a lib.struct
        bincfg = get(bincfg);
        bincfg = renameStructField2(bincfg, 'dataclass', 'class');
    end
    bincfg = clearbincfg(bincfg, S(1));
end
% to double
bincfg = everything2single(bincfg, 'any', 'double');

S = S(:).';
Slength = length(S);
data = zeros(bincfg.size, Slength, 'uint8');
data(:) = uint8(255);
data = recursepack(S, bincfg, data);

if nargin < 4
    wora = 'w';
end

if nargin>2 && ~isempty(outputfile)
    fid = fopen(outputfile, wora);
    fwrite(fid, data, 'uint8');
    fclose(fid);
end
    
end


function data = recursepack(S, bincfg, data, offsetshift)
if nargin<4
    offsetshift = 0;
end

cfgfields = fieldnames(bincfg);
[repS, numS] = size(S);
for ii = 1:repS
%     offsetcount = bincfg.offset;
    offsetcount = 0;
    for ifield = 1:length(cfgfields)
        % get field configure
        field_ii = cfgfields{ifield};
        cfg_ii = bincfg.(field_ii);
        if ~isstruct(cfg_ii)
            continue
        end
        % number
        if isempty(cfg_ii.number) || isempty(cfg_ii.size)
            % get the cfg_ii by structbincfg from S(ii,1)
            % WARN: it is not a good idea!
            cfg_tmp = structbincfg(S(ii,1).(field_ii));
            % what? field_ii is not in S? It's a stupid mistake.
            if isempty(cfg_ii.number)
                cfg_ii.number = cfg_tmp.number;
            end
            if isempty(cfg_ii.size)
                cfg_ii.size = cfg_tmp.size;
            end
        end
        fieldsize = cfg_ii.size * cfg_ii.number;
        if ~isavail(cfg_ii.offset)
            cfg_ii.offset = offsetcount;
        else
            offsetcount = cfg_ii.offset;
        end
        % not in S?
        if ~isfield(S, field_ii)
            % offsetcount +=fieldsize；
            offsetcount = offsetcount + fieldsize;
            continue
        end
        % empty field?
        if fieldsize==0
            continue
        end
        % pack the field
        if isfield(cfg_ii, 'class')
            theclass = cfg_ii.class;
        else
            theclass = cfg_ii.dataclass;
        end
        % % check the size with dataclass % not always true
        % if classsize(lower(theclass)) ~= cfg_ii.size
        %     error('The data class of %s is not coupled with the configured size %d.', field_ii, cfg_ii.size);
        % end

        switch lower(theclass)
            case 'struct'
                % to recurse
                Srec = cat(2, S(ii,:).(field_ii));
                if ~isempty(Srec)
                    Srec = reshape(Srec, [], numS);
                    data = recursepack(Srec, cfg_ii, data, offsetshift + offsetcount);
                end
            case 'cell'
                % TBC
                1;
            otherwise
                fielddata = cat(2, S(ii,:).(field_ii));
                if ~isempty(fielddata)
                    % cast to unit8
                    fielddata = castuint8(fielddata(:), theclass);
                    fielddata = reshape(fielddata, [], numS);

                    % auto fix the fielddata size
                    if size(fielddata, 1) < fieldsize
                        % auto fill zeros
                        fielddata(fieldsize, :) = uint8(0);
                    elseif size(fielddata, 1) > fieldsize
                        % auto cut
                        fielddata = fielddata(1 : fieldsize, :);
                    end
                    
                    % write the fielddata to data
                    if size(data, 2) == size(fielddata, 2)
                        data(offsetshift+cfg_ii.offset+(1:fieldsize), :) = fielddata;
                    else
                        minsize = min(size(data, 2), size(fielddata, 2));
                        data(offsetshift+cfg_ii.offset+(1:fieldsize), 1:minsize) = fielddata(:, 1:minsize);
                    end
                end
        end
        % offsetcount +=fieldsize；
        offsetcount = offsetcount + fieldsize;
    end
    offsetshift = offsetshift + bincfg.size;
end
       
end
