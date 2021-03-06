function [data, bincfg] = packstruct(S, bincfg, outputfile)
% transform a struct to bin data and/or write in a file
% [data, bincfg] = packstruct(S, bincfg, outputfile);
% or, packstruct(S, bincfg, outputfile);
% INPUT:
%   S               the structure to pack
%   bincfg          the configure of the data format
%   outputfile      output file name (path)
% OUTPUT:
%   data            the bin data to (or did) save in file
%   bincfg          the returned configure of the data format
%   write the outputfile
% NOTE: If you don't know how to get the bincfg, try this: bincfg = readcfgfile(cfgmatchrule(outputfile, cfgpath));
% where the cfgpath is the folder of the format configure files, ~/IO/standard/;

if nargin<2
    bincfg = structbincfg(S);
else
    bincfg = clearbincfg(S(1), bincfg);
end
% to double
bincfg = everything2single(bincfg, 'any', 'double');

S = S(:).';
Slength = length(S);
data = zeros(bincfg.size, Slength, 'uint8');
data(:) = uint8(255);
data = recursepack(S, bincfg, data);

if nargin>2
    fid = fopen(outputfile, 'w');
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
    offsetcount = bincfg.offset;
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
        if isnan(cfg_ii.offset)
            cfg_ii.offset = offsetcount;
        else
            offsetcount = cfg_ii.offset;
        end
        offsetcount = offsetcount + fieldsize;
        % not in S?
        if ~isfield(S, field_ii)
            continue
        end
        % empty field?
        if fieldsize==0
            continue
        end
        % pack the field
        switch lower(cfg_ii.class)
            case 'struct'
                % to recurse
                Srec = cat(2, S(ii,:).(field_ii));
                if ~isempty(Srec)
                    Srec = reshape(Srec, [], numS);
                    data = recursepack(Srec, cfg_ii, data, offsetshift);
                end
            case 'cell'
                % TBC
                1;
            otherwise
                fielddata = cat(2, S(ii,:).(field_ii));
                if ~isempty(fielddata)
                    % fielddata = typecast(fielddata(:), 'uint8');
                    fielddata = castuint8(fielddata(:), cfg_ii.class);
                    fielddata = reshape(fielddata, fieldsize, []);
                    data(offsetshift+cfg_ii.offset+(1:fieldsize), :) = fielddata;
                end
        end
    end
    offsetshift = offsetshift + bincfg.size;
end
       
end
