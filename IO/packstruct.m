function [data, bincfg] = packstruct(S, bincfg)
% transform a struct to bin data
% [data, bincfg] = packstruct(S, bincfg);

if nargin<2
    bincfg = structbincfg(S);
else
    bincfg = clearbincfg(S(1), bincfg);
end

S = S(:).';
Slength = length(S);
data = zeros(bincfg.size, Slength, 'uint8');
data = recursepack(S, bincfg, data);
    
end

function data = recursepack(S, bincfg, data, offsetshift)
if nargin<4
    offsetshift = 0;
end

cfgfields = fieldnames(bincfg);
[repS, numS] = size(S);
for ii = 1:repS
    for ifield = 1:length(cfgfields)
        field_ii = cfgfields{ifield};
        cfg_ii = bincfg.(field_ii);
        if ~isstruct(cfg_ii)
            continue
        elseif ~isfield(S, field_ii)
            continue
        end
        switch cfg_ii.class
            case 'struct'
                % to recurse
                Srec = reshape(cat(2, S(ii,:).(field_ii)), [], numS);
                data = recursepack(Srec, cfg_ii, data, offsetshift);
            case 'cell'
                % TBC
                1;
            otherwise
                fielddata = cat(2, S(ii,:).(field_ii));
                %                     fielddata = typecast(fielddata(:), 'uint8');
                fielddata = castuint8(fielddata(:), cfg_ii.class);
                fieldsize = cfg_ii.size * cfg_ii.number;
                fielddata = reshape(fielddata, fieldsize, []);
                data(offsetshift+cfg_ii.offset+(1:fieldsize), :) = fielddata;
        end
    end
    offsetshift = offsetshift + bincfg.size;
end
       
end
