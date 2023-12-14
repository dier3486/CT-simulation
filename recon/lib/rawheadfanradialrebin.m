function rawhead = rawheadfanradialrebin(rawhead, rebin)
% viewangle/2 for DFS in reconnode_FanRadialrebin

headfields = fieldnames(rawhead);
Nfields = length(headfields);

for ifield = 1:Nfields
    datasize = size(rawhead.(headfields{ifield}), 1);
    switch headfields{ifield}
        case 'refblock'
            rawhead.refblock = squeeze(any(reshape(rawhead.refblock, datasize, rebin.Nfocal, []), 2));
        case 'viewangle'
            rawhead.viewangle = rawhead.viewangle(1:rebin.Nfocal:end) + rebin.DFSviewshift;
        case {'mA', 'KV'}
            rawhead.(headfields{ifield}) = mean(reshape(rawhead.(headfields{ifield}), rebin.Nfocal, []), 1);
        case 'Shot_Number'
            rawhead.Shot_Number = rawhead.Shot_Number(1:rebin.Nfocal:end);
        case 'Reading_Number'
            rawhead.Reading_Number = rawhead.Reading_Number(1:rebin.Nfocal:end);
        otherwise
            % remove the useless fields
            rawhead = rmfield(rawhead, headfields{ifield});
    end

end

end