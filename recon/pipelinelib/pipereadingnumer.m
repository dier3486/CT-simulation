function outdata = pipereadingnumer(outdata, outpool, plconsol, headfield)

if nargin < 4
    headfield = 'rawhead';
end

index_out = plconsol.Index_out(1) : plconsol.Index_out(2);
if outpool.circulatemode
    index_out = mod(index_out - 1, outpool.poolsize) + 1;
end

if isfield(outdata, headfield)
    % Reading Number
    if isfield(outdata.(headfield), 'Reading_Number')
        readingindex = (plconsol.Index_out(1) : plconsol.Index_out(2)) - outpool.ReadStart + 1;
        if outpool.circulatemode
            readingindex = mod(readingindex-1, outpool.poolsize) + 1;
        end
        outdata.(headfield).Reading_Number(index_out) = readingindex;
    end
    
    if isfield(outdata.(headfield), 'Shot_Start')
        % Shot Start
        outdata.(headfield).Shot_Start(index_out) = 0;
        ReadStart = outpool.ReadStart;
        if ReadStart >= plconsol.Index_out(1)
            if outpool.circulatemode
                ReadStart = mod(ReadStart-1, outpool.poolsize) + 1;
            end
            outdata.(headfield).Shot_Start(ReadStart) = 1;
        end

        % Shot End
        if plconsol.isshotend
            ReadEnd = outpool.ReadEnd;
            if outpool.circulatemode
                ReadEnd = mod(ReadEnd-1, outpool.poolsize) + 1;
            end
            outdata.(headfield).Shot_Start(ReadEnd) = -1;
        end
    end
end

end