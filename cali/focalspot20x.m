function focalspot_0x = focalspot20x(focalspot)
% to explain the focalspot in protocol
% focalspot_0x = focalspot20x(focalspot);
% QFS:  0x10000000   (left end)
% DFS:  0x01100000

if isnumeric(focalspot)
    focalspot_0x = focalspot;
else
    switch upper(focalspot)
        case 'QFS'
            % the 1st focal spot
            focalspot_0x = 1;       % 0x10000000
        case {'DFS', 'XDFS'}
            % flying between 2nd and 3rd focal sopt
            focalspot_0x = 6;       % 0x01100000
        case {'ZDFS'}
            % flying between 4th and 5th focal sopt
            focalspot_0x = 24;      % 0x00011000
        otherwise
            % unknown focal spot
            focalspot_0x = 0;
    end     
end

focalspot_0x = uint16(focalspot_0x);

end