function H = loadfilter(filter, Npixel, delta_d)
% load filter
% H = loadfilter(filter, Npixel, delta_d);

% design filter
if isfield(filter, 'name')  || ischar(filter)
    % design filter by filter's name
    if isfield(filter, 'freqscale')
        freqscale = filter.freqscale;
    else
        freqscale = 1.0;
    end
    if ischar(filter)
        filtname = filter;
    else
        filtname = filter.name;
    end
    H = filterdesign(filtname, Npixel, delta_d, freqscale);
elseif isfield(filter, 'file')
    % load filter from a file
    if ~exist(filter.file, 'file')
        error('Can not find filter file: %s', filter.file);
    end
    fid = fopen(filter.file, 'r');
    H = fread(fid, inf, 'single');
    fclose(fid);
    H = H(:);
else
    H = struct([]);
    fields = fieldnames(filter);
    for ii = 1:length(fields)
        if isstruct(filter.(fields{ii}))
            Hi = loadfilter(filter.(fields{ii}), Npixel, delta_d);
            if ~isempty(Hi)
                H(1).(fields{ii}) = Hi;
            end
        end
    end
end

end