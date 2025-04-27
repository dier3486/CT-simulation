function pipeOut = combinenodes(pipeIn)
% to 'combine' the pipeline configure with loopings
%   pipeOut = combinenodes(pipeIn);
% e.g. prmflow.pipe = combinenodes(prmflow.pipe);

pipenodes = fieldnames(pipeIn);
pipeOut = struct();

for ifield = pipenodes(:)'
    nodename_slip = regexp(ifield{1}, '_', 'split');
    currfields = fieldnames(pipeOut);
    Ncf = length(currfields);
    if strcmpi(nodename_slip{1}, 'loop')
        gotonode = pipeIn.(ifield{1}).goto;
        gotoindex = find(strcmp(gotonode, currfields), 1);
        if isempty(gotoindex)
            continue;
        end
        for ii = 1:pipeIn.(ifield{1}).times-1
            for jj = 0:Ncf-gotoindex
                newfname = [currfields{gotoindex+jj} sprintf('_%02d', ii)];
                pipeOut.(newfname) = pipeOut.(currfields{gotoindex+jj});
            end
        end
    else
        pipeOut.(ifield{1}) = pipeIn.(ifield{1});
        if isfield(pipeIn.(ifield{1}), 'loop')
            for ii = 1:pipeIn.(ifield{1}).loop-1
                newfname = [ifield{1} sprintf('_%02d', ii)];
                pipeOut.(newfname) = pipeIn.(ifield{1});
            end
        end
    end
end

end