function [outpool, prmflow, status] = allpurposefakekernel(outpool, inpool, prmflow, status)
% a fake node kernel function


% parameters set in pipe
nodename = status.nodename;
if isfield(prmflow.pipe, nodename)
    nodeprm = prmflow.pipe.(nodename);
else
    nodeprm = struct();
end
% nextnode
nextnode = status.pipeline.(nodename).nextnode;

% object type (input)
objecttype = lower(nodeprm.pipeline.inputobjecttype);

% operator
if isfield(nodeprm, 'mainopt')
    mainopt = prmflow.pipe.(nodename).mainopt;
else
    mainopt = @(x) x;
end

if status.currentjob.pipeline_onoff
    % pipe-line mode

    % viewrescale and other
    viewrescale = nodeprm.pipeline.viewrescale;
    viewrely = nodeprm.pipeline.viewrely_out;
    relystrategy = nodeprm.pipeline.relystrategy;
    viewextra = nodeprm.pipeline.viewextra;
    plconsol = status.currentjob.pipeline;
    
    kernellevel = ~isempty(inpool);
    % step 0, pool and index
    % last node
    if isempty(outpool)
        % no where to output
        status.jobdone = 1;
        return;
    end
    if kernellevel==0
        inpool = outpool;
    end
    % Index_in/Index_out
    index_in = poolindex(inpool, plconsol.Index_in);
    index_out = poolindex(outpool, plconsol.Index_out);
    n_inres = length(index_in)*viewrescale(1)/viewrescale(2);
    n_out = length(index_out);
    
    % step 1
    if kernellevel
        % special A-H.1.S
        if inpool.circulatemode && ~outpool.circulatemode && relystrategy==2
            viewextra = -viewrely;
        end
        % pick rescaled data
        indata_res = datarescale(inpool.data, index_in, viewrescale);
        % data copy
        in_start = max(1, 1 + plconsol.Do2i);
        out_start = max(1, 1 - plconsol.Do2i);
        in_end = min(n_inres, n_out + plconsol.Do2i);
        out_end = in_end - plconsol.Do2i;
        if relystrategy==1
            % G
            outpool.data = kerneldatacum(outpool.data, indata_res, index_out(out_start : out_end), in_start : in_end);
        else
            outpool.data = kerneldatacopy(outpool.data, indata_res, index_out(out_start : out_end), in_start : in_end);
        end
    end
    
    % step 2
    % main opt
    if relystrategy==0 || relystrategy==1
        % N or G
        out_start = max(1, 1 - plconsol.Do2i);
        out_end = min(n_inres - plconsol.Do2i, n_out);
    else
        % S
        out_start = 1;
        out_end = n_out;
    end
    if isfield(outpool.data, objecttype)
        outpool.data.(objecttype)(:, index_out(out_start:out_end)) = ...
            mainopt(outpool.data.(objecttype)(:, index_out(out_start:out_end)));
    end

    % step 3
    % fill up viewextra
    if ~outpool.circulatemode
        if plconsol.isshotstart && viewextra(1) > 0
            outpool.data = kerneldatacopy(outpool.data, outpool.data, ...
                index_out(1 : viewextra(1)), index_out(viewextra(1)+1) );
        end
        if plconsol.isshotend && viewextra(2) > 0
            outpool.data = kerneldatacopy(outpool.data, outpool.data, ...
                index_out(n_out-viewextra(2)+1 : n_out), index_out(n_out-viewextra(2)) );
        end
    end

    % step 4
    % view rely
    if relystrategy==1
        % G only
        % left
        out_start = 1;
        if plconsol.isshotend && ~outpool.circulatemode
            out_end = n_out;
        else
            % out_end = n_out - viewrely(1) - viewrely(2);
            out_end = n_inres - plconsol.Do2i - viewrely(2);
        end
        if isfield(outpool.data, objecttype)
            outpool.data.(objecttype)(:, index_out(out_start:out_end)) = ...
                outpool.data.(objecttype)(:, index_out(out_start:out_end)) + 0.5;
        end
        % right
        if plconsol.isshotstart && ~outpool.circulatemode
            out_start = 1;
        else
            % out_start = 1 + viewrely(1) + viewrely(2);
            out_start = 1 - plconsol.Do2i + viewrely(1);
        end
        out_end = n_out;
        if isfield(outpool.data, objecttype)
            outpool.data.(objecttype)(:, index_out(out_start:out_end)) = ...
                outpool.data.(objecttype)(:, index_out(out_start:out_end)) - 0.5;
        end
    end

    % step 5, reset Reading_Number and Shot_Start (select)
    outpool.data = pipereadingnumer(outpool.data, outpool, plconsol, 'rawhead');

else
    % non-pipeline, I know the outpool = inpool = dataflow.
    outpool = inpool;
    if isfield(outpool, objecttype)
        outpool.(objecttype) = mainopt(outpool.(objecttype));
    end
end

end


function resdata = datarescale(indata, index_in, Q)

if nargin<3
    Q = [1 1];
end

n = length(index_in);
m = n*Q(1)/Q(2);
fields = fieldnames(indata);
resdata = struct();
for ii = 1:length(fields)
    if isstruct(indata.(fields{ii}))
        resdata.(fields{ii}) = datarescale(indata.(fields{ii}), index_in, Q);
    else
        sizein = size(indata.(fields{ii}));
        sizein(2) = n;
        sizeout = sizein;
        sizeout(2) = m;
        A = cast(mean(reshape(indata.(fields{ii})(:, index_in), sizein(1), Q(2), []), 2), 'like', indata.(fields{ii}));
        resdata.(fields{ii}) = reshape(repelem(A, 1, Q(1)), sizeout);
    end
end

end

function outdata = kerneldatacopy(outdata, indata, index_out, index_in, copyfields)

if nargin < 5 || isempty(copyfields)
    copyfields = fieldnames(indata);
end
for ii = 1:length(copyfields)
    if isempty(copyfields{ii})
        % pass the {''}
        continue
    end
    if ~isfield(indata, copyfields{ii})
        % pass the not exist fields
        continue;
    end
    if isstruct(indata.(copyfields{ii}))
        if ~isfield(outdata, copyfields{ii})
            outdata.(copyfields{ii}) = struct();
        end
        outdata.(copyfields{ii}) = kerneldatacopy(outdata.(copyfields{ii}), indata.(copyfields{ii}), index_out, index_in);
    else
        if ~isfield(outdata, copyfields{ii})
            outdata.(copyfields{ii}) = cast([], 'like', indata.(copyfields{ii}));
        end
        if length(index_in) ~= 1
            outdata.(copyfields{ii})(:, index_out) = indata.(copyfields{ii})(:, index_in);
        else
            outdata.(copyfields{ii})(:, index_out) = repmat(indata.(copyfields{ii})(:, index_in), 1, length(index_out));
        end
    end
end

end

function outdata = kerneldatacum(outdata, indata, index_out, index_in, copyfields)

if nargin < 5 || isempty(copyfields)
    copyfields = fieldnames(indata);
end
for ii = 1:length(copyfields)
    if isempty(copyfields{ii})
        % pass the {''}
        continue
    end
    if ~isfield(indata, copyfields{ii})
        % pass the not exist fields
        continue;
    end
    if isstruct(indata.(copyfields{ii}))
        if ~isfield(outdata, copyfields{ii})
            outdata.(copyfields{ii}) = struct();
        end
        outdata.(copyfields{ii}) = kerneldatacum(outdata.(copyfields{ii}), indata.(copyfields{ii}), index_out, index_in);
    else
        if ~isfield(outdata, copyfields{ii})
            outdata.(copyfields{ii}) = cast([], 'like', indata.(copyfields{ii}));
            outdata.(copyfields{ii})(:, index_out) = indata.(copyfields{ii})(:, index_in);
        else
            outdata.(copyfields{ii})(:, index_out) = ...
                outdata.(copyfields{ii})(:, index_out) + indata.(copyfields{ii})(:, index_in);
        end
        
    end
end

end