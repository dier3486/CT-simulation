function [dataflow, prmflow, status] = reconnode_badchannelprepare(dataflow, prmflow, status)
% corr node, bad channel correction prepare
% [dataflow, prmflow, status] = reconnode_badchannelprepare(dataflow, prmflow, status);

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% calibration table
if isfield(prmflow.corrtable, nodename)
    badcorr = prmflow.corrtable.(nodename);
    badindex_corr = badcorr.badindex(:);
else
    badcorr = [];
    badindex_corr = [];
end

% input bad channel index    
if isfield(nodeprm, 'badindex')
    if ischar(nodeprm.badindex)
        prmflow.pipe.(nodename).badindex = cellfun(@str2num, regexp(nodeprm.badindex, '\d+', 'match'));
    end
else
    prmflow.pipe.(nodename).badindex = [];
end
% unique badindex
prmflow.pipe.(nodename).badindex = unique([prmflow.pipe.(nodename).badindex; badindex_corr(:)]);

% badmax
if ~isfield(nodeprm, 'badmax')
    if isfield(badcorr, 'badmax')
        badmax = badcorr.badmax;
    else
        badmax = 10;
    end
else
    badmax = nodeprm.badmax;
end

% prepare the bad channel correction
badindex = prmflow.pipe.(nodename).badindex;
Nbadchennel = length(badindex);
prmflow.pipe.(nodename).Nbadchennel = Nbadchennel;
prmflow.pipe.(nodename).interpindex = zeros(Nbadchennel, 2);
prmflow.pipe.(nodename).interpalpha = zeros(Nbadchennel, 2, 'single');
Npixel = prmflow.raw.Npixel;
% loop the bad channels
for ii = 1:length(badindex)
    ibad = badindex(ii);
    % to find the left good one
    for jj = 1:badmax
        i_left = ibad-jj;
        if any(i_left==badindex)
            % bad again
            i_left = [];
        elseif mod(i_left, Npixel) == 0
            % touch the edge
            i_left = [];
            a_left = inf;
            break;
        else
            % got
            a_left = jj;
            break;
        end
    end

    % to find the right good one
    for jj = 1:badmax
        i_right = ibad+jj;
        if any(i_right==badindex)
            % bad again
            i_right = [];
        elseif mod(i_right, Npixel) == 1
            % touch the edge
            i_right = [];
            a_right = inf;
            break;
        else
            % got
            a_right = jj;
            break;
        end
    end
    % interp weight
    alpha = [a_right, a_left]./(a_right+a_left);
    alpha(isnan(alpha)) = 0;

    % save in prmflow
    prmflow.pipe.(nodename).interpindex(ii, :) = [i_left i_right];
    prmflow.pipe.(nodename).interpalpha(ii, :) = single(alpha);
end

% pipe line
if pipeline_onoff
    % pipeline console paramters, default
    prmflow.pipe.(nodename).pipeline = struct();
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end