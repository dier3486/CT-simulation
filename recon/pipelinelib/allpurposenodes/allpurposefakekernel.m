function [dataflow, Rinfo] = allpurposefakekernel(dataflow, prmflow, status, Ininfo)
% a fake node kernel function

% parameters set in pipe
nodename = status.nodename;
if isfield(prmflow.pipe, nodename)
    nodeprm = prmflow.pipe.(nodename);
else
    nodeprm = struct();
end

% operator
if isfield(nodeprm, 'mainopt')
    mainopt = prmflow.pipe.(nodename).mainopt;
else
    mainopt = @(x) x;
end

if nargin<4
    if isfield(dataflow, 'rawdata')
        dataflow.rawdata = mainopt(dataflow.rawdata);
    end
    if isfield(dataflow, 'image')
        dataflow.image = mainopt(dataflow.image);
    end
    Rinfo = [];
else
    datastart = Ininfo.StartPoint;
    dataend = Ininfo.StartPoint + Ininfo.AvailNumber - 1;

    if isfield(dataflow, 'rawdata')
        dataflow.rawdata(:, datastart:dataend) = mainopt(dataflow.rawdata(:, datastart:dataend));
        % not support circulation
    end
%     if isfield(dataflow, 'image') && (nargin>2) && ~isempty(imagestartend)
%         dataflow.image(:, imagestartend(1):imagestartend(2)) = dataflow.image(:, imagestartend(1):imagestartend(2)) + 1;
%     end
    % return Rinfo, ummm seems nothing to return
    Rinfo = Ininfo;
end

end