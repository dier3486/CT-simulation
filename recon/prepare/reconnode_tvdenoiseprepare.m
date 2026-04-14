function [dataflow, prmflow, status] = reconnode_tvdenoiseprepare(dataflow, prmflow, status)
% post recon node, TV denoise prepare
% [dataflow, prmflow, status] = reconnode_tvdenoiseprepare(dataflow, prmflow, status);

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% I konw the prmflow.image was prepared by recon nodes, but if missed call this
if isempty(fieldnames(prmflow.image))
    prmflow.image = protocol2poster([], prmflow.protocol);
end

% mu
if isfield(nodeprm, 'mu') && ~isempty(nodeprm.mu)
    prmflow.image.TV_mu = nodeprm.mu;
else
    prmflow.image.TV_mu = 0.2;
end
% Cl
if isfield(nodeprm, 'Cl') && ~isempty(nodeprm.Cl)
    prmflow.image.TV_Cl = nodeprm.Cl;
else
    if isfield(nodeprm, 'lambda') && ~isempty(nodeprm.lambda)
        prmflow.image.TV_Cl = nodeprm.lambda / prmflow.image.TV_mu;
    else
        prmflow.image.TV_Cl = 0.1;
    end
end

% default denoise parameters
if prmflow.pipe.(nodename).pipeline.CUDAonoff
    % tmply
    nodeprm.DIM = "2D";
    nodeprm.MaxIter = 1;
else

end

% TV opstions
defaultTVopts = TVoptions();
nodeTVopts = rmfield(nodeprm, setdiff(fieldnames(nodeprm), fieldnames(defaultTVopts)));
nodeTVopts = everything2single(nodeTVopts, 'char', 'string');
nodeTVopts.LeftDim = 3;
prmflow.image.TVopts = structmerge(nodeTVopts, defaultTVopts);
% plz use namedargs2cell(TVopts) in calling BregmanTV.

% image rely
if isfield(prmflow.image.TVopts, 'DIM') && strcmp(prmflow.image.TVopts.DIM, '2D')
    prmflow.pipe.(nodename).imagerely = [0 0];
else
    if ~isfield(nodeprm, 'imagerely') || isempty(nodeprm.imagerely)
        prmflow.pipe.(nodename).imagerely = [1 1];
    end
end

% pipe line
if pipeline_onoff
    % pipeline console paramters
    % the denoise is H-H.0.S
    prmflow.pipe.(nodename).pipeline.kernellevel = 0;
    prmflow.pipe.(nodename).pipeline.viewrely = prmflow.pipe.(nodename).imagerely;
    if any(prmflow.pipe.(nodename).imagerely~=0)
        prmflow.pipe.(nodename).pipeline.relystrategy = 'stingy';
    else
        prmflow.pipe.(nodename).pipeline.relystrategy = 0;
    end
    % inputminlimit
    if isfield(nodeprm, 'inputminlimit')
        prmflow.pipe.(nodename).pipeline.inputminlimit = nodeprm.inputminlimit;
    else
        prmflow.pipe.(nodename).pipeline.inputminlimit = 8;
    end
    % ask next pool to keep some slices of the images
    prmflow.pipe.(nodename).pipeline.nextkeepbottom = prmflow.pipe.(nodename).imagerely(1);

    % % carried (should be)
    % if prmflow.protocol.tocarrythepools     % default is true
    %     prmflow.pipe.(nodename).pipeline.iscarried = true;
    %     % default was false
    % end
end

% GPU on/off
prmflow = defaultGPUonoff(prmflow, status, nodename);
% while GPU on
if prmflow.pipe.(nodename).pipeline.GPUonoff > 0
    GPUfields = {'TV_mu', 'TV_Cl'};
    prmflow.image = putfieldsinGPU(prmflow.image, GPUfields);
end

% CUDA on/off
if prmflow.pipe.(nodename).pipeline.CUDAonoff
    % call cuda
    [dataflow, prmflow, status] = cudatvdenpiseprepare(dataflow, prmflow, status);
end


% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end