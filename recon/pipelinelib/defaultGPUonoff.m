function prmflow = defaultGPUonoff(prmflow, status, nodename)
% to set the GPUonoff related values
% [prmflow, status] = defaultGPUonoff(prmflow, status, nodename).
% call me in node prepare or defaultpipelineprepare.

if strcmpi(nodename, 'NULL')
    return;
end

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

if ~prmflow.protocol.gpuonoff
    % GPU closed
    GPUonoff = false;
elseif prmflow.pipe.(nodename).pipeline.GPUonoff == -1
    % default GPUonoff
    if pipeline_onoff
        % to follow the previous node
        prevnode = status.pipeline.(nodename).prevnode;
        Isprevexist = ~strcmpi(prevnode, 'NULL') && status.pipeline.(prevnode).pipeline_onoff;
        if Isprevexist && prmflow.pipe.(prevnode).pipeline.nextgpuDevice
            % GPU on
            GPUonoff = true;
        else
            GPUonoff = false;
        end
    else
        % GPU on
        GPUonoff = true;
    end
else
    GPUonoff = logical(prmflow.pipe.(nodename).pipeline.GPUonoff);   % 0 or 1
end

prmflow.pipe.(nodename).pipeline.GPUonoff = GPUonoff;

% if GPUonoff
%     prmflow.pipe.(nodename).pipeline.GPUonoff = 1;
%     % only 1 GPU devices can be employed now, not opened to configurable yet.
%     prmflow.pipe.(nodename).pipeline.nextgpuDevice = prmflow.system.GPUdevices(1);
% else
%     prmflow.pipe.(nodename).pipeline.GPUonoff = 0;
%     prmflow.pipe.(nodename).pipeline.nextgpuDevice = 0;
% end

end