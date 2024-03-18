function [dataflow, prmflow, status] = fakenodereadrawdata(dataflow, prmflow, status)

% parameters set in pipe
nodename = status.nodename;

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% node main function
if pipeline_onoff
    [dataflow.buffer.(nodename), dataflow.buffer.(nodename).outputpool.data, prmflow, status] = ...
        fakekernelreadrawdata(dataflow.buffer.(nodename), dataflow.buffer.(nodename).outputpool.data, ...
        prmflow, status);
else
    [~, dataflow, prmflow, status] = fakekernelreadrawdata([], dataflow, prmflow, status);
end
    
end

