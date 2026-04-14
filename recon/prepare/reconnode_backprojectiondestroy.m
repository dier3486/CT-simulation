function [dataflow, prmflow, status] = reconnode_backprojectiondestroy(dataflow, prmflow, status)
% destroy node if backprojection
% [dataflow, prmflow, status] = reconnode_backprojectiondestroy(dataflow, prmflow, status);

% Copyright Dier Zhang
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% parameters set in pipe
nodename = status.nodename;
nextnode = status.pipeline.(nodename).nextnode;

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% echo onoff
echo_onoff = status.debug.echo_onoff;

% ini status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

% CUDA on/off
if ~pipeline_onoff || ~prmflow.pipe.(nodename).pipeline.CUDAonoff
    % nothing to destroy
    return;
end
% % call cuda
% if ~ismember(lower(prmflow.recon.method), {'helicalpiline'})
%     status.errorcode = 710;
%     status.errormsg = 'The CUDA BP (destroy) function only support helical-piline yet.';
%     return;
% end

% cuda detroy
if libisloaded('HelicalBackprojection')
    % CUDA Destroy
    calllib('HelicalBackprojection', 'OuterCudaDestroy', dataflow.buffer.(nodename).CUDA.Wrecon, ...
        dataflow.buffer.(nodename).CUDA.devPointers.Value);
end
if libisloaded('CUpipeline')
    % free the private buffers
    calllib('CUpipeline', 'cudaworkspaceFree', dataflow.buffer.(nodename).CUDA.Wreconstruct);
    calllib('CUpipeline', 'cudaworkspaceFree', dataflow.buffer.(nodename).CUDA.rawdata);
    calllib('CUpipeline', 'cudaworkspaceFree', dataflow.buffer.(nodename).CUDA.viewangle);
    % calllib('CUpipeline', 'cudaworkspaceFree', dataflow.buffer.(nodename).CUDA.imageout);
    1;
    % free the (next) pools (hard code)
    if isa(dataflow.pipepool.(nextnode)(1).data.image, 'lib.pointer')
        calllib('CUpipeline', 'cudaworkspaceFree', dataflow.pipepool.(nextnode)(1).data.image);
    end
    % We shall move the freeing of pools to reconnode_pipelinedestroy.m
end

if echo_onoff
    fprintf(' (%.2fms)', dataflow.buffer.(nodename).CUDA.timecost);
end

% clear buffer.CUDA
dataflow.buffer.(nodename) = rmfield(dataflow.buffer.(nodename), 'CUDA');


end

