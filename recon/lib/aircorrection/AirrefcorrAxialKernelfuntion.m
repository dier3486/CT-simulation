function [dataflow, reflast, Nrenew] = AirrefcorrAxialKernelfuntion(dataflow, prmflow, status, shotindex, buffer)
% recon node, Air reference correction function
% [dataflow, reflast, Nrenew] = AirrefcorrAxialKernelfuntion(dataflow, prmflow, status, buffer);
% where the buffer mostly is dataflow.buffer.(nodename)

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

% nodename
nodename = status.nodename;
% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% parameters set in pipe
airprm = prmflow.pipe.(nodename);
% slice independent refernece?
if isfield(airprm, 'sliceindependent')
    sliceindependent = airprm.sliceindependent;
else
    sliceindependent = false;
end

Nviewprot = prmflow.raw.Nviewprot;
Nfocal = prmflow.raw.Nfocal;
blockwindow = prmflow.raw.air.blockwindow;
refpixel = prmflow.raw.air.refpixelindex;
Npixel = prmflow.raw.Npixel;
Nslice = prmflow.raw.Nslice;
% airKVmA = prmflow.raw.air.airKVmA;

if pipeline_onoff
    % 'reflast' is the last reference value of the previous datablock, be used as a left boundary condition.
    reflast = buffer.reflast;
    % The 'refblock' is a logic array to sign if the referece(s) are blocked;   
    refblock = dataflow.rawhead.refblock;
    % I know the refblock is initialed in prepare.
    % Nrenew to return
    dataendpoint = mod(buffer.WritePoint + buffer.Nrenew - 2, Nviewprot) + 1;
    if dataendpoint > Nviewprot - blockwindow
        Nrenew = 0;
    elseif dataendpoint < Nviewprot - blockwindow
        Nrenew = max(0, dataendpoint - buffer.AvailPoint - blockwindow);
    else % dataendpoint == Nviewprot - blockwindow
        Nrenew = dataendpoint - buffer.AvailPoint + blockwindow;
    end
    % viewindex to calculate the reference and blocked reference
    if buffer.WritePoint > Nviewprot - blockwindow
        viewindex = mod((buffer.WritePoint : buffer.WritePoint + max(buffer.Nrenew, Nrenew)) - 1, Nviewprot) + 1;
        viewindex_rawref = (1 : Nrenew) + (Nviewprot - buffer.WritePoint + 1);
        datastartpoint = buffer.WritePoint;
    else
        viewindex = buffer.AvailPoint+1 : buffer.AvailPoint + Nrenew;
        viewindex_rawref = 1 : Nrenew;
        datastartpoint = buffer.AvailPoint+1;
    end
    
else
    % ini for normal mode
    reflast = cell(1, Nfocal);
    refblock = [];
    Nrenew = Nviewprot;
    viewindex = dataflow.rawhead.Shot_Number == shotindex;
    datastartpoint = find(viewindex, 1, 'first');
    
end
    
% % KVmA
% KVmA_fix = -log2(dataflow.rawhead.KV(:, viewindex).* ...
%     dataflow.rawhead.mA(:, viewindex)) - airKVmA;
% refernce by raw the reference error
[rawref, referr] = airreference2(dataflow.rawdata(:, viewindex), refpixel, Npixel, Nslice, ...
    sliceindependent);

refblock = isrefblocked(refblock, prmflow, referr, datastartpoint);

if pipeline_onoff
    viewindex_refblk = buffer.AvailPoint+1 : buffer.AvailPoint+Nrenew;
    viewindex_renew = viewindex_refblk;
else
    % boundary circshift
    viewindex_refblk = [blockwindow+1:Nrenew  1:blockwindow];
    viewindex_rawref = viewindex_refblk;
    viewindex_renew = viewindex_refblk + datastartpoint - 1;
end

% first reading number
reading1 = dataflow.rawhead.Reading_Number(viewindex_renew(1));
% loop the focals (DFS) to correct the rawdata
for ifocal = 1:Nfocal
    % the first view is focal 1 or focal 2? depending on the first reading number
    ifocal_index = (mod(reading1(1), 1) ~= (ifocal-1)) + 1;
    % view index to work
    viewindex_rawref_ifc = viewindex_rawref(ifocal_index:Nfocal:Nrenew);
    viewindex_ref_ifc = viewindex_refblk(ifocal_index:Nfocal:Nrenew);
    viewindex_renew_ifc = viewindex_renew(ifocal_index:Nfocal:Nrenew);
    % call reference correction funtion
    [rawref_ifc, reflast{ifocal}] = referencecorr(rawref(:, viewindex_rawref_ifc), refblock(:, viewindex_ref_ifc), ...
        reflast{ifocal});
    if sliceindependent
        dataflow.rawdata(:, viewindex_renew_ifc) = dataflow.rawdata(:, viewindex_renew_ifc) - repelem(rawref_ifc, Npixel, 1);
    else
        dataflow.rawdata(:, viewindex_renew_ifc) = dataflow.rawdata(:, viewindex_renew_ifc) - rawref_ifc;
    end
end
% to return the refblock
if pipeline_onoff
    dataflow.rawhead.refblock = refblock;
else
    dataflow.rawhead.refblock = [dataflow.rawhead.refblock refblock];
end

end