function [dataflow, reflast, Nrenew] = AirrefcorrHelicalKernelfuntion(dataflow, prmflow, status, buffer)
% recon node, Air reference correction function
% [dataflow, reflast, Nrenew] = AirrefcorrHelicalKernelfuntion(dataflow, prmflow, status, buffer);
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

Nview = prmflow.raw.Nview;
% Nviewprot = prmflow.raw.Nviewprot;
Nfocal = prmflow.raw.Nfocal;
% Nshot = prmflow.raw.Nshot;
blockwindow = prmflow.raw.air.blockwindow;
refpixel = prmflow.raw.air.refpixelindex;
Npixel = prmflow.raw.Npixel;
Nslice = prmflow.raw.Nslice;

if pipeline_onoff
    datastartpoint = buffer.AvailPoint + 1;
    dataendpoint = buffer.WritePoint + buffer.Nrenew - 1;
else
    Nrenew = size(dataflow.rawdata, 2);
    datastartpoint = 1;
    dataendpoint = datastartpoint + Nrenew - 1;
end
    
% KVmA
% KVmA = -log2(dataflow.rawhead.KV(:, datastartpoint:dataendpoint).* ...
%     dataflow.rawhead.mA(:, datastartpoint:dataendpoint));
% refernce by raw the reference error
% [rawref, referr] = airreference2(dataflow.rawdata(:, datastartpoint:dataendpoint), refpixel, Npixel, Nslice, ...
%     sliceindependent);

[rawref, referr] = airreference2(dataflow.rawdata(:, datastartpoint:dataendpoint), refpixel, Npixel, Nslice, ...
    sliceindependent);

if pipeline_onoff
    % 'reflast' is the last reference value of the previous datablock, be used as a left boundary condition.
    reflast = buffer.reflast;
    % The 'refblock' is a logic array to sign if the referece(s) are blocked;   
    % 'new' refblock
    Nref = size(dataflow.rawhead.refblock, 1);
    Nrefblock = buffer.WritePoint+blockwindow-1;
    refblock = [dataflow.rawhead.refblock(:, 1:Nrefblock) false(Nref, buffer.Nrenew)];
    
    % new avialviews
    switch prmflow.raw.scan
        case 'static'
            Nrenew = buffer.Nrenew;
        case {'helical', 'halfaxial'}
            if buffer.AvailViewindex + buffer.Nrenew + blockwindow>= Nview
                % reaching the end
                Nrenew = Nview - buffer.AvailViewindex;
            else
                Nrenew =  buffer.WritePoint + buffer.Nrenew - blockwindow - buffer.AvailPoint - 1;
            end
        otherwise
            error('Can not run the air correction for %s in the Helical Kernelfuntion!', prmflow.raw.scan);
    end
else
    % ini for normal mode
    reflast = cell(1, Nfocal);
    refblock = [];
end

refblock = isrefblocked(refblock, prmflow, referr, datastartpoint);

% first reading number
reading1 = dataflow.rawhead.Reading_Number(datastartpoint);
for ifocal = 1:Nfocal
    % the first view is focal 1 or focal 2?
    ifocal_index = (mod(reading1, 1) ~= (ifocal-1)) + 1;
    viewindex_ifc = (ifocal_index:Nfocal:Nrenew);
    viewindex_renew = viewindex_ifc + datastartpoint - 1;
%     [rawref_ifc, reflast{ifocal}] = referencecorr(rawref(:, viewindex_ifc), refblock(:, viewindex_renew), ...
%         KVmA(:, viewindex_ifc), reflast{ifocal});
    [rawref_ifc, reflast{ifocal}] = referencecorr(rawref(:, viewindex_ifc), refblock(:, viewindex_renew), ...
        reflast{ifocal});
    if sliceindependent
        dataflow.rawdata(:, viewindex_renew) = dataflow.rawdata(:, viewindex_renew) - repelem(rawref_ifc, Npixel, 1);
    else
        dataflow.rawdata(:, viewindex_renew) = dataflow.rawdata(:, viewindex_renew) - rawref_ifc;
    end
end

% over write the refblock
dataflow.rawhead.refblock = refblock;

end