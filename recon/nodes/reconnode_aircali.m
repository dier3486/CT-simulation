function [dataflow, prmflow, status] = reconnode_aircali(dataflow, prmflow, status)
% air calibration
% [dataflow, prmflow, status] = reconnode_aircali(dataflow, prmflow, status)

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

% not prepared?
if ~status.pipeline.(status.nodename).prepared
    [dataflow, prmflow, status] = reconnode_aircaliprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% pipeline_onoff
pipeline_onoff = status.currentjob.pipeline_onoff;

% prio
if pipeline_onoff
    % node prio-step
    [dataflow, prmflow, status] = nodepriostep(dataflow, prmflow, status);
    if status.currentjob.topass
        % error or pass
        return;
    end
end

% main
aircaliKernelfunction();

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done

% Kernel funtion
    function aircaliKernelfunction()
        % The anonymous function is static
        debug = [];
        
        nodename = status.nodename;
        nodeprm = prmflow.pipe.(nodename);
        if pipeline_onoff
            nextnode = status.currentjob.nextnode;
            carrynode = status.currentjob.carrynode;
            if isempty(nextnode) || strcmpi(nextnode, 'NULL')
                return;
            end
        end

        if pipeline_onoff
            plconsol = status.currentjob.pipeline;
            index_out = poolindex(dataflow.pipepool.(nextnode), plconsol.Index_out);
            isshotend = plconsol.isshotend;
        else
            isshotend = true;
        end
        
        % parameters to use in prmflow
        Npixel = prmflow.raw.Npixel;
        Nslice = prmflow.raw.Nslice;
        Nfocal = prmflow.raw.Nfocal;
        Nshot = prmflow.raw.Nshot;
        Nviewpershot = prmflow.raw.viewpershot(1);
        Nviewprot = prmflow.raw.Nviewprot;
        Nmulti = Nviewpershot/Nviewprot;
        angulationcode = prmflow.system.angulationcode;
        angulationzero = prmflow.system.angulationzero;     % double
        anglepercode = (pi*2) / double(angulationcode);
        
        % parameters to use
        caliprm = prmflow.pipe.(status.nodename);
        if isfield(caliprm, 'Nsection')
            Nsection = caliprm.Nsection;
        else
            Nsection = 24;
        end
        if isfield(caliprm, 'refpixel')
            refpixel = caliprm.refpixel;
        else
            refpixel = 16;
        end
        if isfield(caliprm, 'refpixelskip')
            refpixelskip = caliprm.refpixelskip;
        else
            refpixelskip = 0;
        end
        if isfield(caliprm, 'corrversion')
            corrversion = caliprm.corrversion;
        else
            % default version is v1.0
            corrversion = 'v1.0';
        end
        if isfield(caliprm, 'referrcutscale')
            referrcutscale = caliprm.referrcutscale;
        else
            referrcutscale = 1.2;
        end
        if isfield(caliprm, 'stabletol')
            stabletol = caliprm.stabletol;
        else
            stabletol = 0.05;
        end
        
        % viewangle and KVmA
        if pipeline_onoff
            viewangle = double(dataflow.pipepool.(carrynode).data.rawhead.AngleEncoder(index_out)).*anglepercode + angulationzero;
            KVmA = dataflow.pipepool.(carrynode).data.rawhead.mA(index_out).*dataflow.pipepool.(carrynode).data.rawhead.KV(index_out);
            % check if the raw is stable
            if stablecheck(reshape(dataflow.pipepool.(carrynode).data.rawdata(:,index_out), Npixel, Nslice, []), stabletol)
                status.errormsg = 'The air calibration is banned due to the unstable air data! please redo the data scan.';
            %     warning(status.errormsg);
                status.jobdone = true;
                status.errorcode = 2;
                return;
            end
        else
            viewangle = double(dataflow.rawhead.AngleEncoder).*anglepercode + angulationzero;
            KVmA = dataflow.rawhead.mA.*dataflow.rawhead.KV;
            % check if the raw is stable
            if stablecheck(reshape(dataflow.rawdata, Npixel, Nslice, []), stabletol)
                status.errormsg = 'The air calibration is banned due to the unstable air data! please redo the data scan.';
            %     warning(status.errormsg);
                status.jobdone = true;
                status.errorcode = 2;
                return;
            end
        end
        

        % refpixel index
        refpixelindex =  [(1:refpixel) + refpixelskip; (Npixel-refpixel+1:Npixel) - refpixelskip];
        % airmain and reference
        aircorr = caliprmforcorr(prmflow, corrversion);
        if pipeline_onoff
            [aircorr.main, aircorr.referrcut, aircorr.referenceKVmA] = ...
                aircalibration(reshape(dataflow.pipepool.(carrynode).data.rawdata(:,index_out), Npixel, Nslice, []), viewangle, refpixelindex, Nsection, Nfocal, KVmA);
        else
            [aircorr.main, aircorr.referrcut, aircorr.referenceKVmA] = ...
                aircalibration(reshape(dataflow.rawdata, Npixel, Nslice, []), viewangle, refpixelindex, Nsection, Nfocal, KVmA);
        end
        
        
        % create aircorr
        if isshotend
            % paramters for aircorr            
            aircorr.Nslice = Nslice;
            aircorr.Nsection = Nsection;
            aircorr.firstangle = 0;     % useless now
            aircorr.refpixel = refpixel;
            aircorr.refnumber = 2;
            % mainsize
            aircorr.mainsize = length(aircorr.main(:));

            % to scale the referrcut
            aircorr.referrcut = aircorr.referrcut.*referrcutscale.*sqrt(Nmulti*Nshot);
            
            dataflow.calibration.aircorr = aircorr;
        end

        if ~pipeline_onoff
            % jobdone
            status.jobdone = true;
        end
    end
end


function r = stablecheck(rawdata, tol)
rawmean = 2.^mean(reshape(rawdata, [], length(rawdata)));
r = any( abs(rawmean./mean(rawmean) - 1) > tol);

end
