function system = systemconfigure(system_cfg)
% configure CT system

system = struct();

% path
if isfield(system_cfg, 'path')
    system.path = system_cfg.path;
else
    system.path = struct();
end
% IO format govening
if isfield(system.path, 'IOstandard')
    cfggov = system.path.IOstandard;
else
    cfggov = '';
end

% world 
if isfield(system_cfg, 'world')
    system.world = system_cfg.world;
end

% detector
if isfield(system_cfg, 'detector')
    system.detector = struct();
    % copy
    system.detector.detector_corr = system_cfg.detector;
    % detector frame
    frame_base = system_cfg.detector.frame_base;
    if ~isempty(frame_base)
        % load frame_base
        detector_base = loaddata(frame_base, cfggov);
        system.detector.detector_corr = ...
            structmerge(system.detector.detector_corr, detector_base, true);
    end
    % extra frame
    if isfield(system_cfg.detector, 'frame_extra') && ~isempty(system_cfg.detector.frame_extra)
        % load frame_extra
        detector_extra = loaddata(system_cfg.detector.frame_extra, cfggov);
        system.detector.detector_corr = ...
            structmerge(system.detector.detector_corr, detector_extra, true);
        % I will put the ASG in frame_extra
    end
    % response
    if isfield(system_cfg.detector, 'spectresponse') && ~isempty(system_cfg.detector.spectresponse)
        % load detector response
        detector_resp = loaddata(system_cfg.detector.spectresponse, cfggov);
        system.detector.detector_corr = ...
            structmerge(system.detector.detector_corr, detector_resp, true); 
    end
    % crosstalk
    if isfield(system_cfg.detector, 'crosstalk') && ~isempty(system_cfg.detector.crosstalk)
        % load crosstalk
        detector_crosstalk =  loaddata(system_cfg.detector.crosstalk, cfggov);
        system.detector.detector_corr = ...
            structmerge(system.detector.detector_corr, detector_crosstalk, true);
    end
    
    % Note: the settings in system_cfg.detector will not be overwritten by those loaded data.
end

% source
if isfield(system_cfg, 'source')
    system.source = struct();
    % copy
    system.source.tube_corr = system_cfg.source;
    % tube data
    tube_corr = loaddata(system_cfg.source.tubedata, cfggov);
    system.source.tube_corr = structmerge(system.source.tube_corr, tube_corr, true);
    % fill up values
    if ~isfield(system.source.tube_corr, 'focalposition')
        system.source.tube_corr.focalposition = system.detector.detector_corr.focalposition;
    elseif ischar(system.source.tube_corr.focalposition)
        tmp = loaddata(system.source.tube_corr.focalposition);
        system.source.tube_corr.focalposition = tmp.focalposition;
    end
    % reshape to (n, 3)
    system.source.tube_corr.focalposition = reshape(system.source.tube_corr.focalposition, [], 3);
    % off-focal (campatible with previous setup)
    if ~isfield(system.source.tube_corr, 'offfocal') && isfield(system.source.tube_corr, 'offfocalintensity') ...
            && isfield(system.source.tube_corr, 'offfocalwidth')
        system.source.tube_corr.offfocal.width = system.source.tube_corr.offfocalwidth;
        system.source.tube_corr.offfocal.intensity = system.source.tube_corr.offfocalintensity;
    end
end

% collimation
if isfield(system_cfg, 'collimation')
    system.collimation = struct();
    % bowtie
    if isfield(system_cfg.collimation, 'bowtie')
        % multi-bowtie
        if iscell(system_cfg.collimation.bowtie)
            bowtie = system_cfg.collimation.bowtie;
        else
            bowtie = num2cell(system_cfg.collimation.bowtie);
        end
        Nbowtie = length(bowtie(:));
        system.collimation.bowtie = cell(1, Nbowtie);
        system.collimation.bowtie(:) = {struct()};
        for ii = 1:Nbowtie
            % copy
            system.collimation.bowtie{ii}.bowtie_corr = bowtie{ii};
            % bowtie curve
            bowtie_corr = loaddata(bowtie{ii}.bowtiedata, cfggov);
            system.collimation.bowtie{ii}.bowtie_corr = ...
                structmerge(system.collimation.bowtie{ii}.bowtie_corr, bowtie_corr, true);
        end
    end
    % filter
    if isfield(system_cfg.collimation, 'filter')
%         system.collimation.filter = system_cfg.collimation.filter;
        % multi-filter
        if iscell(system_cfg.collimation.filter)
            system.collimation.filter = system_cfg.collimation.filter;
        else
            system.collimation.filter = num2cell(system_cfg.collimation.filter);
        end
    end
    % blades
    if isfield(system_cfg.collimation, 'blade')
        % TBC
        1;
    end
end

% data collector (DCB)
if isfield(system_cfg, 'datacollector')
    system.datacollector = system_cfg.datacollector;
    % log2 (default false)
    if ~isfield(system.datacollector, 'islog2')
        % is the rawdata in log2 compressed
        system.datacollector.islog2 = false;
    end
    % time norm (default false)
    if ~isfield(system.datacollector, 'istimenorm')
        % is to do the time normiliztion of log2
        system.datacollector.istimenorm = false;
    end
    % log2 prm
    if ~isfield(system.datacollector, 'log2scale')
        system.datacollector.log2scale = 11;
    end
    if ~isfield(system.datacollector, 'log2shift')
        system.datacollector.log2shift = 1024;
    end
    % gain
    if ~isfield(system.datacollector, 'DBBzero')
        % zero point of rawdata, defualt is '0x4000'
        system.datacollector.DBBzero = 16384;
    end
    if ~isfield(system.datacollector, 'DBBgain')
        % electronic gain, almost equal to the intensity measurement of one photon in energy of reference keV (e.g. 60keV).
        system.datacollector.DBBgain = 20;
    end
    if ~isfield(system.datacollector, 'Quantumgain')
        % quantum effectivity of bremsstrahlung and X-ray photoelectric effect, which will determine the quantum noise.
        system.datacollector.Quantumgain = 0.5;
        % the total effectivity will be 0.5*0.01
    end
    % rotor
    if ~isfield(system.datacollector, 'angulationcode')
        % the scale number of angulation coder
        system.datacollector.angulationcode = 69120;
    end
    if ~isfield(system.datacollector, 'angulationzero')
        % the zero position of the angulationcode
        system.datacollector.angulationcode = 0;
    end
    % mover (table/couch/conveyor)
    if ~isfield(system.datacollector, 'movercode')
        % the scale number of mover's coder
        system.datacollector.movercode = 2^16;
    end
    if ~isfield(system.datacollector, 'moverlength')
        % the length (in mm) of that movercode covered
        system.datacollector.moverlength = 2000;
    end
    if ~isfield(system.datacollector, 'moveruppersample')
        % upper-sampling of the movercode (to write in rawdata)
        system.datacollector.moveruppersample = 2^16;
    end

end

% simulation method
if isfield(system_cfg, 'simulation')
    system.simulation = system_cfg.simulation;
    % memory
    if ~isfield(system.simulation, 'memorylimit')
        system.simulation.memorylimit = 2.0;        % 2.0 = 2G
    end
    % extra artifacts:
    % quantumnoise
    if ~isfield(system.simulation, 'quantumnoise')
        system.simulation.quantumnoise = false;
    end
    % crosstalk
    if ~isfield(system.simulation, 'crosstalk')
        system.simulation.crosstalk = false;
    end
    % offfocal
    if ~isfield(system.simulation, 'offfocal')
        system.simulation.offfocal = false;
    end
    % GPU
    if ~isfield(system.simulation, 'GPUonoff')
        system.simulation.GPUonoff = 0;
    elseif system.simulation.GPUonoff
        system.simulation.GPUinfo = gpuDevice();
    end
end

% scatter parameters
% TBC

% calibration parameters
% TBC

% console
if isfield(system_cfg, 'console')
    system.console = system_cfg.console;
    % protocoltrans
    if ~isfield(system.console, 'protocoltrans')
        if isfield(system.console, 'protocaltrans')
            % It was a bug of spelling mistake
            % rename to fix a spelling mistake
            system.console.protocoltrans = system.console.protocaltrans;
            system.console = rmfield(system.console, 'protocaltrans');
            % some old configure files could maintain this bug.
        else
            system.console.protocoltrans = struct();
        end
    end
    % protocoltrans.collimationexplain
    if isfield(system.console.protocoltrans, 'collimatorexplain') && ...
            ~isempty(system.console.protocoltrans.collimatorexplain)
        if ischar(system.console.protocoltrans.collimatorexplain)
            system.console.protocoltrans.collimatorexplain_file = system.console.protocoltrans.collimatorexplain;
            system.console.protocoltrans.collimatorexplain = readcfgfile(system.console.protocoltrans.collimatorexplain);
            if ~iscell(system.console.protocoltrans.collimatorexplain.collimator)
                system.console.protocoltrans.collimatorexplain.collimator = ...
                    num2cell(system.console.protocoltrans.collimatorexplain.collimator);
            end
        else
            % do nothing
            1;
        end
    end
    % filetagsrule
    if isfield(system.console.protocoltrans, 'filetagsrule') && ...
            ~isempty(system.console.protocoltrans.filetagsrule)
        if ischar(system.console.protocoltrans.filetagsrule)
            system.console.protocoltrans.filetagsrule_file = system.console.protocoltrans.filetagsrule;
            system.console.protocoltrans.filetagsrule = readcfgfile(system.console.protocoltrans.filetagsrule);
        else
            1;
        end
    end
%     % corrcouplerule
%     if isfield(system.console.protocoltrans, 'corrcouplerule')
%         % nothing to do?
%         1;
%     end
end

% output
if isfield(system_cfg, 'output')
    system.output = system_cfg.output;
    % namerule
    if ~isfield(system.output, 'namerule')
        system.output.namerule = '';
    end
    % dicom dictionary
    if isfield(system_cfg, 'dicomdictionary') && ~isempty(system_cfg.dicomdictionary)
        dicomdict('set', system_cfg.dicomdictionary);
    end
end

return
