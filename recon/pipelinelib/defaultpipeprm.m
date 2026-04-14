function pipeline = defaultpipeprm()
% call me in prepare,
% prmflow.pipe.(nodename).pipeline = defaultpipeprm();

% I know it was called in pipelineinitial.m
pipeline = struct();

pipeline.priostep = true;
pipeline.poststep = true;
pipeline.kernellevel = 0;       % 0, 1 (or 2 not yet)
pipeline.viewrely = [0 0];
pipeline.viewrely_out = [0 0];
pipeline.relystrategy = 'NULL';
pipeline.inputminlimit = 1;
pipeline.inputmaxlimit = inf;
pipeline.outputminlimit = 0;
pipeline.outputmaxlimit = inf;
pipeline.viewdelay = 0;
pipeline.viewexpand = 0;
pipeline.viewextra = [0 0];
pipeline.viewrescale = [1 1];
pipeline.viewcommon = nan;   
pipeline.viewcommonfactor = 1;
pipeline.outputmethod = 'NULL';
pipeline.iscarried = false;
pipeline.currcirculte = [];
pipeline.nextcirculte = [];
pipeline.objecttype = [];       % do not set this in recon configure file, user can only set inputobjecttype/nextobjecttype.
pipeline.inputobjecttype = 'NULL';
pipeline.nextobjecttype = 'NULL';
pipeline.currdatasize = nan; 
pipeline.nextdatasize = nan;
pipeline.currpoolsize = nan;     
pipeline.nextpoolsize = nan;
pipeline.keepbottom = [];       % not used, user can not set it in recon configure file
pipeline.nextkeepbottom = nan;  % do not set this in recon configure file, it is remained to the prepare function
pipeline.currgpuDevice = 0;     % uint, 0: CPU, 1: GPU#1, other: not now (will could be a vector)
pipeline.nextgpuDevice = 0;     % uint, 0: CPU, 1: GPU#1, other: not now (will could be a vector)
pipeline.currResource = 0;      % uint, 0: default, 1: Matlab, 2: sharedC (lib.pointer)
pipeline.nextResource = 0;      % uint, 0: default, 1: Matlab, 2: sharedC (lib.pointer)
pipeline.GPUonoff = -1;         % GPUonoff flag, -1: default, 0: close, 1: on.
pipeline.nextdatabasis = 0;     % real/complex flag, 0: default, 1: real, 2: complex, 3: 3-basis, 4: quaternion.
                                % we don't have 'currdatabasis' because user can not set it
pipeline.CUDAonoff = false;     % call cuda dll, very limited

% all those nan and inf are special int.
end