function SYS = systemprepare(SYS)
% prepare before the simulation

% samplekeV
if ~isfield(SYS.world, 'samplekeV')
    samplekeV = SYS.world.samplekeV_range(1) : SYS.world.samplekeV_step : SYS.world.samplekeV_range(2);
    SYS.world.samplekeV = samplekeV;
else
    samplekeV = SYS.world.samplekeV;
end

% load matreials
SYS = materialconfigure(SYS, samplekeV, SYS.world.elementsdata, SYS.world.materialdata);

end