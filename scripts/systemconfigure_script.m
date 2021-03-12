cfgfile = '.\system\mod\sample_configure.xml';
configure = readcfgfile(cfgfile);
configure = configureclean(configure);
SYS = systemconfigure(configure.system);
if isfield(configure, 'phantom')
    SYS.phantom = phantomconfigure(configure.phantom);
end
SYS = systemprepare(SYS);