function configure = configureclean(configure)
% clean the configure

% relative pathes
if isfield(configure, 'system')
    % main path
    if isfield(configure.system.path, 'main')
        mainpath = configure.system.path.main;
    else
        mainpath = '';
    end
    if isempty(mainpath)
        % set current path as mainpath
        mainpath = [pwd() '/'];
        configure.system.path.main = mainpath;
    end
    % replace the pathes
    configure = cleanpath(configure, mainpath);
end

% sub configure
if isfield(configure, 'system')
    fields = fieldnames(configure.system);
    for ii = 1:length(fields)
        if ischar(configure.system)
        end
    end
end

% fillup protocal
if isfield(configure, 'protocal')
    if ~isfield(configure.protocal, 'seriesnumber')
        configure.protocal.seriesnumber = 1;
    end
    if ~iscell(configure.protocal.series)
        % to cell
        configure.protocal.series = {configure.protocal.series};
    end
    if configure.protocal.seriesnumber>1
        for ii = 2:configure.protocal.seriesnumber
            configure.protocal.series{ii} = ...
                structmerge(configure.protocal.series{ii}, configure.protocal.series{ii-1});
        end
    end
end

end


function configure = cleanpath(configure, mainpath)
% clean the pathes

% to replace relative pathes ~/ and ~\ by mainpath
configure.system = structregexprep(configure.system, '~(/|\\)', regexptranslate('escape', mainpath));
% to replace relative pathes $pathname$ by configure.system.path.(pathname)
configure = structregexprep(configure, '\$(\w*)(\\|/)', '${regexptranslate(''escape'', root.system.path.($1))}');
% kown bug: \ -> \\
configure = structregexprep(configure, '\\+', '\\');
end


function cfg = reload(cfg)
% load sub-configure files
fields = fieldnames(configure.system);
    for ii = 1:length(fields)
        if ischar(configure.system)
        end
    end

end

