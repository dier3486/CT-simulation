function configure = configureclean(configure)
% clean the configure
% configure = configureclean(configure)

% load sub confgirue files
configure = reload(configure, true);

% sub configure
if isfield(configure, 'system')
    configure.system = reload(configure.system);
end

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
    if mainpath(end)~='/' && mainpath(end)~='\'
        mainpath = [mainpath '/'];
    end
    % replace the pathes
    configure = cleanpath(configure, mainpath);
    % double check
    pathes = fieldnames(configure.system.path);
    for ipth = 1:length(pathes)
        path_i = configure.system.path.(pathes{ipth});
        if ~exist(path_i, 'dir')
            errmsg = ['Not exist dir ''' path_i ''' in configure.system.path.' pathes{ipth}];
            error(errmsg);
        end
    end
end

% fillup protocol
if isfield(configure, 'protocol')
    if ~isfield(configure.protocol, 'seriesnumber')
        configure.protocol.seriesnumber = 1;
    end
    if ~iscell(configure.protocol.series)
        % to cell
        configure.protocol.series = {configure.protocol.series};
    end
    if configure.protocol.seriesnumber>1
        for ii = 2:configure.protocol.seriesnumber
            configure.protocol.series{ii} = ...
                structmerge(configure.protocol.series{ii}, configure.protocol.series{ii-1});
        end
    end
end

end


% function configure = cleanpath(configure, mainpath)
% % clean the pathes
% 
% % to replace relative pathes ~/ and ~\ by mainpath
% configure.system = structregexprep(configure.system, '~(/|\\)', regexptranslate('escape', mainpath));
% % to replace relative pathes $pathname$ by configure.system.path.(pathname)
% configure = structregexprep(configure, '\$(\w*)(\\|/)', '${regexptranslate(''escape'', root.system.path.($1))}');
% % kown bug: \. -> . , \- -> -
% configure = structregexprep(configure, '\\([\-,\.])', '$1');
% % kown bug: \ -> \\
% configure = structregexprep(configure, '\\+', '\\');
% end


function cfg = reload(cfg, catcherror)
% load sub-configure files

if nargin < 2
    catcherror = false;
end

fields = fieldnames(cfg);
for ii = 1:length(fields)
    % loop the fields to find out the files like 'xxx.xml' to read
    field_ii = cfg.(fields{ii});
    if ~ischar(field_ii)
        continue
    end
    % delete blank
    field_ii = strtrim(field_ii);
    if isempty(field_ii)
        continue
    end
    % try to read
    if exist(field_ii, 'file')
        [~, ~, fileext] = fileparts(field_ii);
        % only xml and json files are supported
        if strcmp(fileext, '.xml') || strcmp(fileext, '.json')
            % a sub configure file
            cfg.(fields{ii}) = readcfgfile(field_ii);
        end
        % NOTE: do not use .mat or bin file as a configure file
    elseif catcherror
        error(['Can not find ' field_ii]);
    end
end
% NOTE: I don't think a recursed configure file structure is a good idea.
end

