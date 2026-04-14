function hout = cttool(varargin)

try
    hFig = imcttool(varargin{:});
catch ME
    if (strcmp(ME.identifier,'MATLAB:UndefinedFunction'))
        firstrun;
    else
        rethrow(ME);
    end
end

if (nargout > 0)
    % Only return handle if caller requested it.
    hout = hFig;
end
end


function firstrun()

YorN = input('It will copy files to MATLAB/toolbox/ (don''t try this if u don''t know how to uninstall them) continue? (y/n) ','s');
if strcmpi(YorN, 'y')
    % come on
    1;
else
    disp('abort...');
    return;
end

imtoolpath = fileparts(which('imtool'));
cttoolpath = fileparts(which(mfilename));
if ~isempty(imtoolpath)
    v = ver('MATLAB');
    version = cellfun(@str2num, regexp(v.Version, '\.', 'split'));
    if version(1)<9 || (version(1)==9 && version(2)<=1)
        % old version (<2019b)
        copyfile(fullfile(cttoolpath, 'cttool\+cttool_2016\*'), imtoolpath, 'f');
    elseif version(1) < 24
        % version 2021
        copyfile(fullfile(cttoolpath, 'cttool\+cttool_2021\*'), imtoolpath, 'f');
    else
        % new version (2024)
        copyfile(fullfile(cttoolpath, 'cttool\+cttool\*'), imtoolpath, 'f');
    end
end

warnmsg = sprintf(['First time in running cttool?\n', ...
    'Do this: in MATLAB panel click Preferences -> General -> Toolbox path cache -> Update Toolbox Path Cache\n', ...
    '  or 预设 -> 常规 -> 工具箱路径缓存 -> 更新工具箱路径缓存\n', ...
    'and try again.']);
warning(warnmsg);

end