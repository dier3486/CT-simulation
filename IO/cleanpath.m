function configure = cleanpath(configure, mainpath, escapepath)
% clean the pathes

if nargin < 3
    escapepath = 'system.path';
end

% to replace relative pathes ~/ and ~\ by mainpath
configure = structregexprep(configure, '~(/|\\)', regexptranslate('escape', mainpath));
% to replace relative pathes $pathname$ by configure.system.path.(pathname)
configure = structregexprep(configure, '\$(\w*)(\\|/)', ['${regexptranslate(''escape'', root.' escapepath '.($1))}']);
% kown bug: \. -> . , \- -> -
configure = structregexprep(configure, '\\([\-,\.])', '$1');
% kown bug: \ -> \\
configure = structregexprep(configure, '\\+', '\\');

end