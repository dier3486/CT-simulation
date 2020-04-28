function offcorrprm = offfocalloadkernel(offfocalkernel, protocol)
% find out the coupled offfocal kernel
% a hard code based on an external offfocal-kernel defination xml file

offcorrprm = struct();
if ~isempty(offfocalkernel)
    % cell to list
    offfocalkernel = structcellcurse(offfocalkernel);
    % check KV & bowtie
    isKV = [offfocalkernel.offfocalpara(:).KVp] == protocol.KV;
    isbowtie = strcmpi({offfocalkernel.offfocalpara(:).Bowtietype}, protocol.bowtie);
    % found any?
    index_cp = find(isKV & isbowtie, 1);
    if ~isempty(index_cp)
        % check collimator
        iscollimator = strcmp({offfocalkernel.offfocalpara(index_cp).collimation(:).collimationwidth}, ...
            protocol.collimator);
        index_colli = find(iscollimator, 1);
        if ~isempty(index_colli)
            % found
            offcorrprm = offfocalkernel.offfocalpara(index_cp).collimation(index_colli);
        end
    end
end

end