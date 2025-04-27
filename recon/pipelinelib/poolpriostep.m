function [currpool, nextpool, jobstatus] = poolpriostep(currpool, nextpool, pipeprm, jobstatus)
% real prio-step of a pipeline node

% exception
if isempty(currpool)
    % need not input
    return;
    % we will analysis these case later
end

% for multi pools, only employ the 1st
if length(currpool) > 1 || length(nextpool) > 1
    [currpool(1), nextpool(1), jobstatus] = poolpriostep(currpool(1), nextpool(1), pipeprm, jobstatus);
    return;
end

% make them easy
% nodetype
if ~strcmpi(pipeprm.nodetype, 'NULL')
    nodetype = pipeprm.nodetype(3:end);
else
    nodetype = 'NULL';
end
% I know the pipeprm.nodetype is [HA]-[HA].[01].[NGS] or NULL
% and
Lcurr = pipeprm.viewrely;
L = pipeprm.viewrely_out;
E = pipeprm.viewextra;
Q = pipeprm.viewrescale;
Nexp = pipeprm.viewexpand;
% only use in this function

% pass due to nextpool is stucked
if ~isempty(nextpool) && nextpool.WriteStuck
    jobstatus.readnumber = 0;
    jobstatus.writenumber = 0;
    jobstatus.newAvail = 0;
    jobstatus.jobdone = 6;
    % done and keep waking
    return;
    % while a node without output data, we should never lock the nextpool.
end

% check A-H
if (~isempty(nextpool) && ~nextpool.circulatemode) && currpool.circulatemode
    % for A-H.*, specially setting the E to -L
    E = -L;
    flag_AH = true;
else
    flag_AH = false;
end

% check shot start 1
if ~isempty(nextpool)
    jobstatus.isshotstart = nextpool.isshotstart;
    % close it
    nextpool.isshotstart = false;
else
    jobstatus.isshotstart = currpool.ReadPoint == currpool.ReadStart;
end
% jobstatus.isshotstart = currpool.ReadPoint == currpool.ReadStart;
% we use the ReadPoint/ReadStart to judge the shot start
if jobstatus.isshotstart
    % #0
    % views of the current shot
    Ns_curr = currpool.ReadEnd - currpool.ReadStart + 1;
    Ns_in = floor(Ns_curr * Q(1) / Q(2));
    % Note: I know the / is before * that means the remainder views will be abandoned.
    if flag_AH
        % for A-H
        Ns_in = Ns_in + L(1) + L(2);
        % I know the A-H.* are equivalent to E = -L but Ns_in += L(1) + L(2).
    end
    % the planed output view number (of this shot)
    Ns_out = Ns_in + E(1) + E(2);

    % check poolsize
    if ~isempty(nextpool) && nextpool.circulatemode
        % type A.*
        if nextpool.poolsize ~= Ns_out
            % should be some buffer re-malloc here.
            nextpool.poolsize = Ns_out;
        end
    end
    % Here is a check point for C++ codes to new or re-malloc the buffers.

    % #1
    % initial nextpool
    if ~isempty(nextpool)
        % Wstart, Rstart
        if nextpool.circulatemode
            % type A.*
            if pipeprm.iscarried
                nextpool.WriteStart = currpool.ReadStart;
            elseif ~isavail(nextpool.WriteStart)
                nextpool.WriteStart = 1;
            end
            switch nodetype
                case {'A.0.N', 'A.0.G', 'A.0.S'}
                    nextpool.ReadStart = nextpool.WriteStart + L(1);
                case {'A.1.N', 'A.1.G'}
                    nextpool.ReadStart = nextpool.WriteStart + L(1) + L(2);
                case 'A.1.S'
                    nextpool.ReadStart = nextpool.WriteStart;
                case 'NULL'
                    0;
                otherwise
                    nodetypeerror(nodetype);
                    return;
            end
            % period
            nextpool.WriteStart = mod(nextpool.WriteStart - 1, nextpool.poolsize) + 1;
            nextpool.ReadStart = mod(nextpool.ReadStart - 1, nextpool.poolsize) + 1;
        else
            % type H.*
            switch nodetype
                case {'H.0.N', 'H.0.G', 'H.0.S'}
                    if pipeprm.iscarried
                        nextpool.WriteStart = currpool.ReadStart;
                        if nextpool.WriteStart < max(1, 1 + E(1))
                            carryerror();
                            return;
                        end
                    else
                        nextpool.WriteStart = max(1, 1 + E(1));
                    end
                    nextpool.ReadStart = nextpool.WriteStart - E(1);
                case {'H.1.N', 'H.1.G'}
                    nextpool.WriteStart = 1 + E(1) - L(2);
                    nextpool.ReadStart = 1;
                case 'H.1.S'
                    nextpool.WriteStart = 1 + E(1) + L(1);
                    nextpool.ReadStart = 1;
                case 'NULL'
                    0;
                otherwise
                    nodetypeerror(nodetype);
                    return;
            end
        end

        % Wend, Rend
        if nextpool.circulatemode
            % type A.*
            nextpool.WriteEnd = nextpool.WriteStart + Ns_in - 1;
            nextpool.ReadEnd = nextpool.ReadStart + Ns_out - 1;
        else
            switch nodetype
                case {'H.0.N', 'H.0.G', 'H.0.S', 'H.1.N', 'H.1.G'}
                    nextpool.WriteEnd = nextpool.WriteStart + Ns_in - 1;
                case 'H.1.S'
                    nextpool.WriteEnd = nextpool.WriteStart + Ns_in - 1 - L(1) - L(2);
                case 'NULL'
                    0;
                otherwise
                    nodetypeerror(nodetype);
                    return;
            end
            nextpool.ReadEnd = nextpool.ReadStart + Ns_out - 1;
        end

        % W0, R0
        nextpool.WritePoint = nextpool.WriteStart;
        nextpool.ReadPoint = nextpool.ReadStart;

        % PA0
        nextpool.AvailPoint = nextpool.ReadPoint - 1;
        switch nodetype
            case {'H.0.N', 'H.0.G', 'H.0.S', 'H.1.N', 'H.1.G'}
                nextpool.AvailPoint = nextpool.AvailPoint + E(1) - L(2);
            case 'H.1.S'
                nextpool.AvailPoint = nextpool.AvailPoint + E(1) + L(1);
            case {'A.0.N', 'A.0.G', 'A.0.S', 'A.1.N', 'A.1.G'}
                nextpool.AvailPoint = nextpool.AvailPoint - L(1) - L(2);
            case 'A.1.S'
                1;
            case 'NULL'
                0;
            otherwise
                nodetypeerror(nodetype);
                return;
        end
    end
    % ini Index_in/Index_out
    jobstatus.Index_in = [0, 0];
    jobstatus.Index_out = [0, 0];
end

% #2
% minlimit
switch nodetype
    case {'H.0.N', 'H.0.G', 'H.0.S', 'H.1.N', 'H.1.G', 'A.0.N', 'A.0.G', 'A.0.S', 'A.1.N', 'A.1.G'}
        jobstatus.minlimit = pipeprm.inputminlimit;
    case {'H.1.S', 'A.1.S'}
        jobstatus.minlimit = pipeprm.inputminlimit + Lcurr(1) + Lcurr(2);
    case 'NULL'
        jobstatus.minlimit = pipeprm.inputminlimit;
    otherwise
        nodetypeerror(nodetype);
        return;
end
% The minlimit is not less than the prepared pipeprm.inputminlimit. We suggest (while no speical requirements) to prepare the
% inputminlimit=1 for normal nodes and inputminlimit=0 for the nodes need not inputs.

% Index_in.L
if pipeprm.kernellevel==0
    if ~isempty(nextpool)
        jobstatus.Index_in(1) = nextpool.WritePoint;
    else
        jobstatus.Index_in(1) = currpool.ReadPoint;
    end
elseif pipeprm.kernellevel==1
    jobstatus.Index_in(1) = currpool.ReadPoint;
else
    error('%d', pipeprm.kernellevel);
end

% Intdex_out.L
if ~isempty(nextpool)
    AvailNumber = nextpool.AvailPoint - nextpool.ReadPoint + 1;
    switch nodetype
        case {'H.0.N', 'H.0.G', 'H.0.S', 'A.0.S'}
            jobstatus.Index_out(1) = max(nextpool.WritePoint - L(2), nextpool.WritePoint - L(2) - AvailNumber);
            % in starting of a shot the AvailNumber could be negative
        case {'H.1.N', 'H.1.G'}
            jobstatus.Index_out(1) = max(nextpool.WritePoint, nextpool.WritePoint - AvailNumber);
        case {'H.1.S', 'A.1.N', 'A.1.G', 'A.1.S'}
            jobstatus.Index_out(1) = nextpool.WritePoint;
        case {'A.0.N', 'A.0.G'}
            jobstatus.Index_out(1) = nextpool.WritePoint - L(2);
        case 'NULL'
            0;
        otherwise
            nodetypeerror(nodetype);
            return;
    end
end

% D of out to in
if ~isempty(nextpool)
    switch nodetype
        case {'H.0.N', 'H.0.G', 'H.1.N', 'H.1.G', 'H.0.S', 'A.0.S'}
            jobstatus.Do2i = max(-L(2), -L(2)-AvailNumber);
        case {'H.1.S', 'A.1.S'}
            jobstatus.Do2i = L(1);
        case {'A.0.N', 'A.0.G', 'A.1.N', 'A.1.G'}
            jobstatus.Do2i = -L(2);
        case 'NULL'
            jobstatus.Do2i = nan;
        otherwise
            nodetypeerror(nodetype);
            return;
    end
end

% Nexpand
if isempty(nextpool) || nextpool.circulatemode
    % type A.*
    jobstatus.Nexpand = 0;
else
    switch nodetype
        case {'H.0.N', 'H.0.G'}
            jobstatus.Nexpand = L(1) * (nextpool.WriteEnd + E(2) > nextpool.poolsize);
        case {'H.1.N', 'H.1.G'}
            jobstatus.Nexpand = (L(1)+L(2)) * (nextpool.WriteEnd + L(2) + E(2) > nextpool.poolsize);
        case 'H.0.S'
            jobstatus.Nexpand = 0;
        case 'H.1.S'
            jobstatus.Nexpand = -L(1) - L(2);
%             if nextpool.WriteEnd > nextpool.poolsize
%                 jobstatus.Nexpand = -L(1) - L(2);
%             else
%                 jobstatus.Nexpand = E(2) - L(1);
%             end
        case 'NULL'
            0;
        otherwise
            nodetypeerror(nodetype);
            return;
    end
end

% #3
% check shot start 2
if jobstatus.isshotstart
    % fix minlimit
    switch nodetype
        case {'H.0.N', 'H.0.G', 'H.1.N', 'H.1.G'}
            E1 = ceil(E(1) * Q(2)/Q(1));
            jobstatus.minlimit = max(jobstatus.minlimit, 1 - E1);
        case 'H.0.S'
            E1 = ceil(E(1) * Q(2)/Q(1));
            jobstatus.minlimit = max(jobstatus.minlimit, 1 - E1 + Lcurr(2));
        case {'H.1.S', 'A.0.N', 'A.0.G', 'A.0.S', 'A.1.N', 'A.1.G', 'A.1.S'}
            1;
        case 'NULL'
            0;
        otherwise
            nodetypeerror(nodetype);
            return;
    end

    % Intdex_out.L, D
    if strcmpi(nodetype, 'H.1.S') && ~isempty(nextpool)
        jobstatus.Index_out(1) = 1;
        jobstatus.Do2i = -E(1);
    end
end

% #4
% check shot end 1 (the input data has reached the end)
Ishotend1 = currpool.WritePoint == currpool.WriteEnd + 1;
if Ishotend1
    % special Nexpand
    if strcmpi(nodetype, 'H.0.S') && ~currpool.circulatemode
        % H-H.0.S
        jobstatus.Nexpand = max(0, E(2));
    elseif strcmpi(nodetype, 'H.1.S')
        % H.1.S
        jobstatus.Nexpand = E(2) - L(1);
    end
    % special curr AvailPoint, to move after the ReadEnd
    if ~isempty(nextpool)
        if (currpool.circulatemode && ~nextpool.circulatemode) || strcmpi(nodetype, 'A.1.S')
            % A-H or A.1.S
            currpool.AvailPoint = currpool.ReadEnd + Lcurr(1) + Lcurr(2);
            %             AvailNumber = AvailNumber + Lcurr(1) + Lcurr(2);
        end
    end
end
% add set Nexp to Nexpand
if ~Ishotend1
    jobstatus.Nexpand = jobstatus.Nexpand + Nexp;
end

% #5
% n, m
currAvailNumber = currpool.AvailPoint - currpool.ReadPoint + 1;
checkNumber = min(pipeprm.inputmaxlimit, currAvailNumber);
[m, n] = checkreadingnum(nextpool, checkNumber, pipeprm.viewrescale, pipeprm.viewcommon, jobstatus.Nexpand, Ishotend1);
% m = n*Q(1)/Q(2);
jobstatus.Ninres = m;

% pass due to not enough input data
if n < jobstatus.minlimit && ~Ishotend1
    jobstatus.readnumber = 0;
    jobstatus.writenumber = 0;
    jobstatus.newAvail = 0;
    currAvailNumber_common = floor(currAvailNumber/pipeprm.viewcommon)*pipeprm.viewcommon;
    if currAvailNumber_common < jobstatus.minlimit    
        jobstatus.jobdone = 3;
        % nothing to read, shall pass
    else
        % seems nextpool is stucked
        jobstatus.jobdone = 6;
    end
    % withdraw the initial of next pool
    if jobstatus.isshotstart && ~isempty(nextpool)
        nextpool.isshotstart = true;
    end
    return;
    % while a node needs not input data, we should set the jobstatus.minlimit=0 to avoid the pass.
end

% check error about n
if ~isempty(nextpool)
    W2WE = nextpool.WriteEnd - nextpool.WritePoint + 1;
    switch nodetype
        case {'H.0.N', 'H.0.G', 'H.0.S', 'H.1.N', 'H.1.G', 'A.0.N', 'A.0.G', 'A.0.S', 'A.1.N', 'A.1.G'}
            if W2WE < m
                datarangeerror();
                return;
            end
        case {'H.1.S', 'A.1.S'}
            if W2WE < m - L(1) - L(2)
                datarangeerror();
                return;
            end
        case 'NULL'
            0;
        otherwise
            nodetypeerror(nodetype);
            return;
    end
end

% readnumber, writenumber and newAvail
switch nodetype
    case {'H.0.N', 'H.0.G', 'H.0.S', 'H.1.N', 'H.1.G', 'A.0.N', 'A.0.G', 'A.0.S', 'A.1.N', 'A.1.G'}
        jobstatus.readnumber = n;
        jobstatus.writenumber = m;
        jobstatus.newAvail = m;
    case {'H.1.S', 'A.1.S'}
        jobstatus.readnumber = n - Lcurr(1) - Lcurr(2);
        jobstatus.writenumber = m - L(1) - L(2);
        jobstatus.newAvail = m - L(1) - L(2);
    case 'NULL'
        jobstatus.readnumber = n;
        jobstatus.writenumber = 0;
        jobstatus.newAvail = 0;
    otherwise
        nodetypeerror(nodetype);
        return;
end

% #6
% Index_in.U
jobstatus.Index_in(2) = jobstatus.Index_in(1) + n - 1;

% Index_out.U
if ~isempty(nextpool)
    switch nodetype
        case {'H.0.N', 'H.0.G'}
            jobstatus.Index_out(2) = min(nextpool.WritePoint + m + L(1) - 1, nextpool.WriteEnd + E(2));
        case {'H.1.N', 'H.1.G'}
            jobstatus.Index_out(2) = min(nextpool.WritePoint + m + L(1) + L(2) - 1, nextpool.WriteEnd + L(2) + E(2));
        case {'H.0.S', 'A.0.S'}
            jobstatus.Index_out(2) = nextpool.WritePoint + m - L(2) - 1;
        case {'H.1.S', 'A.1.S'}
            jobstatus.Index_out(2) = nextpool.WritePoint + m - L(1) - L(2) - 1;
        case {'A.0.N', 'A.0.G'}
            jobstatus.Index_out(2) = nextpool.WritePoint + m + L(1) - 1;
        case {'A.1.N', 'A.1.G'}
            jobstatus.Index_out(2) = nextpool.WritePoint + m + L(1) + L(2) - 1;
        case 'NULL'
            0;
        otherwise
            nodetypeerror(nodetype);
            return;
    end
end

% #7
% check shot end 2 (will output the shot end)
if ~isempty(nextpool)
    jobstatus.isshotend = (W2WE == jobstatus.writenumber) && Ishotend1;
else
    jobstatus.isshotend = currpool.WritePoint == currpool.WriteEnd;
end

if ~isempty(nextpool) && jobstatus.isshotend
    % I know the W2WE is distance of WritePoint to WriteEnd (+1)
    % newAvail
    %     if ~nextpool.circulatemode
    %         % H.*
    %         jobstatus.newAvail = jobstatus.newAvail + L(2) + E(2);
    %     else
    switch nodetype
        case {'H.0.N', 'H.0.G', 'H.0.S', 'H.1.N', 'H.1.G'}
            jobstatus.newAvail = jobstatus.newAvail + L(2) + E(2);
        case 'H.1.S'
            jobstatus.newAvail = jobstatus.newAvail + L(2) + E(2);
        case {'A.0.N', 'A.0.G', 'A.0.S', 'A.1.N', 'A.1.G'}
            jobstatus.newAvail = jobstatus.newAvail + L(1) + L(2);
        case 'A.1.S'
            2;
        case 'NULL'
            0;
        otherwise
            nodetypeerror(nodetype);
            return;
    end
    % Index_out.U
    if pipeprm.relystrategy==2
        % All S
        switch nodetype
            case 'H.0.S'
                jobstatus.Index_out(2) = nextpool.WriteEnd + E(2);
            case 'H.1.S'
                jobstatus.Index_out(2) = nextpool.WriteEnd + E(2) + L(2);
            case 'A.0.S'
                jobstatus.Index_out(2) = nextpool.WriteEnd + L(1);
            case 'A.1.S'
                jobstatus.Index_out(2) = nextpool.WriteEnd;
            case 'NULL'
                0;
            otherwise
                nodetypeerror(nodetype);
                return;
        end
        % we don't use the nextpool.ReadEnd that means to move the ReadEnd is safety.
    end
    % readnumber
    if strcmpi(nodetype, 'H.1.S') || strcmpi(nodetype, 'A.1.S')
        jobstatus.readnumber = n;
    end
    if n~=currAvailNumber
        jobstatus.readnumber = currAvailNumber;
    end
end

% #99
% % check error
% if ~isempty(nextpool)
%     if W2WE < jobstatus.writenumber
%         error('WritePoint error!');
%     end
%     if jobstatus.readnumber > currAvailNumber
%         error('readnumber error!')
%     end
% end

% #done
if jobstatus.readnumber == currAvailNumber
    % all done
    jobstatus.jobdone = 1;
    % Note: while new available view <= 0 the jobdone should be 4,
    if ~isempty(nextpool) && (nextpool.AvailPoint + jobstatus.newAvail < nextpool.ReadPoint)
        % technically passed
        jobstatus.jobdone = 4;
    end
else
    % partly done
    jobstatus.jobdone = 2;
    % Note: while new available view <= 0 the jobdone should be 7,
    if ~isempty(nextpool) && (nextpool.AvailPoint + jobstatus.newAvail < nextpool.ReadPoint)
        % call me again
        jobstatus.jobdone = 7;
    end
end

    function nodetypeerror(nodetype)
        jobstatus.jobdone = 0;
        jobstatus.errorcode = -309;
        jobstatus.errormsg = sprintf('illeagal nodetype %s.', nodetype);
        return;
    end

    function datarangeerror()
        jobstatus.jobdone = 0;
        jobstatus.errorcode = -302;
        jobstatus.errormsg = sprintf('outof writable range.');
        return;
    end

    function carryerror()
        jobstatus.jobdone = 0;
        jobstatus.errorcode = -303;
        jobstatus.errormsg = sprintf('can not carry the node due to viewextra.');
        return;
    end
end


function [m, n] = checkreadingnum(nextpool, Availnum, viewrescale, viewcommon, viewexpand, isshotend)
% to check how many data to read in pipeline node's prio step

% if nargin<2 || isempty(Availnum)
%     Availnum = inf;
% end

if nargin<3 || isempty(viewrescale)
    viewrescale = [1 1];
end

if nargin<4 || isempty(viewcommon)
    viewcommon = viewrescale(2);
end

if nargin<5 || isempty(viewexpand)
    viewexpand = 0;
end

if nargin<6 || isempty(isshotend)
    isshotend = false;
end

Availnum = max(0, Availnum);
% qc = lcm(viewrescale(2), viewcommon);
qc = viewcommon;

if isempty(nextpool) || nextpool.circulatemode || viewrescale(1) == 0
    % a null pool or a circulate pool
    if ~isshotend
        n = floor(Availnum/qc)*qc;
        m = n*viewrescale(1)/viewrescale(2);
    else
        n = Availnum;
        m = floor(Availnum*viewrescale(1)/viewrescale(2));
    end
    % A null pool have infinite space, a circulate pool is also.
elseif isfield(nextpool, 'WriteStuck') && nextpool.WriteStuck
    % writing locked
    m = 0;
    n = 0;
else
    % non-circulate
    d = nextpool.poolsize - nextpool.WritePoint + 1 - viewexpand;
    if ~isshotend
        m1 = floor(Availnum/qc)*qc*viewrescale(1)/viewrescale(2);
        m2 = floor(d/viewrescale(1))*viewrescale(1);
        m = min(m1, m2);
        n = m*viewrescale(2)/viewrescale(1);
    else
        m3 = floor(Availnum*viewrescale(1)/viewrescale(2));
        if m3 < d
            n = Availnum;
            m = m3;
        else
            m1 = floor(Availnum/qc)*qc*viewrescale(1)/viewrescale(2);
            m2 = floor(d/viewrescale(1))*viewrescale(1);
            m = min(m1, m2);
            n = m*viewrescale(2)/viewrescale(1);
        end
    end
end

end
