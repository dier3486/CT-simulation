function [jobdone, errorcode] = jobdonetable(prio1, prio2, main, post)

if main == 0
    % main function returned error
    jobdone = 0;
    errorcode = '0';
    return
end

errorcode = [];

if prio1==2
    if post==2
        % typical partly output
        jobdone = 2;
        return
    else
        % error, output pool error
        jobdone = 0;
        errorcode = sprintf('%d%d%d', post, prio1, prio2);
        return
    end
end

if prio2==2
    % partly input
    job_mp = 2;
else
    job_mp = main;
end

switch job_mp
    case {1, 2, 5}
        if post==2
            jobdone = 2;
        else
            jobdone = job_mp;
        end
    case {3, 4}
        if post~=0
            jobdone = 0;
            errorcode = sprintf('%d%d%d', post, prio1, main);
        else
            if prio1==0
                jobdone = job_mp;
            else % prio1==1
                jobdone = 1;
            end
        end
    otherwise
        jobdone = 0;
        errorcode = sprintf('%d', main);
end

end