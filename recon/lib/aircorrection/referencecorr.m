function [ref, reflast] = referencecorr(rawref, refblock, ref0)
% reference correction

if nargin < 3
    ref0 = [];
end

% idx_both = ~blk1 & ~blk2;
idx_both = ~refblock(1, :) & ~refblock(2, :);
% idx_both are the views not blocked

% refblk = [blk1; blk2];

Nview = size(rawref, 2);
rawref = reshape(rawref, [], 2, Nview);
Nref = size(rawref, 1);
ref = zeros(Nref, Nview, 'single');
ref(:, idx_both) = squeeze((rawref(:, 1, idx_both) + rawref(:, 2, idx_both))./2);

if any(~idx_both)
    % some views are blocked
    blkidx_1 = ~refblock(1, :) & refblock(2, :);
    blkidx_2 = refblock(1, :) & ~refblock(2, :);
    %
    view_bkl1 = find(~idx_both, 1, 'first');
    view0 = find(idx_both, 1, 'first');
    if isempty(view0) && isempty(ref0)
        % all views are blocked ??
        ref = zeros(Nref, Nview, 'single');
        reflast = ref0;
        return;
    end
    if view_bkl1==1
        % the first view is blocked
        if ~isempty(ref0)
            % ref on '0'-th view is inputed
            if blkidx_1(1)
                % ref1
                ref(:, 1) = ref0.ref-squeeze(ref0.rawref(:, 1)-rawref(:, 1, 1));
            elseif blkidx_2(1)
                % ref2
                ref(:, 1) = ref0.ref-squeeze(ref0.rawref(:, 2)-rawref(:, 2, 1));
            else
                % mA
                ref(:, 1) = ref0.ref;
            end
            % move on
            view_bkl1 = 2;
        else
            %1 go back
            for ii = view0-1:-1:1
                if blkidx_1(ii)
                    % ref1
                    ref(:, ii) = ref(:, ii+1)-squeeze(rawref(:, 1, ii+1)-rawref(:, 1, ii));
                elseif blkidx_2(ii)
                    % ref2
                    ref(:, ii) = ref(:, ii+1)-squeeze(rawref(:, 2, ii+1)-rawref(:, 2, ii));
                else
                    % mA
                    ref(:, ii) = ref(:, ii+1);
                end
            end
            %2 forward
            view_bkl1 = find(~idx_both(view0:end), 1, 'first');
        end
    end
    for ii = view_bkl1:Nview
        if idx_both(ii)
            continue
        end
        if blkidx_1(ii)
            % ref1
            ref(:, ii) = ref(:, ii-1)-squeeze(rawref(:, 1, ii-1)-rawref(:, 1, ii));
        elseif blkidx_2(ii)
            % ref2
            ref(:, ii) = ref(:, ii-1)-squeeze(rawref(:, 2, ii-1)-rawref(:, 2, ii));
        else
            % mA
            ref(:, ii) = ref(:, ii-1);
        end
    end
end

% last view
reflast = struct();
reflast.ref = ref(:, end);
reflast.rawref = rawref(:, :, end);

end