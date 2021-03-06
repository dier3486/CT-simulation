function [u, w] = cylinderconecut(u, object, cone, isinverse)
% mapping of a cylinder cut by collimator cones (and its inverse)

if nargin<4
    isinverse = false;
end

if ~strcmpi(object.type, 'cylinder')
    error('Only Cylinder!');
end

Np = size(u, 1);
u1 = u;
r1 = sqrt(sum(u1(:,[1 2]).^2, 2));
s_r1 = r1>1;
u1(s_r1, [1 2]) = u1(s_r1, [1 2])./r1(s_r1);
v0 = u1*object.vector+object.O;

% Z direction
zsign = sign(object.vector(3,3));
% if zsign<0
%     cone = cone([2 1]);
% end

% tcone
tcone = zeros(Np, 2);
for ic = 1:2
    u_ic = (v0-cone(ic).O)/cone(ic).vector;
    nz1 = [0 0 1]*object.vector/cone(ic).vector;
    
    a = nz1(1)^2+nz1(2)^2-nz1(3)^2;
    b = u_ic*[nz1(1); nz1(2); -nz1(3)];
    c = u_ic(:,1).^2 + u_ic(:,2).^2 - u_ic(:,3).^2;
    d = real(sqrt(b.^2-a.*c));
    t1 = (-d-b)./a;
    t2 = (d-b)./a;
    k1 = u_ic + t1*nz1;
    k2 = u_ic + t2*nz1;
    s1 = k1(:, 3)>=0 & k1(:, 2)>=0 & d>0;
    s2 = k2(:, 3)>=0 & k2(:, 2)>=0 & d>0;
    t1(~s1) = nan;
    t2(~s2) = nan;
    tcone(:, ic) = min(t1, t2);
%     s_rng = tcone(:, ic)>0;
    s_rng = tcone(:, ic).*zsign>0;
    s_nan = isnan(tcone(:, ic));
    if ic==1
        s_zsign = (a>=0) & (nz1(3)>0);
    end
    
    if ic==1
        s_inc = (c<=0 | u_ic(:,2)<0) & u_ic(:, 3)>0;
%         t1(~s1 & ~s_inc) = inf;
%         t1(~s1 & s_inc) = -inf;
%         t2(~s2 & ~s_inc) = inf;
%         t2(~s2 & s_inc) = -inf;
%         tcone(:, ic) = min(t1, t2);
%         tcone(~s_inc, ic) = min(t1(~s_inc), t2(~s_inc));
%         tcone(s_inc, ic) = max(t1(s_inc), t2(s_inc));
        

        tcone((s_nan | s_rng) & s_inc, ic) = -inf*(s_zsign-0.5);
        tcone((s_nan | ~s_rng) & ~s_inc, ic) = inf*(s_zsign-0.5);

        
        
    else
        % ic==2
        s_inc = (c>0 | u_ic(:, 3)<0) & u_ic(:,2)>0;
        tcone((s_nan | ~s_rng) & s_inc, ic) = inf*(s_zsign-0.5);
        tcone((s_nan | s_rng) & ~s_inc, ic) = -inf*(s_zsign-0.5);
        
%         t1(~s1 & ~s_inc) = -inf;
%         t1(~s1 & s_inc) = inf;
%         t2(~s2 & ~s_inc) = -inf;
%         t2(~s2 & s_inc) = inf;
%         tcone(~s_inc, ic) = max(t1(~s_inc), t2(~s_inc));
%         tcone(s_inc, ic) = min(t1(s_inc), t2(s_inc));
    end
%     if zsign>=0
%         
%     else
%         t1(~s1) = -inf;
%         t2(~s2) = -inf;
%         tcone(:,ic) = max(t1, t2);
%     end
end
tcone = tcone + u(:, 3);
tcone(tcone>1) = 1;
tcone(tcone<-1) = -1;

if ~s_zsign
    tcone = tcone(:, [2 1]);
end
dtc = tcone(:,2)-tcone(:,1);
dtc(dtc<=0) = nan;
% u
if isinverse
    u(:, 3) = (u(:, 3).*2 - sum(tcone,2))./dtc;
else
    u(:, 3) = (dtc.*u(:, 3) + sum(tcone,2))./2;
end

% w
w = dtc./2;

end
