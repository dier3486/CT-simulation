function mid_DFS = DFSmidchannel(mid_U, longorshot)
% hard code of DFS mode, should be configurable

if longorshot
    a = 3/8;
    b = -1/8;
else
    a = 1/8;
    b = -3/8;
end

La = floor(mid_U+a);
Da = mid_U+a - La;
Lb = floor(mid_U+b);
Db = mid_U+b - Lb;

mid_DFS = La + Lb + min(Da, Db)*2;

end