img0 = head3{1};

Nimg = size(img0, 3);
Lb = 950;
Ub = 1150;

img1 = img0;
A1 = img0;
A1(A1>Ub) = Ub;
A1(A1<Lb) = Lb;
A1 = A1 - 1000;

Nv = 180;
v = (0:Nv-1).*(180/Nv);
ii = 187;
% for ii = 1:Nimg
    P1 = radon(A1(:,:,ii));
    P2 = P1(2:end-1, :).*2 - P1(3:end, :) - P1(1:end-2, :);
    P3 = mean(P2, 2);
%     P3  = P3 + flipud(P3);
%     P3((size(P3,1)+1)/2) = P3((size(P3,1)+1)/2)./2;
    P3 = repmat(P3, 1, Nv);
    B = iradon(P3, 0:179);
    B = B(2:end-1, 2:end-1);
    img1(:,:,ii) = img1(:,:,ii) - B;
% end
1;
% A1 = head3{1}(:,:,187)';
% A1 = A1-1000;
% A1(A1<0) = 0;
% A1(A1>100) = 0;
% imtool(A1)
A1 = head3{1}(:,:,187)'-1000;
A1(A1<0) = 0;
A1(A1>100) = 100;
% imtool(A1)
P1 = radon(A1);
% size(P1)
% imtool(P1)
P2 = P1(2:end-1, :).*2 - P1(1:end-2, :) - P1(3:end, :);
% imtool(P2)
% figure
% sum(P2,2)
% figure
% plot(sum(P2,2))
% 727/2
% grid on
% size(P2)
% B1=iradon(P2./2, 512);
% help iradon
% B1=iradon(P2./2, 0:179);
% size(B1)
% imtool(B1)
% imtool(A1-B1(2:end-1,2:end-1))
% imtool(A1+B1(2:end-1,2:end-1))
% imtool(A1-B1(2:end-1,2:end-1).*2)
% imtool(A1-B1(2:end-1,2:end-1).*1.5)
% close all hidden
% imtool(A1, [0 100])
% imtool(A1-B1(2:end-1,2:end-1).*1.8, [0 100])
P3 = repmat(mean(P2, 2), 1, 180);
size(P3)
B3=iradon(P3, 0:179);
% imtool(B3)
% imtool(A1-B3, [0 100])
% imtool(A1-B3(2:end, 2:end), [0 100])
% size(B3)
% imtool(A1-B3(2:end-, 2:end-), [0 100])
% imtool(A1-B3(2:end-1, 2:end-1), [0 100])
% imtool(A1)
% imtool(P1)
% imtool(P2)
% figure
% plot(mean(P2, 2))
% P3 = repmat(mean(P2, 2), 1, 180);