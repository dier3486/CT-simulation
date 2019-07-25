function D = linesinobject(A, B, objecttype)

N = size(A,1);
switch objecttype
    case 'sphere'
        Lab = sqrt(sum((A - B).^2, 2));
        d0 = pararea(A, B)./Lab;
        d1 = sqrt(sum(A.^2, 2) - d0.^2);
        d2 = 1 - d0.^2;
        d2(d2<0) = 0;
        d2 = sqrt(d2);

        L = [(d1 - d2)./Lab, zeros(N, 1)];
        R = [(d1 + d2)./Lab, ones(N, 1)];

    case 'cylinder'
        Lxy = sqrt(sum((A(:,1:2) - B(:,1:2)).^2, 2));
        S = abs(A(:,1).*B(:,2) - A(:,2).*B(:,1));
        d0 = S./Lxy;
        d1 = sqrt(sum(A(:,1:2).^2, 2) - d0.^2);
        d2 = 1 - d0.^2;
        d2(d2<0) = 0;
        d2 = sqrt(d2);
        L1 = (d1-d2)./Lxy;
        R1 = (d1+d2)./Lxy;

        LR2 = [-A(:,3), 1 - A(:,3)]./repmat(B(:,3) - A(:,3), 1, 2);
        sn = B(:,3)<A(:,3);
        LR2(sn, :) = fliplr(LR2(sn, :)); 

        L = [L1, LR2(:,1), zeros(N, 1)];
        R = [R1, LR2(:,2), ones(N, 1)];
        
    case 'blade'
        LR = [-A(:,3), 1 - A(:,3)]./repmat(B(:,3) - A(:,3), 1, 2);
        sn = B(:,3)<A(:,3);
        LR(sn, :) = fliplr(LR(sn, :)); 

        L = [LR(:,1), zeros(N, 1)];
        R = [LR(:,2), ones(N, 1)];
        
    otherwise
        D = zeros(N, 1);
        return 
end

D = multiinsect(L, R);
return