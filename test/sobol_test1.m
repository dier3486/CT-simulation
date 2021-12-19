P = sobolset(3);

m = 100;
N = 500;
t1 = zeros(1, 50);
t2 = zeros(1, 50);
for ii = 1:N
    x1 = P(1:m*ii, :);
    t1(ii) = sum(x1(:,1)>x1(:,2))./(m*ii);
    
    x2 = rand(m*ii, 3);
    t2(ii) = sum(x2(:,1)>x2(:,2))./(m*ii);  
end
