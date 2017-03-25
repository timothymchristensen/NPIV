function [XX, DX, DDX, kts] = bsplinemat_quantile(x,m,r)

if m == 0
    kts = [zeros(1,r),ones(1,r)];
else
    kts = [zeros(1,r),quantile(x,(1:1:m)/(m+1)),ones(1,r)];
end
[XX, DX, DDX] = bsplinemat(x,m,r,kts);