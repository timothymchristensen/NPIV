function H0 = H0fun(x,design)

if design == 1
    H0 = -2+4*x;
elseif design == 2
    H0 = log(abs(16*x-8)+1).*sign(16*x-8);
end
